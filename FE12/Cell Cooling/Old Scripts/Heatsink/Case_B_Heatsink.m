current = 13:0.26:52;
IR = 0.0135; % Ohms
radius = 10.775/1000; % meters
height = 70.15/1000; % meters
Tmax = 60 + 273.15; % Kelvin, set this at the center of the axial axis, so z(max) = 1/2*cell height
Tambient = 30 + 273.15; % Kelvin
Rpaste = 0.001; %K/W
axialthermalconductivity = 20.61; % W/mK

heat_flux = zeros(length(current), 1);
h = zeros(length(current), 1);
Tsurface = zeros(length(current), 1);
Tsurface_withconv = zeros(length(current), 1);
qconv = zeros(length(current), 1);
Rcell = zeros(length(current), 1);
Rfins = zeros(length(current), 1);

for i = 1:length(current)
    heat_flux(i) = ((current(i))^2*IR)/(pi*(radius)^2); % is acutally a volume so need to multiply by height/2
    %heat_flux(i) = ((current(i))^2*IR)/(height/2*radius); %area as a rectangular region
    
    Tsurface(i) = Tmax - heat_flux(i)/(2*axialthermalconductivity)*((height/2)^2); % Kelvin
    h(i) = (heat_flux(i)*radius)/(Tsurface(i)-Tambient);
    %qconv(i) = (h(i)*2*pi*radius^2+2*pi*radius)*(Tsurface(i)-Tambient);
    %Tsurface_withconv(i) = Tmax - (heat_flux(i)-qconv(i))/(2*axialthermalconductivity)*((height/2)^2); % Kelvin

    %Rfins = (n/Rfin + 1/Rbase)^-1; %where n is number of fins
    Rfins(i) = (Tmax - Tambient)/(current(i)^2*IR)-Rpaste; % necessary thermal resistance of heatsink and base
end
% display(Rcell);
% format long;

figure
plot(Tsurface-273.15,Rfins);
xlabel('Surface Temperature (degC)');
ylabel('Fin Thermal Resistance');
title('Rfin vs Surface Temperature');

figure
plot(current,Rfins);
xlabel('Curent (A)');
ylabel('Fin Thermal Resistance');
title('Rfin vs Current');

figure
plot(Tsurface-273.15,h);
xlabel('Surface Temperature (degC)');
ylabel('Heat Transfer Coefficient (W/m2K)'); 
title('Heat Transfer Coefficient vs Surface Temperature');

% t = thickness;
% w = width;
% L = length;
% Lc = L + t/2;
% P = perimeter;
% k = thermal conductivity of heatsink material;
% N = number of fins;
% A_c = fin cross sectional area;
% A_b = base area exposed to convection;
% A_t = A_b + N*A_f; total surface area exposed to convection
% theta = Tsurface - Tambient; % Tsurface at set cell temp, Tsurface above
% 
% % for a straight, rectangular fin
% A_p = t*L;
% A_f = 2*w*Lc;
% mLc = sqrt((2*h)/(k*A_p))*Lc^(3/2);
% 
% q_f = (sqrt(h*P*k*A_c)*theta*tanh(mLc)); % fin heat transfer rate
% fin_efficiency = q_f / (h*A_f*theta);
% 
% Rfins = theta/qf;
% Rbase = 1/(h*A_c);
% fin_effectiveness = Rbase/Rfins;
% 
% q_t = N*(fin_efficiency*h*A_f*theta)+(h*A_b*theta);
% overall_fin_efficiency = q_t/(h*A_t*theta);


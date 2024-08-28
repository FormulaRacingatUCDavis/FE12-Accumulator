current = 6.5:0.135:26;
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
    %heat_flux(i) = ((current(i))^2*IR)/(pi*(radius)^2);
    heat_flux(i) = ((current(i))^2*IR)/(height/2*radius); %area as a rectangular region
    
    Tsurface(i) = Tmax - heat_flux(i)/(2*axialthermalconductivity)*((height/2)^2); % Kelvin
    h(i) = (heat_flux(i)*radius)/(Tsurface(i)-Tambient);
    %qconv(i) = (h(i)*2*pi*radius^2+2*pi*radius)*(Tsurface(i)-Tambient);
    %Tsurface_withconv(i) = Tmax - (heat_flux(i)-qconv(i))/(2*axialthermalconductivity)*((height/2)^2); % Kelvin

    Rcell(i) = (Tmax - Tsurface(i))/((current(i))^2*IR);
    Rfins(i) = (Tmax - Tambient)/(current(i)^2*IR)-Rcell(i)-Rpaste;
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
title('Heat Transfer Coefficient vs Temperature');


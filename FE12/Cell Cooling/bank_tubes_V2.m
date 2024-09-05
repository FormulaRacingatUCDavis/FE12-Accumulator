%spacing is about 1.45 mm (random spacing right now), 1-5 mm should be fine 
% we want the transverse spacing (plus diameter/2-which is variable A1) to be greater than the diagonal spacing so
% that Vmax is at A1

diameter = 21.55/1000; % meters
axialthermalconductivity = 2.21; % W/mK
radius = 10.775/1000; % meters
height = 70.15/1000; % meters
SurfaceArea = 2*pi*radius^2+2*pi*radius*height; % m^2

Tambient = 35 + 273.15; % Kelvin
Tmax = 50 + 273.15; % Kelvin

current = 13:0.26:52; 
resistance = 0.0135; % Ohms
Ohmic_heat = zeros(length(current), 1);
heatflux = zeros(length(current), 1);
Tsurface = zeros(length(current), 1);
Rconv = zeros(length(current), 1);
heattransfer_coefficient = zeros(length(current), 1);
Nu = zeros(length(current), 1);
max_Reynolds = zeros(length(current), 1);
V_max = zeros(length(current), 1);
V_inlet = zeros(length(current), 1);

%Air Properties at 30 C
rho = 1.164; % [kg/m^3]
Cp = 1006; % [J/kg K]
mu = 18.72*10^-6; % [N s/m^2]
k = 26.62*10^-3; % [W/m K] 
Pr = mu*Cp/k; %

%Air Properties at 50 C (Surface Temp)
rho_s = 1.093; % [kg/m^3]
Cp_s = 1007; % [J/kg K]
mu_s = 19.53*10^-6; % [N s/m^2]
k_s = 28.08*10^-3; % [W/m K]
Pr_s = mu_s*Cp_s/k_s; %


%spacing 
spacing_transverse = 24/1000; % meters
spacing_longitudinal = 23/1000; % meters
spacing_diagonal = sqrt(spacing_longitudinal^2 + spacing_longitudinal^2); %[m] diagonal spacing
spacing_ratio = spacing_transverse/spacing_longitudinal;

% if spacing_ratio < 2
%     C1 = 0.35*(spacing_ratio)^(1/5);
%     m = 0.6;
% else 
%     C1 = 0.4;
%     m = 0.6;
% end
C1 = 0.9;
C2 = 0.935; % 6 rows
m = 0.4;


for i = 1:length(current)
    Ohmic_heat(i) = current(i)^2*resistance;
    heatflux(i) = Ohmic_heat(i)/(pi*radius^2);

    Tsurface(i) = Tmax - heatflux(i)*radius^2/4*axialthermalconductivity; % Kelvin
    Rconv(i) = (Tsurface(i)-Tambient) / Ohmic_heat(i);
    heattransfer_coefficient(i) = 1/(SurfaceArea*Rconv(i)); %average

    Nu(i) = heattransfer_coefficient(i) * diameter / axialthermalconductivity;
    max_Reynolds(i) = (Nu(i) / ((C1*C2)*(Pr^(0.36))*(Pr/Pr_s)^(0.25)) )^(1/m);
    
    %Checking where max velocity occurs
    A1 = (spacing_transverse-diameter);
    A2 = 2*(spacing_diagonal-diameter);
    if A2<A1
        V_max(i) = max_Reynolds(i)*mu/(rho*2*(spacing_diagonal-diameter)); %Vmax is at A2
        V_inlet(i) = V_max(i)*2*(spacing_diagonal-diameter)/spacing_transverse;
    else
        V_max(i) = max_Reynolds(i)*mu/(rho*(spacing_transverse-diameter)); %Vmax is at A1
        V_inlet(i) = V_max(i)*(spacing_transverse-diameter)/spacing_transverse;
    end

    %Pressure Drop
    P = current*432/1000;
    N_rows = 6;
    P_l = spacing_longitudinal/diameter;
    P_t = spacing_transverse/diameter;
    P_ratio = P_t/P_l;
    %Rough estimates of pressure drop coefficient
    X=1.1; 
    f = 0.9;
    P_delta(i) = N_rows*X*(rho*V_max(i)^2/2)*f; %[Pa]
end

figure()
subplot(2,1,1)
plot(P,V_inlet)
ylabel('Inlet Velocity [m/s]')
xlabel('Average Power Output [kW]')

subplot(2,1,2)
plot(P,P_delta)
ylabel('Delta P [Pa]')
xlabel('Average Power Output [kW]')

diameter = 21.55; %mm
velocity = V_inlet; %initial air velocity m/s

spacing_transverse = 22:0.5:28; %mm, spacing is from center of one cell to center of another cell
spacing_longitudinal = 22:0.5:28; %mm
spacing_diagonal = zeros(length(spacing_longitudinal), 1);

Vmax_table = cell(length(spacing_longitudinal), length(spacing_transverse));
Vmax_location_table = cell(length(spacing_longitudinal), length(spacing_transverse));

%Check Vmax direction
for i = 1:length(spacing_longitudinal)
    for j = 1:length(spacing_transverse)
        A1 = (spacing_transverse(j)+diameter)/2;
        spacing_diagonal(j) = sqrt((spacing_transverse(j)/2)^2 + spacing_longitudinal(i)^2);

        if spacing_diagonal<A1
            Vmax = spacing_transverse(j)/(2*(spacing_diagonal(j)-diameter)) * velocity; %Vmax is at A2
            %fprintf('Vmax is at A2\n');
            Vmax_table{i,j} = Vmax;
            Vmax_location_table{i, j} = 'A2';
        else
            Vmax = spacing_transverse(j)/(spacing_transverse(j)-diameter) * velocity; %Vmax is at A1
            %fprintf('Vmax is at A1\n');
            Vmax_table{i,j} = Vmax;
            Vmax_location_table{i, j} = 'A1';
        end
    end
end





% Pressure drop 
% Tsurface = 50;
% Tinitial = 30;
% rows = 6;
% tubes_per_row = 12;
% totaltubes = rows*tubes_per_row;
% 
% Ts_minus_To = (Tsurface-Tinitial)*exp((-pi*diameter/1000*totaltubes*h)/(rho*velocity*tubes_per_row *spacing_transverse/1000*cp*1000));
% 
% logmeantemp = (Tsurface-Tinitial-Ts_minus_To)/log((Tsurface-Tinitial)/Ts_minus_To); %degC, difference in inlet and outlet air temp
% 
% heattransferrate= totaltubes*h*pi*diameter/1000*logmeantemp/1000; %kW/m

% %pressure drop
% PT = spacing_transverse/diameter;
% PL = spacing_longitudinal/diameter;
% Pratio = PT/PL;
% correction = 1; % This changes based on Pratio and PT
% frictionfactor = 1; % This changes based on Pratio and PT
% 
% deltaP = rows*correction*(rho*Vmax^(2)/2)*frictionfactor; %Pascals


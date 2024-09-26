%% Bank of Tubes Analysis on Subpack Cooling
clear all 
% close all
%% Assumptions
%Nusselt Number correlation does not technically hold true since our flow
%regieme changes as a function of our velocity, however it is a decent
%enough of an approximation

%Average internal resistance of cell used
%Entropic heat not accounted for

%Air properties only evaluated at 30C and 50C
%% Parameters
P_avg = 18.6; % [kW] Average power 
P_max = 20; % [kW] Max power 


%Cell Spacing (from center of one cell to another)
% s_t = 23*10^-3; % [m] traverse spacing
s_l = 24*10^-3; % [m] longitudinal spacing
s_t = 23e-3; 
%s_t = 23e-3: 0.5e-3 : 24e-3; % [m] longitudinal spacing
%s_d = sqrt(s_l^2 + s_t^2); %[m] diagonal spacing

%Constants Based on Cell Config 
%C1= 0.35*(s_t/s_l)^(1/5);
%C1 = 0.9;
C2 = 0.95;
m = 0.6;

%Cell Properties
T_cell_max = 55; %[C] 
T_ambient = 35; %[C]
N_rows = 12; % Number of rows
N_cells = 6; % Number of cells per row
totalcells = N_rows*N_cells;
D_cell = 21.55 *10^-3; %[m] Cell diameter
R_cell = D_cell/2; % [m] Cell radius
k_cell = 2.21; %[W /m K] 
L_cell = 70.15*10^-3; % [m] Cell length/height
SA_cell = 0.5*(L_cell-12.7*10^-3)*2*pi*R_cell; %[m^2] Surface Area of cell
R_internal = 0.0087; % [Ohms]

% %Air Properties at 30 C
% rho = 1.164; % [kg/m^3]
% Cp = 1006; % [J/kg K]
% mu = 18.72*10^-6; % [N s/m^2]
% k = 26.62*10^-3; % [W/m K] 
% Pr = mu*Cp/k; %

%Air Properties at 35 C
rho = 1.1459; % [kg/m^3]
Cp = 1006.7; % [J/kg K]
mu = 18.915*10^-6; % [N s/m^2]
k = 26.71*10^-3; % [W/m K] 
Pr = mu*Cp/k; %

% %Air Properties at 50 C (Surface Temp)
% rho_s = 1.093; % [kg/m^3]
% Cp_s = 1007; % [J/kg K]
% mu_s = 19.53*10^-6; % [N s/m^2]
% k_s = 28.08*10^-3; % [W/m K]
% Pr_s = mu_s*Cp_s/k_s; 

%Air Properties at 55 C (Surface Temp)
rho_s = 1.0759; % [kg/m^3]
Cp_s = 1007.7; % [J/kg K]
mu_s = 19.835*10^-6; % [N s/m^2]
k_s = 28.16*10^-3; % [W/m K]
Pr_s = mu_s*Cp_s/k_s; %


%Pack Properties
V_nom = 432; % [V]



%% Calculations

% For s_l = 24 mm

%Matrix Initialization
n = 40; % number of intervals
P = 0:P_avg/n:P_max;
V_inlet = size(P);
P_delta = size(P);

for j = 1:length(s_t)
    for i = 1:length(P)
        P_avg = P(i);
        %Heat generation and temperature calcs
        I = (P_avg*1000/V_nom)/3; %[A]
        q_total = I^2*R_internal; % [W]
        q_vol_gen = q_total/(pi*R_cell^2*L_cell); % [W/m^3]
        T_s = T_cell_max-(q_vol_gen*R_cell^2)/(4*k_cell); % [K]
        
        %Heat transfer coefficient 
        h= q_total/((T_s-T_ambient)*SA_cell);
        
        %Flow characteristics
        C1 = 0.35*(s_t(j)/s_l)^(1/5);
        Nu = h*D_cell/k_s;
        Re_max = (Nu/(C1*C2*Pr^0.36*(Pr/Pr_s)^0.25))^(1/m);
        
        %Checking where max velocity occurs

        s_d(j) = sqrt(s_l^2 + (s_t(j)/2)^2); %[m] diagonal spacing

        A1 = (s_t(j)-D_cell);
        A2 = 2*(s_d(j)-D_cell);
        if A2<A1
            V_max = Re_max*mu/(rho*2*(s_d(j)-D_cell)); %Vmax is at A2, denominator is just diameter
            V_inlet(j,i) = V_max*2*(s_d(j)-D_cell)/s_t(j);
        else
            V_max = Re_max*mu/(D_cell); %Vmax is at A1
            V_inlet(j,i) = V_max*(s_t(j)-D_cell)/s_t(j);
        end
        
        %Pressure Drop
        
        P_l = s_l/D_cell;
        P_t = s_t(j)/D_cell;
        Pans = (P_t-1)/(P_l-1);
        P_ratio = P_t/P_l;

        %Rough estimates of pressure drop coefficient
        X= 1.1; % staggered graph
        f = 0.9; % staggered graph
        P_delta(j,i) = N_rows*X*(rho*V_max^2/2)*f; %[Pa]
    
    end
end

%% Plotting
figure()
subplot(2,1,1)
plot(V_inlet*0.01018*60,P);
xlabel('Inlet Velocity [m3/min]')
ylabel('Average Power Output [kW]')
legend('23','23.5','24','24.5','25','25.5','26','26.5','27')

subplot(2,1,2)
plot(P,P_delta)
ylabel('Pressure Drop [Pa]')
xlabel('Average Power Output [kW]')
legend('23','23.5','24','24.5','25','25.5','26','26.5','27')
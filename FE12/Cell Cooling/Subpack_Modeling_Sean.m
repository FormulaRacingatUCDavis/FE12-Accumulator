%% Bank of Tubes Analysis on Subpack Cooling
clear all 
close all
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
s_l = 11.5*10^-3; % [m] longitudinal spacing
%s_t = 23e-3: 0.5e-3 : 24e-3; % [m] longitudinal spacing
% s_l = 23*10^-3 : 0.5e-3 : 24e-3; % [m] traverse spacing
s_t = 24e-3;%:0.5e-3:24e-3; % [m] longitudinal spacing
%s_l = 23*10^-3 : 0.5e-3 : 27e-3; % [m] longitudinal spacing
%s_d = sqrt(s_l^2 + s_t^2); %[m] diagonal spacing

%Constants Based on Cell Config 
% C1= 0.35*(s_t/s_l)^(1/5);
% C1 = 0.9;
% C2 = 0.99;
% m = 0.4;

%Single Cylinder in Crossflow configs
C= 0.683;
m = 0.466;
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
Cp_cell = 1360; %[J/kgK]  

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

%%
%Assume lumped thermal mass 
%From CFD assume we know air flow characteristics 
V_max = 3.5; %[m/s]
D = s_t-D_cell;
Re = rho*V_max*D/mu;
Nu = C*(Re^m) * (Pr^(1/3));
% Nu = (C1*(Re^m) * (Pr^0.36) * (Pr/Pr_s)^.25)/C2;
h = Nu*k_cell/D_cell;

R = 6; %number of cells per row scaling factor
%Heat generation and temperature calcs
I = (P_avg*1000/V_nom)/3; %[A]
q_gen = I^2*R_internal; % [W] 


T_surface = zeros(1,N_rows);
T_inlet_air = zeros(1,N_rows+1); %First index is at N=0 for no 
T_inlet_air(1) = 35; %[C]

for N = 1:N_rows
    T_surface(N) = R*q_gen/(h*SA_cell) + T_inlet_air(N);
    x = (-pi * D_cell * N * h)/(rho*V_max*R*N_rows*s_t*Cp_cell);
    T_outlet = -(q_gen/(h*SA_cell)+T_inlet_air(N)-T_inlet_air(1))*exp(x) + T_inlet_air(1) +q_gen/(h*SA_cell);
    T_inlet_air(N+1) = T_outlet;
end

N_plot = 0:N_rows;
figure;
hold on
plot(N_plot(2:end),T_surface,'DisplayName','Cell Surface Temp [°C]')
plot(N_plot,T_inlet_air,'DisplayName','Inlet Air Temp [°C]')
% plot(N_plot(2:end),T_inlet_air(2:end),'DisplayName','Outlet Air Temp [°C]')
legend('Location','best')
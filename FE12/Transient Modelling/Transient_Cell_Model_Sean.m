%% Transient Cell Modelling Script
% clear all
% close all
% clc
%% Assumptions
%Lumped Capacitance for all cells
%Constant properties (i.e. Internal Resistance or Air properties)
%Convective heat transfer between cells and air is approximated as hA(T_cell-T_air_inlet)
%Only 50% of SA of cell is exposed to convective heat transfer

%% Parameters

%Air Properties at 35 C
rho_air = 1.1459; % [kg/m^3]
Cp_air = 1006.7; % [J/kg K]
mu_air = 18.915*10^-6; % [N s/m^2]
k_air = 26.71*10^-3; % [W/m K] 
Pr = mu_air*Cp_air/k_air; %

%Cell Properties
T_ambient = 35; %[C]
N_cells = 12; % Number of cells per row
D_cell = 21.55 *10^-3; %[m] Cell diameter
R_cell = D_cell/2; % [m] Cell radius
k_cell = 2.21; %[W /m K] Radial thermal conductivity
L_cell = 70.15*10^-3; % [m] Cell length/height
SA_cell = 0.5*(L_cell-12.7*10^-3)*2*pi*R_cell; %[m^2] Surface Area of cell
R_internal = 0.0087; % [Ohms]

%Cell Spacing (from center of one cell to another)
s_l = 11.5*10^-3; % [m] longitudinal spacing
s_t = 24e-3;

%Pack Properties
V_nom = 432; % [V]

%% Calculations 

u_air = 4.5; %[m/s] Air velocity (set by fans)

N = 1000; %Number of steps
time = linspace(0,30*6,N); %[s]
P = 18.6*ones(size(time)); %[kW]

%Finding convective htc [h]
s_d = sqrt(s_l^2 + (s_t)^2); %[m] diagonal spacing
A1 = (s_t-D_cell);
A2 = 2*(s_d-D_cell);
if A2<A1
    Re_max = rho_air*u_air*A2/mu_air;
else
    Re_max = rho_air*u_air*A1/mu_air;
end
%Re_max appears to fall within 10^2 - 10^3 range, -> able to approximate as
%single cylinder in cross flow
C = 0.683;
m = 0.466;
Nu_average = C*(Re_max^m)*(Pr^(1/3));
h = Nu_average*k_air/D_cell; %[W/m^2 K]

Bi = h*D_cell/k_cell; 
%Lumped Capacitance assumption gets less and less accurate the higher the
%inlet air velocity

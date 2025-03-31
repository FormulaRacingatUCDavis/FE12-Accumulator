s_l = 11.5*10^-3; % [m] longitudinal spacing
s_t = 24e-3;

%Cell Properties
T_cell_max = 55; %[C] 
T_ambient = 35; %[C]
N_rows = [0,1,2,3,5,6,7,8,9,10,11,12]; % Number of rows
N_cells = 6; 
D_cell = 21.55 *10^-3; %[m] Cell diameter
R_cell = D_cell/2; % [m] Cell radius
k_cell = 2.21; %[W /m K] 
L_cell = 70.15*10^-3; % [m] Cell length/height
SA_cell = 0.5*(L_cell-12.7*10^-3)*2*pi*R_cell; %[m^2] Surface Area of cell
R_internal = 0.0087; % [Ohms]

P = 18.6;
Vnom = 432;

% 50 C
rho = 1.1459; % [kg/m^3]
Cp = 1006.7; % [J/kg K]
mu = 18.915*10^-6; % [N s/m^2]
k = 26.71*10^-3; % [W/m K] 
Pr = mu*Cp/k;

h = 1000
current = (P_avg*1000/V_nom)/3;
qgen = 6*current**2*R_internal/(h*A)

rows = [1,2,3,4,5,6,7,8,9]
Tinlet(1) = 30

for i = 1:length(rows):
    Tsurface(i) = qgen + Tinlet(i);
    Tout(i) 
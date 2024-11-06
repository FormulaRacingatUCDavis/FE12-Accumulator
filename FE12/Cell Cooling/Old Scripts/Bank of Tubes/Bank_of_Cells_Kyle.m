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

%Cell Spacing (from center of one cell to another)
s_t = 23*10^-3; % [m] traverse spacing
s_l = 24*10^-3; % [m] longitudinal spacing
s_d = sqrt(s_l^2 + s_t^2); %[m] diagonal spacing

%Constants Based on Cell Config (assuming we are 
C1= 0.35*(s_t/s_l)^(1/5);
%C1 = 0.9;
C2 = 0.95;
m = 0.6;
%m = 0.4

%Cell Properties
T_cell_max = 50; %[C] 
T_ambient = 35; %[C]
N_rows = 12; % Number of rows
N_cells = 6; % Number of cells per rw
D_cell = 21.55 *10^-3; %[m] Cell diameter
R_cell = D_cell/2; % [m] Cell radius
%R_cell = 12.93e-3; % this is gerater than radius of the cell, larger
%surface area = lower pressure drop
k_cell = 2.21; %[W /m K] 
L_cell = 70.15*10^-3; % [m] Cell length/height
SA_cell = 0.5*(L_cell-12.7*10^-3)*2*pi*R_cell; %[m^2] Surface Area of cell
R_internal = 0.0135061; % [Ohms]

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


%Pack Properties
V_nom = 432; % [V]



%% Calculations

%Matrix Initialization
n = 20; % number of intervals
P = 0:P_avg/n:P_avg;
V_inlet = size(P);
P_delta = size(P);

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
    Nu = h*D_cell/k_s;
    Re_max = (Nu/(C1*C2*Pr^0.36*(Pr/Pr_s)^0.25))^(1/m);
    
    %Checking where max velocity occurs
    A1 = (s_t-D_cell);
    A2 = 2*(s_d-D_cell);
    if A2<A1
        V_max = Re_max*mu/(rho*2*(s_d-D_cell)); %Vmax is at A2, denominator is just diameter
        V_inlet(i) = V_max*2*(s_d-D_cell)/s_t;
    else
        V_max = Re_max*mu/(D_cell); %Vmax is at A1
        V_inlet(i) = V_max*(s_t-D_cell)/s_t;
    end
    
    %Pressure Drop
    
    P_l = s_l/D_cell;
    P_t = s_t/D_cell;
    Pans = (P_t-1)/(P_l-1);
    P_ratio = P_t/P_l;
    %Rough estimates of pressure drop coefficient
    X=1.1; % staggered graph
    f = 0.9; % staggered graph
    %X = 3; % aligned graph; hypothesis was that spacing was almost the
    %same and pratically was like the aligned cells format, but these gave
    %different results
    %f = 0.7; % aligned graph
    f = 0.9;
    P_delta(i) = N_rows*X*(rho*V_max^2/2)*f; %[Pa]
end

totalcells = N_rows*N_cells;

Ts_minus_To(i) = (Tsurface-Tinitial)*exp((-pi*D_cell/1000*totalcells*h)/(rho*V_inlet(i)*N_cells*s_t/1000*Cp*1000));

logmeantemp(i) = (Tsurface-Tinitial-Ts_minus_To(i))/log((Tsurface-Tinitial)/Ts_minus_To); %degC, difference in inlet and outlet air temp

heattransferrate= totaltubes*h*pi*diameter/1000*logmeantemp/1000; %kW/m

%% Plotting
figure()
subplot(2,1,1)
plot(V_inlet*0.01018*60,P);
xlabel('Inlet Velocity [m3/min]')
ylabel('Average Power Output [kW]')

subplot(2,1,2)
plot(P,P_delta)
ylabel('Delta P [Pa]')
xlabel('Average Power Output [kW]')
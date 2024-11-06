%spacing is about 1.45 mm (random spacing right now), 1-5 mm should be fine 
% we want the transverse spacing (plus diameter/2-which is variable A1) to be greater than the diagonal spacing so
% that Vmax is at A1

spacing_transverse = 24; %mm, spacing is from center of one cell to center of another cell
spacing_longitudinal = 23; %mm
spacing_diagonal = sqrt((spacing_transverse/2)^2 + spacing_longitudinal^2);
diameter = 21.55; %m
velocity = 0.3; %initial air velocity m/s
C2 = 0.935; % 6 rows

%Check Vmax direction
A1 = (spacing_transverse+diameter)/2;
if spacing_diagonal<A1
    Vmax = spacing_transverse/(2*(spacing_diagonal-diameter)) * velocity; %Vmax is at A2
    fprintf('Vmax is at A2');
else
    Vmax = spacing_transverse/(spacing_transverse-diameter) * velocity; %Vmax is at A1
    fprintf('Vmax is at A1');
end

%air at 30 C
rho = 1.1614; %kg/m3
cp = 1.007; %kJ/kg K
viscosity = 15.89e-6; %m2/s
k = 26.3e3; %W/mK
Pr = 0.707; %at 30 C
Prs = 0.702; %at 60 C

diameter = 21.55/1000;

max_Reynolds = (Vmax*diameter)/(viscosity);
spacing_ratio = spacing_transverse/spacing_longitudinal;

C1 = 0.35*(spacing_transverse/spacing_longitudinal)^(1/5);
m = 0.60;

Nu = C2*C1*(max_Reynolds^(m))*(Pr^(0.36))*(Pr/Prs)^(0.25);
h = Nu*k/diameter;

Tsurface = 60;
Tinitial = 30;
rows = 6;
tubes_per_row = 12;
totaltubes = rows*tubes_per_row;

Ts_minus_To = (Tsurface-Tinitial)*exp((-pi*diameter*totaltubes*h)/(rho*velocity*tubes_per_row *spacing_transverse/1000*cp*1000));

logmeantemp = (Tsurface-Tinitial-Ts_minus_To)/log((Tsurface-Tinitial)/Ts_minus_To); %degC, difference in inlet and outlet air temp

heattransferrate= totaltubes*h*pi*diameter*logmeantemp/1000; %kW/m
%heat transfer rate per unit length of the cell, match this to power per
%cell for cooling q' = Q/L*deltaT where Q is power in kW

%pressure drop
PT = spacing_transverse/diameter*1000;
PL = spacing_longitudinal/diameter*1000;
Pratio = PT/PL;
correction = 1;
frictionfactor = 1; 

deltaP = rows*correction*(rho*Vmax^(2)/2)*frictionfactor; %Pascals


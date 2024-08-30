% 298 K, 303 K, 313 K; NEXT data point
% 20 30 40 50 60 70 90 100 - 80% soc failed to fit for 298 K and 303 K

OCV = [3.388,  3.396, 3.405; 3.495, 3.502, 3.508; 3.578, 3.586, 3.600; 3.659, 3.667, 3.678; 3.752, 3.760, 3.771; 3.846, 3.854, 3.863; 4.006, 4.013, 4.023; 4.190, 4.119, 4.129];
%R0 = [1.015e-2, 8.735e-3, 7.051e-3; 8.854e-3,8.148e-3, 7.010e-3;8.138e-3,7.448e-3, 6.391e-3; 7.944e-3,7.183e-3,6.458e-3; 6.761e-3, 6.167e-3, 5.510e-3; 8.941e-3, 8.448e-03,7.606e-3; 5.416e-3, 4.632e-3, 3.615e-3; 2.252e-2, 1.067e-2,9.876e-3];
R1 = [1.553e-3, 1.042e-3, 7.357e-4; 1.062e-3, 6.767, 3.742e-4; 2.269e-3, 2.303e-3,1.291e-3; 1.84e-3,2.082e-3,1.348e-3; 1.079e-3, 8.499e-4, 8.860e-4; 5.908e-4, 5.991e-4,3.182e-4; 2.786e-3, 1.953e-3, 9.038e-4; 8.237e-4, 5.795e-3, 2.273e-4];
R2 = [8.115e-3, 7.461e-3, 7.171e-3; 1.548e-3, 1.485e-3, 1.687e-3; 4.280e-3, 2.185e-3, 9.197e-4; 4.646e-3, 3.258e-3, 3.128e-3; 2.078e-3, 1.937e-3, 1.433e-3; 3.664e-3, 3.806e-3,3.628e-3; 1.527e-3, 1.281e-3, 9.251e-4; 1.050e-3, 8.934e-4, 6.163e-4];
T1 = [4.488e1, 2.941e1, 2.314e1; 6.333e1, 5.184e1, 3.609e1; 9.440e1, 9.948e1, 8.308e1; 8.794e1, 1.035e2, 5.770e1; 6.438e1, 5.881e1, 9.427e1; 4.297e1, 5.109e1, 2.136e1; 9.776e1, 6.897e1, 6.713e1; 6.504e1, 7.342e1, 3.640e1];
T2 = [5.608e2, 5.613e2, 5.924e2; 5.611e2, 5.187e2, 5.904e2; 3.910e2, 4.241e2, 4.190e2; 4.377e2, 4.738e2, 4.967e2; 4.217e2, 4.547e2, 5.356e2; 5.026e2, 5.757e2, 5.195e2; 4.093e2, 4.284e2, 5.288e2; 5.518e2, 6.160e2,3.445e2];
RevHeat_graph = [0.147601536, -0.096130371,0.265421186, 0.334903172,0.381578718, 0.219563075,0.14163426, 0.092015948];
%  RevHeat is the coefficients without 10 and 80 SOC since ansys didn't
%  give those, size of arrays have to match for the battery block in
%  simulink

%R0 = [0.016,0.016,0.016;0.016,0.016,0.016;0.016,0.016,0.016;0.016,0.016,0.016;0.016,0.016,0.016;0.016,0.016,0.016;0.016,0.016,0.016;0.016,0.016,0.016];

%R0 manually found is below
R0 = [0.014784,0.013498,0.011892;0.014551,0.013293,0.011714;0.014481,0.013214,0.011704;0.014439,0.013259,0.01174;0.014484,0.013295,0.011795;0.014594,0.013376,0.011867;0.014876,0.01353,0.011846;0.016217,0.014558,0.012716];


%entropic heat coefficients
% 313K = 40 C, 303K = 30 C, 298K = 25 C

coefficients = [1.025390625, 0.147601536, -0.096130371, 0.265421186, 0.334903172,0.381578718, 0.219563075, 0.213922773, 0.14163426, 0.092015948];
% mV/K

figure
x = linspace(10,100,10);
xfit = linspace(10,101,92);
pp = spline(x,coefficients);
interpolated = fnval(pp,xfit);
plot(x,coefficients,'o',xfit,interpolated);
xlim([0,100]);
ylabel('Entropic Coefficient (mV/K)');
xlabel('SOC (%)');
title('Entropic Coefficient vs SOC');

coeffdata = [];
for i = 1:91
    coeffdata(i+9,1) = round(xfit(1,i),0);
    coeffdata(i+9,2) = interpolated(1,i);
end

fulldata = xlsread("FE11 Endurance Full Data V2.xlsx");
revheat = zeros(length(fulldata), 1);
soc_temporary = fulldata(:,7);
j = 0;

new_current = zeros(length(fulldata), 1);
new_temp = zeros(length(fulldata), 1);
current_temporary = fulldata(:,3);
temp_temporary = fulldata(:,4);

for i = 1:2:length(current_temporary)
    if i < 9238
        new_current(i,1) = current_temporary(i) + current_temporary(i+1);
    end
end

for i = 1:length(fulldata)
    value = soc_temporary(i);
    value = int64(value);
    if value == 0
        continue;
    end
    if value > 100
        continue;
    end
    if value == coeffdata(value,1) 
        revheat_pack(i,1) = new_current(i,1) * (fulldata(i,4)+273.15) * (coeffdata(value,2)/1000);
        revheat_cell(i,1) = new_current(i,1)/4 * (fulldata(i,4)+273.15) * (coeffdata(value,2)/1000); %current*temp*dOCV/dT, 
        % current divided by 4 to get the reversible heat of one cell
    end
end

revheat_time = [];
for i = 1:length(revheat_cell)
    revheat_time(i,1) = fulldata(i,1);
    revheat_time(i,2) = revheat_cell(i,1);
end
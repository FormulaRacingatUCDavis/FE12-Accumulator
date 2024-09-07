% Temperatures: 303 K, 308 K, 313 K, 318 K, 323 K

params2C = readmatrix("P45B_Parameters_2C_spreadsheet.xlsx");
params4C = readmatrix("P45B_Parameters_4C_spreadsheet.xlsx");

VOC_2C = zeros(5);
R1_2C = zeros(5);
R2_2C = zeros(5);
C1_2C = zeros(5);
C2_2C = zeros(5);
R0_2C = zeros(5);

VOC_4C = zeros(5);
R1_4C = zeros(5);
R2_4C = zeros(5);
C1_4C = zeros(5);
C2_4C = zeros(5);
R0_4C = zeros(5);

a = 1;
b = 2; 
c = 3; 
d = 4;
e = 5;

for k = 1:5
    if k >= 2
        a = a + 5;
        b = b + 5;
        c = c + 5;
        d = d + 5;
        e = e + 5;
    end

    for i = 2:9
        p = i - 1;
        VOC_2C(p,k) = params2C(i,a);
        R1_2C(p,k) = params2C(i,b);
        C1_2C(p,k) = params2C(i,c);
        R2_2C(p,k) = params2C(i,d);
        C2_2C(p,k) = params2C(i,e);
    end

    for i = 2:11
        p = i - 1;
        VOC_4C(p,k) = params4C(i,a);
        R1_4C(p,k) = params4C(i,b);
        C1_4C(p,k) = params4C(i,c);
        R2_4C(p,k) = params4C(i,d);
        C2_4C(p,k) = params4C(i,e);
    end
end

R0data_2C = xlsread("HPPC_2C_IR_spreadsheet.xlsx")
R0data_4C = xlsread("HPPC_4C_IR_spreadsheet.xlsx")

for i = 1:8
    for j = 1:5
        R0_2C(i,j) = R0data_2C(i+2,j)/1000;
    end
end

for i = 1:10
    for j = 1:5
        R0_4C(i,j) = R0data_4C(i,j)/1000;
    end
end

T1_2C = R1_2C.*C1_2C;
T2_2C = R2_2C.*C2_2C;

T1_4C = R1_4C.*C1_4C;
T2_4C = R2_4C.*C2_4C;

avgIR_2C = 8.19999695/1000;
avgIR_4C = 8.75918494/1000;


%2C starts at 30% SOC
coefficients2C =[-0.08430478,-0.5011749,0.19809724,0.1816177,0.19474032,0.10135652,0.11390688,0.01502992];

figure
x = linspace(30,100,8);
xfit = linspace(30,101,72);
pp = spline(x,coefficients2C);
interpolated = fnval(pp,xfit);
plot(x,coefficients2C,'o',xfit,interpolated);
xlim([0,100]);
ylabel('Entropic Coefficient (mV/K)');
xlabel('SOC (%)');
title('Entropic Coefficient vs SOC, 2C discharge pulses');

coeffdata2C = [];
for i = 1:71
    coeffdata2C(i+29,1) = round(xfit(1,i),0);
    coeffdata2C(i+29,2) = interpolated(1,i);
end

fulldata = xlsread("FE11 Endurance Full Data V2.xlsx");
soc_temporary = fulldata(:,7);
current_temporary = fulldata(:,3);
temp_temporary = fulldata(:,4);

for i = 1:2:size(current_temporary)
    if i < 9238
        new_current(i,1) = current_temporary(i) + current_temporary(i+1);
    end
end

for i = 1:9237
    value = soc_temporary(i);
    value = int64(value);
    if value == 0
        continue;
    end
    if value > 100
        continue;
    end
    if value == coeffdata2C(value,1) 
        revheat_cell2C(i,1) = new_current(i,1)/4 * (fulldata(i,4)+273.15) * (coeffdata2C(value,2)/1000); %current*temp*dOCV/dT, 
        % current divided by 4 to get the reversible heat of one cell
    end
end

revheat_time2C = [];
for i = 1:size(revheat_cell2C)
    revheat_time2C(i,1) = fulldata(i,1);
    revheat_time2C(i,2) = revheat_cell2C(i,1);
end


%4C
coefficients4C = [-0.11108398,-0.13713834,-0.11940002,-0.5365753,0.16471862,0.10856632,0.11627198,0.0710678,0.10696412,-0.00057222];
% mV/K

figure
x = linspace(10,100,10);
xfit = linspace(10,101,92);
pp = spline(x,coefficients4C);
interpolated = fnval(pp,xfit);
plot(x,coefficients4C,'o',xfit,interpolated);
xlim([0,100]);
ylabel('Entropic Coefficient (mV/K)');
xlabel('SOC (%)');
title('Entropic Coefficient vs SOC, 4C discharge pulses');

coeffdata4C = [];
for i = 1:91
    coeffdata4C(i+9,1) = round(xfit(1,i),0);
    coeffdata4C(i+9,2) = interpolated(1,i);
end

for i = 1:9237
    value = soc_temporary(i);
    value = int64(value);
    if value == 0
        continue;
    end
    if value > 100
        continue;
    end
    if value == coeffdata4C(value,1) 
        revheat_cell4C(i,1) = new_current(i,1)/4 * (fulldata(i,4)+273.15) * (coeffdata4C(value,2)/1000); %current*temp*dOCV/dT, 
        % current divided by 4 to get the reversible heat of one cell
    end
end

revheat_time4C = [];
for i = 1:size(revheat_cell4C)
    revheat_time4C(i,1) = fulldata(i,1);
    revheat_time4C(i,2) = revheat_cell4C(i,1);
end
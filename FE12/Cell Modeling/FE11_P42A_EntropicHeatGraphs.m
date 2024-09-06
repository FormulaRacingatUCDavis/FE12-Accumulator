%entropic heat coefficients
% 313K = 40 C, 303K = 30 C, 298K = 25 C
%see spreadsheet "FE12_P42A_EntropicCoefficient_notavgOCV for full data for code below 

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
%Updated upstream
        revheat_pack(i,1) = new_current(i,1) * (fulldata(i,4)+273.15) * (coeffdata(value,2)/1000);
        revheat_cell(i,1) = new_current(i,1)/4 * (fulldata(i,4)+273.15) * (coeffdata(value,2)/1000); %current*temp*dOCV/dT, 
        % current divided by 4 to get the reversible heat of one cell

        revheat_pack(i,1) = new_current(i,1) * fulldata(i,4) * (coeffdata(value,2)/1000);
        revheat_cell(i,1) = new_current(i,1)/4 * (new_temp(i,1)+273.15) * (coeffdata(value,2)/1000); %current*temp*dOCV/dT, 
        % current divided by 4 to get the reversible heat of one cell

 %Stashed changes
    end
end

revheat_time = [];
for i = 1:length(revheat_cell)
    revheat_time(i,1) = fulldata(i,1);
    revheat_time(i,2) = revheat_cell(i,1);
end





% radial_thermalconductivity = 2.21 % W/m*K
% axial_thermalconductivity = 20.61 % W/m*K
% density = 2682 % kg/m^3
% specificheatcapacity = 1360 % J/kg*K
% emissivity = 0.24
% 
% alum_3003_thermalconductivity = 193 % W/m*K

%entropic heat coefficient
% 313K = 40 C, 303K = 30 C, 298K = 25 C

% coefficients using method of average OCV but OCV changes with time after
% pulse, data looks wrong

% time = [];
% voltage313 = [];
% voltage303 = [];
% voltage298 = [];
% 
% pulsedata313 = xlsread("FE11_HPPC_313_DataOnly.xlsx");
% for i = 1:length(pulsedata313)
%     voltage313(i,1) = pulsedata313(i,2);
% end
% 
% pulse100_313 = mean(voltage313(92:110));
% pulse90_313 = mean(voltage313(175:205));
% pulse80_313 = mean(voltage313(270:300));
% pulse70_313 = mean(voltage313(365:395));
% pulse60_313 = mean(voltage313(460:490));
% pulse50_313 = mean(voltage313(555:585));
% pulse40_313 = mean(voltage313(650:680));
% pulse30_313 = mean(voltage313(745:775));
% pulse20_313 = mean(voltage313(841:870));
% pulse10_313 = mean(voltage313(935:965));
% 
% pulsedata303 = xlsread("FE11_HPPC_303_DataOnly.xlsx");
% for i = 1:length(pulsedata303)
%     voltage303(i,1) = pulsedata303(i,2);
% end
% 
% pulse100_303 = mean(voltage303(87:114));
% pulse90_303 = mean(voltage303(184:209));
% pulse80_303 = mean(voltage303(277:304));
% pulse70_303 = mean(voltage303(369:399));
% pulse60_303 = mean(voltage303(464:494));
% pulse50_303 = mean(voltage303(559:589));
% pulse40_303 = mean(voltage303(654:684));
% pulse30_303 = mean(voltage303(749:779));
% pulse20_303 = mean(voltage303(844:874));
% pulse10_303 = mean(voltage303(944:969));
% 
% pulsedata298 = xlsread("FE11_HPPC_298_DataOnly.xlsx");
% for i = 1:length(pulsedata298)
%     voltage298(i,1) = pulsedata298(i,2);
% end
% 
% pulse100_298 = mean(voltage298(87:117));
% pulse90_298 = mean(voltage298(182:212));
% pulse80_298 = mean(voltage298(279:307));
% pulse70_298 = mean(voltage298(373:402));
% pulse60_298 = mean(voltage298(467:497));
% pulse50_298 = mean(voltage298(562:592));
% pulse40_298 = mean(voltage298(658:687));
% pulse30_298 = mean(voltage298(752:782));
% pulse20_298 = mean(voltage298(847:877));
% pulse10_298 = mean(voltage298(942:972));
% 
% temp = [298, 303, 313];
% pulse100 = [pulse100_298,pulse100_303,pulse100_313];
% pulse90 = [pulse90_298,pulse90_303,pulse90_313];
% pulse80 = [pulse80_298,pulse80_303,pulse80_313];
% pulse70 = [pulse70_298,pulse70_303,pulse70_313];
% pulse60 = [pulse60_298,pulse60_303,pulse60_313];
% pulse50 = [pulse50_298,pulse50_303,pulse50_313];
% pulse40 = [pulse40_298,pulse40_303,pulse40_313];
% pulse30 = [pulse30_298,pulse30_303,pulse30_313];
% pulse20 = [pulse20_298,pulse20_303,pulse20_313];
% pulse10 = [pulse10_298,pulse10_303,pulse10_313];
% 
% 
% figure
% subplot(2,5,1)
% P = polyfit(temp,pulse100,1);
% xfit = linspace(298,313,100);
% yfit100 = polyval(P,xfit);
% plot(temp,pulse100,'r');
% title('100% SOC')
% slope100 = P(1);
% 
% subplot(2,5,2)
% P = polyfit(temp,pulse90,1);
% yfit90 = polyval(P,xfit);
% plot(temp,pulse90,'r-');
% title('90% SOC')
% slope90 = P(1);
% 
% subplot(2,5,3)
% P = polyfit(temp,pulse80,1);
% yfit80 = polyval(P,xfit);
% plot(temp,pulse80,'r-');
% title('80% SOC')
% slope80 = P(1);
% 
% subplot(2,5,4)
% P = polyfit(temp,pulse70,1);
% yfit70 = polyval(P,xfit);
% plot(temp,pulse70,'r-');
% title('70% SOC')
% slope70 = P(1);
% 
% subplot(2,5,5)
% P = polyfit(temp,pulse60,1);
% yfit60 = polyval(P,xfit);
% plot(temp,pulse60,'r-');
% title('60% SOC')
% slope60 = P(1);
% 
% subplot(2,5,6)
% P = polyfit(temp,pulse50,1);
% yfit50 = polyval(P,xfit);
% plot(temp,pulse50,'r-');
% title('50% SOC')
% slope50 = P(1);
% 
% subplot(2,5,7)
% P = polyfit(temp,pulse40,1);
% yfit40 = polyval(P,xfit);
% plot(temp,pulse40,'r-');
% title('40% SOC')
% slope40 = P(1);
% 
% subplot(2,5,8)
% P = polyfit(temp,pulse30,1);
% yfit30 = polyval(P,xfit);
% plot(temp,pulse30,'r-');
% title('30% SOC')
% slope30 = P(1);
% 
% subplot(2,5,9)
% P = polyfit(temp,pulse20,1);
% yfit20 = polyval(P,xfit);
% plot(temp,pulse20,'r-');
% title('20% SOC')
% slope20 = P(1);
% 
% subplot(2,5,10)
% P = polyfit(temp,pulse10,1);
% yfit10 = polyval(P,xfit);
% plot(temp,pulse10,'r-');
% title('10% SOC')
% slope10 = P(1);
% 
% 
% figure
% subplot(2,5,1)
% plot(xfit,yfit100,'b-');
% title('100% SOC')
% subplot(2,5,2)
% plot(xfit,yfit90,'b-');
% title('90% SOC')
% subplot(2,5,3)
% plot(xfit,yfit80,'b-');
% title('80% SOC')
% subplot(2,5,4)
% plot(xfit,yfit70,'b-');
% title('70% SOC')
% subplot(2,5,5)
% plot(xfit,yfit60,'b-');
% title('60% SOC')
% subplot(2,5,6)
% plot(xfit,yfit50,'b-');
% title('50% SOC')
% subplot(2,5,7)
% plot(xfit,yfit40,'b-');
% title('40% SOC')
% subplot(2,5,8)
% plot(xfit,yfit30,'b-');
% title('30% SOC')
% subplot(2,5,9)
% plot(xfit,yfit20,'b-');
% title('20% SOC')
% subplot(2,5,10)
% plot(xfit,yfit10,'b-');
% title('10% SOC')
% 
% % format long
% % P100 = polyfit(xfit,yfit100,1);
% % P90 = polyfit(xfit,yfit90,1);
% % P80 = polyfit(xfit,yfit80,1);
% % P70 = polyfit(xfit,yfit70,1);
% % P60 = polyfit(xfit,yfit60,1);
% % P50 = polyfit(xfit,yfit50,1);
% % P40 = polyfit(xfit,yfit40,1);
% % P30 = polyfit(xfit,yfit30,1);
% % P20 = polyfit(xfit,yfit20,1);
% % P10 = polyfit(xfit,yfit10,1);
% 
% cumulative = [slope10,slope20,slope30,slope40,slope50,slope60,slope70,slope80,slope90,slope100];
% x = linspace(10,100,10);
% figure
% plot(x,cumulative,'r-');
% 
% xfit = linspace(10,101,92);
% pp = spline(x,cumulative);
% interpolated = fnval(pp,xfit);
% plot(x,cumulative,'o',xfit,interpolated);
% xlim([0,100]);
% ylabel('Entropic Coefficient (V/K)');
% xlabel('SOC (%)');
% title('Entropic Coefficient vs SOC');



%see spreadsheet "FE12_P42A_EntropicCoefficient_notavgOCV for full data for
%code below 

coefficients = [1.025390625, 0.147601536, -0.096130371,0.265421186, 0.334903172,0.381578718, 0.219563075, 0.213922773, 0.14163426, 0.092015948];
% mV/K
% 
figure
x = linspace(10,100,10);
xfit = linspace(10,101,92);
pp = spline(x,coefficients);
interpolated = fnval(pp,xfit);
plot(x,coefficients,'o',xfit,interpolated);
xlim([0,100]);
ylabel('Entropic Coefficient (V/K)');
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
        new_temp(i,1) = (temp_temporary(i) + temp_temporary(i+1))/2;
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
        revheat_pack(i,1) = new_current(i,1) * fulldata(i,4) * (coeffdata(value,2)/1000);
        revheat_cell(i,1) = new_current(i,1)/4 * new_temp(i,1) * (coeffdata(value,2)/1000); %current*temp*dOCV/dT, 
        % current divided by 4 to get the reversible heat of one cell

        revheat_pack(i,1) = revheat_pack(i,1);
        revheat_cell(i,1) = revheat_cell(i,1);
    end
end


revheat_time = [];
for i = 1:length(revheat_cell)
    revheat_time(i,1) = fulldata(i,1);
    revheat_time(i,2) = revheat_cell(i,1);
end

% bms = out.BMS_discrepency;
% plot(bms);
% bms_error = sum(bms)/26200; %average error, but there is a big outlier,
% 2.3507, I estimate error to be around 4 degC



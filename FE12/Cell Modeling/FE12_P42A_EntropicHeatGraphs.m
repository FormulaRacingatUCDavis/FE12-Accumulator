% radial_thermalconductivity = 2.21 % W/m*K
% axial_thermalconductivity = 20.61 % W/m*K
% density = 2682 % kg/m^3
% specificheatcapacity = 1360 % J/kg*K
% emissivity = 0.24
% 
% alum_3003_thermalconductivity = 193 % W/m*K

%entropic heat coefficient
% 313K = 40 C, 303K = 30 C, 298K = 25 C

time = [];
voltage313 = [];
voltage303 = [];
voltage298 = [];

pulsedata313 = xlsread("FE11_HPPC_313_DataOnly.xlsx");
for i = 1:length(pulsedata313)
    voltage313(i,1) = pulsedata313(i,2);
end

pulse100_313 = mean(voltage313(59:110));
pulse90_313 = mean(voltage313(175:237));
pulse80_313 = mean(voltage313(270:300));
pulse70_313 = mean(voltage313(365:395));
pulse60_313 = mean(voltage313(460:490));
pulse50_313 = mean(voltage313(555:585));
pulse40_313 = mean(voltage313(650:680));
pulse30_313 = mean(voltage313(745:775));
pulse20_313 = mean(voltage313(841:870));
pulse10_313 = mean(voltage313(935:965));

pulsedata303 = xlsread("FE11_HPPC_303_DataOnly.xlsx");
for i = 1:length(pulsedata303)
    voltage303(i,1) = pulsedata303(i,2);
end

pulse100_303 = mean(voltage313(57:114));
pulse90_303 = mean(voltage313(179:241));
pulse80_303 = mean(voltage313(275:304));
pulse70_303 = mean(voltage313(369:399));
pulse60_303 = mean(voltage313(464:494));
pulse50_303 = mean(voltage313(559:589));
pulse40_303 = mean(voltage313(654:684));
pulse30_303 = mean(voltage313(749:779));
pulse20_303 = mean(voltage313(844:874));
pulse10_303 = mean(voltage313(944:969));

pulsedata298 = xlsread("FE11_HPPC_298_DataOnly.xlsx");
for i = 1:length(pulsedata298)
    voltage298(i,1) = pulsedata298(i,2);
end

pulse100_298 = mean(voltage298(56:117));
pulse90_298 = mean(voltage298(182:212));
pulse80_298 = mean(voltage298(279:307));
pulse70_298 = mean(voltage298(373:402));
pulse60_298 = mean(voltage298(467:555));
pulse50_298 = mean(voltage298(562:592));
pulse40_298 = mean(voltage298(658:698));
pulse30_298 = mean(voltage298(752:814));
pulse20_298 = mean(voltage298(847:877));
pulse10_298 = mean(voltage298(942:972));

temp = [298, 303, 313];
pulse100 = [pulse100_298,pulse100_303,pulse100_313];
pulse90 = [pulse90_298,pulse90_303,pulse90_313];
pulse80 = [pulse80_298,pulse80_303,pulse80_313];
pulse70 = [pulse70_298,pulse70_303,pulse70_313];
pulse60 = [pulse60_298,pulse60_303,pulse60_313];
pulse50 = [pulse50_298,pulse50_303,pulse50_313];
pulse40 = [pulse40_298,pulse40_303,pulse40_313];
pulse30 = [pulse30_298,pulse30_303,pulse30_313];
pulse20 = [pulse20_298,pulse20_303,pulse20_313];
pulse10 = [pulse10_298,pulse10_303,pulse10_313];


figure
subplot(2,5,1)
P = polyfit(temp,pulse100,1);
xfit = linspace(298,313,100);
yfit100 = polyval(P,xfit);
plot(temp,pulse100,'r');
title('100% SOC')
slope100 = P(1);

subplot(2,5,2)
P = polyfit(temp,pulse90,1);
yfit90 = polyval(P,xfit);
plot(temp,pulse90,'r-');
title('90% SOC')
slope90 = P(1);

subplot(2,5,3)
P = polyfit(temp,pulse80,1);
yfit80 = polyval(P,xfit);
plot(temp,pulse80,'r-');
title('80% SOC')
slope80 = P(1);

subplot(2,5,4)
P = polyfit(temp,pulse70,1);
yfit70 = polyval(P,xfit);
plot(temp,pulse70,'r-');
title('70% SOC')
slope70 = P(1);

subplot(2,5,5)
P = polyfit(temp,pulse60,1);
yfit60 = polyval(P,xfit);
plot(temp,pulse60,'r-');
title('60% SOC')
slope60 = P(1);

subplot(2,5,6)
P = polyfit(temp,pulse50,1);
yfit50 = polyval(P,xfit);
plot(temp,pulse50,'r-');
title('50% SOC')
slope50 = P(1);

subplot(2,5,7)
P = polyfit(temp,pulse40,1);
yfit40 = polyval(P,xfit);
plot(temp,pulse40,'r-');
title('40% SOC')
slope40 = P(1);

subplot(2,5,8)
P = polyfit(temp,pulse30,1);
yfit30 = polyval(P,xfit);
plot(temp,pulse30,'r-');
title('30% SOC')
slope30 = P(1);

subplot(2,5,9)
P = polyfit(temp,pulse20,1);
yfit20 = polyval(P,xfit);
plot(temp,pulse20,'r-');
title('20% SOC')
slope20 = P(1);

subplot(2,5,10)
P = polyfit(temp,pulse10,1);
yfit10 = polyval(P,xfit);
plot(temp,pulse10,'r-');
title('10% SOC')
slope10 = P(1);


figure
subplot(2,5,1)
plot(xfit,yfit100,'b-');
title('100% SOC')
subplot(2,5,2)
plot(xfit,yfit90,'b-');
title('90% SOC')
subplot(2,5,3)
plot(xfit,yfit80,'b-');
title('80% SOC')
subplot(2,5,4)
plot(xfit,yfit70,'b-');
title('70% SOC')
subplot(2,5,5)
plot(xfit,yfit60,'b-');
title('60% SOC')
subplot(2,5,6)
plot(xfit,yfit50,'b-');
title('50% SOC')
subplot(2,5,7)
plot(xfit,yfit40,'b-');
title('40% SOC')
subplot(2,5,8)
plot(xfit,yfit30,'b-');
title('30% SOC')
subplot(2,5,9)
plot(xfit,yfit20,'b-');
title('20% SOC')
subplot(2,5,10)
plot(xfit,yfit10,'b-');
title('10% SOC')

% format long
% P100 = polyfit(xfit,yfit100,1);
% P90 = polyfit(xfit,yfit90,1);
% P80 = polyfit(xfit,yfit80,1);
% P70 = polyfit(xfit,yfit70,1);
% P60 = polyfit(xfit,yfit60,1);
% P50 = polyfit(xfit,yfit50,1);
% P40 = polyfit(xfit,yfit40,1);
% P30 = polyfit(xfit,yfit30,1);
% P20 = polyfit(xfit,yfit20,1);
% P10 = polyfit(xfit,yfit10,1);

cumulative = [slope10,slope20,slope30,slope40,slope50,slope60,slope70,slope80,slope90,slope100];
x = linspace(10,100,10);
figure
plot(x,cumulative,'r-');

xfit = linspace(10,101,92);
pp = spline(x,cumulative);
interpolated = fnval(pp,xfit);
plot(x,cumulative,'o',xfit,interpolated);
xlim([0,100]);
ylabel('Entropic Coefficient (V/K)');
xlabel('SOC (%)');
title('Entropic Coefficient vs SOC');

pulldata = [];
for i = 1:91
    pulldata(i+9,1) = round(xfit(1,i),0);
    pulldata(i+9,2) = interpolated(1,i);
end


%Testing Parameters:
current = 40; %Amps(A)
initialTemp = 30; %C
crossSectionalArea =  40 * 0.01; %cm^2 [mm^2 -> cm^2] %%THIS IS THE VALUE YOU WANT TO CHANGE TO TEST DIFFERENT CROSS SECTIONS
time = 1800; %seconds 
Length = 57.5e-2;%cm

%Material Specs: Aluminum 3003
density = 2.73; %g/cm^3
specificHeatCapacity = 0.893; %J/g*C
resistivity = 4.16E-06; %ohm-cm at 20 C
tempCoefficientResistivity = 0.0038;
initialResistivity = resistivity*(1 + tempCoefficientResistivity*(initialTemp-20)); %Resistivity at Initial Temperature
thermalconductivty = 190/100; % W/mK --> W/cmK
heattransfercoeff = thermalconductivty/Length; 
emissitivty = 0.3;

Stefanboltzman_constant = 5.67e-8 * 1e4; % W / m^2*K^4 --> W / cm^2*K^4


%Variables Initialized
busbarTemp = initialTemp;
busbarTempData = [];
tempChange = 0;
totalTime = 0:time;

function F = tempChangefunc(x)

F(1) = (current^2*finalresistivity) / (1/((-density*specificHeatCapacity*x*(heattransfercoeff*x+emissitivty*Stefanboltzman_constant*x^4)/time^2)-1/Length)) + crossSectionalArea^2;

end

for i = 1:31 % loop through possible cross sectional areas
    crossSectionalArea = (i + 40)*0.01;
    tempChange = 0;
    busbarTemp = initialTemp;
    
    for t = 1:time 
        tempChange = (current^2*initialResistivity)/(density*crossSectionalArea^2*specificHeatCapacity);
        busbarTemp = busbarTemp + tempChange;
        initialResistivity = resistivity*(1 + tempCoefficientResistivity*(busbarTemp-20));
        busbarTempData(t+1) = busbarTemp;

        busbarTempData = transpose(busbarTempData); 
        finalresistivity = initialResistivity;
        busbarMaxTemp = busbarTempData(end-1);

        
    
        % busbarTemp = busbarTemp + tempChange;
        % initialResistivity = resistivity*(1 + tempCoefficientResistivity*(busbarTemp-20));
        % 
        % busbarTempData_wrt_area(i,t) = busbarTemp;
        % busbarMaxTemp_wrt_area(i) = max(busbarTemp);
    
    end

    fun = @tempChangefunc;
    x0 = busbarMaxTemp - initialTemp;
    x = fsolve(fun,x0);

end

% figure
% for i = 1:31
%     plot(1:length(busbarTempData_wrt_area),busbarTempData_wrt_area(i,:))
%     hold on
% end
% 
% figure
% plot(1:length(busbarTempData_wrt_area),busbarTempData_wrt_area(1,:))
% 
% xlim([1,1800])
% ylim([1,175])
%busbarMaxTemp is the maximum temperature the busbar will reach according to
%your parameters

% function F = tempChangefunc(x)
% 
% F(1) = time*(current^2*initialResistivity*Length/x - heattransfercoeff*crossSectionalArea*x-emissitivty*Stefanboltzman_constant*crossSectionalArea*x^4) - Length*density*specificHeatCapacity*x;
% end
% fun = @tempChangefunc;
% x0 = 20;
% x = fsolve(fun,x0)
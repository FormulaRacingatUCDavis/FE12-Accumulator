
%Testing Parameters:
current = 40; %Amps(A)
initialTemp = 30; %C
crossSectionalArea =  40 * 0.01; %cm^2 [mm^2 -> cm^2] %%THIS IS THE VALUE YOU WANT TO CHANGE TO TEST DIFFERENT CROSS SECTIONS
time = 1800; %seconds 
voltage = ;

%Material Specs: Aluminum 3003
density = 2.73; %g/cm^3
specificHeatCapacity = 0.893; %J/g*C
resistivity = 2.65E-06; %ohm-cm at 20 C
tempCoefficientResistivity = 0.0038;
initialResistivity = resistivity*(1 + tempCoefficientResistivity*(initialTemp-20)); %Resistivity at Initial Temperature
thermalconductivty = 190; % W/mK
heattransfercoeff = thermalconductivty/Length; % length?
emissitivty = 0.3;

Stefanboltzman_constant = 5.67e-8; % W / m^2*K^4

busbarTemp = initialTemp;
busbarTempData = [];
tempChange = 0;
for i = 40:71 % loop through possible cross sectional areas
    for t = 0:time 
    
        %tempChange = ( (current^2*initialResistivity*Length*time) / (Length*density*specificHeatCapacity + heattransfercoeff + ...
         %   emissitivty*Stefanboltzman_constant*busbarTemp^3) ) / crossSectionalArea(i)^2;
    
        tempChange = ( (current^2*initialResistivity)/(crossSectionalArea(i)^2*density*specificHeatCapacity) + (current*voltage)/(crossSectionalArea(i)*(heattransfercoeff+ ...
            emissitivty*Stefanboltzman_constant*busbarTemp^3)) ) * time;
    
        busbarTemp = busbarTemp + tempChange;
        initialResistivity = resistivity*(1 + tempCoefficientResistivity*(busbarTemp-20));
        busbarTempData(t+1) = busbarTemp;

        busbarTempData_wrt_area(i,t+1) = busbarTemp;
        busbarMaxTemp_wrt_area(i) = busbarTempData(end-1);
    
    end
end

busbarTempData = transpose(busbarTempData); 
busbarMaxTemp = busbarTempData(end-1);

%busbarMaxTemp is the maximum temperature the busbar will reach according to
%your parameters
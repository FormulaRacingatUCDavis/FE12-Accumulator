
%Testing Parameters:
current = 65; %Amps(A)
initialTemp = 30; %C
crossSectionalArea = 52 * 0.01; %cm^2 [mm^2 -> cm^2] %%THIS IS THE VALUE YOU WANT TO CHANGE TO TEST DIFFERENT CROSS SECTIONS
time = 1800; %seconds 

%Material Specs: Copper
density = 8.93; %g/cm^3
specificHeatCapacity = 0.385; %J/g*C
resistivity = 1.70E-06; %ohm-cm at 20 C
tempCoefficientResistivity = 0.0039;
initialResistivity = resistivity*(1 + tempCoefficientResistivity*(initialTemp-20)); %Resistivity at Initial Temperature

%Variables Initialized
busbarTemp = initialTemp;
busbarTempData = [];
tempChange = 0;
totalTime = 0:time;

for t = 0:time 

    tempChange = (current^2*initialResistivity)/(density*crossSectionalArea^2*specificHeatCapacity);
    busbarTemp = busbarTemp + tempChange;
    initialResistivity = resistivity*(1 + tempCoefficientResistivity*(busbarTemp-20));
    busbarTempData(t+1) = busbarTemp;

end

busbarTempData = transpose(busbarTempData); 
busbarMaxTemp = busbarTempData(end-1);

plot(totalTime,busbarTempData);

%busbarMaxTemp is the maximum temperature the busbar will reach according to
%your parameters

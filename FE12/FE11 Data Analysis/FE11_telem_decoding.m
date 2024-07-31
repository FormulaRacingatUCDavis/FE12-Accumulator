TelemData  = xlsread("EnduranceBlueMax_parsed.csv");
%[num,text,raw] = xlsread("EnduranceBlueMax_parsed.csv"); 

VoltageData = [ ]; %packvoltage
CellTemperatureData = [ ]; %HI temp
CurrentData = [ ];
SOC = [ ]; %percentage
Torque = [ ];
Speed = [ ]; %angular velocity in RPM
timeCount = 1;
timeCount2 = 1;
timeCount3 = 1;
timeCount5 = 1;

%keyword = 'c0';
%first_column = text(:,1);
%matches = strcmp(first_column,keyword);

FullData = [ ];
AutoData = [ ];

for i = 1:length(TelemData)

    if(TelemData(i,1) == 380)

        VoltageData(timeCount,1) = round(TelemData(i,10),1); %Voltage Reading Time
        VoltageData(timeCount,2) = TelemData(i,6); %Voltage At Such Time

        CellTemperatureData(timeCount,1) = round(TelemData(i,10),1); %Cell Temp Reading Time
        CellTemperatureData(timeCount,2) = TelemData(i,2); %Cell Temp At Such Time
        
        SOC(timeCount,1) = round(TelemData(i,10),1); %SOC Reading Time
        SOC(timeCount,2) = TelemData(i,3); %SOC At Such Time

        timeCount = timeCount + 1;

    end

    if(TelemData(i,1) == 387)

        CurrentData(timeCount2,1)= round(TelemData(i,10),1); %Current Reading Tim
        CurrentData(timeCount2,2) = TelemData(i,2); %Current At Such Time

        timeCount2 = timeCount2 + 1;

    end


    %if any(matches(:))
    if(TelemData(i,1) == 0)

         Torque(timeCount3,1) = round(TelemData(i,10),1); %Torque Reading Time
         Torque(timeCount3,2) = TelemData(i,2); %Torque At Such Time

         timeCount3 = timeCount3 + 1;
    end

    if(TelemData(i,1) == 5)

         Speed(timeCount5,1) = round(TelemData(i,10),1); %Speed Reading Time
         Speed(timeCount5,2) = TelemData(i,4); %Speed At Such Time

         timeCount5 = timeCount5 + 1;

    end %15.5in or 16in

end

timeCount4 = 1;

for i = 1:length(VoltageData)
    for j = 1:length(CurrentData)

        if VoltageData(i,1) == CurrentData(j,1)

            FullData(timeCount4,1) = VoltageData(i,1);
            FullData(timeCount4,2) = VoltageData(i,2);
            FullData(timeCount4,3) = CurrentData(j,2);

            %FullData(timeCount4,4) = CellTemperatureData(i,2);
            %FullData(timeCount4,5) = SOC(i,2);
            %FullData(timeCount4,6) = Torque(i,2);
            %FullData(timeCount4,7) = Speed(i,2);

            timeCount4 = timeCount4 + 1;

        end

    end

end

for i = 1:length(FullData)
    for j = 1:length(CellTemperatureData)

        if FullData(i,1) == CellTemperatureData(j,1)

            FullData(i,4) = CellTemperatureData(j,2);

        end
    end
end

for i = 1:length(FullData)
    for j = 1:length(Torque)

        if FullData(i,1) == Torque(j,1)

            FullData(i,5) = Torque(j,2);

        end
    end
end

for i = 1:length(FullData)
    for j = 1:length(Speed)

        if FullData(i,1) == Speed(j,1)

            FullData(i,6) = Speed(j,2);

        end
    end
end

for i = 1:length(FullData)
    for j = 1:length(SOC)

        if FullData(i,1) == SOC(j,1)

            FullData(i,7) = SOC(j,2);

        end
    end
end

%writematrix(FullData,'Full Data_Run16.xls');

% %get Data for Autocross runs without time matching
% timeCount5 = 1
% 
% for i = 1:64           
%     AutoData(timeCount5,1) = VoltageData(i,1);
%     AutoData(timeCount5,2) = VoltageData(i,2);
%     AutoData(timeCount5,3) = CurrentData(i,2);
%     AutoData(timeCount5,4) = CellTemperatureData(i,2);
%     AutoData(timeCount5,5) = SOC(i,2);
%     AutoData(timeCount5,6) = Torque(i,2);
%     AutoData(timeCount5,7) = Speed(i,2);
% 
%     timeCount5 = timeCount5 + 1;
% 
% end
% writematrix(AutoData,'Auto Data_Run16.xls');

revheat = zeros(length(SOC), 1);
soc_temp = SOC(:,2);

for i = 1:length(soc_temp)
    value = soc_temp(i);
    value = int64(value);
    if value == 0
        continue;
    end
    if value > 100
        continue;
    end
    %if mod(i,2) == 1 %since there are two readings for each temperature
    if soc_temp(i) == pulldata(value,1) 
        revheat(i) = CellTemperatureData(i,2) * pulldata(value,2); %multiply this by current in simulink
    end
    %else
        continue;
    %end
end

total = sum(revheat);
deltaT = total/(0.007*1360); %heat released by convection


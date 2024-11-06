data = xlsread("cooling curve P45B.xlsx")

time = data(:,1)-9868.61;
Tcell = data(:,2)+273.15;
deltaT = data(:,3)+273.15;

area = 5478.73/1e6;
Tamb = 30 + 273.15;
cp = 1360;
mass = 0.07;

% x0 = [1];
% h_opt = zeros(size(Tcell)); 
% for i = 2:length(Tcell)
%     fun = @(h) mass*cp*deltaT(i)/(time(i)-time(i-1)) - h*area*(Tamb-Tcell(i));
%     [h_opt(i),fval(i)] = fminsearch(fun, x0);
% end
% 
% % Print results
% fprintf('Optimal h values:\n');
% disp(h_opt);

x0 = [1];
lb = -100;
ub = 0;
h_opt = zeros(size(Tcell)); 
for i = 2:length(Tcell)
    fun = @(h) mass*cp*deltaT(i)/(time(i)-time(i-1)) - h*area*(Tamb-Tcell(i));
    [h_opt(i),fval(i)] = fmincon(fun, x0,lb,ub);
end

% Print results
fprintf('Optimal h values:\n');
disp(h_opt);
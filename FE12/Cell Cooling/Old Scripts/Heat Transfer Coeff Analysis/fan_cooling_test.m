N_rows = 12;
rho = 1.1640;
X=1.1; % staggered graph
f = 0.9; % staggered graph
V_max = 0:.36/9:.36;
P_delta = N_rows.*X.*(rho.*V_max.^2./2).*f; %[Pa]
figure;
plot(V_max,P_delta)
xlabel("V_max")
ylabel("P [Pa]")
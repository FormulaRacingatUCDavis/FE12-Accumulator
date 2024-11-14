%%Heat Sink OptimizationData
clear all
close all
clc

% Constants
H = 0.14; %[m]
w = .32; %[m]
k = 237; %[W/m K]
h = 30; %[W/m^2 K]

S_min = 0.01; S_max = 0.02;  % [m]
t_f_min = 0.01; t_f_max = 0.02;% [m]
L_f_min = 0.001; L_f_max = 0.0158; % [m]

%N = 100;
N = 20;
S = S_min:(S_max-S_min)/(N-1):S_max;
t_f = t_f_min:(t_f_max-t_f_min)/(N-1):t_f_max;
L_f = L_f_min:(L_f_max-L_f_min)/(N-1):L_f_max;


Data = zeros(N,N,N);
N_f = zeros(size(Data));

for i = 1:length(L_f) %row
    for j = 1:length(t_f) % column
        for k = 1:length(S) %matrix number (ex. first matrix out of 20 total)
            [Data(k,j,i),N_f(i,j,k)] = effectiveness(S(k),t_f(j),L_f(i),w,H,h,k);
            % if Data(i,j,k) == inf
            %     Data(i,j,k) = 0;
            % end
        end
    end
end

% Data = zeros(N,N);
% 
% for i = 1:length(L_f)
%     for j = 1:length(t_f)
%         for k = 1:1
%             Data(i,j,k) = effectiveness(S(k),t_f(j),L_f(i),w,H,h,k);
%         end
%     end
% end

max = max(max(max(Data)));

x = L_f(1:19);
x = t_f(1:19);
y = S(1:19);

%[X,Y,Z] = meshgrid(x,y,z);
[X,Y] = meshgrid(x,y);
% surf(X,Y,Z,Data)
% figure()
% scatter3(X,Y,Z,50,Data,'filled');

for i = 1:19
    figure;
    surf(X*100,Y*100,Data(1:19,1:19,i), 'EdgeColor', 'none');
    %max_effectiveness(i) = max(Data(19,19,i));
    %min_effectiveness(i) = min(Data(1,1,1));
    colormap(spring); % Adjust colormap as needed
    colorbar;
    xlabel('Fin Length (cm)');
    ylabel('Fin Thickness (cm)');
    zlabel('Effectiveness');
    title('Effectiveness, Spacing = 2 cm')
view(3); % Adjust viewpoint as desired
end


%%
function [eff,N_f] = effectiveness(S,t_f,L_f, w, H, h, k)

    A_f = 2 * w * (L_f + t_f / 2);
    m = sqrt(h * (2 * w + 2 * t_f) / (k * w * t_f));
    N_f = (H + S) / (t_f + S);
    % fprintf('A_f: %.5f, m: %.5f, N_f: %.5f\n', A_f, m, N_f);
    eff =  A_f * (tanh(m*(L_f+t_f/2)) / (m*(L_f+t_f/2))) / (w*H - N_f*t_f*w);
end
% Define known constants
w = 0.32;  % Example width in meters (adjust to your problem)
H = 0.14;   % Example height in meters
h = 20;   % Example convective heat transfer coefficient (W/m^2-K)
k = 25;   % Example thermal conductivity of fin material (W/m-K)

% Set bounds for the variables
S_min = 0.01; S_max = H;  % Bounds for S (spacing between fins) in meters
t_f_min = 0.01; t_f_max = H; % Bounds for t_f (fin thickness) in meters
L_f_min = 0.001; L_f_max = 0.0158;  % Bounds for L_f (fin length) in meters

% Set initial guess for the variables
S_initial = 0.005;  % Initial guess for spacing between fins
t_f_initial = 0.002;  % Initial guess for fin thickness
L_f_initial = 0.001;  % Initial guess for fin length

% Combine initial guesses into a vector
initial_guess = [S_initial, t_f_initial, L_f_initial];

% Define bounds as matrices
lb = [S_min, t_f_min, L_f_min];  % Lower bounds for S, t_f, L_f
ub = [S_max, t_f_max, L_f_max];  % Upper bounds for S, t_f, L_f

% Define the anonymous objective function to pass extra arguments
objective = @(vars) effectiveness(vars, w, H, h, k);

% Call fmincon for constrained optimization
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[optimal_vars, fval] = fmincon(objective, initial_guess, [], [], [], [], lb, ub, [], options);

% Display the results
optimal_S = optimal_vars(1);
optimal_t_f = optimal_vars(2);
optimal_L_f = optimal_vars(3);

fprintf('Optimal Spacing (S): %.5f m\n', optimal_S);
fprintf('Optimal Fin Thickness (t_f): %.5f m\n', optimal_t_f);
fprintf('Optimal Fin Length (L_f): %.5f m\n', optimal_L_f);
fprintf('Maximum Effectiveness: %.5f\n', -fval);  % Effectiveness is the negative of fval

%%
function eff = effectiveness(vars, w, H, h, k)
    % Extract variables from the input
    S = vars(1);   % Spacing between fins
    t_f = vars(2); % Fin thickness
    L_f = vars(3); % Fin length

    % Calculate intermediate variables
    A_f = 2 * w * (L_f + t_f / 2);
    m = sqrt(h * (2 * w + 2 * t_f) / (k * w * t_f));
    N_f = (H + S) / (t_f + S);
    fprintf('A_f: %.5f, m: %.5f, N_f: %.5f\n', A_f, m, N_f);
    
    % Calculate effectiveness
    numerator = A_f * (tanh(m * L_f) / (m * L_f));
    denominator = w * H - N_f * t_f * w;
    
    % Output negative effectiveness for minimization
    eff = - numerator / denominator;
end

% Define symbolic variable
syms x

% Define the equation
eq = ((sqrt(3) - 1i)/20) + ((1i * sqrt(3))/10) * x == 0;

% Solve for x
sol = solve(eq, x);

% Display solution
disp('Solution for x:')
disp(sol)

%% Part 3
% Define symbolic variable
syms m2

% Define the coefficients
a1 = (sqrt(2)-3 + (-sqrt(2)-4)*1i) / 100;
a2 = (3/2 - sqrt(3)/2 + ((-sqrt(3))/2 - 11/2)*1i) / 100;

% Assume m1 = 0.5 and set up the equation
eq = a1 * 0.5 + a2 * m2 == 0;

% Solve for m2
sol_m2 = solve(eq, m2);

% Convert to numeric and truncate to 5 decimal places
sol_m2_truncated = vpa(sol_m2, 5);

% Display the solution
disp('Solution for m2 when m1 = 0.5 (truncated to 5 decimal places):')
disp(sol_m2_truncated)

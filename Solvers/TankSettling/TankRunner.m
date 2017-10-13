% Solve alpha equation for settling tanks problems
% Alpha is the continuous phase

% Physical and numerical paramaters
% see case file

% Load experimental data
T55 = load('T55First.dat');
T65 = load('T65First.dat');

% Discard null values in experimental data
thr = 1E-6;
i55 = max(find(T55(:, 2) < thr));
T55 = T55(i55:end, 1:2);

i65 = max(find(T65(:, 2) < thr));
T65 = T65(i65:end, 1:2);

% Shift times
T55(:, 1) = T55(:, 1) - T55(1, 1);
T65(:, 1) = T65(:, 1) - T65(1, 1);

% Plot data
if 1
    figure(3)
    plot(T55(:, 1), T55(:, 2), 'b-')
    hold on
    plot(T65(:, 1), T65(:, 2), 'r-')
end
% Run solver
alphaEqnIsolatedUm

% Plot solution
figure(1)
plot(xC,u.internal,'r*-')

% Plot front advancement
figure(3)
plot(frontData(:, 1), frontData(:, 2))


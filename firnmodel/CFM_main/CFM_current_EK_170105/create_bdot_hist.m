%% Create accumulation history files for CFM

% Create a .csv file that has bdot history for input into CFM

% Length of temperature history
L = 1000;

% Create year vector
year = 0:10:L;

% Create bdot vector
bdot = .08*ones(1,length(year));

% Combine into matrix
bdot_hist = [year; bdot];

% Save as a .csv file
dlmwrite('bdot08kyr1stp10.csv', bdot_hist, 'delimiter', ',');


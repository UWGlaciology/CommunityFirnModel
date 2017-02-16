%% Create temperature history files for CFM

% Create a .csv file that has temperature history for input into CFM

% Length of temperature history in years
L = 1000;

% Create year vector
year = 0:10:L;

% Create temperature vector
temp = 223.15*ones(1,length(year));

% Combine into matrix
temp_hist = [year; temp];

% Save as a .csv file
dlmwrite('T50kyr1stp10.csv', temp_hist, 'delimiter', ',');


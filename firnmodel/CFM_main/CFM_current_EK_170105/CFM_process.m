function [] = CFM_process(experimentname)

%% CFM results processing
% Process and plot results from a given results folder

% Import temperature and accumulation histories used
% in = strfind(experimentname, '_');
% Ts = importdata([experimentname(1:in(1)-1) '.csv']);
% bdot = importdata([experimentname(in(1)+1:in(2)-1) '.csv']);
% Ts = importdata('FM_temp2_4000yr.csv');
% bdot = importdata('FM_bdot13_4000yr.csv');

% CD into results folder for the chosen experiment. Save main directory as
% oldfolder
oldfolder = cd(experimentname);

% Import results
iso_spin1 = importdata('isoSpin.csv');
iso_mat = importdata('iso.csv');
% temp_spin = importdata('tempSpin.csv');
% temp_mat = importdata('temp.csv');
depth_mat = importdata('depth.csv');
% density = importdata('density.csv');

% Crop data
%Ts = Ts(2,:);
%bdot = bdot(2,:);
iso_spin = iso_spin1(2:end);
iso_mat = iso_mat(:,2:end);
% temp_spin = temp_spin(2:end);
% temp_mat = temp_mat(:,2:end);
depth_mat = depth_mat(:,2:end);
% density = density(1,2:end);

% Find total time of run T
T = size(iso_mat,1);    % total time run for to create time steps

% Plot only final isotope profile
h1 = figure;
plot(depth_mat(end,:), iso_mat(end,:))
xlabel('Depth (m)')
ylabel('Isotope (permil)')
title('Final isotope profile')
% axis([0 150 -53 -46])

% Isotope profile through time
h2 = figure;
plot(depth_mat(T/T,:), iso_mat(T/T,:))
xlabel('Depth (m)')
ylabel('isotope')
title('Plot isotope profiles through time versus depth vector')
hold on; % next plot steps through time
plot(depth_mat(floor(T/4),:), iso_mat(floor(T/4),:))
plot(depth_mat(floor(T/2),:), iso_mat(floor(T/2),:))
plot(depth_mat(floor(3*T/4),:), iso_mat(floor(3*T/4),:))
plot(depth_mat(T,:), iso_mat(T,:))
% axis([0 150 -53 -46])
legend('1 timestep', '1/4 timesteps', '1/2 timesteps', '3/4 timesteps', 'all timesteps')

% % DENSITY versus depth
% figure;
% plot(depth_mat(T/T,:), density)
% xlabel('Depth (m)')
% ylabel('Density')
% title('Plot density through time versus depth vector')

% % Plot input histories of temperature and accumulation
% h3 = figure;
% plot(Ts)
% xlabel('Time')
% ylabel('Temp')
% title('Surface temperature history')
% 
% h4 = figure;
% plot(bdot)
% xlabel('Time')
% ylabel('Accumulation rate')
% title('Surface accumulation rate history')

% Save results for this experiment
cd(oldfolder)      % go back to main folder
cd('Matlab_CFM_results')  % go into results folder
mkdir(experimentname)   % make new directory for this experiment
cd(experimentname)      % go into new directory

% Choose variables to save for later analysis
%save resultvariables.mat bdot Ts depth_mat iso_spin iso_mat 
save resultvariables.mat depth_mat iso_spin iso_mat 


% Save figures
savefig(h1, 'finaliso.fig')
savefig(h2, 'stepiso.fig')
%savefig(h3, 'temphist.fig')
%savefig(h4, 'bdothist.fig')


cd(oldfolder)      % go back to main folder
end
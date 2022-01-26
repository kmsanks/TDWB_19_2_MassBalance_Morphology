clear all 
close all
%The purpose of this script is to compare the Google Earth Engine
%delta hypsometry for four field-scale deltas: %1. Mississippi River Delta (MRD), Ganges Brahmaputra Meghna River Delta
%(GBMD), 3. Mekong River Delta (mekong), and 4. Rio Grande River Delta (riogrande)
%downloaded from hand-drawn polygons in Google Earth Engine (GEE) and
%to the elevation hypsometry of the control and treatment
%experiments
%This script produces Figure 5

%%First lets load the GEE delta data
cd('../data');
load('MRD.mat');
load('GBMD.mat');
load('mekong.mat');
load('riogrande.mat');

%%Now we need to load and process elevation data from the experiments
%In order to compare experimental data to field data, we are going to
%scale by channel depth. We will crop the experimental data to include
%only data between -1 and 3 channel depths reative to sea level (rsl),
%which is equivalent to -15 to 45 mm rsl

%%%Control
%%Load control data
load('TDB_18_data.mat')
cd '..\code'

%We only need the elevation data (ZD)
ZD_18 = Z_maps(:,:,1:560); %remove last 5 hours as that was post-experiment data collection, elevation data is in mm 
clear('C_maps','B_maps','G_maps','R_maps','Z_maps');
ZD_18(ZD_18(:,:,:)==0.)=NaN; %make data outside the basin NaN

%%Parameters needed for analysis
nx_18 = size(ZD_18, 1); %number of x locations on map
ny_18 = size(ZD_18,2); %number of y locations on map
nt_18 = size(ZD_18,3); %number of time steps in data set
dx = 5; %5 mm grid cells in x
dy = 5; %5 mm grid cells in y
dt_18 = 1; %delta t of time steps (hr)
baselevel_rr = 0.25; %base level rise rate (mm/hr)
ocean_zero = 25; %ocean elevation at beginning of experiment (mm)

%Loop through time steps to determine elevation reltive to sea level
elev_bin = -15:2:45; %-1 to 3 channel depths relative to sea level using 2 mm bin
k = 0; %start counter
elevations_18_1 = []; %initialize matrix of elevations
N_18 = nan(281,31); %initialize matrix of number of elevations, 31 is the number of bins
edges_18 = nan(281,31); %initialize matrix to store bin eges, 31 is the number of bins
for i=1:2:size(ZD_18,3); %we will use every other timestep for control to compare to treatment 
    k = k+1; %set counter
    %What is sea level at time i
    sl = (i*baselevel_rr*dt_18)+ocean_zero; 
    elevationmask = ZD_18(:,:,i); 
    elevationmask_rslr = elevationmask - sl;
    %clip data to remove data below and above -1 to 3 channel depths
    %relative to sea level
    elevationmask_rslr(elevationmask_rslr > 45) = NaN; %3 channel depths rsl
    elevationmask_rslr(elevationmask_rslr < -15) = NaN; %-1 channel depth rsl
    %remove nan to not bias distribution
    elevationmask_rslr_rmnan = elevationmask_rslr(~isnan(elevationmask_rslr));
    elevations_18_1(1:length(elevationmask_rslr_rmnan),k) = elevationmask_rslr_rmnan;
    %counts and bin edges
    [N, edges] = histcounts(elevationmask_rslr_rmnan(~isnan(elevationmask_rslr_rmnan)), 'normalization', 'probability','binwidth', 2);
    N_18(k, 1:length(N)) = N;
    edges_18(k, 1:length(edges)) = edges;
end 

%reshape elevation data
elev_18_1_vals = reshape(elevations_18_1, 1, []);

%%%Treatment
%%%Treatment experiment
%%Load treatment data
cd '..\data'
load('ZD_19_2_dry.mat'); %load topography array. In here should be a 3D topo array called ZD, oriented space x space x time
cd '..\code'
%change topography data to match the name of control experiment
ZD_19 = ZD_19_2_dry(:,:,1:280); %same time step as control
clear('ZD_19_2_dry');

%%Parameters needed for analysis
nx_19 = size(ZD_19,1); %number of x locations on map
ny_19 = size(ZD_19,2); %number of y locations on map
nt_19 = size(ZD_19,3); %number of time steps in data set
dt_19 = 2; %delta t of time steps (hr), dry scans only every 2 horus as compared to every hour for the control

%Loop through times to determine elevation distribution
elevations_19_2 = []; %initializ matrix for elevations
N_19 = nan(280,31); %initialize matrix for number, 31 is the number of bins 
edges_19 = nan(280,31); %initialzie matrix for edges, 31 is the number of bins
for i=1:nt_19;
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr*dt_19)+ocean_zero; %subtract 1 from i because scans start at hour 0 for te treatment
    elevationmask = ZD_19(:,:,i); 
    elevationmask_rslr = elevationmask - sl;
    %remove data below and above -1 and 3 channel depths rsl
    elevationmask_rslr(elevationmask_rslr > 45) = NaN;
    elevationmask_rslr(elevationmask_rslr < -15) = NaN;
    %remove nan to not bias distribution
    elevationmask_rslr_rmnan = elevationmask_rslr(~isnan(elevationmask_rslr));
    elevations_19_2(1:length(elevationmask_rslr_rmnan),i) = elevationmask_rslr_rmnan;
    [N, edges] = histcounts(elevationmask_rslr_rmnan(~isnan(elevationmask_rslr_rmnan)), 'normalization', 'probability','binwidth', 2);
    N_19(i, 1:length(N)) = N;
    edges_19(i, 1:length(edges)) = edges;
end 
%reshape the data
elev_19_2_vals = reshape(elevations_19_2, 1, []);

%%Now lets analyze and plot the dat
%Statistics
%mean and standard deviation of each bin
%control
mean_prob_18 = [];
std_prob_18 = [];
for i = 1:size(N_18,2);
    mean_prob_18(i) = nanmean(N_18(:,i)); 
    std_prob_18(i) = nanstd(N_18(:,i));
end 
%treatment
mean_prob_19 = [];
std_prob_19 = [];
for i = 1:size(N_19,2);
    mean_prob_19(i) = nanmedian(N_19(:,i)); 
    std_prob_19(i) = nanstd(N_19(:,i));
end 

%%Plot elevation histograms
%We need to normalize the data by channel depth to compare to field-scale
%deltas
%One channel depth in the control and treatment experiments is ~15 mm 
%Create array of mean bin, mean, and standard devition
%Control
array18 = [elev_bin; mean_prob_18; std_prob_18];
cols = any(isnan(array18),1);
array18(:,cols) = [];
%Treatment
array19 = [elev_bin; mean_prob_19; std_prob_19];
cols = any(isnan(array19),1);
array19(:,cols) = [];

%fill standard deviation
%control
y18 = array18(2,:); % your mean vector;
x18 = array18(1,:)/15; %normalize data by channel depth (elevation/channel depth (-))
std18 = array18(3,:);
curve1_18 = y18 + std18;
curve2_18 = y18 - std18;
%treatment
y19 = array19(2,:); % your mean vector;
x19 = array19(1,:)/15; %normalize data by channel depth (elevation/channel depth (-))
std19 = array19(3,:);
curve1_19 = y19 + std19;
curve2_19 = y19 - std19;

fig = figure()
plot(x18, y18, 'b:', 'LineWidth', 2)
hold on
plot(x19, y19, 'g:', 'LineWidth', 2)
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
alpha(0.25)
plot(MRD(:,3), MRD(:,4), 'linewidth', 2)
plot(GBMD(:,3), GBMD(:,4), 'linewidth', 2)
plot(mekong(:,3), mekong(:,4), 'linewidth', 2)
plot(riogrande(:,3), riogrande(:,4), 'linewidth', 2)
xlim([-1 3])
ylim([0 0.2])
xlabel('elevation/channel depth (-)')
ylabel('probability (-)')
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev', 'MRD', 'GBMD', 'Mekong', 'Rio Grande')
saveas(fig, '..\figures\Figure5_hypsometry.pdf')

%%Calculate the sum of probability between 0 and 0.5 for the treatment and
%%global
%GBMD
prob_GBMD = sum(GBMD(find(GBMD(:,3)>=0 & GBMD(:,3)<=0.5),4));
prob_mekong = sum(mekong(find(mekong(:,3)>=0 & mekong(:,3)<=0.5),4));
prob_MRD = sum(MRD(find(MRD(:,3)>=0 & MRD(:,3)<=0.5),4));
prob_riogrande = sum(riogrande(find(riogrande(:,3)>=0 & riogrande(:,3)<=0.5),4));
array19 = array19';
array19(:,1) = array19(:,1)/15;
prob_19 = sum(array19(find(array19(:,1)>=0 & array19(:,1)<=0.5),2));

array18 = array18';
array18(:,1) = array18(:,1)/15;
prob_18 = sum(array18(find(array18(:,1)>=0 & array18(:,1)<=0.5),2));

%%Find peak for control
[~,idx] = max(array18(:,2));%max percent to find peak
peak_18 = array18(idx,1:2); %first value is peak and second is percent

%%Test to make sure dominant shape of histgram shows in Figure 5 for each
%tmie step
figure()
for i=1:size(N_18,1)
    plot(N_18(i,:))
    hold on
end
hold off

figure()
for i=1:size(N_19,1)
    plot(N_19(i,:))
    hold on
end
hold off
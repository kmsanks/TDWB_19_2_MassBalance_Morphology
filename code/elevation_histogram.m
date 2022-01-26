clear all; close all
%This script calculates elevation distributions for both the control and
%treatment experiment. This data is used to create Figure 2b 

%%First, lets set our relative path, we should have opened this file from 
cd '..\' %We should have opened this code from 'C:\Users\XXXXXXX\Documents\GitHub\TDWB_19_2_MassBalance_Morphology\code'

%%%Control experiment
%%Load control data
cd '.\data'
load('ZD_18.mat') %topography; elevation data (mm)
cd '..\code'

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
elev_bin = -50:5:50; % the edges of the bins for this analysis
k = 0; %start counter
elevations_18_1 = nan(1000000,280); %initialize matrix of elevations
N_18 = nan(280,21); %initialize matrix of number of elevations, 51 is the number of bins
edges_18 = nan(280,21); %initialize matrix to store bin eges, 51 is the number of bins
for i=1:2:size(ZD_18,3); %we will use every other timestep for control to compare to treatment 
    k = k+1; %set counter
    %What is sea level at time i
    sl = (i*baselevel_rr*dt_18)+ocean_zero; 
    elevationmask = ZD_18(:,:,i); 
    elevationmask_rslr = elevationmask - sl;
    %clip data to remove extraneous data above and below about
    %50 mm
    elevationmask_rslr(elevationmask_rslr > 50) = NaN; %45
    elevationmask_rslr(elevationmask_rslr < -50) = NaN; %-50
    %remove nan to not bias distribution
    elevationmask_rslr(elevationmask_rslr == 0.) = NaN;
    elevationmask_rslr_rmnan = elevationmask_rslr(~isnan(elevationmask_rslr));
    elevations_18_1(1:length(elevationmask_rslr_rmnan),k) = elevationmask_rslr_rmnan;
    %counts and bin edges
    [N, edges] = histcounts(elevationmask_rslr_rmnan(~isnan(elevationmask_rslr_rmnan)), 'normalization', 'probability','binwidth', 5);
    N_18(k, 1:length(N)) = N;
    edges_18(k, 1:length(edges)) = edges;
end 

%reshape elevation data
elev_18_1_vals = reshape(elevations_18_1, 1, []);

%%%Treatment experiment
%%Load treatment data
cd '..\data'
load('ZD_19.mat'); %load topography array (mm). In here should be a 3D topo array called ZD, oriented space x space x time
ZD_19 = ZD_19(:,:,1:280); %same size as control
cd '..\code'

%%Parameters needed for analysis
nx_19 = size(ZD_19,1); %number of x locations on map
ny_19 = size(ZD_19,2); %number of y locations on map
nt_19 = size(ZD_19,3); %number of time steps in data set
dt_19 = 2; %delta t of time steps (hr), dry scans only every 2 horus as compared to every hour for the control

%Loop through times to determine elevation distribution
elevations_19_2 = nan(1000000,280); %initializ matrix for elevations
N_19 = nan(280,21); %initialize matrix for number, 51 is the number of bins 
edges_19 = nan(280,21); %initialzie matrix for edges, 51 is the number of bins
for i=1:nt_19;
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr*dt_19)+ocean_zero; 
    elevationmask = ZD_19(:,:,i); 
    elevationmask_rslr = elevationmask - sl;
    %remove extraneous data
    elevationmask_rslr(elevationmask_rslr > 50) = NaN;
    elevationmask_rslr(elevationmask_rslr < -50) = NaN;
    elevationmask_rslr_rmnan = elevationmask_rslr(~isnan(elevationmask_rslr));
    %remove nan to not bias distribution
    elevations_19_2(1:length(elevationmask_rslr_rmnan),i) = elevationmask_rslr_rmnan;
    [N, edges] = histcounts(elevationmask_rslr_rmnan(~isnan(elevationmask_rslr_rmnan)), 'normalization', 'probability','binwidth', 5);
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

%%plot the mean distribution with standard deviation about the mean
%control
array18 = [elev_bin; mean_prob_18; std_prob_18];
cols = any(isnan(array18),1);
array18(:,cols) = [];
%treatment
array19 = [elev_bin; mean_prob_19; std_prob_19];
cols = any(isnan(array19),1);
array19(:,cols) = [];

%fill standard deviation
y18 = array18(2,:); % your mean vector;
x18 = array18(1,:);
std18 = array18(3,:);
curve1_18 = y18 + std18;
curve2_18 = y18 - std18;

y19 = array19(2,:); % your mean vector;
x19 = array19(1,:);
std19 = array19(3,:);
curve1_19 = y19 + std19;
curve2_19 = y19 - std19;

figure()
plot(x18, y18, 'b-', 'LineWidth', 2)
hold on
plot(x19, y19, 'g-', 'LineWidth', 2)
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
alpha(0.25)
xlabel('elevation relative to sea level (mm)')
xticks(-50:10:50)
ylabel('probability (-)')
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev')
saveas(gcf, '..\figures\Figure2b_elevation_pdf.pdf')

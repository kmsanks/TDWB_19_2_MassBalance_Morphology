clear all; close all
%This script calculates slopes for both the control and
%treatment experiment subaerial marsh window -9 to 5 mm relative to sea level (rsl) and above marsh (> 5 mm rsl). 
%This data is used to create Figure 2c 

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

%%Control experiment slope
%subaerial marsh window
FX_18_marsh = [];
FY_18_marsh = [];
FX_mean_18_marsh = [];
FY_mean_18_marsh = [];
for i=1:nt_18;
    %What is sea level at time i
    sl = (i*baselevel_rr*dt_18)+ocean_zero; 
    elevationmask = ZD_18(:,:,i);
    elevationmask_marsh = elevationmask - sl;
    elevationmask_marsh(elevationmask_marsh > 5) = [NaN];
    elevationmask_marsh(elevationmask_marsh < 0) = [NaN];
    %Slope in x and y direction for each pixel
    [FX_time_18, FY_time_18] = gradient(elevationmask_marsh, 15);
    FX_18_marsh(:,:,i) = FX_time_18;
    FY_18_marsh(:,:,i) = FY_time_18;
    %mean of slopes in x and y direction
    FX_mean_18_marsh(i) = nanmean(nanmean(abs(FX_time_18)));
    FY_mean_18_marsh(i) = nanmean(nanmean(abs(FY_time_18)));
end 

%slope at a pixel is = sqrt(x^2 + y^2)
marsh_slope_18 = sqrt((FX_mean_18_marsh.^2) + (FY_mean_18_marsh.^2));
marsh_mean_slope_18 = nanmean(marsh_slope_18);
marsh_std_slope_18 = nanstd(marsh_slope_18);

%above marsh 
FX_18_am = [];
FY_18_am = [];
FX_mean_18_am = [];
FY_mean_18_am = [];
%loop through time steps
for i=1:nt_18;
    %What is sea level at time i
    sl = (i*baselevel_rr*dt_18)+ocean_zero; 
    elevationmask = ZD_18(:,:,i);
    elevationmask_marsh = elevationmask - sl;
    elevationmask_marsh(elevationmask_marsh <= 5) = [NaN];
    %calculate gradient in x and y for each pixel
    [FX_time_18, FY_time_18] = gradient(elevationmask_marsh, 15);
    FX_18_am(:,:,i) = FX_time_18;
    FY_18_am(:,:,i) = FY_time_18;
    FX_mean_18_am(i) = nanmean(nanmean(abs(FX_time_18)));
    FY_mean_18_am(i) = nanmean(nanmean(abs(FY_time_18)));
end
%calculate slope from x and y gradient data
am_slope_18 = sqrt((FX_mean_18_am.^2) + (FY_mean_18_am.^2));
am_mean_slope_18 = nanmean(am_slope_18);
am_std_slope_18 = nanstd(am_slope_18);

%%%Treatment experiment
%%Load treatment data
cd '..\data'
load('ZD_19.mat'); %load topography array (mm). In here should be a 3D topo array called ZD, oriented space x space x time
cd '..\code'

%%Parameters needed for analysis
nx_19 = size(ZD_19,1); %number of x locations on map
ny_19 = size(ZD_19,2); %number of y locations on map
nt_19 = size(ZD_19,3); %number of time steps in data set
dt_19 = 2; %delta t of time steps (hr), dry scans only every 2 horus as compared to every hour for the control

%loop through timesteps to calculate slope
%subaerial marsh
FX_19_marsh = [];
FY_19_marsh = [];
FX_mean_19_marsh = [];
FY_mean_19_marsh = [];
for i=1:281;
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr*dt_19)+ocean_zero; %scans start at hour 0 so we need to do i-1
    elevationmask = ZD_19(:,:,i);
    elevationmask_marsh = elevationmask - sl;
    elevationmask_marsh(elevationmask_marsh > 5) = [NaN];
    elevationmask_marsh(elevationmask_marsh < 0) = [NaN];
    %calculate x and y gradient
    [FX_time_19, FY_time_19] = gradient(elevationmask_marsh, 15);
    FX_19_marsh(:,:,i) = FX_time_19;
    FY_19_marsh(:,:,i) = FY_time_19;
    FX_mean_19_marsh(i) = nanmean(nanmean(abs(FX_time_19)));
    FY_mean_19_marsh(i) = nanmean(nanmean(abs(FY_time_19)));
end 
%calculate slope
marsh_slope_19 = sqrt((FX_mean_19_marsh.^2) + (FY_mean_19_marsh.^2));
marsh_mean_slope_19 = nanmean(marsh_slope_19);
marsh_std_slope_19 = nanstd(marsh_slope_19);

%above marsh
FX_19_am = [];
FY_19_am = [];
FX_mean_19_am = [];
FY_mean_19_am = [];
%loop through time steps
for i=1:nt_19;
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr*dt_19)+ocean_zero; 
    elevationmask = ZD_19(:,:,i);
    elevationmask_marsh = elevationmask - sl;
    elevationmask_marsh(elevationmask_marsh <= 5) = [NaN];
    %[ASPECT, SLOPE, gradN, gradE] = gradientm(elevationmask_marsh, gridrv);
    [FX_time_19, FY_time_19] = gradient(elevationmask_marsh, 15);
    FX_19_am(:,:,i) = FX_time_19;
    FY_19_am(:,:,i) = FY_time_19;
    FX_mean_19_am(i) = nanmean(nanmean(abs(FX_time_19)));
    FY_mean_19_am(i) = nanmean(nanmean(abs(FY_time_19)));
end 
am_slope_19 = sqrt((FX_mean_19_am.^2) + (FY_mean_19_am.^2));
am_mean_slope_19 = nanmean(am_slope_19);
am_std_slope_19 = nanstd(am_slope_19);

%%Now lets plot the data 
%Boxplot will show difference in slope between the experiments
fig = figure();
slope_grp = [am_slope_18, am_slope_19,...
    marsh_slope_18, marsh_slope_19];
grp = [zeros(1,length(am_slope_18)), ones(1,length(am_slope_19)),...
    2*ones(1,length(marsh_slope_18)), 3*ones(1,length(marsh_slope_19))];
boxplot(slope_grp,grp, 'labels',{'above marsh control', 'above marsh treatment',...
    'subaerial marsh control', 'subaerial marsh treatment'});
ylabel('Delta Top Slope (-)');
saveas(fig, '../figures/Figure2c_slope.pdf')

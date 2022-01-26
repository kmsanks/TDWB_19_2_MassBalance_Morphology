clear all
close all
%%The purpose of this script is to calculate an average radial profile of
%%the experimental delta tops. This sript creates Figure 2d
%% Definitions
%%First, lets set our relative path, we should have opened this file from 
cd '..\' %We should have opened this code from 'C:\Users\XXXXXXX\Documents\GitHub\TDWB_19_2_MassBalance_Morphology\code'

%%Load control data
cd '.\data'
load('TDB_18_data.mat')
cd '..\code'

%We only need the elevation data (ZD)
ZD_18 = Z_maps(:,:,1:560); %remove last 5 hours as that was post-experiment data collection, elevation data is in mm 
ZD_18(ZD_18(:,:,:)==0.)=NaN; %make data outside the basin NaN
chanMaps_18 = C_maps;
chanMaps_18(:,:,298) = ~chanMaps_18(:,:,298); %inverse map
chanMaps_18(:,:,333) = ~chanMaps_18(:,:,333); %inverse map
clear('C_maps','B_maps','G_maps','R_maps','Z_maps');

%%Parameters needed for analysis
nx_18 = size(ZD_18, 1); %number of x locations on map
ny_18 = size(ZD_18,2); %number of y locations on map
nt_18 = size(ZD_18,3); %number of time steps in data set
dx = 5; %5 mm grid cells in x
dy = 5; %5 mm grid cells in y
dt_18 = 1; %delta t of time steps (hr)
baselevel_rr = 0.25; %base level rise rate (mm/hr)
ocean_zero = 25; %ocean elevation at beginning of experiment (mm)
%channel entrance
xentrance_18 = 109; %x grid node location of the entrance channel
yentrance_18 = 271; %y grid node location of the entrance channel

%%Load treatment data
%%Load treatment data
cd '..\data'
load('ZD_19_2_dry.mat'); %load topography array. In here should be a 3D topo array called ZD, oriented space x space x time
load('TDWB_19_2_chanMaps.mat');
cd '..\code'
%change topography data to match the name of control experiment
ZD_19 = ZD_19_2_dry(:,:,1:280);
chanMaps_19 = TDWB_19_2_chanMaps;
chanMaps_19 = chanMaps_19(:,:,2:2:end); %every other to match elevation data 
clear('ZD_19_2_dry','TDWB_19_2_chanMaps');

%%Parameters needed for analysis
nx_19 = size(ZD_19,1); %number of x locations on map
ny_19 = size(ZD_19,2); %number of y locations on map
nt_19 = size(ZD_19,3); %number of time steps in data set
dt_19 = 2; %delta t of time steps (hr), dry scans only every 2 horus as compared to every hour for the control
%channel entrance
xentrance_19 = 214; %x grid node location of the entrance channel (x is down dip)
yentrance_19 = 397; %y grid node location of the entrance channel (y is strike)

%% Remove Channels from topo
%Make elevation screen, so flow can be clipped to the areas above sea
%level
RSLR18_screen = [];
for i= 1:nt_18;
    %What is sea level at time i
    sl = (i*baselevel_rr*dt_18)+25; %first wet scan starts 48 minutes into the hour so i = i is a good approximation
    elevationmask = ZD_18(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    %elevationmask_rslr(elevationmask_rslr < -30) = NaN; %
    RSLR18_screen(:,:,i) = elevationmask_rslr;
end

RSLR19_screen = [];
for i= 1:nt_19;
    %What is sea level at time i
    sl = ((i-1)*baselevel_rr*dt_19)+25; %first wet scan starts 48 minutes into the hour so i = i is a good approximation
    elevationmask = ZD_19(:,:,i); 
    elevationmask(elevationmask == 0) = NaN;
    elevationmask_rslr = elevationmask - sl;
    %elevationmask_rslr(elevationmask_rslr < -30) = NaN;
    RSLR19_screen(:,:,i) = elevationmask_rslr;
end

%Need channel pixels to be 0 but currently they are 1
chanMaps_18 = ~chanMaps_18;
chanMaps_19 = ~chanMaps_19;

% Remove channels 
z18 = [];
for i =1:nt_18
  z18_RSL=(chanMaps_18(:,:,i)).*(RSLR18_screen(:,:,i));
  z18_RSL(z18_RSL == 0) = NaN; %0 is flow area 
  z18(:,:,i) = z18_RSL;
end

z19 = [];
for i =1:nt_19
  z19_RSL=(chanMaps_19(:,:,i)).*(RSLR19_screen(:,:,i));
  z19_RSL(z19_RSL == 0) = NaN; %0 is flow area 
  z19(:,:,i) = z19_RSL;
end

%% Main loop
elev_18 = [];%empty array to be filled with elevation from radial transects every 5 mm from source 
elev_std_18 = [];%empty array to be filled with std elevation from radial transects every 5 mm from source 
elev_25prct_18 = [];%fill loop for 1.5 iqr
elev_75prct_18 = [];
%Loop through radial transects
for i = 1:700;%loop to run through different radial distances from the end of the entrance channel.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    th = 0:1/i:2*pi;
    xunit = i * cos(th) + xentrance_18;
    yunit = i * sin(th) + yentrance_18;
    xunit = round(xunit);
    yunit = round(yunit);
    xunit18 = [];
    yunit18 = [];
    for j = 1:max(size(xunit))
        xc = xunit(j);
        yc = yunit(j);
        if xc >= 1;
            if yc >= 1;
                if xc <= nx_18
                    if yc <= ny_18
                        xunit18 = [xunit18;xc];
                        yunit18 = [yunit18;yc];
                    end
                end
            end
        end
    end
    xunit = xunit18;
    yunit = yunit18;
    XY18 = [yunit xunit];
    XY18 = sortrows(XY18);
    xunit = XY18(:,2);
    yunit = XY18(:,1);
    
    zs18 = []; 
    for j = 1:nanmax(size(xunit));
        xs_shot = z18(xunit(j),yunit(j),:);
        zs18 = [zs18;xs_shot];
    end
    zs18 = squeeze(zs18);%matrix of topo along radial transect
    %fill matrices with data for each radial transect 
    elev_18(:,i) = nanmedian(zs18);%empty array to be filled with elevation from radial transects
    elev_25prct_18(:,i) = prctile(zs18,25);
    elev_25prct_18(:,i) = prctile(zs18,25);
    elev_std_18(:,i) = nanstd(zs18);
    radial_dist_18(:,i) = i*5; %mm because each pixel is 5mm in length
end

%Statistics of mean and standard deviation for each radial transect
elev_18_mean = nanmean(elev_18);
elev_18_stdmean = nanstd(elev_18);

%%Treatment experiment
elev_19 = [];%empty array to be filled with elevation from radial transects every 5 mm from source 
elev_std_19 = [];%empty array to be filled with std elevation from radial transects every 5 mm from source 
elev_25prct_19 = [];
elev_75prct_19 = [];
%Loop through radial transects
for i = 1:666;%loop to run through different radial distances from the end of the entrance channel.
    %section to find x,y nodes for radial transect and generate matrix of
    %cross section topo and channel (yes/no) data
    th = 0:1/i:2*pi;
    xunit = i * cos(th) + xentrance_19;
    yunit = i * sin(th) + yentrance_19;
    xunit = round(xunit);
    yunit = round(yunit);
    xunit19 = [];
    yunit19 = [];
    for j = 1:max(size(xunit))
        xc = xunit(j);
        yc = yunit(j);
        if xc >= 1;
            if yc >= 1;
                if xc <= nx_19
                    if yc <= ny_19
                        xunit19 = [xunit19;xc];
                        yunit19 = [yunit19;yc];
                    end
                end
            end
        end
    end
    xunit = xunit19;
    yunit = yunit19;
    XY19 = [yunit xunit];
    XY19 = sortrows(XY19);
    xunit = XY19(:,2);
    yunit = XY19(:,1);
    
    zs19 = []; 
    for j = 1:nanmax(size(xunit));
        xs_shot = z19(xunit(j),yunit(j),:);
        zs19 = [zs19;xs_shot];
    end
    zs19 = squeeze(zs19);%matrix of topo along radial transect

    elev_19(:,i) = nanmedian(zs19);%empty array to be filled with elevation from radial transects
    elev_25prct_19(:,i) = prctile(zs19,25);
    elev_75prct_19(:,i) = prctile(zs19,75);
    elev_std_19(:,i) = nanstd(zs19);
    radial_dist_19(:,i) = i*5; %mm
end 

elev_19_mean = nanmean(elev_19);
elev_19_stdmean = nanstd(elev_19);


%% Plot the data
ybars = [-9 5]; %marsh window
%create errorbar about mean
%array of distance, mean, stdev
%control
array18 = [radial_dist_18; elev_18_mean; elev_18_stdmean];
cols = any(isnan(array18),1);
array18(:,cols) = [];
%treatment
array19 = [radial_dist_19; elev_19_mean; elev_19_stdmean];
cols = any(isnan(array19),1);
array19(:,cols) = [];

%fill standard deviation
y18 = array18(2,:); % your mean vector;
x18 = array18(1,:); %distance
std18 = array18(3,:); %standard deviation
curve1_18 = y18 + std18;
curve2_18 = y18 - std18;

y19 = array19(2,:); % your mean vector;
x19 = array19(1,:); %distance
std19 = array19(3,:); %standard deviation
curve1_19 = y19 + std19;
curve2_19 = y19 - std19;

%Plot the data
fig = figure()
plot(x18, y18, 'b', 'LineWidth', 2)
hold on
plot(x19, y19, 'g', 'LineWidth', 2)
patch([x18 fliplr(x18)], [curve1_18 fliplr(curve2_18)], 'b')
patch([x19 fliplr(x19)], [curve1_19 fliplr(curve2_19)], 'g')
patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], 'k')
alpha(0.15)
yline(0, 'k-', 'linewidth', 2)
plot(x18, y18, 'b', 'LineWidth', 2)
plot(x19, y19, 'g', 'LineWidth', 2)
ylim([-40 50])
xlim([0 2000])
ylabel('elevation relative to sea level (mm)')
xlabel('radial distance from channel entrance (mm)') 
legend('control mean', 'treatment mean', 'control stdev', 'treatment stdev', 'marsh window', 'sea level')
set(gcf, 'PaperUnits', 'inches');
y_width=7.25;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '..\figures\Figure2d_radialelevation.pdf')

clear all; close all
%This script calculates delta top area over time for four different
%categorys (1. above sea level (>= 0 mm relative to sea level (rsl)), 2. total delta top (>= -9 mm rsl), 3. marsh window area (-9 mm to 5
%mm rsl), and 4. above marsh window (> 5 mm rsl). For the manuscript, we
%plot only 2. total delta top and 3. marsh window area in Fiugre 2a.

%%First, lets set our relative path
cd '../' %We should have opened this code from 'C:\Users\XXXXXXX\Documents\GitHub\TDWB_19_2_MassBalance_Morphology\code'

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

%%Area Analysis
%Create empty data frames
subaerial_area_18 = []; % >= 0 mm rsl
total_area_18 = []; % >= -9 rsl
marsh_area_18 = []; %-9 to 5 mm rsl
am_area_18 = []; %> 5 mm rsl
elevationmask_contour_18 = []; %elevations relative to sea level
%Loop through all 560 hours in control
for i=1:nt_18;
    %What is sea level at time i
    sl_18 = ((i*baselevel_rr*dt_18)+ocean_zero); %18 scans start at hour 1 so subtract 1 from i 
    elevationmask_18 = ZD_18(:,:,i); 
    elevationmask_18(elevationmask_18 == 0) = NaN;
    
    %subaerial delta
    elevationmask_rslr_18 = elevationmask_18 - sl_18;
    elevationmask_rslr_18(elevationmask_rslr_18 < 0) = NaN;
    subaerial_area_18(:,i) = length(elevationmask_rslr_18(~isnan(elevationmask_rslr_18)));
    
    %total delta top 
    elevationmask_rslr_am_18 = elevationmask_18 - sl_18;
    elevationmask_rslr_am_18(elevationmask_rslr_am_18 < -9) = NaN;
    total_area_18(:,i) = length(elevationmask_rslr_am_18(~isnan(elevationmask_rslr_am_18)));
    
    %marsh window
    elevationmask_rslr_marsh_18 = elevationmask_18 - sl_18;
    elevationmask_rslr_marsh_18(elevationmask_rslr_marsh_18 < -9) = NaN;
    elevationmask_rslr_marsh_18(elevationmask_rslr_marsh_18 >= 5) = NaN;
    marsh_area_18(:,i) = length(elevationmask_rslr_marsh_18(~isnan(elevationmask_rslr_marsh_18)));
    
    %above marsh 
    elevationmask_rslr_marsh_18 = elevationmask_18 - sl_18;
    elevationmask_rslr_marsh_18(elevationmask_rslr_marsh_18 < 5) = NaN;
    am_area_18(:,i) = length(elevationmask_rslr_marsh_18(~isnan(elevationmask_rslr_marsh_18)));
    
    %elevation mask relative to sea level
    elevationmask_contour_18(:,:,i) = elevationmask_18 - sl_18;
end 

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

%%Area Analysis
%Create empty data frames
subaerial_area_19 = []; % >= 0 mm rsl
total_area_19 = []; % >= -9 rsl
marsh_area_19 = []; %-9 to 5 mm rsl
am_area_19 = []; %> 5 mm rsl;
elevationmask_contour_18 = []; %elevations relative to sea level
%Loop through all hours in the treatment
for i=1:nt_19
    %What is sea level at time i
    sl_19 = ((i-1)*baselevel_rr*dt_19)+ocean_zero; %19 scans start at hour 0 so need to subtract 1 from hour
    elevationmask_19 = ZD_19(:,:,i); 
    elevationmask_19(elevationmask_19 == 0) = NaN;
    
    %subaerial 
    elevationmask_rslr_19 = elevationmask_19 - sl_19;
    elevationmask_rslr_19(elevationmask_rslr_19 < 0) = NaN;
    subaerial_area_19(:,i) = length(elevationmask_rslr_19(~isnan(elevationmask_rslr_19)));
    
    %total delta top
    elevationmask_rslr_am_19 = elevationmask_19 - sl_19;
    elevationmask_rslr_am_19(elevationmask_rslr_am_19 < -9) = NaN;
    total_area_19(:,i) = length(elevationmask_rslr_am_19(~isnan(elevationmask_rslr_am_19)));
    
    %marsh window
    elevationmask_rslr_marsh_19 = elevationmask_19 - sl_19;
    elevationmask_rslr_marsh_19(elevationmask_rslr_marsh_19 < -9) = NaN;
    elevationmask_rslr_marsh_19(elevationmask_rslr_marsh_19 >= 5) = NaN;
    marsh_area_19(:,i) = length(elevationmask_rslr_marsh_19(~isnan(elevationmask_rslr_marsh_19)));
    
    %above marsh 
    elevationmask_rslr_marsh_19 = elevationmask_19 - sl_19;
    elevationmask_rslr_marsh_19(elevationmask_rslr_marsh_19 < 5) = NaN;
    am_area_19(:,i) = length(elevationmask_rslr_marsh_19(~isnan(elevationmask_rslr_marsh_19)));

    %mask for contour
    elevationmask_contour_19(:,:,i) = elevationmask_19 - sl_19;
end

%Convert to meters
%subaerial
AerialAreaMeter_18 = subaerial_area_18*(2.5*10^-5); %z data is stored in mm
AerialAreaMeter_19 = subaerial_area_19*(2.5*10^-5);
%total delta top
TotalAreaMeter_18 = total_area_18*(2.5*10^-5);
TotalAreaMeter_19 = total_area_19*(2.5*10^-5);
%marshwindow
MarshAreaMeter_18 = marsh_area_18*(2.5*10^-5);
MarshAreaMeter_19 = marsh_area_19*(2.5*10^-5);
%above marsh
AmAreaMeter_18 = am_area_18*(2.5*10^-5);
AmAreaMeter_19 = am_area_19*(2.5*10^-5);

%%Statistics
%Calculate mean and standard deviation
%subaerial 
mean_area_subaerial_18 = nanmean(AerialAreaMeter_18);
mean_area_subaerial_19 = nanmean(AerialAreaMeter_19);
std_area_subaerial_18 = nanstd(AerialAreaMeter_18);
std_area_subaerial_19 = nanstd(AerialAreaMeter_19);
%total delta top
mean_area_total_18 = nanmean(TotalAreaMeter_18);
mean_area_total_19 = nanmean(TotalAreaMeter_19);
std_area_total_18 = nanstd(TotalAreaMeter_18);
std_area_total_19 = nanstd(TotalAreaMeter_19);
%marsh window
mean_area_marsh_18 = nanmean(MarshAreaMeter_18);
mean_area_marsh_19 = nanmean(MarshAreaMeter_19);
std_area_marsh_18 = nanstd(MarshAreaMeter_18);
std_area_marsh_19 = nanstd(MarshAreaMeter_19);
%above marsh
mean_area_am_18 = nanmean(AmAreaMeter_18);
mean_area_am_19 = nanmean(AmAreaMeter_19);
std_area_am_18 = nanstd(AmAreaMeter_18);
std_area_am_19 = nanstd(AmAreaMeter_19);

%calculate the percent land gained in treatment versus control
percent_increase_total = ((mean_area_total_19 - mean_area_total_18)/mean_area_total_18)*100;
percent_increase_marsh = ((mean_area_marsh_19 - mean_area_marsh_18)/mean_area_marsh_18)*100;
percent_increase_subaerial = ((mean_area_subaerial_19 - mean_area_subaerial_18)/mean_area_subaerial_18)*100;
percent_increase_am = ((mean_area_am_19 - mean_area_am_18)/mean_area_am_18)*100;

%%marsh fraction
%control
marshfrac18 = []; %initiate empty matrix
for i = 1:nt_18;
    marshfrac = MarshAreaMeter_18(i)/TotalAreaMeter_18(i); %marsh frac = marsh area/delta area
    marshfrac18 = [marshfrac18, marshfrac];
end
%treatment
marshfrac19 = [];
for i = 1:nt_19;
    marshfrac = MarshAreaMeter_19(i)/TotalAreaMeter_19(i);
    marshfrac19 = [marshfrac19, marshfrac];
end

%mean and standard deviation for marsh fraction to report in manuscript
%control
meanmarshfrac18 = mean(marshfrac18);
stdmarshfrac18 = std(marshfrac18);
%treatment
meanmarshfrac19 = mean(marshfrac19);
stdmarshfrac19 = std(marshfrac19);

%%Now, lets plot the data
t_19 = linspace(0,562,281); %x-axis, hours
t_18 = 1:nt_18; %hours

%Plot area figure
fig = figure;
%marsh window control
plot(t_18,MarshAreaMeter_18, 'c-');
yline(mean_area_marsh_18, 'c--');
hold on
%marsh window treatment
plot(t_19,MarshAreaMeter_19, 'g-');
yline(mean_area_marsh_19, 'g--');
hold on
%total delta top control
plot(t_18, TotalAreaMeter_18, 'b-');
yline(mean_area_total_18, 'b--');
%total delta top treatment
plot(t_19, TotalAreaMeter_19, 'color',[0 0.5 0], 'linestyle', '-');
yline(mean_area_total_19, 'color',[0 0.5 0], 'linestyle', '-');
ylabel('Area (m^2)', 'FontSize',12,...
       'FontWeight','bold','Color','k')
xlabel('Time (hours)', 'FontSize',12,...
       'FontWeight','bold','Color','k')
ylim([0,4.5])
legend({'marsh window (control)', 'mean (control)','marsh window (treatment)', 'mean (treatment)',...
    'delta top (control)', 'mean (control)', 'delta top (treatment)', 'mean (treatment)'},...
    'Orientation','Horizontal', 'NumColumns', 2,'Location', 'north');
set(gcf, 'PaperUnits', 'inches');
y_width=7.25 ;x_width=9.125
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'..\figures\Figure2a_area.pdf')
hold off 

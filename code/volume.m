clear all 
close all
%%The purpose of this script is to calculate the clastic and marsh volume
%%accumulated on the delta top (>= -9 mm relative to sea level (rsl)), marsh window (-9 to 5 mm
%%rsl), above the marsh window (>5 mm rsl), and off shore (< -9 mm rsl)
%%throught the experiments. This script will produce the data for Table 1
%%For trapping efficiencies in Table 1 and Figure 3 (a and b), please run the code: "tapping_efficiency.m" 

%%%Control experiment
%%Load control and treatment data
cd '..\data'
load('ZD_18.mat') %topography; elevation data (mm) control
load('ZD_19.mat') %topography; elevation data (mm) treatment

%Load stratigraphric interpolated marsh fraction
load('interpolated_marsh_strat_frac.mat');
cd '..\code'

%%Parameters needed for analysis
nt_18 = size(ZD_18,3); %number of scans control
nt_19 = size(ZD_19,3); %number of scans treatment
dt_18 = 1; %time between scans control
dt_19 = 2; %time between scans treatment
baselevel_rr = 0.25; %mm/hr; background sea level rise rate
ocean_zero = 25; %25 mm; starting ocean level

%%Calculate area above marsh window (> 5 mm rsl)
%control
pland_am_18 = []; %empty matrix for binary
for i = 1:nt_18;
    z = ZD_18(:,:,i); %topo
    z = z - ((baselevel_rr*i*dt_18)+ocean_zero); %topo relative to sea level
    z(z < 0) = NaN;
    z(z <= 5) = NaN;
    z(z > 5) = 1; %only area above 5 mm rsl
    pland_am_18(:,:,i) = z; %above marsh binary
end 

%we need to apply a 90% filter to each image
%we only want pixels that were above the marsh for at least 90% of the
%experiment to minimize marsh effects in this area
pland_am_18_sum = nansum(pland_am_18,3); %sum binary along time
pland_am_18_sum(pland_am_18_sum == 0) = NaN; %remove data outside basin
pland_am_18_frac = pland_am_18_sum/560; %fraction of time pixel is above marsh

%binary of 1s and NaNs
pland_am_18_frac(pland_am_18_frac >= 0.9) = 1.; %filter for >=90% of time
pland_am_18_frac(pland_am_18_frac < 1) = NaN; %remove all fractions

%Now lets multiply the binary matrix for pixels above marsh by the
%elevation data; this will eventually give us volume
z0_18 = ZD_18(:,:,1).*pland_am_18_frac;
z_last_18 = ZD_18(:,:,560).*pland_am_18_frac; %last - initial = volume (dz; mm)
dz_am_18 = (z_last_18 - z0_18)/1000; %dz in m
v_am_18 = nansum(dz_am_18(:))*(2.5E-5); %multiply by area (one pixel is 5mmx5mm = 2.5E-5m)

%%Calculate area above marsh window (> 5 mm rsl)
%treatment
pland_am_19 = []; %empty matrix for land fraction
for i = 1:nt_19;
    z = ZD_19(:,:,i); %topo
    z = z - (baselevel_rr*(i-1)*dt_19+ocean_zero); %subtract 1 because scans start at hour 0
    z(z < 0) = NaN;
    z(z <= 5) = NaN;
    z(z > 5) = 1;
    pland_am_19(:,:,i) = z;
end 

%we need to apply a 90% filter to each image
%we only want pixels that were above the marsh for at least 90% of the
%experiment to minimize marsh effects in this area
pland_am_19_sum = nansum(pland_am_19,3);
pland_am_19_sum(pland_am_19_sum == 0) = NaN;
pland_am_19_frac = pland_am_19_sum/281; %281 hours for control

%binary of 1s and NaNs
pland_am_19_frac(pland_am_19_frac >= 0.9) = 1.;
pland_am_19_frac(pland_am_19_frac < 1) = NaN;

%Now lets multiply the binary matrix for pixels above marsh by the
%elevation data; this will eventually give us volume
z0_19 = ZD_19(:,:, 1).*pland_am_19_frac;
z_last_19 = ZD_19(:,:,281).*pland_am_19_frac;
dz_am_19 = (z_last_19 - z0_19)/1000; %total dz including marsh
dz_am_19_clastic = ((z_last_19 - z0_19).*(1-interpolated_marsh_strat_frac))/1000; %dz in m for clastic material only
v_am_19 = nansum(dz_am_19(:))*(2.5E-5);
v_am_19_clastic = nansum(dz_am_19_clastic(:))*(2.5E-5);

%Now we will do the same but for the marsh window (-9 to 5 mm rsl)
%To create an area that doesn't over lap the above marsh window, we will
%use any pixel in the marsh window for at least 10% of the time
%control
pland_marsh_18 = [];
for i = 1:nt_18;
    z = ZD_18(:,:,i);
    z = z - (baselevel_rr*i*dt_18+ocean_zero);
    z(z >= -9) = 1;
    z(z < -9) = NaN;
    pland_marsh_18(:,:,i) = z;
end

%we need to apply the 10% filter to each image
pland_marsh_18_sum = nansum(pland_marsh_18,3);
pland_marsh_18_sum(pland_marsh_18_sum == 0) = NaN;
pland_marsh_18_frac = pland_marsh_18_sum/560;
pland_marsh_18_frac(pland_marsh_18_frac > 0.1) = 1.;
pland_marsh_18_frac(pland_marsh_18_frac < 1) = NaN;

%remove above marsh from marsh window, so that the above marsh continues
%right into the marsh and we have no overlap 
tmp = cat(3,pland_am_18_frac, pland_marsh_18_frac);
marsh18 = nansum(tmp, 3);
marsh18(marsh18 > 1) = NaN;
marsh18(marsh18 < 1) = NaN;
%rename
pland_marsh_18_frac = marsh18;

z0_18 = ZD_18(:,:,1).*pland_marsh_18_frac;
z_last_18 = ZD_18(:,:,560).*pland_marsh_18_frac;
dz_marsh_18 = (z_last_18 - z0_18)/1000; %dz in m
v_marsh_18 = nansum(dz_marsh_18(:))*(2.5E-5); %volume in m

%Treatment marsh window
pland_marsh_19 = [];
for i = 1:nt_19;
    z = ZD_19(:,:,i);
    z = z - (baselevel_rr*(i-1)*dt_19+ocean_zero);%
    z(z >= -9) = 1;
    z(z < -9) = NaN;
    pland_marsh_19(:,:,i) = z;
end 

%we need to apply the 50% filter to each image
pland_marsh_19_sum = nansum(pland_marsh_19,3);
pland_marsh_19_sum(pland_marsh_19_sum == 0) = NaN;
pland_marsh_19_frac = pland_marsh_19_sum/281;
pland_marsh_19_frac(pland_marsh_19_frac > 0.1) = 1.;
pland_marsh_19_frac(pland_marsh_19_frac < 1) = NaN;

%remove spot that is clearly not in the marsh window
pland_marsh_19_frac(680:750,360:515) = NaN;

%remove above marsh from marsh window, so that the above marsh continues
%right into the marsh and we have no overlap 
tmp = cat(3,pland_am_19_frac, pland_marsh_19_frac);
marsh19 = nansum(tmp, 3);
marsh19(marsh19 > 1) = NaN;
marsh19(marsh19 < 1) = NaN;
%rename
pland_marsh_19_frac = marsh19;

z0_19 = ZD_19(:,:,1).*pland_marsh_19_frac;
z_last_19 = ZD_19(:,:,281).*pland_marsh_19_frac;
dz_marsh_19 = (z_last_19 - z0_19)/1000; %dz in m
dz_marsh_19_clastic = ((z_last_19 - z0_19).*(1-interpolated_marsh_strat_frac))/1000; %dz in m
v_marsh_19 = nansum(dz_marsh_19(:))*(2.5E-5);
v_marsh_19_clastic = nansum(dz_marsh_19_clastic(:))*(2.5E-5);

%Now we will calculate volume accumulated on the delta top (>=-9 mm rsl)
%for area that is in this region for at least 50% of the experiment to
%compare average conditions between the two experiments
%control
pland_dtop_18 = [];
for i = 1:nt_18;
    z = ZD_18(:,:,i);
    z = z - (baselevel_rr*i*dt_18+ocean_zero);
    z(z >= -9) = 1;
    z(z < -9) = NaN;
    pland_dtop_18(:,:,i) = z;
end 

%we need to apply the 50% filter to each image
pland_dtop_18_sum = nansum(pland_dtop_18,3);
pland_dtop_18_sum(pland_dtop_18_sum == 0) = NaN;
pland_dtop_18_frac = pland_dtop_18_sum/560;
pland_dtop_18_frac(pland_dtop_18_frac >= 0.5) = 1.;
pland_dtop_18_frac(pland_dtop_18_frac < 1) = NaN;

%volume
z0_18 = ZD_18(:,:,1).*pland_dtop_18_frac;
z_last_18 = ZD_18(:,:,560).*pland_dtop_18_frac;
dz_dtop_18 = (z_last_18 - z0_18)/1000; %dz in m
v_dtop_18 = nansum(dz_dtop_18(:))*(2.5E-5);

%Now we will calculate volume accumulated on the delta top (>=-9 mm rsl)
%for area that is in this region for at least 50% of the experiment to
%compare average conditions between the two experiments
%treatment
pland_dtop_19 = [];
for i = 1:nt_19;
    z = ZD_19(:,:,i);
    z = z - (baselevel_rr*(i-1)*dt_19+ocean_zero);
    z(z >= -9) = 1;
    z(z < -9) = NaN;
    pland_dtop_19(:,:,i) = z;
end 

%we need to apply the 50% filter to each image
pland_dtop_19_sum = nansum(pland_dtop_19,3);
pland_dtop_19_sum(pland_dtop_19_sum == 0) = NaN;
pland_dtop_19_frac = pland_dtop_19_sum/281;
pland_dtop_19_frac(pland_dtop_19_frac >= 0.5) = 1.;
pland_dtop_19_frac(pland_dtop_19_frac < 1) = NaN;

%calculate volume balance
z0_19 = ZD_19(:,:,1).*pland_dtop_19_frac;
z_last_19 = ZD_19(:,:,281).*pland_dtop_19_frac;
dz_dtop_19 = (z_last_19 - z0_19)/1000; %dz in m
dz_dtop_19_clastic = ((z_last_19 - z0_19).*(1-interpolated_marsh_strat_frac))/1000; %dz in m
v_dtop_19 = nansum(dz_dtop_19(:))*(2.5E-5);
v_dtop_19_clastic = nansum(dz_dtop_19_clastic(:))*(2.5E-5);

%%Now lets analyze the data
%turn matrix into vector
dz_air_18 = dz_dtop_18(:)*1000/560;;
dz_air_19 = dz_dtop_19_clastic(:)*1000/559;
dz_air_18 = dz_air_18(~isnan(dz_air_18));
dz_air_19 = dz_air_19(~isnan(dz_air_19));

dz_above_18 = dz_am_18(:)*1000/560;;
dz_above_19 = dz_am_19_clastic(:)*1000/559;
dz_above_18 = dz_above_18(~isnan(dz_above_18));
dz_above_19 = dz_above_19(~isnan(dz_above_19));

dz_m_18 = dz_marsh_18(:)*1000/560;;
dz_m_19 = dz_marsh_19_clastic(:)*1000/559;
dz_m_18 = dz_m_18(~isnan(dz_m_18));
dz_m_19 = dz_m_19(~isnan(dz_m_19));

%%THIS IS FOR A PLOT IN "trapping_efficiency.m" SCRIPT
%Rerun marsh window calculation
%control
pland_marshplot_18 = [];
for i = 1:nt_18;
    z = ZD_18(:,:,i);
    z = z - (baselevel_rr*i*dt_18+ocean_zero);
    z(z >= -9) = 1;
    z(z < -9) = NaN;
    pland_marshplot_18(:,:,i) = z;
end 

%apply 10% filter
pland_marshplot_18_sum = nansum(pland_marshplot_18,3);
pland_marshplot_18_sum(pland_marshplot_18_sum == 0) = NaN;
pland_marshplot_18_frac = pland_marshplot_18_sum/560;
pland_marshplot_18_frac(pland_marshplot_18_frac > 0.1) = 1.;
pland_marshplot_18_frac(pland_marshplot_18_frac < 1) = NaN;

%treatment
pland_marshplot_19 = [];
for i = 1:nt_19;
    z = ZD_19(:,:,i);
    z = z - (baselevel_rr*(i-1)*dt_19+ocean_zero);
    z(z >= -9) = 1;
    z(z < -9) = NaN;
    pland_marshplot_19(:,:,i) = z;
end 

%apply 10% filter
pland_marshplot_19_sum = nansum(pland_marshplot_19,3);
pland_marshplot_19_sum(pland_marshplot_19_sum == 0) = NaN;
pland_marshplot_19_frac = pland_marshplot_19_sum/281;
pland_marshplot_19_frac(pland_marshplot_19_frac > 0.1) = 1.;
pland_marshplot_19_frac(pland_marshplot_19_frac < 1) = NaN;
%remove spot that is clearly not in the marsh window
pland_marshplot_19_frac(680:750,360:515) = NaN;

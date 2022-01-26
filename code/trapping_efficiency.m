%The purpose of this script is to calculate the off-shore volume
%accumulated and the trapping efficiency for each region in the volume.m
%script
%This script produces figure 3 (a and b) and trapping efficiency data in
%Table 1

run('volume.m');

%convert river sed to volume
%3.1 lbs per hour is ~1406.14 g/hr
river_flux = 1406.14; %g/hr
time = 560; %hours
bulk_density = (1-0.55)*2.65; %density of sediment assuming 55% porosity and particle density of 2.65 g/cm^3
mt = ((river_flux/bulk_density)*time)*10^-6; %mass in m^3

%sediment past above masrsh area (above marsh is anything above marsh >=
%90% of the time)
river_flux_marsh_18 = mt - v_am_18;
river_flux_marsh_19 = mt - v_am_19_clastic; %we only trap clastic sediment

%volume accumulated off-shore
v_offshore_18 = mt - v_dtop_18; %m^3
v_offshore_19 = mt - v_dtop_19_clastic; %m^3

%trapping efficiency (clastic volume accumulated/volume delivered); see supporting information for more detail 
te_am_18 = v_am_18/mt;
te_am_19 = v_am_19_clastic/mt;
te_dtop_18 = v_dtop_18/mt;
te_dtop_19 = v_dtop_19_clastic/mt;
te_marsh_18 = v_marsh_18/mt;
te_marsh_19 = v_marsh_19_clastic/mt;
te_basin_18 = v_offshore_18/mt;
te_basin_19 = v_offshore_19/mt;

%Now lets calculate the marsh trapping efficiency using only bypassed
%sediment; this is the second calculation in Table 1 for the marsh window
%and is represented in the footnote a
te_marsh_18_deliver = v_marsh_18/river_flux_marsh_18;
te_marsh_19_deliver = v_marsh_19_clastic/river_flux_marsh_19;

%fraction of volume of marsh on total delta top 
frac_marsh = 1 - (v_dtop_19_clastic/v_dtop_19);

%%area for volume calculation of above marsh and marsh window figure 3
%add binarys together
tmp = cat(3,pland_am_18_frac, pland_marshplot_18_frac);
area_sep18 = nansum(tmp, 3);
tmp = cat(3,pland_am_19_frac, pland_marshplot_19_frac);
area_sep19 = nansum(tmp, 3);

%figure of volume and areas contoured
fig = figure();
subplot(1,2,1)
imagesc(area_sep18(109:645,1:522))
axis equal
colorbar
caxis([0 2])
subplot(1,2,2)
imagesc(area_sep19(214:750,130:651))
axis equal
colorbar
caxis([0 2])
set(gcf, 'PaperUnits', 'inches');
y_width=11 ;x_width=8.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(fig, '..\figures\Figure3_partitioning.pdf')


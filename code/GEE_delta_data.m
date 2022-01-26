clear all
close all
%The purpose of this script is to compile a global (ish) data set of 
%deltaic elevations. Here we will process elevation data for four gobal deltas:
%1. Mississippi River Delta (MRD), Ganges Brahmaputra Meghna River Delta
%(GBMD), 3. Mekong River Delta (mekong), and 4. Rio Grande River Delta (riogrande)
%downloaded from hand-drawn polygons in Google Earth Engine (GEE)

%Note: this script will not run out of the box, please read through code
%and instructions. However, the matrices made here are already stored in
%the data folder, so you don't need to run this script if you are not
%interested in how those matrices were made

%%Process elevation data from GEE
%MRD
elev_MRD = -31:2:89; %-1 to 3 channel depths relative to sea level, using a 2 m bin size
MRD = zeros(61,4); %empty data
MRD(:,1) = elev_MRD'; %transpose elevation bins
%copy and paste counts from excel into column 2 (raw excel data can be
%found in the data folder
MRD(:,3) = MRD(:,1)/30; %normalize by chnnel depth (30m) 
MRD(:,4) = MRD(:,2)/nansum(MRD(:,2));
save '..\data\MRD.mat' MRD;

%GBMD
elev_GBMD = -31:2:89;
GBMD = zeros(61,4);
GBMD(:,1) = elev_GBMD';
%copy and paste counts from excel into column 2
GBMD(:,3) = GBMD(:,1)/30; %normalize by chnnel depth (30 m) 
GBMD(:,4) = GBMD(:,2)/nansum(GBMD(:,2)); 
save '..\data\GBMD.mat' GBMD;

%Mekong
elev_mekong = -11:2:31;
mekong = zeros(22,4);
mekong(:,1) = elev_mekong';
%copy and paste counts from excel into column 2
mekong(:,3) = mekong(:,1)/10; %normalize by chnnel depth (10 m)
mekong(:,4) = mekong(:,2)/nansum(mekong(:,2)); 
save '..\data\mekong.mat' mekong;

%Rio Grande
elev_riogrande = -15:2:45;
riogrande = zeros(31,4);
riogrande(:,1) = elev_riogrande';
%copy and paste counts from excel into column 2
riogrande(:,3) = riogrande(:,1)/15; %normalize by chnnel depth (60 m)
riogrande(:,4) = riogrande(:,2)/nansum(riogrande(:,2));
save '..\data\riogrande.mat' riogrande;


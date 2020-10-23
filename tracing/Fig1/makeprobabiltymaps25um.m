%%% script to make Med anterior, Med posterior, Interposed and Lateral CN
%%% projection probability maps.

%%
clear all
close all

%% make projection probability maps
datadir='./';
load('regions.mat')
DCNnames={'FN','FN_round10','IN','DN'};

% for d=1:numel(DCNnames)
for d=3
load([DCNnames{d},'_Matlaboutput/TI_rawcounts.mat'])
if d==1
   tmp=P;
   P=[];
   P(1).rawdatafull=tmp(1).rawdatafull;
   P(2).rawdatafull=tmp(2).rawdatafull;
   P(3).rawdatafull=tmp(4).rawdatafull;
   P(4).rawdatafull=tmp(5).rawdatafull;
   P(5).rawdatafull=tmp(6).rawdatafull;
end
if d==2
   tmp=P;
   P=[];
   P(1).rawdatafull=tmp(1).rawdatafull;
   P(2).rawdatafull=tmp(3).rawdatafull;
   P(3).rawdatafull=tmp(4).rawdatafull;
   P(4).rawdatafull=tmp(5).rawdatafull;
   P(5).rawdatafull=tmp(7).rawdatafull;
end
if d==3
   tmp=P;
   P=[];
   P(1).rawdatafull=tmp(1).rawdatafull;
   P(2).rawdatafull=tmp(2).rawdatafull;
   P(3).rawdatafull=tmp(4).rawdatafull;
   P(4).rawdatafull=tmp(5).rawdatafull;
   P(5).rawdatafull=tmp(6).rawdatafull;
   P(6).rawdatafull=tmp(7).rawdatafull;
end


% set cutoff parameter
cutoff=20000;

% make ventricle mask
ventricles=findchildren(73,annotated,pathId);
mask=ismember(regions,ventricles);

% make probability heatmaps


[prob,smoothedprob]=calculateprojectionprobability_norm(P,cutoff,['smooth_test_',DCNnames{d}],mask);
end

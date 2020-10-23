%%%%% analyse the output of regionQuant

%%
clear all
close all

%% parameters
datadir='./';

load([datadir,'ZI_Matlaboutput\NormalizedRegionalDensity.mat']);

dataL=NormalizedRegionalDensityL;
dataR=NormalizedRegionalDensityR;
DCN=ones(1,size(NormalizedInnervationL,2));


load([datadir,'Ret_Matlaboutput\NormalizedRegionalDensity.mat']);
dataL=[dataL,NormalizedRegionalDensityL];
dataR=[dataR,NormalizedRegionalDensityR];
DCN=[DCN,2.*ones(1,size(NormalizedInnervationL,2))];

%% filter data

%for analysis, I'd like to mask out cortex, corpus callosum, layer 1 of
%superior colliculus, and hypothalamus. possibly the striatum too.
% 
% 
for i = 1:height(annotated)
    a = str2double(strsplit(annotated.structure_id_path{i},'/'));
    pathId{i} = a(~isnan(a));
end

[isoctx]=findchildren(315,annotated,pathId);
[hypothalamus]=findchildren(1097,annotated,pathId);
[corpuscallosum]=findchildren(776,annotated,pathId);
[SC_sensory]=findchildren(302,annotated,pathId); % might be able to reduce this to only the zonal layer, or superfical grey.
fibertracts=findchildren(1009,annotated,pathId);
ventricles=findchildren(73,annotated,pathId);
striatum=findchildren(477,annotated,pathId);

root=997;
DCNids=findchildren(519,annotated,pathId);
% % combine all these areas tp be ignored into one mask
mask=[isoctx,corpuscallosum,fibertracts,ventricles,root,striatum];

maskedareasL=ismember(annotated_collapsed.id,mask);
maskedareasR=ismember(annotated_collapsed.id,[mask,DCNids]);


annotated_filteredL=annotated_collapsed(~isnan(sum(dataR,2)) & ~maskedareasL,:);
annotated_filteredR=annotated_filteredL(~ismember(annotated_filteredL.id,DCNids),:);

%split by hemisphere and remove injection site
data_l=dataL(~isnan(sum(dataL,2)) & ~maskedareasL,:);
data_r=dataR(~isnan(sum(dataR,2)) & ~maskedareasL & ~ismember(annotated_collapsed.id,DCNids),:);


DCN_r=DCN;
DCN_l=DCN;


%% do pca

[coeff,score,latent,tsquared,explained,mu] = pca([data_r;data_l]);
figure;
plot(coeff(DCN_r==1,1),coeff(DCN_r==1,2),'o');hold on
plot(coeff(DCN_r==2,1),coeff(DCN_r==2,2),'o');
xlabel('PC 1')
ylabel('PC 2')
title('PCA 1v2')
legend('ZI','Ret');




%% run stats.

%right hemisphere
[pvalue_R,names_R,names_Fs_R,sortedmean_R,sorteddata_R]=regionttests(data_r,DCN_r,annotated_filteredR);

%left hemisphere
[pvalue_L,names_L,names_Fs_L,sortedmean_L,sorteddata_L]=regionttests(data_l,DCN_l,annotated_filteredL);

%save to csv

tableR=table(names_R,pvalue_R,sortedmean_R(:,1),sortedmean_R(:,2),...
    sorteddata_R(:,1),sorteddata_R(:,2),sorteddata_R(:,3),sorteddata_R(:,4),...
    sorteddata_R(:,5),sorteddata_R(:,6),sorteddata_R(:,7));
tableR.Properties.VariableNames = {'Area','pvalue','ZI_mean','Ret_mean',...
    'ZI_mouse1','ZI_mouse2','ZI_mouse3','ZI_mouse4',...
    'Ret_mouse1','Ret_mouse2','Ret_mouse3'};
writetable(tableR,'ZI_Ret_densities.xls','Sheet','ipsi')  

tableL=table(names_L,pvalue_L,sortedmean_L(:,1),sortedmean_L(:,2),...
    sorteddata_L(:,1),sorteddata_L(:,2),sorteddata_L(:,3),sorteddata_L(:,4),...
    sorteddata_L(:,5),sorteddata_L(:,6),sorteddata_L(:,7));
tableL.Properties.VariableNames = {'Area','pvalue','ZI_mean','Ret_mean',...
    'ZI_mouse1','ZI_mouse2','ZI_mouse3','ZI_mouse4',...
    'Ret_mouse1','Ret_mouse2','Ret_mouse3'};
writetable(tableL,'ZI_Ret_densities.xls','Sheet','contra')  


%% sort by average projection strength
sortandplot(data_l,DCN_l,annotated_filteredL.name,50,{'ZI contra','Ret contra'})
sortandplot(data_r,DCN_r,annotated_filteredR.name,50,{'ZI ipsi','Ret ipsi'})

plotsortedheatmap(data_l,data_r,DCN_l,[50],annotated_filteredL,{'ZI','Ret'})

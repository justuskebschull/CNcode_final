%%%%% analyse the output of regionQuant

%%
clear all
close all

%% parameters
datadir='./';
DCNnames={'FN','FN_round10','IN','DN'};
dataL=[];
dataR=[];
DCN=[];

%load data and drop out brains with poor axon detection performance.

for i=1:numel(DCNnames)
    load([datadir,DCNnames{i},'_Matlaboutput/NormalizedRegionalDensity.mat']);
    if i==2
        dataL=[dataL,NormalizedRegionalDensityL(:,[1,3:5,7])];
        dataR=[dataR,NormalizedRegionalDensityR(:,[1,3:5,7])];
        DCN=[DCN,i.*ones(1,size(NormalizedInnervationL(:,[1,3:5,7]),2))];
    elseif i==1
        dataL=[dataL,NormalizedRegionalDensityL(:,[1,2,4:6])];
        dataR=[dataR,NormalizedRegionalDensityR(:,[1,2,4:6])];
        DCN=[DCN,i.*ones(1,size(NormalizedInnervationL(:,[1,2,4:6]),2))];
    elseif i==3
        dataL=[dataL,NormalizedRegionalDensityL(:,[1,2,4:7])];
        dataR=[dataR,NormalizedRegionalDensityR(:,[1,2,4:7])];
        DCN=[DCN,i.*ones(1,size(NormalizedInnervationL(:,[1,2,4:7]),2))];
    else
        dataL=[dataL,NormalizedRegionalDensityL];
        dataR=[dataR,NormalizedRegionalDensityR];
        DCN=[DCN,i.*ones(1,size(NormalizedInnervationL,2))];
    end
end


%% filter data


for i = 1:height(annotated)
    a = str2double(strsplit(annotated.structure_id_path{i},'/'));
    pathId{i} = a(~isnan(a));
end

[isoctx]=findchildren(315,annotated,pathId);
[corpuscallosum]=findchildren(776,annotated,pathId);
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
for i=1:numel(unique(DCN_r))
plot(coeff(DCN_r==i,1),coeff(DCN_r==i,2),'o');hold on
end
xlabel('PC 1')
ylabel('PC 2')
title('PCA 1v2')
legend(DCNnames);
print(gcf,['figure_exports/PCA_12'],'-dpdf','-r1000','-fillpage')

figure;
for i=1:numel(unique(DCN_r))
plot(coeff(DCN_r==i,2),coeff(DCN_r==i,3),'o');hold on
end
xlabel('PC 2')
ylabel('PC 3')
title('PCA 2v3')
legend(DCNnames);
print(gcf,['figure_exports/PCA_23'],'-dpdf','-r1000','-fillpage')

% interestingly brain 215 and 219 are quite different, both previously labeled as very
% sparse.. remove from analysis.

% plot in 3D
figure;
for i=1:numel(unique(DCN_r))
scatter3(coeff(DCN_r==i,1),coeff(DCN_r==i,2),coeff(DCN_r==i,3),'o');hold on
end
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
title('PCA')
legend(DCNnames);
print(gcf,['figure_exports/PCA_123'],'-dpdf','-r1000','-fillpage')

%% do tsne
Y=tsne([data_r;data_l]');
gscatter(Y(:,1),Y(:,2),DCN_r)


%% run stats.
[p_r,p_multcomp_r,sorteddata_r,sortednames_r,sortedp_r]=region_anova(data_r,DCN_r,annotated_filteredR,0.05,'ipsi hemi 0.05');
ax=gca;
ax.CLim=[0 8e-7];
colorbar;
print(gcf,['figure_exports/anovasorted_ipsi'],'-dpdf','-r1000','-fillpage')


[p_l,p_multcomp_l,sorteddata_l,sortednames_l,sortedp_l]=region_anova(data_l,DCN_l,annotated_filteredL,0.05,'contra hemi 0.05');
ax=gca;
ax.CLim=[0 8e-7];
colorbar;
print(gcf,['figure_exports/anovasorted_contra'],'-dpdf','-r1000','-fillpage')

tableR=table(sortednames_r,sortedp_r',...
    sorteddata_r(:,DCN_r==1),sorteddata_r(:,DCN_r==2),sorteddata_r(:,DCN_r==3),sorteddata_r(:,DCN_r==4));
tableR.Properties.VariableNames = {'Area','Anovapvalue','Medial_anterior','Medial_posterior','Interposed','Lateral'};
writetable(tableR,'anterograde_densities.xls','Sheet','ipsi_pcutoff')  


tableL=table(sortednames_l,sortedp_l',...
    sorteddata_l(:,DCN_l==1),sorteddata_l(:,DCN_l==2),sorteddata_l(:,DCN_l==3),sorteddata_l(:,DCN_l==4));
tableR.Properties.VariableNames = {'Area','Anovapvalue','Medial_anterior','Medial_posterior','Interposed','Lateral'};
writetable(tableL,'anterograde_densities.xls','Sheet','contra_pcutoff')  


%save to csv
[p_r,p_multcomp_r,sorteddata_r,sortednames_r,sortedp_r]=region_anova(data_r,DCN_r,annotated_filteredR,1,'ipsi hemi 1');

[p_l,p_multcomp_l,sorteddata_l,sortednames_l,sortedp_l]=region_anova(data_l,DCN_l,annotated_filteredL,1,'contra hemi 1');



tableR=table(sortednames_r,sortedp_r',...
    sorteddata_r(:,DCN_r==1),sorteddata_r(:,DCN_r==2),sorteddata_r(:,DCN_r==3),sorteddata_r(:,DCN_r==4));
tableR.Properties.VariableNames = {'Area','Anovapvalue','Medial_anterior','Medial_posterior','Interposed','Lateral'};
writetable(tableR,'anterograde_densities.xls','Sheet','ipsi_all')  


tableL=table(sortednames_l,sortedp_l',...
    sorteddata_l(:,DCN_l==1),sorteddata_l(:,DCN_l==2),sorteddata_l(:,DCN_l==3),sorteddata_l(:,DCN_l==4));
tableR.Properties.VariableNames = {'Area','Anovapvalue','Medial_anterior','Medial_posterior','Interposed','Lateral'};
writetable(tableL,'anterograde_densities.xls','Sheet','contra_all')  



%% sort by average projection strength
sortandplot(data_l,DCN_l,annotated_filteredL.name,50,DCNnames)
sortandplot(data_r,DCN_r,annotated_filteredR.name,50,DCNnames)
%% plot sorted heatmaps of top 50 targetes per hemisphere.
numareas=plotsortedheatmap(data_l,data_r,DCN_l,[50],annotated_filteredL,DCNnames,1);
%% figure out the number of projection targets per hemisphere.
numareas=plotsortedheatmap(data_l,data_r,DCN_l,[],annotated_filteredL,DCNnames,0)

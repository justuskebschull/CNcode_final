%%%%% analyse the output of regionQuant

%%
clear all
close all

%% parameters
datadir='./';

load([datadir,'adipo_cre_Matlaboutput/NormalizedRegionalDensity.mat']);

dataL=NormalizedRegionalDensityL;
dataR=NormalizedRegionalDensityR;
DCN=ones(1,size(NormalizedInnervationL,2));

% %annotate by injection site.
% %1=SC;2=spinalcord;3=crus1;4=RN;5=PN;6=vest.;7=vta,8=cm thal;
% %9=vermis,10=VLthal
% DCN=[1,1,2,3,4,5,4,4,5,6,3,5,7,8,9,9,1,10,3,10];


%annotate by injection site.
%1=VLthal;2=cm thal;3=SC;4=RN;5=PN;6=vta;7=vest.;8=vermis,9=crus1;10=spinalcord;
%
DCN=[3,3,10,9,4,5,4,4,5,7,9,5,6,2,8,8,3,1,9,1];

%% filter out bad brains
badbrains=[1,13,17];
indx=ones(size(DCN));
indx(badbrains)=0;
dataL=dataL(:,logical(indx));
dataR=dataR(:,logical(indx));
DCN=DCN(logical(indx));

%sort data.
[DCN,i]=sort(DCN);
dataL=dataL(:,i);
dataR=dataR(:,i);


%% filter data


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
for i=unique(DCN)
plot(coeff(DCN_r==i,1),coeff(DCN_r==i,2),'o');hold on
end
xlabel('PC 1')
ylabel('PC 2')
title('PCA 1v2')
legend('VLthal','CMthal','SC','RN','PN','VestNuc','Vermis','Crus1','spinalcord');

print(gcf,['figure_exports/PCA'],'-dpdf','-r1000')


figure;
for i=unique(DCN)
plot(coeff(DCN_r==i,2),coeff(DCN_r==i,3),'o');hold on
end
xlabel('PC 2')
ylabel('PC 3')
title('PCA 2v3')
legend('VLthal','CMthal','SC','RN','PN','VestNuc','Vermis','Crus1','spinalcord');

%%
[p_r,p_multcomp_r,sorteddata_r,sortednames_r,sortedp_r]=region_anova(data_r,DCN_r,annotated_filteredR,1,'ipsi hemi all');
ax=gca;
ax.CLim=[0 8e-7];
colorbar;
print(gcf,['figure_exports/anovasorted_ipsi'],'-dpdf','-r1000','-fillpage')


[p_l,p_multcomp_l,sorteddata_l,sortednames_l,sortedp_l]=region_anova(data_l,DCN_l,annotated_filteredL,1,'contra hemi all');
ax=gca;
ax.CLim=[0 8e-7];
colorbar;
print(gcf,['figure_exports/anovasorted_contra'],'-dpdf','-r1000','-fillpage')



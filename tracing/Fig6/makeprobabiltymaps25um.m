%% Script that takes the output of regionQuantwholebrain_collapsed_round9.m 
% and derives Class A and Class B probability maps.

%%
clear all
close all


%% load data
load('round9_rawdataZI.mat')
PZI=P;

load('round9_rawdataRet.mat')
PRet=P;

%% set cutoff parameter
cutoff=20000;

%% make ventricle mask
for i = 1:height(annotated)
    a = str2double(strsplit(annotated.structure_id_path{i},'/'));
    pathId{i} = a(~isnan(a));
end

ventricles=findchildren(73,annotated,pathId);
mask=ismember(regions,ventricles);

%% make smoothed heatmaps
% ZI=makesmoothedoutput(PZI,cutoff,'smoothZI',mask);
% Ret=makesmoothedoutput(Pret,cutoff,'smoothRet',mask);

%% lets solve the system of linear equations of P(ZI) and P(Ret) for P(A) and P(B)
%first calculate proper projection probability maps

[ZIprob,ZIsmoothedprob]=calculateprojectionprobability(PZI,cutoff,'smoothZI',mask);
[Retprob,Retsmoothedprob]=calculateprojectionprobability(PRet,cutoff,'smoothRet',mask);

% variables from retro starmap experiment. P(ZI)=xP(A)+(1-x)P(B) and
% P(Ret)=yP(A)+(1-y)P(B)
load('C:\Users\justus\Documents\postdoc2\DCN_sequencing\scRNAseq\STARmap\RetroStarmap\Round2\infected.mat')
infected2=infected(4:7,:,20:28);
infected2(1,:,:)=infected2(1,:,[4,3,6,5,8,7,1,9,2]);
infected2(2,:,:)=infected2(2,:,[4,3,6,5,8,7,1,9,2]);
infected2(3,:,:)=infected2(3,:,[9,7,6,8,5,3,2,4,1]);
infected2(4,:,:)=infected2(4,:,[9,7,6,8,5,3,2,4,1]);

infected3=infected2(:,[13,14,15],[4,7]);
infected4(:,1,:)=infected3(:,1,:);
infected4(:,2,:)= infected3(:,2,:)+infected3(:,3,:);

for t=1:2
    ratio(t,:)=squeeze(infected4(:,1,t))./(squeeze(infected4(:,1,t))+squeeze(infected4(:,2,t)));
end
x=mean(ratio(1,:));
y=mean(ratio(2,:),'omitnan');


[A,B]=getprobability(ZIsmoothedprob,Retsmoothedprob,x,y);
Ascaled=uint16(A.*(~mask)*2^16);
Bscaled=uint16(B.*(~mask)*2^16);

mkdir('classA')
mkdir('classB')

for i=1:size(A,3)
    imwrite(squeeze(Ascaled(:,:,i)),['classA/classA',int2str(i),'.tif']);
    imwrite(squeeze(Bscaled(:,:,i)),['classB/classB',int2str(i),'.tif']);
end


%% let's define areas where B>>A to calculate second order projections from thalamus.
cutofflo=0.3;
cutoffhi=0.5;
Bdiff=B-A;
Bspecificlo=Bdiff>cutofflo;
Bspecifichi=B>3.*A & B>0.3;

Bspecificlo_scaled=uint16(Bspecificlo.*~mask*2^16);
Bspecifichi_scaled=uint16(Bspecifichi.*~mask*2^16);

mkdir 'classBscificlo'
mkdir 'classBspecifichi'

for i=1:size(A,3)
    imwrite(squeeze(Bspecificlo_scaled(:,:,i)),['classBscificlo/classBscificlo',int2str(i),'.tif']);
    imwrite(squeeze(Bspecifichi_scaled(:,:,i)),['classBspecifichi/classBspecifichi',int2str(i),'.tif']);
end

%% let's define areas where B<<A to calculate second order projections from thalamus.
cutofflo=0.3;
cutoffhi=0.5;
Adiff=A-B;
Aspecificlo=Adiff>cutofflo;
Aspecifichi=A>3.*B & A>cutofflo;

Aspecificlo_scaled=uint16(Aspecificlo.*~mask*2^16);
Aspecifichi_scaled=uint16(Aspecifichi.*~mask*2^16);

mkdir 'classAscificlo'
mkdir 'classAspecifichi'

for i=1:size(A,3)
    imwrite(squeeze(Aspecificlo_scaled(:,:,i)),['classAscificlo/classAscificlo',int2str(i),'.tif']);
    imwrite(squeeze(Aspecifichi_scaled(:,:,i)),['classAspecifichi/classAspecifichi',int2str(i),'.tif']);
end




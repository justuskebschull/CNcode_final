%% collect all data into a usefull datastructure

allbrains=[];
DCN=[];
k=1;
datadir='Y:\justus\iDisco\round8\brainsforanalysis\FN';
load([datadir,'\Matlab_output_full\ilastik_rawcounts.mat']);
DCN=[DCN,ones(1,numel(P))];
for i=1:numel(P)
allbrains(:,:,:,k)=P(i).rawdatafull;
k=k+1;
end


datadir='Y:\justus\iDisco\round8\brainsforanalysis\IN';
load([datadir,'\Matlab_output_full\ilastik_rawcounts.mat']);
DCN=[DCN,2*ones(1,numel(P))];
for i=1:numel(P)
allbrains(:,:,:,k)=P(i).rawdatafull;
k=k+1;
end

datadir='Y:\justus\iDisco\round8\brainsforanalysis\DN';
load([datadir,'\Matlab_output_full\ilastik_rawcounts.mat']);
DCN=[DCN,3*ones(1,numel(P))];
for i=1:numel(P)
allbrains(:,:,:,k)=P(i).rawdatafull;
k=k+1;
end


%% treshold
cutoff=1000;
allbrains2=allbrains;
allbrains2(allbrains<cutoff)=0;
%% normalize and smoothe

for i=1:size(allbrains,4)
    allbrains_norm(:,:,:,i)=squeeze(allbrains2(:,:,:,i))./sum(sum(sum(squeeze(allbrains2(:,:,:,i)))));
    allbrains3(:,:,:,i)=imgaussfilt3(squeeze(allbrains_norm(:,:,:,i)),2);
end


%% save mean values per DCN.
FNtmp=mean(allbrains_norm(:,:,:,DCN==1),4,'omitnan');
INtmp=mean(allbrains_norm(:,:,:,DCN==2),4,'omitnan');
DNtmp=mean(allbrains_norm(:,:,:,DCN==3),4,'omitnan');

FN=imgaussfilt3(FNtmp,2);
IN=imgaussfilt3(INtmp,2);
DN=imgaussfilt3(DNtmp,2);


mkdir FN_heatmap
n=max(max(max(FN)));
for i=1:size(FN,3)
imwrite(uint16(FN(:,:,i)./n.*2^16),['FN_heatmap/FN_',int2str(i),'.tif']);
end
mkdir IN_heatmap
n=max(max(max(IN)));
for i=1:size(IN,3)
imwrite(uint16(IN(:,:,i)./n.*2^16),['IN_heatmap/IN_',int2str(i),'.tif']);
end
mkdir DN_heatmap
n=max(max(max(DN)));
for i=1:size(DN,3)
imwrite(uint16(DN(:,:,i)./n.*2^16),['DN_heatmap/DN_',int2str(i),'.tif']);
end


%% do ttests


[~,p23]=ttest2(permute(double(allbrains3(:,:,:,DCN==2)),[4,1,2,3]),permute(double(allbrains3(:,:,:,DCN==3)),[4,1,2,3]),'Tail','right');
p23=squeeze(p23);
p23(isnan(p23))=1; %remove nans
% p23(p23>0.05)=1;
p23=1-p23;
p23_scaled=uint16(p23.*2^16);
mkdir p23
for i=1:size(p23,3)
imwrite(p23_scaled(:,:,i),['p23/p23test_',int2str(i),'.tif']);
end


[~,p32]=ttest2(permute(double(allbrains3(:,:,:,DCN==2)),[4,1,2,3]),permute(double(allbrains3(:,:,:,DCN==3)),[4,1,2,3]),'Tail','left');
p32=squeeze(p32);
p32(isnan(p32))=1; %remove nans
% p32(p32>0.05)=1;
p32=1-p32;
p32_scaled=uint16(p32.*2^16);
mkdir p32
for i=1:size(p32,3)
imwrite(p32_scaled(:,:,i),['p32/p32test_',int2str(i),'.tif']);
end


%%

[~,p13]=ttest2(permute(double(allbrains3(:,:,:,DCN==1)),[4,1,2,3]),permute(double(allbrains3(:,:,:,DCN==3)),[4,1,2,3]),'Tail','right');
p13=squeeze(p13);
p13(isnan(p13))=1; %remove nans
% p23(p23>0.05)=1;
p13=1-p13;
p13_scaled=uint16(p13.*2^16);
mkdir p13
for i=1:size(p13,3)
imwrite(p13_scaled(:,:,i),['p13/p13test_',int2str(i),'.tif']);
end


[~,p31]=ttest2(permute(double(allbrains3(:,:,:,DCN==1)),[4,1,2,3]),permute(double(allbrains3(:,:,:,DCN==3)),[4,1,2,3]),'Tail','left');
p31=squeeze(p31);
p31(isnan(p31))=1; %remove nans
% p32(p32>0.05)=1;
p31=1-p31;
p31_scaled=uint16(p31.*2^16);
mkdir p31
for i=1:size(p31,3)
imwrite(p31_scaled(:,:,i),['p31/p31test_',int2str(i),'.tif']);
end

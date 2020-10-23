%%
clear all
close all



%% collect all data into a usefull datastructure

allbrains=[];
DCN=[];
k=1;

load('round9_rawdataZI.mat');
DCN=[DCN,ones(1,numel(P))];
for i=1:numel(P)
allbrains(:,:,:,k)=P(i).rawdatafull;
k=k+1;
end


load('round9_rawdataRet.mat');
DCN=[DCN,2*ones(1,numel(P))];
for i=1:numel(P)
allbrains(:,:,:,k)=P(i).rawdatafull;
k=k+1;
end



%% treshold
cutoff=20000;
allbrains2=allbrains;
allbrains2(allbrains<cutoff)=0;
%% normalize and smoothe

for i=1:size(allbrains,4)
    allbrains_norm(:,:,:,i)=squeeze(allbrains2(:,:,:,i))./sum(sum(sum(squeeze(allbrains2(:,:,:,i)))));
    allbrains3(:,:,:,i)=imgaussfilt3(squeeze(allbrains_norm(:,:,:,i)),2);
end



%% do ttests


[~,p21]=ttest2(permute(double(allbrains3(:,:,:,DCN==1)),[4,1,2,3]),permute(double(allbrains3(:,:,:,DCN==2)),[4,1,2,3]),'Tail','right');
p21=squeeze(p21);
p21(isnan(p21))=1; %remove nans
% p23(p23>0.05)=1;
p21=1-p21;
p21_scaled=uint16(p21.*2^16);
mkdir p21
for i=1:size(p21,3)
imwrite(p21_scaled(:,:,i),['p21/p21test_',int2str(i),'.tif']);
end

p21_thresh=uint16(double(p21>0.99).*2^16);
mkdir p21_binary
for i=1:size(p21,3)
imwrite(p21_thresh(:,:,i),['p21_binary/p21test_',int2str(i),'.tif']);
end




[~,p12]=ttest2(permute(double(allbrains3(:,:,:,DCN==2)),[4,1,2,3]),permute(double(allbrains3(:,:,:,DCN==1)),[4,1,2,3]),'Tail','right');
p12=squeeze(p12);
p12(isnan(p12))=1; %remove nans
% p32(p32>0.05)=1;
p12=1-p12;
p12_scaled=uint16(p12.*2^16);
mkdir p12
for i=1:size(p12,3)
imwrite(p12_scaled(:,:,i),['p12/p12test_',int2str(i),'.tif']);
end



p12_thresh=uint16(double(p12>0.99).*2^16);
mkdir p12_binary
for i=1:size(p12,3)
imwrite(p12_thresh(:,:,i),['p12_binary/p12test_',int2str(i),'.tif']);
end
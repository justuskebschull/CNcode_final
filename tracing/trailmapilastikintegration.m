%% script to go through folder and max combine TRAILMAP and Ilastik segmentation before any analysis

%% parameters
datadir='/home/drew/Documents/synology/iDisco/round8/brainsforanalysis/DN';
brains = dir([datadir,'/brain*']);
for i = 1:size(brains,1)
    ILASTIKdatafolder{i} = [datadir,'/',brains(i).name,'/ilastik_oct2019/'];
    tmp=dir([datadir,'/',brains(i).name,'/seg-647_*']);
    TRAILMAPdatafolder{i} = [datadir,'/',brains(i).name,'/',tmp(1).name,'/'];
    outputfolder{i}= [datadir,'/',brains(i).name,'/seg-TI'];
end

%% multiply all brains in folder by 3 to match background levels to those seen in TRAILMAP training set.

for i = [1:3,5:size(brains,1)]
    tmp=dir([datadir,'/',brains(i).name,'/*647*'])
    raw=dir([datadir,'/',brains(i).name,'/',tmp(1).name,'/*.tif']);
    mkdir([datadir,'/',brains(i).name,'/647_X3']);
    parfor s=1:size(raw,1)
        tmp=(imread([raw(s).folder,'/',raw(s).name]).*3);
        imwrite(tmp,[datadir,'/',brains(i).name,'/647_X3/scaled647_',num2str(s,'%04.f'),'.tif']);
    end
end


%% loop over brains and max combine each slice.
%for i = 1:size(brains,1)
    for i=3
    ilastik=dir([datadir,'/',brains(i).name,'/ilastik_oct2019/*.tif']);
    trailmap=dir([char(TRAILMAPdatafolder(i)),'*.tif']);
    mkdir(char(outputfolder(i)));
    parfor s=1:size(ilastik,1)
%         hInfo = imfinfo([trailmap(s).folder,'/',trailmap(s).name]);
%         ImageHeight = hInfo(1).Height; 
%         ImageWidth = hInfo(1).Width;        
%         tiffobject=Tiff([trailmap(s).folder,'/',trailmap(s).name],'r');
%         ImageMatrix = zeros(ImageHeight,ImageWidth,'uint32');
%         ImageMatrix(:,:)=tiffobject.read();

        tmp1=(imread([ilastik(s).folder,'/',ilastik(s).name]));
        tmp2=double(imread([trailmap(s).folder,'/',trailmap(s).name])).*65535;
        tmp3=zeros([2,size(tmp1)]);
        tmp3(1,:,:)=tmp1;
        tmp3(2,:,:)=tmp2;
        tmp4=uint16(squeeze(max(tmp3,[],1)));
        imwrite(tmp4,[char(outputfolder(i)),'/TI_',num2str(s,'%04.f'),'.tif'])
    end
end
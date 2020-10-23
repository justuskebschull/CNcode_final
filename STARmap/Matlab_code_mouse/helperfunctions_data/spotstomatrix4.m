function [outputmatrix,genelabels,cells]=spotstomatrix4(nissl,spotfile,rawfile,section,cutoff,ploting)
%find spotfiles
files=dir([spotfile,section,'*.csv']);
genelabels=[];
cells=bwlabel(nissl);
outputmatrix=zeros(length(files),max(max(cells))+1);

%make corresponding list of rawfiles
rawfiles=dir([rawfile,'*.tif']);


%plot spots on nissl image;
if ploting==1
figure;imshow(nissl);
axis on;
hold on;
end

%loop over spotfiles
for i=1:length(files)
%load spot file
spots=readtable([spotfile,files(i).name]);
spotcoords=spots{:,[3,4]};
gene=char(spots{1,10});

%filter spots by intensity
rawspots=imread([rawfile,rawfiles(i).name]);
files(i).name
rawfiles(i).name
spotintensities=[];
for k=1:size(spotcoords,1)
    spotintensities(k)=rawspots(spotcoords(k,2)+1,spotcoords(k,1)+1);
end

channel_tmp=split(files(i).name,'_');
channel_tmp2=split(channel_tmp(3),'');
channel=str2double(cell2mat(channel_tmp2(3)));

spotcoords=spotcoords(spotintensities>cutoff(channel),:);

%plot spots on nissl image;
if ploting==1
plot(spotcoords(:,1),spotcoords(:,2),'.','Color',rand(1,3));
end
%associate spots with cells;
celloforigin=[];
for j=1:size(spotcoords,1)
    if spotcoords(j,2)<size(cells,1) 
        celloforigin(j)=cells(spotcoords(j,2)+1,spotcoords(j,1)+1); %is this correct?? holds up to santity check.
    else
        celloforigin(j)=0;
    end
end

%put into matrix. first column is outofcell
outputmatrix(i,:)=hist(celloforigin,0:1:max(max(cells)));
genelabels(i).name=gene;
end

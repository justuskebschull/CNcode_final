%% run separately for every section chopped into NxN pieces

function [outputmatrix,genelabels,cells]=spotstomatrix_chopped(nissl,spotfile,rawfile,N,section,cutoff)
%find spotfiles
files=dir([spotfile,'/*.csv']);
genelabels=[];
cells=bwlabel(nissl);


%make corresponding list of rawfiles
rawfiles=dir([rawfile,'*.tif']);

%chop nissl/cell image into same chunks as the data was separated in for
%spotfinding
cells_dim=round(size(cells)./N);
chopinfo=imfinfo([rawfile,rawfiles(1).name]);

choppedcells=[];
k=1; %counter
for r=1:N
    for c=1:N
        choppedcells(k).c=cells(((r-1)*chopinfo.Height+1) : (r*chopinfo.Height), ...
                                    ((c-1)*chopinfo.Width+1) : (c*chopinfo.Width));
        k=k+1;
    end
end
    

outputmatrix=zeros(length(files)./(N^2),max(max(cells))+1);




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



%associate spots with cells;

%determine correct nissl fragment to compare.
name_tmp=char(files(i).name);
frag=str2double(name_tmp(2));

celloforigin=[];
for j=1:size(spotcoords,1)
    if spotcoords(j,2)<size(choppedcells(frag).c,1) 
        celloforigin(j)=choppedcells(frag).c(spotcoords(j,2)+1,spotcoords(j,1)+1); %is this correct?? holds up to santity check.
    else
        celloforigin(j)=0;
    end
end

%put into matrix. first column is outofcell
cellspresent=unique(choppedcells(frag).c);
counts=hist(celloforigin,cellspresent);

%which gene are we looking at?
G=strsplit(gene);
gi=str2double(G(2));

outputmatrix(gi,cellspresent+1)=counts;
genelabels(gi).name=gene;
end

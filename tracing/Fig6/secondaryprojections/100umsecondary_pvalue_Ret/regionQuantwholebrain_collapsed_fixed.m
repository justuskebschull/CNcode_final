%%% analyse the output of second order projection tracing using the te
%%% allen voxel level connectivity model

clear all
close all

%% load annotation files
Regionsfull=nrrdread('../../../../atlasfiles/100um_atlas/annotation_100_coronal.nrrd');
annotated = readtable('../../annotation_info_0118_1327_justusedits.csv');
sz = size(Regionsfull);

%% collapse regions according to annotation table.
for i = 1:height(annotated)
    a = str2double(strsplit(annotated.structure_id_path{i},'/'));
    pathId{i} = a(~isnan(a));
end

%combineids = a list of region ids that have all children collapsed
%for the "full" annotation file, only 73 and 1009 are collapsed
%for the "collapsed" annotation file, see the 0108_1327.csv
combineids = annotated.id(logical(annotated.collapse_logical));
regions=Regionsfull;
for i = 1:numel(pathId)
    if any(ismember(pathId{i},combineids))
        
        tmp=pathId{i}(ismember(pathId{i},combineids)); 
        
        regions(regions == annotated.id(i)) = tmp(end); %this allows to include areas like nucleus of darkschewitsch and fields of forrel that are embedded within their parent strucutre.
        i
    end
end

annotated_collapsed=annotated(annotated.collapse_logical==1,:);



%% read in brain
% load brain image, and save intensity values of pixels in each annotation
%region 
% the midline is at 177 in cropped7_525
PixelbyRegionR = cell(height(annotated_collapsed),1);
PixelbyRegionL = cell(height(annotated_collapsed),1);
H=height(annotated_collapsed);

s = imfinfo('../p12_binary-Thalamus_projections.tif');
rawdataL = (zeros(s(1).Height,s(1).Width/2,numel(s)));
rawdataR = (zeros(s(1).Height,s(1).Width/2,numel(s)));
rawdatafull = (zeros(s(1).Height,s(1).Width,numel(s)));

    for i = 1:numel(s) 
        tmp=imread('../p12_binary-Thalamus_projections.tif',i);
        rawdataL(:,:,i) = tmp(:,1:57);
        rawdataR(:,:,i) = tmp(:,58:end);
        rawdatafull(:,:,i) = tmp;
    end

%normalize data by starting voxels
startingvoxels=69;
rawdatafull=rawdatafull./startingvoxels;
rawdataL=rawdataL./startingvoxels;
rawdataR=rawdataR./startingvoxels;
    
    
   for i = 1:H
           PixelbyRegionL{i} = rawdataL(regions(:,1:57,:) == annotated_collapsed.id(i)); %each saves intensity value for each pixel in a region
           PixelbyRegionR{i} = rawdataR(regions(:,58:end,:) == annotated_collapsed.id(i)); %each saves intensity value for each pixel in a region
   end


%%
% calculate region by region axon innervations



%sum up pixelsVALUES in each region and normalize by region size
  

for i = 1:length(PixelbyRegionL)
    RegionalDensityL(i) = sum(PixelbyRegionL{i})/numel(PixelbyRegionL{i}); 
    RegionalDensityR(i) = sum(PixelbyRegionR{i})/numel(PixelbyRegionR{i}); 
end

%%%%%%alternatively do not normalize by region size

for i=1:length(PixelbyRegionL)
    AxonsByRegionL(i) = sum(PixelbyRegionL{i});
    AxonsByRegionR(i) = sum(PixelbyRegionR{i});
    end



% save
save((['NormalizedRegionalDensity.mat']),...
    'RegionalDensityL','RegionalDensityR',...
    'AxonsByRegionL','AxonsByRegionR',...
    'annotated','PixelbyRegionL','PixelbyRegionR');%save output




%% remove "thalamic" areas. kind of self projections and plot bargraph.
thal=findchildren(549,annotated,pathId);
i=ismember(annotated_collapsed.id,thal);
A=AxonsByRegionR(~i);
B=RegionalDensityR(~i);
Nt=annotated_collapsed.name(~i);

n=20;

[dens,ix]=sort(A,'descend');
nms=Nt(ix);

figure;
bar(dens(1:n))
set(gca, 'XTick', 1:n, 'XTickLabel', nms(1:n));
xtickangle(45);
    title('AxonsByRegion')



tmp=isnan(B);
D=B(~tmp);
N=Nt(~tmp);
[dens,ix]=sort(D,'descend');
nms=N(ix);

figure;
bar(dens(1:n))
set(gca, 'XTick', 1:n, 'XTickLabel', nms(1:n));
xtickangle(45);
    title('RegionalDensity')




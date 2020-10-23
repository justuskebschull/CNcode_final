

clear all
close all

% apparently Fiji messes up the annotation file if its cropped and then
% saved in there. so do the cropping in matlab instead.

%% load annotation files
Regionsfull=nrrdread('../../../atlasfiles/horizontalsections/annotation_25_full.nrrd');
Regionsfull=Regionsfull(end:-1:1,1:end,end:-1:1);

annotated = readtable('../analysis/annotation_info_0118_1327_justusedits.csv');
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



atlasCrop=regions;


%% make an outlines file.
outlines_horizontal=zeros(size(atlasCrop));
for i=1:size(atlasCrop,3)
    outlines_horizontal(:,:,i)=edge(squeeze(atlasCrop(:,:,i)),'Sobel',0);
end

outlines_coronal=zeros(size(atlasCrop));
for i=1:size(atlasCrop,1)
    outlines_coronal(i,:,:)=edge(squeeze(atlasCrop(i,:,:)),'Sobel',0);
end
outlines_coronal=permute(outlines_coronal,[3,2,1]);
outlines_coronal=outlines_coronal(end:-1:1,:,end:-1:1);

outlines_sagital=zeros(size(atlasCrop));
for i=1:size(atlasCrop,2)
    outlines_sagital(:,i,:)=edge(squeeze(atlasCrop(:,i,:)),'Sobel',0);
end
outlines_sagital=permute(outlines_sagital,[3,1,2]);
outlines_sagital=outlines_sagital(end:-1:1,end:-1:1,:);

nrrdWriter('../../../atlasfiles/100um_atlas/outlines_25_horizontal_collapsed.nrrd', outlines_horizontal, [25,25,25], [0,0,0], 'raw');
nrrdWriter('../../../atlasfiles/100um_atlas/outlines_25_sagital_collapsed.nrrd', outlines_sagital, [25,25,25], [0,0,0], 'raw');
nrrdWriter('../../../atlasfiles/100um_atlas/outlines_25_coronal_collapsed.nrrd', outlines_coronal, [25,25,25], [0,0,0], 'raw');





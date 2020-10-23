clear all
close all

% apparently Fiji messes up the annotation file if its cropped and then
% saved in there. so do the cropping in matlab instead.
%% make cropped atlases
atlas = nrrdread('../../../atlasfiles/100um_atlas/annotation_100.nrrd');
atlashorizontal=permute(atlas,[2,3,1]);
atlashorizontal=atlashorizontal(132:-1:1,1:114,80:-1:1);
nrrdWriter('../../../atlasfiles/100um_atlas/annotation_100_horizontal.nrrd', atlashorizontal, [100,100,100], [0,0,0], 'raw');

atlascoronal=permute(atlas,[1,3,2]);
atlashorizontal=atlashorizontal(132:-1:1,1:114,80:-1:1);
nrrdWriter('../../../atlasfiles/100um_atlas/annotation_100_coronal.nrrd', atlascoronal, [100,100,100], [0,0,0], 'raw');

% make an outlines file 100um.
outlines_horizontal=zeros(size(atlashorizontal));
for i=1:size(atlashorizontal,3)
    outlines_horizontal(:,:,i)=edge(squeeze(atlashorizontal(:,:,i)),'Sobel',0);
end

outlines_coronal=zeros(size(atlashorizontal));
for i=1:size(atlashorizontal,1)
    outlines_coronal(i,:,:)=edge(squeeze(atlashorizontal(i,:,:)),'Sobel',0);
end
outlines_coronal=permute(outlines_coronal,[3,2,1]);
outlines_coronal=outlines_coronal(end:-1:1,:,end:-1:1);

outlines_sagital=zeros(size(atlashorizontal));
for i=1:size(atlashorizontal,2)
    outlines_sagital(:,i,:)=edge(squeeze(atlashorizontal(:,i,:)),'Sobel',0);
end
outlines_sagital=permute(outlines_sagital,[3,1,2]);
outlines_sagital=outlines_sagital(end:-1:1,end:-1:1,:);

nrrdWriter('../../../atlasfiles/100um_atlas/outlines_100_horizontal.nrrd', outlines_horizontal, [100,100,100], [0,0,0], 'raw');
nrrdWriter('../../../atlasfiles/100um_atlas/outlines_100_sagital.nrrd', outlines_sagital, [100,100,100], [0,0,0], 'raw');
nrrdWriter('../../../atlasfiles/100um_atlas/outlines_100_coronal.nrrd', outlines_coronal, [100,100,100], [0,0,0], 'raw');





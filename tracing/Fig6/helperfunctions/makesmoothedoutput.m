function[smoothed_scaled]=makesmoothedoutput(P,cutoff,output,mask)
overlay=zeros(size(P(1).rawdatafull));
for i=1:numel(P)
    if isempty(mask)==0
        binarized=(P(i).rawdatafull.*uint16(~mask))>cutoff;
    else
    binarized=P(i).rawdatafull>cutoff;
    end
    
    pixelnum=sum(sum(sum(binarized)));
%     norm=binarized;
    norm=binarized./pixelnum;

    overlay=overlay+norm;
end

smoothed=smooth3(overlay,'box',9);
smoothed_scaled=uint16(smoothed./max(max(max(smoothed)))*2^16);

for i=1:size(overlay,3)
    imwrite(squeeze(smoothed_scaled(:,:,i)),[output,'/smoothed_',int2str(i),'.tif']);
end
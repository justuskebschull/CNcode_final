function [prob,smoothedprob]=calculateprojectionprobability(P,cutoff,output,mask)
mkdir(output);
overlay=zeros(size(P(1).rawdatafull));
for i=1:numel(P)
    if isempty(mask)==0
        binarized=(P(i).rawdatafull.*uint16(~mask))>cutoff;
    else
        binarized=P(i).rawdatafull>cutoff;
    end
    norm=binarized;
    
%     pixelnum=sum(sum(sum(binarized)));
%     norm=binarized./pixelnum;

    overlay=overlay+norm;
end

prob=overlay./numel(P);
smoothedprob=smooth3(prob,'box',9);
smoothed_scaled=uint16(smoothedprob.*2^16);

for i=1:size(overlay,3)
    imwrite(squeeze(smoothed_scaled(:,:,i)),[output,'/smoothed_',int2str(i),'.tif']);
end
function []=plotsortedheatmap(data_l,data_r,dcn,top)
DCN={'FN','IN','DN'};
for i=1:3
    %left hemi
    dl=data_l(:,dcn==i);
    [~,ix]=sort(mean(dl,2),'descend');
    dl=dl(ix,:);
    if isempty(top)==0
        dl=dl(1:top,:);
    else
%         t=ttest(dl',0,'Alpha',0.05/size(data_l,1));
        t=ttest(dl',0,'Alpha',0.05);
        t(isnan(t))=0;
        dl=dl(logical(t),:);
    end
    
    figure;imagesc(dl);
    ax=gca;
    ax.CLim=[0 5e-7];
    title([DCN{i},' projections to left hemisphere']);
    
    %right hemi
    dr=data_r(:,dcn==i);
    [~,ix]=sort(mean(dr,2),'descend');
    dr=dr(ix,:);
    if isempty(top)==0
        dr=dr(1:top,:);
    else
%         t=ttest(dr',0,'Alpha',0.05/size(data_r,1));
        t=ttest(dr',0,'Alpha',0.05);
        t(isnan(t))=(0);
        dr=dr(logical(t),:);
    end
    figure;imagesc(dr);
    ax=gca;
    ax.CLim=[0 5e-7];
    title([DCN{i},' projections to right hemisphere']);
end
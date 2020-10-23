function [pvalue,names,namesFs,sortedmean,sorteddata]=regionttests(data_r,DCN_r,annotated_filtered)

p=[];
for i=1:size(data_r,1)
    [~,p(i,:)]=ttest2(data_r(i,DCN_r==1),data_r(i,DCN_r==2));
end

[pvalue,idx1]=sort(p);
names=annotated_filtered.name(idx1);

%make heatmap for each comparison, sorting areas by mean expression
%difference.

%calculate means
m=zeros(size(data_r,1),numel(unique(DCN_r)));

for i=1:numel(unique(DCN_r))
    m(:,i)=mean(data_r(:,DCN_r==i),2);
end
sortedmean=m(idx1,:);
sorteddata=data_r(idx1,:);

d=m(idx1,1)-m(idx1,2);
d=d(pvalue<0.05);
[~,idxd]=sort(d);
dat=data_r(idx1,DCN_r==1 | DCN_r==2);
dat=dat(pvalue<0.05,:);
namesF=names(pvalue<0.05);
namesFs=namesF(idxd);
figure;imagesc(dat(idxd,:));title('ZIvsRet')
set(gca, 'YTick', 1:numel(namesFs), 'YTickLabel', namesFs);

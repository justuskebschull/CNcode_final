function [FNvsIN,FNvsDN,INvsDN]=region_anova(data_r,DCN_r,annotated_filtered)

p_multcomp=[];
for i=1:size(data_r,1)    
[p,tbl,stats]=anova1(data_r(i,:),DCN_r,'off');
[c,~,~,gnames] = multcompare(stats,[],'off');
p_multcomp(:,i)=c(:,6);
end

[FNvsIN.pvalue,idx1]=sort(p_multcomp(1,:));
FNvsIN.names=annotated_filtered.name(idx1);


[FNvsDN.pvalue,idx2]=sort(p_multcomp(2,:));
FNvsDN.names=annotated_filtered.name(idx2);


[INvsDN.pvalue,idx3]=sort(p_multcomp(3,:));
INvsDN.names=annotated_filtered.name(idx3);

%make heatmap for each comparison, sorting areas by mena expression
%difference.

%calculate means
m=zeros(size(data_r,1),numel(unique(DCN_r)));

for i=1:numel(unique(DCN_r))
    m(:,i)=mean(data_r(:,DCN_r==i),2);
end

%FNvsIN
d=m(idx1,1)-m(idx1,2);
d=d(FNvsIN.pvalue<0.05);
[~,idxd]=sort(d);
dat=data_r(idx1,DCN_r==1 | DCN_r==2);
dat=dat(FNvsIN.pvalue<0.05,:);
figure;imagesc(dat(idxd,:));title('FNvsIN')


%FNvsDN
d=m(idx2,1)-m(idx2,3);
d=d(FNvsDN.pvalue<0.05);
[~,idxd]=sort(d);
dat=data_r(idx2,DCN_r==1 | DCN_r==3);
dat=dat(FNvsDN.pvalue<0.05,:);
figure;imagesc(dat(idxd,:));title('FNvsDN')


%INvsDN
d=m(idx3,2)-m(idx3,3);
d=d(INvsDN.pvalue<0.05);
[~,idxd]=sort(d);
dat=data_r(idx3,DCN_r==2 | DCN_r==3);
dat=dat(INvsDN.pvalue<0.05,:);
figure;imagesc(dat(idxd,:));title('INvsDN')


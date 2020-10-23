function [p,p_multcomp,sorteddata,sortednames,sortedp]=region_anova(data_r,DCN_r,annotated_filtered,pcutoff,T)
p=[];
p_multcomp=[];
for i=1:size(data_r,1)    
[p(i),tbl,stats]=anova1(data_r(i,:),DCN_r,'off');
[c,~,~,gnames] = multcompare(stats,[],'off');
p_multcomp(:,i)=c(:,6);
end



%calculate mean innervations per area and injection site.

m=zeros(size(data_r,1),numel(unique(DCN_r)));

for i=1:numel(unique(DCN_r))
    m(:,i)=mean(data_r(:,DCN_r==i),2);
end

[m_max,m_max_i]=max(m,[],2,'omitnan');

m_max_f=m_max(p<pcutoff);
m_max_if=m_max_i(p<pcutoff);
anno_ff=annotated_filtered(p<pcutoff,:);
data_f=data_r(p<pcutoff,:);
p_f=p(p<pcutoff);

sortednames=[];
sorteddata=[];
sortedp=[];
for i=1:numel(unique(DCN_r))
    tmp_names=anno_ff.name(m_max_if==i);
    tmp_data=data_f(m_max_if==i,:);
    tmp_mean=m_max_f(m_max_if==i);
    tmp_p=p_f(m_max_if==i);
    
    [~,idx]=sort(tmp_mean,'descend');
    
    sortednames=[sortednames;tmp_names(idx)];
    sorteddata=[sorteddata;tmp_data(idx,:)];
    sortedp=[sortedp,tmp_p(idx)];
end

%plot figure
figure;imagesc(sorteddata)

ax = gca;

% Set where ticks will be
ax.YTick = 1:numel(sortednames);
% Set TickLabels;
ax.YTickLabel = sortednames;


title(T)



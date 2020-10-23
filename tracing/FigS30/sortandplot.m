function []=sortandplot(data,DCN,names,top2)

% calculate mean
DCNlabel={'FN','IN','DN'};
for i=1:3
    [mean_projections{i},idx]=sort(median(data(:,DCN==i),2),'descend');
    names_projections{i}=names(idx);
    std_projections_tmp=std(data(:,DCN==i),[],2);
    std_projections{i}=std_projections_tmp(idx);
    
    if isempty(top2)~=1
        top=top2;
    else
        top=find(mean_projections{i}==0,1,'first')-1;
    end
    
    figure;
    bar((mean_projections{i}(1:top)));
    xticks(1:top);
    xticklabels(names_projections{i}(1:top));
    xtickangle(45);
    hold on;
    er=errorbar(1:top,mean_projections{i}(1:top),std_projections{i}(1:top)./sqrt(sum(DCN==i)));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    title([DCNlabel{i},' top projections']);
    hold off;

    %make heatmap
    figure;
    data_s=data(idx,:);
    imagesc(data_s);
    title(['sorted by',DCNlabel{i}]);





end


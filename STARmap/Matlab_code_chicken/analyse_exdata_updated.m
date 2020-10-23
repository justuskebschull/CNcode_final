%% 
clear all
close all

%% load data
load('data.mat');

%% drop empty channels
for s=1:7
    data(s).outputmatrix=data(s).outputmatrix(sum(data(s).outputmatrix,2)~=0,:);
end

%% display the number of discovered signals per cell
% s=1;
% tmp=zeros(size(data(s).cells));
% for i=1:size(data(s).outputmatrix,2)-1
%     [row,col] = find(data(s).cells==i);
%     tmp(row,col)=data(s).outputmatrix(3,i+1);
% end
% figure;
% imagesc(tmp);


%%
s=5;
excells=ismember(data(s).cells,find(data(s).outputmatrix(2,2:end)>3));
figure;imshow(excells)


%% precalculate binarized expression images
genes={"SLC17A6","TAC1","KCNAB1","ADCYAP1","GAD1","RTN4R","PKIB","XDH"};
    
    
% updated thresholds    
thresholds=5.*ones(length(genes),1);
thresholds(1)=10;
thresholds(2)=3;
thresholds(5)=7;
thresholds(7)=7;
thresholds(8)=3; %only for chicken2, otherwise 5.
thresholds(3)=10;


for s=1:7
    for i=1:numel(genes)
        data(s).expression{i}=ismember(data(s).cells,find(data(s).outputmatrix(i,2:end)>thresholds(i)));
    end
end


% %% get nissl images
% L=[];
% sections=dir('nisslsegmentation\raw\*.tif');
% for s=1:5
%     L(s).nissl=imread(['nisslsegmentation\raw\',sections(s).name]);
% end


%% make custom colormap

% custommap=[cmap([1 0 0],6,10,10);cmap([0 1 0],6,10,10);cmap([1 1 0],3,10,10);repmat([0.4 0.4 0.4],20,1);[0 0 1]];
% map2=[[1 0 0];[0 1 0];[1 0 1];[0 0 1]; [1,0,1];[0 1 1];repmat([0.4 0.4 0.4],50,1);[0.8, 0.8, 0.8];[1,0.55,0]];
map2=[[1 0 0];[0 1 0];[0 0 1];[1 0 1]; [1,1,0];[0 1 1];repmat([0.4 0.4 0.4],50,1);[0.8, 0.8, 0.8];[1,0.55,0]];

%% find subnuclei.

for s=1:7

expression=data(s).expression;

Med=expression{1}.*expression{3}.*(1-expression{4});
MedL=expression{1}.*expression{2};%.*(1-expression{7});
IntX=expression{1}.*expression{4}; 
%IntP=expression{7}.*(expression{6} | expression{4}).*(1-expression{2});
IntP=expression{1}.*(expression{7}).*(1-expression{3}).*(1-expression{4});

m=makemask(cat(3,Med,MedL,IntX,IntP),0);
L(s).subnuclei=m; 
end


clear Med MedL IntX IntP
%% quick subnuclei from saved
%subnucleipal=[hex2rgb('#EA9292');hex2rgb('#FF7F00');hex2rgb('#33a02');hex2rgb('#1F78B4');repmat([0.95 0.95 0.95],60,1)];
%subnucleipal=[hex2rgb('#EA9292');hex2rgb('#ED1C24');hex2rgb('#33a02');hex2rgb('#1F78B4');repmat([0.95 0.95 0.95],60,1)];

subnucleipal=[hex2rgb('#EA9292');hex2rgb('#F7F008');hex2rgb('#33a02');hex2rgb('#1F78B4');repmat([0.95 0.95 0.95],60,1)];
for s=1:7
    L(s).ex=data(s).expression{1};
    
    boundaries = bwboundaries(L(s).ex);
    
    n=L(s).subnuclei;
    n(L(s).ex>0 & n==0)=57;
    figure;
    h=imshow(label2rgb(n,subnucleipal,'white'));
    title(['section ',int2str(s)]);
    ax = gca;
    if(s>4)
    ax.YDir = 'normal'
    end
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'black');
    end
   % pause();
    print(gcf,['figure_exports/subnuclei',int2str(s)],'-dpdf','-fillpage','-r1000')
    close gcf
end


%% look at Class-A vs Class-B.

for s=1:7

expression=data(s).expression;

ClassA=expression{1}.*(1-expression{8});
ClassB=expression{1}.*expression{8};

m=makemask(cat(3,ClassA,ClassB),0);
L(s).classes=m;    
end

%% quick Classes from saved
Classpal=[hex2rgb('#ff6347');hex2rgb('#00688b');repmat([0.95 0.95 0.95],60,1)];

for s=1:7
    L(s).ex=data(s).expression{1};
    
    boundaries = bwboundaries(L(s).ex);
    
    n=L(s).classes;
    n(L(s).ex>0 & n==0)=57;
    figure;
    h=imshow(label2rgb(n,Classpal,'white'));
    title(['section ',int2str(s)]);
    ax = gca;
    if s>4
    ax.YDir = 'normal'
    end
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'black');
    end
      print(gcf,['figure_exports/classes_ch2',int2str(s)],'-dpdf','-fillpage','-r1000')
    close gcf
end


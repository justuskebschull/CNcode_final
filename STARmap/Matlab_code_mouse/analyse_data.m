%% 
clear all
close all

%% load data
load('data.mat');
% load('datavirus.mat');
load('ROIS3.mat');

%% precalculate binarized expression images
genesold={'slc17a6','sorcs3','trpm3','acan','sv2c','fign','ankfn1','cd24a',...
    'cpne4','man1a','stac2','nxph1','rik','penk','pcdh11x','zfp8004b',...
    'calb2','cacna1e','asic2','spp1'};
genes={'slc17a6','slc24a4','gad1','slc6a5',...
    'acan','dscaml1','sv2c','slc6a1', ...
    'ankfn1','penk','stac2','sorcs2',...
    'pcdh11x','epha5','calb2','kitl',...
    'cdh22','nckap5','sez6','flp0',...
    'cre','rh2','rh1','rh3',...
    'rh5','rh10','rh6','rh13'};
    
    
% updated thresholds    
thresholds=2.*ones(length(genes),1);
% thresholds(5)=4; %acan
thresholds(5)=6; %acan


thresholds(7)=4;
thresholds(11)=2; %stac
thresholds(9)=0; %ankfn1
thresholds(16)=0; %kitl
thresholds(19)=2; %sez6
threshold(6)=1; %dscaml1
threshold(17)=0; %csh22
threshold(10)=10; %penk
threshold(1)=2; %slc17a6

threshold(4)=0;%slc6a5
threshold(3)=4;%gad1

% for s=1:numel(exdata)
for s=[1:2:11]
    for i=1:numel(genes)
        exdata(s).expression{i}=ismember(exdata(s).cells,find(exdata(s).outputmatrix(i,2:end)>thresholds(i)));
    end
end

%% get nissl images
L=[];
sections=dir('nisslsegmentation\raw\*.tif');
for s=1:numel(exdata)
    L(s).nissl=imread(['nisslsegmentation\raw\',sections(s).name]);
end

%% memory management
% look only at left (contra) side of brain
L=[];

L.nissl=[];
for s=[2:2:12]
    exdata(s).outputmatrix=[];
    exdata(s).genelabels=[];
    exdata(s).cells=[];
    exROIs(s).FN=[];
    exROIs(s).IN=[];
    exROIs(s).DN=[];
end


%% make custom colormap

custommap=[cmap([1 0 0],6,10,10);cmap([0 1 0],6,10,10);cmap([1 1 0],3,10,10);repmat([0.4 0.4 0.4],20,1);[0 0 1]];
map2=[[1 0 0];[0 1 0];[0 0 1];[1 1 0]; [1,0,1];[0 1 1];repmat([0.4 0.4 0.4],50,1);[1,0.55,0]];
RLpal=readtable('RLpal.csv');
RLpal_sorted=RLpal([3 2 1 8 10 9 11 13 12 6 4 5 7 14 15],2);
RLpal_sorted_rgb=hex2rgb(table2array(RLpal_sorted));
RLpal_sorted_rgb_padded=[RLpal_sorted_rgb;repmat([0.8 0.8 0.8],60,1)];
VZpal=readtable('VZpal.csv');
VZpal_sorted=VZpal([1,3,4,5],2);
% VZpal_sorted=VZpal([5,1,3,4],2);
VZpal_sorted_rgb=hex2rgb(table2array(VZpal_sorted));
VZpal_sorted_rgb2=[VZpal_sorted_rgb;RLpal_sorted_rgb(6,:)];
VZpal_sorted_rgb_padded=[VZpal_sorted_rgb2;repmat([0.8 0.8 0.8],60,1)];

%% look at DN

for s=[1:2:11]
    expression=exdata(s).expression;

    class2=expression{1}.*(expression{5} & expression{7}); %slc17a6+acan or Sv2c
    class1=expression{1}.*(1-(expression{5} & expression{7}));%.*expression{8}; %slc17a6-(acan&sv2c)+slc6a1 
    ex7=exROIs(s).DN.*class1;
    ex0=exROIs(s).DN.*class2.*(1-expression{10});%class2 -penk
    ex14=exROIs(s).DN.*class2.*expression{10}; %class2 +penk
    
    m=makemask(cat(3,ex7,ex0,ex14),12);
    L(s).DN=m;
end
clear ex0 ex7 ex14

%% look at FN
close all;
for s=[1:2:11]
    expression=exdata(s).expression;
    exdata(s).expression{4}=ismember(exdata(s).cells,find(exdata(s).outputmatrix(4,2:end)>10));
    
    
    class2=expression{1}.*(expression{5} & expression{7}); %slc17a6+acan and/or Sv2c
    class1=(expression{1}.*(1-(expression{5} & expression{7}))); 

            
            
    ex6=exROIs(s).FN.*class1.*expression{15}; %slc17a6+calb2
    ex2=exROIs(s).FN.*class1.*(1-expression{15}).*(1-expression{10}); %class1-calb2-penk
    ex15=exROIs(s).FN.*class1.*expression{10}; %class1+penk
            
    ex9=exROIs(s).FN.*class2.*expression{10}; %class2+penk
    ex11=exROIs(s).FN.*class2.*(1-expression{10});%class2-penk

    gly4=(exROIs(s).FN.*expression{3}.*expression{4}.*(expression{5}).*(1-expression{19})) | ...
        (exROIs(s).FN.*(1-expression{3}).*expression{4}.*expression{5});        


            m=makemask(cat(3,ex6,ex2,ex15,ex9,ex11,gly4),0);
            L(s).FN=m;

    exdata(s).expression{4}=ismember(exdata(s).cells,find(exdata(s).outputmatrix(4,2:end)>thresholds(4)));
end
clear ex6 ex2 ex15 ex9 ex11 gly4 

%% look at IN
% close all;

for s=[1:2:11]
    expression=exdata(s).expression;
    class2=expression{1}.*(expression{5} & expression{7}); %slc17a6+acan or Sv2c
    class1=(expression{1}.*(1-(expression{5} & expression{7}))); 
     
    ex8=exROIs(s).IN.*class2.*expression{9};%.*expression{11}; %class2+Ankfn1+stac2
    
    ex1=exROIs(s).IN.*class2.*(1-expression{9}).*expression{11}.*(1-expression{6}).*expression{16}; %class2-ankfn1+Stac2-Dscaml1
    
    ex4=exROIs(s).IN.*class2.*(1-expression{9}).*(1-expression{11}).*(1-expression{19}); %class2-ankfn1-Stac2

    ex10=exROIs(s).IN.*class1.*(1-expression{9}).*(1-expression{8});%.*expression{17};%class1-Ankfn1
  
    ex12=exROIs(s).IN.*class1.*expression{9}.*expression{16};%class1+Ankfn1+kitl
    ex5=exROIs(s).IN.*class1.*expression{9}.*(1-expression{16});%class1+Ankfn1-kitl
    
    clear class1 class2 expression
    m=makemask(cat(3,ex8,ex1,ex4,ex10,ex12,ex5),6);
    clear ex8 ex1 ex4 ex10 ex12 ex5
    L(s).IN=m;
    clear m


end
clear ex8 ex1 ex4 ex10 ex12 ex5
%% what am I missing?
for s=[1:2:11]
    expression=exdata(s).expression;
  L(s).all=(exROIs(s).FN | exROIs(s).IN | exROIs(s).DN).*expression{1};
end

        



%% now, let's look at inhibitory cells.
genesinhold={'Gad1',"Ptprz1",'Slc6a5','Sez6',...
    'Chsy3','Nckap5','Cdh22','Grik1',...
    'Spon1','Slc24a4','Kirrel3','Zfhx4'};

for s=[1:2:11]
    expression=exdata(s).expression;
    expression{19}=ismember(exdata(s).cells,find(exdata(s).outputmatrix(19,2:end)>2));
    
    inhROIs=(exROIs(s).FN | exROIs(s).IN | exROIs(s).DN).*(1-expression{1});
    
    inh1=inhROIs.*expression{3}.*(1-expression{4}).*(1-expression{17}); %gad1-slc6a5 no slc24a4 so far..
    gly2=inhROIs.*expression{3}.*(1-expression{4}).*expression{17};%gad1-slc6a5+cdh22
    
    gly0=inhROIs.*expression{3}.*expression{4}.*(1-expression{19}) .* (1- expression{5});
    gly3=inhROIs.*expression{3}.*expression{4}.*(expression{19}).*(1-expression{5});
    gly4=(inhROIs.*expression{3}.*expression{4}.*(expression{5}).*(1-expression{19})) | ...
        (inhROIs.*(1-expression{3}).*expression{4});
 
    m=makemask(cat(3,inh1,gly0,gly2,gly3,gly4),15);
    L(s).inh=m;

end


%% quick inh from saved
close all;

for s=[1:2:11]
% for s=6
    expression=exdata(s).expression;
  L(s).allinh=(exROIs(s).FN | exROIs(s).IN | exROIs(s).DN).*(expression{3} | expression{4}).*(1-expression{1});



    n=L(s).inh;
    n(n>0)=n(n>0)-15;
    n(L(s).allinh>0 & n==0)=57;

    figure;
    h=imshow(label2rgb(n,VZpal_sorted_rgb_padded,'white'));
    boundaries = bwboundaries(n);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'black');
    end   
    hold off;    title(['section ',int2str(s)]);
    ax = gca;
    ax.YDir = 'normal'
end


%% plot one CN at a time, splitting class A and class B into individual plots with same colors for same subarea.
% this time around, have non-plotted class very faint.
subareapal=readtable('subareapal.csv');
subareapal1=table2cell(subareapal(:,2));
subareapal1{8}='#BF93E2';
subareapal1{9}='#B8D3F4';
subareapal1{10}='#961528';

subareapal_sorted=subareapal1([2,10,1,1,3,2,5,4,8,4,5,9,6,6,7]);
subareapal_sorted_rgb=hex2rgb(subareapal_sorted);
subareapal_sorted_rgb=[subareapal_sorted_rgb;repmat([0.95 0.95 0.95],60,1)];

classA=[1,2,3,10,11,12,13];
classB=[4,5,6 7,8,9,14,15];

subareapal_B=subareapal_sorted_rgb;
for i=classA
subareapal_B(i,:)=[0.9 0.9 0.9];
end

subareapal_A=subareapal_sorted_rgb;
for i=classB
subareapal_A(i,:)=[0.9 0.9 0.9];
end

%%plot class A
for s=[1:2:11]
    figure;
    %class A
    m=L(s).FN+(L(s).IN)+L(s).DN;
    m(L(s).all>0 & m==0)=36;
    cA=m;
    cB=m;
    cA(ismember(m,classB))=0;
    cB(ismember(m,classA))=0;

    h=imshow(label2rgb(m,subareapal_A,'white'));
    boundaries = bwboundaries(cA);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'black');
    end
  
    boundaries = bwboundaries(cB);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'Color',[0.7 0.7 0.7]);
    end
 
    hold off;
    title(['A section',int2str(s)])

    ax = gca;
    ax.YDir = 'normal';
    print(gcf,['figure_exports_fixed9/A_updated_sectioncontra',int2str(s)],'-dpdf','-fillpage','-r1000')
    close gcf
end

for s=[1:2:11]
    figure;
    m=L(s).FN+(L(s).IN)+L(s).DN;
    m(L(s).all>0 & m==0)=36;
    cA=m;
    cB=m;
    cA(ismember(m,classB))=0;
    cB(ismember(m,classA))=0;

    h=imshow(label2rgb(m,subareapal_B,'white'));
    boundaries = bwboundaries(cB);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'black');
    end
  
    boundaries = bwboundaries(cA);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'Color',[0.7 0.7 0.7]);
    end
 
    hold off;
    title(['B section',int2str(s)])

%     set(h,'AlphaData',(m>0).*0.5);
    ax = gca;
    ax.YDir = 'normal';
    print(gcf,['figure_exports_fixed9/B_updated_sectioncontra',int2str(s)],'-dpdf','-fillpage','-r1000')
    close gcf
end

classI=[16,17,18,19];
for s=[1:2:11]
    figure;
    m=L(s).inh;
    m(L(s).allinh>0 & m==0)=36;
    cA=m;
    cB=m;
    cA(ismember(m,[classI,20]))=0;
    cB(ismember(m,[20]))=0;
    
    h=imshow(label2rgb(cB,[repmat([0 0 0],15,1);VZpal_sorted_rgb_padded],'white'));
    boundaries = bwboundaries(cB);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'black');
    end
  
    boundaries = bwboundaries(cA);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'Color',[0.7 0.7 0.7]);
    end
 
    hold off;
    title(['B section',int2str(s)])

%     set(h,'AlphaData',(m>0).*0.5);
    ax = gca;
    ax.YDir = 'normal';
    print(gcf,['figure_exports_fixed9/inh_updated_section_contra',int2str(s)],'-dpdf','-fillpage','-r1000')
    close gcf
end


%% plot overlap of class A and class b
Classpal=[hex2rgb('#ff6347');hex2rgb('#00688b');repmat([0.95 0.95 0.95],60,1)];
classA=[1,2,3,10,11,12,13];
classB=[4,5,6 7,8,9,14,15];

%%plot class A
for s=[1:2:11]
    figure;
    %class A
    m=L(s).FN+(L(s).IN)+L(s).DN;
    m(L(s).all>0 & m==0)=36;
    cA=m;
    cB=m;

    cA(ismember(m,classB))=2;
    cA(ismember(m,classA))=1;
    
    cB(ismember(m,[classB,classA]))=0;
    
    h=imshow(label2rgb(cA,Classpal,'white'));
    boundaries = bwboundaries(cA);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'black');
    end
  
    boundaries = bwboundaries(cB);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'Color',[0.7 0.7 0.7]);
    end
 
    hold off;
    title(['A section',int2str(s)])

    ax = gca;
    ax.YDir = 'normal';
    print(gcf,['figure_exports_fixed9/overlay',int2str(s)],'-dpdf','-fillpage','-r1000')
    close gcf
end

%% plot all clusters of class A and class b
classA=[1,2,3,10,11,12,13];
classB=[4,5,6 7,8,9,14,15];
inh=[16,17,18,19];
clusterpal=[RLpal_sorted_rgb;VZpal_sorted_rgb;repmat([0.95 0.95 0.95],60,1)];
%%plot class A
for s=[1:2:11]
    figure;
    %class A
    L(s).inh(L(s).inh==20)=0;
    
    m=L(s).FN+(L(s).IN)+L(s).DN+L(s).inh;
    m((L(s).FN>0+L(s).IN>0+L(s).DN>0)>1)=36;
    m((L(s).all>0 | L(s).allinh>0) & m==0)=36;
    cA=m;
    cB=m;
    
    cB(ismember(m,[classB,classA,inh]))=0;
    
    h=imshow(label2rgb(cA,clusterpal,'white'));
    boundaries = bwboundaries(cA);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'black');
    end
  
    boundaries = bwboundaries(cB);
    hold on;
    for k=1:numel(boundaries)
        b=boundaries{k};
        plot(b(:,2),b(:,1),'Color',[0.7 0.7 0.7]);
    end
 
    hold off;
    title(['A section',int2str(s)])

    ax = gca;
    ax.YDir = 'normal';
    print(gcf,['figure_exports_fixed9/allclusters',int2str(s)],'-dpdf','-fillpage','-r1000')
    close gcf
end



%% can I make a nice heatmap of gene expression?
for c=1:19
mat(c).data=[];
end
for s=[1:2:11]
        L(s).inh(L(s).inh==20)=0;

    m=L(s).FN+(L(s).IN)+L(s).DN+L(s).inh;
    for c=1:19
        cells=unique(exdata(s).cells(m==c));
        mat(c).data=[mat(c).data,exdata(s).outputmatrix(:,cells+1)];
    end
end
bgmat=[];
for c=1:19
   bgmat=[bgmat,mat(c).data]; 
end
%%
clustergram(bgmat([1,3,4,5,7,8,9,10,11,15,16,19],:)','Standardize','none','Cluster','column','Colormap','parula','Symmetric',0)


figure;imagesc(zscore(bgmat([1,3,4,5,7,8,9,10,11,15,16,19],:))');colorbar;

figure;imagesc((bgmat([1,3,4,5,7,8,9,10,11,15,16,19],:))');
colorbar;caxis([0 4])

%%
meanexpression=[];
sd=[];
CIu=[];
CIl=[];
SEM=[];



for c=1:19
    meanexpression(:,c)=mean(mat(c).data,2);
    sd(:,c)=std(mat(c).data,0,2);
    SEM(:,c) = std(mat(c).data,0,2)./sqrt(size(mat(c).data,2));
end

genenames_retro={'slc17a6','slc24a4','gad1','slc6a5',...
    'acan','dscaml1','sv2c','slc6a1', ...
    'ankfn1','penk','stac2','sorcs2',...
    'pcdh11x','epha5','calb2','kitl',...
    'cdh22','nckap5','sez6','flp0',...
    'cre','rh2','rh1','rh3',...
    'rh5','rh10','rh6','rh13'};

colors=[RLpal_sorted_rgb;VZpal_sorted_rgb2];
figure;
k=1;
genes=[1,3,4,5,7,8,9,10,11,15,16,19];
ordering=[3 2 1 11 12 10 13 4 6 5 7 9 8 14 15 16:19];
colors_s=colors(ordering,:);

for i=genes
    subplot(numel(genes),1,k)
    b=bar(1:19,meanexpression(i,ordering),'FaceColor','flat');
    for z = 1:19
    b.CData(z,:) = colors_s(z,:);
    end
    
    set(gca,'xtick',[])

    k=k+1;    
    hold on
%     er = errorbar(1:19,meanexpression(i,:),CIl(i,:),CIu(i,:));    

    er = errorbar(1:19,meanexpression(i,ordering),SEM(i,ordering));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';      
    hold off
        title(genenames_retro{i})

end
    print(gcf,'figure_exports_fixed9/exbargraph','-dpdf','-r1000','-fillpage')

 

%% load cell type classification, define subnuclei, and quantify number of each cell type per subnucleus

%%
load L_contraplotting_fixed_snapshot.mat


subnuclei_ROIs=[];

%% set ROIS
% 
% % this time around, have non-plotted class very faint.
% subareapal=readtable('subareapal.csv');
% subareapal1=table2cell(subareapal(:,2));
% subareapal1{8}='#BF93E2';
% subareapal1{9}='#B8D3F4';
% subareapal1{10}='#961528';
% 
% subareapal_sorted=subareapal1([2,10,1,1,3,2,5,4,8,4,5,9,6,6,7]);
% subareapal_sorted_rgb=hex2rgb(subareapal_sorted);
% subareapal_sorted_rgb=[subareapal_sorted_rgb;repmat([0.95 0.95 0.95],60,1)];
% 
% classA=[1,2,3,10,11,12,13];
% classB=[4,5,6 7,8,9,14,15];
% 
% subareapal_B=subareapal_sorted_rgb;
% for i=classA
% subareapal_B(i,:)=[0.9 0.9 0.9];
% end
% 
% subareapal_A=subareapal_sorted_rgb;
% for i=classB
% subareapal_A(i,:)=[0.9 0.9 0.9];
% end
% 
% for s=[1:2:11]
%     figure;
%     m=L(s).FN+(L(s).IN)+L(s).DN;
%     m(L(s).all>0 & m==0)=36;
%     cA=m;
%     cB=m;
%     cA(ismember(m,classB))=0;
%     cB(ismember(m,classA))=0;
% 
%     h=imshow(label2rgb(m,subareapal_B,'white'));
%     boundaries = bwboundaries(cB);
%     hold on;
%     for k=1:numel(boundaries)
%         b=boundaries{k};
%         plot(b(:,2),b(:,1),'black');
%     end
%   
%     boundaries = bwboundaries(cA);
%     hold on;
%     for k=1:numel(boundaries)
%         b=boundaries{k};
%         plot(b(:,2),b(:,1),'Color',[0.7 0.7 0.7]);
%     end
%  
%     hold off;
%     title(['B section',int2str(s)])
% 
% %     set(h,'AlphaData',(m>0).*0.5);
%     ax = gca;
%     ax.YDir = 'normal';
%     [int2str(s),'Lat']
%     subnucleiROIs(s).Lat=roipoly;
%     [int2str(s),'IntA']
%     subnucleiROIs(s).IntA=roipoly;
%     [int2str(s),'IntP']
%     subnucleiROIs(s).IntP=roipoly;
%     [int2str(s),'MedL']
%     subnucleiROIs(s).MedL=roipoly;
%     [int2str(s),'Med']
%     subnucleiROIs(s).Med=roipoly;
%     [int2str(s),'MedDL']
%     subnucleiROIs(s).MedDL=roipoly;
%     close all;
% end
% 
% % save('contra_subnuclei_ROIS1.mat','subnucleiROIs')

%% quantify number of cells per section per subnucleus according to above defined Rois
classA=[1,2,3,10,11,12,13];
classB=[4,5,6 7,8,9,14,15];
inh=[16,17,18,19];

quantification=zeros(12,20,6);
ROIarea=zeros(12,6);

for s=[1:2:11]
    L(s).inh(L(s).inh==20)=0;
    
    m=L(s).FN+(L(s).IN)+L(s).DN+L(s).inh;
    
    m((L(s).FN>0+L(s).IN>0+L(s).DN>0)>1)=36;
    m((L(s).all>0 | L(s).allinh>0) & m==0)=36;
    
    %record area taken up by each ROI
   
    area(s,1)=sum(sum(subnucleiROIs(s).Lat));
    area(s,2)=sum(sum(subnucleiROIs(s).IntA));
    area(s,3)=sum(sum(subnucleiROIs(s).IntP));
    area(s,4)=sum(sum(subnucleiROIs(s).MedL));
    area(s,5)=sum(sum(subnucleiROIs(s).Med));
    area(s,6)=sum(sum(subnucleiROIs(s).MedDL));
    
    for C=1:20 %loop through cell types
        tmp=bwconncomp(subnucleiROIs(s).Lat.*m==C);
        quantification(s,C,1)=tmp.NumObjects;
        tmp=bwconncomp(subnucleiROIs(s).IntA.*m==C);
        quantification(s,C,2)=tmp.NumObjects;
        tmp=bwconncomp(subnucleiROIs(s).IntP.*m==C);
        quantification(s,C,3)=tmp.NumObjects;
        tmp=bwconncomp(subnucleiROIs(s).MedL.*m==C);
        quantification(s,C,4)=tmp.NumObjects;
        tmp=bwconncomp(subnucleiROIs(s).Med.*m==C);
        quantification(s,C,5)=tmp.NumObjects;
        tmp=bwconncomp(subnucleiROIs(s).MedDL.*m==C);
        quantification(s,C,6)=tmp.NumObjects;
    end
end


%%

load('contra_quantification.mat','quantification','area');

celltypequant=squeeze(sum(quantification,1));
%zero out cell recognitions outsid their home nuclei;
celltypequant(1:12,1)=0;
celltypequant(1:6,2)=0;celltypequant(13:15,2)=0;
celltypequant(1:6,3)=0;celltypequant(13:15,3)=0;
celltypequant(7:15,4)=0;
celltypequant(7:15,5)=0;
celltypequant(7:15,6)=0;



celltypequant_fraction=celltypequant./repmat(sum(celltypequant,2),1,size(celltypequant,2));
tmp=celltypequant./repmat(sum(area,1),size(celltypequant,1),1);
celltypequant_fraction_norm=tmp./repmat(sum(tmp,2),1,size(tmp,2));


celltypequant=celltypequant(:,[5,4,6,2,3,1]);
celltypequant_fraction=celltypequant_fraction(:,[5,4,6,2,3,1]);
celltypequant_fraction_norm=celltypequant_fraction_norm(:,[5,4,6,2,3,1]);


% get subarea pal
subareapal=readtable('subareapal.csv');
subareapal1=table2cell(subareapal(:,2));
subareapal_sorted_rgb=hex2rgb(subareapal1);



figure;
imagesc(celltypequant_fraction)
title('Retro1_contra_quant')
print(gcf,['figure_exports_fixed9/Retro1_contra_quant_heatmap'],'-dpdf','-r1000')
close gcf


ordering=[3 2 1 11 12 10 13 4 6 5 7 9 8 14 15 16:19];
figure;
b=bar(celltypequant_fraction(ordering,:),'FaceColor','flat');
for z = 1:6
    b(z).FaceColor = subareapal_sorted_rgb(z,:);
end
title('Retro1_contra_quant')
print(gcf,['figure_exports_fixed9/Retro1_contra_quant_bars'],'-dpdf','-r1000')
close gcf

figure;
imagesc(celltypequant_fraction_norm)
title('Retro1_contra_quant_norm')
print(gcf,['figure_exports_fixed9/Retro1_contra_quant_norm_heatmap'],'-dpdf','-r1000')
close gcf

ordering=[3 2 1 11 12 10 13 4 6 5 7 9 8 14 15 16:19];
figure;
b=bar(celltypequant_fraction_norm(ordering,:),'FaceColor','flat');
for z = 1:6
    b(z).FaceColor = subareapal_sorted_rgb(z,:);
end
title('Retro1_contra_quant_norm')
print(gcf,['figure_exports_fixed9/Retro1_contra_quant_norm_bars'],'-dpdf','-r1000')
close gcf


% save('contra_quantification.mat','quantification','area');



Fnorm_ex_contra=celltypequant_fraction_norm;
F_ex_contra=celltypequant_fraction;

save('quantification_output_contra.mat','Fnorm_ex_contra','F_ex_contra')

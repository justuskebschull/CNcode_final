%% integrate quanitification info from ipsi and contra retro mouse 1 and from mouse oct into one plot
clear all
close all

%% load data
load('quantification_output_contra.mat')
load('quantification_output_ipsi.mat')
load('../../../StarmapOct2019/round2_retakes/quantification_output.mat')
load('../../../StarmapOct2019/round2_retakes/quantification_output_inh.mat')

%% combine

% get subarea pal
subareapal=readtable('subareapal.csv');
subareapal1=table2cell(subareapal(:,2));
subareapal_sorted_rgb=hex2rgb(subareapal1);




F_ex_oct(16:19,:)=F_inh_oct(1:4,:);
Fnorm_ex_oct(16:19,:)=Fnorm_inh_oct(1:4,:);

%make tensor
combined(:,:,1)=F_ex_ipsi;
combined(:,:,2)=F_ex_contra;
combined(:,:,3)=F_ex_oct;

combined_norm(:,:,1)=Fnorm_ex_ipsi;
combined_norm(:,:,2)=Fnorm_ex_contra;
combined_norm(:,:,3)=Fnorm_ex_oct;


mean_F=mean(combined,3);
std_F=std(combined,[],3);

mean_Fnorm=mean(combined_norm,3);
std_Fnorm=std(combined_norm,[],3);


ordering=[3 2 1 11 12 10 13 4 6 5 7 9 8 14 15 16:19];
y=mean_Fnorm(ordering,:);
err=std_Fnorm(ordering,:);

figure;
b=bar(y,'FaceColor','flat');
for z = 1:6
    b(z).FaceColor = subareapal_sorted_rgb(z,:);
end
hold on;
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er=errorbar(x, y(:,i), err(:,i), '.');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
end
hold off
ylim([0 1])
legend

print(gcf,'figure_exports_fixed9/quant_norm','-dpdf','-r1000')

y=mean_F(ordering,:);
err=std_F(ordering,:);
figure;
b=bar(y,'FaceColor','flat');
for z = 1:6
    b(z).FaceColor = subareapal_sorted_rgb(z,:);
end
hold on;
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er=errorbar(x, y(:,i), err(:,i), '.');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
end
hold off
ylim([0 1])
legend

print(gcf,'figure_exports_fixed9/quant','-dpdf','-r1000')




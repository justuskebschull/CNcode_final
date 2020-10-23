clear all;
close all;
%%
cd('C:\Users\justus\Documents\postdoc2\DCN_sequencing\chicken\Starmap\round2-2_sept2020\analysis\');

sections=dir('nisslsegmentation\binary\*.tif');
N=3;
data=[];
for s=[1:numel(sections)]
    tmp=split(sections(s).name,'_');
    sectionname=tmp{1};
    tmp=char(sectionname);    
    short=string(tmp(1:2));
    
    nissl=imread(['nisslsegmentation\binary\',sections(s).name]);
    spotfile=char(['./spots/',sectionname,'/']);
    rawfile=char(['./rawdata/',sectionname,'/']);
    section=char(short);
%      cutoff=[400,500,600,150];
    cutoff=[100,200,600,200];
    
    [data(s).outputmatrix,data(s).genelabels,data(s).cells]=...
        spotstomatrix_chopped(nissl,spotfile,rawfile,N,sectionname,cutoff);
end
save('data.mat','data','-v7.3');

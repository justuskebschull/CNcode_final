clear all;
close all;
%%
%get cell segmentation
sections=dir('nisslsegmentation\binary\*.tif');

exdata=[];
for s=1:numel(sections)
    tmp=split(sections(s).name,'_');
    sectionname=tmp{1};
    tmp=char(sectionname);    
    short=string(tmp(1:2));
    %get starfish output, check for flurescencein te raw images, adn if
    %above threshold count it towards the tally for every cell.
    nissl=imread(['nisslsegmentation\binary\',sections(s).name]);
    spotfile=char(['./spots/',sectionname,'/']);
    rawfile=char(['./rawdata/',sectionname,'/']);
    section=char(short);
     cutoff=[500,500,600,150];
%     cutoff=[600,600,600,200];

    ploting=0;
    
    [exdata(s).outputmatrix,exdata(s).genelabels,exdata(s).cells]=...
        spotstomatrix4(nissl,spotfile,rawfile,section,cutoff,ploting);
end
save('data2.mat','exdata','-v7.3');



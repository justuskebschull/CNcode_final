%% function to define FN,IN,DN by hand on the nissl images.
cd('C:\Users\justus\Documents\postdoc\DCN_sequencing\scRNAseq\STARmap\RetroStarmap\Round1\RetroMouse1_Feb2020\');

sections=dir('nisslsegmentation\raw\*.tif');
exROIs=[];
for s=1:numel(sections)

    nissl=imread(['nisslsegmentation\raw\',sections(s).name]);
    [int2str(s),'FN']
    imshow(nissl,[100,1500])
    exROIs(s).FN=roipoly;
    [int2str(s),'IN']
    exROIs(s).IN=roipoly;
    [int2str(s),'DN']
    exROIs(s).DN=roipoly;
    close all;
end

save('ROIS.mat','exROIs')


%% function to define FN,IN,DN by hand on the nissl images with help of slc17a6channel.
cd('C:\Users\justus\Documents\postdoc\DCN_sequencing\scRNAseq\STARmap\RetroStarmap\Round1\RetroMouse1_Feb2020\');

sections=dir('nisslsegmentation\raw\*.tif');
exROIs=[];
for s=1:numel(sections)
    sectionname=split(sections(s).name,'_');
    S=char(sectionname(1));
    nissl=imread(['nisslsegmentation\raw\',sections(s).name]);
    slc17a6=readtable(['spots/',S,'/',S(1:2),'_R1_C1.csv']);
    penk=readtable(['spots/',S,'/',S(1:2),'_R3_C2.csv']);
    gad1=readtable(['spots/',S,'/',S(1:2),'_R1_C3.csv']);
    imshow(nissl,[100,1500])
    hold on;
    spotcoords=slc17a6{:,[3,4]};
    plot(spotcoords(:,1),spotcoords(:,2),'.','Color','r');
    spotcoords=penk{:,[3,4]};
    plot(spotcoords(:,1),spotcoords(:,2),'.','Color','g');
    spotcoords=gad1{:,[3,4]};
    plot(spotcoords(:,1),spotcoords(:,2),'.','Color','m')
    hold off;
    
    [int2str(s),'FN']
    exROIs(s).FN=roipoly;
    [int2str(s),'IN']
    exROIs(s).IN=roipoly;
    [int2str(s),'DN']
    exROIs(s).DN=roipoly;
    close all;
end

save('ROIS2.mat','exROIs')


%% function to define FN,IN,DN by hand on the nissl images with help of slc17a6channel.
%fix slice 9 annotation.
cd('C:\Users\justus\Documents\postdoc2\DCN_sequencing\scRNAseq\STARmap\RetroStarmap\Round1\RetroMouse1_Feb2020\');
load('ROIS2.mat');
sections=dir('nisslsegmentation\raw\*.tif');
for s=9
    sectionname=split(sections(s).name,'_');
    S=char(sectionname(1));
    nissl=imread(['nisslsegmentation\raw\',sections(s).name]);
    slc17a6=readtable(['spots/',S,'/',S(1:2),'_R1_C1.csv']);
    penk=readtable(['spots/',S,'/',S(1:2),'_R3_C2.csv']);
    gad1=readtable(['spots/',S,'/',S(1:2),'_R1_C3.csv']);
    calb2=readtable(['spots/',S,'/',S(1:2),'_R4_C3.csv']);
    imshow(nissl,[100,1500])
    hold on;
    spotcoords=slc17a6{:,[3,4]};
    plot(spotcoords(:,1),spotcoords(:,2),'.','Color','r');
    spotcoords=penk{:,[3,4]};
    plot(spotcoords(:,1),spotcoords(:,2),'.','Color','g');
    spotcoords=gad1{:,[3,4]};
    plot(spotcoords(:,1),spotcoords(:,2),'.','Color','m')
    spotcoords=calb2{:,[3,4]};
    plot(spotcoords(:,1),spotcoords(:,2),'.','Color','b')
    hold off;
    
    [int2str(s),'FN']
    exROIs(s).FN=roipoly;
    [int2str(s),'IN']
    exROIs(s).IN=roipoly;
    [int2str(s),'DN']
    exROIs(s).DN=roipoly;
    close all;
end

save('ROIS3.mat','exROIs')


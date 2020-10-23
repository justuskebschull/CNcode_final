filelist = getArgument();
file = split(filelist,'#');
open(file[0]);
run("Properties...", "channels=1 slices="+nSlices+" frames=1 unit=pixel pixel_width=0.4 pixel_height=0.4 voxel_depth=0.4");
run("MHD/MHA ...", "save="+file[1]);
close();
run("Quit");

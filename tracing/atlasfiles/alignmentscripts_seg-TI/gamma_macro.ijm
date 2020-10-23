filelist = getArgument();
file = split(filelist,'#');
open(file[0]);
run("Gamma...", "value=0.75 stack");
saveAs("Tiff", file[1]);
close();
run("Quit");


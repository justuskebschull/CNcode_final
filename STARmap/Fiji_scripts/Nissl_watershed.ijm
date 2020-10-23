parentpath=getDirectory("Choose a Directory");
output=getDirectory("_Choose an output Directory");

images=getFileList(parentpath);

for (i=0;i <images.length;i++){
	open(parentpath + images[i]);
	run("Gaussian Blur...", "sigma=3");
	setAutoThreshold("Default");
	//run("Threshold...");
	run("Convert to Mask");
	run("Invert");
	run("Fill Holes");
	run("Watershed");
	name=getTitle;
	saveAs("Tiff", output + name);
	close();
}

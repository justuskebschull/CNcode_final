parentpath=getDirectory("Choose a Directory");
output=getDirectory("_Choose an output Directory");
croppedoutput=getDirectory("Directory for cropped output")
channels=getFileList(parentpath);
numC=5;

//set dapi contrast
rounds=getFileList(parentpath + "/C0");

for (s=0; s < rounds.length; s++){
	open(parentpath + "/C0/" + rounds[s]);
	setMinAndMax(90, 200);
	saveAs("Tiff",parentpath + "/C0/" + rounds[s]);
	close();
}


// registration on the dapi channel
run("Register Virtual Stack Slices", 
	"source=" + parentpath + "/C0 output=" + output +"/C0 feature=Similarity registration=[Similarity           -- translate + rotate + isotropic scale] save");
// "source=" + parentpath + "/C0 output=" + output +"/C0 feature=Rigid registration=[Rigid                -- translate + rotate                  ] save");
run("Duplicate...", "duplicate");
// find shared area and save as roi
setThreshold(0, 1);
setOption("BlackBackground", false);
run("Convert to Mask", "method=Default background=Light");
run("Invert", "stack");
run("Invert LUT");
run("Z Project...", "projection=[Min Intensity]");
run("Create Selection");
run("Make Inverse");
roiManager("Add");
close("Registered C0-1")
close("MIN_Registered C0-1")

// crop to ROI area
selectWindow("Registered C0");
roiManager("Select", 0);
run("Crop");


name=getTitle;
 
//define a function to separately save the different rounds.
function saveRounds(croppedoutput,name,parentpath,channel) {
	run("Stack to Images");
	while (nImages()>0){		
		windowname=getTitle;
		saveAs("Tiff", croppedoutput + "/" + channel + "/" + windowname + ".tif");
		close();
		}
}

//save cropped output
saveRounds(croppedoutput,name,parentpath,"C0");


// now apply same transformation to other channels
for (c=1;c<numC;c++){
run("Transform Virtual Stack Slices", 
	"source=" + parentpath + "/C" + c + " output=" + output + "/C" + c + 
	" transforms=" + parentpath +"/C0 interpolate");

roiManager("Select", 0);
run("Crop");
currentchannel= "C" + c;
name2=getTitle;
saveRounds(croppedoutput,name2,parentpath,currentchannel);

}
run("Close");


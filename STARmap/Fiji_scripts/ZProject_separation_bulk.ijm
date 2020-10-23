parentpath=getDirectory("Choose a Directory");
output=getDirectory("_Choose an output Directory");
slices=getFileList(parentpath);

for (r=0; r<slices.length;r++){
rounds=getFileList(parentpath + slices[r]);

for (s=0; s < rounds.length; s++){
	
	if(rounds[s]=='Nissl/'){numChannels =2;}
	numC=5;
	open(parentpath + slices[r] + rounds[s]);
//	if (s==1)
//	{
//	run("Flip Horizontally");
//	}
//if (s==3)
//	{
//	run("Flip Horizontally");
//	}
name=getTitle;
print(name);
newname = substring(name, 0, lengthOf(name)-4);
run("Split Channels");


if (s==0)
{

selectWindow("C1-" + name);
saveAs("Tiff", output + slices[r] + "/C0/" + newname + "_Cnissl.tif");
close();
selectWindow("C2-" + name);
saveAs("Tiff", output + slices[r] + "/C1/" + newname + "_Cnissl.tif");
saveAs("Tiff", output + slices[r] + "/C2/" + newname + "_Cnissl.tif");
saveAs("Tiff", output + slices[r] + "/C3/" + newname + "_Cnissl.tif");
saveAs("Tiff", output + slices[r] + "/C4/" + newname + "_Cnissl.tif");
close();
}
if (s>0){
for (c=1; c<=numC;c++){
selectWindow("C" + c +"-" + name);
saveAs("Tiff", output + slices[r] + "/C" + c-1 +"/" + newname + "_C" + c-1 + ".tif");
close();
}

}
}
}
function zProjChannels(parentpath,outpath,filepath,filelist,numC, numZ,R) {
	zname = newArray(numC+1);
	for (c = 1; c <= numC; c++) {
		r=(c-1)*numZSlices+1;
		
		input = "open=" + filepath + " number=" + numZ + " starting=" + r + " sort";
		
		run("Image Sequence...", input);
		name_tmp1=getTitle; 

		run("Z Project...", "projection=[Max Intensity]");
		zname[c]=getTitle;
		close(name_tmp1);
		print(name_tmp1);
	}
	
	if (numC==5){
	run("Merge Channels...", "c1=" + zname[1] + 
		" c2=" + zname[2] +
		" c3=" + zname[3] +
		" c4=" + zname[4] +
		" c5=" + zname[5] +
		" create");
	newname1=substring(name_tmp1,7,lengthOf(name_tmp1));
//	newname1=name_tmp1;
	}
	if(numC==2){
	run("Merge Channels...", "c1=" + zname[1] + 
		" c2=" + zname[2] +
		" create");
	newname1=substring(name_tmp1,6,lengthOf(name_tmp1));
//	newname1=name_tmp1;
	}
	newname2 = substring(R, 0, lengthOf(R)-1);
	saveAs("Tiff", outpath + newname1 + '/' + newname1+ '_' + newname2 + ".tif");
	close();
}


parentpath=getDirectory("Choose a Directory");
outpath=getDirectory("Output");
rounds=getFileList(parentpath);

for (r=0; r <rounds.length; r++){
numChannels = 5;
print(rounds[r]);
if(rounds[r]=='Nissl/'){numChannels =2;}
sections=getFileList(parentpath+rounds[r]);

for (s=0; s < sections.length; s++){

filepath = parentpath + rounds[r] + sections[s];
list = getFileList(filepath);

numZSlices = list.length / numChannels;


zProjChannels(parentpath,outpath,filepath,list, numChannels, numZSlices,rounds[r]);
}
}
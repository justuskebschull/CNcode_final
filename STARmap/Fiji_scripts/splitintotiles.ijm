parentpath=getDirectory("Choose a Directory");
output=getDirectory("_Choose an output Directory");
channels=getFileList(parentpath);


for (s=1; s < channels.length; s++){

	files=getFileList(parentpath + channels[s]);
	
	for (f=0; f < files.length; f++){
		open(parentpath + channels[s] + files[f]);
		name=getTitle;
		run("Montage to Stack...", "columns=3 rows=3 border=0");
		run("Stack to Images");
		close(name);
		newname = substring(name, 4, lengthOf(name));
			for (c=1; c<=9;c++){
				if(c<10){
					selectWindow("Stack-000" + c);
				} else {
			selectWindow("Stack-00" + c);
			}
			saveAs("Tiff", output + "/S" + c +"L_" + newname);
			close();
}
}
}

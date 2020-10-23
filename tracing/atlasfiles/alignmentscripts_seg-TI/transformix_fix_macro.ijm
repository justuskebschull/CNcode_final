filelist = getArgument();
file = split(filelist,'#');
open(file[0]);
run("Duplicate...", "title=maskA duplicate");

setThreshold(-32766, -1);

run("Convert to Mask","method=Defaut background");


run("Duplicate...", "title=maskB duplicate");
run("Invert", "stack");
run("Invert LUT");
selectWindow("maskA");
run("Invert LUT");
selectWindow("maskB");
run("Divide...", "value=255.000 stack");
selectWindow("maskA");
run("Divide...", "value=255.000 stack");
imageCalculator("Multiply create 32-bit stack", "maskA","result.raw");
selectWindow("Result of maskA");
setMinAndMax(-32766, -1);
run("Abs", "stack");
setMinAndMax(0,32766);
run("Invert", "stack");
run("Add...", "value=32767 stack");

imageCalculator("Multiply create 32-bit stack", "Result of maskA","maskA");
imageCalculator("Multiply create 32-bit stack", "maskB", "result.raw");


imageCalculator("Add create 32-bit stack", "Result of maskB","Result of Result of maskA");
selectWindow("Result of Result of maskB");
setMinAndMax(0, 65535);
run("16-bit");
setMinAndMax(0, 65535);
saveAs("Tiff", file[1]);
close();
run("Quit");

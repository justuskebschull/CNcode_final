#!/usr/bin/env bash
brains=($1/brain*)
atlas=$2


aligner () {
    a=$1
    b=$2
    autofluofile="${a}/autofluo_resampled_25um_gamma.tif"
    signalfile="${a}/signal_resampled_25um.tif"
	signalfile2="${a}/seg-TI_resampled_25um.tif"
    elastix_param1="/home/justus/Documents/justus_synology/atlasfiles/elastix_parameters/Par0000affine.txt"
    elastix_param2="/home/justus/Documents/justus_synology/atlasfiles/elastix_parameters/Par0000bspline.txt"
    elastixoutput="${a}/elastix_output"
   
    mkdir $elastixoutput

#run elastix
    elastix -f $2 -m $autofluofile -out $elastixoutput -p $elastix_param1 -p $elastix_param2

#change parameters
    elastixfile1="${elastixoutput}/TransformParameters.0.txt"
    elastixfile2="${elastixoutput}/TransformParameters.1.txt"
    sed -i 's/^.*(FinalBSplineInterpolationOrder 3)$/(FinalBSplineInterpolationOrder 1)/' $elastixfile1
    sed -i 's/^.*(FinalBSplineInterpolationOrder 3)$/(FinalBSplineInterpolationOrder 1)/' $elastixfile2
    #run transformix
    transformixoutput="${a}/transformix_output"
    mkdir $transformixoutput
    transformix -out $transformixoutput -in $signalfile -tp $elastixfile2
	
	#run transformix for seg-TI
    transformixoutput2="${a}/transformix_output_seg-TI"
    mkdir $transformixoutput2
    transformix -out $transformixoutput2 -in $signalfile2 -tp $elastixfile2
	
	
}

fullsize () {
    directory=$1
    #change parameters
    sed -i 's/^.*(Spacing 1.0000000000 1.0000000000 1.0000000000)$/(Spacing 0.4000000000 0.4000000000 0.4000000000)/' $directory/elastix_output/TransformParameters.1.txt
    sed -i 's/^.*(Spacing 1.0000000000 1.0000000000 1.0000000000)$/(Spacing 0.4000000000 0.4000000000 0.4000000000)/' $directory/elastix_output/TransformParameters.0.txt
    sed -i 's/^.*(Size 354 354 260)$/(Size 885 885 650)/' $directory/elastix_output/TransformParameters.1.txt
    sed -i 's/^.*(Size 354 354 260)$/(Size 885 885 650)/' $directory/elastix_output/TransformParameters.0.txt
    
    #run transformix
    transformixoutput="${directory}/transformix_output_10um_647"
    mkdir $transformixoutput
    transformix -out $transformixoutput -in $directory/signal_resampled_10um.mha -tp $directory/elastix_output/TransformParameters.1.txt

    #run transformix for seg-TI output
    transformixoutput="${directory}/transformix_output_10um_seg-TI"
    mkdir $transformixoutput
    transformix -out $transformixoutput -in $directory/seg-TI_resampled_10um.mha -tp $directory/elastix_output/TransformParameters.1.txt


    
}

fileformat () {
    /home/justus/Documents/Fiji.app/ImageJ-linux64 --headless -macro /home/justus/Documents/justus_synology/iDisco/alignmentscripts/mha_maker.ijm $1#$2
}

gamma () {
    /home/justus/Documents/Fiji.app/ImageJ-linux64 --headless -macro /home/justus/Documents/justus_synology/iDisco/alignmentscripts/gamma_macro.ijm $1#$2
}

holes () {
    /home/justus/Documents/Fiji.app/ImageJ-linux64 --headless -macro /home/justus/Documents/justus_synology/iDisco/alignmentscripts/transformix_fix_macro.ijm $1#$2
}

echo $brains
for ((i=0;i<${#brains[@]};i++))

do
     gamma ${brains[i]}/autofluo_resampled_25um.tif ${brains[i]}/autofluo_resampled_25um_gamma.tif


     aligner ${brains[i]} ${atlas}
     

     holes ${brains[i]}/transformix_output/result.mhd ${brains[i]}/transformix_output/result_fixed.tif
     holes ${brains[i]}/transformix_output_seg-TI/result.mhd ${brains[i]}/transformix_output_seg-TI/result_fixed.tif


     fileformat ${brains[i]}/signal_resampled_10um.tif ${brains[i]}/signal_resampled_10um.mha

     fileformat ${brains[i]}/seg-TI_resampled_10um.tif ${brains[i]}/seg-TI_resampled_10um.mha
     

     fullsize ${brains[i]}

     holes ${brains[i]}/transformix_output_10um_647/result.mhd ${brains[i]}/transformix_output_10um_647/result_fixed.tif

     holes ${brains[i]}/transformix_output_10um_seg-TI/result.mhd ${brains[i]}/transformix_output_10um_seg-TI/result_fixed.tif

done

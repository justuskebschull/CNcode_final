

######################### Import modules

import os
import scipy.stats

import ClearMap.Settings as settings
import ClearMap.IO as io

from ClearMap.Alignment.Resampling import resampleData;


    
    #Resolution of the Raw Data (in um / pixel)
OriginalResolution = (4.0625, 4.0625, 3);

    #Orientation: 1,2,3 means the same orientation as the reference and atlas files.
    #Flip axis with - sign (eg. (-1,2,3) flips x). 3D Rotate by swapping numbers. (eg. (2,1,3) swaps x and y)
FinalOrientation = (1,2,3);

    #Resolution of the Atlas (in um/ pixel)
AtlasResolution = (25, 25, 25);
    
    ########################### Config parameters

    #Processes to use for Resampling (usually twice the number of physical processors)
ResamplingParameter = {
        "processes": 24 
        };


    #Files for Auto-fluorescence (Atlas Registration)
RegistrationResamplingParameter488 = ResamplingParameter.copy();
RegistrationResamplingParameter488["source"] = AutofluoFile;
RegistrationResamplingParameter488["sink"] =  os.path.join(BaseDirectory, 'autofluo_resampled_25um.tif');
RegistrationResamplingParameter488["resolutionSink"]  = AtlasResolution;   
RegistrationResamplingParameter488["resolutionSource"] = OriginalResolution;
RegistrationResamplingParameter488["orientation"] = FinalOrientation;
    
RegistrationResamplingParameter647 = RegistrationResamplingParameter488.copy();
RegistrationResamplingParameter647["source"]= sigFile;
RegistrationResamplingParameter647["sink"]=os.path.join(BaseDirectory,'signal_resampled_25um.tif');
    
RegistrationResamplingParameterred = RegistrationResamplingParameter488.copy();
RegistrationResamplingParameterred["source"]= sigFile2;
RegistrationResamplingParameterred["sink"]=os.path.join(BaseDirectory,'seg-TI_resampled_25um.tif');

    
    ##################### run resampling
resampleData(**RegistrationResamplingParameter488);
resampleData(**RegistrationResamplingParameter647);
resampleData(**RegistrationResamplingParameterred);






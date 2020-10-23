#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 12:56:11 2018
Edited May21 2020 
@author: Justus
"""

import os
import subprocess
parentDir = '/home/justus/Documents/justus_synology/iDisco/round8/brainsforanalysis/DN'
Atlas='/home/justus/Documents/justus_synology/atlasfiles/horizontalsections/brighttracts.tif'

for i in os.walk(os.path.join(parentDir,'.')).next()[1]:
        BaseDirectory = parentDir+'/'+i;
        print 'Brain '+ i
        folders = os.walk(os.path.join(BaseDirectory,'.')).next()[1]
        auto = [f for f in folders if ('488' in f)][0]
        signal = [f for f in folders if ('647' in f)][0]
        signal2 = [f for f in folders if ('seg-TI' in f)][0]
       
#        print signal2
        
        affilename = os.listdir(BaseDirectory+'/'+auto)[0].split(' ')[0]
        sffilename = os.listdir(BaseDirectory+'/'+signal)[0].split(' ')[0]
        sf2filename = os.listdir(BaseDirectory+'/'+signal2)[0].split(' ')[0]
       
        print affilename
        print sffilename
        print sf2filename
        
        sigFile = os.path.join(BaseDirectory, signal+'/'+sffilename+' Z\d{4}.ome.tif');
        sigFile2 = os.path.join(BaseDirectory, signal2+'/TI_\d{4}.tif');
        AutofluoFile = os.path.join(BaseDirectory, auto+'/'+affilename+' Z\d{4}.ome.tif');
        
#        print sigFile2
#        print 'Brain '+ i
       
        try:
            execfile('/home/justus/Documents/justus_synology/iDisco/alignmentscripts_seg-TI/downsampling_script_25um.py')
        except Exception as e:
            print(e)

        try:
            execfile('/home/justus/Documents/justus_synology/iDisco/alignmentscripts_seg-TI/downsampling_script_10um.py')
        except Exception as e:
            print(e)
            
subprocess.call(['/home/justus/Documents/justus_synology/iDisco/alignmentscripts_seg-TI/alignment_wholebrain_newcap.sh',parentDir,Atlas])




# -*- coding: utf-8 -*-
"""
Created on Wed, Sept 9, 2020

This is a module defining global values for the program

@author: Qingyu.Feng
"""

##########################################################################
# Import modules #########################################################
##########################################################################
import os, sys

##########################################################################
# Define classes #########################################################
##########################################################################

useRay = False
iFlagCopy = False
iRunSWAT = True
# Only used in applyBestParmSet and sa to facilitate testing
iModParm = True

# Folder structure
fdProjSetup = "01projSetupContPara"
fdObs = "02observedData"
fdgisLayers = "03gisLayers"
fdmodelTxtInOut = "04SWATTxtInOut"
fdDMPOTpyFiles = "pyscripts"
fdWorkingDir = "05workingDir"
fdOutputs = "06outputFiles"
fdCalibrated = "07calibratedTIO"
fdSA = "08sensitivityAnalysis"

if not os.path.isdir(fdWorkingDir):
    os.mkdir(fdWorkingDir)
if not os.path.isdir(fdOutputs):
    os.mkdir(fdOutputs)

fdMain = os.getcwd()


platforms = {
    'linux' : 'Linux',
    'darwin' : 'OS X',
    'win32' : 'Windows'
}


osplatform = platforms[sys.platform]

if osplatform == "Linux":
    fnSwatExe = "swat2012.681.gfort.rel"                                
elif osplatform == "Windows":
    fnSwatExe = "swat2012rev670.exe"





# Control File JSON
fnCtrlSetUsr = os.path.join(fdProjSetup, "dmpot_Control.set")
fnCtrlJsonUsr = os.path.join(fdOutputs, "usrControl.json")

# Parameter set names
fnParmSetComb = os.path.join(fdProjSetup, "dmpot_Para_Combined.set")
fnParmUsed = os.path.join(fdOutputs, "usrParmInCal.csv")

# Reach shapefile path-name
fnReachShp = os.path.join(fdgisLayers, "reach.shp")

# Row crop and pasture list
rowCropLst = ["CORN", "SOYB", "WWHT", "AGRL", "AGRR", "AGRC",
                "CCRN", "CSOY", "CWHT", "CYCN", "SYCN", "SYWH", "WHCN"]

pastHayLst = ["PAST", "HAY", "ALFA"]
forestLst = ["FRST", "FRSD", "FRSE"]
urbanLst = ["URHD", "URMD", "URML", "URLD", "UCOM", "UIDU", "UTRN", "UINS"]

# OutVariable header and numbers in the observed files
# In the Control File:
# Output Variables (0-15): '0- Baseflow', '1-Stream Flow',
# '2-Sediment','3-org-N','4-org-P','5-NO_3-N','6-NH_4-N',
# '7-NO_2-N','8-Mineral-P','9-Soluble Pesticide',
# '10-Sorbed Pesticide','11-Total Phosphorus',
# '12-Total Nitrogen','13-Total Pesticide','14-TKN',
#  '15-NO2+NO3'
# Headers in the Observed data:
# sf(m3/s)	sed(t/ha)
# orgn(kg/ha)	orgp(kg/ha)	no3n(kg/ha)
# 	nh4n(kg/ha)	no2n(kg/ha)	minp(kg/ha)
# 	solpst(mg/ha)	sorpst(mg/ha)	tp(kg/ha)
# 	tn(kg/ha)	tpst(ppb)
varIDObsHdrPair = {
1 : "sf(m3/s)",
2 : "sed(t/ha)",
3 : "orgn(kg/ha)",
4 : "orgp(kg/ha)",
5 : "no3n(kg/ha)",
6 : "nh4n(kg/ha)",
7 : "no2n(kg/ha)",
8 : "minp(kg/ha)",
9 : "solpst(mg/ha)",
10 : "sorpst(mg/ha)",
11 : "tp(kg/ha)",
12 : "tn(kg/ha)",
13 : "tpst(ppb)"
}


# Output variables determined in the file.cio file
# 1bfr,2sf,3sed,4orgn,5orgp,6no3n,7nh4n,8no2n,9minp,
# 10solpst,11sorpst,12tp,13tn,14tpst
# rchVarSelLst = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]
rchVarLst = [2, 6, 9, 11, 13, 15, 17, 19, 27, 29]

# if calibration of sediements,nutrient, or pesticides is based on
# concentrations, out_id will be 1 since SF will be needed to compute conc.
# set all to zero to get the loads: Sediment (ton), P & N (kg), Pesticide (mg) 
# cliabrateConc = [0, 0, 0]

bsnLvlFlExt = [".bsn", "crop.dat", ".wwq"]
subLvlFlExtLst = [".sub", ".rte", ".swq", ".res"]
hruLvlFlExtLst = [".gw", ".hru", ".mgt", ".sol", ".chm"] 


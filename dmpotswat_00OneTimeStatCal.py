#!/usr/bin/python3

"""
This program was developed to calculate the statistics of the default
simulation without any modification and calibration.

"""

##########################################################################
# Import modules #########################################################
##########################################################################

import sys
import math, numpy
import datetime
import pytz

from pyscripts.globVars import *
from pyscripts.DMPOTUtil import *


##########################################################################
# Program start  #########################################################
##########################################################################

print("=================================================================")
print("          Distributed Model Parameter Optimization Tool          ")
print("                   Developed by Qingyu Feng                      ")
print("          Research Center for Eco-Enviormental Sciences          ")
print("                  Chinese Academy of Science                     ")
print("                   Last Updated: 2021-04-22                      ")
print("=================================================================")


timeZone = pytz.timezone('Asia/Shanghai')
startTime = datetime.datetime.now(timeZone)
todayDate = datetime.date.today()

nCalVal = "Calibration"

##########################################################################
# Read Settings ##########################################################
print(".....Read in DMPOT setting from the control file.....")
ctrlSetting = ctrlSetToJSON(fnCtrlSetUsr, fnCtrlJsonUsr)

if (ctrlSetting == "Error"):
    print("Your control Setting was not corrected prepared. Please double check!")
    sys.exit(1)

##########################################################################
# Read observed data and preprocessing ###################################
print(".....Read in Observed data from file.....")
obsDataLst = readOBSSet(ctrlSetting, fdObs, varIDObsHdrPair)

##########################################################################
# Copy TxtInOut contents to working Dir ##################################
if iFlagCopy:
    copySWATTxtIO(timeZone, startTime, fdmodelTxtInOut, fdWorkingDir)

# After copying modify the file.cio to make sure the simulated output contains
# interested outlet Variables.
iPrintForCio = modFileCio(ctrlSetting, fdWorkingDir, rchVarLst)

##########################################################################
# Read Parameter and preprocessing #######################################
print(".....Read in DMPOT parameter from the parameter file.....")
parmBsnLvl, parmBsnLvlFExtLst, parmSubLvl, parmSubLvlFExtLst = getParmSets()

##########################################################################
# Process subarea groups based on user setting ###########################
# Determine whether user select to group subareas based
# on outlet.
# The output will be a subGroups.
# Even if the users select not to group, the structure will be
# the same to maintain a consistent structure.
# The subGroups and subParSetBest has the same key.
# The design of the key for the dictionary is important, since it will be 
# used across the running process to the groups of subarea list, the parameter
# set and the calculation of statistics.
# Sometimes, the user might set groupSubareaIndex == 1 even they
# only have one outlet. This situation will be corrected.
if (len(ctrlSetting["outLetList"]) == 1):
    ctrlSetting["groupSubareaIdx"] = 0

subNoGroups, parmObjFnKeys, rcvRchLst = getSubGroupsRchList(ctrlSetting)

subParaFns, subObjFunFns, subParaSel01Fns = initOutFileParmObjSublvl(
    ctrlSetting, parmSubLvl, parmObjFnKeys, rcvRchLst)
    
fnBsnPara, fnBsnParaSel = initOutFileParmObjBsnlvl(parmBsnLvl)

subObfBestDict,subObfTestDict, bsnObfBest = initOFValDict(parmObjFnKeys)

subParGroups, swatSubFnGroups, swatHruFnGroups = initParmInFilenameSubLvl(
    subNoGroups, parmSubLvl, fdWorkingDir)

##########################################################################
# Start optimization procedure for all runs ##############################
##########################################################################
print(".....DMPOT opitmization procedure start.....")

##########################################################################
##########################################################################
runIdx = 999999
# After modifying, run the SWAT model
if iRunSWAT:
    runCode = runSWATModel(get_osplatform(), fdWorkingDir,
                    fdMain, fdDMPOTpyFiles, fnSwatExe)
    runEndTime = datetime.datetime.now(timeZone)
    runTimeTotal = runEndTime - startTime
    print(".....Time for running SWAT: {}; Total Time: {}.....". format(
        runTimeTotal,
        runEndTime - startTime)) 
else:
    runEndTime = datetime.datetime.now(timeZone)
    runTimeTotal = runEndTime - startTime

# Calculate Statistics 
subObfTestDict, obfValNonOther, subAllStats = calObjFuncValues(iPrintForCio, 
            rcvRchLst, 
            obsDataLst, 
            parmObjFnKeys,
            subObfTestDict, 
            fdWorkingDir, True
            , nCalVal, ctrlSetting["groupSubareaIdx"])

totalTimeThisrun = runEndTime - startTime
# Update the best parameter values and objective functions
# based on the objective Functions
probVal = 1
subParGroups, subObfBestDict, bsnObfBest, parmBsnLvl = updateBestParm(ctrlSetting, 
            subParGroups, 
            subObfBestDict, 
            subObfTestDict,
            subParaFns,
            runIdx,
            obfValNonOther,
            parmObjFnKeys,
            subAllStats,
            subObjFunFns, 
            totalTimeThisrun,
            probVal,
            subParaSel01Fns,
            parmBsnLvl,
            bsnObfBest,
            fnBsnPara,
            fnBsnParaSel
            )

calStatTime = datetime.datetime.now(timeZone)
print("Time for calculating statistics: {}; Total Time: {}". format(
        calStatTime - runEndTime,
        calStatTime - startTime)) 
print("=================================================================")


# End of for loop for total runs
print("--------------------------------------")
print("Congratulations!!! It's done nicely~~~")
print("--------------------------------------")

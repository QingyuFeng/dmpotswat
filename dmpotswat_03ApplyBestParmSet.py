#!/usr/bin/python3

"""
This program was developed to apply the best parameter from the user
selected run numbers

"""

##########################################################################
# Import modules #########################################################
##########################################################################

import sys
import datetime
import pytz

from pyscripts.globVars import *
from pyscripts.DMPOTUtil import *


nCalVal = "Validation"
fdCalibrated = fdCalibrated + nCalVal
if not os.path.isdir(fdCalibrated):
    os.mkdir(fdCalibrated)
    
##########################################################################
# Program start  #########################################################
##########################################################################

print("=================================================================")
print("          Distributed Model Parameter Optimization Tool          ")
print("                   Developed by Qingyu Feng                      ")
print("          Research Center for Eco-Enviormental Sciences          ")
print("                  Chinese Academy of Science                     ")
print("=================================================================")


timeZone = pytz.timezone('Asia/Shanghai')
startTime = datetime.datetime.now(timeZone)
todayDate = datetime.date.today()
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
    copySWATTxtIO(timeZone, startTime, fdmodelTxtInOut, fdCalibrated)

# After copying modify the file.cio to make sure the simulated output contains
# interested outlet Variables.
iPrintForCio = modFileCio(ctrlSetting, fdCalibrated, rchVarLst)

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

subParaFns, subParaCalFns, subObjFunCalFns, subParaSel01Fns = initOutFileParmObjSublvlCal(
    ctrlSetting, parmSubLvl, parmObjFnKeys, rcvRchLst, nCalVal)
    
fnBsnPara, fnBsnParaCal, fnBsnParaSel = initOutFileParmObjBsnlvlCal(parmBsnLvl, nCalVal)

subObfBestDict,subObfTestDict, bsnObfBest = initOFValDict(parmObjFnKeys)

subParGroups, swatSubFnGroups, swatHruFnGroups = initParmInFilenameSubLvl(
    subNoGroups, parmSubLvl, fdCalibrated)

# Get the user best bsn lvl and subarea level parameters
print(".....Getting best paramter values from the user selected run No.....")

usrBestParmBsnDf, usrBestParGroups = getUsrBestParmSet(
            fnBsnPara, 
            ctrlSetting["runNoBestPar"], subParGroups, 
            subParaFns, subNoGroups, parmSubLvl)

##########################################################################
# Start optimization procedure for all runs ##############################
##########################################################################
print(".....DMPOT opitmization procedure start.....")

##########################################################################
##########################################################################
runIdx = 9999999
print(".....Modifying subarea level with best parameters.....")
modifyStartTime = datetime.datetime.now(timeZone)
if len(parmSubLvl.index) > 0:
    for subGPKey, subGPL in subParGroups.items():
        # Update parameter with best parameter values 
        for parIdx in subGPL.index:
            parSymb = subGPL.loc[parIdx,"Symbol"]
            subGPL.loc[parIdx,"TestVal"] = usrBestParGroups[subGPKey][parSymb]  

        if iModParm:  
            # Update parameter values in file, run model, and calculate
            # objective function
            modifyParInFileSub(subGPL, 
                            parmSubLvlFExtLst, 
                            swatSubFnGroups, 
                            subGPKey, 
                            swatHruFnGroups,
                            fdCalibrated)
if len(parmBsnLvl.index) > 0:
    print(".....Modifying basin level with best parameters.....")
    for parBIdx in parmBsnLvl.index:
        parSymb2 = parmBsnLvl.loc[parBIdx,"Symbol"]
        parmBsnLvl.loc[parBIdx,"TestVal"] = usrBestParmBsnDf[parSymb2]
    # After modifying parameter values in file, 
    modifyParInFileBsn(parmBsnLvl, 
                        parmBsnLvlFExtLst,
                        fdCalibrated)
modifyEndTime = datetime.datetime.now(timeZone)
modifyTimeTotal = modifyEndTime - modifyStartTime
print(".....Time for modifying parameter values: {}; Total Time: {}.....". format(
    modifyTimeTotal,
    modifyEndTime - startTime))
    
# After modifying, run the SWAT model
if iRunSWAT:
    runCode = runSWATModel(get_osplatform(), fdCalibrated,
                    fdMain, fdDMPOTpyFiles, fnSwatExe)
    runEndTime = datetime.datetime.now(timeZone)
    runTimeTotal = runEndTime - modifyEndTime
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
            fdCalibrated, True, nCalVal, ctrlSetting["groupSubareaIdx"])

totalTimeThisrun = runEndTime - startTime
# Update the best parameter values and objective functions
# based on the objective Functions
probVal = 999999
subParGroups, subObfBestDict, bsnObfBest, parmBsnLvl =updateBestParm(ctrlSetting, 
            subParGroups, 
            subObfBestDict, 
            subObfTestDict,
            subParaCalFns,
            runIdx,
            obfValNonOther,
            parmObjFnKeys,
            subAllStats,
            subObjFunCalFns, 
            totalTimeThisrun,
            probVal,
            subParaSel01Fns,
            parmBsnLvl,
            bsnObfBest,
            fnBsnPara,
            fnBsnParaSel
            )

# Write the basin level parameter
writeBsnParm(parmBsnLvl,
            runIdx,
            fnBsnParaCal, 
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

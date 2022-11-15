#!/usr/bin/python3

"""
This program was developed to calibrate the SWAT model for large watershed.
The general steps include:
1. Read in dmpot settings
2. Read in observed data
3. Read in model parameters
4. Optimization routines
5. Compute statistics
6. Ending criteria

Updated Nov 26, 2020 by Qingyu Feng
1. Added the making plot of simulated vs observed after each run.

Updated Oct 27, 2021
1. Due to the poorer output by Dist mode, the criteria to update
best parameter set was changed from individual statistic to the sum
of statistics of all outlets 



Developed by Qingyu Feng
RCEES
Last Update: Oct 30, 2020

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

nCalVal = "Calibration"
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

subObfBestDict, subObfTestDict, bsnObfBest = initOFValDict(parmObjFnKeys)

subParGroups, swatSubFnGroups, swatHruFnGroups = initParmInFilenameSubLvl(
    subNoGroups, parmSubLvl, fdWorkingDir)

##########################################################################
# Start optimization procedure for all runs ##############################
##########################################################################
print(".....DMPOT opitmization procedure start.....")

##########################################################################
##########################################################################
# Start running for the main running loops

# Two situations:
# User select to use random parameter initParmIdx == 0
# User select to use initial parameter initParmIdx == 1
initRunNo = math.ceil(0.005 * ctrlSetting["totalModelRuns"])
totalRuns = ctrlSetting["totalModelRuns"] - initRunNo

# Initial a counter to record the runs of random
iCall = 0
# Using Random parameter
if (ctrlSetting["initParmIdx"] == 0):
    for runIdx in range(initRunNo):
        print(".....DMPOT simulation NO: {}.....".format(runIdx+1))
        # Mainly update the testValue in parameterSelected
        # If the user selected subarea parameters
        modifyStartTime = datetime.datetime.now(timeZone)
        iCall = runIdx + 1

        if len(parmSubLvl.index) > 0:
            # print(".....Modifying subarea level parameter values Randomly.....")
            # Generate random values for parameter updating
            totalNoSelParSub = parmSubLvl.shape[0]
            
            for subGPKey, subGPL in subParGroups.items():
                # Update parameter completely random. 
                ranNumSub = numpy.random.rand(1, totalNoSelParSub)
                subGPL = generateRandomParVal(subGPL, ranNumSub)
                # Update parameter values in file, run model, and calculate
                # objective function
                modifyParInFileSub(subGPL, 
                                parmSubLvlFExtLst, 
                                swatSubFnGroups, 
                                subGPKey, 
                                swatHruFnGroups,
                                fdWorkingDir)
        if len(parmBsnLvl.index) > 0:
            # print(".....Modifying basin level parameter values Randomly.....")
            # Generate random values for parameter updating
            totalNoSelParBsn = parmBsnLvl.shape[0]
            ranNumBsn = numpy.random.rand(1, totalNoSelParBsn)

            parmBsnLvl = generateRandomParVal(parmBsnLvl, ranNumBsn)
            # After modifying parameter values in file, 
            modifyParInFileBsn(parmBsnLvl, 
                                parmBsnLvlFExtLst,
                                fdWorkingDir)
        modifyEndTime = datetime.datetime.now(timeZone)
        modifyTimeTotal = modifyEndTime - modifyStartTime
        print(".....Time for modifying parameter values: {}; Total Time: {}.....". format(
            modifyTimeTotal,
            modifyEndTime - startTime))
    
        # After modifying, run the SWAT model
        if iRunSWAT:
            runCode = runSWATModel(get_osplatform(), fdWorkingDir,
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
                    fdWorkingDir, False, nCalVal, ctrlSetting["groupSubareaIdx"])
        
        totalTimeThisrun = runEndTime - startTime
        # Update the best parameter values and objective functions
        # based on the objective Functions
        probVal = 1
        subParGroups, subObfBestDict, bsnObfBest, parmBsnLvl = updateBestParm(
                    ctrlSetting, 
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


# Using Initial parameter
elif (ctrlSetting["initParmIdx"] == 1):
    runIdx = 0 # Only one time
    print(".....DMPOT simulation NO: {}.....".format(runIdx+1))
    # Mainly update the testValue in parameterSelected
    # If the user selected subarea parameters
    print(".....Modifying parameter values using Default Parameter.....")
    modifyStartTime = datetime.datetime.now(timeZone)
    iCall = runIdx + 1
    totalRuns = ctrlSetting["totalModelRuns"] - iCall
    if len(parmSubLvl.index) > 0:
        for subGPKey, subGPL in subParGroups.items():
            # Do not Update parameter and directly using the default values. 
            # Update parameter values in file, run model, and calculate
            # objective function
            modifyParInFileSub(subGPL, 
                            parmSubLvlFExtLst, 
                            swatSubFnGroups, 
                            subGPKey, 
                            swatHruFnGroups,
                            fdWorkingDir)
    if len(parmBsnLvl.index) > 0:
        # After modifying parameter values in file, 
        modifyParInFileBsn(parmBsnLvl, 
                            parmBsnLvlFExtLst,
                                fdWorkingDir)
    modifyEndTime = datetime.datetime.now(timeZone)
    modifyTimeTotal = modifyEndTime - modifyStartTime
    print(".....Time for modifying parameter values: {}; Total Time: {}.....". format(
            modifyTimeTotal,
            modifyEndTime - startTime))

    # After modifying, run the SWAT model
    if iRunSWAT:
        runCode = runSWATModel(get_osplatform(), fdWorkingDir, fdMain, fdDMPOTpyFiles, fnSwatExe)
        runEndTime = datetime.datetime.now(timeZone)
        runTimeTotal = runEndTime - modifyEndTime
        print(".....Time for running SWAT: {}; Total Time: {}.....". format(
                runTimeTotal,
                runEndTime - startTime))  
    else:
        runEndTime = datetime.datetime.now(timeZone)

    # Calculate Statistics 
    subObfTestDict, obfValNonOther, subAllStats = calObjFuncValues(iPrintForCio, 
                rcvRchLst, 
                obsDataLst, 
                parmObjFnKeys,
                subObfTestDict,
                fdWorkingDir, False, nCalVal, ctrlSetting["groupSubareaIdx"])

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


# Start the main loop
for runIdx in range(1, totalRuns+1):
    print(".....DMPOT simulation NO: {}.....".format(runIdx + initRunNo))
    # Update parameter using DDS
    # Mainly update the testValue in parameterSelected
    # If the user selected subarea parameters
    modifyStartTime = datetime.datetime.now(timeZone)

    # Calculate the probability value of each run over total runs
    probVal = 1.0-(numpy.log(runIdx)/numpy.log(totalRuns))
    if len(parmSubLvl.index) > 0:
        print(".....Modifying subarea level parameter values with DDS.....")
        for subGPKey, subGPL in subParGroups.items():
            # Update parameter using DDS. 
            subGPL = generateDDSParVal(subGPL, 
                                        probVal, 
                                        ctrlSetting["perturbFactor"])
            # Update parameter values in file, run model, and calculate
            # objective function
            modifyParInFileSub(subGPL, 
                            parmSubLvlFExtLst, 
                            swatSubFnGroups, 
                            subGPKey, 
                            swatHruFnGroups,
                            fdWorkingDir)

    if len(parmBsnLvl.index) > 0:
        print(".....Modifying basin level parameter values with DDS.....")
        parmBsnLvl = generateDDSParVal(parmBsnLvl, 
                                        probVal, 
                                        ctrlSetting["perturbFactor"])
        # After modifying parameter values in file, 
        modifyParInFileBsn(parmBsnLvl, 
                        parmBsnLvlFExtLst,
                        fdWorkingDir)
    
    modifyEndTime = datetime.datetime.now(timeZone)
    modifyTimeTotal = modifyEndTime - modifyStartTime
    print(".....Time for modifying parameter values: {}; Total Time: {}.....". format(
            modifyTimeTotal,
            modifyEndTime - startTime))
        
    # After modifying, run the SWAT model
    if iRunSWAT:
        runStartTime = datetime.datetime.now(timeZone)
        runCode = runSWATModel(get_osplatform(), fdWorkingDir, fdMain, fdDMPOTpyFiles, fnSwatExe)
        runEndTime = datetime.datetime.now(timeZone)
        runTimeTotal = runEndTime - modifyEndTime
        print(".....Time for running SWAT: {}; Total Time: {}.....". format(
                runTimeTotal,
                runEndTime - startTime))  
    else:
        runEndTime = datetime.datetime.now(timeZone)

    # Calculate Statistics 
    subObfTestDict, obfValNonOther, subAllStats = calObjFuncValues(iPrintForCio, 
                rcvRchLst, 
                obsDataLst, 
                parmObjFnKeys,
                subObfTestDict,
                fdWorkingDir, False, nCalVal, ctrlSetting["groupSubareaIdx"])

    totalTimeThisrun = runEndTime - startTime
    # Update the best parameter values and objective functions
    # based on the objective Functions
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
print(datetime.datetime.now(timeZone))
print("--------------------------------------")
print("Congratulations!!! It's done nicely~~~")
print("--------------------------------------")

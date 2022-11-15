# -*- coding: utf-8 -*-
"""
Created on Tue May 12 2020 

This class is designed to be a collection of functions dealing with 
DDS procedure.

@author: Qingyu.Feng
"""

##########################################################################
# Import modules #########################################################
##########################################################################
import json
import os, sys
import math
import numpy, glob
import subprocess
import pandas
import datetime
import fortranformat as ff
import shutil

from .globVars import *
from pyscripts.GRAPHUtil import *
from pyscripts.PLOTUtil import *

if useRay: 
    from .RaySWATUtil import *
else:
    from .SWATUtil import *

# Set up the random seed
numpy.random.seed(1)

# Turn the future warning of pandas off
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

##########################################################################
# Define functions #######################################################
##########################################################################
##########################################################################

##########################################################################
def getUsrBestParmSet(fnBsnPara, runNoBestPar, subParGroups, 
                        subParaFns, subNoGroups, parmSubLvl):

    calParmBsnDf = pandas.read_csv(fnBsnPara, sep = ",")
    calParmBsnDf = calParmBsnDf.set_index("RunNO")

    usrBestParmBsnDf = calParmBsnDf.loc[runNoBestPar,]

    usrBestParGroups = {}
    for subGI in subNoGroups.keys():
        subParGroups[subGI] = copy.deepcopy(parmSubLvl)
        calParmDf = pandas.read_csv(subParaFns[subGI], sep = ",")
        calParmDf = calParmDf.set_index("RunNO")
        usrBestParGroups[subGI] = calParmDf.loc[runNoBestPar,]

    return usrBestParmBsnDf, usrBestParGroups





##########################################################################
def getSubGroupsRchList(ctrlSetting):

    # Get reach connection information from shapefile
    if not os.path.isfile(fnReachShp):
        print("Reach shapefile was not found, please check and \
            put it under the 03gisLayers folder")
        return
    
    print(".....Reading reach Shapefile.....")
    field_names, subStrAtt = readShapeAttributes(fnReachShp)

    # Get a graph of the watershed 
    print(".....Generating watershed graphs.....")
    wsGraph, rcvRchLst = getGraph(field_names, subStrAtt)

    if (ctrlSetting["groupSubareaIdx"] == 1):
        print(".....Grouping subareas based on outlet.....")
        # Get Groups of the watershed
        subNoGroupsOrig = {}
        subNoGroupsOrig = groupSubForOutlet(ctrlSetting["outLetList"], wsGraph, rcvRchLst)
        subNoGroups = dealWithOverlayInSubGroups(subNoGroupsOrig)
        # Initiate the objective function values and out files
        parmObjFnKeys = list(subNoGroups.keys())

    elif (ctrlSetting["groupSubareaIdx"] == 0):
        subNoGroups = {"NotGrouping": rcvRchLst}
        # Initiate the objective function values and out files
        parmObjFnKeys = list(map(str, ctrlSetting["outLetList"]))
    
    return subNoGroups, parmObjFnKeys, rcvRchLst


##########################################################################
def initOutFileParmObjSublvl(ctrlSetting, parmSubLvl, parmObjFnKeys, rcvRchLst):

    subParaFns = {}
    subObjFunFns = {}

    # Added April 19, 2021 for recording parameter selection
    subParaSel01Fns = {}

    # Headers of obf value and param value files.
    obfValueHeader = "RunNO,Outlet_Var_Freq_Stat_Weight,PBIAS,NSE,RMSE,R2,MSE,TestOF,BestOF,probVal,TimeThisRun\n"
    paraValHeader = "RunNO,Outlet_Var_Freq_Stat_Weight," + ",".join(parmSubLvl["Symbol"].to_list()) + "\n"

    if (ctrlSetting["groupSubareaIdx"] == 1):
    
        # Files recording the autocalibration processes
        for opKeys in parmObjFnKeys:
            # Initialize the parameter files
            fnParaEachRun = os.path.join(fdOutputs,
                "DMPOT_Para_{}.out".format(opKeys))
            if os.path.isfile(fnParaEachRun):
                os.remove(fnParaEachRun)
            subParaFns[opKeys] = fnParaEachRun
            with open(subParaFns[opKeys], 'w') as parvFile:
                parvFile.writelines(paraValHeader)

            # Initialize the parameter selection files
            fnParaSelEachRun = os.path.join(fdOutputs,
                "DMPOT_ParaSel_{}.out".format(opKeys))
            if os.path.isfile(fnParaSelEachRun):
                os.remove(fnParaSelEachRun)
            subParaSel01Fns[opKeys] = fnParaSelEachRun
            with open(subParaSel01Fns[opKeys], 'w') as parselFile:
                parselFile.writelines(paraValHeader)


            # Initialize the objective value files
            fnObjFunEachRun =  os.path.join(fdOutputs,
                "DMPOT_ObjFun{}.out".format(opKeys))
            if os.path.isfile(fnObjFunEachRun):
                os.remove(fnObjFunEachRun)
            subObjFunFns[opKeys] = fnObjFunEachRun
            with open(subObjFunFns[opKeys], 'w') as obfFile:
                obfFile.writelines(obfValueHeader)

    elif (ctrlSetting["groupSubareaIdx"] == 0):
        subNoGroups = {"NotGrouping": rcvRchLst}
        # Initiate the objective function values and out files
        parmObjFnKeys = list(map(str, ctrlSetting["outLetList"]))
        
        # Only one subparam file need to be initialized
        fnParaEachRun = os.path.join(fdOutputs,
                "DMPOT_Para_NotGrouping.out")
        if os.path.isfile(fnParaEachRun):
            os.remove(fnParaEachRun)
        subParaFns["NotGrouping"] = fnParaEachRun
        with open(subParaFns["NotGrouping"], 'w') as parvFile:
            parvFile.writelines(paraValHeader)

        # Only one subparam file need to be initialized
        fnParaSelEachRun = os.path.join(fdOutputs,
                "DMPOT_ParaSel_NotGrouping.out")
        if os.path.isfile(fnParaSelEachRun):
            os.remove(fnParaSelEachRun)
        subParaSel01Fns["NotGrouping"] = fnParaSelEachRun
        with open(subParaSel01Fns["NotGrouping"], 'w') as parsFile:
            parsFile.writelines(paraValHeader)


        # Files recording the autocalibration processes
        for opKeys in parmObjFnKeys:
            # Initialize the objective value files
            fnObjFunEachRun =  os.path.join(fdOutputs,
                "DMPOT_ObjFun{}.out".format(opKeys))
            if os.path.isfile(fnObjFunEachRun):
                os.remove(fnObjFunEachRun)
            subObjFunFns[opKeys] = fnObjFunEachRun
            with open(subObjFunFns[opKeys], 'w') as obfFile:
                obfFile.writelines(obfValueHeader)

    return subParaFns, subObjFunFns, subParaSel01Fns


##########################################################################
def initOutFileParmObjSublvlCal(ctrlSetting, parmSubLvl, 
                                parmObjFnKeys, rcvRchLst, nCalVal):

    # For reading
    subParaFns = {}
    # Added April 19, 2021 for recording parameter selection
    subParaSel01Fns = {}

    # For recording
    subParaCalFns = {}
    subObjFunCalFns = {}

    # Headers of obf value and param value files.
    obfValueHeader = "RunNO,Outlet_Var_Freq_Stat_Weight,PBIAS,NSE,RMSE,R2,MSE,TestOF,BestOF,probVal,TimeThisRun\n"
    paraValHeader = "RunNO,Outlet_Var_Freq_Stat_Weight," + ",".join(parmSubLvl["Symbol"].to_list()) + "\n"

    if (ctrlSetting["groupSubareaIdx"] == 1):
        
        # Files recording the autocalibration processes
        for opKeys in parmObjFnKeys:
            # Initialize the parameter files
            fnParaEachRun = os.path.join(fdOutputs,
                "DMPOT_Para_{}.out".format(opKeys))
            subParaFns[opKeys] = fnParaEachRun

            fnParaEachRunCal = os.path.join(fdOutputs,
            "DMPOT_Para_{}_{}.out".format(opKeys, nCalVal))
            if os.path.isfile(fnParaEachRunCal):
                os.remove(fnParaEachRunCal)
            subParaCalFns[opKeys] = fnParaEachRunCal
            with open(subParaCalFns[opKeys], 'w') as parvFile:
                parvFile.writelines(paraValHeader)

            fnObjFunEachRunCal =  os.path.join(fdOutputs,
            "DMPOT_ObjFun{}_{}.out".format(opKeys, nCalVal))
            if os.path.isfile(fnObjFunEachRunCal):
                os.remove(fnObjFunEachRunCal)
            subObjFunCalFns[opKeys] = fnObjFunEachRunCal
            with open(subObjFunCalFns[opKeys], 'w') as obfFile:
                obfFile.writelines(obfValueHeader)

            # Initialize the parameter selection files
            fnParaSelEachRun = os.path.join(fdOutputs,
                "DMPOT_ParaSel_{}.out".format(opKeys))
            subParaSel01Fns[opKeys] = fnParaSelEachRun


    elif (ctrlSetting["groupSubareaIdx"] == 0):
        subNoGroups = {"NotGrouping": rcvRchLst}
        # Initiate the objective function values and out files
        parmObjFnKeys = list(map(str, ctrlSetting["outLetList"]))
        
        # Only one subparam file need to be initialized
        fnParaEachRun = os.path.join(fdOutputs,
                "DMPOT_Para_NotGrouping.out")
        subParaFns["NotGrouping"] = fnParaEachRun

        fnParaEachRunCal = os.path.join(fdOutputs,
                "DMPOT_Para_NotGrouping_{}.out".format(nCalVal))
        if os.path.isfile(fnParaEachRunCal):
            os.remove(fnParaEachRunCal)
        subParaCalFns["NotGrouping"] = fnParaEachRunCal
        with open(subParaCalFns["NotGrouping"], 'w') as parvFile:
                    parvFile.writelines(paraValHeader)

        # Files recording the autocalibration processes
        for opKeys in parmObjFnKeys:
            # Initialize the objective value files
            fnObjFunEachRunCal =  os.path.join(fdOutputs,
            "DMPOT_ObjFun{}_{}.out".format(opKeys, nCalVal))
            if os.path.isfile(fnObjFunEachRunCal):
                os.remove(fnObjFunEachRunCal)
            subObjFunCalFns[opKeys] = fnObjFunEachRunCal
            with open(subObjFunCalFns[opKeys], 'w') as obfFile:
                obfFile.writelines(obfValueHeader)

        # Only one subparam file need to be initialized
        fnParaSelEachRun = os.path.join(fdOutputs,
                "DMPOT_ParaSel_NotGrouping.out")
        subParaSel01Fns["NotGrouping"] = fnParaSelEachRun

    return subParaFns, subParaCalFns, subObjFunCalFns, subParaSel01Fns



##########################################################################
def initOFValDict(parmObjFnKeys):

    subObfBestDict = {}
    subObfTestDict = {}
    # Files recording the autocalibration processes
    for opKeys in parmObjFnKeys:
        # Initialize the objective values for each group
        subObfBestDict[opKeys] = 10e2
        subObfTestDict[opKeys] = 10e2
            
    bsnObfBest = 10e2

    return subObfBestDict,subObfTestDict, bsnObfBest


##########################################################################
def initOutFileParmObjBsnlvl(parmBsnLvl):

    # Initialize the basin level parameter file as a record
    # First Deal with basin level parameter
    fnBsnPara = os.path.join(fdOutputs, "DMPOT_Para_Bsn.out")
    if os.path.isfile(fnBsnPara):
        os.remove(fnBsnPara)
    paraValBSNHeader = "RunNO," + ",".join(parmBsnLvl["Symbol"].to_list()) + "\n"
    with open(fnBsnPara, 'w') as bsnParmFile:
        bsnParmFile.writelines(paraValBSNHeader)


    fnBsnParaSel = os.path.join(fdOutputs, "DMPOT_ParaSel_Bsn.out")
    if os.path.isfile(fnBsnParaSel):
        os.remove(fnBsnParaSel)
    with open(fnBsnParaSel, 'w') as bsnParmSFile:
        bsnParmSFile.writelines(paraValBSNHeader)

    return fnBsnPara, fnBsnParaSel



##########################################################################
def initOutFileParmObjBsnlvlCal(parmBsnLvl, nCalVal):

    # Initialize the basin level parameter file as a record
    # First Deal with basin level parameter
    fnBsnPara = os.path.join(fdOutputs, "DMPOT_Para_Bsn.out")
    fnBsnParaCal = os.path.join(fdOutputs,
        "DMPOT_Para_Bsn_{}.out".format(nCalVal))
    # Commented by Qingyu Feng 20211115
    # This was commented because the old file was overwritten when 
    # applyBest Parameter function is conducted.
    # if os.path.isfile(fnBsnParaCal):
    #     os.remove(fnBsnParaCal)
    # paraValBSNHeader = "RunNO," + ",".join(parmBsnLvl["Symbol"].to_list()) + "\n"
    # Commented by Qingyu Feng 20211115
    # with open(fnBsnParaCal, 'w') as bsnParmFile:
    #     bsnParmFile.writelines(paraValBSNHeader)

    fnBsnParaSel = os.path.join(fdOutputs, "DMPOT_ParaSel_Bsn.out")
    # Commented by Qingyu Feng 20211115
    # if os.path.isfile(fnBsnParaSel):
    #     os.remove(fnBsnParaSel)
    # with open(fnBsnParaSel, 'w') as bsnParmSFile:
    #     bsnParmSFile.writelines(paraValBSNHeader)


    return fnBsnPara, fnBsnParaCal, fnBsnParaSel


##########################################################################
def initParmInFilenameSubLvl(subNoGroups, parmSubLvl, runningDir):

    # Initialize the parameter set and values of obj func for each subarea group
    subParGroups = {}
    # Create the subarea and hru file names
    swatSubFnGroups = {}
    swatHruFnGroups = {}

    for subGI in subNoGroups.keys():
        # Initialize the parameter set and values of obj func for each subarea group
        subParGroups[subGI] = copy.deepcopy(parmSubLvl)

        # Geneate sub level file name list for modifying
        fnSWATSubFlLst = []
        fnSWATSubFlLst = list(map(buildSWATSubFn, subNoGroups[subGI]))
        swatSubFnGroups[subGI] = fnSWATSubFlLst

        # Generate hru level file name list for modifying
        fnSWATHruFlLst = []
        for subidx in fnSWATSubFlLst:
            fnSWATHruFlLst = fnSWATHruFlLst + buildSWATHruFn(subidx, runningDir)
        swatHruFnGroups[subGI] = fnSWATHruFlLst

    return subParGroups, swatSubFnGroups, swatHruFnGroups



##########################################################################
def copySWATTxtIO(timeZone, startTime, srcDr, DestDr):
    print(".....Copy TxtInOut contents to working directory.....")

    copyTioStartTime = datetime.datetime.now(timeZone)
    
    flTioSource = glob.glob("{}/*".format(srcDr))
    if useRay:
        copyStatus = []
        copyStatus = [copySWATFileToWDRay.remote(fnSWAT, DestDr)
                            for fnSWAT in flTioSource]
        ray.get(copyStatus)
    else:
        for fnSWAT in flTioSource:
            copySWATFileToWD(fnSWAT, DestDr)

    copyTioEndTime = datetime.datetime.now(timeZone)
    print(".....Time for copy: {}; Total Time: {}.....". format(
        copyTioEndTime - copyTioStartTime,
        copyTioEndTime-startTime))

##########################################################################
def generateDDSParVal(parDF, probVal, perturbFactor):
    """
    This function update parameter values with DDS procedure.
    This is conducted individually for each group
    """
    
    # Counter of how many DV seleceted for perturbation
    dvn_select = 0.0 
    # Initialize a parameter value for updating
    parUpdated = 0.0
    # Initialize a random variable for updating
    # Uniformly distributed random number for para selection
    uniRand = 0.0001

    # Added by Qingyu Feng April 6, 2021
    # One very important part is to replace the Test Value
    # with the best value. The logic here is we allways want to
    # try based on the best parameter values. If not, 
    # DDS will try always like the first time and without improvement.
    parDF["TestVal"] = parDF["BestVal"]
    # Reset the modThisRun to 0 everytime a new modification is going
    # to be made since we are recording the status of each parameter
    # in each run.
    parDF["ModThisRun"] = [0] * len(parDF["TestVal"])
    # Added by Qingyu Feng April 6, 2021

    for parIdx in parDF.index:
        uniRand = numpy.random.uniform(0.0, 1.0, 1)
        # print("Random value for uniRand: {}".format(uniRand))
        if (uniRand < probVal): 
            dvn_select = dvn_select + 1
            parUpdated = updateParValDDS(parDF.loc[parIdx,:],
                perturbFactor)
            parDF.loc[parIdx,"TestVal"] = parUpdated
            parDF.loc[parIdx,"ModThisRun"] = 1

    # After updating, deal with a special case for each group
    # This special case is very important. This is because
    # after about 1/3 of total runs, the program will fall
    # in local cycle and can not update the value. Thus, we pick
    # up one paramter randomly to be updated.
    if dvn_select == 0:
        
        uniRandIdx = numpy.random.randint(1, len(parDF.index), 1)
        # print("Random value for uniRandIdx: {}".format(uniRandIdx))
        parUpdated = updateParValDDS(parDF.loc[parDF.index[uniRandIdx],:],
                perturbFactor)
        # Update the value in dataframe
        parDF.loc[parDF.index[uniRandIdx],"TestVal"] = parUpdated
        parDF.loc[parDF.index[uniRandIdx],"ModThisRun"] = 1
        # print("dvn == 0, dealing with special case: updating parameter {} at index{}".format(
        #     parDF.loc[parDF.index[uniRandIdx],"Symbol"],  parDF.index[uniRandIdx]
        # ))

    return parDF
    
##########################################################################
def writeBsnParm(parmBsnLvl,
                runIdx,
                fnBsnPara,
                fnBsnParaSel):
    """
    This function update the best parameter values based on objective
    function value and write the test parameter values into a file.
    """
    # Also write the value of basin level parameters to a 
    # corresponding file.
    subParValLstBsn = ["{:.3f}".format(parVl) for parVl in parmBsnLvl["TestVal"]]
    lfwParaBsn = "{},".format(runIdx) + ",".join(subParValLstBsn) + "\n"
    with open(fnBsnPara, 'a') as parmFileBsn:
        parmFileBsn.writelines(lfwParaBsn) 


    subParSelLstBsn = ["{:.3f}".format(parVl) for parVl in parmBsnLvl["ModThisRun"]]
    lfwParaSBsn = "{},".format(runIdx) + ",".join(subParSelLstBsn) + "\n"
    with open(fnBsnParaSel, 'a') as parmSFileBsn:
        parmSFileBsn.writelines(lfwParaSBsn) 

##########################################################################
def updateBestParm(ctrlSetting, 
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
                    ):
    """
    This function update the best parameter values based on objective
    function value and write the test parameter values into a file.
    """
    # After logging and plotting, update the best and current parameter function,
    # update the parameter values to the best value or test value.  
    # Then, compare the stat between best and test of each group while grouping,
    # based on which, the value of parameters is used.

    # Deal with best value for other basin
    # Using sum of objective functions
    bsnObfTest = sum(obfValNonOther)
    # Using Average
    # bsnObfTest = sum(obfValNonOther)/len(obfValNonOther)
    # subObfTestDict["Other"] = bsnObfTest

    # Deal with printing and parameter updating based on the comparison between the overall
    # performance of this run. At the same time, the parameters and statistics need to be recorded.
    # For the grouping strategy:
    if (ctrlSetting["groupSubareaIdx"] == 1):
        # Loop over each strategy
        for subGPKey, subGPL in subParGroups.items():
            # First calculate the statistics with outlets
            if not subGPKey == "Other":
                # Print the values before updating
                print("Current and best objective function value for outlet {} are {:.3f}, {:.3f}, respectively!".format(
                        subGPKey, subObfTestDict[subGPKey], subObfBestDict[subGPKey]))
                
                # Update the best statistics based on overall performance:
                # if bsnObfTest < bsnObfBest:
                #     # First update the statistics with outlets
                #     if subObfTestDict[subGPKey] < subObfBestDict[subGPKey]:
                #         subObfBestDict[subGPKey] = subObfTestDict[subGPKey]
                #     # Then update the subarea level and basin level parameters
                #     # When overall performance is good, the best parameter is kept.
                #     # even those not improved is not good.
                #     subGPL["BestVal"] = subGPL["TestVal"]
                # # When the overall performance is not good, the improvement in individual
                # # outlets with improvements will need to be encouraged.
                # elif bsnObfTest >= bsnObfBest:
                #     # First update the statistics with outlets
                if subObfTestDict[subGPKey] < subObfBestDict[subGPKey]:
                    subObfBestDict[subGPKey] = subObfTestDict[subGPKey]
                    # Then update the subarea level and basin level parameters
                    subGPL["BestVal"] = subGPL["TestVal"]

                # After updating the parameter values, write them into the file
                # Write the parameter values and objective functions into 
                # corresponding files for recording. 
                # As a record, all parameters tried will be recorded
                subParValLst = ["{:.3f}".format(parVl) for parVl in subParGroups[subGPKey]["TestVal"]]
                lfwParaSub = "{},{},".format(runIdx, subGPKey) + ",".join(subParValLst) + "\n"
                with open(subParaFns[subGPKey], 'a') as parmFile:
                    parmFile.writelines(lfwParaSub) 

                subParSelLst = ["{}".format(parVl) for parVl in subParGroups[subGPKey]["ModThisRun"]]
                lfwParaSelSub = "{},{},".format(runIdx, subGPKey) + ",".join(subParSelLst) + "\n"
                with open(subParaSel01Fns[subGPKey], 'a') as parmsFile:
                    parmsFile.writelines(lfwParaSelSub) 

                # Need to write these variables into a file 
                #lfwAllStat = "RunNo\t{}\tOutLetVar\t{}\tPBIAS\t{}\tNSE\t{}\tRMSE\t{}\tR2\t{}\tMSE\t{}\n"
                for statK, statV in subAllStats[subGPKey].items():
                    lfwAllStat = "{},{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{},{}\n".format(
                        runIdx, statK,
                        statV["PBIAS"], statV["NSE"], statV["RMSE"], 
                        statV["R2"], statV["MSE"], 
                        subObfTestDict[subGPKey], subObfBestDict[subGPKey],
                        probVal, totalTimeThisrun 
                    )
                    with open(subObjFunFns[subGPKey], 'a') as obfFile:
                        obfFile.writelines(lfwAllStat)

            elif subGPKey == "Other":
                subObfTestDict["Other"] = bsnObfTest
                # Update the best statistics based on overall performance:
                if bsnObfTest < bsnObfBest:
                    subObfBestDict["Other"] = subObfTestDict["Other"]
                    subParGroups["Other"]["BestVal"] = subParGroups["Other"]["TestVal"]

                # Record the test values, status, and statistics
                subParValLst = ["{:.3f}".format(parVl) for parVl in subParGroups["Other"]["TestVal"]]
                lfwParaSub = "{},{},".format(runIdx, "Other") + ",".join(subParValLst) + "\n"
                with open(subParaFns["Other"], 'a') as parmFile:
                    parmFile.writelines(lfwParaSub) 

                subParSelLst = ["{}".format(parVl) for parVl in subParGroups["Other"]["ModThisRun"]]
                lfwParaSelSub = "{},{},".format(runIdx, "Other") + ",".join(subParSelLst) + "\n"
                with open(subParaSel01Fns["Other"], 'a') as parmsFile:
                    parmsFile.writelines(lfwParaSelSub) 

                # Write the objective functions for other groups
                lfwAllStat = "{},Others,{:.3f},{:.3f},{},{}\n".format(
                        runIdx,
                        subObfTestDict["Other"], subObfBestDict["Other"],
                        probVal, totalTimeThisrun 
                    )
                with open(subObjFunFns["Other"], 'a') as obfFile:
                        obfFile.writelines(lfwAllStat)
    
    # Dealing with the non grouping strategies
    elif (ctrlSetting["groupSubareaIdx"] == 0):
        # First calculate the mean of best and test of objective functions
        for subGPKey in parmObjFnKeys:
            print("Current and best objective function value for outlet {} are {:.3f}, {:.3f}, respectively!".format(
                    subGPKey, subObfTestDict[subGPKey], subObfBestDict[subGPKey]))

            # First calculate the statistics with outlets
            # Update the best statistics based on overall performance:
            if bsnObfTest < bsnObfBest:
                if subObfTestDict[subGPKey] < subObfBestDict[subGPKey]:
                    subObfBestDict[subGPKey] = subObfTestDict[subGPKey]
            
            # Need to write these variables into a file 
            #lfwAllStat = "RunNo\t{}\tOutLetVar\t{}\tPBIAS\t{}\tNSE\t{}\tRMSE\t{}\tR2\t{}\tMSE\t{}\n"
            for statK, statV in subAllStats[subGPKey].items():
                lfwAllStat = "{},{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{},{}\n".format(
                        runIdx, statK,
                        statV["PBIAS"], statV["NSE"], statV["RMSE"], 
                        statV["R2"], statV["MSE"], 
                        subObfTestDict[subGPKey], subObfBestDict[subGPKey],
                        probVal, totalTimeThisrun 
                    )
                with open(subObjFunFns[subGPKey], 'a') as obfFile:
                    obfFile.writelines(lfwAllStat)

        # Update the parameter values for recording
        if bsnObfTest < bsnObfBest:
            subParGroups["NotGrouping"]["BestVal"] = subParGroups["NotGrouping"]["TestVal"]
        
        # Write the parameter values and objective functions into 
        # corresponding files for recording.
        subParValLst = ["{:.3f}".format(parVl) for parVl in subParGroups["NotGrouping"]["TestVal"]]
        lfwParaSub = "{},{},".format(runIdx, "9999") + ",".join(subParValLst) + "\n"
        with open(subParaFns["NotGrouping"], 'a') as parmFile:
            parmFile.writelines(lfwParaSub) 

        subParSelLst = ["{}".format(parVl) for parVl in subParGroups["NotGrouping"]["ModThisRun"]]
        lfwParaSelSub = "{},{},".format(runIdx, "9999") + ",".join(subParSelLst) + "\n"
        with open(subParaSel01Fns["NotGrouping"], 'a') as parmsFile:
            parmsFile.writelines(lfwParaSelSub) 

    # Deal with the basin level parameters
    if bsnObfTest < bsnObfBest:
        parmBsnLvl["BestVal"] = parmBsnLvl["TestVal"]

    # Also write the value of basin level parameters to a 
    # corresponding file.
    subParValLstBsn = ["{:.3f}".format(parVl) for parVl in parmBsnLvl["TestVal"]]
    lfwParaBsn = "{},".format(runIdx) + ",".join(subParValLstBsn) + "\n"
    with open(fnBsnPara, 'a') as parmFileBsn:
        parmFileBsn.writelines(lfwParaBsn) 

    subParSelLstBsn = ["{:.3f}".format(parVl) for parVl in parmBsnLvl["ModThisRun"]]
    lfwParaSBsn = "{},".format(runIdx) + ",".join(subParSelLstBsn) + "\n"
    with open(fnBsnParaSel, 'a') as parmSFileBsn:
        parmSFileBsn.writelines(lfwParaSBsn) 
            
    # The best statistics will be updated later since the comparison is needed in dealing
    # with updating parameters. If the best is updated here, the comparison between test and
    # best will not be valid since they will be equal to each other.
    # Update the best statistics:
    # The printing here is to check whether the obj value of this run is better than the
    # history best runs. So, they were printed before being updated.
    print("""Current and best sum values of objective functions are {:.3f}, {:.3f}, respectively!""".format(
                bsnObfTest, bsnObfBest))
    if bsnObfTest < bsnObfBest:
        bsnObfBest = bsnObfTest 

    #if (ctrlSetting["groupSubareaIdx"] == 1):
        # for subGPKey, subGPL in subParGroups.items():
            # Commented by Qingyu Feng on Oct 27, 2021
            # The old version below has one issue:
            # When overall (sum) value was not know, the best statistics values was updated.
            # Thus, it can not reflect the actual progress of updating.
            # Since sum of statistics was used to update parameter values, it
            # should also be used in updating statistic values.
            ##########################################
            #     if not subGPKey == "Other":
            #         print("Current and best objective function value for outlet {} are {:.3f}, {:.3f}, respectively!".format(
            #             subGPKey, subObfTestDict[subGPKey], subObfBestDict[subGPKey]))
                    
            #         # commented by Qingyu Feng 20211027
            #         # This block of code was updating the best parameter values based on the statistics
            #         # of incifidual outlet. This was slowing down the updating since some outlet
            #         # might not have improvements at all and will not be updated for many runs.
            #         # Instead, the overall performance was used to try and see whether it converges faster.
            #         # # Update the best parameter to test parameter if the new values are better
            #         # if subObfTestDict[subGPKey] < subObfBestDict[subGPKey]:
            #         #     subObfBestDict[subGPKey] = subObfTestDict[subGPKey]
            #         #     subGPL["BestVal"] = subGPL["TestVal"]
            #         # commented by Qingyu Feng 20211027

            #         # Write the parameter values and objective functions into 
            #         # corresponding files for recording. 
            #         # As a record, all parameters tried will be recorded
            #         subParValLst = ["{:.3f}".format(parVl) for parVl in subParGroups[subGPKey]["TestVal"]]
            #         lfwParaSub = "{},{},".format(runIdx, subGPKey) + ",".join(subParValLst) + "\n"
            #         with open(subParaFns[subGPKey], 'a') as parmFile:
            #             parmFile.writelines(lfwParaSub) 

            #         subParSelLst = ["{}".format(parVl) for parVl in subParGroups[subGPKey]["ModThisRun"]]
            #         lfwParaSelSub = "{},{},".format(runIdx, subGPKey) + ",".join(subParSelLst) + "\n"
            #         with open(subParaSel01Fns[subGPKey], 'a') as parmsFile:
            #             parmsFile.writelines(lfwParaSelSub) 

            #         # Need to write these statistics into a file 
            #         #lfwAllStat = "RunNo\t{}\tOutLetVar\t{}\tPBIAS\t{}\tNSE\t{}\tRMSE\t{}\tR2\t{}\tMSE\t{}\n"
            #         for statK, statV in subAllStats[subGPKey].items():
            #             lfwAllStat = "{},{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{},{}\n".format(
            #                 runIdx, statK,
            #                 statV["PBIAS"], statV["NSE"], statV["RMSE"], 
            #                 statV["R2"], statV["MSE"], 
            #                 subObfTestDict[subGPKey], subObfBestDict[subGPKey],
            #                 probVal, totalTimeThisrun 
            #             )
            #             with open(subObjFunFns[subGPKey], 'a') as obfFile:
            #                 obfFile.writelines(lfwAllStat)

            # # Deal with best value for other basin
            # # Using sum of objective functions
            # bsnObfTest = sum(obfValNonOther)
            # # Using Average
            # # bsnObfTest = sum(obfValNonOther)/len(obfValNonOther)
            # # subObfTestDict["Other"] = bsnObfTest
            # print("""Current and best sum values of objective functions are {:.3f}, {:.3f}, respectively!""".format(
            #             bsnObfTest, bsnObfBest))
            
            # if bsnObfTest < bsnObfBest:
            #     bsnObfBest = bsnObfTest 
            #     # print("Testing objective function is better than best objective function, updating parameter values to the test.")
            #     # Update the subarea level and basin level parameters
            #     subGPL["BestVal"] = subGPL["TestVal"]
            #     parmBsnLvl["BestVal"] = parmBsnLvl["TestVal"]
            #     # Also update the paramter value for others if they are in the group
            #     if "Other" in list(subParGroups.keys()):
            #         subObfBestDict["Other"] = subObfTestDict["Other"]
            #         subParGroups["Other"]["BestVal"] = subParGroups["Other"]["TestVal"]
            # # else:
            # #     print("Testing objective function is larger than best objective function, not updating parameter values to the test.")
            # # After updating the parameter values, write them into the file
            # # Write the parameter values and objective functions into 
            # # corresponding files for recording. 
            # # As a record, all parameters tried will be recorded
            # if "Other" in list(subParGroups.keys()):
            #     subParValLst = ["{:.3f}".format(parVl) for parVl in subParGroups["Other"]["TestVal"]]
            #     lfwParaSub = "{},{},".format(runIdx, "Other") + ",".join(subParValLst) + "\n"
            #     with open(subParaFns["Other"], 'a') as parmFile:
            #         parmFile.writelines(lfwParaSub) 

            #     subParSelLst = ["{}".format(parVl) for parVl in subParGroups["Other"]["ModThisRun"]]
            #     lfwParaSelSub = "{},{},".format(runIdx, "Other") + ",".join(subParSelLst) + "\n"
            #     with open(subParaSel01Fns["Other"], 'a') as parmsFile:
            #         parmsFile.writelines(lfwParaSelSub) 

            #     # Write the objective functions for other groups
            #     lfwAllStat = "{},Others,{:.3f},{:.3f},{},{}\n".format(
            #             runIdx,
            #             subObfTestDict["Other"], subObfBestDict["Other"],
            #             probVal, totalTimeThisrun 
            #         )
            #     with open(subObjFunFns["Other"], 'a') as obfFile:
            #             obfFile.writelines(lfwAllStat)
            ##########################################
            # Commented by Qingyu Feng on Oct 27, 2021

        # In this version, for the grouping strategy, the first thing to do is to judge whether the 
        # current sum of statistics is good or not for the best sum of statistics.
        # If so, update the best values. This will need to be done before printing.
        

    return subParGroups, subObfBestDict, bsnObfBest, parmBsnLvl


##########################################################################
def calObjFuncValues(iPrintForCio, 
                    rcvRchLst, 
                    obsDataLst, 
                    parmObjFnKeys,
                    subObfTestDict, 
                    runningDir, iGenLineChart, nCalVal, groupSubareaIdx):
    """
    This function modify parameter values in files at the basin level.
    """
    
    # There are still some preparation to be done
    # Initialize the objective function Value
    # This is a dictionary containing the 5 statistics for each outlet_variable
    # list. 
    subAllStats = {}

    fnRch = os.path.join(runningDir, "output.rch")
    try:
        rchDFWhole = getRch2DF(fnRch, iPrintForCio, len(rcvRchLst))
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnRch, e))
        exit(1)

    # Then construct series of observed and simulated pairs
    # for stat calculation
    obsSimPair = buildObsSimPair(obsDataLst,
                                rchDFWhole,
                                varIDObsHdrPair)
    
    obsSimPairKeys = list(obsSimPair.keys())
    obsSimPairKeys = [oltVF.split("_") for oltVF in obsSimPairKeys]
    
    obsSimPairKeysSplit = {}
    for obpk in obsSimPairKeys:
        obsSimPairKeysSplit[obpk[0]] = "_".join(obpk)
    # Then calculate statistics and objective function. 
    # This need to be done based on the sub groups.
    # In the subarea grouping function, there are two options:
    # 1. If the user selected to use subarea grouping, the subarea will 
    # be grouped based on their upstream-downstream relationships.
    # Under this situation, there might be some other subareas not 
    # included, which means observed values are not provided for these
    # areas. Then, the objective function can not be calculated to 
    # infer how the parameter should be updated.
    # Here, an assumption was made. The average value of the objective
    # function from user specified outlets will be used as the reference
    # for other groups.  
    # 2. If the user does not select to use subarea grouping, there will
    # be only one group and the average value of the objective function
    # for all outlets and variables will be used to infer how the 
    # parameter will be modified.
    # In the code, the objective function values for each outlet will
    # be calculated, and the "Other" group will be dealt as a special calse.

    # subObfBestDict = {} # 10e2 = 10 * 10^2
    # subObfTestDict = {} # 10e2 = 10 * 10^2
    # if (ctrlSetting["groupSubareaIdx"] == 1):
    # A var to store the value of statistics of all non other sub groups.
    # In order to calculate the average value of the average OBJ and 
    # provide reference for the other groups.
    obfValNonOther = []

    # Calculating objective functions and making plots
    for subGPKey in parmObjFnKeys:
        # First calculate the statistics with outlets
        if not subGPKey == "Other":
            subAllStats[subGPKey] = {}
            # Calculating objective function value for outlet
            subAllStats[subGPKey] = calAllStatEachOlt(obsSimPair, 
                                                    subGPKey)

            # After calculating the statistics, calculate the objective functions
            subObfTestDict[subGPKey] = calObjFunValue(subAllStats[subGPKey])
            obfValNonOther.append(subObfTestDict[subGPKey])

    # Generate Plots for the interested outlets
    # TODO: To stop making figures if the model is not run successfully
    # because there is no simulated data. Thus, the program will crash
    # for this run.
    if iGenLineChart:
        genPDFLinePlotObsSim(fdOutputs, 
                         parmObjFnKeys,
                         obsSimPairKeysSplit,
                         obsSimPair, nCalVal, groupSubareaIdx)
        # # Generate combined figures
        # fnpPlotPdf = "{}/ObsVsSim_UserBest.pdf".format(fdOutputs)
        # if os.path.isfile(fnpPlotPdf):
        #     os.remove(fnpPlotPdf)
        # outPlotPdf = PdfPages(fnpPlotPdf)
        
        # # Calculating objective functions and making plots
        # for subGPKey in parmObjFnKeys:
        #     # First calculate the statistics with outlets
        #     if not subGPKey == "Other":
        #         # Generateing plots for each outlet
        #         obsSimThisKey = obsSimPairKeysSplit[subGPKey]
        #         plt.plot(obsSimPair[obsSimThisKey][0], label="Observation")
        #         plt.plot(obsSimPair[obsSimThisKey][1], label="Simulation")
        #         plt.legend(loc="upper right",fontsize='x-large')
        #         plt.title("Station {} Calibration".format(subGPKey))
        #         outPlotPdf.savefig()
        #         plt.close('all')

        # outPlotPdf.close()

        # Generate single figure for each outlet
        genPNGLinePlotObsSim(fdOutputs, 
                         parmObjFnKeys,
                         obsSimPairKeysSplit,
                         obsSimPair, nCalVal, groupSubareaIdx)
        # # Calculating objective functions and making plots
        # for subGPKey in parmObjFnKeys:
        #     # First calculate the statistics with outlets
        #     if not subGPKey == "Other":
        #         # Generate combined figures
        #         fnpPlotPng = "{}/ObsVsSim_UserBest{}.png".format(fdOutputs, subGPKey)
        #         if os.path.isfile(fnpPlotPng):
        #             os.remove(fnpPlotPng)
        #         obsSimThisKey = obsSimPairKeysSplit[subGPKey]
        #         genFigSingleOlt(fnpPlotPng, obsSimThisKey, obsSimPair)

    return subObfTestDict, obfValNonOther, subAllStats



##########################################################################
def modifyParInFileBsn(parmBsnLvl, 
                    parmBsnLvlFExtLst,
                    runningDir):
    """
    This function modify parameter values in files at the basin level.
    """
    # They use different file names.
    for flExtBLvl in parmBsnLvlFExtLst:
        # Get the list of parameters in a certain file
        selParInFile = parmBsnLvl.loc[parmBsnLvl["File"] == flExtBLvl]
        
        if flExtBLvl == ".bsn":
            if useRay:
                updateParInBsnRay(selParInFile, runningDir)
            else:
                updateParInBsn(selParInFile, runningDir)

        if flExtBLvl == "crop.dat":
            if useRay:
                updateParInCropRay(selParInFile, runningDir, fdmodelTxtInOut)
            else:
                updateParInCrop(selParInFile, runningDir, fdmodelTxtInOut)

        if flExtBLvl == ".wwq":
            if useRay:
                updateParInWwqRay(selParInFile, runningDir)
            else:
                updateParInWwq(selParInFile, runningDir)


##########################################################################
def modifyParInFileSub(subGPL, 
                    parmSubLvlFExtLst, 
                    swatSubFnGroups, 
                    subGPKey, 
                    swatHruFnGroups,
                    runningDir
                    ):
    """
    This function modify parameter values in files.
    """
    # They use different file names.
    for flExtSLvl in parmSubLvlFExtLst:
        # Get the list of parameters in a certain file
        selParInFile = subGPL.loc[subGPL["File"] == flExtSLvl]
                
        # Update subarea level files
        # subLvlFlExtLst = [".sub", ".rte", ".swq"]
        if flExtSLvl == ".sub":
            if useRay:
                updateStatus = []
                updateStatus = [updateParInSubRay.remote(fnSwatSubLvl, selParInFile, runningDir)
                                    for fnSwatSubLvl in swatSubFnGroups[subGPKey]]
                ray.get(updateStatus)
            # Code to run the function in a for loop
            else:
                for fnSwatSubLvl in swatSubFnGroups[subGPKey]:
                    updateParInSub(fnSwatSubLvl, selParInFile, runningDir)

        elif flExtSLvl == ".rte":
            if useRay:
                updateStatus = []
                updateStatus = [updateParInRteRay.remote(fnSwatSubLvl, selParInFile,
                                        runningDir, fdmodelTxtInOut)
                                    for fnSwatSubLvl in swatSubFnGroups[subGPKey]]
                ray.get(updateStatus)
            # Code to run the function in a for loop
            else:
                for fnSwatSubLvl in swatSubFnGroups[subGPKey]:
                    updateParInRte(fnSwatSubLvl, selParInFile, runningDir, fdmodelTxtInOut)

        elif flExtSLvl == ".swq":
            if useRay:
                updateStatus = []
                updateStatus = [updateParInSwqRay.remote(fnSwatSubLvl, selParInFile, runningDir)
                                    for fnSwatSubLvl in swatSubFnGroups[subGPKey]]
                ray.get(updateStatus)
            # Code to run the function in a for loop
            else:
                for fnSwatSubLvl in swatSubFnGroups[subGPKey]:
                    updateParInSwq(fnSwatSubLvl, selParInFile, runningDir)

        # TODO: Add parameters for reservoir
        # elif flExtSLvl == ".res":
        #     if useRay:
        #         updateParInResRay(selParInFile, fdWorkingDir, subNosWithRes, swatSubFnGroups[subGPKey], fdmodelTxtInOut)
        #     else:
        #         updateParInRes(selParInFile, fdWorkingDir, subNosWithRes, swatSubFnGroups[subGPKey], fdmodelTxtInOut)

        # Start processing HRU level files
        # hruLvlFlExtLst = [".gw", ".hru", ".mgt", ".sol", ".chm"] 
        elif flExtSLvl == ".gw":
            if useRay:
                updateStatus = []
                updateStatus = [updateParInGwRay.remote(fnSwatSubLvl, selParInFile,
                                    runningDir, fdmodelTxtInOut)
                                    for fnSwatSubLvl in swatHruFnGroups[subGPKey]]
                ray.get(updateStatus)
            # Code to run the function in a for loop
            else:
                for fnSwatSubLvl in swatHruFnGroups[subGPKey]:
                    updateParInGw(fnSwatSubLvl, selParInFile, runningDir, fdmodelTxtInOut)
        
        elif flExtSLvl == ".hru":
            if useRay:
                updateStatus = []
                updateStatus = [updateParInHruRay.remote(fnSwatSubLvl, selParInFile, 
                                    runningDir, fdmodelTxtInOut)
                                    for fnSwatSubLvl in swatHruFnGroups[subGPKey]]
                ray.get(updateStatus)
            # Code to run the function in a for loop
            else:
                for fnSwatSubLvl in swatHruFnGroups[subGPKey]:
                    updateParInHru(fnSwatSubLvl, selParInFile, runningDir, fdmodelTxtInOut)

        elif flExtSLvl == ".mgt":
            if useRay:
                updateStatus = []
                updateStatus = [updateParInMgtRay.remote(fnSwatSubLvl, selParInFile, 
                                    runningDir, rowCropLst, fdmodelTxtInOut)
                                    for fnSwatSubLvl in swatHruFnGroups[subGPKey]]
                ray.get(updateStatus)
            # Code to run the function in a for loop
            else:
                for fnSwatSubLvl in swatHruFnGroups[subGPKey]:
                    updateParInMgt(fnSwatSubLvl, selParInFile, 
                    runningDir, rowCropLst, fdmodelTxtInOut)

        elif flExtSLvl == ".chm":
            if useRay:
                updateStatus = []
                updateStatus = [updateParInChmRay.remote(fnSwatSubLvl, selParInFile, 
                                    runningDir)
                                    for fnSwatSubLvl in swatHruFnGroups[subGPKey]]
                ray.get(updateStatus)
            # Code to run the function in a for loop
            else:
                for fnSwatSubLvl in swatHruFnGroups[subGPKey]:
                    updateParInChm(fnSwatSubLvl, selParInFile, 
                    runningDir)

        elif flExtSLvl == ".sol":
            if useRay:
                updateStatus = []
                updateStatus = [updateParInSolRay.remote(fnSwatSubLvl, selParInFile, 
                                    runningDir, fdmodelTxtInOut)
                                    for fnSwatSubLvl in swatHruFnGroups[subGPKey]]
                ray.get(updateStatus)
            # Code to run the function in a for loop
            else:
                for fnSwatSubLvl in swatHruFnGroups[subGPKey]:
                    updateParInSol(fnSwatSubLvl, selParInFile, 
                    runningDir, fdmodelTxtInOut)

        # End of for loop updating subarea and hru level parameters parameters            





##########################################################################
def generateRandomParVal(parDF, ranNum):
    """
    This function calculates the random values for initial runs.
    """
    # print("Random value for ranNum in updating parameters:　{}".format(ranNum))
    # ranNum = [random.random() for i in range(len(parDF.index))]
    parDF["ModThisRun"] = [0] * len(parDF["TestVal"])
    parDF["TestVal"] = parDF["LowerBound"] + ranNum[0] * (parDF["UpperBound"] - parDF["LowerBound"])
    parDF["ModThisRun"] = [1] * len(parDF["TestVal"])

    return parDF
    




##########################################################################
def calObjFunValue(subAllStatsOneOlt):
    """
    This function calculates the objective function value for each outlet
    based on the user specified and weights.
    """
    # StatDict is a dictionary containing 5 statistics
    # statDict["PBIAS"] = PBIAS
    # statDict["NSE"] = NSE
    # statDict["RMSE"] = RMSE
    # statDict["R2"] = R2
    # statDict["MSE"] = MSE
    # Only one stat is used in calculating the objective 
    # funcitons.
    # '1=1-NS', '2=BIAS','3=RMSE','4=R2','5=MSE'

    statKey = list(subAllStatsOneOlt.keys())[0]
    statVals = subAllStatsOneOlt[statKey]
    keySplit = statKey.split("_")
    statIdx = keySplit[3]
    statWgt = keySplit[4]
    
    # Determine the stat value used as objective function
    if statIdx == "1":
        statVal = 1 - statVals["NSE"]
    elif statIdx == "2":
        statVal = statVals["PBIAS"]
    elif statIdx == "3":
        statVal = statVals["RMSE"]
    elif statIdx == "4":
        statVal = 1 - statVals["R2"]
    elif statIdx == "5":
        statVal = statVals["MSE"]

    return statVal * float(statWgt)



##########################################################################
def calStatistics(obsTS, simTS):

    """
    PBIAS = 100*(meanQ_Out-meanQ_In)/meanQ_In
    NSE = 1 - (errorSum2/errorSumI2)
    RMSE = (errorSum2/TimeStatsTotal)**0.5            
    R2 = (errorSumR2**2)/(errorSumI2*errorSumO2) 
    MSE = errorSum2/TimeStatsTotal
        
    """
    # Notes for the change of variable names
    # between IPEAT and dmpot
    # meanQ_Out = simMean
    # meanQ_IN = obsMean
    

    obsTSMean = sum(obsTS)/len(obsTS)
    if len(simTS) == 0:
        print("Your SWAT might not run successfully this time, trying next run")
        simTSMean = 9999
        PBIAS = 9999 
        NSE = -9998 
        RMSE = 9999 
        R2 = -9998 
        MSE = 9999 

    else:
        simTS = list(map(float, simTS))
        simTSMean = sum(simTS)/len(simTS)

        # Sum of errors above
        sumErrObSm = 0.0
        sumSqErrObSm = 0.0
        sumSqErrObs = 0.0
        sumSqErrSim = 0.0
        sumProdErrObSm = 0.0

        PBIAS = 999.99
        NSE = 999.99
        RMSE = 999.99           
        R2 = 999.99 
        MSE = 999.99

        for tsidx in range(len(obsTS)):
            # IPEATPlus Error: error between OBS(i) and SIM(i)
            errObSm = 0.0
            # IPEATPlus error2: square of error between OBS(i) and SIM(i)
            sqErrObSm = 0.0
            # IPEATPlus errorI2: square of the error between OBS(i) and OBS(mean)
            sqErrObs = 0.0
            # IPEATPlus errorO2: square of the error between SIM(i) and SIM(mean)
            sqErrSim = 0.0
            # IPEATPlus errorR2: product of OBS(i) error -mean and SIM(i) error - mean
            prodErrObSm = 0.0

            errObSm = simTS[tsidx] - obsTS[tsidx]
            sqErrObSm = errObSm**2
            sqErrObs = (obsTS[tsidx] - obsTSMean)**2
            sqErrSim = (simTS[tsidx] - simTSMean)**2
            prodErrObSm = (obsTS[tsidx] - obsTSMean) * (simTS[tsidx] - simTSMean)

            sumErrObSm = sumErrObSm + errObSm
            sumSqErrObSm = sumSqErrObSm + sqErrObSm
            sumSqErrObs = sumSqErrObs + sqErrObs
            sumSqErrSim = sumSqErrSim + sqErrSim
            sumProdErrObSm = sumProdErrObSm + prodErrObSm

        # Added by Qingyu Feng Mar 26, 2021
        # This modification was done to prevent potential errors of 
        # division by zero. It was found that simulated flow were all 0s
        # and this will cause error. It is very rare but still happened.
        # Calculate R2
        if ((sumSqErrObs == 0.0) or (sumSqErrSim == 0.0)):
            R2 = -9998
        else:
            R2 = (sumProdErrObSm**2)/(sumSqErrObs*sumSqErrSim) 

        # Calculate PBIAS = 100 * (sum(sim-obs)/sum(obs))
        if obsTSMean == 0.0:
            PBIAS = 9999
        else:
            PBIAS = 100*(simTSMean-obsTSMean)/obsTSMean

        if sumSqErrObs == 0.0:
            NSE = -9998
        else:
            NSE = 1 - (sumSqErrObSm/sumSqErrObs)
        
        if len(obsTS) == 0:
            RMSE = 999.99
            MSE = 9999
        else:
            RMSE = (sumSqErrObSm/len(obsTS))**0.5  
            MSE = sumProdErrObSm/len(obsTS)
        # Added by Qingyu Feng Mar 26, 2021
        
        # Restrict the values for printing
        if PBIAS > 9999: 
            PBIAS = 9999
        elif PBIAS < -999:
            PBIAS = 9999

        if NSE > 9999: 
            NSE = -9998
        elif NSE < -999:
            NSE = -9998

        if RMSE > 9999: 
            RMSE = 9999
        elif RMSE < -999:
            RMSE = 9999

        if R2 > 9999: 
            R2 = -9998
        elif R2 < -999:
            R2 = -9998
            
        if MSE > 9999: 
            MSE = 9999
        elif MSE < -999:
            MSE = 9999

    statDict = {}
    statDict["PBIAS"] = PBIAS
    statDict["NSE"] = NSE
    statDict["RMSE"] = RMSE
    statDict["R2"] = R2
    statDict["MSE"] = MSE

    # print("Stat self cal")
    # for key, val in statDict.items():
    #     print(key, val)

    # print("Stat numpy")
    # obsTSAry = numpy.array(obsTS)
    # simTSAry = numpy.array(simTS)
    # difObsSim = obsTSAry - simTSAry
    # numpyRE = (numpy.mean(obsTSAry) - numpy.mean(obsTSAry))/numpy.mean(obsTSAry) * 100
    # numpyBIAS = numpy.sum(difObsSim)/obsTSAry.size
    # numpySSE = numpy.sum(numpy.power(difObsSim, 2))
    # numpyRMSE = numpy.sqrt(numpySSE/obsTSAry.size)
    # if numpyRE == 100:
    #     numpyR2 = 0.0
    # else:
    #     numpyCorr1 = numpy.corrcoef(simTSAry, obsTSAry)
    #     numpyR = numpyCorr1[0, 1]
    #     numpyR2 = numpy.power(numpyR, 2)
    # # Compute Nash Sutcliff coefficient
    # numpyObsMean = numpy.mean(obsTSAry)
    # numpyErrorSQObsSim = numpy.power(difObsSim, 2)
    # numpyErrorSQObs = numpy.power((obsTSAry - numpyObsMean), 2)
    # numpyNS = 1 - numpy.sum(numpyErrorSQObsSim)/numpy.sum(numpyErrorSQObs)

    # statDictNumpy = {}
    # statDictNumpy["PBIASnumpy"] = numpyBIAS
    # statDictNumpy["NSEnumpy"] = numpyNS
    # statDictNumpy["RMSEnumpy"] = numpyRMSE
    # statDictNumpy["R2numpy"] = numpyR2
    # statDictNumpy["MSEnumpy"] = numpySSE

    # for key, val in statDictNumpy.items():
    #     print(key, val)

    return statDict




##########################################################################
def calAllStatEachOlt(obsSimPair, subGPKey):

    """
    This function calculate the average objective function values over
    the objective function for all variables across all time frequencies.
    Users might specify different objective function values for 
    different variables. This tool offer this feature.
    While calculating, the selected objective function for each 
    variable will be calculated.
    """
    # First construct a key to extract the corresponding obsSimPar DF
    oltVarFreqAll = list(obsSimPair.keys())

    oltVarFreqAll = [oltVF.split("_") for oltVF in oltVarFreqAll]
    
    oltVarFreq = []
    
    if subGPKey == "NotGrouping": 
        oltVarFreq = oltVarFreqAll
    else: 
        for oltVF in oltVarFreqAll:
            if (oltVF[0] == subGPKey):
                oltVarFreq.append(oltVF)

    outStat = {}

    for ovf in oltVarFreq:
        obsKey = "_".join(ovf)
        obsSimPairThis = obsSimPair[obsKey]
        # TODO: Remove missing values in the observed values.
        # This can be done after combining the obs and simulated 
        # series to avoid messing up the dates.

        # The variable headers are added _x and _y for 
        # observed and simulated columns in the dataframe.
        obsTS = obsSimPairThis[0]
        simTS = obsSimPairThis[1]
        
        # StatDict is a dictionary containing 5 statistics
        # statDict["PBIAS"] = PBIAS
        # statDict["NSE"] = NSE
        # statDict["RMSE"] = RMSE
        # statDict["R2"] = R2
        # statDict["MSE"] = MSE
        # Only one stat is used in calculating the objective 
        # funcitons.
        statDict = calStatistics(obsTS, simTS)
        
        outStat[obsKey] = statDict


    return outStat
        

##########################################################################
def getOltVarDict(outLetList, outputVarList):

    """
    This function find the groups of variables belonging to the same
    outlet. Here the frequency are not included.
    """

    outLetVarDict = {}
    
    for iOlt in range(len(outLetList)):
        if not outLetList[iOlt] in outLetVarDict.keys():
            outLetVarDict[outLetList[iOlt]] = []
        outLetVarDict[outLetList[iOlt]].append(outputVarList[iOlt])

    return outLetVarDict


##########################################################################
def buildObsSimPair(obsDict, simOutRchDF, varIDObsHdrPair):

    """
    This function read the data from simDF and add corresponding columns
    into the obsDict.
    """

    # # Construct daily, monthly and annual time series.
    # startSimDate = pandas.Timestamp("{}/{}/{}".format(
    #         ctrlSetting["startSimDate"][0] + ctrlSetting["warmUpPeriod"],
    #         ctrlSetting["startSimDate"][1],
    #         ctrlSetting["startSimDate"][2]),
    #         freq="D")

    # endSimDate = pandas.Timestamp("{}/{}/{}".format(
    #         ctrlSetting["endSimDate"][0],
    #         ctrlSetting["endSimDate"][1],
    #         ctrlSetting["endSimDate"][2]),
    #         freq="D")
    
    # simMonthRange = pandas.date_range(startSimDate,
    #                                     endSimDate,
    #                                     freq='M')

    # simDateRange = pandas.date_range(startSimDate,
    #                                     endSimDate,
    #                                     freq='D')

    # simYearRange = pandas.date_range(startSimDate,
    #                                     endSimDate,
    #                                     freq='Y')

    obsSimPair = {}


    for obsKey, obsDF in obsDict.items():
        
        # Retrive outlet-IPrint-VarHdrNo information from the key
        obsKeySP = obsKey.split("_")
        outLetNo = int(obsKeySP[0])
        outVarFreq =  int(obsKeySP[1])
        outVarIdinCtrl =  int(obsKeySP[2])
        outVarHdrNm = varIDObsHdrPair[outVarIdinCtrl]

        # Get the correspoinding frequency from the simDF, which 
        # was read from the output.rch.
        # The iPrint in file.cio was determined before based on IPRINT
        # Thus, it will meet the requirement here.
        # However, the extraction of values will need corresponding 
        # aggregation process. The date here is marked as column "MON"
        # in the dataframe.
        # For daily, no aggregation is needed.
        # Frist, extract the corresponding outlet and variable
        # Then, aggregate as needed.
        simOutVarHeader = ["RCH","GIS", "MON", "AREAkm2"] + [outVarHdrNm]
        simOutVarDF = simOutRchDF.loc[simOutRchDF["RCH"] == outLetNo][simOutVarHeader]

        # Add a time series for better matching with observed data
        # if outVarFreq == 1: 
        #     simOutVarDF["Date"] = simDateRange
        # elif outVarFreq == 2:
        #     simOutVarDF["Date"] = simMonthRange
        # elif outVarFreq == 3:
        #     simOutVarDF["Date"] = simYearRange

        # Combine the two pair based on dates
        # obsSimPair[obsKey] = obsDF.merge(simOutVarDF, how="left", left_on="Date", right_on="Date", copy=True)    

        # Updated on Dec 1, 2020
        # The older version tried to match the obs and sim by date.
        # Actually, this might not be necessary as we have checked observed data
        # to be the length of simulateion dates based on the start and end date set
        # in the control.set file. 
        # The simulation was conducted also for the same period.
        # The output was extracted for the same period.
        # So, two lists will be the easist way.   
        obsSimPair[obsKey] = []
        obsSimPair[obsKey].append(list(obsDF[outVarHdrNm]))
        obsSimPair[obsKey].append(list(simOutVarDF[outVarHdrNm]))

    return obsSimPair


##########################################################################
def getRch2DF(fnRch, iPrintForCio, totalRchNum):
    # Read output.rch into pandas dataframe
    # Original header
    # hedr = ["RCH","GIS", "MON","AREAkm2",
    #         "FLOW_OUTcms", " SED_OUTtons", "ORGN_OUTkg",
    #         "ORGP_OUTkg", "NO3_OUTkg", "NH4_OUTkg",
    #         "NO2_OUTkg", "MINP_OUTkg", "SOLPST_OUTmg",
    #         "SORPST_OUTmg"]

    # Corresponding heads used in the observed data
    # facilitating the extraction of output variables.
    hedr = ["RCH","GIS", "MON","AREAkm2",
            "sf(m3/s)", " sed(t/ha)", "orgn(kg/ha)",
            "orgp(kg/ha)", "no3n(kg/ha)", "nh4n(kg/ha)",
            "no2n(kg/ha)", "minp(kg/ha)", "solpst(mg/ha)",
            "sorpst(mg/ha)"]

    # ffRchDataLnRdr = ff.FortranRecordReader('(A5,1X,I5,1X,I8,1X,I3,I3,I5,3X,11E12.4)')
    # fidRch = open(fnRch, "r")
    # lifRch = fidRch.readlines()
    # fidRch.close()
    # del(lifRch[:9])

    # Version 1: normal serial for loop
    # for idx in range(len(lifRch)):
    #     # print(lifRch[idx])
    #     lifRch[idx] = ffRchDataLnRdr.read(lifRch[idx])[1:]

    # Version 2: use pandas read_fwf
    # rchDF = pandas.read_fwf(fnRch, colspecs='infer', skiprows=8)
    

    # Based on the value of iPrintForCio, the format are different
    # When iPrintForCio == 0, for month and iPrintForCio == 2 for annual,
    # The output.rch will include lines for annual sum and year average for
    # each rch.
    # When iPrintForCio == 1 for daily, these does not exist
    colWidth = [5, 6, 9, 6] + [12] * 11
    if (iPrintForCio == 0) or (iPrintForCio == 2):
        rchDF = pandas.read_fwf(fnRch, widths=colWidth, skiprows=9, names=hedr, skipfooter=totalRchNum)
        rchDF = rchDF.loc[rchDF["MON"] < 13]
    else:
        rchDF = pandas.read_fwf(fnRch, widths=colWidth, skiprows=9, names=hedr)
   
    return rchDF



##########################################################################
def determineIprintVal(iPrint):

    """
    This function modify the file.cio file based on the user input.
    These include:
    1. IPRINT need to specified for output.rch file
    2. start and end date
    3. number of skip years.
    4. output variable in the output.rch file to meet the requirement of
    different variables specified in for calibration.
    """
    
    # Determine values
    # IPRINT for file cio
    iPrintForCio = 1
    iPrintForVars = iPrint # ctrlSetting["iPrint"]
    uniqIPrint = list(set(iPrintForVars))
    # As long as 1 is in the unique Iprint values, 1 should 
    # be set for the iPrintForCio
    # If 1 is not, check whether we need to get monthly values
    # Else, the annual output will be printed.
    # Daily in the Control file
    if 1 in uniqIPrint:
        iPrintForCio = 1
    # Monthly in the control file
    elif 2 in uniqIPrint:
        iPrintForCio = 0
    # Annual in the control file
    else:
        iPrintForCio = 2
        
    return iPrintForCio



##########################################################################
def modFileCio(ctrlSetting, runningDir, rchVarLst):

    """
    This function modify the file.cio file based on the user input.
    These include:
    1. IPRINT need to specified for output.rch file
    2. start and end date
    3. number of skip years.
    4. output variable in the output.rch file to meet the requirement of
    different variables specified in for calibration.
    """

    iPrintForCio = determineIprintVal(ctrlSetting["iPrint"])

    # Determin NBYR, IYR, IDFA, and IDAL
    usrStartJD = datetime.date(ctrlSetting["startSimDate"][0], 
                            ctrlSetting["startSimDate"][1], 
                            ctrlSetting["startSimDate"][2]).timetuple().tm_yday

    usrEndJD = datetime.date(ctrlSetting["endSimDate"][0], 
                            ctrlSetting["endSimDate"][1], 
                            ctrlSetting["endSimDate"][2]).timetuple().tm_yday
    firstJDStartYr = datetime.date(ctrlSetting["startSimDate"][0], 
                            1, 1).timetuple().tm_yday
    lastJDEndYr = datetime.date(ctrlSetting["endSimDate"][0], 
                            1, 1).timetuple().tm_yday

    NBYR = ctrlSetting["endSimDate"][0] - ctrlSetting["startSimDate"][0] + 1
    IYR = ctrlSetting["startSimDate"][0]
    IDAF = usrStartJD - firstJDStartYr + 1
    IDAL = usrEndJD - lastJDEndYr + 1

    # Readin file.cio, modify and write new values.
    fnp = os.path.join(runningDir, "file.cio")

    try:
        with open(fnp, 'r') as swatFile:
            lif = swatFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnp, e))
        exit(1)

    for lidx in range(len(lif)):
        # Line 8 for parameter NBYR
        if (lidx == 7):
            lif[lidx] = """{:16d}    | NBYR : Number of years simulated\n""".format(
                    NBYR)
        
        # Line 9 for parameter IYR
        if (lidx == 8):
            lif[lidx] = """{:16d}    | IYR : Beginning year of simulation\n""".format(
                    IYR)
        
        # Line 10 for parameter IDAF
        if (lidx == 9):
            lif[lidx] = """{:16d}    | IDAF : Beginning julian day of simulation\n""".format(
                    IDAF)
        
        # Line 11 for parameter IDAL
        if (lidx == 10):
            lif[lidx] = """{:16d}    | IDAL : Ending julian day of simulation\n""".format(
                    IDAL)
        
        # Line 59 for parameter IPRINT
        if (lidx == 58):
            lif[lidx] = """{:16d}    | IPRINT : print code (month, day, year)\n""".format(
                    iPrintForCio)
        
        # Line 60 for parameter NYSKIP
        if (lidx == 59):
            lif[lidx] = """{:16d}    | NYSKIP : number of years to skip output printing\n""".format(
                    ctrlSetting["warmUpPeriod"])
 
        # Line 65 for parameter out variable in Rch
        if (lidx == 64):
            # Construct out variable lines 4 spacex * 20 var
            lineForWrite = rchVarLst + [0] * (20 - len(rchVarLst))
            lineForWrite = "".join(["{:4d}".format(varRch) for varRch in lineForWrite])
            lif[lidx] = """{}\n""".format(lineForWrite)
 
        # Line 67 for parameter out variable in subbasin
        if (lidx == 66):
            # Construct out variable lines 4 spacex * 20 var
            lif[lidx] = """   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0\n"""
 
        # Line 69 for parameter out variable in HRU
        if (lidx == 68):
            # Construct out variable lines 4 spacex * 20 var
            lif[lidx] = """   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0\n"""
 
        # Line 71 for parameter out variable in HRU
        if (lidx == 70):
            # Construct out variable lines 4 spacex * 20 var
            lif[lidx] = """   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0\n"""
 
        # Line 79 for parameter IA_B
        if (lidx == 78):
            lif[lidx] = """{:16d}    | IA_B : Code for binary output of files (.rch, .sub, .hru files only)\n""".format(
                    1)
 
        # Line 85 for parameter ICALEN, 1 for y/m/d
        if (lidx == 84):
            lif[lidx] = """{:16d}    | ICALEN: Code for printing out calendar or julian dates to .rch, .sub and .hru files\n""".format(
                    0)

    # Then write the contents into the same file
    with open(fnp, 'w') as swatFile:
        swatFile.writelines(lif)

    return iPrintForCio


##########################################################################
def getParmSets():

    """
    This function takes the subarea no from GIS and convert
    them to the SWAT subarea name convention.
    """
    parmSetRead = pandas.read_table(fnParmSetComb, header=0, index_col=0)

    # Remove those non selected parameters
    parmSetSelected = parmSetRead.loc[parmSetRead['selectFlag'] == 1].copy(deep=True)

    # Add two columns to parameter test and best values
    parmSetSelected["BestVal"] = parmSetSelected["InitVal"]
    parmSetSelected["TestVal"] = parmSetSelected["InitVal"]

    # Updated April 19, 2021
    # Add one column to store whether a parameter is selected in one run
    parmSetSelected["ModThisRun"] = [0] * len(parmSetSelected["InitVal"])
    ##########################################################################
    # Generate file extension list and file names ############################
    # Creating a set of parameters that will need to be updated 
    # at the basin level and the subarea and hru level.
    # Then, create an individual copy of parameter at the subarea
    # and hru level for each subarea group.
    # The updating of basin level will be based on the overall
    # objective function value, and those at the subarea/hru leve
    # will be conducted based on the subarea group objective 
    # function value.
    parmBsnLvl = copy.deepcopy(parmSetSelected)
    parmBsnLvl = parmBsnLvl.loc[parmBsnLvl["File"].isin(bsnLvlFlExt)]
    parmBsnLvlFExtLst = parmBsnLvl["File"].unique()

    subruLvlFlExt = subLvlFlExtLst + hruLvlFlExtLst
    parmSubLvl = copy.deepcopy(parmSetSelected)
    parmSubLvl = parmSubLvl.loc[parmSubLvl["File"].isin(subruLvlFlExt)]
    parmSubLvlFExtLst = parmSubLvl["File"].unique()

    return parmBsnLvl, parmBsnLvlFExtLst, parmSubLvl, parmSubLvlFExtLst

##########################################################################
def buildSWATSubFn(subGISNo):

    """
    This function takes the subarea no from GIS and convert
    them to the SWAT subarea name convention.
    """
    
    if (len(subGISNo)) < 2:
        sub_zeros = '0000'
    elif (len(subGISNo) >= 2) and (len(subGISNo) < 3):
        sub_zeros = '000'
    elif (len(subGISNo) >= 3) and (len(subGISNo) < 4):
        sub_zeros = '00'
    elif (len(subGISNo) >= 4) and (len(subGISNo) < 5):
        sub_zeros = '0'
    elif (len(subGISNo) >= 5):
        sub_zeros = ''

    fnSWATSub = "{}{}0000".format(sub_zeros, subGISNo)

    return fnSWATSub
    
    
##########################################################################
def buildSWATHruFn(fnSwatSub, runningDir):

    """
    This function takes the .sub file name, read the file, 
    get the total number of hrus in the subarea, and create
    a list of hru files for each subarea.
    """
    
    fnSWATHruLst = []

    fnpSwatSub = os.path.join(runningDir,
        "{}.sub".format(fnSwatSub))

    try: 
        with open(fnpSwatSub, 'r', encoding="ISO-8859-1") as f:
            lif = f.readlines()
    except IOError as e:
        print("File {} does not exist: {}. The \
            sub GIS no is obtained from the reach.shp\
            please make sure you are using the corresponding \
            shapefile and the files in the txtinout folder. \
            ".format(fnpSwatSub, e))
        exit(1)

    totalHruNo = int(lif[52].split("|")[0])

    subPrefix = os.path.split(fnpSwatSub)[1][:5]
    for hruIdx in range(1, totalHruNo+1):
        # Determint the number of 0 to be added before hru no in the 
        # hru file name.
        if (hruIdx < 10):
            hruZeros = "000"
        elif ((hruIdx >= 10) and (hruIdx < 100)):
            hruZeros = "00"
        elif ((hruIdx >= 100) and (hruIdx < 1000)):
            hruZeros = "0"
        elif (hruIdx > 1000):
            hruZeros = ""
            
        fnSWATHru = "{}{}{}".format(subPrefix, hruZeros, hruIdx) 
        fnSWATHruLst.append(fnSWATHru)

    return fnSWATHruLst


##########################################################################
def updateParValDDS(parRow, perturbFactor):

    """
    Purpose is to generate a neighboring decision variable
    value for a single decision variable value being 
    perturbed by the DDS optimization algorithm.
    New DV (Decision Variable) value respects the 
    upper and lower DV bounds.
    Coded by Bryan Tolson, Nov 2005.

    I/O variable definitions:
    x_cur - current decision variable (DV) value
    x_min - min DV value
    x_max - max DV value
    r  - the neighborhood perturbation factor
    new_value - new DV variable value 
    (within specified min and max)

    Qingyu Feng 20201008
    The code was originally in fortran and is translated into
    python here. The input variables are provided as a 
    row of dataframe.
    """

    # Get a range of parameter value
    parValRange = float(parRow["UpperBound"]) - float(parRow["LowerBound"])

    # Generate a standard normal random variate (zvalue) 
    # Below returns a standard Gaussian random number based
    # upon Numerical recipes gasdev and 
    # Marsagalia-Bray Algorithm

    ranVal = 0.0001
    work1 = 0.0
    work2 = 0.0
    work3 = 2.0

    while ((work3 >= 1.0) or (work3 == 0.0)):
        ranVal1 = numpy.random.uniform(0, 1)
        ranVal2 = numpy.random.uniform(0, 1)
        work1 = 2.0 * float(ranVal1) - 1.0
        work2 = 2.0 * float(ranVal2) - 1.0
        work3 = work1 * work1 + work2 * work2

    work3Base = (-2.0 * math.log(work3))/work3
    work3final = pow(work3Base, 0.5)

	# pick one of two deviates at random
    # (don't worry about trying to use both):
    ranVal3 = numpy.random.uniform(0, 1)
    # print("Random value for updating values:{}, {}, {}".format(ranVal1, ranVal2, ranVal3))
    zvalue = 0.000
    if (ranVal3 < 0.5):
        zvalue = work3final * work1
    else:
        zvalue = work3final * work2

    # done standard normal random variate generation
    
    # Calculate new decision variable value
    parValUpdated = float(parRow["TestVal"]) + zvalue * perturbFactor * parValRange

    # Check new value is within DV bounds.  
    # If not, bounds are reflecting.
    if (parValUpdated < float(parRow["LowerBound"])):
        parValUpdated = float(parRow["LowerBound"]) + (float(parRow["LowerBound"]) - parValUpdated)
        if (parValUpdated > float(parRow["UpperBound"])):
            parValUpdated = float(parRow["LowerBound"])
    elif (parValUpdated > float(parRow["UpperBound"])):
        parValUpdated = float(parRow["UpperBound"]) - (parValUpdated - float(parRow["UpperBound"]))
        if (parValUpdated < float(parRow["LowerBound"])):
            parValUpdated = float(parRow["UpperBound"])

    return parValUpdated


##########################################################################
def initParaObjFile(fn):

    if os.path.isfile(fn):
        os.remove(fn)

    fid = open(fn, 'w')
    fid.close()


##########################################################################
# This three functions were called by pandas.DataFrame.apply function
# to convert the year month and date information into the pandas.
# Timestamp object.
def y2TimeStamp(rowVal):
    newdate = pandas.Timestamp("{}".format(int(rowVal["yyyy"])), freq = "Y")
    return newdate


def ym2TimeStamp(rowVal):
    newdate = pandas.Timestamp("{}-{}".format(
        int(rowVal["yyyy"]), 
        int(rowVal["mm"])),
        freq = "M")
    return newdate


def ymd2TimeStamp(rowVal):
    newdate = pandas.Timestamp("{}-{}-{}".format(
        int(rowVal["yyyy"]), 
        int(rowVal["mm"]),
        int(rowVal["dd"])),
        freq = "D")
    return newdate

##########################################################################
def readOBSSet(ctrlSetting, fdObs, varIDObsHdrPair):
        
    outLetList = ctrlSetting["outLetList"] 
    outputVarList = ctrlSetting["outputVarList"] 
    iPrint = ctrlSetting["iPrint"] 
    usrStartDate = ctrlSetting["startSimDate"] 
    usrEndDate = ctrlSetting["endSimDate"] 
    warmUpPeriod = ctrlSetting["warmUpPeriod"] 
    objFuncLst =  ctrlSetting["objFuncLst"] 
    objFuncWgtLst =  ctrlSetting["objFuncWgtLst"] 

    # Observation are provided in the following format
    # Obs_FREQOLTNO.set
    # Each outlet_variable pair will be considered
    # as one because they will be used to calculate
    # statistics. 
    # Thus, the input will looks like:
    # Total variables: 3
    # Outlet: 6 6 8
    # iPrint: 1 2 3 (daily, monthly, annual)
    # output variables: 1 2 1 (stream flow, sediment, streamflow)
    # First create a dict pair to to get 
    oltVarDFDict = {}
    
    # Construct daily, monthly and annual time series.
    startSimDate = pandas.Timestamp("{}/{}/{}".format(
            usrStartDate[0] + warmUpPeriod,
            usrStartDate[1],
            usrStartDate[2]
            ),
            freq="D")

    endSimDate = pandas.Timestamp("{}/{}/{}".format(
            usrEndDate[0],
            usrEndDate[1], 
            usrEndDate[2]
            ),
            freq="D")

    obsHeader = ["yyyy", "mm", "dd", "sf(m3/s)","sed(t/ha)","orgn(kg/ha)",
                "orgp(kg/ha)","no3n(kg/ha)","nh4n(kg/ha)",
                "no2n(kg/ha)", "minp(kg/ha)", "solpst(mg/ha)",
                "sorpst(mg/ha)", "tp(kg/ha)", "tn(kg/ha)", "tpst(ppb)"]

    for iOlt in range(len(outLetList)):
        # A predefined output header for each outlet variable pair
        dateHeader = ["yyyy", "mm", "dd", "Date"]

        # Daily 
        if iPrint[iOlt] == 1:
            fnObs = "obs_daily{}.prn".format(outLetList[iOlt])
            fdnObs = os.path.join(fdObs, fnObs)

            # Check file existence:
            if not os.path.isfile(fdnObs):
                print("File does {} not exist!".format(fdnObs))
                sys.exit()
            else:
                obsDF = pandas.read_table(fdnObs, names=obsHeader, skiprows=1)
                obsDF["Date"] = obsDF.apply(ymd2TimeStamp, axis=1)

        # Monthly
        elif iPrint[iOlt] == 2:
            fnObs = "obs_monthly{}.prn".format(outLetList[iOlt])
            fdnObs = os.path.join(fdObs, fnObs)

            # Check file existence:
            if not os.path.isfile(fdnObs):
                print("File does {} not exist!".format(fdnObs))
                sys.exit()
            else:
                obsDF = pandas.read_table(fdnObs, names=obsHeader, skiprows=1)
                obsDF["Date"] = obsDF.apply(ym2TimeStamp, axis=1)

        # Annual
        elif iPrint[iOlt] == 3:
            fnObs = "obs_yearly{}.prn".format(outLetList[iOlt])
            fdnObs = os.path.join(fdObs, fnObs)

            # Check file existence:
            if not os.path.isfile(fdnObs):
                print("File does {} not exist!".format(fdnObs))
                sys.exit()
            else:
                obsDF = pandas.read_table(fdnObs, names=obsHeader, skiprows=1)
                obsDF["Date"] = obsDF.apply(y2TimeStamp, axis=1)

        # Check the start and end date by User in the control file and
        # in the observed Data.
        obsStartDate = obsDF.iloc[[0], [0, 1, 2]]
        obsEndDate = obsDF.iloc[[-1], [0, 1, 2]]

        # Construct daily, monthly and annual time series.
        obsStartDateDT = pandas.Timestamp("{}/{}/{}".format(
                int(obsStartDate["dd"]),
                int(obsStartDate["mm"]), 
                int(obsStartDate["yyyy"])
                ),
                freq="D")

        obsEndDateDT = pandas.Timestamp("{}/{}/{}".format(
                int(obsEndDate["dd"]),
                int(obsEndDate["mm"]), 
                int(obsEndDate["yyyy"])
                ),
                freq="D")

        if startSimDate > obsEndDateDT:
            print("The started date in the control file id later than the observed file in outlet no {} with IPRINT value {}".format(
                outLetList[iOlt], iPrint[iOlt]))
            exit(1)
        elif obsStartDateDT > endSimDate:
            print("The end date of the observed data for outlet no {} with IPRINT value {} is earlier than the start date specified in the control file".format(
                outLetList[iOlt], iPrint[iOlt]))
            exit(1)

        # Trim obsDF to the user specified start and end date
        obsDF = obsDF.loc[obsDF["Date"] >= startSimDate]
        obsDF = obsDF.loc[obsDF["Date"] <= endSimDate]

        # Get the interested columns from the whole data frame
        outVarIdinCtrl = outputVarList[iOlt]
        varColNms = dateHeader + [varIDObsHdrPair[outVarIdinCtrl]]

        # Key: outletno, iprint, outVarID, statID, statWeight
        oVKey = "{}_{}_{}_{}_{:.1f}".format(outLetList[iOlt],
                                        iPrint[iOlt], 
                                        outVarIdinCtrl,
                                        objFuncLst[iOlt],
                                        objFuncWgtLst[iOlt],
                                        )
        oltVarDFDict[oVKey] = obsDF.loc[:, varColNms]

    return oltVarDFDict


##########################################################################
def ctrlSetToJSON(fnCtrlSet, fnCtrlJSON):
    
    """
    This function read in the para_set file and convert it into a json file,
    which is used later for grouping.
    """

    # First step: define parameter set
    fid = open(fnCtrlSet, "r")
    lif = fid.readlines()
    fid.close()

    outdict = {}

    # Control file has a strict format now.
    # I will get the corresponding lines for results
    # Line 2 for perturbation factor
    outdict["perturbFactor"] = float(lif[1].split("|")[0])
    
    # Line 3 for total number of model evaluation
    outdict["totalModelRuns"] = int(lif[2].split("|")[0])
    
    # Line 5 for Initial parameter index
    # 0: Use Random Para
    # 1: Use Initial Para
    outdict["initParmIdx"] = int(lif[4].split("|")[0])

    # Line 6 for restart mechanism
    # 0: Do Not Restart
    # 1: Restart
    outdict["restartIdx"] = int(lif[5].split("|")[0])

    # Line 9 for total variable Nos
    # Based on the later code that uses this variable, 
    # it is the total outLet-Variable No. 
    # Meaning that if you are calibrating flow for two
    # outlets, or flow and sediment for the same outlet,
    # this number is two.
    outdict["totalVarNo"] = int(lif[8])

    # Line 11 for soft data contraints
    outdict["totalSoftNo"] = int(lif[10])

    # Line 13 for list of outlets separated by space
    # If you have two variable for one outlet, then,
    # the number of the oultet should be listed twice.
    # Example: 6 6
    # Or if you have the same variables for two outlet,
    # the number of the two outlets will be listed.
    lif[12] = lif[12].split(" ")
    lif[12][-1] = lif[12][-1][:-1]
    outdict["outLetList"] = list(map(int, lif[12]))

    # Line 15 for print space
    lif[14] = lif[14].split(" ")
    lif[14][-1] = lif[14][-1][:-1]
    outdict["iPrint"] = list(map(int, lif[14]))

    # Line 17 for list of variables
    # Output Variables (0-15): 
    # '0- Baseflow', 
    # '1-Stream Flow',
    # '2-Sediment',
    # '3-org-N',
    # '4-org-P',
    # '5-NO_3-N',
    # '6-NH_4-N',
    # '7-NO_2-N',
    # '8-Mineral-P',
    # '9-Soluble Pesticide',
    # '10-Sorbed Pesticide',
    # '11-Total Phosphorus',
    # '12-Total Nitrogen',
    # '13-Total Pesticide',
    # '14-TKN',
    # '15-NO2+NO3'
    lif[16] = lif[16].split(" ")
    lif[16][-1] = lif[16][-1][:-1]
    outdict["outputVarList"] = list(map(int, lif[16]))

    # Line 19 for statistics
    # Statistics (as Objective Functions):
    # '1=1-NS', 
    # '2=BIAS',
    # '3=RMSE',
    # '4=R2',
    # '5=MSE'
    # These statistics are list of statistics for
    # each outlet-variable pair.
    lif[18] = lif[18].split(" ")
    lif[18][-1] = lif[18][-1][:-1]
    outdict["objFuncLst"] = list(map(int, lif[18]))

    # Line 21 for weights of objective functions
    # Originally, this is used combine the variables 
    # into one to get the minimum of objective functions
    # as the stop criteria.
    # TODO: In this distributed calibrated, this needs further
    # thoughts on how to use this wisely.
    lif[20] = lif[20].split(" ")
    lif[20][-1] = lif[20][-1][:-1]
    outdict["objFuncWgtLst"] = list(map(float, lif[20]))

    # Line 23 for warm up period in file.fig
    outdict["warmUpPeriod"] = int(lif[22])

    # Line 25, 26 for total simulation period
    lif[24] = lif[24].split(" ")
    lif[24][-1] = lif[24][-1][:-1]
    outdict["startSimDate"] = list(map(int, lif[24]))
    lif[25] = lif[25].split(" ")
    lif[25][-1] = lif[25][-1][:-1]
    outdict["endSimDate"] = list(map(int, lif[25]))

    # Line 28 for constraints applied
    # Apply Constraints (1 = YES, 0 = NO) for:
    # (1)DENITRIFICATION, 
    # (2) NO3 LEACHED, 
    # (3) P LEACHED, and 
    # (4) SSQ/(SQ+SSQ) NO3 Yield
    # example input: 1 1 1 1 or 0 0 0 0
    lif[27] = lif[27].split(" ")
    lif[27][-1] = lif[27][-1][:-1]
    outdict["constraintLst"] = list(map(int, lif[27]))

    # Line 30, 31 for constraints lower and upper bounds
    # Apply Constraints (1 = YES, 0 = NO) for:
    # (1)DENITRIFICATION, 
    # (2) NO3 LEACHED, 
    # (3) P LEACHED, and 
    # (4) SSQ/(SQ+SSQ) NO3 Yield
    # example input: 1 1 1 1 or 0 0 0 0
    lif[29] = lif[29].split(" ")
    lif[29][-1] = lif[29][-1][:-1]
    outdict["constraintLBLst"] = list(map(float, lif[29]))
    lif[30] = lif[30].split(" ")
    lif[30][-1] = lif[30][-1][:-1]
    outdict["constraintUBLst"] = list(map(float, lif[30]))

    # Line 33 for input uncertainty index
    # Apply Constraints (1 = YES, 0 = NO) for:
    # 0 = exclude Input uncertainty ,
    # 1 = include Input uncertainty)
    outdict["inputUncertIdx"] = int(lif[32])

    # Line 35 for measurement uncertainty index
    #  0 = w/o measurement uncertainty ,
    #  1 = Probable Error Range , 
    #  2 = Probability Distribution)
    outdict["measureUncertIdx"] = int(lif[34])

    # Line 37 for measurement error (%)
    outdict["measureError"] = float(lif[36])

    # Line 39 for whether grouping subarea based on outlet
    outdict["groupSubareaIdx"] = int(lif[38])

    # Line 41 for run no of best parameter sets for outlets
    # lif[40] = lif[40].split(" ")
    # lif[40][-1] = lif[40][-1][:-1]
    # outdict["runNoBestPar"] = list(map(int, lif[40]))
    outdict["runNoBestPar"] = int(lif[40])

    # Updated Oct 28, 2021 
    # Added 4 variables for sensitivity analysis
    # Line 43 for the sensitivity analysis sampling argument 
    # with the stallite sample method
    outdict["saMethod"] = int(lif[42])

    # Line 45 for the sensitivity analysis sampling argument 
    # with the stallite sample method
    outdict["saltelliArg"] = int(lif[44])

    # Line 45 for the sensitivity analysis sampling argument 
    # with the morris sample method
    outdict["morrisRes"] = int(lif[46])

    # Line 45 for the sensitivity analysis sampling argument 
    # with the fast sample method
    outdict["fastRes"] = int(lif[48])
    # Updated Oct 28, 2021 
    # TODO: Add input check program, either here or later
    # as a separate function
    # Check the settings for potential errors
    if ((len(outdict["outLetList"]) != outdict["totalVarNo"])
     or (len(outdict["iPrint"]) != outdict["totalVarNo"])
     or (len(outdict["outputVarList"]) != outdict["totalVarNo"])
     or (len(outdict["outputVarList"]) != outdict["totalVarNo"])
    ):
        print("Please check the number of total outlet-variable No, \
            outlet list, iprint, output Varables. They should match each \
            other!")

        return "Error"
    else:
        with open(fnCtrlJSON, 'w') as outfile:
            json.dump(outdict, outfile)

        return outdict




##########################################################################
def parmSetToJSON(fnParmSet, fnParmJSON):
    
    """
    This function read in the para_set file and convert it into a json file,
    which is used later for grouping.
    """

    # First step: define parameter set
    fid = open(fnParmSet, "r")
    lif = fid.readlines()
    fid.close()

    del(lif[0])

    outdict = {}

    for idx in range(len(lif)):
        lif[idx] = lif[idx].split("\t")
        lif[idx][-1] = lif[idx][-1][:-1]

        outdict[lif[idx][1]] = {}
        outdict[lif[idx][1]]["orderNo"] = lif[idx][0]
        outdict[lif[idx][1]]["inFile"] = lif[idx][2]
        outdict[lif[idx][1]]["unit"] = lif[idx][3]
        outdict[lif[idx][1]]["initVal"] = lif[idx][4]
        outdict[lif[idx][1]]["select"] = lif[idx][5]
        outdict[lif[idx][1]]["lowerBd"] = lif[idx][6]
        outdict[lif[idx][1]]["upperBd"] = lif[idx][7]

    with open(fnParmJSON, 'w') as outfile:
        json.dump(outdict, outfile)

    return outdict





##########################################################################
def get_osplatform():

    platforms = {
        'linux1' : 'Linux',
        'linux2' : 'Linux',
        'darwin' : 'OS X',
        'win32' : 'Windows'
    }
    if sys.platform not in platforms:
        return sys.platform

    return platforms[sys.platform]



##########################################################################
def runSWATModel(osplatform, runningDir, fdMain, fdDMPOTpyFiles, fnSwatExe):
    """
    This function change the directory to the apex run folder,
    run apex, and change back.
    TODO: Add the swat2012 linux into the package folder
    and copy it to the working directory.
    """

    # Copy the swat exe to the working Dir
    fnpSWATExe = os.path.join(fdMain, 
                    fdDMPOTpyFiles, 
                    fnSwatExe)
    fnpSWATExeWD = os.path.join(runningDir, 
                    fnSwatExe)
    try:
        shutil.copy(fnpSWATExe, fnpSWATExeWD)
    except: 
        print("could not copy SWAT executable to the working directory")
        return

    os.chdir(runningDir)

    if osplatform == "linux":
        procCommand = "./{}".format(fnSwatExe)                                  
    elif osplatform == "Windows":
        procCommand = "{}".format(fnSwatExe)

    ## call date command ##
    try: 
        p = subprocess.run(procCommand, shell=True, stderr=subprocess.PIPE)
        # # Change back to the fdworking folder
        os.chdir(fdMain)

        return p.returncode()
    except:
        return "Error"




##########################################################################
def determineFileLoopFlag(
    parmFileExtLst,
    subLvlFlExtLst, 
    hruLvlFlExtLst):

    """
    This function determines wheter the sub and hru level should be
    started.
    Input: list of file extentions in the parameter dataframe
    Output: two flags, one for subarealevel and one for hru level.
    """

    flagSubLoop = 0
    flagHruLoop = 0

    for fExt in parmFileExtLst:
        if fExt in subLvlFlExtLst:
            flagSubLoop = 1
            break

    for fExt2 in parmFileExtLst:
        if fExt2 in hruLvlFlExtLst:
            flagHruLoop = 1
            break        

    return flagSubLoop, flagHruLoop


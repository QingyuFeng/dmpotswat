#!/usr/bin/python3

"""
This program was developed to conduct sensitivity analysis for SWAT model.
For sensitivity analysis, the division of watershed into subarea groups
are not required since we are using the same number of input parameter lists.

The sensitivity analysis is conducted for same variables in one run.
For example, for flow, only varaibles for flow can be included.
And the variables in the control setting need to be the same for all outlets.
The program will check and see whether this requirement is meet.

If the users select variables for an outlet that are not for the target
variables, the output will not make sense.

"""

##########################################################################
# Import modules #########################################################
##########################################################################

import sys, os
import datetime
import pytz

from pyscripts.globVars import *
from pyscripts.DMPOTUtil import *
from pyscripts.SAUtil import *

from SALib.analyze import sobol, fast
from SALib.sample import saltelli, fast_sampler
import SALib.sample.morris as morris_sample
import SALib.analyze.morris as morris_analyze

# import SALib.analyze as SALib_Analyze
# import SALib.sample as SALib_Sample

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

nCalVal = "Sensitivity Analysis"
if not os.path.isdir(fdSA):
    os.mkdir(fdSA)

##########################################################################
# Read Settings ##########################################################
print(".....Read in DMPOT setting from the control file.....")
ctrlSetting = ctrlSetToJSON(fnCtrlSetUsr, fnCtrlJsonUsr)

if (ctrlSetting == "Error"):
    print("Your control Setting was not corrected prepared. Please double check!")
    sys.exit(1)

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
# All selected parameters were included.
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
ctrlSetting["groupSubareaIdx"] = 0

subNoGroups, parmObjFnKeys, rcvRchLst = getSubGroupsRchList(ctrlSetting)

subParaFns, subObjFunFns, subParaSel01Fns = initOutFileParmObjSublvl(
    ctrlSetting, parmSubLvl, parmObjFnKeys, rcvRchLst)
    
fnBsnPara, fnBsnParaSel = initOutFileParmObjBsnlvl(parmBsnLvl)

subObfBestDict, subObfTestDict, bsnObfBest = initOFValDict(parmObjFnKeys)

subParGroups, swatSubFnGroups, swatHruFnGroups = initParmInFilenameSubLvl(
    subNoGroups, parmSubLvl, fdWorkingDir)

##########################################################################
# Start Sensitivity analysis  for all runs ##############################
##########################################################################
print(".....DMPOT Sensitivity analysis procedure start.....")

##########################################################################
##########################################################################

runEndTime = datetime.datetime.now(timeZone)
runTimeTotal = runEndTime - startTime
# Step 1: Define model parameters for SALib input
parmForSA = geneParmSetsForSA(parmBsnLvl, parmSubLvl)

# Step 2: # Generate samples
# The Saltelli sampler generates N*(2D+2) samples, where in this example
# N is the argument we supplied and D is the number of model inputs.
if ctrlSetting["saMethod"] == 1: # For Sobol method
    # Convergence properties of the Sobol' sequence is only valid if
    # 'N' (20) is equal to '2^n'
    sampleMethod = "saltelli"
    saltelliArgument = ctrlSetting["saltelliArg"]
    if not checkSaltelliArgument(saltelliArgument):
        print("The argument for saltelli sample need to be a value that is power of 2")
        exit()
    parmSamples = saltelli.sample(parmForSA, saltelliArgument)
elif ctrlSetting["saMethod"] == 2: # For Morris method
    sampleMethod = "Morris"
    parmSamples = morris_sample.sample(parmForSA, ctrlSetting["morrisRes"], num_levels=4)
elif ctrlSetting["saMethod"] == 3: # For FAST method
    sampleMethod = "Fast"
    # SALib.sample.fast_sampler.sample(problem, N, M=4, seed=None)[source]
    parmSamples = fast_sampler.sample(parmForSA, ctrlSetting["fastRes"], M=4)

fnpParmSamples = os.path.join(fdSA, "parmSample_{}.txt".format(sampleMethod))
np.savetxt(fnpParmSamples, parmSamples, fmt='%.3f',delimiter=' ')

# Step 3: Evaluate the model with sampled results and
# store the output into a file
# At this version, the average value will be used.
# Another question to consider is the number of sites used.

# In order to save the time to rewrite the duplicated codes for modifying 
# file values, the parameter values will be updated based on the values
# in each row in the parameter dataframe.

# A list to store the value for allruns, which will be written into a file
# There are multiple oulets, each outlet will have one list representing
# the output for all runs.
outLetAvgAnnList = {}
for oltNoIdx in range(len(ctrlSetting["outLetList"])):
    outLetNo = ctrlSetting["outLetList"][oltNoIdx]
    outLetAvgAnnList[outLetNo] = np.zeros(len(parmSamples))

# for saRunIdx in range(1):#len(parmSamples)):
for saRunIdx in range(len(parmSamples)):
    print("Executing sensitivity analysis run no: {}".format(saRunIdx))

    modifyStartTime = datetime.datetime.now(timeZone)

    # The parameter values are in accordance with the symbols. 
    # Create a dictionary for the values
    parmValDict = generateParmValDict(parmForSA["names"], parmSamples[saRunIdx])

    if iModParm:
        # We are not grouping the subarea. However, the file names of subareas and
        # hrus will need to be generated for modifying. 
        if len(parmSubLvl.index) > 0:
            print(".....Modifying subarea level parameter values with sampled values.....")
            for subGPKey, subGPL in subParGroups.items():
                # Update parameter completely random. 
                subGPL = updateParmInDf(subGPL, parmValDict)
                # Update parameter values in file, run model, and calculate
                # objective function
                modifyParInFileSub(subGPL, 
                                parmSubLvlFExtLst, 
                                swatSubFnGroups, 
                                subGPKey, 
                                swatHruFnGroups,
                                fdWorkingDir)
        if len(parmBsnLvl.index) > 0:
            print(".....Modifying basin level parameter values with sampled values.....")
            parmBsnLvl = updateParmInDf(parmBsnLvl, parmValDict)

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

    # Extract average annual values for each outlet for evaluation
    outLetAvgAnnList = extractAvgAnnEachGroup( 
            fdWorkingDir, 
            iPrintForCio, 
            ctrlSetting["outLetList"], 
            ctrlSetting["outputVarList"],
            outLetAvgAnnList,
            rcvRchLst,
            varIDObsHdrPair, 
            saRunIdx)

# After the run, save the values into a file for later evaluation.
for oltNoIdx in range(len(ctrlSetting["outLetList"])):
    outLetNo = ctrlSetting["outLetList"][oltNoIdx]
    fnOutletVarMeanVal = os.path.join(fdSA, "oltVal_{}_{}.txt".format(outLetNo, sampleMethod))
    np.savetxt(fnOutletVarMeanVal, outLetAvgAnnList[outLetNo])
         
    if ctrlSetting["saMethod"] == 1: # For Sobol method
        # SALib.analyze.sobol.analyze(problem, Y, calc_second_order=True,
        # num_resamples=100, conf_level=0.95, print_to_console=False,
        #  parallel=False, n_processors=None, keep_resamples=False, 
        # seed=None)
        saOutput = sobol.analyze(parmForSA, 
                        outLetAvgAnnList[outLetNo],
                        calc_second_order=True,
                        conf_level=0.95,
                        print_to_console = False)
        # Write the output as a dataframe and into a file
        saOutputDF = saOutput.to_df()
        fnOutletSAOutput_Total = os.path.join(fdSA, "oltSATotal_{}_{}.txt".format(outLetNo, sampleMethod))
        saOutputDF[0].to_csv(fnOutletSAOutput_Total)
        fnOutletSAOutput_First = os.path.join(fdSA, "oltSAFirst_{}_{}.txt".format(outLetNo, sampleMethod))
        saOutputDF[1].to_csv(fnOutletSAOutput_First)

        
    elif ctrlSetting["saMethod"] == 2: # For Morris method
        # SALib.analyze.morris.analyze(problem: Dict, X: numpy.ndarray,
        #  Y: numpy.ndarray, num_resamples: int = 100,
        #  conf_level: float = 0.95, print_to_console: bool = False,
        #  num_levels: int = 4, seed=None) → numpy.ndarray
        # Returns a dictionary with keys ‘mu’, ‘mu_star’, ‘sigma’,
        # and ‘mu_star_conf’, where each entry is a list of parameters
        # containing the indices in the same order as the parameter file.
        saOutput = morris_analyze.analyze(parmForSA, 
                        parmSamples, 
                        outLetAvgAnnList[outLetNo],
                        conf_level=0.95,
                        print_to_console=False, 
                        num_levels=4
                        )

        # Write the output as a dataframe and into a file
        saOutputDF = saOutput.to_df()
        fnOutletSAOutput_Mu = os.path.join(fdSA, "oltSA_{}_{}.txt".format(outLetNo, sampleMethod))
        saOutputDF.to_csv(fnOutletSAOutput_Mu)

    elif ctrlSetting["saMethod"] == 3: # For FAST method
        # SALib.analyze.fast.analyze(problem, 
        # Y, M=4, num_resamples=100,
        #  conf_level=0.95, print_to_console=False, seed=None)
        # Returns a dictionary with keys ‘S1’ and ‘ST’, 
        # where each entry is a list of size D (the number of parameters)
        # containing the indices in the same order as the parameter file.
        saOutput = fast.analyze(parmForSA, 
                    outLetAvgAnnList[outLetNo],
                    M=4, 
                    num_resamples=100,
                    conf_level=0.95,
                    print_to_console = False)
        # Write the output as a dataframe and into a file
        saOutputDF = saOutput.to_df()
        fnOutletSAOutput_Fast = os.path.join(fdSA, "oltSATotal_{}_{}.txt".format(outLetNo, sampleMethod))
        saOutputDF.to_csv(fnOutletSAOutput_Fast)

calStatTime = datetime.datetime.now(timeZone)
print("Time for calculating statistics: {}; Total Time: {}". format(
        calStatTime - runEndTime,
        calStatTime - startTime)) 
print("=================================================================")


# End of for loop for total runs
print("--------------------------------------")
print("Congratulations!!! It's done nicely~~~")
print("--------------------------------------")

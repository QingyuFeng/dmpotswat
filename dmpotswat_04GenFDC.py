#!/usr/bin/python3

"""
This program was developed to generate flow duration curves
for observed and simulated flow

"""

##########################################################################
# Import modules #########################################################
##########################################################################

import sys, os
import datetime
import pytz, math

from pyscripts.globVars import *
from pyscripts.DMPOTUtil import *


# Initial settings of FDC generation
iObsOnlyFDC = False
iSimObsFDC = True

nCalVal = "Cal"
nCalVal2 = "Calibration"
fdCalibrated = fdCalibrated + nCalVal2
if not os.path.isdir(fdCalibrated):
    os.mkdir(fdCalibrated)
fdRchSrc = fdCalibrated

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
# Get frequency
iPrintForCio = determineIprintVal(ctrlSetting["iPrint"])


subNoGroups, parmObjFnKeys, rcvRchLst = getSubGroupsRchList(ctrlSetting)

##########################################################################
# Start optimization procedure for all runs ##############################
##########################################################################
print(".....Generating FDC procedure start.....")


if iSimObsFDC:
    # Check whether output.rch exists in cal folder ##########################
    fnRch = os.path.join(fdRchSrc, "output.rch")

    if not os.path.isfile(fnRch):
        print("""You selected to generate FDC for simulated flow.""")
        print("""But the output.rch file was not found in the specified folder.""")
        print("""Please check your settings to use this feature!""")

    # There are still some preparation to be done
    # Initialize the objective function Value
    # This is a dictionary containing the 5 statistics for each outlet_variable
    # list. 
    try:
        rchDFWhole = getRch2DFPlot(fnRch, iPrintForCio, len(rcvRchLst))
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnRch, e))
        exit(1)

    # Then construct series of observed and simulated pairs
    # for stat calculation
    obsSimPair = buildObsSimPairForPlot(obsDataLst,
                                rchDFWhole,
                                varIDObsHdrPair)

    if ctrlSetting["groupSubareaIdx"] == 0:
        nlumpdist = "Lumped"
    elif ctrlSetting["groupSubareaIdx"] == 1:
        nlumpdist = "Distributed"

    fnOutFig = os.path.join(fdOutputs, "FDCObsVsSimAllOutlets_{}_{}.png".format(nCalVal, nlumpdist))
    # The number of subplots need to be specified by users
    # for customization
    noColFig = 3
    noRowFig = math.ceil(len(list(obsSimPair.keys()))/noColFig)
    figWidth = 21/2.54-2
    figHeight = 29.7/2.54/6 * noRowFig

    fig, axes = plt.subplots(noRowFig, noColFig, 
                        figsize=(figWidth, figHeight), 
                        dpi=300,
                        tight_layout=True)
    axesCtr = 0
    axes = axes.flatten()

    for key, val in obsSimPair.items():
        outLetNo = key.split("_")[0]
        # Obs = val[0], sim = val[1]
        print("Generating FDC of Observed vs Simulated for outlet {}".format(outLetNo))
        flow_duration_curve(val[1], 
                            nCalVal,
                            axes,
                            axesCtr,
                            outLetNo,
                            nlumpdist,
                            comparison=val[0],
                            )
        axesCtr += 1
        
    fig.savefig(fnOutFig, bbox_inches="tight")


    # Generate single FDC plots
    

    for key2, val2 in obsSimPair.items():

        outLetNo = key2.split("_")[0]
        fnOutFigSingle = os.path.join(fdOutputs, "FDCObsVsSimSingleOutlet_{}_{}_{}.png".format(outLetNo, nCalVal, nlumpdist))
        # The number of subplots need to be specified by users
        # for customization
        noColFig = 1
        noRowFig = 1
        figWidth = 4
        figHeight = 4

        fig, axes = plt.subplots(noRowFig, noColFig, 
                            figsize=(figWidth, figHeight), 
                            dpi=300,
                            tight_layout=True)
        
        # Obs = val[0], sim = val[1]
        print("Generating Single FDC of Observed vs Simulated for outlet {}".format(outLetNo))
        flow_duration_curve_single(val2[1],  
                            nCalVal,
                            axes,
                            outLetNo,
                            nlumpdist,
                            comparison=val2[0])
        
        fig.savefig(fnOutFigSingle, bbox_inches="tight")

        
if iObsOnlyFDC:
    fig, ax = plt.subplots(1,1)
 
    for key, value in obsDataLst.items():
        flow_duration_curve(value["sf(m3/s)"].to_list())
        


# End of for loop for total runs
print("--------------------------------------")
print("Congratulations!!! It's done nicely~~~")
print("--------------------------------------")

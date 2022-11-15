# -*- coding: utf-8 -*-
"""
Created on Tue May 12 2020 

This class is designed to be a collection of functions dealing with
the SWAT files.

@author: Qingyu.Feng
"""

##########################################################################
# Import modules #########################################################
##########################################################################

import os
import math, copy
from shutil import copyfile

##########################################################################
# Define functions #######################################################
##########################################################################
# TODO: Double check the read function for swq file, which needs special 
# encoding to process the charactor of degree C.
# TODO: ESCO has two places, the basins.bsn file and hru file. Do I need
# to modify both? or Just the HRU level is good?
# TODO: Check USLE_C to make sure that it is less than 1.0


##########################################################################
def updateParInRes(parInFile, fdWorkingDir, subNoWithRes, swatSubFnGroups, fdmodelTxtInOut):

    for subresNo in subNoWithRes:
        fnResNo = ""

        for subNo in swatSubFnGroups:
            if subresNo in subNo:
                fnResNo = subNo
                break
        
        fnpOrig = os.path.join(fdmodelTxtInOut,
            "{}.res".format(fnResNo))
        try:
            with open(fnpOrig, 'r') as swatFileOrig:
                lifOrig = swatFileOrig.readlines()
        except IOError as e:
            print("File {} does not exist: {}. Please double check your TxtInOut \
                folder and make sure you have a complete set".format(fnpOrig, e))
            exit(1)

        fnpRes = os.path.join(fdWorkingDir,
            "{}.res".format(fnResNo))
        try:
            with open(fnpRes, 'r') as resFile:
                lifRes = resFile.readlines()
        except IOError as e:
            print("File {} does not exist: {}. Please double check your TxtInOut \
                folder and make sure you have a complete set".format(fnpRes, e))
            exit(1)

        # Only values in crop and pasture land are updated.
        # Get a list of variables selected. Since we can get in this function,
        # at least one variable in the sub file was selected. 
        varLst = parInFile["Symbol"].unique()

        for lidx in range(len(lifRes)):

            # Line 5 for parameter RES_ESA 
            # If the original value is 0, the program will still change it, but the new value will
            # still be 0.
            if ((lidx == 4) and ("RES_ESA" in varLst)):
                if (int(parInFile.loc[parInFile["Symbol"] == "RES_ESA"]["selectFlag"]) == 1):
                    origVal = float(lifOrig[lidx].split("|")[0])
                    mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "RES_ESA"]["TestVal"])
                    newVal = origVal * mutiPlier
                    lifRes[lidx] = """{:16.3f}    | RES_ESA : Reservoir surface area when the reservoir is filled to the  emergency spillway [ha]\n""".format(
                        origVal * mutiPlier)
            
            # Line 6 for parameter RES_EVOL 
            # If the original value is 0, the program will still change it, but the new value will
            # still be 0.
            if ((lidx == 5) and ("RES_EVOL" in varLst)):
                if (int(parInFile.loc[parInFile["Symbol"] == "RES_EVOL"]["selectFlag"]) == 1):
                    origVal = float(lifOrig[lidx].split("|")[0])
                    mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "RES_EVOL"]["TestVal"])
                    newVal = origVal * mutiPlier
                    lifRes[lidx] = """{:16.3f}    | RES_EVOL : Volume of water needed to fill the reservoir to the emergency spillway (104 m3)\n""".format(
                        origVal * mutiPlier)

            # Line 7 for parameter RES_PSA 
            if ((lidx == 6) and ("RES_PSA" in varLst)):
                if (int(parInFile.loc[parInFile["Symbol"] == "RES_PSA"]["selectFlag"]) == 1):
                    origVal = float(lifOrig[lidx].split("|")[0])
                    mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "RES_PSA"]["TestVal"])
                    newVal = origVal * mutiPlier
                    lifRes[lidx] = """{:16.3f}    | RES_PSA : Reservoir surface area when the reservoir is filled to the principal spillway [ha]\n""".format(
                        origVal * mutiPlier)

            # Line 8 for parameter RES_PVOL 
            if ((lidx == 7) and ("RES_PVOL" in varLst)):
                if (int(parInFile.loc[parInFile["Symbol"] == "RES_PVOL"]["selectFlag"]) == 1):
                    origVal = float(lifOrig[lidx].split("|")[0])
                    mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "RES_PVOL"]["TestVal"])
                    newVal = origVal * mutiPlier
                    lifRes[lidx] = """{:16.3f}    | RES_PVOL : Volume of water needed to fill the reservoir to the principal spillway [104 m3]\n""".format(
                        origVal * mutiPlier)
     
            # Line 9 for parameter RES_VOL 
            if ((lidx == 8) and ("RES_VOL" in varLst)):
                if (int(parInFile.loc[parInFile["Symbol"] == "RES_VOL"]["selectFlag"]) == 1):
                    origVal = float(lifOrig[lidx].split("|")[0])
                    mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "RES_VOL"]["TestVal"])
                    newVal = origVal * mutiPlier
                    lifRes[lidx] = """{:16.3f}    | RES_VOL : Initial reservoir volume [104 m3]\n""".format(
                        origVal * mutiPlier)

            # Line 10 for parameter RES_SED 
            if ((lidx == 9) and ("RES_SED" in varLst)):
                if (int(parInFile.loc[parInFile["Symbol"] == "RES_SED"]["selectFlag"]) == 1):
                    origVal = float(lifOrig[lidx].split("|")[0])
                    mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "RES_SED"]["TestVal"])
                    newVal = origVal * mutiPlier
                    lifRes[lidx] = """{:16.3f}    | RES_SED : Initial sediment concentration in the reservoir [mg/l]\n""".format(
                        origVal * mutiPlier)
                        
            # Line 11 for parameter RES_NSED 
            if ((lidx == 10) and ("RES_NSED" in varLst)):
                if (int(parInFile.loc[parInFile["Symbol"] == "RES_NSED"]["selectFlag"]) == 1):
                    origVal = float(lifOrig[lidx].split("|")[0])
                    mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "RES_NSED"]["TestVal"])
                    newVal = origVal * mutiPlier
                    lifRes[lidx] = """{:16.3f}    | RES_NSED : Normal sediment concentration in the reservoir [mg/l]\n""".format(
                        origVal * mutiPlier)

            # Line 23 for parameter RES_RR 
            if ((lidx == 22) and ("RES_RR" in varLst)):
                if (int(parInFile.loc[parInFile["Symbol"] == "RES_RR"]["selectFlag"]) == 1):
                    origVal = float(lifOrig[lidx].split("|")[0])
                    mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "RES_RR"]["TestVal"])
                    newVal = origVal * mutiPlier
                    lifRes[lidx] = """{:16.3f}    | RES_RR : Average daily principal spillway release rate [m3/s]\n""".format(
                        origVal * mutiPlier)
     
            # Line 27 for parameter NDTARGR 
            if ((lidx == 26) and ("NDTARGR" in varLst)):
                if (int(parInFile.loc[parInFile["Symbol"] == "NDTARGR"]["selectFlag"]) == 1):
                    origVal = float(lifOrig[lidx].split("|")[0])
                    mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "NDTARGR"]["TestVal"])
                    newVal = origVal * mutiPlier
                    lifRes[lidx] = """{:16.2f}    | NDTARGR : Number of days to reach target storage from current reservoir storage\n""".format(
                        origVal * mutiPlier)

        # Line 29 and 31 for parameter STARG
        # The format is different.
        iresco = int(lifOrig[13].split("|")[0])
        if iresco == 2:
            fldSeasonBegin = int(lifOrig[24].split("|")[0])
            fldSeasonEnd = int(lifOrig[25].split("|")[0])
            if ((fldSeasonBegin > 0) and (fldSeasonEnd > 0)):
                monIdx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
                fldSenIdx = []
                for midx in monIdx:
                    if (midx >= (fldSeasonBegin - 1) and midx <= (fldSeasonEnd - 1)):
                        fldSenIdx.append(midx)
                nonfldSenIdx = []
                for midx2 in monIdx:
                    if not midx2 in fldSenIdx:
                        nonfldSenIdx.append(midx2)
                if (("STARG_nonflood" in varLst) or ("STARG_flood" in varLst)):
                    origVal1 = lifOrig[28].split(" ")
                    origVal2 = lifOrig[30].split(" ")
                    origVal1[-1] = origVal1[-1][:-1]
                    origVal2[-1] = origVal2[-1][:-1]
                    origVal = origVal1 + origVal2
                    while "" in origVal:
                        origVal.remove("")

                    origVal = list(map(float, origVal))    
                    if (int(parInFile.loc[parInFile["Symbol"] == "STARG_nonflood"]["selectFlag"]) == 1):
                        mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "STARG_nonflood"]["TestVal"])
                        for flidx in nonfldSenIdx:
                            origVal[flidx] = origVal[flidx] * mutiPlier
                    if (int(parInFile.loc[parInFile["Symbol"] == "STARG_flood"]["selectFlag"]) == 1):
                        mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "STARG_flood"]["TestVal"])
                        for flidx in fldSenIdx:
                            origVal[flidx] = origVal[flidx] * mutiPlier
                    
                    lifRes[28] = "{:11.1f}  {:11.1f}  {:11.1f}  {:11.1f}  {:11.1f}  {:11.1f}\n".format(
                        origVal[0], origVal[1], origVal[2], origVal[3], origVal[4], origVal[5]
                    )
                    lifRes[30] = "{:11.1f}  {:11.1f}  {:11.1f}  {:11.1f}  {:11.1f}  {:11.1f}\n".format(
                        origVal[6], origVal[7], origVal[8], origVal[9], origVal[10], origVal[11]
                    )    

        # Then write the contents into the same file
        with open(fnpRes, 'w', encoding="ISO-8859-1") as swatResFile:
            swatResFile.writelines(lifRes)

    return "res"




##########################################################################
def updateParInWwq(parInFile, fdWorkingDir):

    """
    Sometimes the crop is named "crop.dat". Other projects
    used "plant.dat".
    I will try both.
    """
    fnp = os.path.join(fdWorkingDir,
        "basins.wwq")

    try:
        with open(fnp, 'r', encoding="ISO-8859-1") as swatFile:
            lif = swatFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnp, e))
        exit(1)

    # Only values in crop and pasture land are updated.
    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()

    for lidx in range(len(lif)):
        # Line 4 for parameter AI0 
        if ((lidx == 3) and ("AI0" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "AI0"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | AI0 : Ratio of chlorophyll-a to algal biomass [chla/mg algae]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "AI0"]["TestVal"]))
        
        # Line 5 for parameter AI1 
        if ((lidx == 4) and ("AI1" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "AI1"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | AI1 : Fraction of algal biomass that is nitrogen [mg N/mg alg]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "AI1"]["TestVal"]))
        
        # Line 6 for parameter AI2 
        if ((lidx == 5) and ("AI2" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "AI2"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | AI2 : Fraction of algal biomass that is phosphorus [mg P/mg alg]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "AI2"]["TestVal"]))

        # Line 11 for parameter MUMAX 
        if ((lidx == 10) and ("MUMAX" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "MUMAX"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | MUMAX : Maximum specific algal growth rate at 20oC [day-1]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "MUMAX"]["TestVal"]))

        # Line 12 for parameter RHOQ 
        if ((lidx == 11) and ("RHOQ" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RHOQ"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RHOQ : Algal respiration rate at 20oC [day-1]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RHOQ"]["TestVal"]))
        
        # Line 15 for parameter K_N 
        if ((lidx == 14) and ("K_N" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "K_N"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | K_N : Michaelis-Menton half-saturation constant for nitrogen [mg N/lL]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "K_N"]["TestVal"]))
        
        # Line 20 for parameter P_N 
        if ((lidx == 19) and ("P_N" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "P_N"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | P_N : Algal preference factor for ammonia\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "P_N"]["TestVal"]))

    # Then write the contents into the same file
    with open(fnp, 'w', encoding="ISO-8859-1") as swatFile:
        swatFile.writelines(lif)

    return "wwq"



##########################################################################
def updateParInCrop(parInFile, fdWorkingDir, fdmodelTxtInOut):

    """
    Sometimes the crop is named "crop.dat". Other projects
    used "plant.dat".
    I will try both.
    """
    fnpOrig = os.path.join(fdmodelTxtInOut,
        "plant.dat")
    if not os.path.isfile(fnpOrig):
        fnpOrig = os.path.join(fdmodelTxtInOut,
            "crop.dat")
    try:
        with open(fnpOrig, 'r') as swatFileOrig:
            lifOrig = swatFileOrig.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpOrig, e))
        exit(1)

    fnp = os.path.join(fdWorkingDir,
        "plant.dat")
    if not os.path.isfile(fnp):
        fnp = os.path.join(fdWorkingDir,
            "crop.dat")
    try:
        with open(fnp, 'r') as swatFile:
            lif = swatFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnp, e))
        exit(1)

    # Only values in crop and pasture land are updated.
    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    if "USLE_C" in varLst:
        if (int(parInFile.loc[parInFile["Symbol"] == "USLE_C"]["selectFlag"]) == 1):
            usleCVal = float(parInFile.loc[parInFile["Symbol"] == "USLE_C"]["TestVal"])
            for lidx in range(len(lif)):
                # USLE_C need to be modified by fraction 
                # Line 4 for AGRL
                if (lidx in [3, 8, 13, 93, 58, 278, 138]):
                    lifOrig[lidx] = lifOrig[lidx][:-1].split(" ")
                    while "" in lifOrig[lidx]:
                        lifOrig[lidx].remove("")
                    lifOrig[lidx] = list(map(float, lifOrig[lidx])) 

                    lif[lidx] = lif[lidx][:-1].split(" ")
                    while "" in lif[lidx]:
                        lif[lidx].remove("")
                    lif[lidx] = list(map(float, lif[lidx]))  
                    lif[lidx][1] = lifOrig[lidx][1] * (1.0 + usleCVal)
                    lif[lidx] = "".join(["{:8.3f}".format(cropVal) for cropVal in lif[lidx]]) + "\n"


    # Then write the contents into the same file
    with open(fnp, 'w') as swatFile:
        swatFile.writelines(lif)

    return "plant.dat"







##########################################################################
def updateParInBsn(parInFile, fdWorkingDir):

    fnp = os.path.join(fdWorkingDir,
        "basins.bsn")

    try:
        with open(fnp, 'r', encoding="ISO-8859-1") as swatFile:
            lif = swatFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnp, e))
        exit(1)

    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    for lidx in range(len(lif)):
        # Line 4 for parameter SFTMP 
        if ((lidx == 3) and ("SFTMP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SFTMP"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SFTMP : Snowfall temperature\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SFTMP"]["TestVal"]))

        # Line 5 for parameter SMTMP 
        if ((lidx == 4) and ("SMTMP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SMTMP"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SMTMP : Snow melt base temperature\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SMTMP"]["TestVal"]))

        # Line 6 for parameter SMFMX 
        if ((lidx == 5) and ("SMFMX" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SMFMX"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SMFMX : Melt factor for snow on June 21 [mm H2O/oC-day]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SMFMX"]["TestVal"]))

        # Line 7 for parameter SMFMN 
        if ((lidx == 6) and ("SMFMN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SMFMN"]["selectFlag"]) == 1):
                smfMnVal = float(parInFile.loc[parInFile["Symbol"] == "SMFMN"]["TestVal"])
                smfMxVal = float(parInFile.loc[parInFile["Symbol"] == "SMFMX"]["TestVal"])
                if smfMnVal> smfMxVal:
                    smfMnVal = smfMxVal
                lif[lidx] = """{:16.3f}    | SMFMN : Melt factor for snow on December 21 [mm H2O/oC-day]\n""".format(
                    smfMnVal)

        # Line 8 for parameter TIMP 
        if ((lidx == 7) and ("TIMP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "TIMP"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | TIMP : Snow pack temperature lag factor\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "TIMP"]["TestVal"]))

        # Line 9 for parameter SNOCOVMX 
        if ((lidx == 8) and ("SNOCOVMX" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SNOCOVMX"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SNOCOVMX : Minimum snow water content that corresponds to 100% snow cover [mm]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SNOCOVMX"]["TestVal"]))

        # Line 10 for parameter SNO50COV 
        if ((lidx == 9) and ("SNO50COV" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SNO50COV"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SNO50COV : Fraction of snow volume represented by SNOCOVMX that corresponds to 50% snow cover\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SNO50COV"]["TestVal"]))

        # # Line 13 for parameter ESCO 
        # if ((lidx == 9) and ("ESCO" in varLst)):
        #     if (int(parInFile.loc[parInFile["Symbol"] == "ESCO"]["selectFlag"]) == 1):
        #         lif[lidx] = """{:16.3f}    | ESCO :  soil evaporation compensation factor\n""".format(
        #             float(parInFile.loc[parInFile["Symbol"] == "ESCO"]["TestVal"]))

        # Line 14 for parameter EPCO 
        if ((lidx == 13) and ("EPCO" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "EPCO"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | EPCO : plant water uptake compensation factor\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "EPCO"]["TestVal"]))

        # Line 15 for parameter EVLAI 
        if ((lidx == 14) and ("EVLAI" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "EVLAI"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | EVLAI : Leaf area index at which no evaporation occurs from water surface [m2/m2]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "EVLAI"]["TestVal"]))

        # Line 20 for parameter SURLAG 
        if ((lidx == 19) and ("SURLAG" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SURLAG"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SURLAG : Surface runoff lag time [days]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SURLAG"]["TestVal"]))

        # Line 21 for parameter ADJ_PKR 
        if ((lidx == 20) and ("ADJ_PKR" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "ADJ_PKR"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | ADJ_PKR : Peak rate adjustment factor for sediment routing in the subbasin (tributary channels)\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "ADJ_PKR"]["TestVal"]))

        # Line 22 for parameter PRF 
        if ((lidx == 21) and ("PRF" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "PRF"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | PRF_BSN : Peak rate adjustment factor for sediment routing in the main channel\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "PRF"]["TestVal"]))

        # Line 23 for parameter SPCON 
        if ((lidx == 22) and ("SPCON" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SPCON"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SPCON : Linear parameter for calculating the maximum amount of sediment that can be reentrained during channel sediment routing\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SPCON"]["TestVal"]))

        # Line 24 for parameter SPEXP 
        if ((lidx == 23) and ("SPEXP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SPEXP"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SPEXP : Exponent parameter for calculating sediment reentrained in channel sediment routing\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SPEXP"]["TestVal"]))

        # Line 26 for parameter RCN 
        if ((lidx == 25) and ("RCN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RCN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RCN : nitrogen in rainfall (ppm)\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RCN"]["TestVal"]))

        # Line 27 for parameter CMN 
        if ((lidx == 26) and ("CMN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CMN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | CMN : Rate factor for humus mineralization of active organic nitrogen\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CMN"]["TestVal"]))

        # Line 28 for parameter N_UPDIS 
        if ((lidx == 27) and ("N_UPDIS" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "N_UPDIS"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | N_UPDIS : Nitrogen uptake distribution parameter\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "N_UPDIS"]["TestVal"]))

        # Line 30 for parameter NPERCO 
        if ((lidx == 29) and ("NPERCO" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "NPERCO"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | NPERCO : Nitrogen percolation coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "NPERCO"]["TestVal"]))

        # Line 31 for parameter PPERCO 
        if ((lidx == 30) and ("PPERCO" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "PPERCO"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | PPERCO : Phosphorus percolation coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "PPERCO"]["TestVal"]))

        # Line 32 for parameter PHOSKD 
        if ((lidx == 31) and ("PHOSKD" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "PHOSKD"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | PHOSKD : Phosphorus soil partitioning coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "PHOSKD"]["TestVal"]))

        # Line 33 for parameter PSP 
        if ((lidx == 32) and ("PSP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "PSP"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | PSP : Phosphorus sorption coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "PSP"]["TestVal"]))

        # Line 34 for parameter PSP 
        if ((lidx == 33) and ("RSDCO" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RSDCO"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RSDCO : Residue decomposition coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RSDCO"]["TestVal"]))

        # Line 36 for parameter PERCOP 
        if ((lidx == 35) and ("PERCOP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "PERCOP"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | PERCOP : Pesticide percolation coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "PERCOP"]["TestVal"]))

        # Line 59 for parameter MSK_CO1 
        if ((lidx == 58) and ("MSK_CO1" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "MSK_CO1"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | MSK_CO1 : Calibration coefficient used to control impact of the storage time constant (Km) for normal flow \n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "MSK_CO1"]["TestVal"]))

        # Line 60 for parameter MSK_CO2 
        if ((lidx == 59) and ("MSK_CO2" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "MSK_CO2"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | MSK_CO2 : Calibration coefficient used to control impact of the storage time constant (Km) for low flow \n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "MSK_CO2"]["TestVal"]))

        # Line 61 for parameter MSK_X 
        if ((lidx == 60) and ("MSK_X" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "MSK_X"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | MSK_X : Weighting factor controlling relative importance of inflow rate and outflow rate in determining water storage in reach segment\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "MSK_X"]["TestVal"]))

        # Line 66 for parameter EVRCH 
        if ((lidx == 65) and ("EVRCH" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "EVRCH"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | EVRCH : Reach evaporation adjustment factor\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "EVRCH"]["TestVal"]))

        # Line 68 for parameter ICN 
        if ((lidx == 67) and ("ICN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "ICN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | ICN  : Daily curve number calculation method\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "ICN"]["TestVal"]))

        # Line 69 for parameter CNCOEF 
        if ((lidx == 68) and ("CNCOEF" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CNCOEF"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | CNCOEF : Plant ET curve number coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CNCOEF"]["TestVal"]))

        # Line 70 for parameter CDN 
        if ((lidx == 69) and ("CDN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CDN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | CDN : Denitrification exponential rate coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CDN"]["TestVal"]))

        # Line 71 for parameter SDNCO 
        if ((lidx == 70) and ("SDNCO" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SDNCO"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SDNCO : Denitrification threshold water content\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SDNCO"]["TestVal"]))

        # Line 81 for parameter DEPIMP_BSN 
        if ((lidx == 80) and ("DEPIMP_BSN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "DEPIMP_BSN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | DEPIMP_BSN : Depth to impervious layer for modeling perched water tables [mm]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "DEPIMP_BSN"]["TestVal"]))

        # Line 88 for parameter FIXCO 
        if ((lidx == 87) and ("FIXCO" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "FIXCO"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | FIXCO : Nitrogen fixation coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "FIXCO"]["TestVal"]))

        # Line 89 for parameter NFIXMX 
        if ((lidx == 88) and ("NFIXMX" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "NFIXMX"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | NFIXMX : Maximum daily-n fixation [kg/ha]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "NFIXMX"]["TestVal"]))

        # Line 90 for parameter ANION_EXCL_BSN 
        if ((lidx == 89) and ("ANION_EXCL_BSN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "ANION_EXCL_BSN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | ANION_EXCL_BSN : Fraction of porosity from which anions are excluded\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "ANION_EXCL_BSN"]["TestVal"]))



    # Then write the contents into the same file
    with open(fnp, 'w', encoding="ISO-8859-1") as swatFile:
        swatFile.writelines(lif)

    return "bsn"






##########################################################################
def updateParInChm(fnSwatHruLvl, parInFile, fdWorkingDir):

    fnpSol = os.path.join(fdWorkingDir,
        "{}.sol".format(fnSwatHruLvl))

    fnpChm = os.path.join(fdWorkingDir,
        "{}.chm".format(fnSwatHruLvl))

    try:
        with open(fnpSol, 'r') as solFile:
            lifSol = solFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpSol, e))
        exit(1)

    try:
        with open(fnpChm, 'r') as chmFile:
            lif = chmFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpChm, e))
        exit(1)

    # Get the soil depth (line 8) and soil Organic N (line 12) from the sol file
    solDepLst = lifSol[7].split(":")[1][:-1].split(" ")
    while "" in solDepLst:
        solDepLst.remove("")
    solDepLst = list(map(float, solDepLst))
    avgSolDep = sum(solDepLst)/float(len(solDepLst))

    solOCLst = lifSol[11].split(":")[1][:-1].split(" ")
    while "" in solOCLst:
        solOCLst.remove("")
    solOCLst = list(map(float, solOCLst))

    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    for lidx in range(len(lif)):

        # Line 4 for parameter SOLN 
        # In the matlab code, the SOLN value was scaled by the following equation:
        # str_SolN = [ones(1,n_layers).*exp(-avg_depth/1000) zeros(1,10-length(sol_depth))];
        # SOLN = x(cellfun(@(x) isequal(x, 'SOLN'), symbol))*str_SolN;
        # We will keep it as it is in the matlab code.
        if ((lidx == 3) and ("SOLN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SOLN"]["selectFlag"]) == 1):
                preText = lif[lidx].split(":")[0]
                solValues = lif[lidx].split(":")[1][:-1].split(" ")
                while "" in solValues:
                    solValues.remove("")
                solValNoDataInLyr = solValues[len(solDepLst):]
                newSolNVal = float(parInFile.loc[parInFile["Symbol"] == "SOLN"]["TestVal"])
                solValHasDataInLyr = [newSolNVal * math.exp(-avgSolDep/1000)] * len(solDepLst)
                wholeLst = solValHasDataInLyr + list(map(float,solValNoDataInLyr))
                newValPartinLine = "".join(["{:12.2f}".format(valLayer) for valLayer in wholeLst])
                lif[lidx] = """{}:{}\n""".format(preText, newValPartinLine)

        # Line 5 for parameter OGRN 
        if ((lidx == 4) and ("ORGN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "ORGN"]["selectFlag"]) == 1):
                preText = lif[lidx].split(":")[0]
                solValues = lif[lidx].split(":")[1][:-1].split(" ")
                while "" in solValues:
                    solValues.remove("")
                solValNoDataInLyr = solValues[len(solDepLst):]
                newSolNVal = float(parInFile.loc[parInFile["Symbol"] == "ORGN"]["TestVal"])
                # Soil Organic N depends on the soil organic carbon content in the first layer
                # print(solOCLst)
                if not solOCLst[0] == 0.00:
                    mutiPlier = [soc/solOCLst[0] for soc in solOCLst]
                else:
                    mutiPlier = [0.00] * len(solOCLst)
                solValHasDataInLyr = [newSolNVal * mutip for mutip in mutiPlier]
                wholeLst = solValHasDataInLyr + list(map(float,solValNoDataInLyr))
                newValPartinLine = "".join(["{:12.2f}".format(valLayer) for valLayer in wholeLst])
                lif[lidx] = """{}:{}\n""".format(preText, newValPartinLine)

        # Line 6 for parameter LABP 
        if ((lidx == 5) and ("LABP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "LABP"]["selectFlag"]) == 1):
                preText = lif[lidx].split(":")[0]
                solValues = lif[lidx].split(":")[1][:-1].split(" ")
                while "" in solValues:
                    solValues.remove("")
                solValNoDataInLyr = solValues[len(solDepLst):]
                newSolNVal = float(parInFile.loc[parInFile["Symbol"] == "LABP"]["TestVal"])
                # Soil Organic N depends on the soil organic carbon content in the first layer
                solValHasDataInLyr =  [newSolNVal] * len(solDepLst)
                wholeLst = solValHasDataInLyr + list(map(float,solValNoDataInLyr))
                newValPartinLine = "".join(["{:12.2f}".format(valLayer) for valLayer in wholeLst])
                lif[lidx] = """{}:{}\n""".format(preText, newValPartinLine)

        # Line 7 for parameter ORGP 
        if ((lidx == 6) and ("ORGP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "ORGP"]["selectFlag"]) == 1):
                preText = lif[lidx].split(":")[0]
                solValues = lif[lidx].split(":")[1][:-1].split(" ")
                while "" in solValues:
                    solValues.remove("")
                solValNoDataInLyr = solValues[len(solDepLst):]
                newSolNVal = float(parInFile.loc[parInFile["Symbol"] == "ORGP"]["TestVal"])
                # Soil Organic N depends on the soil organic carbon content in the first layer
                solValHasDataInLyr =  [newSolNVal] * len(solDepLst)
                wholeLst = solValHasDataInLyr + list(map(float,solValNoDataInLyr))
                newValPartinLine = "".join(["{:12.2f}".format(valLayer) for valLayer in wholeLst])
                lif[lidx] = """{}:{}\n""".format(preText, newValPartinLine)

    # Then write the contents into the same file
    with open(fnpChm, 'w') as swatFile:
        swatFile.writelines(lif)

    return "chm"






##########################################################################
def updateParInSol(fnSwatHruLvl, parInFile, fdWorkingDir, fdmodelTxtInOut):

    fnpOrig = os.path.join(fdmodelTxtInOut,
        "{}.sol".format(fnSwatHruLvl))

    try:
        with open(fnpOrig, 'r') as swatFileOrig:
            lifOrig = swatFileOrig.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpOrig, e))
        exit(1)



    fnp = os.path.join(fdWorkingDir,
        "{}.sol".format(fnSwatHruLvl))

    try:
        with open(fnp, 'r') as swatFile:
            lif = swatFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnp, e))
        exit(1)

    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    for lidx in range(len(lif)):

        # Line 8 for parameter SOL_Z 
        if ((lidx == 7) and ("SOL_Z" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SOL_Z"]["selectFlag"]) == 1):
                preText = lifOrig[lidx].split(":")[0]
                solValues = lifOrig[lidx].split(":")[1][:-1].split(" ")
                while "" in solValues:
                    solValues.remove("")
                origValLst = list(map(float, solValues))
                mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "SOL_Z"]["TestVal"])
                newValLst = "".join(["{:12.3f}".format(origV * mutiPlier) for origV in origValLst])
                lif[lidx] = """{}:{}\n""".format(preText, newValLst)

        # Line 10 for parameter SOL_AWC 
        if ((lidx == 9) and ("SOL_AWC" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SOL_AWC"]["selectFlag"]) == 1):
                preText = lifOrig[lidx].split(":")[0]
                solValues = lifOrig[lidx].split(":")[1][:-1].split(" ")
                while "" in solValues:
                    solValues.remove("")
                origValLst = list(map(float, solValues))
                mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "SOL_AWC"]["TestVal"])
                newValLst = "".join(["{:12.1f}".format(origV * mutiPlier) for origV in origValLst])
                lif[lidx] = """{}:{}\n""".format(preText, newValLst)

        # Line 11 for parameter SOL_K 
        if ((lidx == 10) and ("SOL_K" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SOL_K"]["selectFlag"]) == 1):
                preText = lifOrig[lidx].split(":")[0]
                solValues = lifOrig[lidx].split(":")[1][:-1].split(" ")
                while "" in solValues:
                    solValues.remove("")
                origValLst = list(map(float, solValues))
                mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "SOL_K"]["TestVal"])
                newValLst = "".join(["{:12.1f}".format(origV * mutiPlier) for origV in origValLst])
                lif[lidx] = """{}:{}\n""".format(preText, newValLst)

        # Line 17 for parameter SOL_ALB 
        if ((lidx == 16) and ("SOL_ALB" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SOL_ALB"]["selectFlag"]) == 1):
                preText = lifOrig[lidx].split(":")[0]
                solValues = lifOrig[lidx].split(":")[1][:-1].split(" ")
                while "" in solValues:
                    solValues.remove("")
                origValLst = list(map(float, solValues))
                mutiPlier = 1.00+float(parInFile.loc[parInFile["Symbol"] == "SOL_ALB"]["TestVal"])
                newValLst = "".join(["{:12.3f}".format(origV * mutiPlier) for origV in origValLst])
                lif[lidx] = """{}:{}\n""".format(preText, newValLst)

        # Line 18 for parameter USLE_K 
        if ((lidx == 17) and ("USLE_K" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "USLE_K"]["selectFlag"]) == 1):
                preText = lifOrig[lidx].split(":")[0]
                solValues = lifOrig[lidx].split(":")[1][:-1].split(" ")
                while "" in solValues:
                    solValues.remove("")
                origValLst = list(map(float, solValues))
                mutiPlier = 1+float(parInFile.loc[parInFile["Symbol"] == "USLE_K"]["TestVal"])
                newValLst = "".join(["{:12.3f}".format(origV * mutiPlier) for origV in origValLst])
                lif[lidx] = """{}:{}\n""".format(preText, newValLst)

    # Then write the contents into the same file
    with open(fnp, 'w') as swatFile:
        swatFile.writelines(lif)

    return "sol"




##########################################################################
def updateParInMgt(fnSwatHruLvl, parInFile, fdWorkingDir, rowCropLst, fdmodelTxtInOut):

    """
    For mgt, information from HRU and SOL files are required.
    The required information in the HRU file is land type.
    The required information in the SOL file is the hydrologic soil group. 
    """
    # First readin the contents of the old file
    fnpHru = os.path.join(fdWorkingDir,
        "{}.hru".format(fnSwatHruLvl))

    fnpMgt = os.path.join(fdWorkingDir,
        "{}.mgt".format(fnSwatHruLvl))

    fnpSol = os.path.join(fdWorkingDir,
        "{}.sol".format(fnSwatHruLvl))

    fnpMgtOrig = os.path.join(fdmodelTxtInOut,
        "{}.mgt".format(fnSwatHruLvl))

    try:
        with open(fnpHru, 'r') as hruFile:
            lifHru = hruFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpHru, e))
        exit(1)


    try:
        with open(fnpSol, 'r') as solFile:
            lifSol = solFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpSol, e))
        exit(1)

    try:
        with open(fnpMgtOrig, 'r') as mgtFileOrig:
            lifmgtOrig = mgtFileOrig.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpMgtOrig, e))
        exit(1)

    try:
        with open(fnpMgt, 'r') as mgtFile:
            lif = mgtFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpMgt, e))
        exit(1)

    # First get the required information from the lines in the HRU and SOL file
    isRowCrops = False
    for rcIdx in rowCropLst:
        if lifHru[0].find(rcIdx):
            isRowCrops = True
            break
    
    # Get the soil HSG from soil file
    is_soilHSG_BCD = False
    soilHSG = lifSol[2].split(":")[1].replace(" ", "") 
    if soilHSG in ["B", "C", "D"]:
        is_soilHSG_BCD = True

    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    for lidx in range(len(lif)):
        # Line 10 for parameter BIOMIX 
        if ((lidx == 9) and ("BIOMIX" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "BIOMIX"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | BIOMIX: Biological mixing efficiency\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "BIOMIX"]["TestVal"]))
        
        # Line 11 for parameter CN_F 
        if ((lidx == 10) and ("CN_F" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CN_F"]["selectFlag"]) == 1):
                origVal = float(lifmgtOrig[lidx].split("|")[0])
                newVal = origVal * (1+float(parInFile.loc[parInFile["Symbol"] == "CN_F"]["TestVal"]))
                lif[lidx] = """{:16.3f}    | CN2: Initial SCS CN II value\n""".format(newVal)

        # Line 12 for parameter USLE_P 
        if ((lidx == 11) and ("USLE_P" in varLst) and isRowCrops):
            if (int(parInFile.loc[parInFile["Symbol"] == "USLE_P"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | USLE_P: USLE support practice factor\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "USLE_P"]["TestVal"]))

        # Line 13 for parameter BIOMIN 
        if ((lidx == 12) and ("BIOMIN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "BIOMIN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | BIOMIN: Minimum biomass for grazing (kg/ha)\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "BIOMIN"]["TestVal"]))

        # Line 14 for parameter FILTERW 
        if ((lidx == 13) and ("FILTERW" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "FILTERW"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | FILTERW: width of edge of field filter strip (m)\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "FILTERW"]["TestVal"]))

        # Line 25 for parameter FILTERW 
        if ((lidx == 24) and ("DDRAIN" in varLst) and is_soilHSG_BCD):
            if (int(parInFile.loc[parInFile["Symbol"] == "DDRAIN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | DDRAIN: depth to subsurface tile drain (mm)\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "DDRAIN"]["TestVal"]))

        # Line 26 for parameter TDRAIN 
        if ((lidx == 25) and ("TDRAIN" in varLst) and is_soilHSG_BCD):
            if (int(parInFile.loc[parInFile["Symbol"] == "TDRAIN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | TDRAIN: time to drain soil to field capacity (hr)\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "TDRAIN"]["TestVal"]))

        # Line 27 for parameter GDRAIN 
        if ((lidx == 26) and ("GDRAIN" in varLst) and is_soilHSG_BCD):
            if (int(parInFile.loc[parInFile["Symbol"] == "GDRAIN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | GDRAIN: drain tile lag time (hr)\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "GDRAIN"]["TestVal"]))

    # Then write the contents into the same file
    with open(fnpMgt, 'w') as swatFile:
        swatFile.writelines(lif)

    return "mgt"


##########################################################################
def updateParInHru(fnSwatHruLvl, parInFile, fdWorkingDir, fdmodelTxtInOut):

    # First readin the contents of the old file
    fnpOrig = os.path.join(fdmodelTxtInOut,
        "{}.hru".format(fnSwatHruLvl))
    try:
        with open(fnpOrig, 'r') as swatFileOrig:
            lifOrig = swatFileOrig.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpOrig, e))
        exit(1)

    # First readin the contents of the old file
    fnp = os.path.join(fdWorkingDir,
        "{}.hru".format(fnSwatHruLvl))
    try:
        with open(fnp, 'r') as swatFile:
            lif = swatFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnp, e))
        exit(1)

    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    for lidx in range(len(lif)):
        # Line 3 for parameter SLSUBBSN 
        if ((lidx == 2) and ("SLSUBBSN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SLSUBBSN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SLSUBBSN : Average slope length [m]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SLSUBBSN"]["TestVal"]))

        # Line 4 for parameter SLOPE
        # Slope need to be modified by fraction 
        if ((lidx == 3) and ("SLOPE" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SLOPE"]["selectFlag"]) == 1):
                origVal = float(lifOrig[lidx].split("|")[0])
                newVal = origVal * (1+float(parInFile.loc[parInFile["Symbol"] == "SLOPE"]["TestVal"]))
                lif[lidx] = """{:16.3f}    | SLOPE : Main channel slope [m/m]\n""".format(newVal)
        
        # Line 5 for parameter OV_N 
        if ((lidx == 4) and ("OV_N" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "OV_N"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | OV_N : Manning"s "n" value for overland flow\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "OV_N"]["TestVal"]))

        # Line 9 for parameter CANMX 
        if ((lidx == 8) and ("CANMX" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CANMX"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | CANMX : Maximum canopy storage [mm]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CANMX"]["TestVal"]))

        # Line 10 for parameter ESCO 
        if ((lidx == 9) and ("ESCO" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "ESCO"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | ESCO : Soil evaporation compensation factor\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "ESCO"]["TestVal"]))

        # Line 12 for parameter RSDIN 
        if ((lidx == 11) and ("RSDIN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RSDIN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RSDIN : Initial residue cover [kg/ha]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RSDIN"]["TestVal"]))

        # Line 24 for parameter DEP_IMP 
        if ((lidx == 23) and ("DEP_IMP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "DEP_IMP"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | DEP_IMP : Depth to impervious layer in soil profile [mm]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "DEP_IMP"]["TestVal"]))

    # Then write the contents into the same file
    with open(fnp, 'w') as swatFile:
        swatFile.writelines(lif)

    return "hru"







##########################################################################
def updateParInGw(fnSwatHruLvl, parInFile, fdWorkingDir, fdmodelTxtInOut):

    # First readin the contents of the old file
    fnpOrig = os.path.join(fdmodelTxtInOut,
        "{}.gw".format(fnSwatHruLvl))
    try:
        with open(fnpOrig, 'r') as swatFileOrig:
            lifOrig = swatFileOrig.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpOrig, e))
        exit(1)
    
    # First readin the contents of the old file
    fnp = os.path.join(fdWorkingDir,
        "{}.gw".format(fnSwatHruLvl))
    try:
        with open(fnp, 'r') as swatFile:
            lif = swatFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnp, e))
        exit(1)

    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    for lidx in range(len(lif)):
        #  Line 4 for parameter GW_DELAY 
        if ((lidx == 3) and ("GW_DELAY" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "GW_DELAY"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | GW_DELAY : Groundwater delay [days]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "GW_DELAY"]["TestVal"]))

        # Line 5 for parameter ALPHA_BF 
        if ((lidx == 4) and ("ALPHA_BF" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "ALPHA_BF"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | ALPHA_BF : BAseflow alpha factor [days]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "ALPHA_BF"]["TestVal"]))

        # Line 6 for parameter GWQMN 
        if ((lidx == 5) and ("GWQMN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "GWQMN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | GWQMN : Threshold depth of water in the shallow aquifer required for return flow to occur [mm]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "GWQMN"]["TestVal"]))

        # Line 7 for parameter GW_REVAP 
        if ((lidx == 6) and ("GW_REVAP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "GW_REVAP"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | GW_REVAP : Groundwater "revap" coefficient\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "GW_REVAP"]["TestVal"]))

        # Line 8 for parameter REVEP_MN 
        if ((lidx == 7) and ("REVEP_MN" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "REVEP_MN"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | REVEP_MN : Threshold depth of water in the shallow aquifer for "revap" to occur [mm]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "REVEP_MN"]["TestVal"]))
        
        # Line 9 for parameter RCHRG_DP 
        if ((lidx == 8) and ("RCHRG_DP" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RCHRG_DP"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RCHRG_DP : Deep aquifer percolation fraction\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RCHRG_DP"]["TestVal"]))

        # Line 10 for parameter GWHT 
        if ((lidx == 9) and ("GWHT" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "GWHT"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | GWHT : Initial groundwater height [m]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "GWHT"]["TestVal"]))

        # Line 11 for parameter GW_SPYLD 
        if ((lidx == 10) and ("GW_SPYLD" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "GW_SPYLD"]["selectFlag"]) == 1):
                # Determine whether this is modified by fraction or absolute value
                origVal = float(lifOrig[lidx].split("|")[0])
                newVal = origVal * (1+float(parInFile.loc[parInFile["Symbol"] == "GW_SPYLD"]["TestVal"]))
                lif[lidx] = """{:16.3f}    | GW_SPYLD : Specific yield of the shallow aquifer [m3/m3]\n""".format(newVal)
        
        # Line 12 for parameter SHALLST_N 
        if ((lidx == 11) and ("SHALLST_N" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "SHALLST_N"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | SHALLST_N : Initial concentration of nitrate in shallow aquifer [mg N/l]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "SHALLST_N"]["TestVal"]))

        # Line 14 for parameter HLIFE_NGW 
        if ((lidx == 13) and ("HLIFE_NGW" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "HLIFE_NGW"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | HLIFE_NGW : Half-life of nitrate in the shallow aquifer [day]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "HLIFE_NGW"]["TestVal"]))

    # Then write the contents into the same file
    with open(fnp, 'w') as swatFile:
        swatFile.writelines(lif)

    return "gw"




##########################################################################
def updateParInSwq(fnSwatSubLvl, parInFile, fdWorkingDir):

    # First readin the contents of the old file
    
    fnpSwatSwq = os.path.join(fdWorkingDir,
        "{}.swq".format(fnSwatSubLvl))

    try:
        with open(fnpSwatSwq, 'r', encoding="ISO-8859-1") as swqFile:
            lif = swqFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpSwatSwq, e))
        exit(1)

    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    for lidx in range(len(lif)):
        #  Line 3 for parameter CH_NII 
        if ((lidx == 2) and ("RS1" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RS1"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RS1 : Local algal settling rate in the reach at 20 [m/day]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RS1"]["TestVal"]))

        # Line 4 for parameter RS2 
        if ((lidx == 3) and ("RS2" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RS2"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RS2 : Benthic (sediment) source rate for dissolved phosphorus in the reach at 20 [m/day]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RS2"]["TestVal"]))

        # Line 5 for parameter RS3 
        if ((lidx == 4) and ("RS3" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RS3"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RS3 : Benthic source rate for NH4-N in the reach at 20 [mg NH4-N/[m2ay]]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RS3"]["TestVal"]))

        # Line 6 for parameter RS4 
        if ((lidx == 5) and ("RS4" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RS4"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RS4 : Rate coefficient for organic N settling in the reach at 20 [day-1]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RS4"]["TestVal"]))

        # Line 7 for parameter RS5 
        if ((lidx == 6) and ("RS5" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "RS5"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | RS5 : Organic phosphorus settling rate in the reach at 20 [day-1]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "RS5"]["TestVal"]))
        
        # Line 16 for parameter BC1 
        if ((lidx == 15) and ("BC1" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "BC1"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | BC1 : Rate constant for biological oxidation of NH4 to NO2 in the reach at 20C [day-1]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "BC1"]["TestVal"]))

        # Line 17 for parameter BC2 
        if ((lidx == 16) and ("BC2" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "BC2"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | BC2 : Rate constant for biological oxidation of NO2 to NO3 in the reach at 20C [day-1]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "BC2"]["TestVal"]))

        # Line 18 for parameter BC3 
        if ((lidx == 17) and ("BC3" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "BC3"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | BC3 : Rate constant for hydrolysis of organic N to NH4 in the reach at 20C [day-1]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "BC3"]["TestVal"]))

        # Line 19 for parameter BC4 
        if ((lidx == 18) and ("BC4" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "BC4"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | BC4 : Channel cover factor\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "BC4"]["TestVal"]))

    # Then write the contents into the same file
    with open(fnpSwatSwq, 'w', encoding="ISO-8859-1") as swqFile:
        swqFile.writelines(lif)

    return "swq"







##########################################################################
def updateParInRte(fnSwatSubLvl, parInFile, fdWorkingDir, fdmodelTxtInOut):

    # First readin the contents of the original file in the swattio folder
    fnpSwatRteOrig = os.path.join(fdmodelTxtInOut,
        "{}.rte".format(fnSwatSubLvl))
    try:
        with open(fnpSwatRteOrig, 'r', encoding="ISO-8859-1") as rteFileOrig:
            lifOrig = rteFileOrig.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpSwatRteOrig, e))
        exit(1)

    # First readin the contents of the new file in the working folder
    fnpSwatRte = os.path.join(fdWorkingDir,
        "{}.rte".format(fnSwatSubLvl))
    try:
        with open(fnpSwatRte, 'r', encoding="ISO-8859-1") as rteFile:
            lif = rteFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpSwatRte, e))
        exit(1)

    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    for lidx in range(len(lif)):
        # Line 4 for parameter CH_SII
        # In the matlab code, the CH_S2 was modified using the following code:
        # (in matlab): CH_S2=str2double(strtok(line))*(1+CH_SII);
        # This means, the script get the current value of slope and time it by percent.
        if ((lidx == 3) and ("CH_SII" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CH_SII"]["selectFlag"]) == 1):
                # Determine whether this is modified by fraction or absolute value
                origVal = float(lifOrig[lidx].split("|")[0])
                newVal = origVal * (1+float(parInFile.loc[parInFile["Symbol"] == "CH_SII"]["TestVal"]))
                lif[lidx] = """{:14.3f}    | CH_SII : Main channel slope [m/m]\n""".format(newVal)
        
        #  Line 6 for parameter CH_NII 
        if ((lidx == 5) and ("CH_NII" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CH_NII"]["selectFlag"]) == 1):
                lif[lidx] = """{:14.3f}    | CH_NII : Manning"s nvalue for main channel\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CH_NII"]["TestVal"]))

        # Line 7 for parameter CH_KII 
        if ((lidx == 6) and ("CH_KII" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CH_KII"]["selectFlag"]) == 1):
                lif[lidx] = """{:14.3f}    | CH_KII : Effective hydraulic conductivity [mm/hr]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CH_KII"]["TestVal"]))

        # Line 8 for parameter CH_COV1 
        if ((lidx == 7) and ("CH_COV1" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CH_COV1"]["selectFlag"]) == 1):
                lif[lidx] = """{:14.3f}    | CH_COV1 : Channel erodibility factor\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CH_COV1"]["TestVal"]))

        # Line 9 for parameter CH_COV2 
        if ((lidx == 8) and ("CH_COV2" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CH_COV2"]["selectFlag"]) == 1):
                lif[lidx] = """{:14.3f}    | CH_COV2 : Channel cover factor\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CH_COV2"]["TestVal"]))


    # Then write the contents into the same file
    with open(fnpSwatRte, 'w', encoding="ISO-8859-1") as rteFile:
        rteFile.writelines(lif)

    return "rte"




##########################################################################
def updateParInSub(fnSwatSubLvl, parInFile, fdWorkingDir):

    # First readin the contents of the old file
    fnpSwatSub = os.path.join(fdWorkingDir,
        "{}.sub".format(fnSwatSubLvl))
    try:
        with open(fnpSwatSub, 'r', encoding="ISO-8859-1") as subFile:
            lif = subFile.readlines()
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnpSwatSub, e))
        exit(1)
    
    # Get a list of variables selected. Since we can get in this function,
    # at least one variable in the sub file was selected. 
    varLst = parInFile["Symbol"].unique()
    for lidx in range(len(lif)):
        # Line 26 for parameter CH_S1
        if ((lidx == 25) and ("CH_SI" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CH_SI"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | CH_SI : Average slope of tributary channel [m/m]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CH_SI"]["TestVal"]))
        
        #  Line 28 for parameter CH_K1 
        if ((lidx == 27) and ("CH_KI" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CH_KI"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | CH_KI : Effective hydraulic conductivity in tributary channel [mm/hr]\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CH_KI"]["TestVal"]))

        # Line 29 for parameter CH_N1 
        if ((lidx == 28) and ("CH_NI" in varLst)):
            if (int(parInFile.loc[parInFile["Symbol"] == "CH_NI"]["selectFlag"]) == 1):
                lif[lidx] = """{:16.3f}    | CH_NI : Manning"s "n" value for the tributary channels\n""".format(
                    float(parInFile.loc[parInFile["Symbol"] == "CH_NI"]["TestVal"]))


    # Then write the contents into the same file
    with open(fnpSwatSub, 'w') as subFile:
        subFile.writelines(lif)

    return "sub"





##########################################################################
def copySWATFileToWD(fnSwatFileSrc, destDir):
   
    fnSrc = os.path.split(fnSwatFileSrc)[1]
    ftoDest = os.path.join(destDir, fnSrc)

    try:
        if os.path.isfile(ftoDest):
            os.remove(ftoDest)
        copyfile(fnSwatFileSrc, ftoDest)
        
    except IOError as e:
        print("Unable to copy file. {} due to: {}".format(fnSrc, e))
        exit(1)
    
    return "copySWATFiles"


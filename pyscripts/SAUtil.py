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
from SALib.test_functions import Ishigami
import numpy as np
import math

from .globVars import useRay
from .DMPOTUtil import getRch2DF

if useRay: 
    from .RaySWATUtil import *
else:
    from .SWATUtil import *


##########################################################################
# Define functions #######################################################
##########################################################################

##########################################################################
def geneParmSetsForSA(parmBsnLvl, parmSubLvl):

    """
    This function takes the basin and subarea level parameter selected
    by users and generate the problem varaible required for SALib.
    Input: 
    parmBsnLvl, parmSubLvl (Format: Dataframe)
    Output：
    Dictionary
    problem = {
    'num_vars': 3,
    'names': ['x1', 'x2', 'x3'],
    'bounds': [[-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359]]
    }

    """
    
    parSymbolList = parmBsnLvl["Symbol"].tolist() + parmSubLvl["Symbol"].tolist()
    parBoundsArray = parmBsnLvl[["LowerBound","UpperBound"]].values.tolist() + parmSubLvl[["LowerBound","UpperBound"]].values.tolist()

    problem = {}
    problem['num_vars'] = len(parSymbolList)
    problem['names'] = parSymbolList
    problem['bounds'] = parBoundsArray
    
    return problem


def checkSaltelliArgument(usrSaltelli):
    # Python3 Program to find whether a no is power of two
    # Function to check
    # if x is power of 2
    log2ofN = (math.log10(usrSaltelli) / math.log10(2))
    if (math.ceil(log2ofN) == math.floor(log2ofN)):
        return True
    else:
        return False
    # This code is contributed
    # by mits


def generateParmValDict(parmSymbol, parmValList):
    
    """
    This function combines the symbol and value of each sample run
    into a dictionary.
    """
    outDict = {}

    for sbIdx in range(len(parmSymbol)):
        outDict[parmSymbol[sbIdx]] = parmValList[sbIdx]

    return outDict



def updateParmInDf(parmDF, parmValDict):
    
    """
    This function updates the parameter values to the sample values.
    """
    for key, val in parmValDict.items():
        parmDF.loc[parmDF['Symbol'] == key, "TestVal"] = val

    return parmDF


def extractAvgAnnEachGroup(runningDir, 
            iPrintForCio, 
            outLetList, 
            outputVarList,
            outLetAvgAnnList,
            rcvRchLst,
            varIDObsHdrPair,
            saRunIdx):
    
    """
    This function extracts the rch file and calculate the average annual values.
    """
    fnRch = os.path.join(runningDir, "output.rch")
    try:
        rchDFWhole = getRch2DF(fnRch, iPrintForCio, len(rcvRchLst))
    except IOError as e:
        print("File {} does not exist: {}. Please double check your TxtInOut \
            folder and make sure you have a complete set".format(fnRch, e))
        exit(1)

    # The outlet no and variable list need to be correspond.
    for oltNoIdx in range(len(outLetList)):
        outLetNo = outLetList[oltNoIdx]
        outletVarID = outputVarList[oltNoIdx]
        outVarHdrNm = varIDObsHdrPair[outletVarID]
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
        simOutVarDF = rchDFWhole.loc[rchDFWhole["RCH"] == outLetNo][simOutVarHeader]
        outLetAvgAnnList[outLetNo][saRunIdx] = simOutVarDF[outVarHdrNm].mean()

    return outLetAvgAnnList
    





# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14, 2021
This class is designed to be a collection of functions scenarios.

@author: Qingyu.Feng
"""

##########################################################################
# Import modules #########################################################
##########################################################################

import os
import math
from shutil import copyfile

##########################################################################
# Define functions #######################################################
##########################################################################
##########################################################################

##########################################################################
def createScenList(resParmList):

    scenarioDict = {}
    icount = 1
    for ndidx in resParmList["NDTARGR"]:
        for if1 in resParmList["IFLOD1R"]:
            for if2 in resParmList["IFLOD2R"]:
                scKey = "sce{}".format(icount)
                scenarioDict[scKey] = {}
                if not if1 == if2:
                    scenarioDict[scKey]["NDTARGR"] = ndidx
                    scenarioDict[scKey]["IFLOD1R"] = if1
                    scenarioDict[scKey]["IFLOD2R"] = if2
                    icount = icount + 1

    return scenarioDict



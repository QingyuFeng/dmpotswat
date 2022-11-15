#!/usr/bin/python3

"""
This program was developed to convert the output of obj function
into excel and help the selection of best runs.


"""

##########################################################################
# Import modules #########################################################
##########################################################################

import pandas, numpy
import datetime
import pytz

from pyscripts.globVars import *
from pyscripts.DMPOTUtil import ctrlSetToJSON
import matplotlib.pyplot as plt

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

outLetList = ctrlSetting["outLetList"]

objFunNmDict = dict()
objFunNmDict[1] = "NSE"
objFunNmDict[2] = "BIAS"
objFunNmDict[3] = "RMSE"
objFunNmDict[4] = "R2"
objFunNmDict[5] = "MSE"


fig = None
# The unit of figsize is in units.
# dpi: dots per inch
fig, axes = plt.subplots(figsize=(4, 
                                  4), 
                        dpi=400,
                        tight_layout=True)
fnFigAllObj = os.path.join(fdOutputs, "objValProgressAllOlts.jpg")

for oltidx in range(len(outLetList)):

    oltNo = outLetList[oltidx]
    fnOFOut = os.path.join(fdOutputs, "DMPOT_ObjFun{}.out".format(oltNo))
    oltObjValList = pandas.read_table(fnOFOut, header=0, sep=",") 
    # Sort the dataframe by object function
    oltObjNm = objFunNmDict[ctrlSetting["objFuncLst"][oltidx]]
    # oltObjValList = oltObjValList.sort_values(by=[oltObjNm], 
    #                                     ascending = False)
    totalRunNo = ctrlSetting["totalModelRuns"]
    fnOFOutCsv = os.path.join(fdOutputs, "DMPOT_ObjFun{}_{}runs.csv".format(oltNo, totalRunNo))
    oltObjValList.to_csv(fnOFOutCsv, 
                            encoding="utf-8",
                            index = False)
    
    axes.scatter(range(0, ctrlSetting["totalModelRuns"]),
        oltObjValList["BestOF"],
        marker = ".",
        s = 2,
        label="Outlet {}".format(oltNo))
    # Control legend
    axes.legend(title = "{}".format(oltObjNm),
                fontsize = 8,
                title_fontsize = 8,
                #ncol = 2,
                framealpha = 0.5, 
                loc = 'upper right',
                edgecolor = "white" 
                )
            
    # box = ax1.get_position()
    # # setposition(left, bottom, width, height)
    # ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Control label
    # set: a property batch setter
    # set_xlabel(xlabel, labelpad, **kwargs)
    axes.set_xlabel("Number of Runs",
                fontsize=8)
    axes.set_ylabel("Best Objective Function Value", 
                    fontsize=8)
    
    # Control grids
    axes.grid(which='major',
              linewidth=0.3)
            
    # Control ticks
    axes.set_xlim(left=-5, right=totalRunNo+5)
    axes.set_xticks(numpy.arange(0, totalRunNo+1, totalRunNo/4))
    axes.set_ylim(bottom=0, top=1)
    axes.set_yticks(numpy.arange(0, 1.1, 0.25))
    
    axes.tick_params(
        which='major',
        direction = "in",
        labelsize=8)
    # Updated Nov 5, 2021
    # In order to facilitate the selection of best parameters, a figure
    # displaying the progress of improvement is created based on the 
    # objective function values at all runs for all outlets.

fig.savefig(fnFigAllObj, bbox_inches="tight")


# End of for loop for total runs
print("--------------------------------------")
print("Congratulations!!! It's done nicely~~~")
print("--------------------------------------")

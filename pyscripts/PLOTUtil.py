# -*- coding: utf-8 -*-
"""
Created on Tue May 12 2020 

This module contain functions to generate figures

@author: Qingyu.Feng
"""

##########################################################################
# Import modules #########################################################
##########################################################################
import pandas, numpy
import os

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

##########################################################################
# Define classes #########################################################
##########################################################################
       


##########################################################################
def genPDFLinePlotObsSim(fdOutputs, 
                         parmObjFnKeys,
                         obsSimPairKeysSplit,
                         obsSimPair, nCalVal, groupSubareaIdx):
    """
    This function generate line charts for observed vs simulated flow and 
    save them as pdf files.
"""       
    if groupSubareaIdx == 1:
        nlumpdist = "Distributed"
    elif groupSubareaIdx == 0:
        nlumpdist = "Lumped"
    # Generate combined figures
    fnpPlotPdf = "{}/ObsVsSim_UserBest_{}_{}.pdf".format(fdOutputs, nCalVal, nlumpdist)
    if os.path.isfile(fnpPlotPdf):
        os.remove(fnpPlotPdf)
    outPlotPdf = PdfPages(fnpPlotPdf)
    
    # Calculating objective functions and making plots
    for subGPKey in parmObjFnKeys:
        # First calculate the statistics with outlets
        if not subGPKey == "Other":
            # Generateing plots for each outlet
            obsSimThisKey = obsSimPairKeysSplit[subGPKey]
            plt.plot(obsSimPair[obsSimThisKey][0], label="Observation")
            plt.plot(obsSimPair[obsSimThisKey][1], label="Simulation")
            plt.legend(loc="upper right",fontsize='x-large')


            plt.title("Station {} {} {}".format(subGPKey, nlumpdist, nCalVal))
            outPlotPdf.savefig()
            plt.close('all')

    outPlotPdf.close()
    
    
##########################################################################
def genPNGLinePlotObsSim(fdOutputs, 
                         parmObjFnKeys,
                         obsSimPairKeysSplit,
                         obsSimPair,nCalVal, groupSubareaIdx):
    """
    This function generate line charts for observed vs simulated flow and 
    save them as individual PNG files.
    """

    if groupSubareaIdx == 1:
        nlumpdist = "Distributed"
    elif groupSubareaIdx == 0:
        nlumpdist = "Lumped"
    # Calculating objective functions and making plots
    for subGPKey in parmObjFnKeys:
        # First calculate the statistics with outlets
        if not subGPKey == "Other":
            # Generate combined figures
            fnpPlotPng = "{}/ObsVsSim_UserBest_{}_{}_{}.png".format(fdOutputs, subGPKey, nlumpdist, nCalVal)
            if os.path.isfile(fnpPlotPng):
                os.remove(fnpPlotPng)
            obsSimThisKey = obsSimPairKeysSplit[subGPKey]
            genFigSingleOlt(fnpPlotPng, obsSimThisKey, obsSimPair,nCalVal, groupSubareaIdx)
            
            

##########################################################################
def genFigSingleOlt(fnpPlotPng, subGPKey, obsSimPair,nCalVal, groupSubareaIdx):
    """
    This function generate figure for one outlet.
    """
    fig = None

    xlabelfontsize = 8
    ylabelfontsize = 8
    tickfontsize = 8
    legendfontsize = 7
    scatterMarkersizeNg = 1
    scatterMarkersizeG = 2
    scatterMarkerLineWidth = 0.3

    spineLineWd = 0.2
    gridLineWidth = 0.3
    xtickMin = 0
    xtickMax = 500

    fig, axes = plt.subplots(1, 1,
                        figsize=(5, 5), 
                        dpi=300,
                        tight_layout=True)

    axes.plot(obsSimPair[subGPKey][0], label="Observation")
    axes.plot(obsSimPair[subGPKey][1], label="Simulation")

    if groupSubareaIdx == 0:
        legendTitle = "Lumped"
    elif groupSubareaIdx == 1:
        legendTitle = "Distributed"

    axes.legend(title = "{}".format(legendTitle), 
            loc="upper right",fontsize='x-large')
    oltNo = subGPKey.split("_")[0]
    axes.set_title("Station {} {}".format(oltNo, nCalVal))

    axes.set_xlabel("Number of Runs",
                fontsize=xlabelfontsize)
    axes.set_ylabel("Pertubation Value", 
                    fontsize=ylabelfontsize)

    # Control grids
    axes.grid(which='major',
            linewidth=gridLineWidth)
            
    # Control ticks
    axes.tick_params(
        which='major',
        direction = "in",
        labelsize=tickfontsize)

    fig.savefig(fnpPlotPng, bbox_inches="tight")

    
##########################################################################
def getRch2DFPlot(fnRch, iPrintForCio, totalRchNum):
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
def buildObsSimPairForPlot(obsDict, simOutRchDF, varIDObsHdrPair):

    """
    This function read the data from simDF and add corresponding columns
    into the obsDict.
    """

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
# This function is modified from the function developed by Florian Ulrich Jehn
# https://www.researchgate.net/publication/323800527_An_easy_interface_to_plot_flow_duration_curves_in_Python_with_and_without_uncertainty
def flow_duration_curve(timeSeries, iCalVal, axes, axesCtr, oltNo, nlumpDist, comparison = None):

    """
    Calculates and plots a flow duration curve from timeSeries. 
    
    All observations/simulations are ordered and the empirical probability is
    calculated. This is then plotted as a flow duration curve. 
    
    Additionally a comparison can be given to the function, which is plotted in
    the same ax.
    
    :param timeSeries: list of simulated and/or observed flow
    :param comparison: numpy array or pandas dataframe of discharge that should
    also be plotted in the same ax
    :param axis: int, axis along which x is iterated through
    :param ax: matplotlib subplot object, if not None, will plot in that 
    instance
    :param plot: bool, if False function will not show the plot, but simply
    return the ax object
    :param log: bool, if True plot on loglog axis
    :param percentiles: tuple of int, percentiles that should be used for 
    drawing a range flow duration curve
    :param fdc_kwargs: dict, matplotlib keywords for the normal fdc
    :param fdc_range_kwargs: dict, matplotlib keywords for the range fdc
    :param fdc_comparison_kwargs: dict, matplotlib keywords for the comparison 
    fdc
    
    return: subplot object with the flow duration curve in it
    """
       
    xlabelfontsize = 8
    ylabelfontsize = 8
    tickfontsize = 8
    legendfontsize = 7

    spineLineWd = 0.2
    gridLineWidth = 0.3
    xtickMin = 0
    xtickMax = 500
    
    axes[axesCtr] = plot_single_flow_duration_curve(axes[axesCtr], timeSeries, oltNo)
    
    # Add a comparison to the plot if is present
    if comparison is not None:
        axes[axesCtr] = plot_single_flow_duration_curve(axes[axesCtr], comparison, oltNo, iObs=True) 
        
    # Figure refine
    # Control legend
    axes[axesCtr].legend(title = "{} {} {}".format(nlumpDist, oltNo, iCalVal),
                fontsize = legendfontsize,
                title_fontsize = legendfontsize,
                framealpha = 1, 
                #loc = 'upper right',
                edgecolor = "white" 
                )
            
    # box = ax1.get_position()
    # # setposition(left, bottom, width, height)
    # ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Control label
    # set: a property batch setter
    # set_xlabel(xlabel, labelpad, **kwargs)
    axes[axesCtr].set_xlabel("Exceedance Probability",
                fontsize=xlabelfontsize)
    # Only create y label for the first column of figures
    if axesCtr % 3 == 0:
        axes[axesCtr].set_ylabel("Discharge (cubic meter/second)", 
                    fontsize=ylabelfontsize)
    
    # Control grids
    axes[axesCtr].grid(which='major',
              linewidth=gridLineWidth)
            
    # Control ticks
    axes[axesCtr].set_xlim(left=-5, right=105)
    
    axes[axesCtr].tick_params(
        which='major',
        direction = "in",
        labelsize=tickfontsize)

##########################################################################
# This function is modified from the function developed by Florian Ulrich Jehn
# https://www.researchgate.net/publication/323800527_An_easy_interface_to_plot_flow_duration_curves_in_Python_with_and_without_uncertainty
def flow_duration_curve_single(timeSeries, iCalVal, axes, oltNo, nlumpDist, comparison = None):

    """
    Calculates and plots a flow duration curve from timeSeries. 
    
    All observations/simulations are ordered and the empirical probability is
    calculated. This is then plotted as a flow duration curve. 
    
    Additionally a comparison can be given to the function, which is plotted in
    the same ax.
    
    :param timeSeries: list of simulated and/or observed flow
    :param comparison: numpy array or pandas dataframe of discharge that should
    also be plotted in the same ax
    :param axis: int, axis along which x is iterated through
    :param ax: matplotlib subplot object, if not None, will plot in that 
    instance
    :param plot: bool, if False function will not show the plot, but simply
    return the ax object
    :param log: bool, if True plot on loglog axis
    :param percentiles: tuple of int, percentiles that should be used for 
    drawing a range flow duration curve
    :param fdc_kwargs: dict, matplotlib keywords for the normal fdc
    :param fdc_range_kwargs: dict, matplotlib keywords for the range fdc
    :param fdc_comparison_kwargs: dict, matplotlib keywords for the comparison 
    fdc
    
    return: subplot object with the flow duration curve in it
    """
       
    xlabelfontsize = 8
    ylabelfontsize = 8
    tickfontsize = 8
    legendfontsize = 7

    spineLineWd = 0.2
    gridLineWidth = 0.3
    xtickMin = 0
    xtickMax = 500
    
    axes = plot_single_flow_duration_curve(axes, timeSeries, oltNo)
    
    # Add a comparison to the plot if is present
    if comparison is not None:
        axes = plot_single_flow_duration_curve(axes, comparison, oltNo, iObs=True) 
        
    # Figure refine
    # Control legend
    axes.legend(title = "{} {} {}".format(nlumpDist, oltNo, iCalVal),
                fontsize = legendfontsize,
                title_fontsize = legendfontsize,
                framealpha = 1, 
                #loc = 'upper right',
                edgecolor = "white" 
                )
            
    # box = ax1.get_position()
    # # setposition(left, bottom, width, height)
    # ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Control label
    # set: a property batch setter
    # set_xlabel(xlabel, labelpad, **kwargs)
    axes.set_xlabel("Exceedance Probability",
                fontsize=xlabelfontsize)
    # Only create y label for the first column of figures
    axes.set_ylabel("Discharge (cubic meter/second)", 
                    fontsize=ylabelfontsize)
    
    # Control grids
    axes.grid(which='major',
              linewidth=gridLineWidth)
            
    # Control ticks
    axes.set_xlim(left=-5, right=105)
    
    axes.tick_params(
        which='major',
        direction = "in",
        labelsize=tickfontsize)



def plot_single_flow_duration_curve(ax, timeseries, oltNo, iObs=False):
    """
    Plots a single fdc into an ax.
    
    :param ax: matplotlib subplot object
    :param timeseries: list like iterable

    return: subplot object with a flow duration curve drawn into it
    """
    scatterMarkersizeNg = 1.5
    scatterMarkersizeG = 1.5
    scatterMarkerLineWidth = 0.1

    # Get the probability
    exceedence = numpy.arange(1., len(timeseries) + 1) / len(timeseries)
    exceedence *= 100
    if not iObs:
        ax.scatter(exceedence, 
                sorted(timeseries, reverse=True), 
                label="Simulated",
                marker = "o",
                linewidth=scatterMarkerLineWidth, 
                s = scatterMarkersizeG,
                facecolor = "none",
                edgecolor="blue"
                )
    else:
        ax.scatter(exceedence, 
                sorted(timeseries, reverse=True),
                label="Observed",
                marker = "+",
                s = scatterMarkersizeNg,
                linewidth=scatterMarkerLineWidth,  
                c="red"
                )
    # Figure refine
    
    
        
    return ax


def plot_range_flow_duration_curve(ax, x, percentiles, kwargs):
    """
    Plots a single range fdc into an ax.
    
    :param ax: matplotlib subplot object
    :param x: dataframe of several timeseries
    :param decimal_places: defines how finely grained the range flow duration 
    curve is calculated and drawn. A low values makes it more finely grained.
    A value which is too low might create artefacts.
    :param kwargs: dict, keyword arguments for matplotlib
    
    return: subplot object with a range flow duration curve drawn into it
    """
    # Get the probabilites
    exceedence = numpy.arange(1.,len(numpy.array(x))+1) /len(numpy.array(x))
    exceedence *= 100

    # Sort the data
    sort = numpy.sort(x, axis=0)[::-1]
    
    # Get the percentiles
    low_percentile = numpy.percentile(sort, percentiles[0], axis=1)
    high_percentile = numpy.percentile(sort, percentiles[1], axis=1)
    
    # Plot it, check for empty kwargs
    if kwargs is not None:
        ax.fill_between(exceedence, low_percentile, high_percentile, **kwargs)
    else:
        ax.fill_between(exceedence, low_percentile, high_percentile)
    return ax
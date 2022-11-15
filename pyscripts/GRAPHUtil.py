# -*- coding: utf-8 -*-
"""
Created on Tue May 12 2020 

This class is designed to create watershed and find the contributing 
subareas.


@author: Qingyu.Feng
"""

##########################################################################
# Import modules #########################################################
##########################################################################
import pandas as pd
import os, sys, copy
import math

try:
    from osgeo import ogr, osr, gdal
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')


##########################################################################
# Define classes #########################################################
##########################################################################

# Added on Mar 22, 2021 to take care of the non group situation
##########################################################################
def noGroupSubForOutlet(outLetForGroups, allSubLst):

    """
    This function create the same structure as grouping, with all
    sub group list as empty list. Since when not grouping, all 
    subareas will be modified/updated in a same way.
    """
    subGroups = {}
    # Convert the outlet list from number list to string list
    # The allSubLst has all reach no in str.
    outLetForGroups = list(map(str, outLetForGroups))

    for olIdx in outLetForGroups:
        if not olIdx in allSubLst:
            print("Outlet {} does not exist in this watershed!!".format(olIdx))
            sys.exit()

        subGroups[olIdx] = []

    subGroups["Other"] = []

    return subGroups















##########################################################################
def dealWithOverlayInSubGroups(subNoGroups):
    """
    This function remove the overlays among subarea groups.
    The main principle is to keep the independency of subarea groups, 
    especially those as branches (tributary streams) when they are
    included in the groups for outlets on the main stream.
    """
    
    subNoGroups1 = copy.deepcopy(subNoGroups)
    del(subNoGroups1["Other"])
    subNoGroups2 = copy.deepcopy(subNoGroups1)

    subNoGroupUpdated = {}

    for subOlt1, subLst1 in subNoGroups1.items():
        # print("........processing subgroup for outlet: {}.......".format(subOlt1))
        # subLst1 is for editing
        # subLst2 is for reference.
        tempLst = copy.deepcopy(subLst1)

        for subOlt2, subLst2 in subNoGroups2.items():
            # Exclude comparison between the same outlet
            if subOlt1 != subOlt2:
                # print("comparing between: {} and {}".format(subOlt1, subOlt2))
                # First check whether sublist and sublst2 has overlays.
                if checkCommonInTwoLists(tempLst, subLst2):
                    # If there are common, there will be two cases:
                    # 1. subLst1 > subLst2: subLst1 contains subLst2
                    # 2. subLst1 < subLst2: subLst1 contained by subLst2
                    # The main principle is to maintain the independence
                    # of small lists, due to the facts that it was 
                    # for one individual subarea.
                    # Under situation 1: remove elements of subLst2 from 
                    # subLst1
                    # Under situation 2: keep subList and do not modify.
                    # print("Overlay found between {} and {}".format(subOlt1, subOlt2))
                    # print("len of subLst1 and subLst2: {} and {}".format(len(tempLst), len(subLst2)))
                    if len(tempLst) > len(subLst2):
                        # print("len before remove extra: {}".format(len(tempLst)))
                        tempLst = removeExtra(tempLst, subLst2)
                        # print("len after remove extra: {}".format(len(tempLst)))
                    
        subNoGroupUpdated[subOlt1] = tempLst

    subNoGroupUpdated["Other"] = subNoGroups["Other"]
    # Check whether the remove succeed:
    # print("..................Checking Overlay results..................")
    # subNoGroupsNoOL2 = copy.deepcopy(subNoGroupUpdated)
    # for subOlt11, subLst11 in subNoGroupUpdated.items():
    #     # subLst1 is for editing
    #     # subLst2 is for reference.
    #     for subOlt22, subLst22 in subNoGroupsNoOL2.items():
    #         # Exclude comparison between the same outlet
    #         if subOlt11 != subOlt22:
    #             # First check whether sublist and sublst2 has overlays.
    #             if checkCommonInTwoLists(subLst11, subLst22):
    #                 print("Overlay found between {} and {}".format(subOlt11, subOlt22))
                

    return subNoGroupUpdated               
              

##########################################################################
def removeExtra(list1, list2):
    """
    This function removes the elements in list1 that is also contained in
    list2.
    """
    newLst = []
    for ele1 in list1:
        if not ele1 in list2:
            newLst.append(ele1)

    return newLst




##########################################################################
def checkCommonInTwoLists(list1, list2):
    """
    This function checks whether there are are common elements 
    between list1 and list2.
    Input: two lists
    Output: Boolean:
    """
    one = set(list1)
    two = set(list2)
    if (one & two):
        return True
    else:
        return False









##########################################################################
def getGraph(field_names, subAttSWAT):

    """
    This function read in the attribute table of the shapefile and
    create a graph to represent the watershed.
    The graph will be a dictionary, with the reach no as key and upstream 
    reach no as values.
    For example, wsGraph = {reachNo: [upStrmNo1, upStrmNo2]}

    For the river generated by the SWAT model, the records are stored
    in the logic of downstream, not upstream.
    So, we need to find the reach no and all its up streams. 
    ['OBJECTID', 'ARCID', 'GRID_CODE', 'FROM_NODE', 'TO_NODE', 
    'Subbasin', 'SubbasinR', 'AreaC', 'Len2', 'Slo2', 'Wid2', 
    'Dep2', 'MinEl', 'MaxEl', 'Shape_Leng', 'HydroID', 'OutletID']
    """

    # Get all stream receiving streams, 
    rchAttDF = pd.DataFrame.from_dict(subAttSWAT, orient='index', columns=field_names)

    rcvRchNo = rchAttDF["GRID_CODE"].unique()
    
    wsGraph = {}

    for rchID in rcvRchNo:
        rchAttSubset = rchAttDF[rchAttDF["TO_NODE"] == rchID].dropna()["FROM_NODE"]
        if len(list(rchAttSubset)) > 0:
            wsGraph[rchID] = list(rchAttSubset)
        else:
            wsGraph[rchID] =[]

    return wsGraph, rcvRchNo


##########################################################################
def dfs_iterative(graph, start):
    stack, path, pathminus = [start], [], []

    # Stack as the starting point
    while stack:
        # Signed vertex: for path routing
        signedvertex = stack.pop()
        # remove sign for looping
        vertex = str(abs(int(signedvertex)))
        # Mark vertex as visited.
        if vertex in path:
            continue
        # If not visited, append it.
        path.append(vertex)
        pathminus.append(signedvertex)
        for nbid in range(len(graph[vertex])):
            if nbid > 0:
                neighbor = '-%s' %(graph[vertex][nbid])
            else:
                neighbor = graph[vertex][nbid]
                
            stack.append(neighbor)

    return path, pathminus


##########################################################################
def groupSubForOutlet(outLetForGroups, wsGraph, allSubLst):

    subGroups = {}
    sudProcessed = []

    # Convert the outlet list from number list to string list
    # The allSubLst has all reach no in str.
    outLetForGroups = list(map(str, outLetForGroups))

    for olIdx in outLetForGroups:
        if not olIdx in allSubLst:
            print("Outlet {} does not exist in this watershed!!".format(olIdx))
            sys.exit()

        tempGroup = []
        tempGroup, tempGroup2 = dfs_iterative(wsGraph, olIdx)
        subGroups[olIdx] = tempGroup
        sudProcessed = sudProcessed + tempGroup
    
    subNotProcessed = []
    for subId in allSubLst:
        if not subId in sudProcessed:
            subNotProcessed.append(subId)

    if len(subNotProcessed)>0:
        subGroups["Other"] = subNotProcessed

    return subGroups
    




##########################################################################
def readShapeAttributes(finshp):

    '''
    Read shapefile and return all values in the
    attritube table.
    '''
    # Check file existence:
    if not os.path.isfile(finshp):
        sys.exit('ERROR: Input file does not exist please check!')
    
    subStrAtt = dict()
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(finshp, 0)
    layer = dataSource.GetLayer()

    # Get the field name
    field_names = [field.name for field in layer.schema]

    # Get the value of each field for all layers
    for feature in layer:
        values_list = [str(feature.GetField(j)) for j in field_names]
        subStrAtt[str(values_list[0])] = values_list
    
    return field_names, subStrAtt
import numpy as np
import pandas as pd
#import dask.dataframe as dd
import collections
from hgcalHistHelpers import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import statsmodels.api as sm

#---------------------------------------------------------------------------------------------------
# Analyze SimClusters tree
def analyzeSimClusters(ntuple,tree,maxEvents,outDir,output,verbosityLevel):
    
    #Out variables
    out = collections.defaultdict(list)
    histDict = {}

    #---------------------------------------------------------------------------------------------------
    # start event loop
    for event in ntuple:
        currentevent = event.entry()
        if event.entry() >= maxEvents and maxEvents != -1 : break
        if (verbosityLevel>=1 and currentevent % 1 == 0): print( "\nCurrent event: ", currentevent)

        simClusters = event.simClusters()
        
        for simClusIndex, simClus in enumerate(simClusters):
            for DetId in simClus.hits():
                out["sClusHits"].append(DetId)
                out["EventId"].append(currentevent)
                out["EventIdFromFile"].append(event.event())
                out["sClusIndex"].append(simClusIndex)
                #Fill the per Object variables but per hit to produce the dataframe
                out["pdgId"].append(simClus.pdgId())
                out["id"].append(simClus.id())
                out["particleId"].append(simClus.particleId())
                out["charge"].append(simClus.charge())
                out["simcluster_p"].append(simClus.simcluster_p())
                out["energy"].append(simClus.energy())
                out["et"].append(simClus.et())
                out["mass"].append(simClus.mass())
                out["mt"].append(simClus.mt())
                out["pt"].append(simClus.pt())
                out["phi"].append(simClus.phi())
                out["theta"].append(simClus.theta())
                out["eta"].append(simClus.eta())
                out["rapidity"].append(simClus.rapidity())
                out["status"].append(simClus.status())
                out["longLived"].append(simClus.longLived())
                out["numberOfSimHits"].append(simClus.numberOfSimHits())
                out["simEnergy"].append(simClus.simEnergy())
            for rh in simClus.rechit_eta():
                out["rechit_eta"].append(rh)
            for rh in simClus.rechit_phi():
                out["rechit_phi"].append(rh)
            for rh in simClus.rechit_pt():
                out["rechit_pt"].append(rh)
            for rh in simClus.rechit_energy():
                out["rechit_energy"].append(rh)
            for rh in simClus.rechit_x():
                out["rechit_x"].append(rh)
            for rh in simClus.rechit_y():
                out["rechit_y"].append(rh)
            for rh in simClus.rechit_z():
                out["rechit_z"].append(rh)
            for rh in simClus.rechit_time():
                out["rechit_time"].append(rh)
            for rh in simClus.rechit_simclusterid():
                out["rechit_simclusterid"].append(rh)
            for DetId in simClus.matched_hits():
                out["sClusMatchedHits"].append(DetId)
            for thick in simClus.hits_thickness():
                out["sClusHitsThick"].append(thick)
            for det in simClus.hits_dets():
                out["sClusHitsDet"].append(det)        
            for frac in simClus.fractions():
                out["sClusHitsFractions"].append(frac)
            for layer in simClus.layers():
                out["sClusHitsLayers"].append(layer)
            for wafer_u in simClus.wafers_u():
                out["wafer_u"].append(wafer_u)
            for wafer_v in simClus.wafers_v():
                out["wafer_v"].append(wafer_v)
            for cell_u in simClus.cells_u():
                out["cell_u"].append(cell_u)
            for cell_v in simClus.cells_v():
                out["cell_v"].append(cell_v)
            for cell_type in simClus.cells_type():
                out["cell_type"].append(cell_type)
            for cell_zside in simClus.cells_zside():
                out["cell_zside"].append(cell_zside)
    
    #---------------------------------------------------------------------------------------------------
    #for key, value in out.items(): print(key, len(value))
    #Finished loop over events. Create the per hit dataframe
    df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in out.items() ]))
    df.fillna(-99999,inplace=True)
    #print(df[["sClusHitsThick"]].head())
    print(df.head())

    #histDict = histValue1D(out["sClusHitsDet"], histDict, tag = "LayersFromSimhits", title = "SimCluster Hits Layers",   axunit = "Layers from simhits", binsBoundariesX = [100, 0, 100], ayunit = "# simhits", verbosityLevel=verbosityLevel)

    #histPrintSaveAll(histDict, outDir, output, tree, verbosityLevel)

    return df

#Find weird SimHits in one side when shooting only in the other side
def findSimHitsInOtherSide(df_simCl):
    
    #df =  pd.read_csv(simCl_csv)
    ddf = df_simCl[ df_simCl['sClusHitsLayers'] < 50]
    print(ddf.head())
    #print(ddf.to_string(index=False))
    #print(len(df),len(ddf),float(len(ddf))/float(len(df)))

    return ddf
    
#def recHitCalibration(csv_list): 
def recHitCalibration(df): 
  
    histDict = {}
    #df = pd.concat(map(pd.read_csv, csv_list))
    #Add a column with the rechit_energy x fraction
    df["recHitEneXFraction"] = df.apply(lambda row: (row["rechit_energy"] * row["sClusHitsFractions"]) if (row["sClusMatchedHits"] != 9999 and row["rechit_energy"] > 0.) else 0.,axis=1)
    #We want the sum of the above column per Event per hit thickness and per detector 
    recHitEneXFractionSumPerThick = df.groupby(['EventId','sClusHitsThick','sClusHitsDet']).agg( recHitEneXFractionSum  = ('recHitEneXFraction','sum'))
    print(recHitEneXFractionSumPerThick)
    #The object above has a multiindex so let's save the quantities we want
    #The conventions are: Tracker = 1, Muon = 2,Ecal = 3,Hcal = 4,Calo = 5,Forward = 6,VeryForward = 7,HGCalEE = 8,HGCalHSi = 9,HGCalHSc = 10,HGCalTrigger = 11 
    rHxEnFrSum_CE_E_120um = np.ma.masked_equal(recHitEneXFractionSumPerThick.query("sClusHitsThick == 120 & sClusHitsDet == 8")[['recHitEneXFractionSum']].to_numpy().flatten(),0)
    rHxEnFrSum_CE_E_200um = np.ma.masked_equal(recHitEneXFractionSumPerThick.query("sClusHitsThick == 200 & sClusHitsDet == 8")[['recHitEneXFractionSum']].to_numpy().flatten(),0)
    rHxEnFrSum_CE_E_300um =np.ma.masked_equal( recHitEneXFractionSumPerThick.query("sClusHitsThick == 300 & sClusHitsDet == 8")[['recHitEneXFractionSum']].to_numpy().flatten(),0)
    rHxEnFrSum_CE_H_120um = np.ma.masked_equal(recHitEneXFractionSumPerThick.query("sClusHitsThick == 120 & sClusHitsDet == 9")[['recHitEneXFractionSum']].to_numpy().flatten(),0)
    rHxEnFrSum_CE_H_200um = np.ma.masked_equal(recHitEneXFractionSumPerThick.query("sClusHitsThick == 200 & sClusHitsDet == 9")[['recHitEneXFractionSum']].to_numpy().flatten(),0)
    rHxEnFrSum_CE_H_300um =np.ma.masked_equal( recHitEneXFractionSumPerThick.query("sClusHitsThick == 300 & sClusHitsDet == 9")[['recHitEneXFractionSum']].to_numpy().flatten(),0)
    rHxEnFrSum_CE_H_Scint = np.ma.masked_equal(recHitEneXFractionSumPerThick.query("sClusHitsThick > 400 & sClusHitsDet == 10")[['recHitEneXFractionSum']].to_numpy().flatten(),0)
    #print(rHxEnFrSum_CE_E_120um,rHxEnFrSum_CE_E_200um,rHxEnFrSum_CE_E_300um,rHxEnFrSum_CE_H_120um,rHxEnFrSum_CE_H_200um,rHxEnFrSum_CE_H_300um,rHxEnFrSum_CE_H_Scint)
    #Add the Egen. Yes, hardcoded for the moment. Observe that this is per hit, so we cannot possibly know beforehand
    #which array is of the greatest length to fill the Egen according to that. 
    y = np.full( max(len(rHxEnFrSum_CE_E_120um),len(rHxEnFrSum_CE_E_200um),len(rHxEnFrSum_CE_E_300um),len(rHxEnFrSum_CE_H_120um),len(rHxEnFrSum_CE_H_200um),len(rHxEnFrSum_CE_H_300um),len(rHxEnFrSum_CE_H_Scint)) , 60.)
 
    #Will create the dataframe this way since the length of the columns is different
    df_reg = pd.DataFrame({
        'rHxEnFrSum_CE_E_120um': pd.Series(rHxEnFrSum_CE_E_120um), 
        'rHxEnFrSum_CE_E_200um': pd.Series(rHxEnFrSum_CE_E_200um), 
        'rHxEnFrSum_CE_E_300um': pd.Series(rHxEnFrSum_CE_E_300um), 
        'rHxEnFrSum_CE_H_120um': pd.Series(rHxEnFrSum_CE_H_120um), 
        'rHxEnFrSum_CE_H_200um': pd.Series(rHxEnFrSum_CE_H_200um), 
        'rHxEnFrSum_CE_H_300um': pd.Series(rHxEnFrSum_CE_H_300um), 
        #'rHxEnFrSum_CE_H_Scint': pd.Series(rHxEnFrSum_CE_H_Scint), 
        'theEgen': y})
    df_reg.fillna(0,inplace=True)
    
    print(df_reg.head())
        
    #X = df_reg[['rHxEnFrSum_CE_E_120um','rHxEnFrSum_CE_E_200um','rHxEnFrSum_CE_E_300um','rHxEnFrSum_CE_H_120um','rHxEnFrSum_CE_H_200um','rHxEnFrSum_CE_H_300um','rHxEnFrSum_CE_H_Scint','theEgen']]
    X = df_reg[['rHxEnFrSum_CE_E_120um','rHxEnFrSum_CE_E_200um','rHxEnFrSum_CE_E_300um','rHxEnFrSum_CE_H_120um','rHxEnFrSum_CE_H_200um','rHxEnFrSum_CE_H_300um','theEgen']]
    y = df_reg[['theEgen']]
 
    ols = sm.OLS(y,X)
    results = ols.fit()
    print(results.summary())
    
    print(type(results.params))
    print("Multiple Regression (not inverse): ", results.params.to_numpy()) 
    print("Multiple Regression (inverse): ", np.reciprocal(results.params.to_numpy()))
    print('Linear regression results')
    print('Accuracy (R^2):',results.rsquared)
    
    np.set_printoptions(precision=2)
    print("Multiple Regression (not inverse): ", results.params.to_numpy())
    print("Multiple Regression (inverse): ", np.reciprocal(results.params.to_numpy()))
    
    return df_reg

    #histDict = histValue1D(out["sClusHitsDet"], histDict, tag = "LayersFromSimhits", title = "SimCluster Hits Layers",   axunit = "Layers from simhits", binsBoundariesX = [100, 0, 100], ayunit = "# simhits", verbosityLevel=verbosityLevel)

    #histPrintSaveAll(histDict, outDir, output, tree, verbosityLevel)




#---------------------------------------------------------------------------------------------------
# Analyze LayerClusters from Tracksters 
def analyzeLayerClusters(layerClusters,currentevent,currenteventFromFile):

    #Out variables
    outLC = collections.defaultdict(list)

    for lc in layerClusters: 
        outLC["layerCluster_EventId"].append(currentevent)
        outLC["layerCluster_EventIdFromFile"].append(currenteventFromFile)
        outLC["layerCluster_id"].append(lc.id())
        outLC["layerCluster_trackster_id"].append(lc.trackster_id())
        outLC["layerCluster_eta"].append(lc.eta())
        outLC["layerCluster_phi"].append(lc.phi())
        outLC["layerCluster_pt"].append(lc.pt())
        outLC["layerCluster_energy"].append(lc.energy())
        outLC["layerCluster_x"].append(lc.x())
        outLC["layerCluster_y"].append(lc.y())
        outLC["layerCluster_z"].append(lc.z())
        outLC["layerCluster_layer"].append(lc.layer())
        #outLC["layerCluster_nhitCore"].append(lc.nhitCore)
        #outLC["layerCluster_nhitAll"].append(lc.nhitAll)

    dfl = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in outLC.items() ]))
    dfl.fillna(-99999,inplace=True)

    return dfl

#---------------------------------------------------------------------------------------------------
# Analyze Tracksters tree
def analyzeTracksters(ntuple,tree,maxEvents,outDir,output,verbosityLevel):
    
    #Out variables
    out = collections.defaultdict(list)
    histDict = {}

    dfl = pd.DataFrame()
    #---------------------------------------------------------------------------------------------------
    # start event loop
    for event in ntuple:
        if event.entry() >= maxEvents and maxEvents != -1 : break
        if (verbosityLevel>=1): print( "\nCurrent event: ", event.entry())
        currentevent = event.entry()
        currenteventFromFile = event.event()

        recHits  = event.recHits()
        #These LCs are the ones associated with tracksters. 
        layerClusters = event.layerClusters()
        tracksters = event.tracksters()

        #Create a dataframe with the LayerClusters of the Tracksters
        dflcur = analyzeLayerClusters(layerClusters, currentevent, currenteventFromFile)
        if currentevent == 0: dfl = dflcur
        else: dfl.append(dflcur)
        #dfl.to_csv('TrackstersAndLCs.csv', index=False)
        LC_idtoX = pd.Series(dflcur.layerCluster_x.values,index=dflcur.layerCluster_id).to_dict()
        LC_idtoY = pd.Series(dflcur.layerCluster_y.values,index=dflcur.layerCluster_id).to_dict()
        LC_idtoZ = pd.Series(dflcur.layerCluster_z.values,index=dflcur.layerCluster_id).to_dict()

        #print(len(LC_idtoX), len(LC_idtoY), len(LC_idtoZ))

        for trstIndex, trst in enumerate(tracksters):
            #print(trstIndex)
            for edge in trst.numedges():
                out["EventId"].append(currentevent)
                out["EventIdFromFile"].append(currenteventFromFile)
                out["trstIndex"].append(trstIndex)
                #Per Trackster in per edge
                out["trackster_x"].append(trst.x())
                out["trackster_y"].append(trst.y())
                out["trackster_z"].append(trst.z())
                out["trackster_eta"].append(trst.eta())
                out["trackster_phi"].append(trst.phi())
                out["trackster_firstlayer"].append(trst.firstlayer())
                out["trackster_lastlayer"].append(trst.lastlayer())
                out["trackster_layersnum"].append(trst.layersnum())
                out["trackster_id"].append(trst.id())
                out["trackster_time"].append(trst.time())
                out["trackster_timeError"].append(trst.timeError())
                out["trackster_regr_energy"].append(trst.regr_energy())
                out["trackster_raw_energy"].append(trst.raw_energy())
                out["trackster_raw_em_energy"].append(trst.raw_em_energy())
                out["trackster_raw_pt"].append(trst.raw_pt())
                out["trackster_raw_em_pt"].append(trst.raw_em_pt())
                out["trackster_id_prob_1"].append(trst.id_prob_1())
                out["trackster_id_prob_2"].append(trst.id_prob_2())
                out["trackster_id_prob_3"].append(trst.id_prob_3())
                out["trackster_id_prob_4"].append(trst.id_prob_4())
                out["trackster_id_prob_5"].append(trst.id_prob_5())
                out["trackster_id_prob_6"].append(trst.id_prob_6())
                out["trackster_id_prob_7"].append(trst.id_prob_7())
                out["trackster_id_prob_8"].append(trst.id_prob_8())
                #Per edge
                out["trackster_edge"].append(edge)
            for layer_in in trst.edge_layerin():
                out["trackster_edge_layerin"].append(layer_in)
            for layer_in_id in trst.edge_layerin_id():
                out["trackster_edge_layerin_id"].append(layer_in_id)
                out["trackster_edge_layerin_x"].append(LC_idtoX[layer_in_id])
                out["trackster_edge_layerin_y"].append(LC_idtoY[layer_in_id])
                out["trackster_edge_layerin_z"].append(LC_idtoZ[layer_in_id])
                out["trackster_edge_layerin_R"].append( np.sqrt( LC_idtoX[layer_in_id]**2 + LC_idtoY[layer_in_id]**2 ) )
                #print(layer_in_id, LC_idtoX[layer_in_id])
            for layer_out in trst.edge_layerout():
                out["trackster_edge_layerout"].append(layer_out)
            for layer_out_id in trst.edge_layerout_id():
                out["trackster_edge_layerout_id"].append(layer_out_id)
                out["trackster_edge_layerout_x"].append(LC_idtoX[layer_out_id])
                out["trackster_edge_layerout_y"].append(LC_idtoY[layer_out_id])
                out["trackster_edge_layerout_z"].append(LC_idtoZ[layer_out_id])
                out["trackster_edge_layerout_R"].append( np.sqrt( LC_idtoX[layer_out_id]**2 + LC_idtoY[layer_out_id]**2 ) )     
            for delta_energy in trst.delta_energy():
                out["trackster_delta_energy"].append(delta_energy)
            for delta_energy_relative in trst.delta_energy_relative():
                out["trackster_delta_energy_relative"].append(delta_energy_relative)
            for delta_layer in trst.delta_layer():
                out["trackster_delta_layer"].append(delta_layer)
            for angle_alpha in trst.angle_alpha():
                out["trackster_angle_alpha"].append(angle_alpha)
            for angle_alpha_alternative in trst.angle_alpha_alternative():
                out["trackster_angle_alpha_alternative"].append(angle_alpha_alternative)
            #beta angle has a different number of entries and will be dealt with differently    
            #for angle_beta in trst.angle_beta():
                #out["trackster_angle_beta"].append(angle_beta)

    #print(layer_in)
    #print(dfl[ (dfl['layerCluster_id'] == layer_in) & (dfl['layerCluster_trackster_id'] == trst.id())])
    #ddfl = dfl[ dfl['layerCluster_id'] == layer_in]
    #print(ddfl[['layerCluster_x','layerCluster_y']].values[[0]])
    #out["trackster_edge_layerin_x"].append()
    #out["trackster_edge_layerin_y"].append()
    #print(dfl)
    #---------------------------------------------------------------------------------------------------
    #Finished loop over events. Create the per trackster dataframe
    df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in out.items() ]))
    #df.fillna(-99999,inplace=True)

    #df.to_csv('TrackstersWithEdgesDataFrame.csv', index=False)
    chosenId = 0

    ddf = df[ (df['EventId'] == chosenId)]
    ddfl = dfl[ (dfl['layerCluster_EventId'] == chosenId)]
 
    #Simple 3D scatter plot: Edges_X_Y_Layer_plt.png
    drawEdges3d_plt(ddf, chosenId)
    #Simple 2D : Edges_R_vs_Layer_plt.png
    drawEdges2d_plt(ddf,chosenId)
    #2D plot from DiGraph R vs Layer:  Edges_R_vs_Layer_DiGraph.png
    drawDirectedGraph(ddf,chosenId)
    #Use plotly and save 3D plot as html page 
    drawDirectedGraph3D_plotly(ddf,dfl,chosenId)


    histDict = histValue1D(out["sClusHitsDet"], histDict, tag = "LayersFromSimhits", title = "SimCluster Hits Layers",   axunit = "Layers from simhits", binsBoundariesX = [100, 0, 100], ayunit = "# simhits", verbosityLevel=verbosityLevel)

    #histPrintSaveAll(histDict, outDir, output, tree, verbosityLevel)

    return df
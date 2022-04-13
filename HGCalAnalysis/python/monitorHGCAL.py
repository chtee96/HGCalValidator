# Monitor HGCAL objects

# Example command of a basic setup: 
# python3 Validation/HGCalValidation/python/monitorHGCAL.py --input hgc.root --object SimClusters --maxEvents -1 --verbosityLevel 2 --output output_test --outDir output

#---------------------------------------------------------------------------------------------------
#Necessary imports
from __future__ import print_function
import ROOT
import math
import os,sys
import numpy as np
import pandas as pd
from NtupleDataFormat import HGCalNtuple
from hgcalHelpers import *
import multiprocessing as mp

ROOT.gROOT.SetBatch(True) 

#---------------------------------------------------------------------------------------------------
#Options from command line
from optparse import OptionParser
parser = OptionParser(description="HGCAL Objects Analysis")
parser.add_option("--input",help="The input hgcal ntuple file")
parser.add_option("--object",help="The object we want to monitor")
parser.add_option("--maxEvents",help="The maximum number of events you want to process. -1 for all events. ", type=int)
parser.add_option("--verbosityLevel",help=" 0 - only basic info (default); 1 - additional info; 2 - detailed info printed, histograms produced",type=int)
parser.add_option("--output",default="output.root",help="Name of output root file to keep the histos (default: %default)")
parser.add_option("--genEnergy",help="Generated energy of the sample",type=int)
parser.add_option("--outDir",default="output",help="Output directory with all plots and files (default: %default)")
parser.add_option("--ecut",help="The rechit energy cut relative to noise",type=float)

(options,args)=parser.parse_args()

def main():

    #---------------------------------------------------------------------------------------------------
    # Preparation
    outDir = options.outDir
    if not os.path.exists(outDir): os.makedirs(outDir)

    #The ntuple will be created only when reading the objects
    ntuple = None
    if (options.object == "TrackstersAssociators"):
        ntuple = HGCalNtuple(options.input, "HGCalAnalysis/ticlTrackstersMerge")
    elif options.object != "RecHitCalibration":
        ntuple = HGCalNtuple(options.input, "HGCalAnalysis/%s" %(options.object) )

    #The tree below where we will save all our results apart from 
    #the individual pngs we may need to also plot. 
    tree = ROOT.TTree("ttree_%s_e%s" %(options.object, options.genEnergy), "results")
    
    #Check to see whether the output file is there to avoid overwriting. It 
    #assumes that you are reading from eos, since these are huge files
    if os.path.exists("%s/%s_%s.csv" %(options.outDir,options.input.replace(".root","").split("/")[-1], options.object)): 
        print("The output file exists. Backup if needed and rerun")
        exit()

    if (options.object == "SimClusters"):
        df = analyzeSimClusters(ntuple,tree,options.maxEvents,outDir,options.output,options.verbosityLevel)
        #We will save the df in a csv file, since running in multiple samples will take some time. Then, 
        #we can load that csv in a dataframe and do analysis we want
        #df.to_csv('%s/%s_%s.csv' %(options.outDir,options.input.replace(".root","").split("/")[-1],options.object), index=False)

        #Check for weird SimHits in one side when shooting only in the other side. 
        #simHits = findSimHitsInOtherSide(df)
        #simHits.to_csv('%s/%s_%s_checkOtherSide.csv' %(options.outDir,options.input.replace(".root","").split("/")[-1],options.object), index=False)
        #When simcluster id will be here use the one below
        #ddf = df[ df["rechit_simclusterid"] >= 0 ]
        recHitCalibration(df,tree,options.maxEvents,outDir,options.output,options.genEnergy,options.ecut,options.verbosityLevel)
        #recHitsCalib = recHitCalibration(df[ df["sClusMatchedHits"] != 9999])
       
    if (options.object == "LayerClusters"):
        df = analyzeLayerClusters(ntuple,tree,options.maxEvents,outDir,options.output,options.genEnergy,options.verbosityLevel)

    if (options.object in ["ticlTrackstersEM","ticlTrackstersTrkEM","ticlTrackstersTrk","ticlTrackstersHAD","ticlTrackstersMerge","ticlSimTracksters"] ):
        #df = analyzeTracksters(ntuple,tree,options.maxEvents,outDir,options.output,options.verbosityLevel)
        df = analyzeSimClusters(ntuple,tree,options.maxEvents,outDir,options.output,options.verbosityLevel)
        recHitCalibration(df,tree,options.maxEvents,outDir,options.output,options.genEnergy,options.ecut,options.verbosityLevel)        
        #print(df.head())

    if (options.object == "TrackstersAssociators"):
        df = analyzeTrackstersAssociators(ntuple,tree,options.maxEvents,outDir,options.output,options.verbosityLevel)
        
    if (options.object == "RecHitCalibration"):
        
        csv_list = ["%s/%s" % (options.outDir,cs) for cs in options.input.split(',')]
        print(csv_list)
        recHitsCalib = recHitCalibration(csv_list)
        recHitsCalib.to_csv('%s/%s_RecHitCalib.csv' %(options.outDir,options.object), index=False)
 



if __name__ == '__main__':
    main()

###!/usr/bin/python
import sys
import os
import commands
import optparse
import glob
#import ROOT

# Production based on guidelines by Chris

### parsing input options
def parseOptions():

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-q', '--queue',  dest='QUEUE',  type='string',  default='tomorrow', help='queue to be used with HTCondor, default is tomorrow')
    parser.add_option('-t', '--tag',    dest='TAG',    type='string',  default='', help='tag to be appended to the resulting output dir, default is an empty string')
    parser.add_option('-n', '--nevts',  dest='NEVTS',  type=int,       default=100,  help='total number of events, applicable to runs with GEN stage, default is 100')
    parser.add_option('-e', '--evtsperjob', dest='EVTSPERJOB', type=int, default=-1,   help='number of events per job, if set to -1 it will set to a recommended value (GSD: 4events/1nh, RECO:8events/1nh), default is -1')
    parser.add_option('-d', '--datTier', dest='DTIER',  type='string', default='GSD', help='data tier to run: "GSD" (GEN-SIM-DIGI) or "RECO", default is "GSD"')
    parser.add_option('-p', '--partID', dest='PARTID', type='string',     default='', help='string of particle PDG IDs separated by comma, if empty string - run on all supported (11,12,13,14,15,16,22,111,211,130) and corresponding negative values, default is empty string (all)')
    parser.add_option('', '--nPart',  dest='NPART',  type=int,   default=1,      help='number of times particles of type(s) PARTID will be generated per event, default is 1')
    parser.add_option('', '--nRandomPart',  dest='NRANDOMPART',  type=int,   default=1,      help='This is used together with randomShoot to shoot randomly [1, NParticles +1 ] particles, default is 1')
    parser.add_option('', '--gunMode',   dest='gunMode',   type='string', default='default',    help='default or pythia8 or closeby')
    parser.add_option('', '--gunType',   dest='gunType',   type='string', default='Pt',    help='Pt or E gun')
    parser.add_option('', '--eosArea', dest='eosArea', type='string', default='', help='path to the eos area where the output jobs will be staged out')
    parser.add_option('', '--date', dest='DATE',  type='string', default='nodate', help='If running RECO this variable is the last part of the folder name from the GSD input. If not specified the code will complain.')
    parser.add_option('', '--html-validation-name', dest='HTMLVALNAME', type='string', default='', help='Could be either be hgcalLayerClusters or hgcalMultiClusters')
    parser.add_option('', '--zMin',  dest='zMin',  type=float, default=321.6,  help='min. z value start of EE at V10')
    parser.add_option('', '--zMax',  dest='zMax',  type=float, default=650.0,    help='max. z value')
    parser.add_option('', '--rMin',  dest='rMin',  type=float, default=0.0,  help='min. r value')
    parser.add_option('', '--rMax',  dest='rMax',  type=float, default=300.0,    help='max. r value')
    parser.add_option('', '--pointing',  action='store_false',  dest='pointing',  default=True,    help='pointing to (0,0,0) in case of closeby gun')
    parser.add_option('', '--overlapping',  action='store_true',  dest='overlapping',  default=False,    help='particles will be generated in window [phiMin,phiMax], [rMin,rMax] (true) or with a DeltaPhi=Delta/R (default false) in case of closeby gun')
    parser.add_option('', '--Delta',  dest='Delta',  type=float, default=0.25, help=' arc-distance between two consecutive vertices over the circle of radius R')
    parser.add_option('-y', '--dry-run', action='store_true', dest='DRYRUN', default=False, help='perform a dry run (no jobs are lauched).')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # sanity check for data tiers
    dataTiers = ['GSD', 'RECO','NTUP','AFTERHARVESTING','PLOTTING']
    if opt.DTIER not in dataTiers:
        parser.error('Data tier ' + opt.DTIER + ' is not supported. Exiting...')
        sys.exit()

### processing the external os commands
def processCmd(cmd, quite = 0):
    print(cmd)
    status, output = commands.getstatusoutput(cmd)
    if (status !=0 and not quite):
        print('Error in processing command:\n   ['+cmd+']')
        print('Output:\n   ['+output+'] \n')
    return output

parseOptions()

#If the main source CMSSW area has root files it will cause a problem. Start with a clean one. 

Pt = []
Eta = []
En = []

SampleName = ''

#PHOTONS

'''
if opt.PARTID == '22':
    Pt = [25,35,50,100]
    #Pt = [100]
    Eta = [1.6,2.8]
    SampleName = 'Photons'
'''
#CLOSE BY PHOTONS
#if opt.PARTID == '22,22': 
if opt.PARTID == '22': 
    #En = [20,60,100]
    En = [60]
    Eta = [1.4,4.0]
    #Deltas = [2.5,5.0,7.5,15.0]
    Deltas = [2.5]
    SampleName = 'Close_By_Photons'

#CHARGED PIONS
if opt.PARTID == '211':
    Pt = [5,7,10,15,25,35,50,100]
    Eta = [1.6,2.8]
    SampleName = 'Charged Pions'

if (opt.DTIER == 'GSD'):
    if (opt.gunMode == 'default' or opt.gunMode == 'pythia8' ):
        for thePt in Pt:
            cmd = 'python SubmitHGCalPGun.py --datTier GSD --nevts ' + str(opt.NEVTS) +' --evtsperjob '+ str(opt.EVTSPERJOB) +' --queue '+ opt.QUEUE + ' --partID ' + str(opt.PARTID) + ' --nPart ' + str(opt.NPART) + ' --thresholdMin ' + str(thePt) + ' --thresholdMax ' + str(thePt) + ' --etaMin ' + str(Eta[0]) + ' --etaMax ' + str(Eta[1]) + ' --gunType Pt --eosArea ' + opt.eosArea + ' --tag ${USER}_PDGId%s_nPart%s_Pt%s_eta%sto%s' % (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'))
            if(opt.DRYRUN):
                print('Dry-run: ['+cmd+']')
            else:
                output = processCmd(cmd)
    if (opt.gunMode == 'closeby'):
        for theEn in En:
            for theDelta in Deltas:
                cmd = 'python SubmitHGCalPGun.py --datTier GSD --nevts ' + str(opt.NEVTS) +' --evtsperjob '+ str(opt.EVTSPERJOB) +' --queue '+ opt.QUEUE + ' --partID ' + str(opt.PARTID) + ' --nPart ' + str(opt.NPART) + ' --nRandomPart ' + str(opt.NRANDOMPART)+  ' --thresholdMin ' + str(theEn) + ' --thresholdMax ' + str(theEn) + ' --etaMin ' + str(Eta[0]) + ' --etaMax ' + str(Eta[1]) + ' --zMin ' + str(opt.zMin) + ' --zMax ' + str(opt.zMax) + ' --rMin ' + str(opt.rMin) + ' --rMax ' + str(opt.rMax) + ' --gunType E --gunMode closeby --eosArea ' + opt.eosArea + ' --tag ${USER}_PDGId%s_nPart%s_E%s_eta%sto%s_%s_Delta_%s --Delta %s' % (opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG, str(theDelta).replace('.','p'), str(theDelta))
                if opt.pointing == False: 
                    cmd = cmd + ' --pointing '
                if opt.overlapping == True: 
                    cmd = cmd + ' --overlapping '
                if(opt.DRYRUN):
                    print('Dry-run: ['+cmd+']')
                else:
                    output = processCmd(cmd)

if (opt.DTIER == 'RECO'):
    if (opt.gunMode == 'default' or opt.gunMode == 'pythia8' ):
        for thePt in Pt:
            cmd = 'python SubmitHGCalPGun.py --datTier '+ opt.DTIER + ' --evtsperjob '+ str(opt.EVTSPERJOB) + ' --queue ' + opt.QUEUE +' --inDir FlatRandomPtGunProducer_${USER}_PDGId%s_nPart%s_Pt%s_eta%sto%s_%s_Delta_%s '% (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.DATE) + ' --eosArea ' + opt.eosArea + ' --tag ${USER}_PDGId%s_nPart%s_Pt%s_eta%sto%s ' % (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p')) + ' --keepDQMfile ' #+ ' --keepDQMfile '
            if(opt.DRYRUN):
                print('Dry-run: ['+cmd+']')
            else:
                output = processCmd(cmd)
    if (opt.gunMode == 'closeby'):
        for theEn in En:
            for theDelta in Deltas:
                cmd = 'python SubmitHGCalPGun.py --datTier '+ opt.DTIER + ' --evtsperjob '+ str(opt.EVTSPERJOB) + ' --queue ' + opt.QUEUE +' --inDir CloseByParticleGunProducer_${USER}_PDGId%s_nPart%s_E%s_eta%sto%s_%s_Delta_%s_%s '% (opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG, str(theDelta).replace('.','p'),opt.DATE) + ' --eosArea ' + opt.eosArea + ' --tag ${USER}_PDGId%s_nPart%s_E%s_eta%sto%s_%s_Delta_%s_%s ' % (opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG, str(theDelta).replace('.','p'),opt.DATE) + ' --keepDQMfile ' #+ ' --keepDQMfile '
                if(opt.DRYRUN):
                    print('Dry-run: ['+cmd+']')
                else:
                    output = processCmd(cmd)



if (opt.DTIER == 'NTUP'):
    if (opt.gunMode == 'default' or opt.gunMode == 'pythia8' ):
        for thePt in Pt:
            cmd = 'python SubmitHGCalPGun.py --datTier '+ opt.DTIER + ' --evtsperjob '+ str(opt.EVTSPERJOB) + ' --queue ' + opt.QUEUE +' --inDir FlatRandomPtGunProducer_${USER}_PDGId%s_nPart%s_Pt%s_eta%sto%s_%s_Delta_%s '% (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.DATE) + ' --eosArea ' + opt.eosArea + ' --tag ${USER}_PDGId%s_nPart%s_Pt%s_eta%sto%s ' % (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p')) + ' --keepDQMfile ' #+ ' --keepDQMfile '
            if(opt.DRYRUN):
                print('Dry-run: ['+cmd+']')
            else:
                output = processCmd(cmd)
    if (opt.gunMode == 'closeby'):
        for theEn in En:
            for theDelta in Deltas:
                cmd = 'python SubmitHGCalPGun.py --datTier '+ opt.DTIER + ' --evtsperjob '+ str(opt.EVTSPERJOB) + ' --queue ' + opt.QUEUE +' --inDir CloseByParticleGunProducer_${USER}_PDGId%s_nPart%s_E%s_eta%sto%s_%s_Delta_%s_%s_${USER}_PDGId%s_nPart%s_E%s_eta%sto%s_%s_Delta_%s_%s '% (opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG, str(theDelta).replace('.','p'),opt.DATE,opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG, str(theDelta).replace('.','p'),opt.DATE) + ' --eosArea ' + opt.eosArea + ' --tag ${USER}_PDGId%s_nPart%s_E%s_eta%sto%s_%s_Delta_%s_%s ' % (opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG, str(theDelta).replace('.','p'),opt.DATE) + ' --keepDQMfile ' #+ ' --keepDQMfile '
                if(opt.DRYRUN):
                    print('Dry-run: ['+cmd+']')
                else:
                    output = processCmd(cmd)



if (opt.DTIER == 'AFTERHARVESTING'):

    for thePt in Pt:
        #This is where we want to put the skimmed dqm files
        outDir = opt.eosArea + '/SKIMMED_FlatRandomPtGunProducer_apsallid_PDGId%s_nPart%s_Pt%s_eta%sto%s_%s/DQM'% (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.DATE)
        mergedoutDir = opt.eosArea + '/MERGED_FlatRandomPtGunProducer_apsallid_PDGId%s_nPart%s_Pt%s_eta%sto%s_%s/DQM'% (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.DATE)
        if (not os.path.isdir(outDir)) or (not os.path.isdir(mergedoutDir)) :
            processCmd('mkdir -p '+outDir)
            processCmd('mkdir -p '+mergedoutDir)
        #else:
            #print 'Directory '+outDir+' already exists. Exiting...'
            #sys.exit()
    

        pathtodqmfiles = opt.eosArea + '/FlatRandomPtGunProducer_apsallid_PDGId%s_nPart%s_Pt%s_eta%sto%s_%s/DQM'% (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.DATE)

        #Gather all the DQM files from the RECO step to do the harvesting
        fullinputList = [f for f in processCmd( ' ls ' +  pathtodqmfiles  ).split('\n') ]
        #print fullinputList

        #Check for zombie files
        clearinputList = []
        for infi in fullinputList: 
            fin = ROOT.TFile.Open('root://eoscms.cern.ch/' + pathtodqmfiles + '/' + infi)
            #print infi
            if fin: 
                #print infi
                clearinputList.append(infi)

        print("="*40)    
        #print clearinputList
        if len(clearinputList) != 1: 
            print("Starting Skimming")
            for infi in clearinputList: 
                processCmd( 'root.exe -b -q validationplots.C\(\\"root://eoscms.cern.ch/' + pathtodqmfiles + '/' + infi  +  '\\"\)'  ) 
                processCmd( 'cp newfi.root ' + outDir + '/'+  infi  )
                #sys.exit()

            #Ready to merge
            processCmd( 'hadd -f test.root ' + outDir + '/*.root'   )
            processCmd( ' cp test.root ' + mergedoutDir + '/' + 'closeby_PDGid%s_x10_Pt%s.0To%s.0_DQM.root' % (opt.PARTID,thePt,thePt) )

#For the plotting step we must have only 1 file per pt or E
if (opt.DTIER == 'PLOTTING'):

    #We also need a variable to easily create the html file
    fragments = []
    if (opt.gunMode == 'default' or opt.gunMode == 'pythia8' ):
        for thePt in Pt:
            #Let's first take the input DQM files in plain ROOT format
            pathtodqmfiles = opt.eosArea + '/FlatRandomPtGunProducer_apsallid_PDGId%s_nPart%s_Pt%s_eta%sto%s_%s/DQM'% (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.DATE)
            
            fullinputList = [f for f in processCmd( ' ls ' +  pathtodqmfiles  ).split('\n') ]
            #print fullinputList
            
            #Check for zombie files
            clearinputList = []
            for infi in fullinputList: 
                fin = ROOT.TFile.Open('root://eoscms.cern.ch/' + pathtodqmfiles + '/' + infi)
                #print infi
                if fin: 
                    #print infi
                    #There is no need for the 'root://eoscms.cern.ch/' prefix
                    clearinputList.append(pathtodqmfiles + '/' + infi)

                    print("="*40)  
                    #print clearinputList

            #Ready to create plots for each file
            #We should be at the reco_prodtools directory
            for infi in clearinputList: 
                cmd = 'python ../Validation/HGCalValidation/scripts/makeHGCalValidationPlots.py ' + infi + ' --outputDir HGCValid_%s_Plots --no-ratio --png --separate --html-sample "' %(opt.HTMLVALNAME) + SampleName + ' Pt %s GeV Eta %s to %s" ' % (thePt,str(Eta[0]),str(Eta[1])) + ' --html-validation-name %s --subdirprefix ' %(opt.HTMLVALNAME) + ' plots_PDGId%s_nPart%s_Pt%sEta%sto%s' % (opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'))

                if(opt.DRYRUN):
                    print('Dry-run: ['+cmd+']')
                else:
                    output = processCmd(cmd)
                    #processCmd( ' mv HGCalValidationPlots/index.html HGCalValidationPlots/index_Pt%sEta%sto%s.html'% (thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'))  )
                    processCmd('awk \'NR>=6&&NR<=44\' HGCValid_%s_Plots/index.html > HGCValid_%s_Plots/index_PDGId%s_nPart%s_Pt%sEta%sto%s.html '% (opt.HTMLVALNAME,opt.HTMLVALNAME,opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'))  )
         
                #cmd = cmd + 'sed -i \'s/Tracking/HGCal/g\' HGCalValidationPlots/index_Pt%sEta%sto%s.html'% (thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p')) + '\n '
                #ALWAYS CHECK THE FRAGMENT YOU WANT
                fragments.append( 'HGCValid_%s_Plots/index_PDGId%s_nPart%s_Pt%sEta%sto%s.html'% (opt.HTMLVALNAME,opt.PARTID,opt.NPART,thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p')) )
            
    if (opt.gunMode == 'closeby'):
       #We will create a per tag structure. It seems more readable. 
       processCmd('mkdir -p HGCValid_%s_Plots/%s' %(opt.HTMLVALNAME,opt.TAG) )
       for theEn in En:
           clearinputList = []
           for theDelta in Deltas:
               #Let's first take the input DQM files in plain ROOT format
               pathtodqmfiles = opt.eosArea + '/CloseByParticleGunProducer_apsallid_PDGId%s_nPart%s_E%s_eta%sto%s_%s_Delta_%s_%s/DQM'% (opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG, str(theDelta).replace('.','p'),opt.DATE)
               
               fullinputList = [f for f in processCmd( ' ls ' +  pathtodqmfiles  ).split('\n') ]
               #print fullinputList
               
               #Check for zombie files
               for infi in fullinputList: 
                   fin = ROOT.TFile.Open('root://eoscms.cern.ch/' + pathtodqmfiles + '/' + infi)
                   #print infi
                   if fin: 
                       #print infi
                       #There is no need for the 'root://eoscms.cern.ch/' prefix
                       clearinputList.append(pathtodqmfiles + '/' + infi)
                       
                       print("="*40)  
                       #print clearinputList

           #Ready to create plots for each file
           #We should be at the reco_prodtools directory
           #for infi in clearinputList: 
           #print (' '.join(clearinputList))
           cmd = 'python ../Validation/HGCalValidation/scripts/makeHGCalValidationPlots.py ' + ' '.join(clearinputList) + ' --outputDir HGCValid_%s_Plots --no-ratio --png --separate --html-sample "' %(opt.HTMLVALNAME) + SampleName + ' E %s GeV %s " ' % (theEn,opt.TAG) + ' --html-validation-name %s --subdirprefix ' %(opt.HTMLVALNAME) + ' plots_PDGId%s_nPart%s_E%sEta%sto%s_%s_Delta_%s' % (opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG,"_".join(str(n) for n in Deltas).replace('.','p')) + ' --collection %s' %(opt.HTMLVALNAME)

           if(opt.DRYRUN):
               print('Dry-run: ['+cmd+']')
           else:
               output = processCmd(cmd)
               #processCmd( ' mv HGCalValidationPlots/index.html HGCalValidationPlots/index_Pt%sEta%sto%s.html'% (thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'))  )
               linetocuttheindexhtml = "44" #for layerclusters 44, for multiclusters 22 for the moment
               if "Multi" in opt.HTMLVALNAME: 
                   linetocuttheindexhtml = "28"
               processCmd('awk \'NR>=6&&NR<=%s\' HGCValid_%s_Plots/index.html > HGCValid_%s_Plots/%s/index_PDGId%s_nPart%s_E%sEta%sto%s_%s_Delta_%s.html '% (linetocuttheindexhtml,opt.HTMLVALNAME,opt.HTMLVALNAME,opt.TAG,opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG,"_".join(str(n) for n in Deltas).replace('.','p'))  )
                      
           #cmd = cmd + 'sed -i \'s/Tracking/HGCal/g\' HGCalValidationPlots/index_Pt%sEta%sto%s.html'% (thePt,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p')) + '\n '
           #ALWAYS CHECK THE FRAGMENT YOU WANT
           fragments.append( 'HGCValid_%s_Plots/%s/index_PDGId%s_nPart%s_E%sEta%sto%s_%s_Delta_%s.html'% (opt.HTMLVALNAME,opt.TAG,opt.PARTID.replace(',','_'),opt.NPART,theEn,str(Eta[0]).replace('.','p'),str(Eta[1]).replace('.','p'),opt.TAG,"_".join(str(n) for n in Deltas).replace('.','p')) )
   
    #We need the initial index.html file before creating our own. 
    processCmd('cp /eos/user/a/apsallid/www/index.php HGCValid_%s_Plots/.' %(opt.HTMLVALNAME))
    #The command above created the plots_.html files in the parent directory. We should move them
    processCmd('mv HGCValid_%s_Plots/plots_* HGCValid_%s_Plots/%s/.' %(opt.HTMLVALNAME,opt.HTMLVALNAME,opt.TAG)  )
    #Let's also create the final index xml file. 
    #There is also the index.html file we should move
    processCmd('mv HGCValid_%s_Plots/index.html HGCValid_%s_Plots/test.html' %(opt.HTMLVALNAME,opt.HTMLVALNAME) )
    index_file = open('HGCValid_%s_Plots/%s/index.html'%(opt.HTMLVALNAME,opt.TAG),'w')            
    #Write preamble
    index_file.write('<html>\n')
    index_file.write(' <head>\n')
    index_file.write('  <title>HGCal validation %s </title>\n' %(opt.HTMLVALNAME) )
    index_file.write(' </head>\n')
    index_file.write(' <body>\n')
    
    for frag in fragments:   
        with open(frag,'r') as f:
            lines = f.read().splitlines()
            for line in lines:
                print(line)
                index_file.write(line + '\n')
                #processCmd( 'cat ' + frag + ' >> HGCalValidationPlots/index.html '   )
                #index_file.write(frag)
        
        
    #Writing postamble"
    index_file.write(' </body>\n')
    index_file.write('</html>\n')
    index_file.close()
    
    processCmd( 'cp HGCValid_%s_Plots/%s/index.html HGCValid_%s_Plots/%s/index_%s.html ' %(opt.HTMLVALNAME,opt.TAG,opt.HTMLVALNAME,opt.TAG,opt.PARTID.replace(',','_'))   )
       



import networkx as nx
import ROOT
from array import array
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pylab
import plotly.graph_objects as go
import matplotlib.ticker as ticker

def histValue1D(fValues, histDict, tag="hist1D_", title="hist 1D", axunit="a.u.", binsBoundariesX=[10, -1, 1], ayunit="a.u.", verbosityLevel=0):
    """1D histograming of given list of values."""
    # sanity check for hists
    if histDict is None:
        return
    # sanity check for boundaries
    if (len(binsBoundariesX) != 3 and len(binsBoundariesX) != 2):
        return
    # define event-level hists
    elif len(binsBoundariesX) == 3:  # bondaries in format [nbins, low, high]
        histDict[tag] = ROOT.TH1F(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], binsBoundariesX[1], binsBoundariesX[2])
    elif len(binsBoundariesX) == 2:  # bondaries in format [nbins, list_boundaries]
        histDict[tag] = ROOT.TH1F(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], array('f', binsBoundariesX[1]))
    # set some properties
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset() * 3.0)
    # loop over all values
    if (verbosityLevel >= 3):
        print("tag: ", tag, ", fValues: ", fValues)
    for value in fValues:
        histDict[tag].Fill(value)
    return histDict


def histValues2D(fValues, histDict, tag="hist2D_", title="hist 2D", axunit="a.u.", binsBoundariesX=[10, -1, 1], ayunit="a.u.", binsBoundariesY=[10, -1, 1], weighted2D=False, verbosityLevel=0):
    """2D histograming of given list of values"""
    # sanity check for hists
    if histDict is None:
        return
    # sanity check for boundaries
    if (len(binsBoundariesX) != len(binsBoundariesY)):
        return
    if (len(binsBoundariesX) != 3 and len(binsBoundariesX) != 2):
        return
    # define event-level hists
    elif len(binsBoundariesX) == 3:  # bondaries in format [nbins, low, high]
        histDict[tag] = ROOT.TH2F(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], binsBoundariesX[1], binsBoundariesX[2], binsBoundariesY[0], binsBoundariesY[1], binsBoundariesY[2])
    elif len(binsBoundariesY) == 2:  # bondaries in format [nbins, list_boundaries]
        histDict[tag] = ROOT.TH2F(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], array('f', binsBoundariesX[1]), binsBoundariesY[0], array('f', binsBoundariesY[1]))
    # set some properties
    histDict[tag].GetXaxis().SetTitleOffset(histDict[tag].GetXaxis().GetTitleOffset() * 1.0)
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset() * 3.0)
    # loop over all values
    if (verbosityLevel >= 3):
        print("tag: ", tag, ", fValues: ", fValues)
    if (not weighted2D):
        for (valueX, valueY) in fValues:
            histDict[tag].Fill(valueX, valueY)
    else:
        for (valueX, valueY, valueZ) in fValues:
            print (valueX, valueY, valueZ)
            histDict[tag].Fill(valueX, valueY, valueZ)
    return histDict


def histsPrintSaveSameCanvas(histsAndProps, outDir, tag="hists1D_", latexComment="", funcsAndProps=None, verbosityLevel=0):
    """print/save list of histograms with their properties on one canvas"""
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1
    # set default style values
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.05)
    # create canvas
    canvas = ROOT.TCanvas(outDir + tag, outDir + tag, 500, 500)
    # prepare the legend
    leg = ROOT.TLegend(0.15, 0.90-len(histsAndProps)*0.07, 0.82, 0.9)
    # leg.SetHeader("Energy of the clusters before/after filtering")
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    # prepare latex comment
    ltx = ROOT.TLatex()
    ltx.SetNDC(ROOT.kTRUE)
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.03)
    # set image extensions
    imgTypes = ["pdf", "png", "root"]
    if (verbosityLevel >= 3):
        print("histsAndProps: ", histsAndProps)
        print("funcsAndProps: ", funcsAndProps)
    # loop over all histograms to get max
    y_maxs = [0.01]
    x_maxs = [1.]
    for hist in histsAndProps:
        # do not print/save empty histograms
        if (type(hist) == ROOT.TH1F) or (type(hist) == ROOT.TH2F) or (type(hist) == ROOT.TH3F):
            if hist.GetEntries() == 0:
                continue
            # if (type(hist) == ROOT.TH1F):
            #     hist.Rebin()
        x_maxs.append(hist.GetBinCenter(hist.FindLastBinAbove(1)))
        hist.Scale(1./hist.Integral())
        hist.GetYaxis().SetTitle("a.u.")
        curr_max = hist.GetMaximum()
        if (curr_max < 1./3.): # temp. fix for hists with very different y_max
            y_maxs.append(curr_max)
    # print "y_maxs: ", y_maxs
    # loop over all histograms
    first = True
    for hist in histsAndProps:
        # do not print/save empty histograms
        if (type(hist) == ROOT.TH1F) or (type(hist) == ROOT.TH2F) or (type(hist) == ROOT.TH3F):
            if hist.GetEntries() == 0:
                continue
        # print and save
        hist.SetTitle("")
        if type(hist) == ROOT.TH1F:
            hist.SetLineColor(histsAndProps[hist]["color"])
            hist.SetLineWidth(2)
            leg.AddEntry(hist, histsAndProps[hist]["leg"], "L")
            hist.GetXaxis().SetTitleOffset(hist.GetXaxis().GetTitleOffset() * 1.2)
            hist.GetXaxis().SetTitle("E_{meas}/E_{true}")
            hist.GetYaxis().SetTitleOffset(hist.GetYaxis().GetTitleOffset() * 3.0)
            if (first):
                hist.GetYaxis().SetRangeUser(0, max(y_maxs) * 1.4)
                hist.GetXaxis().SetRangeUser(0, max(x_maxs) * 1.0)
                hist.Draw("hist0 goff")
                first = False
            else:
                hist.Draw("hist0 same goff")
    # check if any function should be drawn
    if funcsAndProps is not None:
        for func in funcsAndProps:
            func.SetLineColor(funcsAndProps[func]["color"])
            leg.AddEntry(func, funcsAndProps[func]["leg"], "L")
            func.Draw("same goff")
    # draw the rest
    leg.Draw("same")
    # print latex comments
    ltx.SetTextColor(ROOT.kBlue)
    for k in range(len(latexComment)):
        ltx.DrawLatex(0.17, 0.86 - len(histsAndProps)*0.07 - k*0.07, latexComment[k])
    # print latex header
    ltx.SetTextColor(ROOT.kBlack)
    ltx.DrawLatex(0.150, 0.935, "CMS Phase-2 Simulation, #sqrt{s} = 14 TeV")
    for imgType in imgTypes:
        canvas.SaveAs("{}/{}.{}".format(outDir, tag, imgType))
    return canvas

#---------------------------------------------------------------------------------------------------
# print/save all histograms
def histPrintSaveAll(histDict, outDir, tree):
    imgType = "png"
    outfile = ROOT.TFile("{}/{}".format(outDir, options.output), "recreate")
    canvas = ROOT.TCanvas(outDir, outDir, 500, 500)
    if (options.verbosityLevel>=3): print( "histDict.items(): ", histDict.items())
    #Write tree
    tree.SetDirectory(outfile)
    tree.Write()
    for key, item in histDict.items():
        # do not save empty histograms
        if (type(item) == ROOT.TH1F) or (type(item) == ROOT.TH2F) or (type(item) == ROOT.TH3F):
            if item.GetEntries() == 0:
                continue
        ROOT.gStyle.SetPalette(ROOT.kBird)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPadTopMargin(0.05)
        ROOT.gStyle.SetPadBottomMargin(0.12)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadRightMargin(0.02)
        if type(item) == ROOT.TH1F:
            item.Draw("hist0")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TH2F:
            item.Draw("colz")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TProfile2D:
            item.Draw("")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            item.Draw("box")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        else:
            continue
    return

#---------------------------------------------------------------------------------------------------
def profValues2D(fValues, histDict, tag="hist2D_", title="hist 2D", axunit="a.u.", binsBoundariesX=[10, -1, 1], ayunit="a.u.", binsBoundariesY=[10, -1, 1], binsBoundariesZ=[-1., 1], weighted2D=False, weight=1, verbosityLevel=0):
    """2D profile of given list of values"""
    # sanity check for hists
    if histDict is None:
        return
    # sanity check for boundaries
    if (len(binsBoundariesX) != len(binsBoundariesY)):
        return
    if (len(binsBoundariesX) != 3 and len(binsBoundariesX) != 2):
        return
    # define event-level hists
    elif len(binsBoundariesX) == 3:  # bondaries in format [nbins, low, high]
        histDict[tag] = ROOT.TProfile2D(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], binsBoundariesX[1], binsBoundariesX[2], binsBoundariesY[0], binsBoundariesY[1], binsBoundariesY[2], binsBoundariesZ[0], binsBoundariesZ[1])
    # set some properties
    histDict[tag].GetXaxis().SetTitleOffset(histDict[tag].GetXaxis().GetTitleOffset() * 1.0)
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset() * 3.0)
    # loop over all values
    if (verbosityLevel >= 3):
        print ("tag: ", tag, ", fValues: ", fValues)
    if (not weighted2D):
        for (valueX, valueY, valueZ) in fValues:
            #print (valueX, valueY, valueZ)
            histDict[tag].Fill(valueX, valueY, valueZ)
    else:
        for (valueX, valueY, valueZ) in fValues:
            #print (valueX, valueY, valueZ)
            histDict[tag].Fill(valueX, valueY, valueZ, weight)
    return histDict

def drawGraphs(graphsAndProps, grOptions, outDir, latexComment=[], tag="graphTest_", verbosityLevel=0):
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1
    # set default style values
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.10)
    ROOT.gStyle.SetPadRightMargin(0.05)
    # create canvas
    canvas = ROOT.TCanvas(tag, tag, 800, 600)
    # prepare the legend (according to the number of entries)
    legLowerBoundary = 0.85-len(graphsAndProps)*0.07
    if (legLowerBoundary<0.45): legLowerBoundary = 0.45
    leg = ROOT.TLegend(0.40, legLowerBoundary, 0.95, 0.9)
    leg.SetHeader(grOptions['title'])
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    # prepare latex comment
    ltx = ROOT.TLatex()
    ltx.SetNDC(ROOT.kTRUE)
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.03)
    # set image extensions
    imgTypes = ["pdf", "png", "root"]
    # prepare latex comment
    ltx = ROOT.TLatex()
    ltx.SetNDC(ROOT.kTRUE)
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.03)
    if (verbosityLevel >= 3):
        print("graphsAndProps: ", graphsAndProps)
    # loop over all graphs to get max
    y_maxs = [gr.GetYaxis().GetXmax() for gr in graphsAndProps]
    y_mins = [gr.GetYaxis().GetXmin() for gr in graphsAndProps]
    if (verbosityLevel >= 3):
        print("y_mins: ", y_mins)
        print("y_maxs: ", y_maxs)
    # loop over all histograms
    first = True
    k = 0
    for gr in graphsAndProps:
        gr.SetTitle("")
        gr.GetXaxis().SetTitle(grOptions['Xaxis'])
        gr.GetXaxis().SetTitleSize(0.05)
        gr.GetXaxis().SetTitleOffset(0.9)
        gr.GetYaxis().SetTitle(grOptions['Yaxis'])
        gr.GetYaxis().SetTitleSize(0.05)
        gr.GetYaxis().SetTitleOffset(0.8)
        colour = graphsAndProps[gr]["color"]
        gr.SetLineColor(colour)
        gr.SetLineWidth(1)
        gr.SetLineStyle(graphsAndProps[gr]["LineStyle"])
        gr.SetMarkerColor(colour)
        gr.SetMarkerStyle(graphsAndProps[gr]["MarkerStyle"])
        gr.SetMarkerSize(0.7)
        gr.SetFillColor(0)
        gr.SetFillStyle(0)
        leg.AddEntry(gr, graphsAndProps[gr]["leg"])
        if (first):
            gr.SetMaximum(max(y_maxs) * 1.5)
            gr.SetMinimum(4.)
            gr.SetTitle("")
            gr.Draw("AP goff")
            first = False
        else:
            gr.Draw("P same goff")
        ltx.SetLineColor(colour)
        ltx.SetTextColor(colour)
        if (len(graphsAndProps)>2): ltx.SetTextSize(0.02)
        if 'latexComment' in graphsAndProps[gr].keys():
            ltxCommentSize = (legLowerBoundary - 0.2)/len(graphsAndProps)
            ltx.DrawLatex(0.45, legLowerBoundary - k*ltxCommentSize, graphsAndProps[gr]['latexComment'])
        k+=1
    # draw the rest
    leg.Draw("same")
    # print common header
    ltx.SetTextColor(ROOT.kBlack)
    ltx.DrawLatex(0.120, 0.935, "CMS Phase-2 Simulation, #sqrt{s} = 14 TeV")
    # save
    for imgType in imgTypes:
        canvas.SaveAs("{}/{}.{}".format(outDir, tag, imgType))
    return canvas

#---------------------------------------------------------------------------------------------------
# print/save all histograms
def histPrintSaveAll(histDict, outDir, output, tree, verbosityLevel = 0):
    imgType = "png"
    outfile = ROOT.TFile("{}/{}".format(outDir, output), "recreate")
    canvas = ROOT.TCanvas(outDir, outDir, 500, 500)
    if (verbosityLevel>=3): print( "histDict.items(): ", histDict.items())
    #Write tree
    tree.SetDirectory(outfile)
    tree.Write()
    for key, item in histDict.items():
        # do not save empty histograms
        if (type(item) == ROOT.TH1F) or (type(item) == ROOT.TH2F) or (type(item) == ROOT.TH3F):
            if item.GetEntries() == 0:
                continue
        ROOT.gStyle.SetPalette(ROOT.kBird)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPadTopMargin(0.05)
        ROOT.gStyle.SetPadBottomMargin(0.12)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadRightMargin(0.02)
        if type(item) == ROOT.TH1F:
            item.Draw("hist0")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TH2F:
            item.Draw("colz")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TProfile2D:
            item.Draw("")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            item.Draw("box")
            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        else:
            continue
    return

def fitGauss(hist, paramRangeFactor=1.8):
    if (hist.GetEntries() == 0):
        return (hist, 0, 0)
    hist.GetXaxis().SetTitleOffset(hist.GetXaxis().GetTitleOffset() * 1.2)
    hist.GetYaxis().SetTitleOffset(hist.GetYaxis().GetTitleOffset() * 3.0)
    # define the range of the fit from the hist mean and RMS
    meanLimitDn = hist.GetMean() - paramRangeFactor * hist.GetRMS()
    meanLimitUp = hist.GetMean() + paramRangeFactor * hist.GetRMS()
    sigmaLimitDn = hist.GetRMS() / paramRangeFactor
    sigmaLimitUp = hist.GetRMS() * paramRangeFactor
    # define the fitting gausian and range of its parameters
    fGauss = ROOT.TF1("f", "[0]*TMath::Gaus(x,[1],[2])", meanLimitDn, meanLimitUp)
    fGauss.SetParLimits(1, meanLimitDn, meanLimitUp)
    fGauss.SetParLimits(2, sigmaLimitDn, sigmaLimitUp)
    # perform fit and extract params
    hist.Fit(fGauss, "Q", "", meanLimitDn, meanLimitUp)
    gaussMean = fGauss.GetParameter(1)
    gaussStd = fGauss.GetParameter(2)
    return (hist, gaussMean, gaussStd)

def fitResolution(graph, fitLineColor = ROOT.kBlue, fitLineStyle = 1, rangeLimitDn = 5., rangeLimitUp = 100.):
    # define the range of the fit from the hist mean and RMS
    stochasticTermLimitDn = 0
    stochasticTermLimitUp = 300
    constantTermLimitDn = 0
    constantTermLimitUp = 100
    noiseTermLimitDn = 0
    noiseTermLimitUp = 500
    # define the fitting gausian and range of its parameters
    fResolution = ROOT.TF1("f", "sqrt([1]*[1] + [0]*[0]/x + [2]*[2]/(x*x))", rangeLimitDn, rangeLimitUp)
    fResolution.SetParLimits(0, stochasticTermLimitDn, stochasticTermLimitUp)  # stochastic term
    fResolution.SetParLimits(1, constantTermLimitDn, constantTermLimitUp)  # constant term
    fResolution.SetParLimits(2, noiseTermLimitDn, noiseTermLimitUp) # noise term
    fResolution.SetLineColor(fitLineColor)
    fResolution.SetLineStyle(fitLineStyle)
    # perform fit and extract params
    graph.Fit(fResolution, "", "", rangeLimitDn, rangeLimitUp)
    stochasticTerm = fResolution.GetParameter(0)
    constantTerm = fResolution.GetParameter(1)
    noiseTerm = fResolution.GetParameter(2)
    return (graph, stochasticTerm, constantTerm, noiseTerm)

def getEffSigma(theHist, wmin=-100, wmax=100, epsilon=0.01):
    """taken from Hgg framework (by Ed)"""
    # initialise
    weight = 0.
    points = []
    thesum = theHist.Integral()
    # return -1 in case of empty histogram
    if (thesum == 0):
        return -1.
    # compute the cumulative distr. points
    for i in range(theHist.GetNbinsX()):
        weight += theHist.GetBinContent(i)
        if weight / thesum > epsilon:
            points.append([theHist.GetBinCenter(i), weight / thesum])
    # initialise
    low = wmin
    high = wmax
    width = wmax - wmin
    # find minimal 0.683 interval
    for i in range(len(points)):
        for j in range(i, len(points)):
            wy = points[j][1] - points[i][1]
            if abs(wy - 0.683) < epsilon:
                wx = points[j][0] - points[i][0]
                if wx < width:
                    low = points[i][0]
                    high = points[j][0]
                    width = wx
    return 0.5 * (high - low)

def getHistMeanStd(histo):
    hEntries = histo.GetEntries()
    gMean = histo.GetMean()
    gMeanError = histo.GetMeanError()
    gStd = histo.GetRMS()
    if (hEntries > 100): # extract mean/error from fit if enough statistics
        (histo, gMean, gStd) = fitGauss(histo)
        gMeanError = gStd/(hEntries**0.5)
    effSigma = getEffSigma(histo)
    return histo, hEntries, gMean, gMeanError, gStd, effSigma

#---------------------------------------------------------------------------------------------------
def acustompalette():
    NRGBs = 7
    NCont = 100
    ncolors = array('i', [])
    ROOT.gStyle.SetNumberContours(NCont);
    stops   = [ 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 ]
    red     = [ 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 ]
    green   = [ 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 ]
    blue    = [ 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 ]
    stopsArray = array('d', stops)
    redArray = array('d', red)
    greenArray = array('d', green)
    blueArray = array('d', blue)
    first_color_number = ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, NCont);
    ROOT.gStyle.SetNumberContours(NCont)

    palsize = NCont
    palette = []
    for i in range(palsize):
        palette.append(first_color_number+i)
        palarray = array('i',palette)

    ROOT.gStyle.SetPalette(palsize,palarray)

#---------------------------------------------------------------------------------------------------
def isMatched(eta,phi, genParticles,minDR = 0.5):
    belong = -1	
    for iGen in genParticles:
        geta,gphi=iGen.eta(),iGen.phi()
        deta=geta-eta
        dphi=TVector2.Phi_mpi_pi(gphi-phi)
        dR=TMath.Sqrt(deta**2+dphi**2)
        if dR<minDR:
            belong=1
            break
    return belong
#---------------------------------------------------------------------------------------------------
def drawEdges3d_plt(df,EventId):
    
    fig = plt.figure(figsize=(12, 9))
    ax = plt.axes(projection ="3d")

    ax.scatter3D(df['trackster_edge_layerin'], df['trackster_edge_layerin_x'], df['trackster_edge_layerin_y'], c=df['trackster_id'])
    ax.scatter3D(df['trackster_edge_layerout'], df['trackster_edge_layerout_x'], df['trackster_edge_layerout_y'], c=df['trackster_id'])

    edgein = df[['trackster_edge_layerin','trackster_edge_layerin_x','trackster_edge_layerin_y']].to_numpy()
    edgeout = df[['trackster_edge_layerout','trackster_edge_layerout_x','trackster_edge_layerout_y']].to_numpy()

    theedges = np.concatenate((edgein, edgeout), axis=1)

    for i in range(0,len(df['trackster_edge_layerin'])):
        theedges_x = []
        theedges_x.append(theedges[i][1])
        theedges_x.append(theedges[i][4])
        theedges_y = []
        theedges_y.append(theedges[i][2])
        theedges_y.append(theedges[i][5])
        theedges_la = []
        theedges_la.append(theedges[i][0])
        theedges_la.append(theedges[i][3])
        theedges_R = []
        theedges_R.append(np.sqrt(theedges[i][1]**2 + theedges[i][2]**2))
        theedges_R.append(np.sqrt(theedges[i][4]**2 + theedges[i][5]**2))
        #print(theedges_x,theedges_y,theedges_la,theedges_R)
        plt.plot(theedges_la, theedges_x,theedges_y, 'ro-')

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.set_title("Tracksters Edges by CA - Event %s " % (EventId))
    ax.set_xlabel('Layer Number', fontweight ='bold')
    ax.tick_params(left=True,bottom=True,labelleft=True,labelbottom=True)
    ax.set_xlabel('Layer Number', fontweight ='bold')
    ax.set_ylabel('X', fontweight ='bold')
    ax.set_zlabel('Y', fontweight ='bold')
    #plt.show()
    ax.figure.savefig('Edges_X_Y_Layer_plt.png') 
 
#---------------------------------------------------------------------------------------------------
def drawEdges2d_plt(df,EventId):
    

    fig = plt.figure(figsize=(12, 9))
    ax = plt.axes()

    trackster_edge_layerin_R = np.sqrt( np.add( 
        np.power(df['trackster_edge_layerin_x'].to_numpy(), 2),
        np.power(df['trackster_edge_layerin_y'].to_numpy(), 2)) 
    )
    trackster_edge_layerout_R = np.sqrt( np.add( 
        np.power(df['trackster_edge_layerout_x'].to_numpy(), 2),
        np.power(df['trackster_edge_layerout_y'].to_numpy(), 2)) 
    )

    ax.scatter(df['trackster_edge_layerin'], trackster_edge_layerin_R)
    ax.scatter(df['trackster_edge_layerout'], trackster_edge_layerout_R)

    edgein = df[['trackster_edge_layerin','trackster_edge_layerin_x','trackster_edge_layerin_y']].to_numpy()
    edgeout = df[['trackster_edge_layerout','trackster_edge_layerout_x','trackster_edge_layerout_y']].to_numpy()

    theedges = np.concatenate((edgein, edgeout), axis=1)

    for i in range(0,len(df['trackster_edge_layerin'])):
        theedges_x = []
        theedges_x.append(theedges[i][1])
        theedges_x.append(theedges[i][4])
        theedges_y = []
        theedges_y.append(theedges[i][2])
        theedges_y.append(theedges[i][5])
        theedges_la = []
        theedges_la.append(theedges[i][0])
        theedges_la.append(theedges[i][3])
        theedges_R = []
        theedges_R.append(np.sqrt(theedges[i][1]**2 + theedges[i][2]**2))
        theedges_R.append(np.sqrt(theedges[i][4]**2 + theedges[i][5]**2))
        #print(theedges_x,theedges_y,theedges_la,theedges_R)
        #plt.plot(theedges_la, theedges_x,theedges_y, 'ro-')
        plt.plot(theedges_la, theedges_R, 'ro-')

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.set_title("Tracksters Edges by CA - Event %s " % (EventId))
    ax.set_xlabel('Layer Number', fontweight ='bold')
    ax.set_ylabel('R (cm)', fontweight ='bold')
    ax.tick_params(left=True,bottom=True,labelleft=True,labelbottom=True)
    #plt.show()
    ax.figure.savefig('Edges_R_vs_Layer_plt.png') 
 
#---------------------------------------------------------------------------------------------------
def drawDirectedGraph(df,EventId):
    
    fig = plt.figure(figsize=(12, 9))
    ax = plt.axes()

    theGraph = nx.from_pandas_edgelist(df,source='trackster_edge_layerin',target='trackster_edge_layerout', edge_attr=['trackster_edge','trackster_angle_alpha'], create_using=nx.DiGraph())

    #print(nx.is_directed(theGraph))
    #print(theGraph.edges())
    #print(theGraph.nodes())

    pos = {}
    df_fornode = df[['trackster_edge','trackster_edge_layerin','trackster_edge_layerin_R']]
    for i,nlrow in df_fornode.iterrows(): 
        tmpdict = nlrow[1:].to_dict()
        tmpname = nlrow['trackster_edge_layerin']
        pos[tmpname] = (tmpdict['trackster_edge_layerin'],tmpdict['trackster_edge_layerin_R'])
        #theGraph.nodes[tmpname]['pos'] = (tmpdict['trackster_edge_layerin'],tmpdict['trackster_edge_layerin_z'])
        #print(nlrow['trackster_edge'], tmpdict['trackster_edge_layerin'],tmpdict['trackster_edge_layerin_z'])
        #g.add_node( tmpname, pos = (tmpdict['trackster_edge_layerin'],tmpdict['trackster_edge_layerin_z']) )


    df_fornode = df[['trackster_edge','trackster_edge_layerout','trackster_edge_layerout_R']]
    for i,nlrow in df_fornode.iterrows(): 
        tmpdict = nlrow[1:].to_dict()
        tmpname = nlrow['trackster_edge_layerout']
        pos[tmpname] = (tmpdict['trackster_edge_layerout'],tmpdict['trackster_edge_layerout_R'])
        #theGraph.nodes[tmpname]['pos'] = (tmpdict['trackster_edge_layerout'],tmpdict['trackster_edge_layerout_z'])
        #print(nlrow['trackster_edge'], tmpdict['trackster_edge_layerin'],tmpdict['trackster_edge_layerin_z'])
        #g.add_node( tmpname, pos = (tmpdict['trackster_edge_layerin'],tmpdict['trackster_edge_layerin_z']) )
    
    #print(theGraph.nodes(data=True))

    #print(pos)

    edge_labels = nx.get_edge_attributes(theGraph,'trackster_edge')
    #edge_labels = nx.get_edge_attributes(theGraph,'trackster_angle_alpha')

    #nx.draw_networkx_edge_labels(theGraph, pos, edge_labels = edge_labels )
    #nx.draw(theGraph, pos, ax = ax )
    #nx.draw_networkx_labels(theGraph, pos)
    #nx.draw_networkx_edges(theGraph, pos, edge_color='r', arrows=True)

    nx.draw(theGraph, pos, node_size = 50, with_labels= False, edge_color='r', arrows=True)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.set_title("Tracksters Edges by CA - Event %s " % (EventId))
    ax.set_xlabel('Layer Number', fontweight ='bold')
    ax.set_ylabel('R (cm)', fontweight ='bold')
    ax.tick_params(left=True,bottom=True,labelleft=True,labelbottom=True)
    limits = plt.axis("on")
    plt.show()
    #pylab.show()
    ax.figure.savefig('Edges_R_vs_Layer_DiGraph.png', dpi=300) 

#---------------------------------------------------------------------------------------------------
def drawDirectedGraph3D_plotly(df,dfl,EventId): #For now EventId, when looping through events will use EventIdFromFile

    #Later create one graph for each trackster and label the nodes by trackster id
    theGraph = nx.from_pandas_edgelist(df,source='trackster_edge_layerin',target='trackster_edge_layerout', edge_attr=['trackster_edge','trackster_angle_alpha'], create_using=nx.DiGraph())

    #print(nx.is_directed(theGraph))
    #print(theGraph.edges())
    #print(theGraph.nodes())
    #print(theGraph.nodes(data=True))

    #Nodes. We want a 3d of x vs y vs layer. 
    x_nodes = dfl[['layerCluster_x']].to_numpy().reshape(-1).tolist()
    y_nodes = dfl[['layerCluster_y']].to_numpy().reshape(-1).tolist()
    z_nodes = dfl[['layerCluster_layer']].to_numpy().reshape(-1).tolist()
    #Color nodes according to trackster id
    color_nodes = dfl[['layerCluster_trackster_id']].to_numpy()
    #Trackster id plus 1 for color to avoid the white
    color_nodes = (color_nodes +1).reshape(-1).tolist()

    #print(x_nodes,y_nodes,z_nodes)

    edge_list = theGraph.edges()

    #Starting and ending coordinates of each edge. 
    x_edges=[]
    y_edges=[]
    z_edges=[]
    
    #We will use the edge dataframe to get these coordinates. Nice also to doublecheck the nodes from LC df.
    df_edges = df[['trackster_edge_layerin','trackster_edge_layerin_x','trackster_edge_layerin_y','trackster_edge_layerout','trackster_edge_layerout_x','trackster_edge_layerout_y']]
    #print(df_edges[ (df_edges['trackster_edge_layerin'] == 53) & (df_edges['trackster_edge_layerout'] == 55)]['trackster_edge_layerin_x'].values )

    #Remember here edge is defined above in the graph with specific source and target. When we will be asking 
    #e.g. for edge[0] we will be getting the trackster_edge_layerin. 
    for edge in edge_list:
        #Choose only the line we are interested
        tmpdf = df_edges[ (df_edges['trackster_edge_layerin'] == edge[0]) & (df_edges['trackster_edge_layerout'] == edge[1])]
        #[beginning,ending,None]
        #print(tmpdf)

        x_coords = [tmpdf['trackster_edge_layerin_x'].to_numpy()[0],tmpdf['trackster_edge_layerout_x'].to_numpy()[0],None]
        x_edges += x_coords

        y_coords = [tmpdf['trackster_edge_layerin_y'].to_numpy()[0],tmpdf['trackster_edge_layerout_y'].to_numpy()[0],None]
        y_edges += y_coords
        
        z_coords = [tmpdf['trackster_edge_layerin'].to_numpy()[0],tmpdf['trackster_edge_layerout'].to_numpy()[0],None]
        z_edges += z_coords
        

    #print(x_edges) 
    #print(y_edges) 
    #print(z_edges)

    #create a trace for the edges
    trace_edges = go.Scatter3d(x=x_edges,
                               y=y_edges,
                               z=z_edges,
                               mode='lines',
                               line=dict(color='black', width=2),
                               hoverinfo='none')
    
    #create a trace for the nodes
    trace_nodes = go.Scatter3d(x=x_nodes,
                               y=y_nodes,
                               z=z_nodes,
                               mode='markers',
                               marker=dict(symbol='circle',
                                           size=10,
                                           color=color_nodes, #color the nodes according to trackster id
                                           line=dict(color='black', width=0.5)),
                               hoverinfo='z')

    #we need to set the axis for the plot 
    x_axis = dict(showbackground=True,
                  showline=True,
                  zeroline=True,
                  showgrid=True,
                  showticklabels=True,
                  range=[-200,200],
                  title='X')

    y_axis = dict(showbackground=True,
                  showline=True,
                  zeroline=True,
                  showgrid=True,
                  showticklabels=True,
                  range=[-200,200],                      
                  title='Y')

    z_axis = dict(showbackground=True,
                  showline=True,
                  zeroline=True,
                  showgrid=True,
                  showticklabels=True,
                  range=[50,100],
                  title='Layer number')

    #also need to create the layout for our plot
    layout = go.Layout(title="Tracksters Edges by CA Event %s " % (EventId),
                       #width=650,
                       #height=625,
                       showlegend=False,
                       scene=dict(xaxis=dict(x_axis),
                                  yaxis=dict(y_axis),
                                  zaxis=dict(z_axis),
                              ),
                       margin=dict(t=100),
                       hovermode='closest')


    #Include the traces we want to plot and create a figure
    data = [trace_edges, trace_nodes]
    fig = go.Figure(data=data, layout=layout)

    fig.write_html('first_figure.html', auto_open=False)
    #fig.show()



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
import math
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import dEdX
import settings

#Mainly for the colors:
#Instead of Green -> 125
#Instead of Blue -> 122
#Instead of Red -> 132
#Instead of Orange -> 131
#Instead of Magenta -> 123
#Instead of Black -> 129

ROOT.gROOT.ProcessLine(".x rootlogon.C")

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


def histsPrintSaveSameCanvas(histsAndProps, outDir, tag="hists1D_", xaxistitle = "", yaxistitle = "", setLogX = False, setLogY = True, latexComment="", funcsAndProps=None, verbosityLevel=0, ratioplot = False, systematics = False):
    """print/save list of histograms with their properties on one canvas"""
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1
    # set default style values
    #ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetPadTopMargin(0.08)
    #ROOT.gStyle.SetPadBottomMargin(0.12)
    #ROOT.gStyle.SetPadLeftMargin(0.12)
    #ROOT.gStyle.SetPadRightMargin(0.05)
    # create canvas
    canvas = None
    if ratioplot:
        canvas = makeRatioPlotCanvas(name = tag)
        canvas.cd(1)
        if setLogX: ROOT.gPad.SetLogx()
        if setLogY: ROOT.gPad.SetLogy()
    else:
        canvas = ROOT.TCanvas(outDir + tag, outDir + tag, settings.canvas_width, settings.canvas_height-150)
        canvas.cd()

    # prepare the legend
    leg = ROOT.TLegend(0.55, 0.90-len(histsAndProps)*0.07, 0.92, 0.9)
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
    #imgTypes = ["pdf", "png", "root"]
    imgTypes = ["png"]
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
        #customizeHisto(hist, ratioplot)
        x_maxs.append(hist.GetBinCenter(hist.FindLastBinAbove(1)))
        if hist.Integral() != 0 : hist.Scale(1./hist.Integral())
        hist.GetYaxis().SetTitle(yaxistitle)
        curr_max = hist.GetMaximum()
        if (curr_max < 1./3.): # temp. fix for hists with very different y_max
            y_maxs.append(curr_max)
    #print ("y_maxs",y_maxs)
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
            hist.GetXaxis().SetTitle(xaxistitle)
            hist.GetYaxis().SetTitleOffset(hist.GetYaxis().GetTitleOffset() * 3.0)
            #hist.GetYaxis().SetLabelSize(0.);
            
            if (first):
                if not setLogY:
                    hist.GetYaxis().SetRangeUser(0, max(y_maxs) * 1.4)
                    hist.GetXaxis().SetRangeUser(0, max(x_maxs) * 1.0)
                #else:
                #    hist.SetMaximum( max(y_maxs) * 10. );
                #    #hist.GetYaxis().SetRangeUser(1, pow(10., (max(y_maxs)) * 1.4))
                #    hist.GetXaxis().SetRangeUser(0, max(x_maxs) * 1.0)
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
    ltx.SetTextColor(122)
    for k in range(len(latexComment)):
        ltx.DrawLatex(0.17, 0.86 - len(histsAndProps)*0.07 - k*0.07, latexComment[k])
    # print latex header
    ltx.SetTextColor(129)
    ltx.DrawLatex(0.150, 0.935, "CMS Phase-2 Simulation, #sqrt{s} = 14 TeV")
    canvas.cd()
    
    if ratioplot:
        canvas.cd(2)
        ROOT.gStyle.SetOptStat(0)
        if settings.ratio_plot_grid :
            ROOT.gPad.SetGridy()
            ROOT.gPad.SetGridx()
        first = True
        icount = 0
        refHist = None
        for hist in histsAndProps:
            # do not print/save empty histograms
            if (type(hist) == ROOT.TH1F) or (type(hist) == ROOT.TH2F) or (type(hist) == ROOT.TH3F):
                if hist.GetEntries() == 0:
                    continue
            # print and save
            hist.SetTitle("")
            if type(hist) == ROOT.TH1F and first == True :
            #if type(hist) == ROOT.TH1F and icount == 0 :
                refHist = hist
                first = False
                icount += 1
                continue
            if type(hist) == ROOT.TH1F and first == False :
            #if type(hist) == ROOT.TH1F and icount != 0 :
                #(errorHist,systHist) = make_stat_progression(hist, systematic_only=systematics)
                #ROOT.SetOwnership(errorHist,0)
                #ROOT.SetOwnership(systHist ,0)
                #errorHist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
                #errorHist.GetYaxis().SetTitle('LC_{2,3,4,5}/LC_{1}')
                #errorHist.GetYaxis().CenterTitle(True)
                #systHist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
                #systHist.GetYaxis().SetTitle('LC_{2,3,4,5}/LC_{1}')
                #systHist.GetYaxis().CenterTitle(True)
                #errorHist.SetLineColor(histsAndProps[hist]["color"])
                #errorHist.SetLineWidth(2)
                #errorHist.SetStats(0)
                print("icount",icount)
                #if icount == 1: errorHist.Draw('E2')
                #else: errorHist.Draw('E2,same')
                #if systematics == True:
                #    systHist.SetLineColor(histsAndProps[hist]["color"])
                #    systHist.SetLineWidth(2)
                #    systHist.SetStats(0)
                #    systHist.Draw('E2,same')

            
                ratioHist = makeRatio(hist1 = hist, hist2 = refHist, isdata = True)
                #ROOT.SetOwnership(ratioHist,0)
                #ROOT.SetOwnership(line,0)
                customizeHisto(ratioHist, ratioplot)
                ratioHist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
                ratioHist.GetYaxis().SetTitle('LC_{2,3,4,5}/LC_{1}')
                ratioHist.GetYaxis().CenterTitle(True)
                ratioHist.SetStats(0)
                
                
                if icount == 1:

                    ratioHist.Draw('')
                    line = ROOT.TLine(ratioHist.GetXaxis().GetXmin(),1,ratioHist.GetXaxis().GetXmax(),1)
                    line.SetLineColor(4)
                    line.SetLineStyle(7)
                    line.Draw('same')

                    (errorHist,systHist) = make_stat_progression(hist, systematic_only=systematics)
                    #ROOT.SetOwnership(errorHist,0)
                    #ROOT.SetOwnership(systHist ,0)
                    systHist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
                    systHist.GetYaxis().SetTitle('LC_{2,3,4,5}/LC_{1}')
                    systHist.GetYaxis().CenterTitle(True)
                    errorHist.SetLineColor(histsAndProps[hist]["color"])
                    errorHist.SetLineWidth(2)
                    errorHist.SetStats(0)
                    errorHist.Draw('E2,same')
                    if systematics == True:
                        systHist.SetLineColor(histsAndProps[hist]["color"])
                        systHist.SetLineWidth(2)
                        systHist.SetStats(0)
                        systHist.Draw('E2,same')                   

                else:
                    ratioHist.Draw('same')
                icount += 1
  
    if setLogX: canvas.SetLogx()
    if setLogY: canvas.SetLogy()
    for imgType in imgTypes:
        canvas.SaveAs("{}/{}.{}".format(outDir, tag, imgType))
    return canvas

#---------------------------------------------------------------------------------------------------
# print/save all histograms
def histPrintSaveAll(histDict, outDir, tree):
    imgType = "png"
    outfile = ROOT.TFile("{}/{}".format(outDir, "test"), "recreate")
    canvas = ROOT.TCanvas(outDir, outDir, 500, 500)
    #canvas.SetLogy()
    if (options.verbosityLevel>=3): print( "histDict.items(): ", histDict.items())
    #Write tree
    tree.SetDirectory(outfile)
    #outfile.cd()
    tree.Write()
    for key, item in histDict.items():
        # do not save empty histograms
        if (type(item) == ROOT.TH1F) or (type(item) == ROOT.TH2F) or (type(item) == ROOT.TH3F):
            if item.GetEntries() == 0:
                continue
        ROOT.gStyle.SetPalette(ROOT.kBird)
        ROOT.gStyle.SetOptStat(1101)
        ROOT.gStyle.SetPadTopMargin(0.05)
        ROOT.gStyle.SetPadBottomMargin(0.12)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadRightMargin(0.02)
        if type(item) == ROOT.TH1F:
            item.Draw("hist0")
            #item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TH2F:
            item.Draw("colz")
            #item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TProfile2D:
            item.Draw("")
            #item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            item.Draw("box")
            #item.Write()
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

def drawGraphs(graphsAndProps, grOptions, outDir, latexComment=[], setLogY = False, tag="graphTest_", verbosityLevel=0):
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
    #leg.SetHeader(grOptions['title'])
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
    #imgTypes = ["pdf", "png", "root"]
    imgTypes = ["png"]
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
    ltx.SetTextColor(129)
    ltx.DrawLatex(0.120, 0.935, "CMS Phase-2 Simulation, #sqrt{s} = 14 TeV")
    if setLogY: canvas.SetLogy()
    # save
    for imgType in imgTypes:
        canvas.SaveAs("{}/{}.{}".format(outDir, tag, imgType))
    return canvas

#---------------------------------------------------------------------------------------------------
# print/save all histograms
def histPrintSaveAll(histDict, outDir, output, tree, verbosityLevel = 0, setLogY = False, removeOptStat = False, setLogz = True, outfilecounter = 0):
    imgType = "png"
    #outfile = ROOT.TFile("{}/{}.root".format(outDir, output), "recreate")
    canvas = ROOT.TCanvas(outDir, outDir, 500, 500)
    if setLogY: canvas.SetLogy()
    if (verbosityLevel>=3): print( "histDict.items(): ", histDict.items())
    #Write tree
    #tree.SetDirectory(outfile)
    #tree.Write()
    for key, item in histDict.items():
        # do not save empty histograms
        if (type(item) == ROOT.TH1F) or (type(item) == ROOT.TH2F) or (type(item) == ROOT.TH3F):
            if item.GetEntries() == 0:
                continue
        if removeOptStat: ROOT.gStyle.SetOptStat(0)
        else : ROOT.gStyle.SetOptStat(1101)
        #ROOT.gStyle.SetPalette(ROOT.kBird)
        #ROOT.gStyle.SetPadTopMargin(0.05)
        #ROOT.gStyle.SetPadBottomMargin(0.12)
        #ROOT.gStyle.SetPadLeftMargin(0.12)
        #ROOT.gStyle.SetPadRightMargin(0.12)
        if type(item) == ROOT.TH1F:
            item.Draw("hist0")
            #item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TH2F:
            acustompalette()
            ex1 = ROOT.TExec("ex1","acustompalette();");
            ex1.Draw();
            item.Draw("colz")
            if setLogz == True: canvas.SetLogz()
            canvas.Update()

            palette = item.GetListOfFunctions().FindObject("palette")
            if palette:
                palette.__class__ = ROOT.TPaletteAxis
                palette.SetX1NDC(0.89)
                palette.SetX2NDC(0.94)
                #palette.SetY1NDC(0.1)
                #palette.SetY2NDC(0.6)
                palette.GetAxis().SetTickSize(.01)
                #palette.GetAxis().SetTitle("Si thick")
                palette.GetAxis().SetTitleOffset(0.82);
                #palette.GetAxis().LabelsOption("v")
                ROOT.gPad.Update()

            item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TProfile2D:
            acustompalette()
            ex1 = ROOT.TExec("ex1","acustompalette();");
            ex1.Draw();
            item.Draw("")

            palette = item.GetListOfFunctions().FindObject("palette")
            if palette:
                palette.__class__ = ROOT.TPaletteAxis
                palette.SetX1NDC(0.89)
                palette.SetX2NDC(0.94)
                #palette.SetY1NDC(0.1)
                #palette.SetY2NDC(0.6)
                palette.GetAxis().SetTickSize(.01)
                #palette.GetAxis().SetTitle("Si thick")
                palette.GetAxis().SetTitleOffset(0.82);
                #palette.GetAxis().LabelsOption("v")
                ROOT.gPad.Update()
            
            #item.Write()
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            item.Draw("box")
            #item.Write()
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

def fitResolution(graph, fitLineColor = 122, fitLineStyle = 1, rangeLimitDn = 5., rangeLimitUp = 100.):
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
def make_stat_progression(myHisto,systematics={},
                          systematic_only=True,
                          combine_with_systematic=True):
    #This function returns a function with the statistical precision in each bin
    

    statPrecision = myHisto.Clone('_ratioErrors_')
    systPrecision = myHisto.Clone('_ratioSysErrors_')
    statPrecision.SetFillColorAlpha(settings.ratio_error_band_color,settings.ratio_error_band_opacity)
    statPrecision.SetFillStyle(settings.ratio_error_band_style)
    statPrecision.SetMarkerColorAlpha(0,0)

    systPrecision.SetFillColorAlpha(settings.ratio_syst_band_color,settings.ratio_error_band_opacity)
    systPrecision.SetFillStyle(settings.ratio_syst_band_style)
    systPrecision.SetMarkerColorAlpha(0,0)


    if len(systematics)==0 : systematic_only = False
    for ibin in range(myHisto.GetNbinsX()+1):
        y    = statPrecision.GetBinContent(ibin)
        stat = statPrecision.GetBinError  (ibin)
        if( y > 0 ):
            statPrecision.SetBinContent(ibin,      1 )
            statPrecision.SetBinError  (ibin, stat/y )
        else:
            statPrecision.SetBinContent(ibin,   1 )
            statPrecision.SetBinError  (ibin,   0 )
        if systematic_only:
            up_err_sum2 = 0
            dw_err_sum2 = 0
            if( y > 0 ):
                up_err_sum2 = (stat/y)*(stat/y)
                dw_err_sum2 = (stat/y)*(stat/y)
                for key,syst in systematics.items():
                    up_diff   = (syst.up_histo.GetBinContent  (ibin)- y)/y
                    dw_diff   = (syst.down_histo.GetBinContent(ibin)- y)/y
                    if( up_diff > 0 ):
                        up_err_sum2  += up_diff*up_diff
                    if( dw_diff < 0 ):
                        dw_err_sum2  += dw_diff*dw_diff
            up_error = math.sqrt(up_err_sum2)
            dw_error = math.sqrt(dw_err_sum2)
            band_max   = 1 + up_error
            band_min   = 1 - dw_error
            
            systPrecision.SetBinContent(ibin, (band_max + band_min)/2.0);
            systPrecision.SetBinError  (ibin, (band_max - band_min)/2.0);
    statPrecision.GetYaxis().SetRangeUser(0.1, 5.0)
    systPrecision.GetYaxis().SetRangeUser(0.1, 5.0)
    return (statPrecision, systPrecision)

#---------------------------------------------------------------------------------------------------
def makeRatio(hist1,hist2,ymax=2.1,ymin=0,norm=False, isdata =False):
    """returns the ratio plot hist2/hist1
    if one of the histograms is a stack put it in as argument 2!"""
    if norm:
        try:
            hist1.Scale(1/hist1.Integral())
            hist2.Scale(1/hist2.Integral())
        except(ZeroDivisionError):
            pass
    retH = hist1.Clone()
    retH.Divide(hist2)
    if isdata:
        for ibin in range(hist2.GetNbinsX()+1):
            ymc  = hist2.GetBinContent(ibin);
            stat = hist1.GetBinError  (ibin);
            if (ymc>0):
                retH.SetBinError  (ibin,stat/ymc);
            else:
                retH.SetBinError  (ibin,0);
    ROOT.SetOwnership(retH,0)
    return retH

#---------------------------------------------------------------------------------------------------
def makeRatioPlotCanvas(name=''):
    """
    returns a divided canvas for ratios
    """
    canv  = ROOT.TCanvas("c_" + name, name,settings.canvas_width,settings.canvas_height)
    canv.cd()
    #padup = ROOT.TPad("padup", "padup", 0, 0.4, 1, 1.0)
    padup = ROOT.TPad("padup", "padup", 0, 0.3, 1, 1.0)
    padup.SetNumber(1)
    #paddw = ROOT.TPad("paddw", "paddw", 0, 0.0, 1, 0.4)
    paddw = ROOT.TPad("paddw", "paddw", 0, 0.0, 1, 0.3)
    paddw.SetNumber(2)
    padup.Draw()
    padup.SetTopMargin(0.08)
    padup.SetBottomMargin(0.00)
    padup.SetLeftMargin(0.14)
    padup.SetRightMargin(0.05)
    padup.SetFrameBorderMode(0)
    paddw.Draw()
    paddw.SetTopMargin(0.00)
    paddw.SetBottomMargin(0.37)
    paddw.SetLeftMargin(0.14)
    paddw.SetRightMargin(0.05)
    paddw.SetFrameBorderMode(0)
    canv.cd()
    #ROOT.SetOwnership(padup,0)
    #ROOT.SetOwnership(paddw,0)
    return canv

#---------------------------------------------------------------------------------------------------
def customizeHisto(hist, ratioplot = True):
    hist.GetYaxis().SetTitleSize  (21)
    hist.GetYaxis().SetTitleFont  (43)
    hist.GetYaxis().SetTitleOffset(1.8)
    hist.GetYaxis().SetLabelFont  (43)
    hist.GetYaxis().SetLabelSize  (18)
    if ratioplot :
        hist.GetXaxis().SetTitleSize  (21)
        hist.GetXaxis().SetTitleFont  (43)
        hist.GetXaxis().SetTitleOffset(3.5)
        hist.GetXaxis().SetLabelOffset(0.02)
        hist.GetXaxis().SetLabelFont  (43)
        hist.GetXaxis().SetLabelSize  (18)
    else:
        hist.GetXaxis().SetTitleSize  (21)
        hist.GetXaxis().SetTitleFont  (43)
        hist.GetXaxis().SetTitleOffset(1.5)
        hist.GetXaxis().SetLabelOffset(0.01)
        hist.GetXaxis().SetLabelFont  (43)
        hist.GetXaxis().SetLabelSize  (18)

#---------------------------------------------------------------------------------------------------
def DataMCratio(histMC,histData,
                log=False,
                xTitle="",
                yTitle="",
                drawMCOpt="",
                drawDataOpt="",
                norm=False,
                ratioMin=0.7,
                ratioMax=1.3):
    
    """Takes two histograms as inputs and returns a canvas with a ratio plot of the two.
    The two optional arguments are for the x Axis and y Axis titles"""

    c = makeRatioPlotCanvas()
    pad1 = c.cd(1)
    pad2 = c.cd(2)
    c.cd(1)

    if log: pad1.SetLogy()
    yMax = max(histData.GetMaximum(),histMC.GetMaximum())
    yMin = 0

    if log: yMin = min(histData.GetMinimum(),histMC.GetMinimum())
    else  : yMin = 0
    if log: yMax = 100*yMax
    else  : yMax = 1.2*yMax

    try:
        histData.GetYaxis().SetRangeUser(0,1.2*yMax)
    except(ReferenceError):
        h = pad1.DrawFrame(histMC.GetXaxis().GetXmin(),yMin,histMC.GetXaxis().GetXmax(),yMax)
        ROOT.SetOwnership(h,0)
    if not norm:
        drawStatErrBand(histMC,drawMCOpt)
        histData.Draw  ('same,'+drawDataOpt)
    else:
        histMC   = histMC.DrawNormalized(drawMCOpt)
        histData = histData.DrawNormalized("same"+ drawDataOpt)

    histData.GetYaxis().SetRangeUser(yMin,yMax)
    c.cd()
    c.cd(2)

    (errorHist,systHist) = make_stat_progression(histMC)
    ROOT.SetOwnership(errorHist,0)
    ROOT.SetOwnership(systHist ,0)
    errorHist.GetXaxis().SetTitle(xTitle)
    errorHist.GetYaxis().SetTitle(yTitle)
    #
    errorHist.Draw('E2')
    sysrHist.Draw('E2,same')
    ratioHist = makeRatio(histData,histMC,ymax= ratioMax,ymin=ratioMin,norm=norm)
    ROOT.SetOwnership(ratioHist,0)
    ratioHist.GetXaxis().SetTitle(xTitle)
    ratioHist.GetYaxis().SetTitle(yTitle)
    
    line = ROOT.TLine(ratioHist.GetXaxis().GetXmin(),1,ratioHist.GetXaxis().GetXmax(),1)
    line.SetLineColor(4)
    line.Draw()
    ROOT.SetOwnership(line,0)
    ratioHist.Draw('same')
    c.cd()
    return c

#---------------------------------------------------------------------------------------------------
def layerClusterPlots(df,dfl,tree,maxEvents,outDir,output,GenEnergy,verbosityLevel = 0):

    #-------------------------------------------------------------------------
    #Sum of LCs energy over generated energy
    ddfl = dfl.loc[ (dfl['nhitAll'] >= 1)].copy()
    sumLCene = ddfl.groupby(['EventId']).agg( LCEneSum  = ('energy','sum'))
    #print(sumLCene.head())
    #print(sumLCene[['LCEneSum']].to_numpy())
    theEgen = np.full(len(sumLCene[['LCEneSum']].to_numpy().flatten()),GenEnergy) 
    sumLCeneOverEgen = sumLCene[['LCEneSum']].to_numpy().flatten() / theEgen
    print(sumLCeneOverEgen,sumLCene[['LCEneSum']],theEgen)

    histDict_sumLCene = {}
    histDict_sumLCene = histValue1D(sumLCeneOverEgen, histDict_sumLCene, tag = "sumLCeneOverEgen", title = "Cumulative LCs energy over generated energy",   axunit = "#sum #left(LCs Energy#right)/Egen", binsBoundariesX = [200, 0, 2], ayunit = "#Events", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_sumLCene, outDir, output, tree, verbosityLevel)

    #-------------------------------------------------------------------------
    #LC energy over LC size: So, average cell energy
    histDict_LCeneOverLCsize = {}
    LCeneOverLCsize = ddfl[['energy']].to_numpy() / ddfl[['nhitAll']].to_numpy()
    LCeneOverLCsize = LCeneOverLCsize.flatten()
    histDict_LCeneOverLCsize = histValue1D(LCeneOverLCsize, histDict_LCeneOverLCsize, tag = "LCeneOverLCsize", title = "LC energy over LC size",   axunit = "Average Cell Energy (GeV)", binsBoundariesX = [200, 0, 2], ayunit = "#Events/0.01 GeV", verbosityLevel=verbosityLevel)
    #Plus a zoomed version
    histDict_LCeneOverLCsize = histValue1D(LCeneOverLCsize, histDict_LCeneOverLCsize, tag = "LCeneOverLCsize_Zoomed", title = "LC energy over LC size",   axunit = "Average Cell Energy (GeV)", binsBoundariesX = [20, 0, 0.2], ayunit = "#Events/0.01 GeV", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_LCeneOverLCsize, outDir, output, tree, verbosityLevel, setLogY = True)

    #-------------------------------------------------------------------------
    #LC size
    histDict_LCsize = {}
    LCsize = ddfl[['nhitAll']].to_numpy()
    LCsize = LCsize.flatten()
    histDict_LCsize = histValue1D(LCsize, histDict_LCsize, tag = "LCsize", title = "LC size",   axunit = "LC Size", binsBoundariesX = [50, 0, 50], ayunit = "#Events x #LCs", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_LCsize, outDir, output, tree, verbosityLevel, setLogY = True)

    #-------------------------------------------------------------------------
    #LC energy over LC size Per size: So, average cell energy per size. 
    histDict_LCeneOverLCsizePerSize = {}
    LCeneOverLCsizePerSize1 = ddfl[ (ddfl['nhitAll'] == 1)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] == 1)][['nhitAll']].to_numpy()
    LCeneOverLCsizePerSize2 = ddfl[ (ddfl['nhitAll'] == 2)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] == 2)][['nhitAll']].to_numpy()
    LCeneOverLCsizePerSize3 = ddfl[ (ddfl['nhitAll'] == 3)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] == 3)][['nhitAll']].to_numpy()
    LCeneOverLCsizePerSize4 = ddfl[ (ddfl['nhitAll'] == 4)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] == 4)][['nhitAll']].to_numpy()
    LCeneOverLCsizePerSizege5 = ddfl[ (ddfl['nhitAll'] >= 5)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] >= 5)][['nhitAll']].to_numpy()

    for i, obj in enumerate([LCeneOverLCsizePerSize1,LCeneOverLCsizePerSize2,LCeneOverLCsizePerSize3,LCeneOverLCsizePerSize4,LCeneOverLCsizePerSizege5]):
        histDict_LCeneOverLCsizePerSize[i] = {}
        thetitle = "LC energy over LC size from LCs of size %s" %(i+1)
        if i == 4: thetitle = "LC energy over LC size from LCs of size >= %s" %(i+1)
        histDict_LCeneOverLCsizePerSize[i] = histValue1D(obj, histDict_LCeneOverLCsizePerSize[i], tag = "LCeneOverLCsizeForLCSize%s"%(i+1), title = thetitle,  axunit = "Average Cell Energy (GeV)", binsBoundariesX = [200, 0, 2], ayunit = "#Events/0.01 GeV ", verbosityLevel=verbosityLevel)
        histDict_LCeneOverLCsizePerSize[i] = histValue1D(obj, histDict_LCeneOverLCsizePerSize[i], tag = "LCeneOverLCsizeForLCSize%s_Zoomed"%(i+1), title = thetitle,  axunit = "Average Cell Energy (GeV)", binsBoundariesX = [20, 0, 0.2], ayunit = "#Events/0.01 GeV ", verbosityLevel=verbosityLevel)
        print(histDict_LCeneOverLCsizePerSize[i])
        histPrintSaveAll(histDict_LCeneOverLCsizePerSize[i], outDir, output, tree, verbosityLevel)

    histsAndProps = {histDict_LCeneOverLCsizePerSize[0]['LCeneOverLCsizeForLCSize1']:{"leg":"LC size == 1","color":132}, histDict_LCeneOverLCsizePerSize[1]['LCeneOverLCsizeForLCSize2']:{"leg":"LC size == 2","color":122}, histDict_LCeneOverLCsizePerSize[2]['LCeneOverLCsizeForLCSize3']:{"leg":"LC size == 3","color":125}, histDict_LCeneOverLCsizePerSize[3]['LCeneOverLCsizeForLCSize4']:{"leg":"LC size == 4","color":131}, histDict_LCeneOverLCsizePerSize[4]['LCeneOverLCsizeForLCSize5']:{"leg":"LC size >= 5","color":129} }
    
    # plot these histograms on top of each other 
    histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "LCeneOverLCsizeNormForLCVariousSizes", xaxistitle = "Average Cell Energy (GeV)", yaxistitle = "a.u./0.01 GeV ", setLogY = True, latexComment = "", funcsAndProps = None, ratioplot = True)

    # make a zoomed version
    histsAndProps = {histDict_LCeneOverLCsizePerSize[0]['LCeneOverLCsizeForLCSize1_Zoomed']:{"leg":"LC size == 1","color":132}, histDict_LCeneOverLCsizePerSize[1]['LCeneOverLCsizeForLCSize2_Zoomed']:{"leg":"LC size == 2","color":122}, histDict_LCeneOverLCsizePerSize[2]['LCeneOverLCsizeForLCSize3_Zoomed']:{"leg":"LC size == 3","color":125}, histDict_LCeneOverLCsizePerSize[3]['LCeneOverLCsizeForLCSize4_Zoomed']:{"leg":"LC size == 4","color":131}, histDict_LCeneOverLCsizePerSize[4]['LCeneOverLCsizeForLCSize5_Zoomed']:{"leg":"LC size >= 5","color":129} }
    histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "LCeneOverLCsizeNormForLCVariousSizes_Zoomed", xaxistitle = "Average Cell Energy (GeV)", yaxistitle = "a.u./0.01 GeV ", setLogY = True, latexComment = "", funcsAndProps = None, ratioplot = True)

    
    #-------------------------------------------------------------------------
    #LC energy over LC size Per size vs Layer: So, average cell energy per size vs Layer. 
    histDict_LCeneOverLCsizePerSizeVsLayer = {}
    #LCeneOverLCsizePerSize1 = ddfl[ (ddfl['nhitAll'] == 1)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] == 1)][['nhitAll']].to_numpy()
    #LCeneOverLCsizePerSize2 = ddfl[ (ddfl['nhitAll'] == 2)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] == 2)][['nhitAll']].to_numpy()
    #LCeneOverLCsizePerSize3 = ddfl[ (ddfl['nhitAll'] == 3)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] == 3)][['nhitAll']].to_numpy()
    #LCeneOverLCsizePerSize4 = ddfl[ (ddfl['nhitAll'] == 4)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] == 4)][['nhitAll']].to_numpy()
    #LCeneOverLCsizePerSizege5 = ddfl[ (ddfl['nhitAll'] >= 5)][['energy']].to_numpy() / ddfl[ (ddfl['nhitAll'] >= 5)][['nhitAll']].to_numpy()

    LayerNumberPerSize1 = ddfl[ (ddfl['nhitAll'] == 1)][['layer']].to_numpy()
    LayerNumberPerSize2 = ddfl[ (ddfl['nhitAll'] == 2)][['layer']].to_numpy()
    LayerNumberPerSize3 = ddfl[ (ddfl['nhitAll'] == 3)][['layer']].to_numpy()
    LayerNumberPerSize4 = ddfl[ (ddfl['nhitAll'] == 4)][['layer']].to_numpy()
    LayerNumberPerSizege5 = ddfl[ (ddfl['nhitAll'] >= 5)][['layer']].to_numpy()

    LCeneOverLCsizePerSize1vsLayer = np.column_stack((LayerNumberPerSize1, LCeneOverLCsizePerSize1))
    LCeneOverLCsizePerSize2vsLayer = np.column_stack((LayerNumberPerSize2, LCeneOverLCsizePerSize2))
    LCeneOverLCsizePerSize3vsLayer = np.column_stack((LayerNumberPerSize3, LCeneOverLCsizePerSize3))
    LCeneOverLCsizePerSize4vsLayer = np.column_stack((LayerNumberPerSize4, LCeneOverLCsizePerSize4))
    LCeneOverLCsizePerSizege5vsLayer = np.column_stack((LayerNumberPerSizege5, LCeneOverLCsizePerSizege5))

    for i, obj in enumerate([LCeneOverLCsizePerSize1vsLayer,LCeneOverLCsizePerSize2vsLayer,LCeneOverLCsizePerSize3vsLayer,LCeneOverLCsizePerSize4vsLayer,LCeneOverLCsizePerSizege5vsLayer]):
        histDict_LCeneOverLCsizePerSizeVsLayer[i] = {}
        thetitle = "LC energy over LC size from LCs of size %s vs Layer" %(i+1)
        if i == 4: thetitle = "LC energy over LC size from LCs of size >= %s vs Layer" %(i+1)
        histDict_LCeneOverLCsizePerSizeVsLayer[i] = histValues2D(obj, histDict_LCeneOverLCsizePerSizeVsLayer[i], tag = "LCeneOverLCsizeForLCSize%sVsLayer"%(i+1), title = thetitle,  axunit = "Layer number", binsBoundariesX = [100, 0,100],  ayunit = "Average Cell Energy (GeV)", binsBoundariesY = [200, 0, 2], verbosityLevel=verbosityLevel)
        #histDict_LCeneOverLCsizePerSizeVsLayer[i] = profValues2D(obj, histDict_LCeneOverLCsizePerSizeVsLayer[i], tag = "LCeneOverLCsizeProfileForLCSize%sVsLayer"%(i+1), title = thetitle,  axunit = "Layer number", binsBoundariesX = [100, 0, 100],  ayunit = "Average Cell Energy (GeV)", binsBoundariesY = [200, 0, 2], verbosityLevel=verbosityLevel)
        print(histDict_LCeneOverLCsizePerSizeVsLayer[i])
        histPrintSaveAll(histDict_LCeneOverLCsizePerSizeVsLayer[i], outDir + "/PerLayer", output, tree, verbosityLevel, removeOptStat = True)

    #histsAndProps = {histDict_LCeneOverLCsizePerSizeVsLayer[0]['LCeneOverLCsizeForLCSize1']:{"leg":"LC size == 1","color":132}, histDict_LCeneOverLCsizePerSizeVsLayer[1]['LCeneOverLCsizeForLCSize2']:{"leg":"LC size == 2","color":122}, histDict_LCeneOverLCsizePerSizeVsLayer[2]['LCeneOverLCsizeForLCSize3']:{"leg":"LC size == 3","color":125}, histDict_LCeneOverLCsizePerSizeVsLayer[3]['LCeneOverLCsizeForLCSize4']:{"leg":"LC size == 4","color":131}, histDict_LCeneOverLCsizePerSizeVsLayer[4]['LCeneOverLCsizeForLCSize5']:{"leg":"LC size >= 5","color":129} }
    
    # plot these histograms on top of each other 
    #histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "LCeneOverLCsizeNormForLCVariousSizes", xaxistitle = "Average Cell Energy (GeV)", yaxistitle = "a.u./0.01 GeV ", setLogY = True, latexComment = "", funcsAndProps = None)

    #-------------------------------------------------------------------------
    #LC energy over LC size Per size: So, average cell energy per size all in one plot for each layer
    #So, total plots = #layers and in each plot #size cases
    #Again, we will start from the per LC dataframe
    histDict_LCenevsLayerPerSize = {}
        
    #LCenevsLayerPerSize = dfl.groupby(['layer','nhitAll'])[['energy']]
    #LCeplps[layer][size] -> energy column to numpy
    LCeplps = {}
    histDict_LCenevsLayerPerSize = {}
    for lay in dfl['layer'].unique():
        LCeplps[lay] = {}
        histDict_LCenevsLayerPerSize[lay] = {}
        for sz in range(0,6): #1->5
            LCeplps[lay][sz] = dfl.query("nhitAll == @sz & layer == @lay ")[['energy']].to_numpy()
            #print(lay,sz,LCeplps[lay][sz])
            histDict_LCenevsLayerPerSize[lay][sz] = {}
            histDict_LCenevsLayerPerSize[lay][sz] = histValue1D(LCeplps[lay][sz], histDict_LCenevsLayerPerSize[lay][sz], tag = "LCeneLayer%sSize%s" %(lay,sz), title = "LC energy for Layer %s and LC size %s" %(lay,sz), axunit = "LC energy (GeV)", binsBoundariesX = [200, 0, 20], ayunit = "#a.u./0.1 GeV", verbosityLevel=verbosityLevel)

        histsAndProps = {
            histDict_LCenevsLayerPerSize[lay][1]['LCeneLayer%sSize1'%(lay)]:{"leg":"LC size == 1","color":132},
            histDict_LCenevsLayerPerSize[lay][2]['LCeneLayer%sSize2'%(lay)]:{"leg":"LC size == 2","color":122},
            histDict_LCenevsLayerPerSize[lay][3]['LCeneLayer%sSize3'%(lay)]:{"leg":"LC size == 3","color":125},
            histDict_LCenevsLayerPerSize[lay][4]['LCeneLayer%sSize4'%(lay)]:{"leg":"LC size == 4","color":131},
            histDict_LCenevsLayerPerSize[lay][5]['LCeneLayer%sSize5'%(lay)]:{"leg":"LC size >= 5","color":129}
        }
    
        # plot these histograms on top of each other 
        histsPrintSaveSameCanvas(histsAndProps, outDir + "/PerLayer", tag = "LCEnergyNormLayerForLCVariousSizes_Layer%s" %(lay), xaxistitle = "LC energy (GeV)", yaxistitle = "a.u./0.01 GeV ", setLogY = True, latexComment = "", funcsAndProps = None)

    '''
    graphsAndProps = {}
    grOptions = {}

    LCenePerSize1 = df[ (df['nhitAll'] == 1)][['energy']].to_numpy().mean()
    LCenePerSize2 = df[ (df['nhitAll'] == 2)][['energy']].to_numpy().mean() 
    LCenePerSize3 = df[ (df['nhitAll'] == 3)][['energy']].to_numpy().mean() 
    LCenePerSize4 = df[ (df['nhitAll'] == 4)][['energy']].to_numpy().mean() 
    LCenePerSizege5 = df[ (df['nhitAll'] >= 5)][['energy']].to_numpy().mean()

    LayerPerSize1 = df[ (df['nhitAll'] == 1)][['rechit_layer']].to_numpy().mean()
    LayerPerSize2 = df[ (df['nhitAll'] == 2)][['rechit_layer']].to_numpy().mean() 
    LayerPerSize3 = df[ (df['nhitAll'] == 3)][['rechit_layer']].to_numpy().mean() 
    LayerPerSize4 = df[ (df['nhitAll'] == 4)][['rechit_layer']].to_numpy().mean() 
    LayerPerSizege5 = df[ (df['nhitAll'] >= 5)][['rechit_layer']].to_numpy().mean()

    graphsAndProps[ROOT.TGraph(len(LayerPerSize1), LayerPerSize1, LCenePerSize1)] = {"leg":"LC size == 1","color":132, "MarkerStyle": 21, "LineStyle": 4},
    graphsAndProps[ROOT.TGraph(len(LayerPerSize2), LayerPerSize2, LCenePerSize2)] = {"leg":"LC size == 2","color":122, "MarkerStyle": 22, "LineStyle": 5}
    graphsAndProps[ROOT.TGraph(len(LayerPerSize3), LayerPerSize3, LCenePerSize3)] = {"leg":"LC size == 3","color":125, "MarkerStyle": 23, "LineStyle": 6}
    graphsAndProps[ROOT.TGraph(len(LayerPerSize4), LayerPerSize4, LCenePerSize4)] = {"leg":"LC size == 4","color":131, "MarkerStyle": 24, "LineStyle": 7}
    graphsAndProps[ROOT.TGraph(len(LayerPerSizege5), LayerPerSizege5, LCenePerSizege5)] = {"leg":"LC size >= 5","color":129, "MarkerStyle": 25, "LineStyle": 8} 

    grOptions['Xaxis'] = "Layer number"
    grOptions['Yaxis'] = "LC Energy (GeV)"

    
    drawGraphs(graphsAndProps, grOptions, outDir, tag=type+"_vs_"+vsDep+"_"+ pidmap[pidSelected] + "_scenarios_" + "_".join(scenarios) + "_" + tag)

    drawGraphs(graphsAndProps, grOptions, outDir, latexComment=[], setLogY = True, tag="LCEneVsLayerForVariousSizes", verbosityLevel=verbosityLevel):
    '''
    #-------------------------------------------------------------------------
    #LC energy Fraction vs Layer per size all in one plot
    #We will start from the per LC dataframe
    histDict_LCeneFractionvsLayerPerSize = {}
        
    #LCeneFractionvsLayerPerSize = dfl.groupby(['layer','nhitAll'])[['energy']]
    #LCefrplps[layer][size] -> energy column to numpy
    LCeplps = {}
    LCefrplps = {}
    histDict_LCeneFractionvsLayerPerSize = {}
    for lay in dfl['layer'].unique():
        LCeplps[lay] = {}
        LCefrplps[lay] = {}
        histDict_LCeneFractionvsLayerPerSize[lay] = {}
        for sz in range(0,6): #1->5
            LCeplps[lay][sz] = dfl.query("nhitAll == @sz & layer == @lay ")[['energy']].to_numpy()
            theEgen = np.full(len(LCeplps[lay][sz].flatten()),GenEnergy) 
            LCefrplps[lay][sz] = LCeplps[lay][sz].flatten() / theEgen
            #print(lay,sz,LCefrplps[lay][sz],LCefrplps[lay][sz])
            histDict_LCeneFractionvsLayerPerSize[lay][sz] = {}
            histDict_LCeneFractionvsLayerPerSize[lay][sz] = histValue1D(LCefrplps[lay][sz], histDict_LCeneFractionvsLayerPerSize[lay][sz], tag = "LCeneFractionLayer%sSize%s" %(lay,sz), title = "LC energy Fraction for Layer %s and LC size %s" %(lay,sz), axunit = "LC energy fraction", binsBoundariesX = [40, 0, 0.4], ayunit = "#a.u.", verbosityLevel=verbosityLevel)

        histsAndProps = {
            histDict_LCeneFractionvsLayerPerSize[lay][1]['LCeneFractionLayer%sSize1'%(lay)]:{"leg":"LC size == 1","color":132},
            histDict_LCeneFractionvsLayerPerSize[lay][2]['LCeneFractionLayer%sSize2'%(lay)]:{"leg":"LC size == 2","color":122},
            histDict_LCeneFractionvsLayerPerSize[lay][3]['LCeneFractionLayer%sSize3'%(lay)]:{"leg":"LC size == 3","color":125},
            histDict_LCeneFractionvsLayerPerSize[lay][4]['LCeneFractionLayer%sSize4'%(lay)]:{"leg":"LC size == 4","color":131},
            histDict_LCeneFractionvsLayerPerSize[lay][5]['LCeneFractionLayer%sSize5'%(lay)]:{"leg":"LC size >= 5","color":129}
        }
    
        # plot these histograms on top of each other 
        histsPrintSaveSameCanvas(histsAndProps, outDir + "/PerLayer", tag = "LCEnergyFractionNormLayerForLCVariousSizes_Layer%s" %(lay), xaxistitle = "LC energy fraction", yaxistitle = "a.u.", setLogY = True, latexComment = "", funcsAndProps = None)

    #-------------------------------------------------------------------------
    histDict_LCenevsLCsize = {}
    LCenevsLCsize = ddfl[['energy','nhitAll']].to_numpy()
    histDict_LCenevsLCsize = histValues2D(LCenevsLCsize, histDict_LCenevsLCsize, tag = "LCenevsLCsize", title = "LC energy vs LC size", axunit = "LC energy (GeV)", binsBoundariesX = [20, 0, 20], ayunit = "LC size", binsBoundariesY=[50, 0., 50.], weighted2D=False, verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_LCenevsLCsize, outDir, output, tree, verbosityLevel, setLogY = True, removeOptStat = True, setLogz = True)

    #-------------------------------------------------------------------------
    #Sum of the cells energy belonging to LCs of specific sizes. Essentially, a double check of the LC.energy() method. 
    cellsFromLCEneSumPerSize = df.groupby(['EventId','nhitAll']).agg( cellsFromLCEneSum  = ('rechit_energy','sum'))
    print(cellsFromLCEneSumPerSize)
    cellEnSum_LC_Size1 = cellsFromLCEneSumPerSize.query("nhitAll == 1")[['cellsFromLCEneSum']].to_numpy().flatten()
    cellEnSum_LC_Size2 = cellsFromLCEneSumPerSize.query("nhitAll == 2")[['cellsFromLCEneSum']].to_numpy().flatten()
    cellEnSum_LC_Size3 = cellsFromLCEneSumPerSize.query("nhitAll == 3")[['cellsFromLCEneSum']].to_numpy().flatten()
    cellEnSum_LC_Size4 = cellsFromLCEneSumPerSize.query("nhitAll == 4")[['cellsFromLCEneSum']].to_numpy().flatten()
    #A little differently for the >=5 case
    ddf = df[ (df['nhitAll'] >= 5)]
    cellsFromLCEneSumPerSize = ddf.groupby(['EventId']).agg( cellsFromLCEneSum  = ('rechit_energy','sum'))
    cellEnSum_LC_Sizege5 = cellsFromLCEneSumPerSize[['cellsFromLCEneSum']].to_numpy().flatten()
    print(cellEnSum_LC_Size1, cellEnSum_LC_Size2, cellEnSum_LC_Size3, cellEnSum_LC_Size4, cellEnSum_LC_Sizege5)
    #histDict[i] = histValue1D(cellEnSum_LC_Size1, histDict, tag = "CellsEnergySumForLCSize1", title = "Cells Energy from LCs of size 1",   axunit = "Cells Energy (GeV)", binsBoundariesX = [GenEnergy, 0, GenEnergy], ayunit = "#Events/1 GeV ", verbosityLevel=verbosityLevel)

    histDict = {}

    for i, obj in enumerate([cellEnSum_LC_Size1,cellEnSum_LC_Size2,cellEnSum_LC_Size3,cellEnSum_LC_Size4,cellEnSum_LC_Sizege5]):
        histDict[i] = {}
        thetitle = "Cells Energy Sum from LCs of size %s" %(i+1)
        if i == 4: thetitle = "Cells Energy Sum from LCs of size >= %s" %(i+1)
        histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergySumForLCSize%s"%(i+1) , title = thetitle,   axunit = "Cells Energy Sum (GeV)", binsBoundariesX = [GenEnergy, 0, GenEnergy], ayunit = "#Events/1 GeV ", verbosityLevel=verbosityLevel)
        histPrintSaveAll(histDict[i], outDir, output, tree, verbosityLevel)

    histsAndProps = {histDict[0]['CellsEnergySumForLCSize1']:{"leg":"LC size == 1","color":132}, histDict[1]['CellsEnergySumForLCSize2']:{"leg":"LC size == 2","color":122}, histDict[2]['CellsEnergySumForLCSize3']:{"leg":"LC size == 3","color":125}, histDict[3]['CellsEnergySumForLCSize4']:{"leg":"LC size == 4","color":131}, histDict[4]['CellsEnergySumForLCSize5']:{"leg":"LC size >= 5","color":123} }
    
    # plot these histograms on top of each other 
    histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "CellsEnergyNormSumForLCVariousSizes", xaxistitle = "Cells Energy (GeV)", yaxistitle = "a.u./1 GeV ", setLogY = False, latexComment = "", funcsAndProps = None)

    #-------------------------------------------------------------------------
    #Sum of the cells energy belonging to LCs of specific sizes. Essentially, a double check of the LC.energy() method. 
    #Cut off at 300 MeV (0.3 GeV) for LC size 1
    df_cut = df[ (df['rechit_energy'] < 0.3) & (df['nhitAll'] == 1) ]
    print(df_cut)
    cellsFromLCEneSumPerSize = df_cut.groupby(['EventId','nhitAll']).agg( cellsFromLCEneSum  = ('rechit_energy','sum'))
    print(cellsFromLCEneSumPerSize)
    cellEnSum_LC_Size1 = cellsFromLCEneSumPerSize.query("nhitAll == 1")[['cellsFromLCEneSum']].to_numpy().flatten()
    #No cut for the rest
    cellsFromLCEneSumPerSize = df.groupby(['EventId','nhitAll']).agg( cellsFromLCEneSum  = ('rechit_energy','sum'))
    cellEnSum_LC_Size2 = cellsFromLCEneSumPerSize.query("nhitAll == 2")[['cellsFromLCEneSum']].to_numpy().flatten()
    cellEnSum_LC_Size3 = cellsFromLCEneSumPerSize.query("nhitAll == 3")[['cellsFromLCEneSum']].to_numpy().flatten()
    cellEnSum_LC_Size4 = cellsFromLCEneSumPerSize.query("nhitAll == 4")[['cellsFromLCEneSum']].to_numpy().flatten()
    #A little differently for the >=5 case
    ddf = df[ (df['nhitAll'] >= 5)]
    cellsFromLCEneSumPerSize = ddf.groupby(['EventId']).agg( cellsFromLCEneSum  = ('rechit_energy','sum'))
    cellEnSum_LC_Sizege5 = cellsFromLCEneSumPerSize[['cellsFromLCEneSum']].to_numpy().flatten()
    print(cellEnSum_LC_Size1, cellEnSum_LC_Size2, cellEnSum_LC_Size3, cellEnSum_LC_Size4, cellEnSum_LC_Sizege5)
    #histDict[i] = histValue1D(cellEnSum_LC_Size1, histDict, tag = "CellsEnergySumForLCSize1", title = "Cells Energy from LCs of size 1",   axunit = "Cells Energy (GeV)", binsBoundariesX = [GenEnergy, 0, GenEnergy], ayunit = "#Events/1 GeV ", verbosityLevel=verbosityLevel)

    histDict = {}

    for i, obj in enumerate([cellEnSum_LC_Size1,cellEnSum_LC_Size2,cellEnSum_LC_Size3,cellEnSum_LC_Size4,cellEnSum_LC_Sizege5]):
        histDict[i] = {}
        thetitle = "Cells Energy Sum from LCs of size %s" %(i+1)
        if i == 4: thetitle = "Cells Energy Sum from LCs of size >= %s" %(i+1)
        histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergySum300MeVcutForLCSize%s"%(i+1) , title = thetitle,   axunit = "Cells Energy Sum (GeV)", binsBoundariesX = [GenEnergy, 0, GenEnergy], ayunit = "#Events/1 GeV ", verbosityLevel=verbosityLevel)
        histPrintSaveAll(histDict[i], outDir, output, tree, verbosityLevel)

    histsAndProps = {histDict[0]['CellsEnergySum300MeVcutForLCSize1']:{"leg":"LC size == 1 with E_{cell} < 300 MeV","color":132}, histDict[1]['CellsEnergySum300MeVcutForLCSize2']:{"leg":"LC size == 2","color":122}, histDict[2]['CellsEnergySum300MeVcutForLCSize3']:{"leg":"LC size == 3","color":125}, histDict[3]['CellsEnergySum300MeVcutForLCSize4']:{"leg":"LC size == 4","color":131}, histDict[4]['CellsEnergySum300MeVcutForLCSize5']:{"leg":"LC size >= 5","color":123} }
    
    # plot these histograms on top of each other 
    histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "CellsEnergyNormSumForLCVariousSizes_300MeV_cut_forLC1", xaxistitle = "Cells Energy (GeV)", yaxistitle = "a.u./1 GeV ", setLogY = False, latexComment = "", funcsAndProps = None)

    #-------------------------------------------------------------------------
    #Sum of the cells energy belonging to LCs of specific sizes in MIPs. In contrary to the next approach
    #approach here uncalib energy is coming from Uncalibrated rechits which is directly in MIPs.
    #So, the results here will be different to the results below since the results here doesn't
    #contain the cce or regional em factors.
    cellsFromLCEneSumPerSize = df.groupby(['EventId','nhitAll']).agg( cellsFromLCEneSum  = ('rechit_uncalib_energy','sum'))
    print(cellsFromLCEneSumPerSize)
    cellEnSum_LC_Size1 = cellsFromLCEneSumPerSize.query("nhitAll == 1")[['cellsFromLCEneSum']].to_numpy().flatten()
    cellEnSum_LC_Size2 = cellsFromLCEneSumPerSize.query("nhitAll == 2")[['cellsFromLCEneSum']].to_numpy().flatten()
    cellEnSum_LC_Size3 = cellsFromLCEneSumPerSize.query("nhitAll == 3")[['cellsFromLCEneSum']].to_numpy().flatten()
    cellEnSum_LC_Size4 = cellsFromLCEneSumPerSize.query("nhitAll == 4")[['cellsFromLCEneSum']].to_numpy().flatten()
    #A little differently for the >=5 case
    ddf = df[ (df['nhitAll'] >= 5)]
    cellsFromLCEneSumPerSize = ddf.groupby(['EventId']).agg( cellsFromLCEneSum  = ('rechit_uncalib_energy','sum'))
    cellEnSum_LC_Sizege5 = cellsFromLCEneSumPerSize[['cellsFromLCEneSum']].to_numpy().flatten()
    print(cellEnSum_LC_Size1, cellEnSum_LC_Size2, cellEnSum_LC_Size3, cellEnSum_LC_Size4, cellEnSum_LC_Sizege5)
    #histDict[i] = histValue1D(cellEnSum_LC_Size1, histDict, tag = "CellsEnergySumForLCSize1", title = "Cells Energy from LCs of size 1",   axunit = "Cells Energy (GeV)", binsBoundariesX = [GenEnergy, 0, GenEnergy], ayunit = "#Events/1 GeV ", verbosityLevel=verbosityLevel)

    histDict = {}

    for i, obj in enumerate([cellEnSum_LC_Size1,cellEnSum_LC_Size2,cellEnSum_LC_Size3,cellEnSum_LC_Size4,cellEnSum_LC_Sizege5]):
        histDict[i] = {}
        thetitle = "Cells Energy Sum from LCs of size %s" %(i+1)
        if i == 4:
            thetitle = "Cells Energy Sum in MIPs from LCs of size >= %s" %(i+1)
            histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergySumInMIPsFromUncalibForLCSize%s"%(i+1) , title = thetitle,   axunit = "Cells Energy Sum (MIPs)", binsBoundariesX = [200, 0, 10000], ayunit = "#Events/50 MIP", verbosityLevel=verbosityLevel)
        else: 
            thetitle = "Cells Energy Sum in MIPs from LCs of size %s" %(i+1)
            histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergySumInMIPsFromUncalibForLCSize%s"%(i+1) , title = thetitle,   axunit = "Cells Energy Sum (MIPs)", binsBoundariesX = [100, 0, 1000], ayunit = "#Events/10 MIP", verbosityLevel=verbosityLevel)
        histPrintSaveAll(histDict[i], outDir, output, tree, verbosityLevel)

    #histsAndProps = {histDict[0]['CellsEnergySumInMIPsFromUncalibForLCSize1']:{"leg":"LC size == 1","color":132}, histDict[1]['CellsEnergySumInMIPsFromUncalibForLCSize2']:{"leg":"LC size == 2","color":122}, histDict[2]['CellsEnergySumInMIPsFromUncalibForLCSize3']:{"leg":"LC size == 3","color":125}, histDict[3]['CellsEnergySumInMIPsFromUncalibForLCSize4']:{"leg":"LC size == 4","color":131}, histDict[4]['CellsEnergySumInMIPsFromUncalibForLCSize5']:{"leg":"LC size >= 5","color":123} }
    
    # plot these histograms on top of each other 
    #histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "CellsEnergyNormSumInMIPsFromUncalibForLCVariousSizes", xaxistitle = "Cells Energy (GeV)", yaxistitle = "a.u./1 GeV ", setLogY = False, latexComment = "", funcsAndProps = None)

    
    #-------------------------------------------------------------------------
    #Sum of the cells energy belonging to LCs of specific sizes in MIPs. So, this
    #cannot be taken from LC.energy(). We should undo the calibration of each cell and then sum
    #The uncalibrated energy will be rechit_energy (GeV) * (1/dedxWeights*0.001) (MIP/GeV), so
    #regional e/m factors and cce corrections are still applied.
    #Observe the manipulation below to get dEdxWeights with the convention 0 to 99. 
    dEdxWeights_allLayers = np.append( np.array(dEdX.weights[1:]), np.array(dEdX.weights[1:]) ) 
    #Layer clusters with at least one hit
    ddf = df.loc[ (df['nhitAll'] >= 1) ].copy()
    ddf["rechit_uncalib_energy"] = ddf.apply(lambda row: row["rechit_energy"] / (dEdxWeights_allLayers[int(row["rechit_layer"])] * 0.001),axis=1)

    cellsFromUncalibLCEneSumPerSize = ddf.groupby(['EventId','nhitAll']).agg( cellsFromUncalibLCEneSum  = ('rechit_uncalib_energy','sum'))
    print(cellsFromUncalibLCEneSumPerSize)
    cellUncalibEnSum_LC_Size1 = cellsFromUncalibLCEneSumPerSize.query("nhitAll == 1")[['cellsFromUncalibLCEneSum']].to_numpy().flatten()
    cellUncalibEnSum_LC_Size2 = cellsFromUncalibLCEneSumPerSize.query("nhitAll == 2")[['cellsFromUncalibLCEneSum']].to_numpy().flatten()
    cellUncalibEnSum_LC_Size3 = cellsFromUncalibLCEneSumPerSize.query("nhitAll == 3")[['cellsFromUncalibLCEneSum']].to_numpy().flatten()
    cellUncalibEnSum_LC_Size4 = cellsFromUncalibLCEneSumPerSize.query("nhitAll == 4")[['cellsFromUncalibLCEneSum']].to_numpy().flatten()
    #A little differently for the >=5 case
    dddf = df.loc[ (df['nhitAll'] >= 5) ].copy()
    dddf["rechit_uncalib_energy"] = dddf.apply(lambda row: row["rechit_energy"] / (dEdxWeights_allLayers[int(row["rechit_layer"])] * 0.001),axis=1)
    cellsFromUncalibLCEneSumPerSize = dddf.groupby(['EventId']).agg( cellsFromUncalibLCEneSum  = ('rechit_uncalib_energy','sum'))
    cellUncalibEnSum_LC_Sizege5 = cellsFromUncalibLCEneSumPerSize[['cellsFromUncalibLCEneSum']].to_numpy().flatten()
    print(cellUncalibEnSum_LC_Size1, cellUncalibEnSum_LC_Size2, cellUncalibEnSum_LC_Size3, cellUncalibEnSum_LC_Size4, cellUncalibEnSum_LC_Sizege5)
    #histDict[i] = histValue1D(cellEnSum_LC_Size1, histDict, tag = "CellsEnergySumForLCSize1", title = "Cells Energy from LCs of size 1",   axunit = "Cells Energy (GeV)", binsBoundariesX = [GenEnergy, 0, GenEnergy], ayunit = "#Events/1 GeV ", verbosityLevel=verbosityLevel)

    histDict = {}

    for i, obj in enumerate([cellUncalibEnSum_LC_Size1,cellUncalibEnSum_LC_Size2,cellUncalibEnSum_LC_Size3,cellUncalibEnSum_LC_Size4,cellUncalibEnSum_LC_Sizege5]):
        histDict[i] = {}
        if i == 4:
            thetitle = "Cells Energy Sum in MIPs from LCs of size >= %s" %(i+1)
            histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergySumInMIPsForLCSize%s"%(i+1) , title = thetitle,   axunit = "Cells Energy Sum (MIPs)", binsBoundariesX = [200, 0, 10000], ayunit = "#Events/50 MIP", verbosityLevel=verbosityLevel)
        else: 
            thetitle = "Cells Energy Sum in MIPs from LCs of size %s" %(i+1)
            histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergySumInMIPsForLCSize%s"%(i+1) , title = thetitle,   axunit = "Cells Energy Sum (MIPs)", binsBoundariesX = [100, 0, 1000], ayunit = "#Events/10 MIP", verbosityLevel=verbosityLevel)
        histPrintSaveAll(histDict[i], outDir, output, tree, verbosityLevel)
   

    #-------------------------------------------------------------------------
    #Max energy belonging to LCs of specific sizes.
    cellsFromLCEneMaxPerSize = df.groupby(['EventId','nhitAll']).agg( cellsFromLCEneMax  = ('rechit_energy','max'))
    print(cellsFromLCEneMaxPerSize)
    cellEnMax_LC_Size1 = cellsFromLCEneMaxPerSize.query("nhitAll == 1")[['cellsFromLCEneMax']].to_numpy().flatten()
    cellEnMax_LC_Size2 = cellsFromLCEneMaxPerSize.query("nhitAll == 2")[['cellsFromLCEneMax']].to_numpy().flatten()
    cellEnMax_LC_Size3 = cellsFromLCEneMaxPerSize.query("nhitAll == 3")[['cellsFromLCEneMax']].to_numpy().flatten()
    cellEnMax_LC_Size4 = cellsFromLCEneMaxPerSize.query("nhitAll == 4")[['cellsFromLCEneMax']].to_numpy().flatten()
    #A little differently for the >=5 case
    ddf = df[ (df['nhitAll'] >= 5)]
    cellsFromLCEneMaxPerSize = ddf.groupby(['EventId']).agg( cellsFromLCEneMax  = ('rechit_energy','max'))
    cellEnMax_LC_Sizege5 = cellsFromLCEneMaxPerSize[['cellsFromLCEneMax']].to_numpy().flatten()
    print(cellEnMax_LC_Size1, cellEnMax_LC_Size2, cellEnMax_LC_Size3, cellEnMax_LC_Size4, cellEnMax_LC_Sizege5)

    for i, obj in enumerate([cellEnMax_LC_Size1,cellEnMax_LC_Size2,cellEnMax_LC_Size3,cellEnMax_LC_Size4,cellEnMax_LC_Sizege5]):
        histDict[i] = {}
        thetitle = "Cells Energy Max from LCs of size %s" %(i+1)
        if i == 4: thetitle = "Cells Energy Max from LCs of size >= %s" %(i+1)
        histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergyMaxForLCSize%s"%(i+1) , title = thetitle,   axunit = "Cells Energy Max (GeV)", binsBoundariesX = [20, 0, 20], ayunit = "#Events/1 GeV ", verbosityLevel=verbosityLevel)
        histPrintSaveAll(histDict[i], outDir, output, tree, verbosityLevel)

    #histsAndProps = {histDict[0]['CellsEnergyMaxForLCSize1']:{"leg":"LC size == 1","color":132}, histDict[1]['CellsEnergyMaxForLCSize2']:{"leg":"LC size == 2","color":122}, histDict[2]['CellsEnergyMaxForLCSize3']:{"leg":"LC size == 3","color":125}, histDict[3]['CellsEnergyMaxForLCSize4']:{"leg":"LC size == 4","color":131}, histDict[4]['CellsEnergyMaxForLCSize5']:{"leg":"LC size >= 5","color":129} }
    
    # plot these histograms on top of each other 
    #histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "CellsEnergyMaxNormForLCVariousSizes", xaxistitle = "Cells Energy (GeV)", yaxistitle = "a.u./0.01 GeV ", latexComment = "", funcsAndProps = None)

    #-------------------------------------------------------------------------
    #Cells energy belonging to LCs of specific sizes.
    cellsFromLCEnePerSize1 = df[ (df['nhitAll'] == 1)][['rechit_energy']].to_numpy().flatten()
    cellsFromLCEnePerSize2 = df[ (df['nhitAll'] == 2)][['rechit_energy']].to_numpy().flatten()
    cellsFromLCEnePerSize3 = df[ (df['nhitAll'] == 3)][['rechit_energy']].to_numpy().flatten()
    cellsFromLCEnePerSize4 = df[ (df['nhitAll'] == 4)][['rechit_energy']].to_numpy().flatten()
    cellsFromLCEnePerSizege5 = df[ (df['nhitAll'] >= 5)][['rechit_energy']].to_numpy().flatten()
    print(cellsFromLCEnePerSize1,cellsFromLCEnePerSize2,cellsFromLCEnePerSize3,cellsFromLCEnePerSize4,cellsFromLCEnePerSizege5)

    for i, obj in enumerate([cellsFromLCEnePerSize1,cellsFromLCEnePerSize2,cellsFromLCEnePerSize3,cellsFromLCEnePerSize4,cellsFromLCEnePerSizege5]):
        histDict[i] = {}
        thetitle = "Cells Energy from LCs of size %s" %(i+1)
        if i == 4: thetitle = "Cells Energy from LCs of size >= %s" %(i+1)
        histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergyForLCSize%s"%(i+1), title = thetitle,   axunit = "Cells Energy (GeV)", binsBoundariesX = [200, 0, 2], ayunit = "#Events/0.01 GeV ", verbosityLevel=verbosityLevel)
        histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergyForLCSize%s_Zoomed"%(i+1), title = thetitle,   axunit = "Cells Energy (GeV)", binsBoundariesX = [20, 0, 0.2], ayunit = "#Events/0.01 GeV ", verbosityLevel=verbosityLevel)
        print(histDict[i])
        histPrintSaveAll(histDict[i], outDir, output, tree, verbosityLevel, setLogY = True)

    histsAndProps = {histDict[0]['CellsEnergyForLCSize1']:{"leg":"LC size == 1","color":132}, histDict[1]['CellsEnergyForLCSize2']:{"leg":"LC size == 2","color":122}, histDict[2]['CellsEnergyForLCSize3']:{"leg":"LC size == 3","color":125}, histDict[3]['CellsEnergyForLCSize4']:{"leg":"LC size == 4","color":131}, histDict[4]['CellsEnergyForLCSize5']:{"leg":"LC size >= 5","color":123}, histDict_LCeneOverLCsize["LCeneOverLCsize"]:{"leg":"Average cell energy","color":129} }
    
    # plot these histograms on top of each other 
    histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "CellsEnergyNormForLCVariousSizes", xaxistitle = "Cells Energy (GeV)", yaxistitle = "a.u./0.01 GeV ", setLogY = True, latexComment = "", funcsAndProps = None)

    # make a version with logX
    histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "CellsEnergyNormForLCVariousSizes_LogX", xaxistitle = "Cells Energy (GeV)", yaxistitle = "a.u./0.01 GeV ", setLogX = True, setLogY = True, latexComment = "", funcsAndProps = None)

    # make a zoomed version
    histsAndProps = {histDict[0]['CellsEnergyForLCSize1_Zoomed']:{"leg":"LC size == 1","color":132}, histDict[1]['CellsEnergyForLCSize2_Zoomed']:{"leg":"LC size == 2","color":122}, histDict[2]['CellsEnergyForLCSize3_Zoomed']:{"leg":"LC size == 3","color":125}, histDict[3]['CellsEnergyForLCSize4_Zoomed']:{"leg":"LC size == 4","color":131}, histDict[4]['CellsEnergyForLCSize5_Zoomed']:{"leg":"LC size >= 5","color":123}, histDict_LCeneOverLCsize["LCeneOverLCsize_Zoomed"]:{"leg":"Average cell energy","color":129} }
    histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "CellsEnergyNormForLCVariousSizes_Zoomed", xaxistitle = "Cells Energy (GeV)", yaxistitle = "a.u./0.01 GeV ", latexComment = "", funcsAndProps = None)

    #-------------------------------------------------------------------------
    #Cells energy in MIPs belonging to LCs of specific sizes.
    df_mips = df[ (df['id'] >= 0)]
    cellsFromLCEnePerSize1 = df_mips[ (df_mips['nhitAll'] == 1)][['rechit_uncalib_energy']].to_numpy().flatten()
    cellsFromLCEnePerSize2 = df_mips[ (df_mips['nhitAll'] == 2)][['rechit_uncalib_energy']].to_numpy().flatten()
    cellsFromLCEnePerSize3 = df_mips[ (df_mips['nhitAll'] == 3)][['rechit_uncalib_energy']].to_numpy().flatten()
    cellsFromLCEnePerSize4 = df_mips[ (df_mips['nhitAll'] == 4)][['rechit_uncalib_energy']].to_numpy().flatten()
    cellsFromLCEnePerSizege5 = df_mips[ (df_mips['nhitAll'] >= 5)][['rechit_uncalib_energy']].to_numpy().flatten()
    #Now, for the average cell energy in MIPs, we need LC energy in MIPs
    LCEneInMIPs = df_mips.groupby(['EventId','id','nhitAll']).agg( LCEneInMIPsFromSum  = ('rechit_uncalib_energy','sum'))
    print(LCEneInMIPs)
    print(LCEneInMIPs.index.get_level_values('nhitAll'))
    LCeneOverLCsizeinMIPs = LCEneInMIPs[['LCEneInMIPsFromSum']].to_numpy().flatten() / LCEneInMIPs.index.get_level_values('nhitAll').to_numpy().flatten()

    histDict_LCeneOverLCsizeinMIPs = {}
    histDict_LCeneOverLCsizeinMIPs = histValue1D(LCeneOverLCsizeinMIPs, histDict_LCeneOverLCsizeinMIPs, tag = "LCeneOverLCsizeinMIPs", title = "LC energy in MIPs over LC size",   axunit = "Average Cell Energy (MIPs)", binsBoundariesX = [100, 0, 100], ayunit = "#Events/1 MIPs", verbosityLevel=verbosityLevel)

    #print(cellsFromLCEnePerSize1,cellsFromLCEnePerSize2,cellsFromLCEnePerSize3,cellsFromLCEnePerSize4,cellsFromLCEnePerSizege5)

    for i, obj in enumerate([cellsFromLCEnePerSize1,cellsFromLCEnePerSize2,cellsFromLCEnePerSize3,cellsFromLCEnePerSize4,cellsFromLCEnePerSizege5]):
        histDict[i] = {}
        thetitle = "Cells Energy in MIPs from LCs of size %s" %(i+1)
        if i == 4: thetitle = "Cells Energy in MIPs from LCs of size >= %s" %(i+1)
        histDict[i] = histValue1D(obj, histDict[i], tag = "CellsEnergyinMIPsForLCSize%s"%(i+1), title = thetitle,   axunit = "Cells Energy (MIPs)", binsBoundariesX = [100, 0, 100], ayunit = "#Events/1 MIP ", verbosityLevel=verbosityLevel)
        print(histDict[i])
        histPrintSaveAll(histDict[i], outDir, output, tree, verbosityLevel, setLogY = True)

    histsAndProps = {histDict[0]['CellsEnergyinMIPsForLCSize1']:{"leg":"LC size == 1","color":132}, histDict[1]['CellsEnergyinMIPsForLCSize2']:{"leg":"LC size == 2","color":122}, histDict[2]['CellsEnergyinMIPsForLCSize3']:{"leg":"LC size == 3","color":125}, histDict[3]['CellsEnergyinMIPsForLCSize4']:{"leg":"LC size == 4","color":131}, histDict[4]['CellsEnergyinMIPsForLCSize5']:{"leg":"LC size >= 5","color":123}, histDict_LCeneOverLCsizeinMIPs["LCeneOverLCsizeinMIPs"]:{"leg":"Average cell energy","color":129} }
    
    # plot these histograms on top of each other 
    histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "CellsEnergyinMIPsNormForLCVariousSizes", xaxistitle = "Cells Energy (MIPs)", yaxistitle = "a.u./1 MIP ", setLogY = True, latexComment = "", funcsAndProps = None)

    # make a version with logX
    histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "CellsEnergyinMIPsNormForLCVariousSizes_LogX", xaxistitle = "Cells Energy (MIPs)", yaxistitle = "a.u./1 MIP ", setLogX = True, setLogY = True, latexComment = "", funcsAndProps = None)

    #histPrintSaveAll(histDict, outDir, output, tree, verbosityLevel)

def maketree(theinarray, addname, tree):

    outbranch = np.empty((1), dtype="float32")
    tree.Branch(addname, outbranch, "%s/F"%(addname))
    
    print(addname)
    for i in theinarray:
        print(i)
        outbranch[0] = i
        tree.Fill()

#---------------------------------------------------------------------------------------------------
def recHitCalibrationPlots(df_m,df_u, tree,maxEvents,outDir,output,GenEnergy,ecut,verbosityLevel = 0):

    histDict = {}
    #For the fit that was previously produced the regional em factors.
    #kept only for crosscheck and historical reasons.
    eratioboundaries = {}
    eratioboundaries["CE_E_Front_120um"] = [0.6 , 1.1]##[0.85 , 1.15]#[0.6 , 0.9][0.85 , 1.15]
    eratioboundaries["CE_E_Front_200um"] = [0.6 , 1.1]#[0.85 , 1.15]#[0.6 , 0.9][0.85 , 1.15]
    eratioboundaries["CE_E_Front_300um"] = [0.6 , 1.1]#[0.8 , 1.2]#[0.6 , 1.1][0.75 , 1.45]
    eratioboundaries["CE_H_Fine_120um"] = [0.5 , 1.2]
    eratioboundaries["CE_H_Fine_200um"] = [0.5 , 1.2]
    eratioboundaries["CE_H_Fine_300um"] = [0.5 , 1.2]
    eratioboundaries["CE_H_Fine_300um_Var1"] = [0.5 , 1.2]
    eratioboundaries["CE_H_Coarse_Scint"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_300um"] = [0.6 , 1.2] #[0.8 , 1.2]#[0.6 , 1.1][0.75 , 1.45]
    eratioboundaries["CE_H_Fine_Scint"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Fine_Scint_Var1"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Fine_Scint_Var2"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_Var1"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_Var2"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_4285"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_4295"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_4305"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_4315"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_4325"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_4335"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_4345"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_4354"] = [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]
    eratioboundaries["CE_H_Coarse_Scint_4364"] =  [0.6 , 1.2]#[0.45 , 0.95]#[0.85 , 1.15]#[0.6 , 1.0]

    EtaBoundariesShower = {}
    EtaBoundariesShower["CE_E_Front_120um"] = [400, 2.3, 2.7]
    EtaBoundariesShower["CE_E_Front_200um"] = [400, 1.8, 2.2]
    EtaBoundariesShower["CE_E_Front_300um"] = [200, 1.5, 1.7]
    EtaBoundariesShower["CE_H_Fine_120um"] = [400, 2.3, 2.7]
    EtaBoundariesShower["CE_H_Fine_200um"] = [400, 1.8, 2.2]
    EtaBoundariesShower["CE_H_Fine_300um"] = [200, 1.5, 1.7]
    EtaBoundariesShower["CE_H_Fine_300um_Var1"] = [200, 1.5, 1.7]
    EtaBoundariesShower["CE_H_Coarse_Scint"] = [400, 1.5, 3.5]
    EtaBoundariesShower["CE_H_Coarse_300um"] = [200, 1.5, 1.7]
    EtaBoundariesShower["CE_H_Fine_Scint"] = [400, 1.5, 3.5]

    # Matched rechits plots
    rHxEnFrSum_m = {}

    #recHitEneXFractionOverEgenSumPerThick = df_m.groupby(['EventId','sClusHitsThick','sClusHitsDet']).agg( recHitEneXFractionOverEgenSum  = ('recHitEneXFractionOvertheEgen','sum'))
    #The object recHitEneXFractionOverEgenSumPerThick has a multiindex so let's save the quantities we want
    #The conventions are: Tracker = 1, Muon = 2,Ecal = 3,Hcal = 4,Calo = 5,Forward = 6,VeryForward = 7,HGCalEE = 8,HGCalHSi = 9,HGCalHSc = 10,HGCalTrigger = 11 
    rHxEnFrSum_m["CE_E_Front_120um"] = df_m[ (df_m["sClusHitsThick"] == 120) & (df_m["sClusHitsDet"] == 8)].groupby(['EventId']).agg( recHitEneXFractionOverEgenSum  = ('recHitEneXFractionOvertheEgen','sum'))[["recHitEneXFractionOverEgenSum"]].to_numpy().flatten()

    rHxEnFrSum_m["CE_E_Front_200um"] = df_m[ (df_m["sClusHitsThick"] == 200) & (df_m["sClusHitsDet"] == 8)].groupby(['EventId']).agg( recHitEneXFractionOverEgenSum  = ('recHitEneXFractionOvertheEgen','sum'))[['recHitEneXFractionOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_m["CE_E_Front_300um"] = df_m[ (df_m["sClusHitsThick"] == 300) & (df_m["sClusHitsDet"] == 8)].groupby(['EventId']).agg( recHitEneXFractionOverEgenSum  = ('recHitEneXFractionOvertheEgen','sum'))[['recHitEneXFractionOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_m["CE_H_Fine_120um"]  = df_m[ (df_m["sClusHitsThick"] == 120) & (df_m["sClusHitsDet"] == 9)].groupby(['EventId']).agg( recHitEneXFractionOverEgenSum  = ('recHitEneXFractionOvertheEgen','sum'))[['recHitEneXFractionOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_m["CE_H_Fine_200um"]  = df_m[ (df_m["sClusHitsThick"] == 200) & (df_m["sClusHitsDet"] == 9)].groupby(['EventId']).agg( recHitEneXFractionOverEgenSum  = ('recHitEneXFractionOvertheEgen','sum'))[['recHitEneXFractionOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_m["CE_H_Fine_300um"]  = df_m[ (df_m["sClusHitsThick"] == 300) & (df_m["sClusHitsDet"] == 9)].groupby(['EventId']).agg( recHitEneXFractionOverEgenSum  = ('recHitEneXFractionOvertheEgen','sum'))[['recHitEneXFractionOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_m["CE_H_Coarse_Scint"] = df_m[ (df_m["sClusHitsThick"] > 400) & (df_m["sClusHitsDet"] == 10)].groupby(['EventId']).agg( recHitEneXFractionOverEgenSum  = ('recHitEneXFractionOvertheEgen','sum'))[['recHitEneXFractionOverEgenSum']].to_numpy().flatten()
    
    #Counting hits of specific thickness
    rHitthick120 = len( df_m[ (df_m['sClusHitsThick'] == 120) ].to_numpy())
    rHitthick200 = len( df_m[ (df_m['sClusHitsThick'] == 200) ].to_numpy())
    rHitthick300 = len( df_m[ (df_m['sClusHitsThick'] == 300) ].to_numpy())
    rHitthickScint = len( df_m[ (df_m['sClusHitsThick'] > 400) ].to_numpy())
    rHitthickallmatched = df_m.shape[0]

    fittree_m = {}
    
    for i, obj in rHxEnFrSum_m.items():
        print(i,obj, "SumEoverEgen_%s"%(i))
        #print(i,len(obj))
        fittree_m[i] = ROOT.TTree("ttree_SumEoverEgen_e%s_%s" %(GenEnergy, i), "results")
        maketree(obj, "SumEoverEgen_%s"%(i), fittree_m[i])
        histDict[i] = {}
        histDict[i] = histValue1D(obj, histDict[i], tag = "SumEoverEgen_%s" %(i), title = "Reconstructed hits energy over generated energy for %s" %(i),   axunit = "#sum E_{i}/E_{gen}",    binsBoundariesX = [400, 0, 2], ayunit = "N(events)", verbosityLevel=verbosityLevel)
        histPrintSaveAll(histDict[i], outDir, output, tree, verbosityLevel)
        #A version of this plot with the fit. 
        ROOT.gStyle.SetOptStat(0);
        mycE = ROOT.TCanvas("P22E%s_thick%s"%(GenEnergy,i), "P22E%s_thick%s"%(GenEnergy,i), 500, 500)
        histDict[i]["SumEoverEgen_%s" %(i)].Draw("")
        histDict[i]["SumEoverEgen_%s" %(i)].GetXaxis().SetRangeUser(eratioboundaries[i][0],eratioboundaries[i][1])
        mycE.Update()
 
        meanE = histDict[i]["SumEoverEgen_%s" %(i)].GetMean()
        rmsE  = histDict[i]["SumEoverEgen_%s" %(i)].GetRMS()

        #histDict[i]["SumEoverEgen_%s" %(i)].Fit("gaus","LR0","", meanE-0.5*rmsE, meanE+1.5*rmsE);
        #histDict[i]["SumEoverEgen_%s" %(i)].Fit("gaus","LR0","");
        histDict[i]["SumEoverEgen_%s" %(i)].Fit("gaus","LR0","",eratioboundaries[i][0],eratioboundaries[i][1])
        fitResult = histDict[i]["SumEoverEgen_%s" %(i)].GetFunction("gaus")

        fitResult.SetLineColor(2)
        meanFitE = fitResult.GetParameter(1)
        rmsFitE = fitResult.GetParameter(2)
        meanFitEerr = fitResult.GetParError(1)
        rmsFitEerr = fitResult.GetParError(2)

        fitResult.Draw("same")
        #latx =  histDict[i]["SumEoverEgen_%s" %(i)].GetXaxis().GetXmin()+(histDict[i]["SumEoverEgen_%s" %(i)].GetXaxis().GetXmax()-histDict[i]["SumEoverEgen_%s" %(i)].GetXaxis().GetXmin())/20.
        #laty =  histDict[i]["SumEoverEgen_%s" %(i)].GetMaximum()
        latx = eratioboundaries[i][0] + (eratioboundaries[i][1]-eratioboundaries[i][0])/20.
        laty =  histDict[i]["SumEoverEgen_%s" %(i)].GetMaximum()

        lat = ROOT.TLatex()
        lat.DrawLatex(latx,laty*0.9, "<#sum E_{i} * frac/E_{gen}> = %3.3f +/- %3.3f"%(fitResult.GetParameter(1),fitResult.GetParError(1))    )
        lat.DrawLatex(latx,laty*0.8, "RMSfit = %3.3f +/- %3.3f"%(fitResult.GetParameter(2),fitResult.GetParError(2))   )
        lat.DrawLatex(latx,laty*0.7, "RMS/meanfit = %3.3f"%(fitResult.GetParameter(2)/fitResult.GetParameter(1))   )
        lat.DrawLatex(latx,laty*0.6, "#chi^{2}/N = %3.3f/%d = %3.3f"%(fitResult.GetChisquare(),fitResult.GetNDF(),fitResult.GetChisquare()/fitResult.GetNDF()) )
        #lat.DrawLatex(latx,laty*0.5, "S/N = %d "% options.ecut    )
        lat.DrawLatex(latx,laty*0.5, "# hits %d #mum = %3.3f %%"%(120,(rHitthick120 * 100.) / rHitthickallmatched ) )
        lat.DrawLatex(latx,laty*0.4, "# hits %d #mum = %3.3f %%"%(200,(rHitthick200 * 100.) / rHitthickallmatched ) )
        lat.DrawLatex(latx,laty*0.3, "# hits %d #mum = %3.3f %%"%(300,(rHitthick300 * 100.) / rHitthickallmatched ) )
        lat.DrawLatex(latx,laty*0.2, "# hits Scint #mum = %3.3f %%"%((rHitthickScint * 100.) / rHitthickallmatched ) )
        mycE.Update()
        mycE.SaveAs("%s/P22E%s_thick%s.png"%(outDir,GenEnergy,i))
 
    
    # Unmatched rechits plots
    rHxEnFrSum_u = {}

    rHxEnFrSum_u["CE_E_Front_120um"] = df_u[ (df_u["sClusHitsThick"] == 120) & (df_u["sClusHitsDet"] == 8)].groupby(['EventId']).agg( recHitEneOverEgenSum  = ('recHitEneOvertheEgen','sum'))[["recHitEneOverEgenSum"]].to_numpy().flatten()

    rHxEnFrSum_u["CE_E_Front_200um"] = df_u[ (df_u["sClusHitsThick"] == 200) & (df_u["sClusHitsDet"] == 8)].groupby(['EventId']).agg( recHitEneOverEgenSum  = ('recHitEneOvertheEgen','sum'))[['recHitEneOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_u["CE_E_Front_300um"] = df_u[ (df_u["sClusHitsThick"] == 300) & (df_u["sClusHitsDet"] == 8)].groupby(['EventId']).agg( recHitEneOverEgenSum  = ('recHitEneOvertheEgen','sum'))[['recHitEneOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_u["CE_H_Fine_120um"]  = df_u[ (df_u["sClusHitsThick"] == 120) & (df_u["sClusHitsDet"] == 9)].groupby(['EventId']).agg( recHitEneOverEgenSum  = ('recHitEneOvertheEgen','sum'))[['recHitEneOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_u["CE_H_Fine_200um"]  = df_u[ (df_u["sClusHitsThick"] == 200) & (df_u["sClusHitsDet"] == 9)].groupby(['EventId']).agg( recHitEneOverEgenSum  = ('recHitEneOvertheEgen','sum'))[['recHitEneOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_u["CE_H_Fine_300um"]  = df_u[ (df_u["sClusHitsThick"] == 300) & (df_u["sClusHitsDet"] == 9)].groupby(['EventId']).agg( recHitEneOverEgenSum  = ('recHitEneOvertheEgen','sum'))[['recHitEneOverEgenSum']].to_numpy().flatten()

    rHxEnFrSum_u["CE_H_Coarse_Scint"] = df_u[ (df_u["sClusHitsThick"] > 400) & (df_u["sClusHitsDet"] == 10)].groupby(['EventId']).agg( recHitEneOverEgenSum  = ('recHitEneOvertheEgen','sum'))[['recHitEneOverEgenSum']].to_numpy().flatten()

    fittree_u = {}

    for i, obj in rHxEnFrSum_u.items():
        print(i,obj)
        fittree_u[i] = ROOT.TTree("ttree_SumEoverEgen_Unmatched_e%s_%s" %(GenEnergy, i), "results")
        maketree(obj, "SumEoverEgen_Unmatched_%s"%(i), fittree_u[i])
        #maketree(obj, "SumEoverEgen_Unmatched_%s"%(i), tree)
        histDict[i] = {}
        histDict[i] = histValue1D(obj, histDict[i], tag = "SumEoverEgenUnmatched_%s" %(i), title = "Unmatched reconstructed hits energy over generated energy for %s" %(i),   axunit = "#sum E_{i}/E_{gen}",    binsBoundariesX = [400, 0, 2], ayunit = "N(events)", verbosityLevel=verbosityLevel)
        histPrintSaveAll(histDict[i], outDir, output, tree, verbosityLevel)

    #----------------------------------------------------------
    #RvsEtavsThickness
    RechitRvsEtavsThickness = {}
    RechitRvsEtavsThickness["CE_E_Front_120um"] = df_m.query("sClusHitsThick == 120 & sClusHitsDet == 8")[['rechit_eta', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsEtavsThickness["CE_E_Front_200um"] = df_m.query("sClusHitsThick == 200 & sClusHitsDet == 8")[['rechit_eta', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsEtavsThickness["CE_E_Front_300um"] = df_m.query("sClusHitsThick == 300 & sClusHitsDet == 8")[['rechit_eta', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsEtavsThickness["CE_H_Fine_120um"] = df_m.query("sClusHitsThick == 120 & sClusHitsDet == 9")[['rechit_eta', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsEtavsThickness["CE_H_Fine_200um"] = df_m.query("sClusHitsThick == 200 & sClusHitsDet == 9")[['rechit_eta', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsEtavsThickness["CE_H_Fine_300um"] = df_m.query("sClusHitsThick == 300 & sClusHitsDet == 9")[['rechit_eta', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsEtavsThickness["CE_H_Coarse_Scint"] = df_m.query("sClusHitsThick > 400 & sClusHitsDet == 10")[['rechit_eta', 'R', 'sClusHitsThick']].to_numpy()

    histDict = {}

    for i, obj in RechitRvsEtavsThickness.items():
        print(i,obj)
        histDict[i] = {}
        histDict[i] = profValues2D(obj, histDict[i], tag = "RvsEtavsThickness_%s" %(i), title = "R vs Eta vs Thickness for %s" %(i) , axunit = "|#eta|", binsBoundariesX = [200, 1.5, 3.5], ayunit = "R (cm)", binsBoundariesY=[100, 0., 300.], binsBoundariesZ=[0.,400.], weighted2D=False, verbosityLevel=verbosityLevel)

        ROOT.gStyle.SetOptStat(0)
        #ROOT.gStyle.SetOptStat(1101)

        mycE1 = ROOT.TCanvas("P22E%s_RvsEtavsThickness_thick%s"%(GenEnergy,i), "P22E%s_RvsEtavsThickness_thick%s"%(GenEnergy,i), 500, 500)
        acustompalette()
        ex1 = ROOT.TExec("ex1","acustompalette();");
        ex1.Draw();

        histDict[i]["RvsEtavsThickness_%s"%(i)].Draw("COLZ")
        #RvsEtavsThickness.GetXaxis().SetTitleOffset(0.95) 
        mycE1.Update()

        palette = histDict[i]["RvsEtavsThickness_%s"%(i)].GetListOfFunctions().FindObject("palette")
        if palette:
            palette.__class__ = ROOT.TPaletteAxis
            palette.SetX1NDC(0.85)
            palette.SetX2NDC(0.9)
            #palette.SetY1NDC(0.1)
            #palette.SetY2NDC(0.6)
            palette.GetAxis().SetTickSize(.01)
            palette.GetAxis().SetTitle("Si thick")
            palette.GetAxis().SetTitleOffset(0.82);
            #palette.GetAxis().LabelsOption("v")
            ROOT.gPad.Update()

        mycE1.SaveAs("%s/P22E%s_RvsEtavsThickness_thick_%s.png"%(outDir,GenEnergy,i))
        #mycE1.Write()

    #histPrintSaveAll(histDict, outDir, output, tree, verbosityLevel)

    #----------------------------------------------------------
    # RechitRvsLayervsThickness
    RechitRvsLayervsThickness = {}
    RechitRvsLayervsThickness["CE_E_Front_120um"] = df_m.query("sClusHitsThick == 120 & sClusHitsDet == 8")[['sClusHitsLayers', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsLayervsThickness["CE_E_Front_200um"] = df_m.query("sClusHitsThick == 200 & sClusHitsDet == 8")[['sClusHitsLayers', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsLayervsThickness["CE_E_Front_300um"] = df_m.query("sClusHitsThick == 300 & sClusHitsDet == 8")[['sClusHitsLayers', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsLayervsThickness["CE_H_Fine_120um"] = df_m.query("sClusHitsThick == 120 & sClusHitsDet == 9")[['sClusHitsLayers', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsLayervsThickness["CE_H_Fine_200um"] = df_m.query("sClusHitsThick == 200 & sClusHitsDet == 9")[['sClusHitsLayers', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsLayervsThickness["CE_H_Fine_300um"] = df_m.query("sClusHitsThick == 300 & sClusHitsDet == 9")[['sClusHitsLayers', 'R', 'sClusHitsThick']].to_numpy()
    RechitRvsLayervsThickness["CE_H_Coarse_Scint"] = df_m.query("sClusHitsThick > 400 & sClusHitsDet == 10")[['sClusHitsLayers', 'R', 'sClusHitsThick']].to_numpy()

    histDict = {}

    for i, obj in RechitRvsLayervsThickness.items():
        print(i,obj)
        histDict[i] = {}
        histDict[i] = profValues2D(obj, histDict[i], tag = "RvsLayervsThickness_%s" %(i), title = "R vs Layer vs Thickness for %s" %(i), axunit = "Layer", binsBoundariesX = [100, 0., 100.], ayunit = "R (cm)", binsBoundariesY=[100, 0., 300.], binsBoundariesZ=[0.,400.], weighted2D=False, verbosityLevel=verbosityLevel)

        ROOT.gStyle.SetOptStat(0)

        mycE2 = ROOT.TCanvas("P22E%s_RvsLayervsThickness_thick%s"%(GenEnergy,i), "P22E%s_RvsLayervsThickness_thick%s"%(GenEnergy,i), 500, 500)

        acustompalette()
        ex1 = ROOT.TExec("ex1","acustompalette();");
        ex1.Draw();

        histDict[i]["RvsLayervsThickness_%s"%(i)].Draw("COLZ")
        mycE2.Update()

        palette = histDict[i]["RvsLayervsThickness_%s"%(i)].GetListOfFunctions().FindObject("palette")
        if palette:
            palette.__class__ = ROOT.TPaletteAxis
            palette.SetX1NDC(0.85)
            palette.SetX2NDC(0.9)
            #palette.SetY1NDC(0.1)
            #palette.SetY2NDC(0.6)
            palette.GetAxis().SetTickSize(.01)
            palette.GetAxis().SetTitle("Si thick")
            palette.GetAxis().SetTitleOffset(0.82);
            #palette.GetAxis().LabelsOption("v")
            ROOT.gPad.Update()

        mycE2.SaveAs("%s/P22E%s_RvsLayervsThickness_thick_%s.png"%(outDir,GenEnergy,i))

    #----------------------------------------------------------
    # RvsEta
    RechitRvsEta = {}
    RechitRvsEta["CE_E_Front_120um"] = df_m.query("sClusHitsThick == 120 & sClusHitsDet == 8")[['rechit_eta', 'R']].to_numpy()
    RechitRvsEta["CE_E_Front_200um"] = df_m.query("sClusHitsThick == 200 & sClusHitsDet == 8")[['rechit_eta', 'R']].to_numpy()
    RechitRvsEta["CE_E_Front_300um"] = df_m.query("sClusHitsThick == 300 & sClusHitsDet == 8")[['rechit_eta', 'R']].to_numpy()
    RechitRvsEta["CE_H_Fine_120um"] = df_m.query("sClusHitsThick == 120 & sClusHitsDet == 9")[['rechit_eta', 'R']].to_numpy()
    RechitRvsEta["CE_H_Fine_200um"] = df_m.query("sClusHitsThick == 200 & sClusHitsDet == 9")[['rechit_eta', 'R']].to_numpy()
    RechitRvsEta["CE_H_Fine_300um"] = df_m.query("sClusHitsThick == 300 & sClusHitsDet == 9")[['rechit_eta', 'R']].to_numpy()
    RechitRvsEta["CE_H_Coarse_Scint"] = df_m.query("sClusHitsThick > 400 & sClusHitsDet == 10")[['rechit_eta', 'R']].to_numpy()

    histDict = {}
    for i, obj in RechitRvsEta.items():
        print(i,obj)
        histDict[i] = {}
        histDict[i] = histValues2D(obj, histDict[i], tag = "RvsEta_%s"%(i), title = "R vs Eta for %s" %(i), axunit = "|#eta|", binsBoundariesX = EtaBoundariesShower[i], ayunit = "R (cm)", binsBoundariesY=[100, 0., 300.], weighted2D=False, verbosityLevel=verbosityLevel)

        ROOT.gStyle.SetOptStat(0)

        mycE3 = ROOT.TCanvas("P22E%s_RvsEta_thick%s_ecut%d"%(GenEnergy,i,ecut), "P22E%s_RvsEta_thick%s_ecut%d"%(GenEnergy,i,ecut), 500, 500)

        acustompalette()
        ex1 = ROOT.TExec("ex1","acustompalette();");
        ex1.Draw();

        histDict[i]["RvsEta_%s"%(i)].Draw("COLZ")
        #histDict[i]["RvsEta_%s"%(i)].GetXaxis().SetTitleOffset(0.95) 
        mycE3.Update()

        palette = histDict[i]["RvsEta_%s"%(i)].GetListOfFunctions().FindObject("palette")
        if palette:
            palette.__class__ = ROOT.TPaletteAxis
            palette.SetX1NDC(0.85)
            palette.SetX2NDC(0.9)
            #palette.SetY1NDC(0.1)
            #palette.SetY2NDC(0.6)
            palette.GetAxis().SetTickSize(.01)
            palette.GetAxis().SetTitle("Si thick")
            palette.GetAxis().SetTitleOffset(0.8);
            #palette.GetAxis().LabelsOption("v")
            ROOT.gPad.Update()

        mycE3.SaveAs("%s/P22E%s_RvsEta_thick%s_ecut%d.png"%(outDir,GenEnergy,i,ecut))



    #----------------------------------------------------------
    # RvsLayer
    RechitRvsLayer = {}
    RechitRvsLayer["CE_E_Front_120um"] = df_m.query("sClusHitsThick == 120 & sClusHitsDet == 8")[['sClusHitsLayers', 'R']].to_numpy()
    RechitRvsLayer["CE_E_Front_200um"] = df_m.query("sClusHitsThick == 200 & sClusHitsDet == 8")[['sClusHitsLayers', 'R']].to_numpy()
    RechitRvsLayer["CE_E_Front_300um"] = df_m.query("sClusHitsThick == 300 & sClusHitsDet == 8")[['sClusHitsLayers', 'R']].to_numpy()
    RechitRvsLayer["CE_H_Fine_120um"] = df_m.query("sClusHitsThick == 120 & sClusHitsDet == 9")[['sClusHitsLayers', 'R']].to_numpy()
    RechitRvsLayer["CE_H_Fine_200um"] = df_m.query("sClusHitsThick == 200 & sClusHitsDet == 9")[['sClusHitsLayers', 'R']].to_numpy()
    RechitRvsLayer["CE_H_Fine_300um"] = df_m.query("sClusHitsThick == 300 & sClusHitsDet == 9")[['sClusHitsLayers', 'R']].to_numpy()
    RechitRvsLayer["CE_H_Coarse_Scint"] = df_m.query("sClusHitsThick > 400 & sClusHitsDet == 10")[['sClusHitsLayers', 'R']].to_numpy()

    histDict = {}
    for i, obj in RechitRvsLayer.items():
        print(i,obj)
        histDict[i] = {}
        histDict[i] = histValues2D(obj, histDict[i], tag = "RvsLayer_%s"%(i), title = "R vs Layer for %s"  %(i), axunit = "Layer", binsBoundariesX = [100, 0., 100.], ayunit = "R (cm)", binsBoundariesY=[100, 0., 300.], weighted2D=False, verbosityLevel=verbosityLevel)

        ROOT.gStyle.SetOptStat(0)

        mycE4 = ROOT.TCanvas("P22E%s_RvsLayer_thick%s_ecut%d"%(GenEnergy,i,ecut), "P22E%s_RvsLayer_thick%s_ecut%d"%(GenEnergy,i,ecut), 500, 500)

        acustompalette()
        ex1 = ROOT.TExec("ex1","acustompalette();");
        ex1.Draw();

        histDict[i]["RvsLayer_%s"%(i)].Draw("COLZ")
        #histDict[i]["RvsLayer_%s"%(i)].GetXaxis().SetTitleOffset(0.95) 
        mycE4.Update()

        palette = histDict[i]["RvsLayer_%s"%(i)].GetListOfFunctions().FindObject("palette")
        if palette:
            palette.__class__ = ROOT.TPaletteAxis
            palette.SetX1NDC(0.85)
            palette.SetX2NDC(0.9)
            #palette.SetY1NDC(0.1)
            #palette.SetY2NDC(0.6)
            palette.GetAxis().SetTickSize(.01)
            palette.GetAxis().SetTitle("Si thick")
            palette.GetAxis().SetTitleOffset(0.8);
            #palette.GetAxis().LabelsOption("v")
            ROOT.gPad.Update()

        mycE4.SaveAs("%s/P22E%s_RvsLayer_thick%s_ecut%d.png"%(outDir,GenEnergy,i,ecut))


    # Reconstructable energy 
    recostr_ene = {}

    recostr_ene["CE_E_Front_120um"] = df_m[ (df_m["sClusHitsThick"] == 120) & (df_m["sClusHitsDet"] == 8)].groupby(['EventId']).agg( recHitEneXFractionSum  = ('rechit_recostructable_energy','sum'))[["recHitEneXFractionSum"]].to_numpy().flatten()

    recostr_ene["CE_E_Front_200um"] = df_m[ (df_m["sClusHitsThick"] == 200) & (df_m["sClusHitsDet"] == 8)].groupby(['EventId']).agg( recHitEneXFractionSum  = ('rechit_recostructable_energy','sum'))[['recHitEneXFractionSum']].to_numpy().flatten()

    recostr_ene["CE_E_Front_300um"] = df_m[ (df_m["sClusHitsThick"] == 300) & (df_m["sClusHitsDet"] == 8)].groupby(['EventId']).agg( recHitEneXFractionSum  = ('rechit_recostructable_energy','sum'))[['recHitEneXFractionSum']].to_numpy().flatten()

    recostr_ene["CE_H_Fine_120um"]  = df_m[ (df_m["sClusHitsThick"] == 120) & (df_m["sClusHitsDet"] == 9)].groupby(['EventId']).agg( recHitEneXFractionSum  = ('rechit_recostructable_energy','sum'))[['recHitEneXFractionSum']].to_numpy().flatten()

    recostr_ene["CE_H_Fine_200um"]  = df_m[ (df_m["sClusHitsThick"] == 200) & (df_m["sClusHitsDet"] == 9)].groupby(['EventId']).agg( recHitEneXFractionSum  = ('rechit_recostructable_energy','sum'))[['recHitEneXFractionSum']].to_numpy().flatten()

    recostr_ene["CE_H_Fine_300um"]  = df_m[ (df_m["sClusHitsThick"] == 300) & (df_m["sClusHitsDet"] == 9)].groupby(['EventId']).agg( recHitEneXFractionSum  = ('rechit_recostructable_energy','sum'))[['recHitEneXFractionSum']].to_numpy().flatten()

    recostr_ene["CE_H_Coarse_Scint"] = df_m[ (df_m["sClusHitsThick"] > 400) & (df_m["sClusHitsDet"] == 10)].groupby(['EventId']).agg( recHitEneXFractionSum  = ('rechit_recostructable_energy','sum'))[['recHitEneXFractionSum']].to_numpy().flatten()

    fittree_recostr_ene = {}

    for i, obj in recostr_ene.items():
        print(i,obj)
        fittree_recostr_ene[i] = ROOT.TTree("ttree_reconstructable_energy_e%s_%s" %(GenEnergy, i), "results")
        maketree(obj, "reconstructable_energy_%s"%(i), fittree_recostr_ene[i])
        #maketree(obj, "reconstructable_energy_%s"%(i), tree)

    # Uncalibrated energy 
    uncalib_ene = {}

    uncalib_ene["CE_E_Front_120um"] = df_m[ (df_m["sClusHitsThick"] == 120) & (df_m["sClusHitsDet"] == 8)].groupby(['EventId']).agg( uncalibRecHitEneSum  = ('rechit_uncalib_energy','sum'))[["uncalibRecHitEneSum"]].to_numpy().flatten()

    uncalib_ene["CE_E_Front_200um"] = df_m[ (df_m["sClusHitsThick"] == 200) & (df_m["sClusHitsDet"] == 8)].groupby(['EventId']).agg( uncalibRecHitEneSum  = ('rechit_uncalib_energy','sum'))[['uncalibRecHitEneSum']].to_numpy().flatten()

    uncalib_ene["CE_E_Front_300um"] = df_m[ (df_m["sClusHitsThick"] == 300) & (df_m["sClusHitsDet"] == 8)].groupby(['EventId']).agg( uncalibRecHitEneSum  = ('rechit_uncalib_energy','sum'))[['uncalibRecHitEneSum']].to_numpy().flatten()

    uncalib_ene["CE_H_Fine_120um"]  = df_m[ (df_m["sClusHitsThick"] == 120) & (df_m["sClusHitsDet"] == 9)].groupby(['EventId']).agg( uncalibRecHitEneSum  = ('rechit_uncalib_energy','sum'))[['uncalibRecHitEneSum']].to_numpy().flatten()

    uncalib_ene["CE_H_Fine_200um"]  = df_m[ (df_m["sClusHitsThick"] == 200) & (df_m["sClusHitsDet"] == 9)].groupby(['EventId']).agg( uncalibRecHitEneSum  = ('rechit_uncalib_energy','sum'))[['uncalibRecHitEneSum']].to_numpy().flatten()

    uncalib_ene["CE_H_Fine_300um"]  = df_m[ (df_m["sClusHitsThick"] == 300) & (df_m["sClusHitsDet"] == 9)].groupby(['EventId']).agg( uncalibRecHitEneSum  = ('rechit_uncalib_energy','sum'))[['uncalibRecHitEneSum']].to_numpy().flatten()

    uncalib_ene["CE_H_Coarse_Scint"] = df_m[ (df_m["sClusHitsThick"] > 400) & (df_m["sClusHitsDet"] == 10)].groupby(['EventId']).agg( uncalibRecHitEneSum  = ('rechit_uncalib_energy','sum'))[['uncalibRecHitEneSum']].to_numpy().flatten()

    fittree_uncalib_ene = {}

    for i, obj in uncalib_ene.items():
        print(i,obj)
        fittree_uncalib_ene[i] = ROOT.TTree("ttree_uncalib_energy_e%s_%s" %(GenEnergy, i), "results")
        maketree(obj, "uncalib_energy_%s"%(i), fittree_uncalib_ene[i])
        #maketree(obj, "uncalib_energy_%s"%(i), tree)

    #Time to save the trees. 
    outfile = ROOT.TFile("{}/{}.root".format(outDir, output), "recreate")

    #The final tree as we want it
    #finaltree = ROOT.TTree("ttree_e%s" %(GenEnergy), "A tree with the needed info for resolution plots")

    #finaldf = {}
    dataout = {}
    maxlen = max(len(recostr_ene["CE_E_Front_120um"]),len(recostr_ene["CE_E_Front_200um"]),len(recostr_ene["CE_E_Front_300um"]),len(recostr_ene["CE_H_Fine_120um"]),len(recostr_ene["CE_H_Fine_200um"]),len(recostr_ene["CE_H_Fine_300um"]) ,len(recostr_ene["CE_H_Coarse_Scint"]))
    for region in ["CE_E_Front_120um","CE_E_Front_200um","CE_E_Front_300um","CE_H_Fine_120um","CE_H_Fine_200um","CE_H_Fine_300um","CE_H_Coarse_Scint"]:
        padtocreatedf = maxlen - len(recostr_ene[region])
        dummy = np.empty(padtocreatedf)
        dummy.fill(99999999)
        if padtocreatedf != 0: 
            dataout['reconstructable_energy_%s' %(region)] = np.append(recostr_ene[region], dummy )      
        else: 
            dataout['reconstructable_energy_%s' %(region)] = recostr_ene[region]

    maxlen = max(len(uncalib_ene["CE_E_Front_120um"]),len(uncalib_ene["CE_E_Front_200um"]),len(uncalib_ene["CE_E_Front_300um"]),len(uncalib_ene["CE_H_Fine_120um"]),len(uncalib_ene["CE_H_Fine_200um"]),len(uncalib_ene["CE_H_Fine_300um"]) ,len(uncalib_ene["CE_H_Coarse_Scint"]))
    for region in ["CE_E_Front_120um","CE_E_Front_200um","CE_E_Front_300um","CE_H_Fine_120um","CE_H_Fine_200um","CE_H_Fine_300um","CE_H_Coarse_Scint"]:
        padtocreatedf = maxlen - len(uncalib_ene[region])
        dummy = np.empty(padtocreatedf)
        dummy.fill(99999999)
        if padtocreatedf != 0: 
            dataout['uncalib_energy_%s' %(region)] = np.append(uncalib_ene[region], dummy )      
        else: 
            dataout['uncalib_energy_%s' %(region)] = uncalib_ene[region]

    maxlen = max(len(rHxEnFrSum_m["CE_E_Front_120um"]),len(rHxEnFrSum_m["CE_E_Front_200um"]),len(rHxEnFrSum_m["CE_E_Front_300um"]),len(rHxEnFrSum_m["CE_H_Fine_120um"]),len(rHxEnFrSum_m["CE_H_Fine_200um"]),len(rHxEnFrSum_m["CE_H_Fine_300um"]) ,len(rHxEnFrSum_m["CE_H_Coarse_Scint"]))
    for region in ["CE_E_Front_120um","CE_E_Front_200um","CE_E_Front_300um","CE_H_Fine_120um","CE_H_Fine_200um","CE_H_Fine_300um","CE_H_Coarse_Scint"]:
        padtocreatedf = maxlen - len(rHxEnFrSum_m[region])
        dummy = np.empty(padtocreatedf)
        dummy.fill(99999999)
        if padtocreatedf != 0: 
            dataout['SumEoverEgen_%s' %(region)] = np.append(rHxEnFrSum_m[region], dummy )      
        else: 
            dataout['SumEoverEgen_%s' %(region)] = rHxEnFrSum_m[region]


    #maxlen = max(len(rHxEnFrSum_u["CE_E_Front_120um"]),len(rHxEnFrSum_u["CE_E_Front_200um"]),len(rHxEnFrSum_u["CE_E_Front_300um"]),len(rHxEnFrSum_u["CE_H_Fine_120um"]),len(rHxEnFrSum_u["CE_H_Fine_200um"]),len(rHxEnFrSum_u["CE_H_Fine_300um"]) ,len(rHxEnFrSum_u["CE_H_Coarse_Scint"]))
    #for region in ["CE_E_Front_120um","CE_E_Front_200um","CE_E_Front_300um","CE_H_Fine_120um","CE_H_Fine_200um","CE_H_Fine_300um","CE_H_Coarse_Scint"]:
    #    padtocreatedf = maxlen - len(rHxEnFrSum_u[region])
    #    dummy = np.empty(padtocreatedf)
    #    dummy.fill(99999999)
    #    if padtocreatedf != 0: 
    #        dataout['SumEoverEgen_Unmatched_%s' %(region)] = np.append(rHxEnFrSum_u[region], dummy )      
    #    else: 
    #        dataout['SumEoverEgen_Unmatched_%s' %(region)] = rHxEnFrSum_u[region]




    finaldf = ROOT.RDF.MakeNumpyDataFrame(dataout)
    #finaldf.Display().Print()
    #for region in ["CE_E_Front_120um","CE_E_Front_200um","CE_E_Front_300um","CE_H_Fine_120um","CE_H_Fine_200um","CE_H_Fine_300um","CE_H_Coarse_Scint"]:
        #finaldf[region] = ROOT.RDataFrame(fittree_m[region])
        #print(finaldf[region])
    

    #for region in ["CE_E_Front_120um","CE_E_Front_200um","CE_E_Front_300um","CE_H_Fine_120um","CE_H_Fine_200um","CE_H_Fine_300um","CE_H_Coarse_Scint"]:
        
    #    for i in fittree[region].GetEntries(): 
    #        print(region, i)
        
        #finaltree = fittree[region].CloneTree(0)
        #fittree[region].GetListoOfClones().Remove(finaltree)
        #finaltree.ResetBranchAddress()
        #finaltree->SetBranchAddress("SumEoverEgen_%s"%(region), "SumEoverEgen_%s"%(region))
        
        #finaltree.GetBranch("SumEoverEgen_%s"%(region))
        #finaltree.CopyEntries(fittree[region])
        
        #fittree[region].GetBranch("SumEoverEgen_%s"%(region))
        

        #finaltree.Branch("SumEoverEgen_%s"%(region), rHxEnFrSum_tree_m[region], "SumEoverEgen_%s/F"%(region))

        #ttree_SumEoverEgen_CE_E_Front_120um_e100->GetBranch("SumEoverEgen_CE_E_Front_120um")

    #Write tree
    #tree.SetDirectory(outfile)
    #tree.Write()
    #for i, obj in fittree_u.items(): fittree_u[i].Write()
    #for i, obj in fittree_m.items(): fittree_m[i].Write()
    #for i, obj in fittree_recostr_ene.items(): fittree_recostr_ene[i].Write()
    #for i, obj in fittree_uncalib_ene.items(): fittree_uncalib_ene[i].Write()
    #finaltree.Write()
    finaldf.Snapshot('finaltree', "{}/{}.root".format(outDir, output))

'''
    maketree(rHxEnFrSum_m, rHxEnFrSum_u, recostr_ene, uncalib_ene)

    rHxEnFrSum_tree_m = {}
    rHxEnFrSum_tree_u = {}
    recostr_tree_ene = {}
    uncalib_tree_ene = {}

    for region in ["CE_E_Front_120um","CE_E_Front_200um","CE_E_Front_300um","CE_H_Fine_120um","CE_H_Fine_200um","CE_H_Fine_300um","CE_H_Coarse_Scint"]: 
        rHxEnFrSum_tree_m[region] = np.empty((1), dtype="float32")
        rHxEnFrSum_tree_u[region] = np.empty((1), dtype="float32")
        recostr_tree_ene[region] = np.empty((1), dtype="float32")
        uncalib_tree_ene[region] = np.empty((1), dtype="float32")

    for region in ["CE_E_Front_120um","CE_E_Front_200um","CE_E_Front_300um","CE_H_Fine_120um","CE_H_Fine_200um","CE_H_Fine_300um","CE_H_Coarse_Scint"]: 
        tree.Branch("SumEoverEgen_%s"%(region), rHxEnFrSum_tree_m[region], "SumEoverEgen_%s/F"%(region))
        tree.Branch("SumEoverEgen_Unmatched_%s"%(region), rHxEnFrSum_tree_u[region], "SumEoverEgen_Unmatched_%s/F"%(region))
        tree.Branch("reconstructable_energy_%s"%(region), recostr_tree_ene[region], "reconstructable_energy_%s/F"%(region))
        tree.Branch("uncalib_energy_%s"%(region), uncalib_tree_ene[region], "uncalib_energy_%s/F"%(region))

    for i in rHxEnFrSum_tree_m


def maketree(rHxEnFrSum_m, rHxEnFrSum_u, recostr_ene, uncalib_ene):

    rHxEnFrSum_tree_m_CE_E_120um = np.empty((1), dtype="float32")
    rHxEnFrSum_tree_m_CE_E_200um = np.empty((1), dtype="float32")
    rHxEnFrSum_tree_m_CE_E_300um = np.empty((1), dtype="float32")
    rHxEnFrSum_tree_m_CE_H_120um = np.empty((1), dtype="float32")
    rHxEnFrSum_tree_m_CE_H_200um = np.empty((1), dtype="float32")
    rHxEnFrSum_tree_m_CE_H_300um = np.empty((1), dtype="float32")
    rHxEnFrSum_tree_m_CE_H_Scint = np.empty((1), dtype="float32")

    tree.Branch("SumEoverEgen_CE_E_120um", rHxEnFrSum_tree_m_CE_E_120um, "SumEoverEgen_CE_E_120um/F")
    tree.Branch("SumEoverEgen_CE_E_200um", rHxEnFrSum_tree_m_CE_E_200um, "SumEoverEgen_CE_E_200um/F")
    tree.Branch("SumEoverEgen_CE_E_300um", rHxEnFrSum_tree_m_CE_E_300um, "SumEoverEgen_CE_E_300um/F")
    tree.Branch("SumEoverEgen_CE_H_120um", rHxEnFrSum_tree_m_CE_H_120um, "SumEoverEgen_CE_H_120um/F")
    tree.Branch("SumEoverEgen_CE_H_200um", rHxEnFrSum_tree_m_CE_H_200um, "SumEoverEgen_CE_H_200um/F")
    tree.Branch("SumEoverEgen_CE_H_300um", rHxEnFrSum_tree_m_CE_H_300um, "SumEoverEgen_CE_H_300um/F")
    tree.Branch("SumEoverEgen_CE_H_Scint", rHxEnFrSum_tree_m_CE_H_Scint, "SumEoverEgen_CE_H_Scint/F")

    for i in rHxEnFrSum_m["CE_E_Front_120um"]: 
       
        outbranch[0] = i
        tree.Fill()
'''
    
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


#---------------------------------------------------------------------------------------------------
def trackstersAssociatorsPlots(df_ts,df_st,tree,maxEvents,outDir,output,verbosityLevel = 0):

    #-------------------------------------------------------------------------
    #Score - TS to SimTS
    #
    histDict_score = {}
    histDict_score = histValue1D(df_ts[['score_trackster2caloparticle']].to_numpy().flatten(), histDict_score, tag = "score_TS2SimTS", title ="Score of Trackster per SimTrackster",  axunit = "score Reco-to-Sim", binsBoundariesX = [200, 0, 2], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_score, outDir, output, tree, verbosityLevel, setLogY = True)

    #-------------------------------------------------------------------------
    #Shared energy - TS to SimTS
    #
    histDict_shene = {}
    histDict_shene = histValue1D(df_ts[['sharedenergy_trackster2caloparticle']].to_numpy().flatten(), histDict_shene, tag = "sharedenergy_TS2SimTS", title ="Shared Energy of Trackster per SimTrackster",  axunit = "shared Reco energy fraction", binsBoundariesX = [200, 0, 2], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_shene, outDir, output, tree, verbosityLevel, setLogY = True)

    #-------------------------------------------------------------------------
    #numberOfHitsInTS: Total number of hits in trackster, that is the sum of the hits of all LCs that
    #                  belongs to the trackster. The histo is filled #tracksters in event x #events
    #
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #tracksters in event x #related SimTracksters x #events. Here we undo the per SimTS entries
    # for the numberOfHitsInTS by grouping and taking the mean. 
    histDict_numhitsints = {}
    df_ts_perEvent_perST = df_ts.groupby(['EventId','scId']).agg( numberOfHitsInTSperTSperEvent  = ('numberOfHitsInTS','mean'))
    print(df_ts_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_ts[ (df_ts['EventId'] == 1099) & (df_ts['scId'] == 0)][['numberOfHitsInTS']])
    #df_ts_perEvent_perST[['numberOfHitsInTSperTSperEvent']]
    
    histDict_numhitsints = histValue1D(df_ts_perEvent_perST[['numberOfHitsInTSperTSperEvent']].to_numpy().flatten(), histDict_numhitsints, tag = "numberOfHitsInTSperTSperEvent", title ="Number of hits in Trackster",  axunit = "#hits in TS", binsBoundariesX = [200, 0, 1000], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_numhitsints, outDir, output, tree, verbosityLevel, setLogY = False)

    #-------------------------------------------------------------------------
    #numberOfNoiseHitsInTS: Total number of noise hits in trackster. In the case of Linking, this is the total number of hits 
    #                       of the trackster not related to a simhit of any SimCluster of a CaloParticle (SimTrackster from
    #                       CPs). In the PR case it counts the number of rechits of a LC of a trackster without any matched
    #                       simhit. This is really stronger than just a simple noise unmatched rechit, since in the
    #                       detIdSimTSId_Map at least one simhit to belong to a LC, the whole LC is considered. 
    #                       The histo is filled #tracksters in event x #events
    #
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #tracksters in event x #related SimTracksters x #events. Here we undo the per SimTS entries
    # for the numberOfNoiseHitsInTS by grouping and taking the mean. 
    histDict_numnoisehitsints = {}
    df_ts_perEvent_perST = df_ts.groupby(['EventId','scId']).agg( numberOfNoiseHitsInTSperTSperEvent  = ('numberOfNoiseHitsInTS','mean'))
    print(df_ts_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_ts[ (df_ts['EventId'] == 1099) & (df_ts['scId'] == 0)][['numberOfNoiseHitsInTS']])
    #df_ts_perEvent_perST[['numberOfNoiseHitsInTSperTSperEvent']]
    
    histDict_numnoisehitsints = histValue1D(df_ts_perEvent_perST[['numberOfNoiseHitsInTSperTSperEvent']].to_numpy().flatten(), histDict_numnoisehitsints, tag = "numberOfNoiseHitsInTSperTSperEvent", title ="Number of noise hits in Trackster",  axunit = "#noise hits in TS", binsBoundariesX = [200, 0, 1000], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_numnoisehitsints, outDir, output, tree, verbosityLevel, setLogY = False)

    #-------------------------------------------------------------------------
    #energyFractionOfTSinCP: 
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #tracksters in event x #related SimTracksters x #events. Here we undo the per SimTS entries
    # for the energyFractionOfTSinCP by grouping and taking the mean. 
    histDict_enefroftsincp = {}
    df_ts_perEvent_perST = df_ts.groupby(['EventId','scId']).agg( energyFractionOfTSinCPperTSperEvent  = ('energyFractionOfTSinCP','mean'))
    print(df_ts_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_ts[ (df_ts['EventId'] == 1099) & (df_ts['scId'] == 0)][['energyFractionOfTSinCP']])
    #df_ts_perEvent_perST[['energyFractionOfTSinCPperTSperEvent']]
    
    histDict_enefroftsincp = histValue1D(df_ts_perEvent_perST[['energyFractionOfTSinCPperTSperEvent']].to_numpy().flatten(), histDict_enefroftsincp, tag = "energyFractionOfTSinCPperTSperEvent", title ="energyFractionOfTSinCP",  axunit = "energyFractionOfTSinCP", binsBoundariesX = [200, 0, 2], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_enefroftsincp, outDir, output, tree, verbosityLevel, setLogY = True)

    #-------------------------------------------------------------------------
    #energyFractionOfCPinTS: 
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #tracksters in event x #related SimTracksters x #events. Here we undo the per SimTS entries
    # for the energyFractionOfCPinTS by grouping and taking the mean. 
    histDict_enefrofcpints = {}
    df_ts_perEvent_perST = df_ts.groupby(['EventId','scId']).agg( energyFractionOfCPinTSperTSperEvent  = ('energyFractionOfCPinTS','mean'))
    print(df_ts_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_ts[ (df_ts['EventId'] == 1099) & (df_ts['scId'] == 0)][['energyFractionOfCPinTS']])
    #df_ts_perEvent_perST[['energyFractionOfCPinTSperTSperEvent']]
    
    histDict_enefrofcpints = histValue1D(df_ts_perEvent_perST[['energyFractionOfCPinTSperTSperEvent']].to_numpy().flatten(), histDict_enefrofcpints, tag = "energyFractionOfCPinTSperTSperEvent", title ="energyFractionOfCPinTS",  axunit = "energyFractionOfCPinTS", binsBoundariesX = [200, 0, 2], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_enefrofcpints, outDir, output, tree, verbosityLevel, setLogY = True)

    #-------------------------------------------------------------------------
    #numofvertices: 
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #tracksters in event x #related SimTracksters x #events. Here we undo the per SimTS entries
    # for the numofvertices by grouping and taking the mean. 
    histDict_numofvert = {}
    df_ts_perEvent_perST = df_ts.groupby(['EventId','scId']).agg( numofverticesperTSperEvent  = ('numofvertices','mean'))
    print(df_ts_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_ts[ (df_ts['EventId'] == 1099) & (df_ts['scId'] == 0)][['numofvertices']])
    #df_ts_perEvent_perST[['numofverticesperTSperEvent']]
    
    histDict_numofvert = histValue1D(df_ts_perEvent_perST[['numofverticesperTSperEvent']].to_numpy().flatten(), histDict_numofvert, tag = "numofverticesperTSperEvent", title ="Number of vertices in Trackster",  axunit = "#vertices", binsBoundariesX = [200, 0, 1000], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_numofvert, outDir, output, tree, verbosityLevel, setLogY = False)

    #-------------------------------------------------------------------------
    #raw_energy: 
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #tracksters in event x #related SimTracksters x #events. Here we undo the per SimTS entries
    # for the raw_energy by grouping and taking the mean. 
    histDict_rawene = {}
    df_ts_perEvent_perST = df_ts.groupby(['EventId','scId']).agg( raw_energyperTSperEvent  = ('raw_energy','mean'))
    print(df_ts_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_ts[ (df_ts['EventId'] == 1099) & (df_ts['scId'] == 0)][['raw_energy']])
    #df_ts_perEvent_perST[['raw_energyperTSperEvent']]
    
    histDict_rawene = histValue1D(df_ts_perEvent_perST[['raw_energyperTSperEvent']].to_numpy().flatten(), histDict_rawene, tag = "raw_energyperTSperEvent", title ="Raw Energy of Trackster",  axunit = "Raw Energy (GeV)", binsBoundariesX = [200, 0, 1000], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_rawene, outDir, output, tree, verbosityLevel, setLogY = False)

    #-------------------------------------------------------------------------
    #Score - SimTS to TS
    #
    histDict_score = {}
    histDict_score = histValue1D(df_st[['sts_score_caloparticle2trackster']].to_numpy().flatten(), histDict_score, tag = "score_SimTS2TS", title ="Score of SimTrackster per Trackster",  axunit = "score Sim-to-Reco", binsBoundariesX = [200, 0, 2], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_score, outDir, output, tree, verbosityLevel, setLogY = True)

    #-------------------------------------------------------------------------
    #Shared energy - SimTS to TS
    #
    histDict_shene = {}
    histDict_shene = histValue1D(df_st[['sts_sharedenergy_caloparticle2trackster']].to_numpy().flatten(), histDict_shene, tag = "sharedenergy_SimTS2TS", title ="Shared Energy of SimTrackster per Trackster",  axunit = "shared Sim energy fraction", binsBoundariesX = [200, 0, 2], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_shene, outDir, output, tree, verbosityLevel, setLogY = True)

    #-------------------------------------------------------------------------
    #sts_SimEnergy:
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #simtracksters in event x #related Tracksters x #events. Here we undo the per TS entries
    # for the sts_SimEnergy by grouping and taking the mean. 
    histDict_sts_SimEnergy = {}
    df_st_perEvent_perST = df_st.groupby(['EventId','sts_ts_id']).agg( sts_SimEnergyperSTSperEvent  = ('sts_SimEnergy','mean'))
    print(df_st_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_st[ (df_st['EventId'] == 1099) & (df_st['sts_ts_id'] == 0)][['sts_SimEnergy']])
    #df_st_perEvent_perST[['sts_SimEnergySTSperEvent']]
    
    histDict_sts_SimEnergy = histValue1D(df_st_perEvent_perST[['sts_SimEnergyperSTSperEvent']].to_numpy().flatten(), histDict_sts_SimEnergy, tag = "SimEnergyperSTSperEvent", title ="SimEnergy of SimTrackster",  axunit = "#simenergy in STS", binsBoundariesX = [200, 0, 1000], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_sts_SimEnergy, outDir, output, tree, verbosityLevel, setLogY = False)

    #-------------------------------------------------------------------------
    #sts_SimEnergyWeight:
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #simtracksters in event x #related Tracksters x #events. Here we undo the per TS entries
    # for the sts_SimEnergyWeight by grouping and taking the mean. 
    histDict_sts_SimEnergyWeight = {}
    df_st_perEvent_perST = df_st.groupby(['EventId','sts_ts_id']).agg( sts_SimEnergyWeightperSTSperEvent  = ('sts_SimEnergyWeight','mean'))
    print(df_st_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_st[ (df_st['EventId'] == 1099) & (df_st['sts_ts_id'] == 0)][['sts_SimEnergyWeight']])
    #df_st_perEvent_perST[['sts_SimEnergyWeightSTSperEvent']]
    
    histDict_sts_SimEnergyWeight = histValue1D(df_st_perEvent_perST[['sts_SimEnergyWeightperSTSperEvent']].to_numpy().flatten(), histDict_sts_SimEnergyWeight, tag = "SimEnergyWeightperSTSperEvent", title ="SimEnergyWeight of SimTrackster",  axunit = "#simenergyweight in STS", binsBoundariesX = [200, 0, 1000], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_sts_SimEnergyWeight, outDir, output, tree, verbosityLevel, setLogY = False)

    #-------------------------------------------------------------------------
    #sts_eta:
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #simtracksters in event x #related Tracksters x #events. Here we undo the per TS entries
    # for the sts_eta by grouping and taking the mean. 
    histDict_sts_eta = {}
    df_st_perEvent_perST = df_st.groupby(['EventId','sts_ts_id']).agg( sts_etaperSTSperEvent  = ('sts_eta','mean'))
    print(df_st_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_st[ (df_st['EventId'] == 1099) & (df_st['sts_ts_id'] == 0)][['sts_eta']])
    #df_st_perEvent_perST[['sts_etaSTSperEvent']]
    
    histDict_sts_eta = histValue1D(df_st_perEvent_perST[['sts_etaperSTSperEvent']].to_numpy().flatten(), histDict_sts_eta, tag = "EtaofSTS", title ="SimTrackster eta",  axunit = "eta", binsBoundariesX = [200, -5.0, 5.0], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_sts_eta, outDir, output, tree, verbosityLevel, setLogY = False)

    #-------------------------------------------------------------------------
    #sts_phi:
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #simtracksters in event x #related Tracksters x #events. Here we undo the per TS entries
    # for the sts_phi by grouping and taking the mean. 
    histDict_sts_phi = {}
    df_st_perEvent_perST = df_st.groupby(['EventId','sts_ts_id']).agg( sts_phiperSTSperEvent  = ('sts_phi','mean'))
    print(df_st_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_st[ (df_st['EventId'] == 1099) & (df_st['sts_ts_id'] == 0)][['sts_phi']])
    #df_st_perEvent_perST[['sts_phiSTSperEvent']]
    
    histDict_sts_phi = histValue1D(df_st_perEvent_perST[['sts_phiperSTSperEvent']].to_numpy().flatten(), histDict_sts_phi, tag = "PhiofSTS", title ="SimTrackster phi",  axunit = "phi", binsBoundariesX = [200, -5.0, 5.0], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_sts_phi, outDir, output, tree, verbosityLevel, setLogY = False)

    #-------------------------------------------------------------------------
    #sts_pt:
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #simtracksters in event x #related Tracksters x #events. Here we undo the per TS entries
    # for the sts_pt by grouping and taking the mean. 
    histDict_sts_pt = {}
    df_st_perEvent_perST = df_st.groupby(['EventId','sts_ts_id']).agg( sts_ptperSTSperEvent  = ('sts_pt','mean'))
    print(df_st_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_st[ (df_st['EventId'] == 1099) & (df_st['sts_ts_id'] == 0)][['sts_pt']])
    #df_st_perEvent_perST[['sts_ptSTSperEvent']]
    
    histDict_sts_pt = histValue1D(df_st_perEvent_perST[['sts_ptperSTSperEvent']].to_numpy().flatten(), histDict_sts_pt, tag = "PtofSTS", title ="SimTrackster pt",  axunit = "pt", binsBoundariesX = [200, 0, 1000], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_sts_pt, outDir, output, tree, verbosityLevel, setLogY = False)

    #-------------------------------------------------------------------------
    #sts_raw_energy:
    #In order to have a dataframe with same number of lines certain variables where saved more times than it was
    #needed, namely #simtracksters in event x #related Tracksters x #events. Here we undo the per TS entries
    # for the sts_raw_energy by grouping and taking the mean. 
    histDict_sts_raw_energy = {}
    df_st_perEvent_perST = df_st.groupby(['EventId','sts_ts_id']).agg( sts_raw_energyperSTSperEvent  = ('sts_raw_energy','mean'))
    print(df_st_perEvent_perST)
    #These two lines below is to double check that I indeed take what I expect. 
    #print(df_st[ (df_st['EventId'] == 1099) & (df_st['sts_ts_id'] == 0)][['sts_raw_energy']])
    #df_st_perEvent_perST[['sts_raw_energySTSperEvent']]
    
    histDict_sts_raw_energy = histValue1D(df_st_perEvent_perST[['sts_raw_energyperSTSperEvent']].to_numpy().flatten(), histDict_sts_raw_energy, tag = "Raw_EnergyofSTS", title ="SimTrackster raw_energy",  axunit = "raw_energy", binsBoundariesX = [200, 0., 1000.], ayunit = "", verbosityLevel=verbosityLevel)
    histPrintSaveAll(histDict_sts_raw_energy, outDir, output, tree, verbosityLevel, setLogY = False)

    







    
 
  

   

import time
import CMS_lumi, tdrstyle
from ROOT import *
from array import *
import math
import numpy as np
from optparse import OptionParser

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = ""
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod=4

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-m','--mshift', action='store_true', dest='doMassShiftFit', default=True, help='Fit (reco-gen)/reco instead of pruned mass!')
parser.add_option('-p','--pruning', action='store_true', dest='doPruning', default=False, help='Fit pruned mass')
parser.add_option('--fitGen', action='store_true', dest='fitGenMass', default=False, help='Fit gen mass')

(options, args) = parser.parse_args()

if options.noX: gROOT.SetBatch(True)

prefix = '/mnt/t3nfs01/data01/shome/thaarres/EXOVVAnalysisRunII/AnalysisOutput/80X/'
gStyle.SetOptFit(1)

def get_line(xmin,xmax,ymin,ymax,style):
   line = TLine(xmin,ymin,xmax,ymax)
   line.SetLineColor(kRed)
   line.SetLineStyle(style)
   line.SetLineWidth(2)
   return line
   
def getCanvas():
  c = TCanvas("c","c",800,800)
  c.GetWindowHeight()
  c.GetWindowWidth()
  c.SetTitle("")
  return c
  
  
def get_palette(mode):

 palette = {}
 palette['gv'] = [] 
 
 colors = ['#40004b','#762a83','#9970ab','#de77ae','#a6dba0','#5aae61','#1b7837','#00441b','#92c5de','#4393c3','#2166ac','#053061']
 colors = ['#762a83','#de77ae','#a6dba0','#92c5de','#4393c3','#2166ac','#053061','#40004b','#762a83','#9970ab','#de77ae']
 colors = ['#FF420E','#80BD9E','#a6dba0','#92c5de','#4393c3','#2166ac','#053061','#40004b','#762a83','#9970ab','#de77ae''#40004b','#762a83','#9970ab']

 for c in colors:
  palette['gv'].append(c)
 
 return palette[mode]
 
palette = get_palette('gv')
col = TColor()

masses = [1000,1200,1400,1600,1800,2000,2500,3000,4000,4500]  
masses = [600,800,1000,1400,1800,2000,2500,3000,3500,4000,4500] #Ops!! 400 masspoint is named for convenience and is actually SM WW, not signal sample!

hCentral = 'massShift_softdrop_eta1v3'
hForward = 'massShift_softdrop_etaUP1v3'
  
lineStyle = [1,1,1,1,3,3,3,3]


signals = ["BulkWW","BulkZZ","ZprimeWW","WprimeWZ"]
signals = ["BulkWW"]


for signal in signals:
 
  filelist = []
  histosCEN = []
  histosFOR = []
  fits = []
  
  l = TLegend(0.7650754,0.7564767,0.8065327,0.876943)
  l.SetTextSize(0.035)
  l.SetLineColor(0)
  l.SetShadowColor(0)
  l.SetLineStyle(1)
  l.SetLineWidth(1)
  l.SetFillColor(0)
  l.SetFillStyle(0)
  l.SetMargin(0.35)
  
  meansCentral = []
  meanErrCentral = []
  sigmasCentral = []
  sigErrCentral = []
  ptsCEN = []
  ptErrCEN = []
  meansForward = []
  meanErrForward = []
  sigmasForward = []
  sigErrForward = []
  ptsFOR = []
  ptErrFOR = []
  masspoints = [] 
  
  for m in masses:
    filename = prefix + 'ExoDiBosonAnalysis.' + signal + '_13TeV_' + "%s"%m + 'GeV.VV.root'
    filetmp = TFile.Open(filename,"READ")
    filetmp.SetName(filename)
    print "opening " ,filename
    filelist.append(filetmp)
    masspoints.append(m)
  i = -1
  for filetmp in filelist:
    i += 1
    print "Masspoint = " ,masspoints[i]
    print "Fit central bin (eta<1.3)" 
    #Central
    histtmpC = TH1F(filetmp.Get(hCentral))
    histtmpC.SetName("Central%i"%masspoints[i])
    # histtmp.Scale(1./histtmp.Integral())
    maxbin=0
    maxcontent=0
    startbin = -0.4
    for b in range(histtmpC.GetXaxis().GetNbins()):
      if histtmpC.GetXaxis().GetBinCenter(b+1) > startbin and histtmpC.GetBinContent(b+1)>maxcontent:
        maxbin = b
        maxcontent = histtmpC.GetBinContent(b+1)   
    tmpmean = histtmpC.GetXaxis().GetBinCenter(maxbin)  
    tmpwidth = 0.8  
    g1 = TF1("g1CEN%i"%masspoints[i],"gaus", tmpmean-tmpwidth,tmpmean+tmpwidth)
    histtmpC.Fit(g1, "SR")
    tmpmean = g1.GetParameter(1)
    tmpwidth = g1.GetParameter(2)
    g1 = TF1("g1CEN%i"%masspoints[i],"gaus", tmpmean-tmpwidth,tmpmean+tmpwidth)
    histtmpC.Fit(g1, "SR")
    tmpmean = g1.GetParameter(1)
    tmpwidth = g1.GetParameter(2)
    g1 = TF1("g1CEN%i"%masspoints[i],"gaus", tmpmean-tmpwidth,tmpmean+tmpwidth)
    histtmpC.Fit(g1, "SR")
    print ""
    histosCEN.append(histtmpC)
    mean    = g1.GetParameter(1)
    meanerr = g1.GetParError(1)
    sigma    = g1.GetParameter(2)
    sigmaerr = g1.GetParError(2)
    meansCentral.append(mean)
    meanErrCentral.append(meanerr)
    sigmasCentral.append(sigma)
    sigErrCentral.append(sigmaerr)
    
    print "Fit forward bin (eta>1.3)" 
    #Forward
    histtmpF = TH1F(filetmp.Get(hForward))
    histtmpF.SetName("Forward%i"%masspoints[i])
    histtmpF.Rebin(4)
    maxbin=0
    maxcontent=0
    startbin = 60.
    if options.doMassShiftFit: startbin = -0.4
    for b in range(histtmpF.GetXaxis().GetNbins()):
      if histtmpF.GetXaxis().GetBinCenter(b+1)>startbin and histtmpF.GetBinContent(b+1)>maxcontent:
        maxbin = b
        maxcontent = histtmpF.GetBinContent(b+1)
    tmpmean = histtmpF.GetXaxis().GetBinCenter(maxbin)
    tmpwidth = 0.8   
    if masspoints[i] >= 3000:
      tmpwidth = tmpwidth*1.5
    g1 = TF1("g1FOR%i"%masspoints[i],"gaus", tmpmean-tmpwidth,tmpmean+tmpwidth)
    histtmpF.Fit(g1, "SR")
    tmpmean = g1.GetParameter(1)
    tmpwidth = g1.GetParameter(2)
    if masspoints[i] >= 3000:
      tmpwidth = tmpwidth*1.5
    g1 = TF1("g1FOR%i"%masspoints[i],"gaus", tmpmean-tmpwidth,tmpmean+tmpwidth)
    histtmpF.Fit(g1, "SR")
    tmpmean = g1.GetParameter(1)
    tmpwidth = g1.GetParameter(2)
    if masspoints[i] >= 3000:
      tmpwidth = tmpwidth*1.5
    g1 = TF1("g1FOR%i"%masspoints[i],"gaus", tmpmean-tmpwidth,tmpmean+tmpwidth)
    histtmpF.Fit(g1, "SR");
    print ""
    histosFOR.append(histtmpF)
    mean    = g1.GetParameter(1)
    meanerr = g1.GetParError(1)
    sigma    = g1.GetParameter(2)
    sigmaerr = g1.GetParError(2)
    meansForward.append(mean)
    meanErrForward.append(meanerr)
    sigmasForward.append(sigma)
    sigErrForward.append(sigmaerr)
    
    ptsCEN.append(TH1F(filetmp.Get("gen_pt_eta1v3")).GetMean())
    ptErrCEN.append(TH1F(filetmp.Get("gen_pt_eta1v3")).GetMeanError())
   
    ptsFOR.append(TH1F(filetmp.Get("gen_pt_etaUP1v3")).GetMean())
    ptErrFOR.append(TH1F(filetmp.Get("gen_pt_etaUP1v3")).GetMeanError())
    

  filename = "JetResolutionFits"
  f = TFile("%s.root"%filename,  "RECREATE")
  
  l1= TLegend(0.7650754,0.5564767,0.8065327,0.876943)
  l1.SetTextSize(0.035)
  l1.SetLineColor(0)
  l1.SetShadowColor(0)
  l1.SetLineStyle(1)
  l1.SetLineWidth(1)
  l1.SetFillColor(0)
  l1.SetFillStyle(0)
  l1.SetMargin(0.35)

  for j in xrange(0,len(histosCEN)):
    histosCEN[j].Write() 
    histosFOR[j].Write()
    histosCEN[j].Scale(1./histosCEN[j].Integral())
    histosFOR[j].Scale(1./histosFOR[j].Integral())
    histosCEN[j].SetLineColor(col.GetColor(palette[j]))
    histosCEN[j].SetLineWidth(2)
    histosCEN[j].Rebin(1) 
    histosFOR[j].SetLineColor(col.GetColor(palette[j]))
    histosFOR[j].SetLineWidth(2)
    histosFOR[j].Rebin(1)
    l1.AddEntry(histosCEN[j], "M = %i"%masspoints[j],"l")
    
  f.Close()
  
  
  for i in range(0,len(masspoints)):
    print "Central:  Mass = %i  GeV pT = %.4f +/- %.4f GeV" %(masspoints[i], ptsCEN[i],  ptErrCEN[i] )
    print "Central:  Softdrop width = %.4f +/- %.4f GeV   " %(sigmasCentral[i], sigErrCentral[i] )
    print ""
    print "Forward:  Mass = %i  GeV pT = %.4f +/- %.4f GeV" %(masspoints[i], ptsFOR[i],  ptErrFOR[i] )
    print "Forward:  Softdrop width = %.4f +/- %.4f GeV   " %(sigmasForward[i], sigErrForward[i] )
    print "";print ""
  
  
  vxCEN = array("f",ptsCEN)
  vxErrCEN = array("f",ptErrCEN)
  vxFOR = array("f",ptsFOR)
  vxErrFOR = array("f",ptErrFOR)
  vyCEN = array("f",sigmasCentral)
  errCEN = array("f",sigErrCentral)
  vyFOR = array("f",sigmasForward)
  errFOR = array("f",sigErrForward)
  # gCEN = TGraphAsymmErrors(len(vxCEN),vxCEN,vyCEN,vxErrCEN,errCEN)
  # gFOR = TGraphAsymmErrors(len(vxFOR),vxFOR,vyFOR,vxErrFOR,errFOR)

  gCEN = TGraphErrors(len(vxCEN),vxCEN,vyCEN,vxErrCEN,errCEN)
  gFOR = TGraphErrors(len(vxFOR),vxFOR,vyFOR,vxErrFOR,errFOR)

  canv = getCanvas()
  canv.cd()
 
  vFrame = canv.DrawFrame(200,0.0,2200,0.2)
  vFrame.SetYTitle("(#sigma_{reco} - #sigma_{gen})/#sigma_{reco}")
  vFrame.SetXTitle("p_{T} (GeV)")
  vFrame.GetXaxis().SetTitleSize(0.06)
  vFrame.GetXaxis().SetTitleOffset(0.95)
  vFrame.GetXaxis().SetLabelSize(0.05)
  vFrame.GetYaxis().SetTitleSize(0.06)
  vFrame.GetYaxis().SetTitleOffset(1.2)
  vFrame.GetYaxis().SetLabelSize(0.05)
  vFrame.GetXaxis().SetNdivisions(809)
  vFrame.GetYaxis().SetNdivisions(703)
  gCEN.SetMarkerSize(1.6)
  gFOR.SetMarkerSize(1.6)
  gCEN.SetMarkerStyle(20)
  gFOR.SetMarkerStyle(20)
  gCEN.SetMarkerColor(col.GetColor(palette[0]))
  gFOR.SetMarkerColor(col.GetColor(palette[1]))

  gCEN.Draw("PLsame")
  gFOR.Draw("PLsame")
  l.AddEntry(gCEN, " |#eta| #leq 1.3","p")
  l.AddEntry(gFOR, " |#eta| > 1.3","p")
  
  addInfo = TPaveText(0.1959799,0.1632124,0.3580402,0.3717617,"NDC")
  addInfo.AddText("PUPPI softdrop mass")  
  if signal.find("BulkWW") != -1: addInfo.AddText("Bulk G #rightarrow WW")
  elif signal.find("BulkZZ") != -1: addInfo.AddText("Bulk G #rightarrow ZZ")
  elif signal.find("WprimeWZ") != -1: addInfo.AddText("W'#rightarrow WZ")
  elif signal.find("ZprimeWW") != -1: addInfo.AddText("Z' #rightarrow WW")     
  addInfo.SetFillColor(0)
  addInfo.SetLineColor(0)
  addInfo.SetFillStyle(0)
  addInfo.SetBorderSize(0)
  addInfo.SetTextFont(42)
  addInfo.SetTextSize(0.040)
  addInfo.SetTextAlign(12)
  addInfo.AddText("AK, R= 0.8")
  addInfo.AddText("Uncorrected")
  addInfo.AddText("p_{T} > 200 GeV, |#eta| < 2.5")
  
  
  l.Draw("same")
  addInfo.Draw("same")

  CMS_lumi.CMS_lumi(canv, iPeriod, iPos)
  canv.Update()

  
  canvname = "massResolution.pdf"
  canv.SaveAs(canvname,"pdf")
  canv.SaveAs(canvname.replace("pdf","root"),"pdf")
  
  # vyCEN.append(vyCEN[0])
#   errCEN.append(errCEN[0])
#   vyFOR.append(vyFOR[0])
#   errFOR.append(errFOR[0])
#   vxCEN.append(200.)
#   vxErrCEN.append(vxErrCEN[0])
#   vxFOR.append(200.)
#   vxErrFOR.append(vxErrFOR[0])
  

  vyCEN.append(vyCEN[-1])
  errCEN.append(errCEN[-1])  
  vyCEN.append(vyCEN[-1])
  errCEN.append(errCEN[-1])
  vyCEN.append(vyCEN[-1])
  errCEN.append(errCEN[-1])
  
  vyFOR.append(vyFOR[-1])
  errFOR.append(errFOR[-1])
  vyFOR.append(vyFOR[-1])
  errFOR.append(errFOR[-1])
  vyFOR.append(vyFOR[-1])
  errFOR.append(errFOR[-1])
  vyFOR.append(vyFOR[-1])
  errFOR.append(errFOR[-1])
 
  
  vxCEN.append(2500.)
  vxErrCEN.append(vxErrCEN[-1])
  vxCEN.append(3000.)
  vxErrCEN.append(vxErrCEN[-1])
  vxCEN.append(3200.)
  vxErrCEN.append(vxErrCEN[-1])
  
  vxFOR.append(1500.)
  vxErrFOR.append(vxErrFOR[-1])
  vxFOR.append(2000.)
  vxErrFOR.append(vxErrFOR[-1])
  vxFOR.append(3000.)
  vxErrFOR.append(vxErrFOR[-1])
  vxFOR.append(3200.)
  vxErrFOR.append(vxErrFOR[-1])
  
  nvyCEN  = np.array(vyCEN)
  nerrCEN = np.array(errCEN)
  nvyCEN  = -1*(nvyCEN - 1)
  nvyCEN  = 1./nvyCEN
  # nerrCEN = nerrCEN
  # nerrCEN = nvyCEN*0.005 #Try 0.5% error on all points
  # nerrCEN[0] = nvyCEN[0]*0.1
  nnvyCEN = array("f",nvyCEN)
  nnerrCEN = array("f",nerrCEN)
  
  nvyFOR  = np.array(vyFOR)
  nerrFOR = np.array(errFOR)
  nvyFOR  = -1*(nvyFOR-1)
  nvyFOR  = 1./nvyFOR
  # nerrFOR = nvyFOR*0.005 #Try 0.5% error on all points
  nnvyFOR = array("f",nvyFOR)
  nnerrFOR = array("f",nerrFOR)

  gC_forCorr = TGraphErrors(len(vxCEN),vxCEN,nnvyCEN,vxErrCEN,nnerrCEN)
  gC_forCorr.SetName("gCentral")
  gF_forCorr = TGraphErrors(len(vxFOR),vxFOR,nnvyFOR,vxErrFOR,nnerrFOR)
  gF_forCorr.SetName("gForward")
  
  filename = "input/massResolution"
  f = TFile("%s.root"%filename,  "RECREATE")
  gC_forCorr.Write()
  gF_forCorr.Write()
  f.Close()
  time.sleep(105)
  del canv

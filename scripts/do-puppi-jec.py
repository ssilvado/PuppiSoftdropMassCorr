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
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod=4

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-m','--mshift', action='store_true', dest='doMassShiftFit', default=False, help='Fit (reco-gen)/reco instead of pruned mass!')
parser.add_option('-p','--pruning', action='store_true', dest='doPruning', default=False, help='Fit pruned mass')
parser.add_option('--fitGen', action='store_true', dest='fitGenMass', default=False, help='Fit gen mass')

(options, args) = parser.parse_args()

if options.noX: gROOT.SetBatch(True)

prefix = '/mnt/t3nfs01/data01/shome/thaarres/EXOVVAnalysisRunII/AnalysisOutput/80X/GEN/'
gStyle.SetOptFit(1)

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
 colors = ['#762a83','#de77ae','#a6dba0','#92c5de','#4393c3','#2166ac','#053061']

 for c in colors:
  palette['gv'].append(c)
 
 return palette[mode]
 
palette = get_palette('gv')
col = TColor()

masses = [1000,1200,1400,1600,1800,2000,2500,3000,4000,4500]  
masses = [400,600,800,1000,1400,2000,2500,3000,4000,4500] #Ops!! 400 masspoint is named for convenience and is actually SM WW, not signal sample!

hCentral = 'gen_SoftdropMass_eta1v3'
hForward = 'gen_SoftdropMass_etaUP1v3'

if options.doMassShiftFit:
  hCentral = 'massShift_softdrop_eta1v3'
  hForward = 'massShift_softdrop_etaUP1v3'
  
if options.doPruning:
  hCentral = 'massShift_pruning_eta1v3'
  hForward = 'massShift_pruning_etaUP1v3'
  options.doMassShiftFit = True

if options.fitGenMass:
  hCentral = 'GenAK8SoftdropMass_eta1v3'
  hForward = 'GenAK8SoftdropMass_etaUP1v3'
  if options.doPruning:
    hCentral = 'GenAK8PrunedMass_eta1v3'
    hForward = 'GenAK8PrunedMass_etaUP1v3'
    options.doMassShiftFit = False

      
  
lineStyle = [1,1,1,1,3,3,3,3]


signals = ["BulkWW","BulkZZ","ZprimeWW","WprimeWZ"]
signals = ["BulkWW"]


for signal in signals:
 
  filelist = []
  histos = []
  fits = []
  l = TLegend(0.4861809,0.7020725,0.6859296,0.9209845)
  if options.fitGenMass or options.doMassShiftFit: l = TLegend(0.6861809,0.7020725,0.7859296,0.9209845)
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
  for filetmp in filelist:
    #Central
    histtmp = TH1F(filetmp.Get(hCentral))
    histtmp.SetName("Central%s"%filetmp.GetName().replace(".","_"))
    histtmp.Scale(1./histtmp.Integral())
    maxbin=0
    maxcontent=0
    startbin = 60.
    if options.doMassShiftFit: startbin = -2.
    for b in range(histtmp.GetXaxis().GetNbins()):
      if histtmp.GetXaxis().GetBinCenter(b+1) > startbin and histtmp.GetBinContent(b+1)>maxcontent:
        maxbin = b
        maxcontent = histtmp.GetBinContent(b+1)   
    tmpmean = histtmp.GetXaxis().GetBinCenter(maxbin)           
    g1 = TF1("g1CEN%i"%m,"gaus", tmpmean-10.,tmpmean+15.)
    if options.doMassShiftFit: g1 = TF1("g1CEN%s"%filetmp.GetName().replace(".","_"),"gaus", tmpmean-0.5,tmpmean+0.5)
    histtmp.Fit(g1, "SR")
    histos.append(histtmp)
    mean    = g1.GetParameter(1)
    meanerr = g1.GetParError(1)
    meansCentral.append(mean)
    meanErrCentral.append(meanerr)
  
    #Forward
    histtmp = TH1F(filetmp.Get(hForward))
    histtmp.SetName("Forward%i"%m)
    histtmp.Scale(1./histtmp.Integral())
    maxbin=0
    maxcontent=0
    startbin = 60.
    if options.doMassShiftFit: startbin = -2.
    for b in range(histtmp.GetXaxis().GetNbins()):
      if histtmp.GetXaxis().GetBinCenter(b+1)>startbin and histtmp.GetBinContent(b+1)>maxcontent:
        maxbin = b
        maxcontent = histtmp.GetBinContent(b+1)
    tmpmean = histtmp.GetXaxis().GetBinCenter(maxbin)
    g1 = TF1("g1FOR%i"%m,"gaus", tmpmean-10.,tmpmean+15.)
    if options.doMassShiftFit: g1 = TF1("g1CEN%i"%m,"gaus", tmpmean-0.5,tmpmean+0.5)
    histtmp.Fit(g1, "SR")
    histos.append(histtmp)
    mean    = g1.GetParameter(1)
    meanerr = g1.GetParError(1)
    meansForward.append(mean)
    meanErrForward.append(meanerr)
 
    ptsCEN.append(TH1F(filetmp.Get("gen_pt_eta1v3")).GetMean())
    ptErrCEN.append(TH1F(filetmp.Get("gen_pt_eta1v3")).GetMeanError())
   
    ptsFOR.append(TH1F(filetmp.Get("gen_pt_etaUP1v3")).GetMean())
    ptErrFOR.append(TH1F(filetmp.Get("gen_pt_etaUP1v3")).GetMeanError())
    
  
  canv = getCanvas()
  canv.cd()
  histos[0].Draw("same")
  histos[0].SetMaximum(histos[0].GetMaximum()*2)
  for i in range(1, len(histos)):
   histos[i].Draw("same")
  canv.Print("JEC_%s_fit.pdf"%hCentral)
  print signal
  for i in range(0,len(masspoints)):
    print "Central:  Mass = %i  GeV pT = %.2f +/- %.2f GeV" %(masspoints[i], ptsCEN[i],  ptErrCEN[i] )
    print "Central:  Softdrop mass = %.2f +/- %.2f GeV" %(meansCentral[i], meanErrCentral[i] )
    print ""
    print "Forward:  Mass = %i  GeV pT = %.2f +/- %.2f GeV" %(masspoints[i], ptsFOR[i],  ptErrFOR[i] )
    print "Forward:  Softdrop mass = %.2f +/- %.2f GeV" %(meansForward[i], meanErrForward[i] )
    print "";print ""
  
  
  vxCEN = array("f",ptsCEN)
  vxErrCEN = array("f",ptErrCEN)
  vxFOR = array("f",ptsFOR)
  vxErrFOR = array("f",ptErrFOR)
  vyCEN = array("f",meansCentral)
  errCEN = array("f",meanErrCentral)
  vyFOR = array("f",meansForward)
  errFOR = array("f",meanErrForward)
  # gCEN = TGraphAsymmErrors(len(vxCEN),vxCEN,vyCEN,vxErrCEN,errCEN)
  # gFOR = TGraphAsymmErrors(len(vxFOR),vxFOR,vyFOR,vxErrFOR,errFOR)

  gCEN = TGraphErrors(len(vxCEN),vxCEN,vyCEN,vxErrCEN,errCEN)
  gFOR = TGraphErrors(len(vxFOR),vxFOR,vyFOR,vxErrFOR,errFOR)

  canv = getCanvas()
  canv.cd()
  if not options.fitGenMass and not options.doMassShiftFit:
    canv.Divide(1,2,0,0,0)
    canv.cd(1)
    p11_1 = canv.GetPad(1)
    p11_1.SetPad(0.01,0.20,0.99,0.98)
    p11_1.SetRightMargin(0.05)
    p11_1.SetTopMargin(0.05)
    p11_1.SetFillColor(0)
    p11_1.SetBorderMode(0)
    p11_1.SetFrameFillStyle(0)
    p11_1.SetFrameBorderMode(0)
  vFrame = canv.DrawFrame(200,65.,2200,85.)
  vFrame.SetYTitle("<m>_{m_{reco}} (GeV)")
  if options.fitGenMass:
    vFrame = canv.DrawFrame(200,70.,2200,90.)
    vFrame.SetYTitle("<m>_{m_{gen}} (GeV)")
  if options.doMassShiftFit:
    vFrame = canv.DrawFrame(200,-0.3,2200,0.1)
    vFrame.SetYTitle("(m_{reco} - m_{gen})/m_{reco}")
  vFrame.SetXTitle("p_{T} (GeV)")
  vFrame.GetXaxis().SetTitleSize(0.06)
  vFrame.GetXaxis().SetTitleOffset(0.95)
  vFrame.GetXaxis().SetLabelSize(0.05)
  vFrame.GetYaxis().SetTitleSize(0.06)
  vFrame.GetYaxis().SetTitleOffset(1.2)
  vFrame.GetYaxis().SetLabelSize(0.05)
  vFrame.GetXaxis().SetNdivisions(809)
  vFrame.GetYaxis().SetNdivisions(908)
  if options.doMassShiftFit:vFrame.GetYaxis().SetNdivisions(707)
  gCEN.SetMarkerSize(1.6)
  gFOR.SetMarkerSize(1.6)
  gCEN.SetMarkerStyle(20)
  gFOR.SetMarkerStyle(20)
  gCEN.SetMarkerColor(col.GetColor(palette[0]))
  gFOR.SetMarkerColor(col.GetColor(palette[1]))
  filetmp = TFile.Open("/mnt/t3nfs01/data01/shome/thaarres/EXOVVAnalysisRunII/AnalysisOutput/80X/ExoDiBosonAnalysis.W_all.root","READ")
  histtmpCEN = TProfile(filetmp.Get("gen_chsJEC_eta1v3"))
  histtmpFOR = TProfile(filetmp.Get("gen_chsJEC_etaUP1v3"))
  histtmpCEN.SetMarkerSize(1.6)
  histtmpFOR.SetMarkerSize(1.6)
  histtmpCEN.SetMarkerStyle(20)
  histtmpFOR.SetMarkerStyle(20)
  histtmpCEN.SetMarkerColor(col.GetColor(palette[0]))
  histtmpFOR.SetMarkerColor(col.GetColor(palette[1]))
  histtmpCEN.Draw("same")
  histtmpCEN.SetMaximum(80.)
  histtmpCEN.SetMinimum(65.)
  if options.doMassShiftFit:
    histtmpCEN.SetMaximum(1.)
    histtmpCEN.SetMinimum(-1.)
     
  histtmpFOR.Draw("same")
  gCEN.Draw("PLsame")
  gFOR.Draw("PLsame")
  l.AddEntry(gCEN, "|#eta|<1.3","p")
  l.AddEntry(gFOR, "|#eta|>1.3","p")

  l1 = TLatex()
  l1.SetNDC()
  l1.SetTextAlign(12)
  l1.SetTextFont(42)
  l1.SetTextSize(0.035)
  l1.DrawLatex(0.20,0.86, "AK, R= 0.8")
  l1.DrawLatex(0.20,0.82, "No L2L3")
  l1.DrawLatex(0.20,0.78, "p_{T} > 200 GeV, |#eta| < 2.5")
  l1.SetTextSize(0.045)
  if not options.fitGenMass:
    if options.doPruning: 
      l1.DrawLatex(0.20,0.91, "Pruned mass")
    else:
      l1.DrawLatex(0.20,0.91, "PUPPI softdrop mass")  
  else:    
    if options.doPruning:  
      l1.DrawLatex(0.20,0.91, "Gen pruned mass") 
    else:
      l1.DrawLatex(0.20,0.91, "Gen softdrop mass") 
  # if signal.find("BulkWW") != -1: l1.DrawLatex(0.20,0.91, " Bulk G #rightarrow WW")
#   elif signal.find("BulkZZ") != -1: l1.DrawLatex(0.20,0.91, " Bulk G #rightarrow ZZ")
#   elif signal.find("WprimeWZ") != -1: l1.DrawLatex(0.20,0.91, " W'#rightarrow WZ")
#   elif signal.find("ZprimeWW") != -1: l1.DrawLatex(0.20,0.91, " Z' #rightarrow WW")
  
  l.Draw("same")
  # l2.Draw("same")
  if not options.fitGenMass and not options.doMassShiftFit: CMS_lumi.CMS_lumi(p11_1, iPeriod, iPos)
  else: CMS_lumi.CMS_lumi(canv, iPeriod, iPos)
  canv.Update()
  
  if not options.fitGenMass and not options.doMassShiftFit:
    canv.cd(2)
    p11_2 = canv.GetPad(2)
    p11_2.SetPad(0.01,0.02,0.99,0.27)
    p11_2.SetBottomMargin(0.35)
    p11_2.SetRightMargin(0.05)
    p11_2.SetGridx()
    p11_2.SetGridy()
    vFrame2 = p11_2.DrawFrame(p11_1.GetUxmin(),1.04, p11_1.GetUxmax(), 1.09)
    vFrame2.SetXTitle("p_{T} (GeV)")
    vFrame2.SetYTitle("CHS L2L3")
    vFrame2.GetXaxis().SetTitleSize(0.06)
    vFrame2.GetYaxis().SetTitleSize(0.15)
    vFrame2.GetYaxis().SetTitleOffset(0.40)
    vFrame2.GetYaxis().SetLabelSize(0.09)
    vFrame2.GetXaxis().SetTitleSize(0.15)
    # vFrame2.GetXaxis().SetTitleOffset(0.90)
    vFrame2.GetXaxis().SetLabelSize(0.12)
    vFrame2.GetXaxis().SetNdivisions(809)
    vFrame2.GetYaxis().SetNdivisions(403)
  
  
  
  histtmpCEN.Draw("same")
  histtmpFOR.Draw("same")
  canv.Update()
  canvname = "JEC_%s_vspt.pdf"%(hCentral)
  canv.SaveAs(canvname,"pdf")
  canv.SaveAs(canvname.replace("pdf","root"),"pdf")
  
  
  if options.fitGenMass or options.doMassShiftFit:
    
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
    if options.fitGenMass:
      nvyCEN = nvyCEN/80.4
      # nerrCEN = nerrCEN/80.4
      nerrCEN = nvyCEN*0.005 #Try 0.5% error on all points
    elif options.doMassShiftFit:
      nvyCEN  = -1*(nvyCEN - 1)
      # nerrCEN = nerrCEN
      nerrCEN = nvyCEN*0.005 #Try 0.5% error on all points
    nnvyCEN = array("f",nvyCEN)
    nnerrCEN = array("f",nerrCEN)
    
    nvyFOR  = np.array(vyFOR)
    nerrFOR = np.array(errFOR)
    if options.fitGenMass:
      nvyFOR = nvyFOR/80.4
      # nerrFOR = nerrFOR/80.4
      nerrFOR = nvyFOR*0.005 #Try 0.5% error on all points
    elif options.doMassShiftFit:
      nvyFOR  = -1*(nvyFOR-1)
      nerrFOR = nvyFOR*0.005 #Try 0.5% error on all points
    nnvyFOR = array("f",nvyFOR)
    nnerrFOR = array("f",nerrFOR)

    gC_forCorr = TGraphErrors(len(vxCEN),vxCEN,nnvyCEN,vxErrCEN,nnerrCEN)
    gC_forCorr.SetName("gC_forCorr")
    gF_forCorr = TGraphErrors(len(vxFOR),vxFOR,nnvyFOR,vxErrFOR,nnerrFOR)
    gF_forCorr.SetName("gF_forCorr")
    
    filename = "input/genCorr"
    if options.doMassShiftFit: filename = "input/recoCorr"
    f = TFile("%s.root"%filename,  "RECREATE")
    gC_forCorr.Write()
    gF_forCorr.Write()
    f.Close()

    
  time.sleep(105)

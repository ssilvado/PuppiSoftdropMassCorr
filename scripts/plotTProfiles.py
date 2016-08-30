import time
import CMS_lumi, tdrstyle
from ROOT import *

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = ""
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod=4

rebin = 1

setYmax = 90.
setYmin = 70.
prefix = '/shome/thaarres/EXOVVAnalysisRunII/AnalysisOutput/80X/'
doFit = False


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
 colors = ['#A43820','#de77ae','#a6dba0','#92c5de','#4393c3','#2166ac','#053061']

 for c in colors:
  palette['gv'].append(c)
 
 return palette[mode]
 
palette = get_palette('gv')

col = TColor()
  
file = 'BulkWW_all.root'



histonames = ['PUPPISDvsNPV','PUPPISDvsETA']
legendname = ['Corrected PUPPI+softdrop','Corrected PUPPI+softdrop']
titles     = ['Number of primary vertices','#eta']
lineStyle = [1,1,1,1,3,3,3,3]
markerStyle = [20,24,22,26,33]

filelist = []



filename = prefix + file
filetmp = TFile.Open(filename,"READ") 
    
ii = -1
for hname in histonames:
  histos = []
  ii += 1
  
  min = 0
  max = 52
  
  if hname.find("ETA") != -1:
    min = -2.52
    max = 2.52
  

  l = TLegend(0.4535176,0.734456,0.5954774,0.8419689)
  l.SetTextSize(0.035)
  l.SetLineColor(0)
  l.SetShadowColor(0)
  l.SetLineStyle(1)
  l.SetLineWidth(1)
  l.SetFillColor(0)
  l.SetFillStyle(0)
  l.SetMargin(0.35)
  

  histtmp = TProfile(filetmp.Get(hname))
  histtmp.SetName(filetmp.GetName())
  histos.append(histtmp)

  for j in xrange(0,len(histos)):
    # histos[j].SetLineColor(col.GetColor(palette[j %4]))
    # histos[j].SetLineStyle(lineStyle[j])
    # histos[j].SetLineWidth(3)
    histos[j].Rebin(rebin)
    histos[j].SetMarkerStyle(markerStyle[j])
    histos[j].SetMarkerColor(col.GetColor(palette[j %4]))
    histos[j].SetLineColor(col.GetColor(palette[j %4]))
    legend = legendname[j]
    l.AddEntry(histos[j],legend,"p")


  fits = []
  for h in histos:
    fittmp = TGraph(h)
    fits.append(fittmp)


  yTitle = "PUPPI softdrop mass (GeV)"

  canv = getCanvas()
  canv.SetName("canv%i"%ii)
  canv.cd()
  vFrame = canv.DrawFrame(min,setYmin,max,setYmax)  
  vFrame.SetXTitle(titles[ii])
  vFrame.SetYTitle(yTitle)
  vFrame.GetXaxis().SetTitleSize(0.06)
  vFrame.GetXaxis().SetTitleOffset(0.95)
  vFrame.GetXaxis().SetLabelSize(0.05)
  vFrame.GetYaxis().SetTitleSize(0.06)
  #vFrame.GetYaxis().SetTitleOffset(1.0)
  vFrame.GetYaxis().SetLabelSize(0.05)
  vFrame.GetXaxis().SetNdivisions(408)
  vFrame.GetYaxis().SetNdivisions(404)


  if doFit:
    for f in fits: f.Draw("Csame")      
  else:
    for h in histos: h.Draw("sameEMP")



  addInfo = TPaveText(0.2085427,0.1761658,0.3454774,0.384715,"NDC")
  addInfo.AddText("PUPPI softdrop mass")  
  if filename.find("BulkWW") != -1: addInfo.AddText("Bulk G #rightarrow WW")
  elif filename.find("BulkZZ") != -1: addInfo.AddText("Bulk G #rightarrow ZZ")
  elif filename.find("WprimeWZ") != -1: addInfo.AddText("W'#rightarrow WZ")
  elif filename.find("ZprimeWW") != -1: addInfo.AddText("Z' #rightarrow WW")     
  addInfo.SetFillColor(0)
  addInfo.SetLineColor(0)
  addInfo.SetFillStyle(0)
  addInfo.SetBorderSize(0)
  addInfo.SetTextFont(42)
  addInfo.SetTextSize(0.040)
  addInfo.SetTextAlign(12)
  addInfo.AddText("AK, R= 0.8")
  addInfo.AddText("Mass corrections applied")
  addInfo.AddText("p_{T} > 200 GeV, |#eta| < 2.5")
  addInfo.Draw("same")
  

  l.Draw("same")
  CMS_lumi.CMS_lumi(canv, iPeriod, iPos)
  canv.Update()
  canvname = "%s.pdf"%hname
  canv.SaveAs(canvname,"pdf")
  canv.SaveAs(canvname.replace("pdf","root"),"pdf")
time.sleep(100)

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

rebin = 5

setmax = 0.39
prefix = '/mnt/t3nfs01/data01/shome/thaarres/EXOVVAnalysisRunII/AnalysisOutput/80X/'
addBKG = False
doFit = True


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
 colors = ['#FF420E','#80BD9E','#336B87','#763626','#003B46','#66A5AD']

 for c in colors:
  palette['gv'].append(c)
 
 return palette[mode]
 
palette = get_palette('gv')

col = TColor()
  
files = [ 'ExoDiBosonAnalysis.BulkWW_13TeV_1000GeV.VV.root',
          'ExoDiBosonAnalysis.BulkZZ_13TeV_1000GeV.VV.root',
          'ExoDiBosonAnalysis.ZprimeWW_13TeV_1000GeV.VV.root',
          'ExoDiBosonAnalysis.WprimeWZ_13TeV_1000GeV.VV.root',
          'ExoDiBosonAnalysis.BulkWW_13TeV_4000GeV.VV.root',
          'ExoDiBosonAnalysis.BulkZZ_13TeV_4000GeV.VV.root',
          'ExoDiBosonAnalysis.ZprimeWW_13TeV_4000GeV.VV.root',
          'ExoDiBosonAnalysis.WprimeWZ_13TeV_4000GeV.VV.root'
          ]


legendname = [ "W (G_{Bulk}#rightarrow WW (Madgraph))",
               "Z (G_{Bulk}#rightarrow ZZ (Madgraph))",
               "W (Z'#rightarrow WW (Madgraph))",       
               "Z (W#rightarrow WZ (Madgraph))",     
             ]
                

histonames = ['gen_SoftdropMass_NEWCORR','gen_SoftdropMass_eta1v3_NEWCORR','gen_SoftdropMass_etaUP1v3_NEWCORR']
titles     = ['PUPPI softdrop mass (GeV)','PUPPI softdrop mass (GeV)','PUPPI softdrop mass (GeV)']


lineStyle = [1,1,1,1,3,3,3,3]



ii = -1
for hname in histonames:
  print "Working on histogram " ,hname
  ii += 1
  
  min = 40
  max = 135
  
  if hname.find("gen_Mass") != -1:
    min = 0
    max = 305
  
  if hname.find("Tau21") != -1 or hname.find("DDT") != -1 :
    min = 0
    max = 1  

  #l = TLegend(.16,.7,.4,.9)
  l = TLegend(0.4861809,0.7020725,0.6859296,0.9209845)
  l.SetTextSize(0.035)
  l.SetLineColor(0)
  l.SetShadowColor(0)
  l.SetLineStyle(1)
  l.SetLineWidth(1)
  l.SetFillColor(0)
  l.SetFillStyle(0)
  l.SetMargin(0.35)

  l2 = TLegend(0.7613065,0.6088083,0.9937186,0.6761658)
  l2.SetTextSize(0.033)
  l2.SetLineColor(0)
  l2.SetShadowColor(0)
  l2.SetLineStyle(1)
  l2.SetLineWidth(1)
  l2.SetFillColor(0)
  l2.SetFillStyle(0)
  l2.SetMargin(0.35)

  filelist = []
  histos = []

  for f in files:
    filename = prefix + f
    filetmp = TFile.Open(filename,"READ") 
    filetmp.SetName(filetmp.GetName())
    filelist.append(filetmp)
    
  for f in filelist:
      histtmp = TH1F(f.Get(hname))
      histtmp.SetName(f.GetName())
      histos.append(histtmp)

  for j in xrange(0,len(histos)):
    print "Working on file " ,histos[j].GetName()
    histos[j].SetLineColor(col.GetColor(palette[j %4]))
    histos[j].SetLineStyle(lineStyle[j])
    histos[j].SetLineWidth(3)
    histos[j].Rebin(rebin)
    histos[j].Scale(1./histos[j].Integral())
    if j < 4:
      legend = legendname[j]
      l.AddEntry(histos[j],legend,"l")
    if j == 0: l2.AddEntry(histos[j],"1.0 TeV","l")
    if j == 4: l2.AddEntry(histos[j],"4.0 TeV","l")
  
  fits = []
  for h in histos:
    fittmp = TGraph(h)
    fits.append(fittmp)
  
  if addBKG:
    bkgfile = TFile.Open('80X/ExoDiBosonAnalysis.QCD_pythia8_cp.root','READ')
    bkghist = TH1F(bkgfile.Get(hname))
    bkghist.SetName("BKG")
    bkghist.SetFillStyle(3002)
    bkghist.SetFillColor(kBlack)
    bkghist.SetLineColor(kBlack)
    bkghist.SetLineStyle(1)
    bkghist.Scale(1./bkghist.Integral())
    bkghist.Rebin(rebin)
    bkghist.SetMaximum(setmax)
    l.AddEntry(bkghist,'QCD (Pythia8)',"f")

  yTitle = "Arbitrary scale"
   
  canv = getCanvas()
  canv.cd()
  setmax = histos[0].GetMaximum()*1.5
  vFrame = canv.DrawFrame(min,0.000005,max,setmax)  
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
    if addBKG: bkghist.Draw('HISTsame')
  else:
    for h in histos: h.Draw("sameHIST")
    if addBKG: bkghist.Draw('HISTsame')
  


  l1 = TLatex()
  l1.SetNDC()
  l1.SetTextAlign(12)
  l1.SetTextFont(42)
  l1.SetTextSize(0.025)
  l1.DrawLatex(0.20,0.80, "AK, R= 0.8")
  if hname.find("pt500") != -1: 
    l1.DrawLatex(0.20,0.74, "p_{T} > 500 GeV, |#eta| < 2.5")
  elif hname.find("eta1v3") != -1: 
    l1.DrawLatex(0.20,0.74, "p_{T} > 200 GeV, |#eta| < 1.3")  
  elif hname.find("etaUP1v3") != -1: 
    l1.DrawLatex(0.20,0.74, "p_{T} > 200 GeV, 1.3 < |#eta| < 2.5")  
  else:
       l1.DrawLatex(0.20,0.74, "p_{T} > 200 GeV, |#eta| < 2.5")
  
  
  if hname.find("CHSCorr") != -1:
    l1.DrawLatex(0.20,0.77, "L2L3 CHS corrected")
  elif hname.find("NEWCORR") != -1:
    l1.DrawLatex(0.20,0.77, "Mass corrections applied")
  elif hname.find("Mass") != -1 and not hname.find("UnCorr") != -1 and not hname.find("1v3") != -1  and not hname.find("Tau21") != -1:
    l1.DrawLatex(0.20,0.77, "L2L3 corrected")
  
  # if not hname.find("Mass") != -1: l1.DrawLatex(0.20,0.76, "65 GeV < M_{G} < 105 GeV")
  l1.DrawLatex(0.20,0.71, "65 GeV < M_{G} < 105 GeV")

  # l1.DrawLatex(0.7,0.42,"65 GeV < M_{p} < 105 GeV")
  

   
  l.Draw("same")
  l2.Draw("same")
  CMS_lumi.CMS_lumi(canv, iPeriod, iPos)
  canv.Update()
  canvname = "%s.pdf"%hname.replace("gen_","")
  canv.SaveAs(canvname,"pdf")
  canv.SaveAs(canvname.replace("pdf","root"),"pdf")
  time.sleep(20)

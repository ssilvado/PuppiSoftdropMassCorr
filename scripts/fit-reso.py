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
parser.add_option('-r','--reco', action='store_true', dest='calculateRecoJEC', default=False, help='calculate reco JEC')

(options, args) = parser.parse_args()

if options.noX: 
  gROOT.SetBatch(True)

gStyle.SetOptFit(0)

def getCanvas():
  c = TCanvas("c","c",800,800)
  c.GetWindowHeight()
  c.GetWindowWidth()
  c.SetTitle("")
  return c
  
def get_palette(mode):
 palette = {}
 palette['gv'] = [] 
 colors = ['#FF420E','#80BD9E','#a6dba0','#92c5de','#4393c3','#2166ac','#053061','#40004b','#762a83','#9970ab','#de77ae']
 for c in colors:
  palette['gv'].append(c)
 return palette[mode]
 
palette = get_palette('gv')
col = TColor()

filetmp = TFile.Open("input/massResolution.root","READ")   
gCentral = filetmp.Get("gCentral")
gForward = filetmp.Get("gForward")
  
canv = getCanvas()
canv.cd()
canv.SetName("fitRECO")
vFrame = canv.DrawFrame(200,1.05,3100,1.2)
vFrame.SetYTitle("#sigma_{reco} / #sigma_{gen}")


vFrame.SetXTitle("p_{T} (GeV)")
vFrame.GetXaxis().SetTitleSize(0.06)
vFrame.GetXaxis().SetTitleOffset(0.95)
vFrame.GetXaxis().SetLabelSize(0.05)
vFrame.GetYaxis().SetTitleSize(0.06)
vFrame.GetYaxis().SetTitleOffset(1.2)
vFrame.GetYaxis().SetLabelSize(0.05)
vFrame.GetXaxis().SetNdivisions(809)
vFrame.GetYaxis().SetNdivisions(703)
gCentral.SetMarkerSize(1.6)
gForward.SetMarkerSize(1.6)
gCentral.SetMarkerStyle(20)
gForward.SetMarkerStyle(20)
gCentral.SetMarkerColor(col.GetColor(palette[0]))
gForward.SetMarkerColor(col.GetColor(palette[1]))
gCentral.SetLineColor(col.GetColor(palette[0]))
gForward.SetLineColor(col.GetColor(palette[1]))
l = TLegend(0.60461809,0.7620725,0.7559296,0.9009845)
l.SetTextSize(0.035)
l.SetLineColor(0)
l.SetShadowColor(0)
l.SetLineStyle(1)
l.SetLineWidth(1)
l.SetFillColor(0)
l.SetFillStyle(0)
l.SetMargin(0.35)

fitmin = TMath.MinElement(gCentral.GetN(),gCentral.GetX())
fitmin = 200.
fitmax = TMath.MaxElement(gCentral.GetN(),gCentral.GetX())
print "Fitting range: [%.2f,%.2f] "%(fitmin,fitmax)

l.AddEntry(gForward, "Mass resol. |#eta| > 1.3","p")
l.AddEntry(gCentral, "Mass resol. |#eta| #leq 1.3","p")

# g = TF1("massResolution_0eta1v3","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",fitmin,fitmax)
g = TF1("massResolution_0eta1v3","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",fitmin,fitmax)
# g = TF1("gJEC","pol4",fitmin,fitmax)
# g.SetParameter(0,      1.21228)
#  g.SetParameter(1, -0.000559701)
#  g.SetParameter(2,  8.96666e-07)
#  g.SetParameter(3, -5.44342e-10)
#  g.SetParameter(4,   1.1131e-13)
#  g.SetParameter(5, -1.95118e-17)

g.SetParameter(0, 1.08382e+00)
g.SetParameter(1, 5.75471e-05)
g.SetParameter(2,-1.46527e-07)
g.SetParameter(3, 1.23074e-10)
g.SetParameter(4,-4.12487e-14)
g.SetParameter(5, 4.82147e-18)
g.SetParameter(6, 1.08382e+00)

print "Fitting for central region" 
for i in range(0,10):
  gCentral.Fit(g, "EX0FRU")

gCentral.Draw("PEsame") 

# print "" ; print "" ; print "" ;
# sC = TSpline3("SPLINE_massResolution_0eta1v3",gForward,"",fitmin,fitmax)
# sC.SetLineColor(kBlue)
# sC.Draw("same")

# g2 = TF1("gJEC2","pol3",fitmin,fitmax)
g2 = TF1("massResolution_1v3eta2v5","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",fitmin,fitmax)
# g2 = TF1("massResolution_1v3eta2v5","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)",fitmin,fitmax)
g2.SetParameter(0,  1.25551e+00)
g2.SetParameter(1, -6.01398e-04)
g2.SetParameter(2,  9.04900e-07)
g2.SetParameter(3, -5.92240e-10)
g2.SetParameter(4,  1.73920e-13)
g2.SetParameter(5, -1.87967e-17)
g2.SetParameter(6,  1.25551e+00)


print "Fitting for forward region" 
for i in range(0,10):
  gForward.Fit(g2, "EX0FRU")


gForward.Draw("PEsame")


# sF = TSpline3("SPLINE_massResolution_1v3eta2v5",gForward,"",fitmin,fitmax)
# sF.SetLineColor(kBlue)
# s2.Draw("same")
# gForward.Draw("PEsame")
massReso_central = gCentral.GetFunction("massResolution_0eta1v3")
massReso_forward = gForward.GetFunction("massResolution_1v3eta2v5")

filename = "weights/puppiSoftdropResol"
f = TFile("%s.root"%filename,  "UPDATE")
print "Writing to file " ,f.GetName()
massReso_central.Write("",TObject.kOverwrite)
massReso_forward.Write("",TObject.kOverwrite)
# sC.Write("",TObject.kOverwrite)
# sF.Write("",TObject.kOverwrite)
canv.Write("",TObject.kOverwrite)
f.Close()

l.Draw("same")
CMS_lumi.CMS_lumi(canv, iPeriod, iPos)
canv.RedrawAxis()
canv.Update()
canvname = "puppiSoftdropResol_fit.pdf"
canv.SaveAs(canvname,"pdf")
canv.SaveAs(canvname.replace("pdf","root"),"pdf")
time.sleep(205)



import sys, math
from optparse import OptionParser, OptionGroup

parser = OptionParser()
parser.add_option("-f", "--FCSensitivity", action="store_false", dest="use5SigDL", default=False, help="Use Feldman-Cousins sensitivity (default)")
parser.add_option("-d", "--5SigmaDL", action="store_true", dest="use5SigDL", help="Use 5-sigma discovery level")
parser.add_option("-q", "--QRPA", action="store_false", dest="useShellModel", default=False, help="Use QRPA NME (default)")
parser.add_option("-s", "--SM", action="store_true", dest="useShellModel", help="Use shell model NME (except for 150Nd)")
parser.add_option("-g", "--set_ga", type="float", dest="quenchGa", default=1.25, help="change g_a from 1.25 to VAL", metavar="VAL")
parser.add_option("-k", "--drawKKDC", action="store_true", dest="drawKKDC", default=False, help="draw KKDC region")
parser.add_option("-m", "--drawMbb", action="store_true", dest="drawMbb", default=False, help="draw Mbb vs exposure")
parser.add_option("-B", "--set_highBG", type="float", dest="highBG", default=10.0, help="set high BG to VAL (default = 4)", metavar="VAL")
parser.add_option("-x", "--xenon-136", action="store_true", dest="xenon136", default=False, help="draw for 136Xe")
parser.add_option("-e", "--cut-efficiency", type="float", dest="cutefficiency", default=0.9*0.92, help="set efficiency of cuts to VAL (default = 0.9*0.92 (MJD PSA*FV))", metavar="VAL")
(options, args) = parser.parse_args()

expMin_ty = 1e-3
expMax_ty = 1000
mbbMin_meV = .99999
mbbMax_meV = 1000
tHalfMin_y = 0.9999e24
tHalfMax_y = 1.e30

from ROOT import gROOT, TCanvas, gStyle, TLegend, TBox, TColor, TLatex, TGaxis, TAxis

gROOT.ProcessLine(".include \"$CLHEP_INCLUDE_DIR\"")
gROOT.ProcessLine(".L SensClass.C+")

from ROOT import SensClass, CLHEP

gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasBorderMode(0)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
lineWidth = 2

canvas = TCanvas("canvas", "canvas", 1024, 768)
canvas.SetFillColor(0)
canvas.SetLogy()
canvas.SetLogx()

statModel = SensClass.kFC90CLSens
if options.use5SigDL: statModel = SensClass.k5SigmaDL
nmeModel = SensClass.kQRPA
if options.useShellModel: nmeModel = SensClass.kShellModel
ga = options.quenchGa
drawKKDC = False
if options.drawKKDC: drawKKDC = True
drawMbb = False
if options.drawMbb: drawMbb = True
highBG = options.highBG
if options.xenon136: 
    isotope = SensClass.kXe136
    enrichment = 0.9
else:
    isotope = SensClass.kGe76
    enrichment = 0.86
efficiency = SensClass.GetROIEfficiency()*options.cutefficiency # ROI x PSA x FV
cPerROIty = 1./CLHEP.tonne/CLHEP.year

bgFreeSens = SensClass(isotope, enrichment, efficiency, 0.*cPerROIty, statModel, nmeModel, ga)
if drawMbb: bgFreeFunc = bgFreeSens.GetMbbVsExposureFunction()
else: bgFreeFunc = bgFreeSens.GetTHalfVsExposureFunction(True, SensClass.kYears)
bgFreeFunc.SetRange(expMin_ty, expMax_ty)

lowBGSens = SensClass(isotope, enrichment, efficiency, 0.1*cPerROIty, statModel, nmeModel, ga)
if drawMbb: lowBGFunc = lowBGSens.GetMbbVsExposureFunction()
else: lowBGFunc = lowBGSens.GetTHalfVsExposureFunction(True, SensClass.kYears)
lowBGFunc.SetRange(expMin_ty, expMax_ty)
lowBGFunc.SetLineStyle(9)
lowBGFunc.SetLineWidth(lineWidth)

tonScaleSens = SensClass(isotope, enrichment, efficiency, 1.*cPerROIty, statModel, nmeModel, ga)
if drawMbb: tonScaleFunc = tonScaleSens.GetMbbVsExposureFunction()
else: tonScaleFunc = tonScaleSens.GetTHalfVsExposureFunction(True, SensClass.kYears)
tonScaleFunc.SetRange(expMin_ty, expMax_ty)
tonScaleFunc.SetLineStyle(7)
tonScaleFunc.SetLineWidth(lineWidth)

hiBGSens = SensClass(isotope, enrichment, efficiency, highBG*cPerROIty, statModel, nmeModel, ga)
if drawMbb: hiBGFunc = hiBGSens.GetMbbVsExposureFunction()
else: hiBGFunc = hiBGSens.GetTHalfVsExposureFunction(True, SensClass.kYears)
hiBGFunc.SetLineStyle(10)
hiBGFunc.SetLineWidth(lineWidth)

hiBGFunc.SetRange(expMin_ty, expMax_ty)
if drawMbb: 
    hiBGFunc.SetMaximum(mbbMax_meV)
    hiBGFunc.SetMinimum(mbbMin_meV)
else:
    hiBGFunc.SetMaximum(tHalfMax_y)
    hiBGFunc.SetMinimum(tHalfMin_y)
hiBGFunc.GetXaxis().SetTitle("Exposure [ton-years]")
hiBGFunc.GetXaxis().CenterTitle()
hiBGFunc.GetXaxis().SetTitleOffset(1.2)
if drawMbb: 
    yTitle = "m_{#beta#beta} "
    if options.use5SigDL: yTitle += "5#sigma DL ("
    else: yTitle += "sensitivity (90% CL, "
    if options.useShellModel: yTitle += "SM, "
    else: yTitle += "QRPA, "
    yTitle += "g_{A}=" + str(ga) + ") [meV]"
else: 
    if options.xenon136: yTitle = "^{136}Xe T_{1/2} "
    else: yTitle = "^{76}Ge T_{1/2} "
    if options.use5SigDL: yTitle += "5#sigma DL "
    else: yTitle += "sensitivity (90% CL) "
    yTitle += "[years]"
hiBGFunc.GetYaxis().SetTitle(yTitle)
hiBGFunc.GetYaxis().CenterTitle()
#hiBGFunc.GetYaxis().SetNoExponent()
if not drawMbb: hiBGFunc.GetYaxis().SetTitleOffset(1.2)
hiBGFunc.Draw()

# eyeballing from Schubert plot m_l -> 0, with PDG 2013 osc parameters
ihMinMinMbb_meV = 16.2
ihMinMbb_meV = 18.3
ihMaxMbb_meV = 48.3
ihMaxMaxMbb_meV = 49.2
ihText = "Inverted Hierarchy (m_{#it{l}} #rightarrow 0 eV) "
if drawMbb:
    bInvHierOut = TBox(expMin_ty, ihMinMinMbb_meV, 0.98*expMax_ty, ihMaxMaxMbb_meV)
    bInvHierIn = TBox(expMin_ty, ihMinMbb_meV, 0.98*expMax_ty, ihMaxMbb_meV)
    lInvHier = TLatex(2.*expMin_ty, math.sqrt(ihMinMinMbb_meV*ihMaxMaxMbb_meV), ihText)
else:
    ihMinMinTHalf_y = SensClass.GetTHalf(ihMaxMaxMbb_meV*CLHEP.meV, isotope, nmeModel, ga)/CLHEP.year
    ihMinTHalf_y = SensClass.GetTHalf(ihMaxMbb_meV*CLHEP.meV, isotope, nmeModel, ga)/CLHEP.year
    ihMaxTHalf_y = SensClass.GetTHalf(ihMinMbb_meV*CLHEP.meV, isotope, nmeModel, ga)/CLHEP.year
    ihMaxMaxTHalf_y = SensClass.GetTHalf(ihMinMinMbb_meV*CLHEP.meV, isotope, nmeModel, ga)/CLHEP.year
    bInvHierOut = TBox(expMin_ty, ihMinMinTHalf_y, 0.98*expMax_ty, ihMaxMaxTHalf_y)
    bInvHierIn = TBox(expMin_ty, ihMinTHalf_y, 0.98*expMax_ty, ihMaxTHalf_y)
    if options.useShellModel: ihText += "(SM"
    else: ihText += "(QRPA"
    ihText += ", g_{A}=" + str(ga) + ")"
    lInvHier = TLatex(2.*expMin_ty, math.sqrt(ihMinTHalf_y*ihMaxTHalf_y), ihText)
bInvHierOut.SetFillColor(TColor.GetColor(0.8, 1.0, 0.8))
bInvHierOut.Draw()
bInvHierIn.SetFillColor(TColor.GetColor(0.7, 0.95, 0.7))
bInvHierIn.Draw()
lInvHier.SetTextAlign(12)
lInvHier.SetTextColor(TColor.GetColor(0.1, 0.4, 0.1))
lInvHier.SetTextSize(0.04)
lInvHier.Draw()

if drawKKDC:
    hvkkTHalfLo_y = (2.23 - 3.*0.31)*1e25
    hvkkTHalfHi_y = (2.23 + 3.*0.44)*1e25
    hvkkMbbLo_meV = SensClass.GetMbb(hvkkTHalfHi_y*CLHEP.year, isotope, nmeModel, ga)/(1.e-3*CLHEP.eV)
    hvkkMbbHi_meV = SensClass.GetMbb(hvkkTHalfLo_y*CLHEP.year, isotope, nmeModel, ga)/(1.e-3*CLHEP.eV)
    hvkkBox = TBox(expMin_ty, hvkkMbbLo_meV, 0.98*expMax_ty, hvkkMbbHi_meV)
    hvkkBox.SetFillColor(TColor.GetColor(0.7, 0.7, 0.95))
    hvkkBox.Draw()
    hvkkLabelText = "Mod. Phys. Lett. A 21 (2006), p. 1547 (3#sigma): (%0.2f-%0.2f) x 10^{25} years" % (hvkkTHalfLo_y/1e25, hvkkTHalfHi_y/1e25)
    lHVKK = TLatex(0.002, math.pow(10, (math.log10(hvkkMbbLo_meV)+math.log10(hvkkMbbHi_meV))/2), hvkkLabelText)
    lHVKK.SetTextAlign(12)
    lHVKK.SetTextSize(0.03)
    lHVKK.SetTextColor(TColor.GetColor(0.1, 0.1, 0.4))
    lHVKK.Draw()

hiBGFunc.GetHistogram().Draw("AXIS SAME")
hiBGFunc.Draw("SAME")
lowBGFunc.Draw("SAME")
tonScaleFunc.Draw("SAME")
bgFreeFunc.Draw("SAME")

canvas.Update()

if drawMbb:
    if not options.useShellModel: legend = TLegend(0.63, 0.73, 0.94, 0.88)
    else: legend = TLegend(0.63, 0.53, 0.94, 0.68)
else: legend = TLegend(0.63, 0.18, 0.94, 0.33)
legend.AddEntry(bgFreeFunc, "Background free", "L")
legend.AddEntry(lowBGFunc, "0.1 counts/ROI/t/y", "L")
legend.AddEntry(tonScaleFunc, "1.0 count/ROI/t/y", "L")
legend.AddEntry(hiBGFunc, str(highBG)+" counts/ROI/t/y", "L")
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.Draw()

if drawMbb:
    tHalfLow_y = SensClass.GetTHalf(mbbMax_meV*1.e-3*CLHEP.eV, isotope, nmeModel, ga)/CLHEP.year
    tHalfHigh_y = SensClass.GetTHalf(mbbMin_meV*1.e-3*CLHEP.eV, isotope, nmeModel, ga)/CLHEP.year
    nLabels = 4
    gAxis = TGaxis(expMax_ty*1.01, mbbMax_meV, expMax_ty - (expMax_ty-expMin_ty)/10000., mbbMin_meV, tHalfLow_y, tHalfHigh_y, nLabels, "-G")
    if(options.xenon136): gAxis.SetTitle("^{136}Xe T^{0#nu}_{1/2} sensitivity (90% CL) [years]")
    else: gAxis.SetTitle("^{76}Ge T^{0#nu}_{1/2} sensitivity (90% CL) [years]")
    gAxis.SetTitleOffset(1.4)
    gAxis.SetBit(TAxis.kRotateTitle)
    gAxis.SetBit(TAxis.kCenterTitle)
    gAxis.Draw()
else:
    mbbHi_meV = SensClass.GetMbb(tHalfMin_y*CLHEP.year, isotope, nmeModel, ga)/CLHEP.meV
    mbbLow_meV = SensClass.GetMbb(tHalfMax_y*CLHEP.year, isotope, nmeModel, ga)/CLHEP.meV
    nLabels = 4
    gAxis = TGaxis(expMax_ty*1.01, tHalfMax_y, expMax_ty - (expMax_ty-expMin_ty)/10000., tHalfMin_y, mbbLow_meV, mbbHi_meV, nLabels, "-G")
    gAxis.SetTitle("m_{#beta#beta} sensitivity (90% CL) [meV]")
    gAxis.SetTitleOffset(1.2)
    gAxis.SetTitleFont(hiBGFunc.GetYaxis().GetTitleFont())
    gAxis.SetTitleSize(hiBGFunc.GetYaxis().GetTitleSize())
    gAxis.SetLabelFont(hiBGFunc.GetYaxis().GetLabelFont())
    gAxis.SetLabelSize(hiBGFunc.GetYaxis().GetLabelSize())
    gAxis.SetBit(TAxis.kRotateTitle)
    gAxis.SetBit(TAxis.kCenterTitle)
    gAxis.Draw()

canvas.Update()
label = "sensVsExp"
if options.use5SigDL: label = "DLVsExp"
if options.xenon136: label += "_Xe"
if options.useShellModel: label += "_SM"
if ga != 1.25: label += "_g"
canvas.Print(label+".pdf")
canvas.Print(label+".png")


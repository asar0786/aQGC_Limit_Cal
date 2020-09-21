import ROOT as r
import sys
import Utils

fileIn = r.TFile("ch1_splitted_TF1_WV.root","read");
#fileIn = r.TFile("Check_bkg_bug.root","read");
#fileIn = r.TFile("Check_bkg_bug_2500.root","read");
#r.gStyle.SetErrorX(0.0001);

c1 = r.TCanvas("c1");
c1, pad1, pad2 = Utils.createCanvasPads();
legend = r.TLegend(0.5,0.7,0.9,0.9);

NominalHist =fileIn.Get("W1+jets");
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("shape_W+jetsUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lDown = fileIn.Get("shape_W+jetsDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

legend.AddEntry(NominalHist,"nominial","lpe");
legend.AddEntry(lUp,"+1#sigma shift","lp");
legend.AddEntry(lDown,"-1#sigma shift","lp");
legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

for i in range(1,5):
  hratio1.SetBinError(i,0.001)
  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("WjetFitSyst_Par0_Wjets.png")
c1.SaveAs("WjetFitSyst_Par0_Wjets.pdf")
c1.SetLogy(1)
c1.SaveAs("WjetFitSyst_Par0_Wjets_log.png")
c1.SaveAs("WjetFitSyst_Par0_Wjets_log.pdf")
c1.SetLogy(0)
c1.Clear()
c1, pad1, pad2 = Utils.createCanvasPads();


lUp = fileIn.Get("shape2_W+jetsUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lDown = fileIn.Get("shape2_W+jetsDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2, 0.4, 1.6)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4, 0.4, 1.6)

for i in range(1,5):
  hratio1.SetBinError(i,0.001)
  hratio2.SetBinError(i,0.001)
hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("WjetFitSyst_Par1_Wjets.png")
c1.SaveAs("WjetFitSyst_Par1_Wjets.pdf")
c1.SetLogy(1)
c1.SaveAs("WjetFitSyst_Par1_Wjets_log.png")
c1.SaveAs("WjetFitSyst_Par1_Wjets_log.pdf")
c1.SetLogy(0)
c1.Clear()
c1, pad1, pad2 = Utils.createCanvasPads();

lUp = fileIn.Get("shape3_W+jetsUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lDown = fileIn.Get("shape3_W+jetsDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2, 0.5, 1.5)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4, 0.5, 1.5)
for i in range(1,5):
  hratio1.SetBinError(i,0.001)
  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("WjetFitSyst_AlternateShape_Wjets.png")
c1.SaveAs("WjetFitSyst_AlternateShape_Wjets.pdf")
c1.SetLogy(1)
c1.SaveAs("WjetFitSyst_AlternateShape_Wjets_log.png")
c1.SaveAs("WjetFitSyst_AlternateShape_Wjets_log.pdf")
c1.SetLogy(0)
c1.Clear()
####################################################################
c1, pad1, pad2 = Utils.createCanvasPads();

lUp = fileIn.Get("shape4_W+jetsUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lDown = fileIn.Get("shape4_W+jetsDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2, 0.5, 1.5)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4, 0.5, 1.5)
for i in range(1,5):
  hratio1.SetBinError(i,0.001)
  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("AlphaFitSyst_Par0_Wjets.png")
c1.SaveAs("AlphaFitSyst_Par0_Wjets.pdf")
c1.SetLogy(1)
c1.SaveAs("AlphaFitSyst_Par0_Wjets_log.png")
c1.SaveAs("AlphaFitSyst_Par0_Wjets_log.pdf")
c1.SetLogy(0)
c1.Clear()
####################################################################
c1, pad1, pad2 = Utils.createCanvasPads();

lUp = fileIn.Get("shape5_W+jetsUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lDown = fileIn.Get("shape5_W+jetsDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);

NominalHist.SetTitle("")

c1.cd()
pad1.cd()
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

legend.Draw()


pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2, 0.5, 1.5)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4, 0.5, 1.5)
for i in range(1,5):
  hratio1.SetBinError(i,0.001)
  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("AlphaFitSyst_Par1_Wjets.png")
c1.SaveAs("AlphaFitSyst_Par1_Wjets.pdf")
c1.SetLogy(1)
c1.SaveAs("AlphaFitSyst_Par1_Wjets_log.png")
c1.SaveAs("AlphaFitSyst_Par1_Wjets_log.pdf")
c1.SetLogy(0)
c1.Clear()
legend.Clear()

####################################################################
c1, pad1, pad2 = Utils.createCanvasPads();

#legend = r.TLegend(0.1,0.0,0.45,0.5);
legend = r.TLegend(0.6,0.5,0.95,0.9);

lUp1 = fileIn.Get("shape2_W+jetsUp");
lUp1.SetMarkerColor(2);
lUp1.SetMarkerSize(1);
lUp1.SetLineColor(2);
lUp1.SetFillStyle(0);
lDown1 = fileIn.Get("shape2_W+jetsDown");
lDown1.SetMarkerColor(4);
lDown1.SetMarkerSize(1);
lDown1.SetLineColor(4);
lDown1.SetFillStyle(0);

lUp0 = fileIn.Get("shape_W+jetsUp");
lUp0.SetMarkerColor(6);
lUp0.SetMarkerSize(1);
lUp0.SetLineColor(6);
lUp0.SetFillStyle(0);
lDown0 = fileIn.Get("shape_W+jetsDown");
lDown0.SetMarkerColor(8);
lDown0.SetMarkerSize(1);
lDown0.SetLineColor(8);
lDown0.SetFillStyle(0);

lUp = fileIn.Get("shape3_W+jetsUp");
lUp.SetMarkerColor(9);
lUp.SetMarkerSize(1);
lUp.SetLineColor(9);
lUp.SetFillStyle(0);
lDown = fileIn.Get("shape3_W+jetsDown");
lDown.SetMarkerColor(49);
lDown.SetMarkerSize(1);
lDown.SetLineColor(49);
lDown.SetFillStyle(0);

NominalHist.SetTitle("")

c1.cd()
pad1.cd()
NominalHist.Draw();
lUp1.Draw("same hist")
lDown1.Draw("same hist")
lUp0.Draw("same hist")
lDown0.Draw("same hist")
#lUp.Draw("same hist")
lDown.Draw("same hist")

legend.AddEntry(NominalHist,"nominial","lpe");
legend.AddEntry(lUp1,"+1#sigma shift (Par-1)","lp");
legend.AddEntry(lDown1,"-1#sigma shift (Par-1)","lp");
legend.AddEntry(lUp0,"+1#sigma shift (Par-0)","lp");
legend.AddEntry(lDown0,"-1#sigma shift (Par-0)","lp");
#legend.AddEntry(lUp,"+1#sigma shift (Alternate)","lp");
#legend.AddEntry(lDown,"-1#sigma shift (Alternate)","lp");
legend.AddEntry(lDown,"Exp(ax) fit","lp");
legend.Draw()


pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 9)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 49)
hratio3 = Utils.createRatio(NominalHist,  lDown0, "M_{WW} (GeV)", 8)
hratio4 = Utils.createRatio(NominalHist,  lUp0, "M_{WW} (GeV)", 6)
hratio5 = Utils.createRatio(NominalHist,  lUp1, "M_{WW} (GeV)", 2)
hratio6 = Utils.createRatio(NominalHist,  lDown1, "M_{WW} (GeV)", 4)
for i in range(1,5):
  hratio1.SetBinError(i,0.001)
  hratio2.SetBinError(i,0.001)
  hratio3.SetBinError(i,0.001)
  hratio4.SetBinError(i,0.001)
  hratio5.SetBinError(i,0.001)
  hratio6.SetBinError(i,0.001)

#hratio1.Draw("")
hratio2.Draw("")
hratio3.Draw("same ")
hratio4.Draw("same ")
hratio5.Draw("same ")
hratio6.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("WjetFit_all.png")
c1.SaveAs("WjetFit_all.pdf")
c1.Clear()


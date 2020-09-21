import ROOT as r
import sys
import Utils

#fileIn = r.TFile("ch1_splitted_TF1_WV.root","read");
fileIn = r.TFile("ch1_splitted_TF1_ZV.root","read");
#fileIn = r.TFile("Check_bkg_bug.root","read");
#fileIn = r.TFile("Check_bkg_bug_2500.root","read");
#r.gStyle.SetErrorX(0.0001);

c1 = r.TCanvas("c1");
c1, pad1, pad2 = Utils.createCanvasPads();
legend = r.TLegend(0.5,0.7,0.9,0.9);

NominalHist =fileIn.Get("diboson");
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.SetStats(0)
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("diboson_diboson_QCDScaleUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lUp.SetStats(0);
lDown = fileIn.Get("diboson_diboson_QCDScaleDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);
lDown.SetStats(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
#NominalHist.SetMaximum(20.0);
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

#for i in range(1,5):
#  hratio1.SetBinError(i,0.001)
#  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("diboson_QCDScale.png")
c1.SaveAs("diboson_QCDScale.pdf")
c1.Clear()

c1, pad1, pad2 = Utils.createCanvasPads();

NominalHist =fileIn.Get("VVjjQCD");
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.SetStats(0)
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("VVjjQCD_pdf_ggbarUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lUp.SetStats(0);
lDown = fileIn.Get("VVjjQCD_pdf_ggbarDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);
lDown.SetStats(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
#NominalHist.SetMaximum(40.0);
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

#legend.AddEntry(NominalHist,"nominial","lpe");
#legend.AddEntry(lUp,"+1#sigma shift","lp");
#legend.AddEntry(lDown,"-1#sigma shift","lp");
legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

#for i in range(1,5):
#  hratio1.SetBinError(i,0.001)
#  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("VVjjQCD_PDFScale.png")
c1.SaveAs("VVjjQCD_PDFScale.pdf")
c1.Clear()

c1, pad1, pad2 = Utils.createCanvasPads();

NominalHist =fileIn.Get("VVjjQCD");
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.SetStats(0)
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("VVjjQCD_VVjjQCD_QCDScaleUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lUp.SetStats(0);
lDown = fileIn.Get("VVjjQCD_VVjjQCD_QCDScaleDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);
lDown.SetStats(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
#NominalHist.SetMaximum(40.0);
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

#legend.AddEntry(NominalHist,"nominial","lpe");
#legend.AddEntry(lUp,"+1#sigma shift","lp");
#legend.AddEntry(lDown,"-1#sigma shift","lp");
legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2, 0.3, 2.0)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4, 0.3, 2.0)

#for i in range(1,5):
#  hratio1.SetBinError(i,0.001)
#  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("VVjjQCD_QCDScale.png")
c1.SaveAs("VVjjQCD_QCDScale.pdf")
c1.Clear()

c1, pad1, pad2 = Utils.createCanvasPads();

NominalHist =fileIn.Get("diboson");
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.SetStats(0)
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("diboson_pdf_qqbarUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lUp.SetStats(0);
lDown = fileIn.Get("diboson_pdf_qqbarDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);
lDown.SetStats(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
#NominalHist.SetMaximum(20.0);
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

#legend.AddEntry(NominalHist,"nominial","lpe");
#legend.AddEntry(lUp,"+1#sigma shift","lp");
#legend.AddEntry(lDown,"-1#sigma shift","lp");
legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

#for i in range(1,5):
#  hratio1.SetBinError(i,0.001)
#  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("diboson_PDFScale.png")
c1.SaveAs("diboson_PDFScale.pdf")
c1.Clear()


c1, pad1, pad2 = Utils.createCanvasPads();

NominalHist =fileIn.Get("CH_WZToLL_M800");
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.SetStats(0)
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("CH_WZToLL_M800_pdf_qqbarUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lUp.SetStats(0);
lDown = fileIn.Get("CH_WZToLL_M800_pdf_qqbarDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);
lDown.SetStats(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
#NominalHist.SetMaximum(20.0);
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

#legend.AddEntry(NominalHist,"nominial","lpe");
#legend.AddEntry(lUp,"+1#sigma shift","lp");
#legend.AddEntry(lDown,"-1#sigma shift","lp");
legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

#for i in range(1,5):
#  hratio1.SetBinError(i,0.001)
#  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("CH_WZToLL_M800_PDFScale.png")
c1.SaveAs("CH_WZToLL_M800_PDFScale.pdf")
c1.Clear()

c1, pad1, pad2 = Utils.createCanvasPads();

NominalHist =fileIn.Get("CH_WZToLL_M900");
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.SetStats(0)
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("CH_WZToLL_M900_pdf_qqbarUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lUp.SetStats(0);
lDown = fileIn.Get("CH_WZToLL_M900_pdf_qqbarDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);
lDown.SetStats(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
NominalHist.SetMaximum(30.0);
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

#legend.AddEntry(NominalHist,"nominial","lpe");
#legend.AddEntry(lUp,"+1#sigma shift","lp");
#legend.AddEntry(lDown,"-1#sigma shift","lp");
legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

#for i in range(1,5):
#  hratio1.SetBinError(i,0.001)
#  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("CH_WZToLL_M900_PDFScale.png")
c1.SaveAs("CH_WZToLL_M900_PDFScale.pdf")
c1.Clear()


c1, pad1, pad2 = Utils.createCanvasPads();

NominalHist =fileIn.Get("CH_WZToLL_M800");
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.SetStats(0)
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("CH_WZToLL_M800_Higgs_QCDScaleUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lUp.SetStats(0);
lDown = fileIn.Get("CH_WZToLL_M800_Higgs_QCDScaleDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);
lDown.SetStats(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
#NominalHist.SetMaximum(20.0);
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

#legend.AddEntry(NominalHist,"nominial","lpe");
#legend.AddEntry(lUp,"+1#sigma shift","lp");
#legend.AddEntry(lDown,"-1#sigma shift","lp");
legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

#for i in range(1,5):
#  hratio1.SetBinError(i,0.001)
#  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("CH_WZToLL_M800_QCDScale.png")
c1.SaveAs("CH_WZToLL_M800_QCDScale.pdf")
c1.Clear()

c1, pad1, pad2 = Utils.createCanvasPads();

NominalHist =fileIn.Get("CH_WZToLL_M900");
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.SetStats(0)
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("CH_WZToLL_M900_Higgs_QCDScaleUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lUp.SetStats(0);
lDown = fileIn.Get("CH_WZToLL_M900_Higgs_QCDScaleDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);
lDown.SetStats(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
#NominalHist.SetMaximum(140.0);
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

#legend.AddEntry(NominalHist,"nominial","lpe");
#legend.AddEntry(lUp,"+1#sigma shift","lp");
#legend.AddEntry(lDown,"-1#sigma shift","lp");
legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

#for i in range(1,5):
#  hratio1.SetBinError(i,0.001)
#  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("CH_WZToLL_M900_QCDScale.png")
c1.SaveAs("CH_WZToLL_M900_QCDScale.pdf")
c1.Clear()


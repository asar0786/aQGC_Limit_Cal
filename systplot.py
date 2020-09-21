import ROOT as r
import sys
import Utils

fileIn = r.TFile("ch1_splitted_TF1_WV.root","read");

r.gStyle.SetOptStat(0)

#c1 = r.TCanvas("c1");
c1, pad1, pad2 = Utils.createCanvasPads(1,0);
legend = r.TLegend(0.5,0.7,0.9,0.9);

l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);

NominalHist =fileIn.Get(sys.argv[1]);
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetTitle("");
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get(sys.argv[1]+"_"+sys.argv[1]+"_bbb_"+sys.argv[1]+"_bin_4Up");
#lUp.SetMarkerColor(2);
#lUp.SetLineColor(2);
lDown = fileIn.Get(sys.argv[1]+"_"+sys.argv[1]+"_bbb_"+sys.argv[1]+"_bin_4Down");
#lDown.SetMarkerColor(3);
#lDown.SetLineColor(4);
#
#c1.cd()
#pad1.cd()
#NominalHist.Draw("P");
#lUp.Draw("same hist")
#lDown.Draw("same hist")
#
#legend.AddEntry(NominalHist,"nominial","lpe");
#legend.AddEntry(lUp,"+1#sigma shift","lp");
#legend.AddEntry(lDown,"-1#sigma shift","lp");
#legend.Draw()
#
#pad2.cd()
#hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
#hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)
#
#hratio1.Draw()
#hratio2.Draw("same")
#l.Draw();
#
#c1.SaveAs(sys.argv[1]+"_bbb_"+sys.argv[1]+"_bin_4Up.png")
#c1.SaveAs(sys.argv[1]+"_bbb_"+sys.argv[1]+"_bin_4Up.pdf")
#c1.SetLogy(1)
##c1.SaveAs(sys.argv[1]+"_bbb_"+sys.argv[1]+"_bin_4Up_log.png")
##c1.SaveAs(sys.argv[1]+"_bbb_"+sys.argv[1]+"_bin_4Up_log.pdf")
#c1.SetLogy(0)
c1.Clear()
#c1, pad1, pad2 = Utils.createCanvasPads(1,0);
#
#NominalHist =fileIn.Get(sys.argv[1]);
#NominalHist.SetMarkerColor(1);
#NominalHist.SetLineColor(1);
#NominalHist.SetTitle("");
#NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")
#
#lUp = fileIn.Get(sys.argv[1]+"_CMS_btagLFUp");
#lUp.SetMarkerColor(2);
#lUp.SetLineColor(2);
#lDown = fileIn.Get(sys.argv[1]+"_CMS_btagLFDown");
#lDown.SetMarkerColor(3);
#lDown.SetLineColor(4);
#
#c1.cd()
#pad1.cd()
#NominalHist.Draw();
#lUp.Draw("same hist")
#lDown.Draw("same hist")
#
#legend.Draw()
#
#pad2.cd()
#hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
#hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)
#
#hratio1.Draw()
#hratio2.Draw("same")
#l.Draw();
#
#c1.SaveAs(sys.argv[1]+"_CMS_btagLFUp.png")
#c1.SaveAs(sys.argv[1]+"_CMS_btagLFUp.pdf")
#c1.SetLogy(1)
##c1.SaveAs(sys.argv[1]+"_CMS_btagLFUp_log.png")
##c1.SaveAs(sys.argv[1]+"_CMS_btagLFUp_log.pdf")
#c1.SetLogy(0)
#c1.Clear()
#c1, pad1, pad2 = Utils.createCanvasPads(1,0);
#
#NominalHist =fileIn.Get(sys.argv[1]);
#NominalHist.SetMarkerColor(1);
#NominalHist.SetLineColor(1);
#NominalHist.SetTitle("");
#NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")
#
#lUp = fileIn.Get(sys.argv[1]+"_CMS_btagHFUp");
#lUp.SetMarkerColor(2);
#lUp.SetLineColor(2);
#lDown = fileIn.Get(sys.argv[1]+"_CMS_btagHFDown");
#lDown.SetMarkerColor(3);
#lDown.SetLineColor(4);
#
#c1.cd()
#pad1.cd()
#NominalHist.Draw();
#lUp.Draw("same hist")
#lDown.Draw("same hist")
#
#legend.Draw()
#
#pad2.cd()
#hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
#hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)
#
#hratio1.Draw()
#hratio2.Draw("same")
#l.Draw();
#
#c1.SaveAs(sys.argv[1]+"_CMS_btagHFUp.png")
#c1.SaveAs(sys.argv[1]+"_CMS_btagHFUp.pdf")
#c1.SetLogy(1)
##c1.SaveAs(sys.argv[1]+"_CMS_btagHFUp_log.png")
##c1.SaveAs(sys.argv[1]+"_CMS_btagHFUp_log.pdf")
#c1.SetLogy(0)
#c1.Clear()
#c1, pad1, pad2 = Utils.createCanvasPads(1,0);
#
#NominalHist =fileIn.Get(sys.argv[1]);
#NominalHist.SetMarkerColor(1);
#NominalHist.SetLineColor(1);
#NominalHist.SetTitle("");
#NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")
#
#lUp = fileIn.Get(sys.argv[1]+"_CMS_puUp");
#lUp.SetMarkerColor(2);
#lUp.SetLineColor(2);
#lDown = fileIn.Get(sys.argv[1]+"_CMS_puDown");
#lDown.SetMarkerColor(3);
#lDown.SetLineColor(4);
#
#c1.cd()
#pad1.cd()
#NominalHist.Draw();
#lUp.Draw("same hist")
#lDown.Draw("same hist")
#
#legend.Draw()
#
#pad2.cd()
#hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
#hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)
#
#hratio1.Draw("P")
#hratio2.Draw("same P")
#l.Draw();
#
#c1.SaveAs(sys.argv[1]+"_CMS_puUp.png")
#c1.SaveAs(sys.argv[1]+"_CMS_puUp.pdf")
#c1.SetLogy(1)
##c1.SaveAs(sys.argv[1]+"_CMS_puUp_log.png")
##c1.SaveAs(sys.argv[1]+"_CMS_puUp_log.pdf")
#c1.SetLogy(0)
#c1.Clear()
#c1, pad1, pad2 = Utils.createCanvasPads(1,0);
#
#NominalHist =fileIn.Get(sys.argv[1]);
#NominalHist.SetMarkerColor(1);
#NominalHist.SetLineColor(1);
#NominalHist.SetTitle("");
#NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")
#
#lUp = fileIn.Get(sys.argv[1]+"_CMS_res_metUp");
#lUp.SetMarkerColor(2);
#lUp.SetLineColor(2);
#lDown = fileIn.Get(sys.argv[1]+"_CMS_res_metDown");
#lDown.SetMarkerColor(3);
#lDown.SetLineColor(4);
#
#c1.cd()
#pad1.cd()
#NominalHist.Draw();
#lUp.Draw("same hist")
#lDown.Draw("same hist")
#
#legend.Draw()
#
#pad2.cd()
#hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
#hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)
#
#hratio1.Draw("P")
#hratio2.Draw("same P")
#l.Draw();
#
#c1.SaveAs(sys.argv[1]+"_CMS_res_metUp.png")
#c1.SaveAs(sys.argv[1]+"_CMS_res_metUp.pdf")
#c1.SetLogy(1)
##c1.SaveAs(sys.argv[1]+"_CMS_res_metUp_log.png")
##c1.SaveAs(sys.argv[1]+"_CMS_res_metUp_log.pdf")
#c1.SetLogy(0)
#c1.Clear()
c1, pad1, pad2 = Utils.createCanvasPads(1,0);

NominalHist =fileIn.Get(sys.argv[1]);
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetTitle("");
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get(sys.argv[1]+"_CMS_scale_jUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lDown = fileIn.Get(sys.argv[1]+"_CMS_scale_jDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);


legend.AddEntry(NominalHist,"nominial","lpe");
legend.AddEntry(lUp,"+1#sigma shift","lp");
legend.AddEntry(lDown,"-1#sigma shift","lp");

c1.cd()
pad1.cd()
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

hratio1.Draw("P")
hratio2.Draw("same P")
l.Draw("same");


c1.SaveAs(sys.argv[1]+"_CMS_scale_jUp.png")
c1.SaveAs(sys.argv[1]+"_CMS_scale_jUp.pdf")
c1.SetLogy(1)
#c1.SaveAs(sys.argv[1]+"_CMS_scale_jUp_log.png")
#c1.SaveAs(sys.argv[1]+"_CMS_scale_jUp_log.pdf")
c1.SetLogy(0)
c1.Clear()
c1, pad1, pad2 = Utils.createCanvasPads(1,0);

NominalHist =fileIn.Get(sys.argv[1]);
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetTitle("");
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get(sys.argv[1]+"_CMS_scale_lUp");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lDown = fileIn.Get(sys.argv[1]+"_CMS_scale_lDown");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);

c1.cd()
pad1.cd()
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

hratio1.Draw("P")
hratio2.Draw("same P")
l.Draw();

c1.SaveAs(sys.argv[1]+"_CMS_scale_lUp.png")
c1.SaveAs(sys.argv[1]+"_CMS_scale_lUp.pdf")
c1.SetLogy(1)
#c1.SaveAs(sys.argv[1]+"_CMS_scale_lUp_log.png")
#c1.SaveAs(sys.argv[1]+"_CMS_scale_lUp_log.pdf")
c1.SetLogy(0)
c1.Clear()

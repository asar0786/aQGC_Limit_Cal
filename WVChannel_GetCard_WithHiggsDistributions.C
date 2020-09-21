#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <assert.h>
#include "TROOT.h"
#include "TLatex.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TSystem.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include <ctime>

#include "tdrstyle.C"
#include "utils.C" // Tokenize
#include "CMS_lumi.C"

#include <Python.h>



typedef struct SampleInfo_t {
  int     index;
  TString samplename;
  TString treefilename;
  double xsecpblumi;
  double otherscale;
  int    nMCevents;
  int	 MCnegEvent;
  int    colorcode;
  int    stackit;
}
SampleInfo_t;
using namespace std;

double intLUMIinvpb;
double fs0[43] = {-50.0, -45.0, -40.0, -35.0, -30.0, -20.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.0, -2.5, -2.0, -1.5, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 20.0, 30.0, 35.0, 40.0, 45.0, 50.0};
double fs1[33] = {-35, -33, -30, -25, -20, -15, -10, -7.5, -5.0, -4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10, 15, 20, 25, 30, 33, 35};
double fm0[41] = {-10, -9, -8, -7, -6, -5, -4, -3, -2.0, -1.5, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3, 4, 5, 6, 7, 8, 9, 10};
double fm1[37] = {-30, -28, -23, -21, -18, -15, -13, -10, -5.0, -3.0, -2.5, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.5, 3.0, 5.0, 10, 13, 15, 18, 21, 23.0, 28, 30};
double fm2[31] = {-60.0, -55.0, -50.0, -45.0, -40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0, -6.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 6.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0};
double fm3[31] = {-105.0, -95.0, - 85.0, -75.0, -65.0, -55.0, -44.0, -31.0, -21.0, -13.0, -8.0, -5.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 31.0, 44.0, 55.0, 65.0, 75.0, 85.0, 95.0, 105.0};
double fm4[37] = {-130.0,-121.0,-115.0,105.0,-95.0, -85.0, -75.0, -65.0, -55.0, -44.0, -31.0, -21.0, -13.0, -8.0, -5.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 31.0, 44.0, 55.0, 65.0, 75.0, 85.0, 95.0, 105.0, 115.0, 121.0, 130.0};
double fm5[45] = {-200.0, -190.0, -170.0, -150.0, -130.0,-121.0,-115.0,105.0,-95.0, -85.0, -75.0, -65.0, -55.0, -44.0, -31.0, -21.0, -13.0, -8.0, -5.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 31.0, 44.0, 55.0, 65.0, 75.0, 85.0, 95.0, 105.0, 115.0, 121.0, 130.0, 150.0, 170.0, 190.0, 200.0};
double fm6[35] = {-20.0, -18.0, -15.0,- -12.0, -10.0, -7.0, -5.0, -3.0, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.7, -0.5, -0.2, 0.0, 0.2, 0.5, 0.7, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 5.0, 7.0, 10.0, 12.0, 15.0, 18.0, 20.0};
double fm7[33] = {-40, -35, -30, -25, -20, -15, -10, -5.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0};
double ft0[35] = {-2.0, -1.8, -1.4, -1.2, -1.0, -0.7, -0.5, -0.3, -0.2, -0.18, -0.14, -0.12, -0.10, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.18, 0.20, 0.30, 0.50, 0.7, 1.0, 1.2, 1.4, 1.8, 2.0};
double ft1[35] = {-2.0, -1.8, -1.4, -1.2, -1.0, -0.7, -0.5, -0.3, -0.2, -0.18, -0.14, -0.12, -0.10, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.18, 0.20, 0.30, 0.50, 0.7, 1.0, 1.2, 1.4, 1.8, 2.0};
double ft2[35] = {-4.5, -3.9, -3.4, -2.9, -2.5, -1.7, -1.2, -0.9, -0.7, -0.5, -0.32, -0.26, -0.20, -0.14, -0.08, -0.04, -0.02, 0, 0.02, 0.04, 0.08, 0.14, 0.20, 0.26, 0.32, 0.5, 0.7, 0.9, 1.2, 1.7, 2.5, 2.9, 3.4, 3.9, 4.5};
double ft5[39] = {-25.0, -22.0, -20.0, -18.0, -15.0, -12.0, -10.0, -7.0, -5.0, -3.0, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.7, -0.5, -0.2, 0.0, 0.2, 0.5, 0.7, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 5.0, 7.0, 10.0, 12.0, 15.0, 18.0, 20.0, 22.0, 25.0};
double ft6[43] = {-29.0, -27.0, -25.0, -22.0, -20.0, -18.0, -15.0, -12.0, -10.0, -7.0, -5.0, -3.0, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.7, -0.5, -0.2, 0.0, 0.2, 0.5, 0.7, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 5.0, 7.0, 10.0, 12.0, 15.0, 18.0, 20.0, 22.0, 25.0, 27.0, 29.0};
double ft7[51] = {-70.0, -65.0, -60.0, -55.0, -50.0, -45.0, -40.0, -35.0, -30.0, -20.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.0, -2.5, -2.0, -1.5, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 20.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0};
double ft8[35] = {-2.0, -1.8, -1.4, -1.2, -1.0, -0.7, -0.5, -0.3, -0.2, -0.18, -0.14, -0.12, -0.10, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.18, 0.20, 0.30, 0.50, 0.7, 1.0, 1.2, 1.4, 1.8, 2.0};
double ft9[41] = {-10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

//======================================================================
class Sample {
public:
  Sample(const SampleInfo_t& sinfo) {
    info_ = sinfo;
    tree_ = 0;
    th1f_ = 0;
    //cout << "sample = " << name_ << endl;
    TFile *f = new TFile (sinfo.treefilename, "READ"); if (!f) { cerr << "Couldn't find file " << sinfo.treefilename << endl; return; }
    //TFile *f =  TFile::Open("root://cmsxrootd.fnal.gov/"+sinfo.treefilename, "READ"); if (!f) { cerr << "Couldn't find file " << sinfo.treefilename << endl; return; }
    tree_ =  (TTree *)f->Get("Events"); if (!tree_) { cerr << "Couldn't find tree Events in file " << sinfo.treefilename << endl; return; }
    th1f_ = (TH1F *)f->Get("TotalEvents"); if (!th1f_) {cerr << "Couldn't find TotalEvents in the file " << sinfo.treefilename << endl; return;}
  }
  ~Sample() { if (tree_) delete tree_; }
  TTree *Tree() const { return tree_; }
//  TH1F *TInputHist() {return th1f_; }
  TH1F *TInputHist() const {return th1f_; }
  TString name() const { return info_.samplename; }
  TString filename() const { return info_.treefilename; }
  bool stackit() const { return (info_.stackit != 0); }
  int colorcode() const { return info_.colorcode; }
  double otherscale() const { return info_.otherscale; }
  double cross() const {return info_.xsecpblumi; }
  int mcevent() const {return info_.nMCevents; }
  int mcevent_neg() const {return info_.MCnegEvent; }
  private:
    SampleInfo_t info_;
    TTree *tree_;
    TH1F *th1f_;
};

//======================================================================
//
void loadSamples(const char *filename,vector<Sample *>& samples)
{
  FILE *fp = fopen(filename,"r");
  if (!fp) {
    cout << "Error, file " << TString(filename) << " not found." << endl;
    exit(-1);
  }

  char line[512];

  intLUMIinvpb=-1; // obvious error condition

  for (int i=0; !feof(fp) && fgets(line,512,fp); i++) {
    if (!strlen(line) || line[0]=='#') continue; // comments are welcome

    string strline(line);
    strline.pop_back();     // shed the \n
    vector<string> fields;

    // expect columns with fields cutname, cutvalue, possible embedded spaces both
    // within and between, so " " or "\t" cannot be used as delimiters. Require quotes
    // instead.
    //
    Tokenize(strline,fields, " \t");

    //for (size_t j=0; j<fields.size(); j++)
    //cout << j << ": \"" << fields[j] << "\"" << endl;

    assert (fields.size()==8);

    SampleInfo_t s;
    s.index        = i;
    s.samplename   = fields[0];
    s.treefilename = fields[1];
    s.xsecpblumi   = str2dbl(fields[2]);
    s.otherscale   = str2dbl(fields[3]);
    s.nMCevents    = str2int(fields[4]);
    s.MCnegEvent   = str2int(fields[5]);
    s.colorcode    = str2int(fields[6]);
    s.stackit      = str2int(fields[7]);
    
    //if (!s.samplename.EqualTo("aQGC")) continue;
    
    cout << "Loading sample " << s.samplename << " -> " << s.treefilename << endl;
    
    
    if (!samples.size()) {
      if (s.samplename.EqualTo("data")) {
	intLUMIinvpb = s.xsecpblumi;
	s.xsecpblumi = 1;
	cout << "intLUMI = " << intLUMIinvpb << " pb^-1" << endl;
      } else {
	cerr << "First sample in the table must be 'data'" << endl;
	//exit(-1);
      }
    } else {
      s.otherscale *= intLUMIinvpb;
    }
    
    samples.push_back(new Sample(s) );
  }
}                                                         // loadSamples
//======================================================================

void loadCutString(const char *filename, TString& cutstring)
{
  FILE *fp = fopen(filename,"r");
  if (!fp) {
    cout << "Error, file " << TString(filename) << " not found." << endl;
    exit(-1);
  }

  char line[512];

  for (int i=0; !feof(fp) && fgets(line,512,fp); i++) {
    if (!strlen(line) || line[0]=='#') continue; // comments are welcome

    if (cutstring.Length()) cutstring += " && ";

    string strline(line);
    strline.pop_back();     // shed the \n
    vector<string> fields;

    // expect columns with fields cutname, cutvalue, possible embedded spaces both
    // within and between, so " " or "\t" cannot be used as delimiters. Require quotes
    // instead.
    //
    Tokenize(strline,fields, "\"");

    //for (size_t j=0; j<fields.size(); j++)
    //cout << j << ": \"" << fields[j] << "\"" << endl;

    assert (fields.size()==3);
    cutstring += TString(fields.at(2));
  }
}                                                       // loadCutString

//======================================================================

void model(const char *samplefilename,
	   const TString OutPutRootFileName = "ch1_splitted_TF1_WV")
{
  cout<< "done..." << endl;

  vector<Sample *> samples;
  
  loadSamples(samplefilename,samples);
  
  // Data
  
  Sample *sdata = samples[0];
  
  if (sdata->Tree())
    cout << "ndata =" << sdata->Tree()->GetEntries() <<endl;

  TFile* wjetBkgSystFile = new TFile("WV_bkg_estimation_4Bins_50GeVLepCut.root","READ");
  
  TH1F* wjet = (TH1F*)wjetBkgSystFile->Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_Nominal");
  TH1F* wjetup = (TH1F*)wjetBkgSystFile->Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_Par0Up");
  TH1F* wjetdown = (TH1F*)wjetBkgSystFile->Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_Par0Down");
  TH1F* wjetup1 = (TH1F*)wjetBkgSystFile->Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_Par1Up");
  TH1F* wjetdown1 = (TH1F*)wjetBkgSystFile->Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_Par1Down");
  TH1F* wjetup2 = (TH1F*)wjetBkgSystFile->Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_AlternateShape_Up");
  TH1F* wjetdown2 = (TH1F*)wjetBkgSystFile->Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_AlternateShape_Down");
  TH1F* wjetup3 = (TH1F*)wjetBkgSystFile->Get("AlphaSyst_Vjet_SR_4bins_Par0Up");
  TH1F* wjetdown3 = (TH1F*)wjetBkgSystFile->Get("AlphaSyst_Vjet_SR_4bins_Par0Down");
  TH1F* wjetup4 = (TH1F*)wjetBkgSystFile->Get("AlphaSyst_Vjet_SR_4bins_Par1Up");
  TH1F* wjetdown4 = (TH1F*)wjetBkgSystFile->Get("AlphaSyst_Vjet_SR_4bins_Par1Down");
  
  //============================================================
  //  VARIABLE LOOP
  //============================================================
  int NBINS = 8;
  double MINRange = 600;
  double MAXRange = 2500;
  //double massLEdges[5] = {150, 300, 500, 1000, 1500};
 // double massLEdges[5] = {600, 1075, 1550, 2025, 2500};
  double massLEdges[9]={600, 700, 800, 900, 1000, 1200, 1500, 2000, 2500};


  TH1 *h1 = new TH1D("EWK", "EWK", 30, 0, 3000);
  TH1 *hists[79];
  TH1 *ChargedHist[561];
  TH1 *ChargedHistQCD[198];
  TH1 *ChargedHistPDF[3300];

  const char* HistName[79] = {	"data_obs",						// 0
  			"diboson", "diboson_CMS_scale_lUp", "diboson_CMS_scale_lDown", "diboson_CMS_scale_jUp", "diboson_CMS_scale_jDown", "diboson_CMS_res_metUp", "diboson_CMS_res_metDown", 
			"diboson_CMS_puUp", "diboson_CMS_puDown",	"diboson_CMS_btagHFUp", "diboson_CMS_btagHFDown", "diboson_CMS_btagLFUp", "diboson_CMS_btagLFDown", // 1	
			"VVjjQCD", "VVjjQCD_CMS_scale_lUp", "VVjjQCD_CMS_scale_lDown", "VVjjQCD_CMS_scale_jUp", "VVjjQCD_CMS_scale_jDown", "VVjjQCD_CMS_res_metUp", "VVjjQCD_CMS_res_metDown", 
			"VVjjQCD_CMS_puUp", "VVjjQCD_CMS_puDown", "VVjjQCD_CMS_btagHFUp", "VVjjQCD_CMS_btagHFDown", "VVjjQCD_CMS_btagLFUp", "VVjjQCD_CMS_btagLFDown", // 14
			"top", "top_CMS_scale_lUp", "top_CMS_scale_lDown", "top_CMS_scale_jUp", "top_CMS_scale_jDown",	"top_CMS_res_metUp", "top_CMS_res_metDown", 
			"top_CMS_puUp", "top_CMS_puDown", "top_CMS_btagHFUp", "top_CMS_btagHFDown", "top_CMS_btagLFUp", "top_CMS_btagLFDown",		// 27
			"Vjets", "Vjets_CMS_scale_lUp", "Vjets_CMS_scale_lDown", "Vjets_CMS_scale_jUp", "Vjets_CMS_scale_jDown", "Vjets_CMS_res_metUp", "Vjets_CMS_res_metDown", 
			"Vjets_CMS_puUp", "Vjets_CMS_puDown", "Vjets_CMS_btagHFUp", "Vjets_CMS_btagHFDown", "Vjets_CMS_btagLFUp", "Vjets_CMS_btagLFDown",	// 40
			"CH_WZ", "CH_WZ_CMS_scale_lUp", "CH_WZ_CMS_scale_lDown", "CH_WZ_CMS_scale_jUp", "CH_WZ_CMS_scale_jDown", "CH_WZ_CMS_res_metUp", "CH_WZ_CMS_res_metDown", 
			"CH_WZ_CMS_puUp", "CH_WZ_CMS_puDown", "CH_WZ_CMS_btagHFUp", "CH_WZ_CMS_btagHFDown", "CH_WZ_CMS_btagLFUp", "CH_WZ_CMS_btagLFDown",	// 53
			"DCH_WW", "DCH_WW_CMS_scale_lUp", "DCH_WW_CMS_scale_lDown", "DCH_WW_CMS_scale_jUp", "DCH_WW_CMS_scale_jDown", "DCH_WW_CMS_res_metUp", "DCH_WW_CMS_res_metDown", 
			"DCH_WW_CMS_puUp", "DCH_WW_CMS_puDown", "DCH_WW_CMS_btagHFUp", "DCH_WW_CMS_btagHFDown", "DCH_WW_CMS_btagLFUp", "DCH_WW_CMS_btagLFDown"	// 66
			};
  

  TString HiggsSampleName[3] = { "CH_WZToLL", "CH_WZToLNu", "DCH_WW"};
  TString MassPoint[11] = { "_M200", "_M300", "_M400", "_M500", "_M600", "_M700", "_M800", "_M900", "_M1000", "_M1500", "_M2000"};
  TString Syst[17] = {"", "_CMS_scale_lUp", "_CMS_scale_lDown", "_CMS_scale_jUp", "_CMS_scale_jDown", "_CMS_res_metUp", "_CMS_res_metDown", "_CMS_puUp", "_CMS_puDown",       "_CMS_btagHFUp", "_CMS_btagHFDown", "_CMS_btagLFUp", "_CMS_btagLFDown", "_Higgs_QCDScaleUp", "_Higgs_QCDScaleDown", "_pdf_qqbarUp", "_pdf_qqbarDown" };

  
  TH1 *histo_aqgc[678];
  for(int j=0;j<678;j++)
    {
      stringstream ss;
      ss << j;
      string temp = ss.str();
      const char* name = temp.c_str();
      histo_aqgc[j] = new TH1D(name, name, NBINS, massLEdges);
      histo_aqgc[j]->Sumw2();
    }
/*
  //let's define few histograms for the uncertainties
  TH1D* histo_diboson_EWK_CMS_QCDScaleBounding[6];
  TH1D* histo_VVjjQCD_EWK_CMS_QCDScaleBounding[6];
  for(int i = 0; i<6; i++)
    {
      stringstream ss;
      ss << i;
      string temp = ss.str()+"_QCD";
      string temp1 = ss.str()+"_diboson";
      const char* name = temp.c_str();
      const char* name1 = temp1.c_str();
      histo_diboson_EWK_CMS_QCDScaleBounding[i] = new TH1D(name1, name1, NBINS,massLEdges);
      histo_diboson_EWK_CMS_QCDScaleBounding[i]->Sumw2();
      histo_VVjjQCD_EWK_CMS_QCDScaleBounding[i] = new TH1D(name, name, NBINS,massLEdges);
      histo_VVjjQCD_EWK_CMS_QCDScaleBounding[i]->Sumw2();
    }
 TH1D* histo_diboson_EWK_CMS_PDFScaleBounding[100];
 TH1D* histo_VVjjQCD_EWK_CMS_PDFScaleBounding[100];
  for(int i = 0; i<100; i++)
    {
      stringstream ss;
      ss << i;
      string temp = ss.str()+"_diboson";
      string temp1 = ss.str()+"_VVJJQCD";
      const char* name = temp.c_str();
      const char* name1 = temp1.c_str();
      histo_diboson_EWK_CMS_PDFScaleBounding[i] = new TH1D(name, name, NBINS,massLEdges);
      histo_diboson_EWK_CMS_PDFScaleBounding[i]->Sumw2();
      histo_VVjjQCD_EWK_CMS_PDFScaleBounding[i] = new TH1D(name1, name1, NBINS,massLEdges);
      histo_VVjjQCD_EWK_CMS_PDFScaleBounding[i]->Sumw2();
    }

 TH1D* histo_VVjjQCD_EWK_CMS_QCDScaleBounding_Up = new TH1D("VVjjQCD_VVjjQCD_QCDScaleUp","", NBINS,massLEdges);
 histo_VVjjQCD_EWK_CMS_QCDScaleBounding_Up->Sumw2();

 TH1D* histo_VVjjQCD_EWK_CMS_QCDScaleBounding_Down = new TH1D("VVjjQCD_VVjjQCD_QCDScaleDown","", NBINS,massLEdges);
 histo_VVjjQCD_EWK_CMS_QCDScaleBounding_Down->Sumw2();

  
 TH1D* histo_VVjjQCD_EWK_CMS_PDFScaleBounding_Up = new TH1D("VVjjQCD_pdf_ggbarUp","", NBINS,massLEdges);
 histo_VVjjQCD_EWK_CMS_PDFScaleBounding_Up->Sumw2();

 TH1D* histo_VVjjQCD_EWK_CMS_PDFScaleBounding_Down = new TH1D("VVjjQCD_pdf_ggbarDown","", NBINS,massLEdges);
 histo_VVjjQCD_EWK_CMS_PDFScaleBounding_Down->Sumw2();
  

 TH1D* histo_diboson_EWK_CMS_QCDScaleBounding_Up = new TH1D("diboson_diboson_QCDScaleUp","", NBINS,massLEdges);
 histo_diboson_EWK_CMS_QCDScaleBounding_Up->Sumw2();

 TH1D* histo_diboson_EWK_CMS_QCDScaleBounding_Down = new TH1D("diboson_diboson_QCDScaleDown","", NBINS,massLEdges);
 histo_diboson_EWK_CMS_QCDScaleBounding_Down->Sumw2();
  
 TH1D* histo_diboson_EWK_CMS_PDFScaleBounding_Up = new TH1D("diboson_pdf_qqbarUp","", NBINS,massLEdges);
 histo_diboson_EWK_CMS_PDFScaleBounding_Up->Sumw2();

 TH1D* histo_diboson_EWK_CMS_PDFScaleBounding_Down = new TH1D("diboson_pdf_qqbarDown","", NBINS,massLEdges);
 histo_diboson_EWK_CMS_PDFScaleBounding_Down->Sumw2();
  
*/
  // histo for JES, JER, UP, Btag, LEP up/down uncertanities
  for (int i=0; i<79; i++)
  {
//    hists[i] = new TH1D(HistName[i],HistName[i], NBINS,MINRange,MAXRange);
    hists[i] = new TH1D(HistName[i],HistName[i], NBINS,massLEdges);
    hists[i]->Sumw2();
  }
/*
  int HistCount = 0;
  for (int i=0; i<3; i++)
     for (int j=0; j<11; j++)
        for (int k=0; k<17; k++)
	{
	   TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
	   cout<< "Hist name will be : " << HistCount << "\t" << name << endl;
	   ChargedHist[HistCount] = new TH1D(name , name, NBINS,massLEdges); 
	   ChargedHist[HistCount]->Sumw2();
	   HistCount++;
	}
  HistCount = 0;
  for (int i=0; i<3; i++)
    for (int j=0; j<11; j++)
      for (int k=0; k<6; k++)
	{
	  TString name = HiggsSampleName[i]+MassPoint[j]+"_CMS_QCDScale";
	  cout<< "Hist name will be : " << HistCount << "\t" << name << endl;
	  ChargedHistQCD[HistCount] = new TH1D(name , name, NBINS, massLEdges); 
	  ChargedHistQCD[HistCount]->Sumw2();
	  HistCount++;
	}
  HistCount = 0;
  for (int i=0; i<3; i++)
    for (int j=0; j<11; j++)
      for (int k=0; k<100; k++)
	{
	  TString name = HiggsSampleName[i]+MassPoint[j]+"_CMS_PDFScale";
	  cout<< "Hist name will be : " << HistCount << "\t" << name << endl;
	  ChargedHistPDF[HistCount] = new TH1D(name , name, NBINS,massLEdges); 
	  ChargedHistPDF[HistCount]->Sumw2();
	  HistCount++;
	}
  HistCount = 0;*/
  //============================================================
  // DRAW THE VARIABLE FOR ALL SAMPLES, CREATE HISTOS
  //============================================================    
  
  for (size_t isamp=0; isamp<samples.size(); isamp++) {
    Sample *s = samples[isamp];
    
    double xsec = s->cross();
    double otherscale = s->otherscale();
    int nmc =  s->mcevent();
    int nneg = s->mcevent_neg();

    cout<< "Running for  ==>  " << s->name() << "\t" << endl;
    if(s->name().EqualTo("data")) { 
    	//intLUMIinvpb = s->cross(); 
      cout<< "Lumi = " << intLUMIinvpb << endl;
    }

  const float MUON_MASS = 0.1056583745;
  const float ELE_MASS  = 0.000511;    

    TTree *mytree = s->Tree();
    TH1F *myth1f = s->TInputHist();

    int run=-1, ls=-1,evt=-1, nJet30=-1, nJet50=-1, nBtag_loose=-1, nBtag_medium=-1,nBtag_tight=-1;
    float nPV=-1, nPU_mean=-1, puWeight=1, puWeight_Up=1, puWeight_Dn=1, L1PFWeight=1, genWeight=1, btagWeight=1;
    float aqgcWeight[1000]={};
    float trig_eff_Weight=1;
    float lep1_pt=-1, lep1_eta=-999, lep1_phi=-1, lep1_m=-1, lep1_q=-1, lep1_iso=-1,lep1_idEffWeight=1, lep1_pt_scaleUp=-1, lep1_pt_scaleDn=-1, lep2_pt=-1, lep2_eta=-999, lep2_phi=-1, lep2_m=-1, lep2_q=-1, lep2_iso=-1,lep2_idEffWeight=1, lep2_pt_scaleUp=-1, lep2_pt_scaleDn=-1;
    float dilep_m=-1, dilep_pt=-1,dilep_eta=-999, dilep_phi=-1, dilep_m_scaleUp=-1, dilep_m_scaleDn=-1, dilep_pt_scaleUp=-1, dilep_pt_scaleDn=-1;
    float MET=-1, MET_phi=-1, MET_2017raw=-1, MET_scaleUp=-1, MET_scaleDn=-1, neu_pz_type0=-999, neu_pz_type0_scaleUp=-999, neu_pz_type0_scaleDn=-999;
    float vbf1_AK4_pt=-1, vbf1_AK4_eta=-999, vbf1_AK4_phi=-1, vbf1_AK4_m=-1, vbf1_AK4_gqid=-1, vbf1_AK4_axis2=-1, vbf1_AK4_ptD =-1, vbf1_AK4_pt_scaleUp=-1, vbf1_AK4_pt_scaleDn=-1, vbf1_AK4_m_scaleUp=-1, vbf1_AK4_m_scaleDn=-1, vbf2_AK4_pt=-1, vbf2_AK4_eta=-999, vbf2_AK4_phi=-1, vbf2_AK4_m=-1, vbf2_AK4_gqid=-1, vbf2_AK4_axis2=-1, vbf2_AK4_ptD =-1, vbf2_AK4_pt_scaleUp=-1, vbf2_AK4_pt_scaleDn=-1, vbf2_AK4_m_scaleUp=-1, vbf2_AK4_m_scaleDn=-1;
    float vbf_pt=-1, vbf_eta=-999, vbf_phi=-1, vbf_m=-1, vbf_pt_scaleUp=-1, vbf_pt_scaleDn=-1, vbf_m_scaleUp=-1, vbf_m_scaleDn=-1;
    float bos_PuppiAK8_m_sd0=-1, bos_PuppiAK8_m_sd0_corr=-1, bos_PuppiAK8_pt=-1, bos_PuppiAK8_eta=-999, bos_PuppiAK8_phi=-1, bos_PuppiAK8_tau2tau1=-999, bos_PuppiAK8_m_sd0_corr_scaleUp=-1, bos_PuppiAK8_m_sd0_corr_scaleDn=-1, bos_PuppiAK8_pt_scaleUp=-1, bos_PuppiAK8_pt_scaleDn=-1, bos_PuppiAK8_e2_sdb1=-999, bos_PuppiAK8_e3_sdb1=-999, bos_PuppiAK8_e3_v1_sdb1=-999, bos_PuppiAK8_e3_v2_sdb1=-999, bos_PuppiAK8_e4_v1_sdb1=-999, bos_PuppiAK8_e4_v2_sdb1=-999, bos_PuppiAK8_e2_sdb2=-999, bos_PuppiAK8_e3_sdb2=-999, bos_PuppiAK8_e3_v1_sdb2=-999, bos_PuppiAK8_e3_v2_sdb2=-999, bos_PuppiAK8_e4_v1_sdb2=-999, bos_PuppiAK8_e4_v2_sdb2=-999;
    float  bos_j1_AK4_pt=-1, bos_j1_AK4_eta=-999, bos_j1_AK4_phi=-1, bos_j1_AK4_m=-1, bos_j1_AK4_pt_scaleUp=-1, bos_j1_AK4_pt_scaleDn=-1, bos_j1_AK4_m_scaleUp=-1, bos_j1_AK4_m_scaleDn=-1, bos_j2_AK4_pt=-1, bos_j2_AK4_eta=-999, bos_j2_AK4_phi=-1, bos_j2_AK4_m=-1, bos_j2_AK4_pt_scaleUp=-1, bos_j2_AK4_pt_scaleDn=-1, bos_j2_AK4_m_scaleUp=-1, bos_j2_AK4_m_scaleDn=-1;
    float bos_AK4AK4_pt=-1, bos_AK4AK4_eta=-999, bos_AK4AK4_phi=-1, bos_AK4AK4_m=-1, bos_AK4AK4_pt_scaleUp=-1, bos_AK4AK4_pt_scaleDn=-1, bos_AK4AK4_m_scaleUp=-1, bos_AK4AK4_m_scaleDn=-1;
    float dibos_m=-1, dibos_pt=-1, dibos_eta=-999, dibos_phi=-1, dibos_m_scaleUp=-1, dibos_m_scaleDn=-1, dibos_pt_scaleUp=-1, dibos_pt_scaleDn=-1, bosCent=-999, zeppLep=-999, zeppHad=-999;

    mytree->SetBranchStatus("*",0);
    mytree->SetBranchStatus("run",1);
    mytree->SetBranchAddress("run",&run);
    mytree->SetBranchStatus("evt",1);
    mytree->SetBranchAddress("evt", & evt);
    mytree->SetBranchStatus("nJet30",1);
    mytree->SetBranchAddress("nJet30", & nJet30);
    mytree->SetBranchStatus("nJet50",1);
    mytree->SetBranchAddress("nJet50", & nJet50);
    mytree->SetBranchStatus("nBtag_loose",1);
    mytree->SetBranchAddress("nBtag_loose", & nBtag_loose);
    mytree->SetBranchStatus("nBtag_medium",1);
    mytree->SetBranchAddress("nBtag_medium", & nBtag_medium);
    mytree->SetBranchStatus("nBtag_tight",1);
    mytree->SetBranchAddress("nBtag_tight", & nBtag_tight);
    mytree->SetBranchStatus("nPV",1);
    mytree->SetBranchAddress("nPV", & nPV);
    mytree->SetBranchStatus("nPU_mean",1);
    mytree->SetBranchAddress("nPU_mean", & nPU_mean);
    mytree->SetBranchStatus("puWeight",1);
    mytree->SetBranchAddress("puWeight", & puWeight);
    mytree->SetBranchStatus("puWeight_Up",1);
    mytree->SetBranchAddress("puWeight_Up", & puWeight_Up);
    mytree->SetBranchStatus("puWeight_Dn",1);
    mytree->SetBranchAddress("puWeight_Dn", & puWeight_Dn);
    mytree->SetBranchStatus("genWeight",1);
    mytree->SetBranchAddress("genWeight", & genWeight);
    mytree->SetBranchStatus("btagWeight",1);
    mytree->SetBranchAddress("btagWeight", & btagWeight);
    mytree->SetBranchStatus("aqgcWeight",1);
    mytree->SetBranchAddress("aqgcWeight", & aqgcWeight);
    mytree->SetBranchStatus("L1PFWeight",1);
    mytree->SetBranchAddress("L1PFWeight", & L1PFWeight);
    mytree->SetBranchStatus("lep1_pt",1);
    mytree->SetBranchAddress("lep1_pt", & lep1_pt);
    mytree->SetBranchStatus("lep1_eta",1);
    mytree->SetBranchAddress("lep1_eta", & lep1_eta);
    mytree->SetBranchStatus("lep1_phi",1);
    mytree->SetBranchAddress("lep1_phi", & lep1_phi);
    mytree->SetBranchStatus("lep1_m",1);
    mytree->SetBranchAddress("lep1_m", & lep1_m);
    mytree->SetBranchStatus("lep1_q",1);
    mytree->SetBranchAddress("lep1_q", & lep1_q);
    mytree->SetBranchStatus("lep1_iso",1);
    mytree->SetBranchAddress("lep1_iso", & lep1_iso);
    mytree->SetBranchStatus("lep1_idEffWeight",1);
    mytree->SetBranchAddress("lep1_idEffWeight", & lep1_idEffWeight);
    mytree->SetBranchStatus("lep1_pt_scaleUp",1);
    mytree->SetBranchAddress("lep1_pt_scaleUp", & lep1_pt_scaleUp);
    mytree->SetBranchStatus("lep1_pt_scaleDn",1);
    mytree->SetBranchAddress("lep1_pt_scaleDn", & lep1_pt_scaleDn);
    mytree->SetBranchStatus("lep2_pt",1);
    mytree->SetBranchAddress("lep2_pt", & lep2_pt);
    mytree->SetBranchStatus("lep2_eta",1);
    mytree->SetBranchAddress("lep2_eta", & lep2_eta);
    mytree->SetBranchStatus("lep2_phi",1);
    mytree->SetBranchAddress("lep2_phi", & lep2_phi);
    mytree->SetBranchStatus("lep2_m",1);
    mytree->SetBranchAddress("lep2_m", & lep2_m);
    mytree->SetBranchStatus("lep2_q",1);
    mytree->SetBranchAddress("lep2_q", & lep2_q);
    mytree->SetBranchStatus("lep2_iso",1);
    mytree->SetBranchAddress("lep2_iso", & lep2_iso);
    mytree->SetBranchStatus("lep2_idEffWeight",1);
    mytree->SetBranchAddress("lep2_idEffWeight", & lep2_idEffWeight);
    mytree->SetBranchStatus("lep2_pt_scaleUp",1);
    mytree->SetBranchAddress("lep2_pt_scaleUp", & lep2_pt_scaleUp);
    mytree->SetBranchStatus("lep2_pt_scaleDn",1);
    mytree->SetBranchAddress("lep2_pt_scaleDn", & lep2_pt_scaleDn);
    mytree->SetBranchStatus("dilep_m",1);
    mytree->SetBranchAddress("dilep_m", & dilep_m);
    mytree->SetBranchStatus("dilep_pt",1);
    mytree->SetBranchAddress("dilep_pt", & dilep_pt);
    mytree->SetBranchStatus("dilep_eta",1);
    mytree->SetBranchAddress("dilep_eta", & dilep_eta);
    mytree->SetBranchStatus("dilep_phi",1);
    mytree->SetBranchAddress("dilep_phi", & dilep_phi);
    mytree->SetBranchStatus("dilep_m_scaleUp",1);
    mytree->SetBranchAddress("dilep_m_scaleUp", & dilep_m_scaleUp);
    mytree->SetBranchStatus("dilep_m_scaleDn",1);
    mytree->SetBranchAddress("dilep_m_scaleDn", & dilep_m_scaleDn);
    mytree->SetBranchStatus("dilep_pt_scaleUp",1);
    mytree->SetBranchAddress("dilep_pt_scaleUp", & dilep_pt_scaleUp);
    mytree->SetBranchStatus("dilep_pt_scaleDn",1);
    mytree->SetBranchAddress("dilep_pt_scaleDn", & dilep_pt_scaleDn);
    mytree->SetBranchStatus("MET",1);
    mytree->SetBranchAddress("MET", & MET);
    mytree->SetBranchStatus("MET_phi",1);
    mytree->SetBranchAddress("MET_phi", & MET_phi);
    mytree->SetBranchStatus("MET_2017raw",1);
    mytree->SetBranchAddress("MET_2017raw", & MET_2017raw);
    mytree->SetBranchStatus("MET_scaleUp",1);
    mytree->SetBranchAddress("MET_scaleUp", & MET_scaleUp);
    mytree->SetBranchStatus("MET_scaleDn",1);
    mytree->SetBranchAddress("MET_scaleDn", & MET_scaleDn);
    mytree->SetBranchStatus("neu_pz_type0",1);
    mytree->SetBranchAddress("neu_pz_type0", & neu_pz_type0);
    mytree->SetBranchStatus("neu_pz_type0_scaleUp",1);
    mytree->SetBranchAddress("neu_pz_type0_scaleUp", & neu_pz_type0_scaleUp);
    mytree->SetBranchStatus("neu_pz_type0_scaleDn",1);
    mytree->SetBranchAddress("neu_pz_type0_scaleDn", & neu_pz_type0_scaleDn);
    mytree->SetBranchStatus("vbf1_AK4_pt",1);
    mytree->SetBranchAddress("vbf1_AK4_pt", & vbf1_AK4_pt);
    mytree->SetBranchStatus("vbf1_AK4_eta",1);
    mytree->SetBranchAddress("vbf1_AK4_eta", & vbf1_AK4_eta);
    mytree->SetBranchStatus("vbf1_AK4_phi",1);
    mytree->SetBranchAddress("vbf1_AK4_phi", & vbf1_AK4_phi);
    mytree->SetBranchStatus("vbf1_AK4_m",1);
    mytree->SetBranchAddress("vbf1_AK4_m", & vbf1_AK4_m);
    mytree->SetBranchStatus("vbf1_AK4_gqid",1);
    mytree->SetBranchAddress("vbf1_AK4_gqid", & vbf1_AK4_gqid);
    mytree->SetBranchStatus("vbf1_AK4_axis2",1);
    mytree->SetBranchAddress("vbf1_AK4_axis2", & vbf1_AK4_axis2);
    mytree->SetBranchStatus("vbf1_AK4_ptD",1);
    mytree->SetBranchAddress("vbf1_AK4_ptD", & vbf1_AK4_ptD);
    mytree->SetBranchStatus("vbf1_AK4_pt_scaleUp",1);
    mytree->SetBranchAddress("vbf1_AK4_pt_scaleUp", & vbf1_AK4_pt_scaleUp);
    mytree->SetBranchStatus("vbf1_AK4_pt_scaleDn",1);
    mytree->SetBranchAddress("vbf1_AK4_pt_scaleDn", & vbf1_AK4_pt_scaleDn);
    mytree->SetBranchStatus("vbf1_AK4_m_scaleUp",1);
    mytree->SetBranchAddress("vbf1_AK4_m_scaleUp", & vbf1_AK4_m_scaleUp);
    mytree->SetBranchStatus("vbf1_AK4_m_scaleDn",1);
    mytree->SetBranchAddress("vbf1_AK4_m_scaleDn", & vbf1_AK4_m_scaleDn);
    mytree->SetBranchStatus("vbf2_AK4_pt",1);
    mytree->SetBranchAddress("vbf2_AK4_pt", & vbf2_AK4_pt);
    mytree->SetBranchStatus("vbf2_AK4_eta",1);
    mytree->SetBranchAddress("vbf2_AK4_eta", & vbf2_AK4_eta);
    mytree->SetBranchStatus("vbf2_AK4_phi",1);
    mytree->SetBranchAddress("vbf2_AK4_phi", & vbf2_AK4_phi);
    mytree->SetBranchStatus("vbf2_AK4_m",1);
    mytree->SetBranchAddress("vbf2_AK4_m", & vbf2_AK4_m);
    mytree->SetBranchStatus("vbf2_AK4_gqid",1);
    mytree->SetBranchAddress("vbf2_AK4_gqid", & vbf2_AK4_gqid);
    mytree->SetBranchStatus("vbf2_AK4_axis2",1);
    mytree->SetBranchAddress("vbf2_AK4_axis2", & vbf2_AK4_axis2);
    mytree->SetBranchStatus("vbf2_AK4_ptD",1);
    mytree->SetBranchAddress("vbf2_AK4_ptD", & vbf2_AK4_ptD);
    mytree->SetBranchStatus("vbf2_AK4_pt_scaleUp",1);
    mytree->SetBranchAddress("vbf2_AK4_pt_scaleUp", & vbf2_AK4_pt_scaleUp);
    mytree->SetBranchStatus("vbf2_AK4_pt_scaleDn",1);
    mytree->SetBranchAddress("vbf2_AK4_pt_scaleDn", & vbf2_AK4_pt_scaleDn);
    mytree->SetBranchStatus("vbf2_AK4_m_scaleUp",1);
    mytree->SetBranchAddress("vbf2_AK4_m_scaleUp", & vbf2_AK4_m_scaleUp);
    mytree->SetBranchStatus("vbf2_AK4_m_scaleDn",1);
    mytree->SetBranchAddress("vbf2_AK4_m_scaleDn", & vbf2_AK4_m_scaleDn);
    mytree->SetBranchStatus("vbf_pt",1);
    mytree->SetBranchAddress("vbf_pt", & vbf_pt);
    mytree->SetBranchStatus("vbf_eta",1);
    mytree->SetBranchAddress("vbf_eta", & vbf_eta);
    mytree->SetBranchStatus("vbf_phi",1);
    mytree->SetBranchAddress("vbf_phi", & vbf_phi);
    mytree->SetBranchStatus("vbf_m",1);
    mytree->SetBranchAddress("vbf_m", & vbf_m);
    mytree->SetBranchStatus("vbf_pt_scaleUp",1);
    mytree->SetBranchAddress("vbf_pt_scaleUp", & vbf_pt_scaleUp);
    mytree->SetBranchStatus("vbf_pt_scaleDn",1);
    mytree->SetBranchAddress("vbf_pt_scaleDn", & vbf_pt_scaleDn);
    mytree->SetBranchStatus("vbf_m_scaleUp",1);
    mytree->SetBranchAddress("vbf_m_scaleUp", & vbf_m_scaleUp);
    mytree->SetBranchStatus("vbf_m_scaleDn",1);
    mytree->SetBranchAddress("vbf_m_scaleDn", & vbf_m_scaleDn);
    mytree->SetBranchStatus("bos_PuppiAK8_m_sd0",1);
    mytree->SetBranchAddress("bos_PuppiAK8_m_sd0", & bos_PuppiAK8_m_sd0);
    mytree->SetBranchStatus("bos_PuppiAK8_m_sd0_corr",1);
    mytree->SetBranchAddress("bos_PuppiAK8_m_sd0_corr", & bos_PuppiAK8_m_sd0_corr);
    mytree->SetBranchStatus("bos_PuppiAK8_pt",1);
    mytree->SetBranchAddress("bos_PuppiAK8_pt", & bos_PuppiAK8_pt);
    mytree->SetBranchStatus("bos_PuppiAK8_eta",1);
    mytree->SetBranchAddress("bos_PuppiAK8_eta", & bos_PuppiAK8_eta);
    mytree->SetBranchStatus("bos_PuppiAK8_phi",1);
    mytree->SetBranchAddress("bos_PuppiAK8_phi", & bos_PuppiAK8_phi);
    mytree->SetBranchStatus("bos_PuppiAK8_tau2tau1",1);
    mytree->SetBranchAddress("bos_PuppiAK8_tau2tau1", & bos_PuppiAK8_tau2tau1);
    mytree->SetBranchStatus("bos_PuppiAK8_m_sd0_corr_scaleUp",1);
    mytree->SetBranchAddress("bos_PuppiAK8_m_sd0_corr_scaleUp", & bos_PuppiAK8_m_sd0_corr_scaleUp);
    mytree->SetBranchStatus("bos_PuppiAK8_m_sd0_corr_scaleDn",1);
    mytree->SetBranchAddress("bos_PuppiAK8_m_sd0_corr_scaleDn", & bos_PuppiAK8_m_sd0_corr_scaleDn);
    mytree->SetBranchStatus("bos_PuppiAK8_pt_scaleUp",1);
    mytree->SetBranchAddress("bos_PuppiAK8_pt_scaleUp", & bos_PuppiAK8_pt_scaleUp);
    mytree->SetBranchStatus("bos_PuppiAK8_pt_scaleDn",1);
    mytree->SetBranchAddress("bos_PuppiAK8_pt_scaleDn", & bos_PuppiAK8_pt_scaleDn);
    mytree->SetBranchStatus("bos_PuppiAK8_e2_sdb1",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e2_sdb1", & bos_PuppiAK8_e2_sdb1);
    mytree->SetBranchStatus("bos_PuppiAK8_e3_sdb1",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e3_sdb1", & bos_PuppiAK8_e3_sdb1);
    mytree->SetBranchStatus("bos_PuppiAK8_e3_v1_sdb1",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e3_v1_sdb1", & bos_PuppiAK8_e3_v1_sdb1);
    mytree->SetBranchStatus("bos_PuppiAK8_e3_v2_sdb1",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e3_v2_sdb1", & bos_PuppiAK8_e3_v2_sdb1);
    mytree->SetBranchStatus("bos_PuppiAK8_e4_v1_sdb1",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e4_v1_sdb1", & bos_PuppiAK8_e4_v1_sdb1);
    mytree->SetBranchStatus("bos_PuppiAK8_e4_v2_sdb1",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e4_v2_sdb1", & bos_PuppiAK8_e4_v2_sdb1);
    mytree->SetBranchStatus("bos_PuppiAK8_e2_sdb2",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e2_sdb2", & bos_PuppiAK8_e2_sdb2);
    mytree->SetBranchStatus("bos_PuppiAK8_e3_sdb2",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e3_sdb2", & bos_PuppiAK8_e3_sdb2);
    mytree->SetBranchStatus("bos_PuppiAK8_e3_v1_sdb2",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e3_v1_sdb2", & bos_PuppiAK8_e3_v1_sdb2);
    mytree->SetBranchStatus("bos_PuppiAK8_e3_v2_sdb2",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e3_v2_sdb2", & bos_PuppiAK8_e3_v2_sdb2);
    mytree->SetBranchStatus("bos_PuppiAK8_e4_v1_sdb2",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e4_v1_sdb2", & bos_PuppiAK8_e4_v1_sdb2);
    mytree->SetBranchStatus("bos_PuppiAK8_e4_v2_sdb2",1);
    mytree->SetBranchAddress("bos_PuppiAK8_e4_v2_sdb2", & bos_PuppiAK8_e4_v2_sdb2);
    mytree->SetBranchStatus("bos_j1_AK4_pt",1);
    mytree->SetBranchAddress("bos_j1_AK4_pt", & bos_j1_AK4_pt);
    mytree->SetBranchStatus("bos_j1_AK4_eta",1);
    mytree->SetBranchAddress("bos_j1_AK4_eta", & bos_j1_AK4_eta);
    mytree->SetBranchStatus("bos_j1_AK4_phi",1);
    mytree->SetBranchAddress("bos_j1_AK4_phi", & bos_j1_AK4_phi);
    mytree->SetBranchStatus("bos_j1_AK4_m",1);
    mytree->SetBranchAddress("bos_j1_AK4_m", & bos_j1_AK4_m);
    mytree->SetBranchStatus("bos_j1_AK4_pt_scaleUp",1);
    mytree->SetBranchAddress("bos_j1_AK4_pt_scaleUp", & bos_j1_AK4_pt_scaleUp);
    mytree->SetBranchStatus("bos_j1_AK4_pt_scaleDn",1);
    mytree->SetBranchAddress("bos_j1_AK4_pt_scaleDn", & bos_j1_AK4_pt_scaleDn);
    mytree->SetBranchStatus("bos_j1_AK4_m_scaleUp",1);
    mytree->SetBranchAddress("bos_j1_AK4_m_scaleUp", & bos_j1_AK4_m_scaleUp);
    mytree->SetBranchStatus("bos_j1_AK4_m_scaleDn",1);
    mytree->SetBranchAddress("bos_j1_AK4_m_scaleDn", & bos_j1_AK4_m_scaleDn);
    mytree->SetBranchStatus("bos_j2_AK4_pt",1);
    mytree->SetBranchAddress("bos_j2_AK4_pt", & bos_j2_AK4_pt);
    mytree->SetBranchStatus("bos_j2_AK4_eta",1);
    mytree->SetBranchAddress("bos_j2_AK4_eta", & bos_j2_AK4_eta);
    mytree->SetBranchStatus("bos_j2_AK4_phi",1);
    mytree->SetBranchAddress("bos_j2_AK4_phi", & bos_j2_AK4_phi);
    mytree->SetBranchStatus("bos_j2_AK4_m",1);
    mytree->SetBranchAddress("bos_j2_AK4_m", & bos_j2_AK4_m);
    mytree->SetBranchStatus("bos_j2_AK4_pt_scaleUp",1);
    mytree->SetBranchAddress("bos_j2_AK4_pt_scaleUp", & bos_j2_AK4_pt_scaleUp);
    mytree->SetBranchStatus("bos_j2_AK4_pt_scaleDn",1);
    mytree->SetBranchAddress("bos_j2_AK4_pt_scaleDn", & bos_j2_AK4_pt_scaleDn);
    mytree->SetBranchStatus("bos_j2_AK4_m_scaleUp",1);
    mytree->SetBranchAddress("bos_j2_AK4_m_scaleUp", & bos_j2_AK4_m_scaleUp);
    mytree->SetBranchStatus("bos_j2_AK4_m_scaleDn",1);
    mytree->SetBranchAddress("bos_j2_AK4_m_scaleDn", & bos_j2_AK4_m_scaleDn);
    mytree->SetBranchStatus("bos_AK4AK4_pt",1);
    mytree->SetBranchAddress("bos_AK4AK4_pt", & bos_AK4AK4_pt);
    mytree->SetBranchStatus("bos_AK4AK4_eta",1);
    mytree->SetBranchAddress("bos_AK4AK4_eta", & bos_AK4AK4_eta);
    mytree->SetBranchStatus("bos_AK4AK4_phi",1);
    mytree->SetBranchAddress("bos_AK4AK4_phi", & bos_AK4AK4_phi);
    mytree->SetBranchStatus("bos_AK4AK4_m",1);
    mytree->SetBranchAddress("bos_AK4AK4_m", & bos_AK4AK4_m);
    mytree->SetBranchStatus("bos_AK4AK4_pt_scaleUp",1);
    mytree->SetBranchAddress("bos_AK4AK4_pt_scaleUp", & bos_AK4AK4_pt_scaleUp);
    mytree->SetBranchStatus("bos_AK4AK4_pt_scaleDn",1);
    mytree->SetBranchAddress("bos_AK4AK4_pt_scaleDn", & bos_AK4AK4_pt_scaleDn);
    mytree->SetBranchStatus("bos_AK4AK4_m_scaleUp",1);
    mytree->SetBranchAddress("bos_AK4AK4_m_scaleUp", & bos_AK4AK4_m_scaleUp);
    mytree->SetBranchStatus("bos_AK4AK4_m_scaleDn",1);
    mytree->SetBranchAddress("bos_AK4AK4_m_scaleDn", & bos_AK4AK4_m_scaleDn);
    mytree->SetBranchStatus("dibos_m",1);
    mytree->SetBranchAddress("dibos_m", & dibos_m);
    mytree->SetBranchStatus("dibos_pt",1);
    mytree->SetBranchAddress("dibos_pt", & dibos_pt);
    mytree->SetBranchStatus("dibos_eta",1);
    mytree->SetBranchAddress("dibos_eta", & dibos_eta);
    mytree->SetBranchStatus("dibos_phi",1);
    mytree->SetBranchAddress("dibos_phi", & dibos_phi);
    mytree->SetBranchStatus("dibos_m_scaleUp",1);
    mytree->SetBranchAddress("dibos_m_scaleUp", & dibos_m_scaleUp);
    mytree->SetBranchStatus("dibos_m_scaleDn",1);
    mytree->SetBranchAddress("dibos_m_scaleDn", & dibos_m_scaleDn);
    mytree->SetBranchStatus("dibos_pt_scaleUp",1);
    mytree->SetBranchAddress("dibos_pt_scaleUp", & dibos_pt_scaleUp);
    mytree->SetBranchStatus("dibos_pt_scaleDn",1);
    mytree->SetBranchAddress("dibos_pt_scaleDn", & dibos_pt_scaleDn);
    mytree->SetBranchStatus("bosCent",1);
    mytree->SetBranchAddress("bosCent", & bosCent);
    mytree->SetBranchStatus("zeppLep",1);
    mytree->SetBranchAddress("zeppLep", & zeppLep);
    mytree->SetBranchStatus("zeppHad",1);
    mytree->SetBranchAddress("zeppHad", & zeppHad);
	///************************************************/////*
    int nEvents=-1, nNegEvents=-1, type=-1, nBTagJet_loose=-1;
   //*************************************************// 
        float nTotal=0, nNeg=0;
	if (!(isamp==0)){
	nTotal = myth1f->GetBinContent(2);
        nNeg = myth1f->GetBinContent(1);
	}
	cout<<nTotal<<"   "<<nNeg<<endl;
    for(int i = 0; i<mytree->GetEntries(); i++)
    {
//if(s->name().EqualTo("data"))      {cout<<i<<" "<<dibos_m<<endl;} 
//	if(i>10000) continue;
      mytree->GetEntry(i);
//	if (i%500==0)cout <<mytree->GetEntry(i)<<"   "<<lep1_m<<"  "<<lep2_m<<endl;
      bool isEle=false, isResolved=false, isZ=false;

      if (bos_PuppiAK8_m_sd0_corr > 0 && bos_AK4AK4_m < 0) { isResolved=false; }
      else if (bos_PuppiAK8_m_sd0_corr < 0 && bos_AK4AK4_m > 0) { isResolved=true; }
      else {
        //cout << "both or neither of resolved and boosted mass is defined" << endl;
                continue;
                      }

      if (lep1_m == ELE_MASS){isEle=true;}
      else if (lep1_m == MUON_MASS){isEle=false;}
      else {
        cout << "lepton is not electron or muon! skipping" << endl;
        continue;
      }
      if (fabs(vbf1_AK4_eta)>2.65 && fabs(vbf1_AK4_eta)<3.139) continue;
      if (fabs(vbf2_AK4_eta)>2.65 && fabs(vbf2_AK4_eta)<3.139) continue;

      if ( vbf_m < 500) continue;
      if ( fabs(vbf1_AK4_eta - vbf2_AK4_eta)<2.5) continue;
//wjets
      if ( !(nBtag_loose==0 && vbf1_AK4_pt>50 && vbf2_AK4_pt>50) ) continue;
      if (isResolved==true && (bos_AK4AK4_m>65 &&bos_AK4AK4_m<105)) continue;
      if (isResolved==false && (bos_PuppiAK8_m_sd0_corr>65 &&bos_PuppiAK8_m_sd0_corr<105)) continue;

      if (isResolved==true && (bos_j1_AK4_pt<50 || bos_j2_AK4_pt<50)) continue;

      if (isResolved==false && bos_PuppiAK8_pt<200) continue;
      if (isResolved==false && bos_PuppiAK8_tau2tau1>0.55) continue;


      if (lep2_pt>0) isZ=true;

      if (isEle==true && (lep1_pt<35 || abs(lep1_eta)>2.5 || (abs(lep1_eta)>1.4442 && abs(lep1_eta)<1.566))) continue;
      if (isEle==false && (lep1_pt<35 || abs(lep1_eta)>2.4)) continue;

      if (isZ==true && (dilep_m < 81 || dilep_m > 101)) continue;
      if (isZ==true && isEle==true && (lep2_pt<20 || abs(lep2_eta)>2.5 || (abs(lep2_eta)>1.4442 && abs(lep2_eta)<1.566))) continue;
      if (isZ==true && isEle==false && (lep2_pt<20 || abs(lep2_eta)>2.4)) continue;
      if (isZ==true && (lep1_q*lep2_q)==1) continue;

      if (isZ==false && MET<30) continue;
      //if(!(type==0||type==1)) continue;
      
      if (1)	//----------------	Nominal, PU up, PU down
        {
	if(i<3){	cout<<"xsec="<<xsec<<": otherscale="<<otherscale<<": genWeight="<<genWeight<<"; trig_eff="<<trig_eff_Weight<<"; id_eff="<<lep1_idEffWeight<<";  puWeight="<<puWeight<<"; btag="<<btagWeight<<endl;
	//double weightC = (xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg));
	//cout<<"dibos_m= "<<dibos_m<<"  weight="<<weightC<<endl;
	}
	if (isResolved==true && isZ==false)
	    {
	double weightC = (xsec*genWeight*otherscale)/(1.0*(nTotal-2*nNeg));
	//cout<<"dibos_m= "<<dibos_m<<"  weight="<<weightC<<endl;
	
	      if(s->name().EqualTo("data"))	 {hists[0]->Fill(dibos_m); h1->Fill(dibos_m);}
	      if(s->name().EqualTo("WV_EWK"))	 {hists[1]->Fill(dibos_m,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));}// cout<<"dibos_m= "<<dibos_m<<"  weight="<<weightC<<endl; }
	      if(s->name().EqualTo("Diboson")) 	 hists[14]->Fill(dibos_m,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("top"))  	 hists[27]->Fill(dibos_m,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[40]->Fill(dibos_m,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("CH_WZ"))	 hists[53]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("DCH_WW"))	 hists[66]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
  	    /*  HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		    { 
		      //TString name = HiggsSampleName[i]+MassPoint[j];
		      TString OrigName = s->name();
		      TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		      if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		      HistCount++;
		//	cout<< OrigName<<"*******************"<<name<<endl;
		    }
	      
	      //------	PU UP
	      if(s->name().EqualTo("WV_EWK"))	 hists[8]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Up*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("Diboson")) 	 hists[21]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Up*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("top"))  	 hists[34]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Up*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[47]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Up*btagWeight)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("CH_WZ"))	 hists[60]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Up*btagWeight)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("DCH_WW"))	 hists[73]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Up*btagWeight)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		    { 
		      //TString name = HiggsSampleName[i]+MassPoint[j];
		      TString OrigName = s->name()+"_CMS_puUp";
		      TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		      if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Up*btagWeight)/(1.0*(nTotal-2*nNeg)));
		      HistCount++;
		    }
	      
	      //------	PU Down
	      if(s->name().EqualTo("WV_EWK"))	 hists[9]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Dn*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("Diboson")) 	 hists[22]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Dn*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("top"))  	 hists[35]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Dn*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[48]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Dn*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("CH_WZ"))	 hists[61]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Dn*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("DCH_WW"))	 hists[74]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Dn*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		    { 
		      //TString name = HiggsSampleName[i]+MassPoint[j];
		      TString OrigName = s->name()+"_CMS_puDown";
		      TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		      if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight_Dn*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
		      HistCount++;
		    }
	      
	      //------	btag HF Up
	      if(s->name().EqualTo("WV_EWK"))	 hists[10]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("Diboson")) 	 hists[23]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("top"))  	 hists[36]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[49]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("CH_WZ"))	 hists[62]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("DCH_WW"))	 hists[75]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		    { 
		      //TString name = HiggsSampleName[i]+MassPoint[j];
		      TString OrigName = s->name()+"_CMS_btagHFUp";
		      TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		      if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpHF)/(1.0*(nTotal-2*nNeg)));
		      HistCount++;
		    }

	      //------	btag HF Down
	      if(s->name().EqualTo("WV_EWK"))	 hists[11]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownHF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("Diboson")) 	 hists[24]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownHF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("top"))  	 hists[37]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownHF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[50]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownHF)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("CH_WZ"))	 hists[63]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownHF)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("DCH_WW"))	 hists[76]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownHF)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		    { 
		      //TString name = HiggsSampleName[i]+MassPoint[j];
		      TString OrigName = s->name()+"_CMS_btagHFDown";
		      TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		      if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownHF)/(1.0*(nTotal-2*nNeg)));
		      HistCount++;
		    }
	      
	      //------	btag LF Up
	      if(s->name().EqualTo("WV_EWK"))	 hists[12]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpLF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("Diboson")) 	 hists[25]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpLF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("top"))  	 hists[38]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpLF)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[51]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpLF)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("CH_WZ"))	 hists[64]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpLF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("DCH_WW"))	 hists[77]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpLF)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		    { 
		      //TString name = HiggsSampleName[i]+MassPoint[j];
		      TString OrigName = s->name()+"_CMS_btagLFUp";
		      TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		      if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtUpLF)/(1.0*(nTotal-2*nNeg)));
		      HistCount++;
		    }

	      //------	btag LF Down
	      if(s->name().EqualTo("WV_EWK"))	 hists[13]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownLF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("Diboson")) 	 hists[26]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownLF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("top"))  	 hists[39]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownLF)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[52]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownLF)/(1.0*(nTotal-2*nNeg)));// cout<<"Vjets works!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
//	      if(s->name().EqualTo("CH_WZ"))	 hists[65]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownLF)/(1.0*(nTotal-2*nNeg)));
//	      if(s->name().EqualTo("DCH_WW"))	 hists[78]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownLF)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		    { 
		      //TString name = HiggsSampleName[i]+MassPoint[j];
		      TString OrigName = s->name()+"_CMS_btagLFDown";
		      TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		      if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(bos_PuppiAK8_m_sd0,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btag0WgtDownLF)/(1.0*(nTotal-2*nNeg)));
		      HistCount++;
		    }
	      
	      // To get QCD scale bounding we need to add QCD scale for all signal and bkg. But except for WV_EWK and Diboson others are taken care of using background estimation. For top there is not QCD scale bounding present in MC.
	      if(s->name().EqualTo("WV_EWK")||s->name().EqualTo("Diboson"))
	  	{
		  if(s->name().EqualTo("WV_EWK"))
		    {
		    histo_diboson_EWK_CMS_QCDScaleBounding[0]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[1]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_diboson_EWK_CMS_QCDScaleBounding[1]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[2]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_diboson_EWK_CMS_QCDScaleBounding[2]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[3]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_diboson_EWK_CMS_QCDScaleBounding[3]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[4]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_diboson_EWK_CMS_QCDScaleBounding[4]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[6]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_diboson_EWK_CMS_QCDScaleBounding[5]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[8]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    for(int npdf=0; npdf<100; npdf++) histo_diboson_EWK_CMS_PDFScaleBounding[npdf]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[9+npdf]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    }
		  if(s->name().EqualTo("Diboson"))
		    {
		    histo_VVjjQCD_EWK_CMS_QCDScaleBounding[0]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[1]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_VVjjQCD_EWK_CMS_QCDScaleBounding[1]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[2]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_VVjjQCD_EWK_CMS_QCDScaleBounding[2]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[3]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_VVjjQCD_EWK_CMS_QCDScaleBounding[3]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[4]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_VVjjQCD_EWK_CMS_QCDScaleBounding[4]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[6]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    histo_VVjjQCD_EWK_CMS_QCDScaleBounding[5]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[8]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    for(int npdf=0; npdf<100; npdf++) histo_VVjjQCD_EWK_CMS_PDFScaleBounding[npdf]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[9+npdf]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    }
		 }
	      HistCount = 0;
	      int LHEWgt[6] = {1, 2, 3, 4, 6, 8};
	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<6; k++)
		    { 
		      //TString name = HiggsSampleName[i]+MassPoint[j];
		      TString OrigName = s->name()+"_CMS_QCDScale";
		      TString name = HiggsSampleName[i]+MassPoint[j]+"_CMS_QCDScale";
		      if (OrigName.EqualTo(name)) 
			{
			  //std::cout << (LHEWeight[LHEWgt[k]]/LHEWeight[0]) << std::endl;
			  ChargedHistQCD[HistCount]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[LHEWgt[k]]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
			}
		      HistCount++;
		    }
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<100; k++)
		    { 
		      //TString name = HiggsSampleName[i]+MassPoint[j];
		      TString OrigName = s->name()+"_CMS_PDFScale";
		      TString name = HiggsSampleName[i]+MassPoint[j]+"_CMS_PDFScale";
		      if (OrigName.EqualTo(name)) ChargedHistPDF[HistCount]->Fill(bos_PuppiAK8_m_sd0,((LHEWeight[9+k]/LHEWeight[0])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		      HistCount++;
		    }*/
	      if(s->name().EqualTo("aQGC"))
	  	{
		  for (int j=0;j<678;j++)
		    {
//		      histo_aqgc[j]->Fill(dibos_m,((aqgcWeight[j])*xsec*otherscale*genWeight)/(1.0*(nTotal-2*nNeg)));
		      histo_aqgc[j]->Fill(dibos_m,((aqgcWeight[j])*xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    }
		}
	    }
	}
/*	  
      if (1)	//--------------------------- LEP up
	    {
	      if (
		  (l_pt2_Up<0 && l_pt1_Up>50) &&
		  (((type==0)&&(abs(l_eta1)<2.4))||((type==1)&&((abs(l_eta1)<2.5)&&!(abs(l_eta1)>1.4442 && abs(l_eta1)<1.566)))) &&
		  (((type==0)&&(pfMET_Corr>50)) || ((type==1)&&(pfMET_Corr>80))) &&
		  ((bos_PuppiAK8_pt>200)&&(abs(ungroomed_PuppiAK8_jet_eta)<2.4)&&(PuppiAK8_jet_tau2tau1<0.55)) &&
		  ((PuppiAK8_jet_mass_so_corr>65) && (PuppiAK8_jet_mass_so_corr<105)) &&
		  (nBTagJet_loose==0) &&
		  (vbf_maxpt_jj_m>800) &&
		  (abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta)>4.0) &&
		  ((vbf_maxpt_j1_pt>30) && (vbf_maxpt_j2_pt>30)) &&
		  (mass_lvj_type0_LEP_Up>600) &&
		  (BosonCentrality_type0_LEP_Up>1.0) &&
		  ((abs(ZeppenfeldWL_type0_LEP_Up)/abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta))<0.3) &&
		  ((abs(ZeppenfeldWH)/abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta))<0.3)
		  ){
		//if(s->name().EqualTo("data"))	 histo_data_LEPUp->Fill(mass_lvj_type0_LEP_Up);
		if(s->name().EqualTo("WV_EWK"))	 hists[2]->Fill(mass_lvj_type0_LEP_Up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		if(s->name().EqualTo("Diboson"))   hists[15]->Fill(mass_lvj_type0_LEP_Up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		if(s->name().EqualTo("top"))  	 hists[28]->Fill(mass_lvj_type0_LEP_Up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[41]->Fill(mass_lvj_type0_LEP_Up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		if(s->name().EqualTo("CH_WZ"))	 hists[54]->Fill(mass_lvj_type0_LEP_Up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		if(s->name().EqualTo("DCH_WW"))	 hists[67]->Fill(mass_lvj_type0_LEP_Up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		HistCount = 0;
		for (int i=0; i<3; i++)
		  for (int j=0; j<11; j++)
		    for (int k=0; k<17; k++)
		      { 
			//TString name = HiggsSampleName[i]+MassPoint[j];
			TString OrigName = s->name()+"_CMS_scale_lUp";
			TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
			if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(mass_lvj_type0_LEP_Up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
			HistCount++;
		  }
	      }
	    }*/
	  
/*	  if (1)	//--------------------------- LEP down
      {
          if (
	      (l_pt2_Down<0 && l_pt1_Down>50) &&
	      (((type==0)&&(abs(l_eta1)<2.4))||((type==1)&&((abs(l_eta1)<2.5)&&!(abs(l_eta1)>1.4442 && abs(l_eta1)<1.566)))) &&
	      (((type==0)&&(pfMET_Corr>50)) || ((type==1)&&(pfMET_Corr>80))) &&
	      ((bos_PuppiAK8_pt>200)&&(abs(ungroomed_PuppiAK8_jet_eta)<2.4)&&(PuppiAK8_jet_tau2tau1<0.55)) &&
	      ((PuppiAK8_jet_mass_so_corr>65) && (PuppiAK8_jet_mass_so_corr<105)) &&
	      (nBTagJet_loose==0) &&
	      (vbf_maxpt_jj_m>800) &&
	      (abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta)>4.0) &&
	      ((vbf_maxpt_j1_pt>30) && (vbf_maxpt_j2_pt>30)) &&
	      (mass_lvj_type0_LEP_Down>600) &&
	      (BosonCentrality_type0_LEP_Down>1.0) &&
	      ((abs(ZeppenfeldWL_type0_LEP_Down)/abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta))<0.3) &&
	      ((abs(ZeppenfeldWH)/abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta))<0.3) 
          ){
	      //if(s->name().EqualTo("data"))	 histo_data_LEPDown->Fill(mass_lvj_type0_LEP_Down);
	      if(s->name().EqualTo("WV_EWK"))	 hists[3]->Fill(mass_lvj_type0_LEP_Down,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("Diboson"))   hists[16]->Fill(mass_lvj_type0_LEP_Down,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("top"))  	 hists[29]->Fill(mass_lvj_type0_LEP_Down,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[42]->Fill(mass_lvj_type0_LEP_Down,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("CH_WZ"))	 hists[55]->Fill(mass_lvj_type0_LEP_Down,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	      if(s->name().EqualTo("DCH_WW"))	 hists[68]->Fill(mass_lvj_type0_LEP_Down,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		  { 
		    //TString name = HiggsSampleName[i]+MassPoint[j];
		    TString OrigName = s->name()+"_CMS_scale_lDown";
		    TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		    if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(mass_lvj_type0_LEP_Down,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    HistCount++;
		  }
          }
      }*/
/*
	if (1)	//-------------------	JES up
	{
	   if ( (l_pt2<0 && l_pt1>50) && (((type==0)&&(abs(l_eta1)<2.4))||((type==1)&&((abs(l_eta1)<2.5)&&!(abs(l_eta1)>1.4442 && abs(l_eta1)<1.566)))) &&
	        (((type==0)&&(pfMET_jes_up>50)) || ((type==1)&&(pfMET_jes_up>80))) &&
		((ungroomed_PuppiAK8_jet_pt_jes_up>200)&&(abs(ungroomed_PuppiAK8_jet_eta_jes_up)<2.4)&&(PuppiAK8_jet_tau2tau1<0.55)) &&
		((PuppiAK8_jet_mass_so_corr>65) && (PuppiAK8_jet_mass_so_corr<105)) &&
		(nBTagJet_loose==0) &&
		(vbf_maxpt_jj_m_jes_up>800) &&
		(abs(vbf_maxpt_j2_eta_jes_up-vbf_maxpt_j1_eta_jes_up)>4.0) &&
		((vbf_maxpt_j1_pt_jes_up>30) && (vbf_maxpt_j2_pt_jes_up>30)) && 
		(mass_lvj_type0_PuppiAK8_jes_up>600) &&
		(BosonCentrality_type0_jes_up>1.0) &&
		((abs(ZeppenfeldWL_type0_jes_up)/abs(vbf_maxpt_j2_eta_jes_up-vbf_maxpt_j1_eta_jes_up))<0.3) &&
		((abs(ZeppenfeldWH_jes_up)/abs(vbf_maxpt_j2_eta_jes_up-vbf_maxpt_j1_eta_jes_up))<0.3)
	   )
	   {
	   //if(s->name().EqualTo("data"))	 histo_data_LEPDown->Fill(mass_lvj_type0_PuppiAK8_jes_up);
	   if(s->name().EqualTo("WV_EWK"))	 hists[4] ->Fill(mass_lvj_type0_PuppiAK8_jes_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("Diboson"))  	 hists[17]->Fill(mass_lvj_type0_PuppiAK8_jes_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("top"))  	 hists[30]->Fill(mass_lvj_type0_PuppiAK8_jes_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[43]->Fill(mass_lvj_type0_PuppiAK8_jes_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("CH_WZ"))	 hists[56]->Fill(mass_lvj_type0_PuppiAK8_jes_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("DCH_WW"))	 hists[69]->Fill(mass_lvj_type0_PuppiAK8_jes_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		  { 
		    //TString name = HiggsSampleName[i]+MassPoint[j];
		    TString OrigName = s->name()+"_CMS_scale_jUp";
		    TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		    if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(mass_lvj_type0_PuppiAK8_jes_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    HistCount++;
		  }
	   }
	}*/
/*
	if (1)	//-------------------	JES down
	{
	   if ( (l_pt2<0 && l_pt1>50) && (((type==0)&&(abs(l_eta1)<2.4))||((type==1)&&((abs(l_eta1)<2.5)&&!(abs(l_eta1)>1.4442 && abs(l_eta1)<1.566)))) &&
	        (((type==0)&&(pfMET_jes_dn>50)) || ((type==1)&&(pfMET_jes_dn>80))) &&
		((ungroomed_PuppiAK8_jet_pt_jes_dn>200)&&(abs(ungroomed_PuppiAK8_jet_eta_jes_dn)<2.4)&&(PuppiAK8_jet_tau2tau1<0.55)) &&
		((PuppiAK8_jet_mass_so_corr>65) && (PuppiAK8_jet_mass_so_corr<105)) &&
		(nBTagJet_loose==0) &&
		(vbf_maxpt_jj_m_jes_dn>800) &&
		(abs(vbf_maxpt_j2_eta_jes_dn-vbf_maxpt_j1_eta_jes_dn)>4.0) &&
		((vbf_maxpt_j1_pt_jes_dn>30) && (vbf_maxpt_j2_pt_jes_dn>30)) && 
		(mass_lvj_type0_PuppiAK8_jes_dn>600) &&
		(BosonCentrality_type0_jes_dn>1.0) &&
		((abs(ZeppenfeldWL_type0_jes_dn)/abs(vbf_maxpt_j2_eta_jes_dn-vbf_maxpt_j1_eta_jes_dn))<0.3) &&
		((abs(ZeppenfeldWH_jes_dn)/abs(vbf_maxpt_j2_eta_jes_dn-vbf_maxpt_j1_eta_jes_dn))<0.3)
	   )
	   {
	   //if(s->name().EqualTo("data"))	 histo_data_LEPDown->Fill(mass_lvj_type0_PuppiAK8_jes_dn);
	   if(s->name().EqualTo("WV_EWK"))	 hists[5]->Fill(mass_lvj_type0_PuppiAK8_jes_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("Diboson"))  	 hists[18]->Fill(mass_lvj_type0_PuppiAK8_jes_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("top"))  	 hists[31]->Fill(mass_lvj_type0_PuppiAK8_jes_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[44]->Fill(mass_lvj_type0_PuppiAK8_jes_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("CH_WZ"))	 hists[57]->Fill(mass_lvj_type0_PuppiAK8_jes_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("DCH_WW"))	 hists[70]->Fill(mass_lvj_type0_PuppiAK8_jes_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		  { 
		    //TString name = HiggsSampleName[i]+MassPoint[j];
		    TString OrigName = s->name()+"_CMS_scale_jDown";
		    TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		    if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(mass_lvj_type0_PuppiAK8_jes_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    HistCount++;
		  }
	   }
	}*/
/*
	if (1)	//-------------------	JER up
	{
          if(
	   (l_pt2<0 && l_pt1>50) &&
	   (((type==0)&&(abs(l_eta1)<2.4))||((type==1)&&((abs(l_eta1)<2.5)&&!(abs(l_eta1)>1.4442 && abs(l_eta1)<1.566)))) &&
	   (((type==0)&&(pfMET_Corr_jerup>50)) || ((type==1)&&(pfMET_Corr_jerup>80))) &&
	   ((ungroomed_PuppiAK8_jet_pt>200)&&(abs(ungroomed_PuppiAK8_jet_eta)<2.4)&&(PuppiAK8_jet_tau2tau1<0.55)) &&
	   ((PuppiAK8_jet_mass_so_corr>65) && (PuppiAK8_jet_mass_so_corr<105)) &&
	   (nBTagJet_loose==0) &&
	   (vbf_maxpt_jj_m>800) &&
	   (abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta)>4.0) &&
	   ((vbf_maxpt_j1_pt>30) && (vbf_maxpt_j2_pt>30)) &&
	   (mass_lvj_type0_PuppiAK8_jer_up>600) &&
	   (BosonCentrality_type0>1.0) &&
	   ((abs(ZeppenfeldWL_type0_jer_up)/abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta))<0.3) &&
	   ((abs(ZeppenfeldWH)/abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta))<0.3)
          ){

	   //if(s->name().EqualTo("data"))	 histo_data_LEPDown->Fill(mass_lvj_type0_PuppiAK8_jer_dn);
	   if(s->name().EqualTo("WV_EWK"))	 hists[6]->Fill(mass_lvj_type0_PuppiAK8_jer_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("Diboson"))  	 hists[19]->Fill(mass_lvj_type0_PuppiAK8_jer_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("top"))  	 hists[32]->Fill(mass_lvj_type0_PuppiAK8_jer_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[45]->Fill(mass_lvj_type0_PuppiAK8_jer_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("CH_WZ"))	 hists[58]->Fill(mass_lvj_type0_PuppiAK8_jer_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("DCH_WW"))	 hists[71]->Fill(mass_lvj_type0_PuppiAK8_jer_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		  { 
		    //TString name = HiggsSampleName[i]+MassPoint[j];
		    TString OrigName = s->name()+"_CMS_res_metUp";
		    TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		    if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(mass_lvj_type0_PuppiAK8_jer_up,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    HistCount++;
		  }
          }
	}*/
/*
	if (1)	//-------------------	JER down
	{
	  if(
	   (l_pt2<0 && l_pt1>50) &&
	   (((type==0)&&(abs(l_eta1)<2.4))||((type==1)&&((abs(l_eta1)<2.5)&&!(abs(l_eta1)>1.4442 && abs(l_eta1)<1.566)))) &&
	   (((type==0)&&(pfMET_Corr_jerdn>50)) || ((type==1)&&(pfMET_Corr_jerdn>80))) &&
	   ((bos_PuppiAK8_pt>200)&&(abs(ungroomed_PuppiAK8_jet_eta)<2.4)&&(PuppiAK8_jet_tau2tau1<0.55)) &&
	   ((PuppiAK8_jet_mass_so_corr>65) && (PuppiAK8_jet_mass_so_corr<105)) &&
	   (nBTagJet_loose==0) &&
	   (vbf_maxpt_jj_m>800) &&
	   (abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta)>4.0) &&
	   ((vbf_maxpt_j1_pt>30) && (vbf_maxpt_j2_pt>30)) &&
	   (mass_lvj_type0_PuppiAK8_jer_dn>600) &&
	   (BosonCentrality_type0>1.0) &&
	   ((abs(ZeppenfeldWL_type0_jer_dn)/abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta))<0.3) &&
	   ((abs(ZeppenfeldWH)/abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta))<0.3)
	  ){

	   if(s->name().EqualTo("WV_EWK"))	 hists[7]->Fill(mass_lvj_type0_PuppiAK8_jer_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("Diboson"))  	 hists[20]->Fill(mass_lvj_type0_PuppiAK8_jer_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("top"))  	 hists[33]->Fill(mass_lvj_type0_PuppiAK8_jer_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("W+jets")||s->name().EqualTo("Z+jets"))	 hists[46]->Fill(mass_lvj_type0_PuppiAK8_jer_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("CH_WZ"))	 hists[59]->Fill(mass_lvj_type0_PuppiAK8_jer_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
	   if(s->name().EqualTo("DCH_WW"))	 hists[72]->Fill(mass_lvj_type0_PuppiAK8_jer_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
  	      HistCount = 0;
  	      for (int i=0; i<3; i++)
	        for (int j=0; j<11; j++)
		  for (int k=0; k<17; k++)
		  { 
		    //TString name = HiggsSampleName[i]+MassPoint[j];
		    TString OrigName = s->name()+"_CMS_res_metDown";
		    TString name = HiggsSampleName[i]+MassPoint[j]+Syst[k];
		    if (OrigName.EqualTo(name)) ChargedHist[HistCount]->Fill(mass_lvj_type0_PuppiAK8_jer_dn,(xsec*otherscale*genWeight*trig_eff_Weight*lep1_idEffWeight*puWeight*btagWeight)/(1.0*(nTotal-2*nNeg)));
		    HistCount++;
		  }
	  }
	}*/
      }
   // cout<<"\n\n ====>  " << hists[1]->Integral() << endl;
  }
  
  
  // include overflow bin in last bin for all histograms
  for (int i=0; i<79; i++)
  {
    hists[i]->SetBinContent(NBINS,hists[i]->GetBinContent(NBINS)+hists[i]->GetBinContent(NBINS+1));
    cout << HistName[i] << " = " << hists[i]->Integral() << endl;
    //hists[i]->Write();
  }
/*  for (int i=0; i<561; i++)
  {
    ChargedHist[i]->SetBinContent(NBINS,ChargedHist[i]->GetBinContent(NBINS)+ChargedHist[i]->GetBinContent(NBINS+1));
    cout << ChargedHist[i]->GetName() << " = " << ChargedHist[i]->Integral() << endl;
    //hists[i]->Write();
  }

  for (int i=0; i<198; i++)
  {
    ChargedHistQCD[i]->SetBinContent(NBINS,ChargedHistQCD[i]->GetBinContent(NBINS)+ChargedHistQCD[i]->GetBinContent(NBINS+1));
  }

  for (int i=0; i<3300; i++)
  {
    ChargedHistPDF[i]->SetBinContent(NBINS,ChargedHistPDF[i]->GetBinContent(NBINS)+ChargedHistPDF[i]->GetBinContent(NBINS+1));
  }

  for (int i=0; i<6; i++)
    {
      histo_diboson_EWK_CMS_QCDScaleBounding[i]->SetBinContent(NBINS,histo_diboson_EWK_CMS_QCDScaleBounding[i]->GetBinContent(NBINS)+histo_diboson_EWK_CMS_QCDScaleBounding[i]->GetBinContent(NBINS+1));
      histo_VVjjQCD_EWK_CMS_QCDScaleBounding[i]->SetBinContent(NBINS,histo_VVjjQCD_EWK_CMS_QCDScaleBounding[i]->GetBinContent(NBINS)+histo_VVjjQCD_EWK_CMS_QCDScaleBounding[i]->GetBinContent(NBINS+1));
    }
  for (int i=0; i<99; i++)
    {
      histo_diboson_EWK_CMS_PDFScaleBounding[i]->SetBinContent(NBINS,histo_diboson_EWK_CMS_PDFScaleBounding[i]->GetBinContent(NBINS)+histo_diboson_EWK_CMS_PDFScaleBounding[i]->GetBinContent(NBINS+1));
      histo_VVjjQCD_EWK_CMS_PDFScaleBounding[i]->SetBinContent(NBINS,histo_VVjjQCD_EWK_CMS_PDFScaleBounding[i]->GetBinContent(NBINS)+histo_VVjjQCD_EWK_CMS_PDFScaleBounding[i]->GetBinContent(NBINS+1));
    }
*/
  for (int j=0;j<678;j++){
    histo_aqgc[j]->SetBinContent(NBINS,histo_aqgc[j]->GetBinContent(NBINS)+histo_aqgc[j]->GetBinContent(NBINS+1));
    //std::cout << "aqgc integral " << histo_aqgc[j]->Integral() << std::endl;
  }

  //ok now we calculate the uncertainty
  std::cout << "EWK Scale uncertainties" << std::endl;
/*  for(int bin=1; bin<NBINS+1; bin++)
    {
      double systQCDScale=0;
      for (int i = 0; i<6; i++)
	{
	  if(TMath::Abs(histo_diboson_EWK_CMS_QCDScaleBounding[i]->GetBinContent(bin)-hists[1]->GetBinContent(bin)) > systQCDScale) systQCDScale = TMath::Abs(histo_diboson_EWK_CMS_QCDScaleBounding[i]->GetBinContent(bin)-hists[1]->GetBinContent(bin));
	}
      std::cout << "bin number " << bin << " " << 1 + systQCDScale/hists[1]->GetBinContent(bin) << std::endl; 
      histo_diboson_EWK_CMS_QCDScaleBounding_Up  ->SetBinContent(bin,hists[1]->GetBinContent(bin) + systQCDScale);
      histo_diboson_EWK_CMS_QCDScaleBounding_Down->SetBinContent(bin,hists[1]->GetBinContent(bin) - systQCDScale);
      systQCDScale=0;
      for (int i = 0; i<6; i++)
	{
	  if(TMath::Abs(histo_VVjjQCD_EWK_CMS_QCDScaleBounding[i]->GetBinContent(bin)-hists[14]->GetBinContent(bin)) > systQCDScale) systQCDScale = TMath::Abs(histo_VVjjQCD_EWK_CMS_QCDScaleBounding[i]->GetBinContent(bin)-hists[14]->GetBinContent(bin));
	}
      histo_VVjjQCD_EWK_CMS_QCDScaleBounding_Up  ->SetBinContent(bin,hists[14]->GetBinContent(bin) + systQCDScale);
      histo_VVjjQCD_EWK_CMS_QCDScaleBounding_Down->SetBinContent(bin,hists[14]->GetBinContent(bin) - systQCDScale);
      std::cout << "bin number " << bin << " " << 1 + systQCDScale/hists[14]->GetBinContent(bin) << std::endl; 
      HistCount = 0; 
      int CountCHhist=0;
      for (int i=0; i<3; i++)
        {
	  for (int j=0; j<11; j++)
	    {
	      systQCDScale=0;
	      for (int k = 0; k<6; k++)
		{
		  //std::cout << "Let's gop " << ChargedHistQCD[HistCount]->GetName() << " " << ChargedHist[CountCHhist*17]->GetName() << std::endl;
		  ////std::cout << "Let's gop " << ChargedHistQCD[HistCount]->GetBinContent(bin) << " " << ChargedHist[CountCHhist*17]->GetBinContent(bin) << std::endl;
		  if(TMath::Abs(ChargedHistQCD[HistCount]->GetBinContent(bin)-ChargedHist[CountCHhist*17]->GetBinContent(bin)) > systQCDScale) systQCDScale = TMath::Abs(ChargedHistQCD[HistCount]->GetBinContent(bin)-ChargedHist[CountCHhist*17]->GetBinContent(bin));
		  HistCount++;
		}
	      ChargedHist[13+CountCHhist*17]->SetBinContent(bin, ChargedHist[CountCHhist*17]->GetBinContent(bin) + systQCDScale);
	      ChargedHist[14+CountCHhist*17]->SetBinContent(bin, ChargedHist[CountCHhist*17]->GetBinContent(bin) - systQCDScale);
	      CountCHhist++;
	    }
	}
    }
*/  
  std::cout << "EWK PDF uncertainties" << std::endl;
/*  for(int bin=1; bin<NBINS+1; bin++)
    {
      double systPDFScale_1=0, systPDFScale_2=0;
      for (int i = 0; i<99; i++)
	{
	  systPDFScale_1 = systPDFScale_1 + (histo_diboson_EWK_CMS_PDFScaleBounding[i]->GetBinContent(bin)-hists[1]->GetBinContent(bin))*(histo_diboson_EWK_CMS_PDFScaleBounding[i]->GetBinContent(bin)-hists[1]->GetBinContent(bin));
	  systPDFScale_2 = systPDFScale_2 + (histo_VVjjQCD_EWK_CMS_PDFScaleBounding[i]->GetBinContent(bin)-hists[14]->GetBinContent(bin))*(histo_VVjjQCD_EWK_CMS_PDFScaleBounding[i]->GetBinContent(bin)-hists[14]->GetBinContent(bin));
	}
      systPDFScale_1 = sqrt(systPDFScale_1/99.);
      systPDFScale_2 = sqrt(systPDFScale_2/99.);
      std::cout << "bin number " << bin << " " << 1 + systPDFScale_1/hists[1]->GetBinContent(bin)  << std::endl; 
      std::cout << "bin number " << bin << " " << 1 + systPDFScale_2/hists[14]->GetBinContent(bin)  << std::endl; 
      histo_diboson_EWK_CMS_PDFScaleBounding_Up->SetBinContent(bin, hists[1]->GetBinContent(bin) +  systPDFScale_1);
      histo_diboson_EWK_CMS_PDFScaleBounding_Down->SetBinContent(bin, hists[1]->GetBinContent(bin) - systPDFScale_1);
      histo_VVjjQCD_EWK_CMS_PDFScaleBounding_Up->SetBinContent(bin, hists[14]->GetBinContent(bin) + systPDFScale_2);
      histo_VVjjQCD_EWK_CMS_PDFScaleBounding_Down->SetBinContent(bin, hists[14]->GetBinContent(bin) - systPDFScale_2);
      HistCount = 0; 
      int CountCHhist=0;
      for (int i=0; i<3; i++)
	{
	  for (int k=0; k<11; k++)
	    {
	      double systPDFScale=0;
	      for (int j = 0; j<100; j++)
		{
		  systPDFScale = systPDFScale + (ChargedHistPDF[HistCount]->GetBinContent(bin)-ChargedHist[CountCHhist*17]->GetBinContent(bin))*(ChargedHistPDF[HistCount]->GetBinContent(bin)-ChargedHist[CountCHhist*17]->GetBinContent(bin));
		  HistCount++;
		}
	      systPDFScale =  sqrt(systPDFScale/99.);
	      ChargedHist[15+CountCHhist*17]->SetBinContent(bin, ChargedHist[CountCHhist*17]->GetBinContent(bin) + systPDFScale);
	      ChargedHist[16+CountCHhist*17]->SetBinContent(bin, ChargedHist[CountCHhist*17]->GetBinContent(bin) - systPDFScale);
	      CountCHhist++;
	    }
	}
    }
 */ 

  
  //map<TString, TH1 *>::iterator mit = m_histos.find("ZV(EWK)");
  TH1 *h = 0;//mit->second;
  std::cout << "Did we ever get here? " << std::endl;
  TFile *outFile = new TFile("ch1_splitted_TF1_hfs0.root","RECREATE"); 
  TFile *outFile1 = new TFile("ch1_splitted_TF1_hfs1.root","RECREATE"); 
  TFile *outFile2 = new TFile("ch1_splitted_TF1_hfm0.root","RECREATE"); 
/*  TFile *outFile3 = new TFile("ch1_splitted_TF1_hfm1.root","RECREATE"); 
  TFile *outFile4 = new TFile("ch1_splitted_TF1_hfm2.root","RECREATE"); 
  TFile *outFile5 = new TFile("ch1_splitted_TF1_hfm3.root","RECREATE"); 
  TFile *outFile6 = new TFile("ch1_splitted_TF1_hfm4.root","RECREATE"); 
  TFile *outFile7 = new TFile("ch1_splitted_TF1_hfm5.root","RECREATE"); 
  TFile *outFile8 = new TFile("ch1_splitted_TF1_hfm6.root","RECREATE");       
  TFile *outFile9 = new TFile("ch1_splitted_TF1_hfm7.root","RECREATE"); 
  TFile *outFile10 = new TFile("ch1_splitted_TF1_hft0.root","RECREATE"); 
  TFile *outFile11 = new TFile("ch1_splitted_TF1_hft1.root","RECREATE"); 
  TFile *outFile12 = new TFile("ch1_splitted_TF1_hft2.root","RECREATE"); 
  TFile *outFile13 = new TFile("ch1_splitted_TF1_hft5.root","RECREATE"); 
  TFile *outFile14 = new TFile("ch1_splitted_TF1_hft6.root","RECREATE"); 
  TFile *outFile15 = new TFile("ch1_splitted_TF1_hft7.root","RECREATE"); 
  TFile *outFile16 = new TFile("ch1_splitted_TF1_hft8.root","RECREATE"); 
  TFile *outFile17 = new TFile("ch1_splitted_TF1_hft9.root","RECREATE");       
  // In this for loop the used histogram `hists[1]` contains information for WV_EWK
*/
  for(int i = 1; i<hists[1]->GetNbinsX()+1; i++)
    {
 gROOT->Reset();  
      stringstream ss;
      ss << i;
      std::string hist_name_temp = "bin_content_par1_"+ss.str();
      const char* hist_name = hist_name_temp.c_str();
      TH1D  *hfs0  = new TH1D(hist_name, hist_name, 42,fs0[0]-0.5,fs0[42]-0.5);
      TH1D  *hfs1  = new TH1D(hist_name, hist_name, 32,fs1[0]-0.25,fs1[32]-0.25);
      TH1D  *hfm0  = new TH1D(hist_name, hist_name, 40,fm0[0]-0.05,fm0[40]-0.05);
/*      TH1D  *hfm1  = new TH1D(hist_name, hist_name, 37,fm1[0]-2.5,fm1[36]+2.5);
      TH1D  *hfm2  = new TH1D(hist_name, hist_name, 31,fm2[0]-2.5,fm2[30]-2.5);
cout<<"+++++++++++++++++++++++++++"<<endl;
      TH1D  *hfm3  = new TH1D(hist_name, hist_name, 30,fm3[0]-2.5,fm3[30]-2.5);
      TH1D  *hfm4  = new TH1D(hist_name, hist_name, 36,fm4[0]-2.5,fm4[36]-2.5);
      TH1D  *hfm5  = new TH1D(hist_name, hist_name, 44,fm5[0]-2.5,fm5[44]-2.5);
      TH1D  *hfm6  = new TH1D(hist_name, hist_name, 34,fm6[0]-1.0,fm6[34]-1.0);
      TH1D  *hfm7  = new TH1D(hist_name, hist_name, 32,fm7[0]-2.5,fm7[32]-2.5);
      TH1D  *hft0  = new TH1D(hist_name, hist_name, 34,ft0[0]-0.1,ft0[34]-0.1);
      TH1D  *hft1  = new TH1D(hist_name, hist_name, 34,ft1[0]-0.25,ft1[34]-0.25);
      TH1D  *hft2  = new TH1D(hist_name, hist_name, 34,ft2[0]-0.25,ft2[34]-0.25);
      TH1D  *hft5  = new TH1D(hist_name, hist_name, 38,ft5[0]-0.25,ft5[38]-0.25);
      TH1D  *hft6  = new TH1D(hist_name, hist_name, 42,ft6[0]-0.25,ft6[42]-0.25);
      TH1D  *hft7  = new TH1D(hist_name, hist_name, 50,ft7[0]-0.25,ft7[50]-0.25);
      TH1D  *hft8  = new TH1D(hist_name, hist_name, 34,ft8[0]-0.25,ft8[34]-0.25);
      TH1D  *hft9  = new TH1D(hist_name, hist_name, 40,ft9[0]-0.25,ft9[40]-0.25);
*/      
      for(int j = 0; j<678; j++)
	{
	  if(j>= 526 && j<569)
	    {
//		cout<<"fs0["<<j-526<<"]: "<<fs0[j-526]<<endl;
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
		cout<<"histo_aqgc[j]: "<<histo_aqgc[j]->GetBinContent(i)<<"   hists[1]:"<<hists[1]->GetBinContent(i)<<endl;
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfs0->SetBinContent(j-527+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
//		cout<<"hfs0: "<<GetBinContent(1);
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));
		cout<<"w: "<<w<<" e1: "<<e1<<" e2:"<<e2<<" err:"<<err<<endl;
	      //std::cout << "fs0 \t " << histo_aqgc[j-446]->GetBinContent(i) << " " << hists[1]->GetBinContent(i) << " " << histo_aqgc[j-446]->GetBinError(i) << " " <<  hists[1]->GetBinError(i) << " " << err << "\t" << hfs0->GetBinContent(j-446+1) << std::endl;
	      hfs0->SetBinError(j-527+1,err);
	    }
	  else if(j> 569 && j<602)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfs1->SetBinContent(j-570+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      //std::cout << "fs1 \t " << histo_aqgc[j-446]->GetBinContent(i) << " " << hists[1]->GetBinContent(i) << " " << histo_aqgc[j-446]->GetBinError(i) << " " <<  hists[1]->GetBinError(i) << " " << err << "\t" << hfs1->GetBinContent(j-537+1) << std::endl;
	      hfs1->SetBinError(j-570+1,err);
	    }
	  else if(j>448 && j<489)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfm0->SetBinContent(j-449+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hfm0->SetBinError(j-449+1,err);
	      //hfm0->SetBinContent(j+1-604,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
/*	  else if(j>= 489 && j<526)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfm1->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hfm1->SetBinError(j+1,err);
	      //hfm1->SetBinContent(j+1-689,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>= 386 && j<417)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfm2->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hfm2->SetBinError(j+1,err);
	      //hfm1->SetBinContent(j+1-689,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>= 417 && j<448)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfm3->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hfm3->SetBinError(j+1,err);
	      //hfm1->SetBinContent(j+1-689,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>= 304 && j<341)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfm4->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hfm4->SetBinError(j+1,err);
	      //hfm1->SetBinContent(j+1-689,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>= 341 && j<386)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfm5->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hfm5->SetBinError(j+1,err);
	      //hfm1->SetBinContent(j+1-689,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>= 238 && j<273)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfm6->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hfm6->SetBinError(j+1,err);
	      //hfm6->SetBinContent(j+1-756,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));	
	    }
	  else if(j>=273 && j<304)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hfm7->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hfm7->SetBinError(j+1,err);
	      //hfm7->SetBinContent(j+1-840,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>=35 && j<70)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hft0->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hft0->SetBinError(j+1,err);
	      //hft0->SetBinContent(j+1-961,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>=0 && j<35)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hft1->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hft1->SetBinError(j+1,err);
	      //hft1->SetBinContent(j+1-1030,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>=70 && j<105)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hft2->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hft2->SetBinError(j+1,err);
	      //hft2->SetBinContent(j+1-1081,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>=105 && j<144)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hft5->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hft5->SetBinError(j+1,err);
	      //hft2->SetBinContent(j+1-1081,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>=195 && j<238)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hft6->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hft6->SetBinError(j+1,err);
	      //hft2->SetBinContent(j+1-1081,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>=144 && j<195)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hft7->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hft7->SetBinError(j+1,err);
	      //hft2->SetBinContent(j+1-1081,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>=602 && j<637)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hft8->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hft8->SetBinError(j+1,err);
	      //hft2->SetBinContent(j+1-1081,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }
	  else if(j>=637 && j<678)
	    {
	      double w = histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i);
	      double e1 = histo_aqgc[j]->GetBinError(i)/histo_aqgc[j]->GetBinContent(i);
	      double e2 = hists[1]->GetBinError(i)/hists[1]->GetBinContent(i);
	      hft9->SetBinContent(j+1,histo_aqgc[j]->GetBinContent(i)/hists[1]->GetBinContent(i));
	      double err = sqrt(e1*e1+e2*e2)*w;//TMath::Abs(((1-2*w)*e1*e1 + w*w*e2*e2 )/(hists[1]->GetBinContent(i)*hists[1]->GetBinContent(i)));  
	      hft9->SetBinError(j+1,err);
	      //hft2->SetBinContent(j+1-1081,histo_aqgc[j-446]->GetBinContent(i)/hists[1]->GetBinContent(i));
	    }*/
	}
      outFile->cd();
      TF1 *fit_1 = new TF1(hist_name,"pol2",fs0[0]-0.5,fs0[42]+0.5);
      hfs0->Fit(hist_name,"R");
   //   hfs0->Fit(fit_1,"R");
     // fit_1->Write();
      hfs0->Write();
// hfs0->Reset();
      outFile1->cd();
      TF1 *fit_2 = new TF1(hist_name,"pol2",fs1[0]-0.25,fs1[32]-0.25);
      //hfs1->Fit(hist_name,"R");
      hfs1->Write(); //hfs1->Reset();
      //fit_2->Write();
      outFile2->cd();
      //hfs0->Write();
      TF1 *fit_3 = new TF1(hist_name,"pol2",fm0[0]-0.05,fm0[40]-0.05);
      //hfm0->Fit(hist_name,"R");
      //fit_3->Write();
      hfm0->Write();
/*      outFile3->cd();
      TF1 *fit_4 = new TF1(hist_name,"pol2",fm1[0]-2.5,fm1[36]-2.5);
      //hfm1->Fit(hist_name,"R");
      //fit_4->Write();
      hfm1->Write();
      outFile4->cd();
      TF1 *fit_5 = new TF1(hist_name,"pol2",fm2[0]-2.5,fm2[30]-2.5);
      //hfm6->Fit(hist_name,"R");
      //fit_5->Write();
      hfm2->Write();
      outFile5->cd();
      TF1 *fit_6 = new TF1(hist_name,"pol2",fm3[0]-2.5,fm3[30]-2.5);
      hfm3->Write();
      outFile6->cd();
      TF1 *fit_7 = new TF1(hist_name,"pol2",fm4[0]-2.5,fm4[36]-2.5);
      hfm4->Write();
      outFile7->cd();
      TF1 *fit_8 = new TF1(hist_name,"pol2",fm5[0]-2.5,fm5[44]-2.5);
      hfm5->Write();
      outFile8->cd();
      TF1 *fit_9 = new TF1(hist_name,"pol2",fm6[0]-1.0,fm6[34]-1.1);
      hfm6->Write();
      outFile9->cd();
      TF1 *fit_10 = new TF1(hist_name,"pol2",fm7[0]-2.5,fm7[32]-2.5);
      hfm7->Write();
      outFile10->cd();
      TF1 *fit_11 = new TF1(hist_name,"pol2",ft0[0]-1.0,ft0[34]-1.0);
      hft0->Write();
      outFile11->cd();
      TF1 *fit_12 = new TF1(hist_name,"pol2",ft1[0]-0.25,ft1[34]-0.25);
      hft1->Write();
      outFile12->cd();
      TF1 *fit_13 = new TF1(hist_name,"pol2",ft2[0]-0.25,ft2[34]-0.25);
      hft2->Write();
      outFile13->cd();
      TF1 *fit_14 = new TF1(hist_name,"pol2",ft5[0]-0.25,ft5[38]-0.25);
      hft5->Write();
      outFile14->cd();
      TF1 *fit_15 = new TF1(hist_name,"pol2",ft6[0]-0.25,ft6[42]-0.25);
      hft6->Write();
      outFile15->cd();
      TF1 *fit_16 = new TF1(hist_name,"pol2",ft7[0]-0.25,ft7[50]-0.25);
      hft7->Write();
      outFile16->cd();
      TF1 *fit_17 = new TF1(hist_name,"pol2",ft8[0]-0.25,ft8[34]-0.25);
      hft8->Write();
      outFile17->cd();
      TF1 *fit_18 = new TF1(hist_name,"pol2",ft9[0]-0.25,ft9[40]-0.25);
      hft9->Write();
      //f.Write();
      delete hft0; delete  hft1; delete hft2; delete hfs0; delete hfs1; delete hfm0;   delete hfm1;   delete hfm6;   delete hfm7; 
*/
	delete hfs0; delete hfs1;
    }
  //outFile->Write();
  outFile->Close();
  outFile1->Close();
//outFile2->Close();outFile3->Close();outFile4->Close();outFile5->Close();outFile6->Close();outFile7->Close();outFile8->Close();
//  outFile9->Close();outFile10->Close();outFile11->Close();outFile12->Close();outFile13->Close();outFile14->Close();outFile15->Close();outFile16->Close();outFile17->Close();



  //TFile f("ch1_splitted_TF1_NoBinbyBin.root", "RECREATE");	// if name change then change this name also in first time where script add_stat_shapes.py appears
  TString OutRootFileSuffix = "_NoBinbyBin";
  TFile f(OutPutRootFileName + OutRootFileSuffix + ".root", "RECREATE");	// if name change then change this name also in first time where script add_stat_shapes.py appears


  // Write all histograms... 
  for (int i=0; i<79; i++)
  {
    hists[i]->Write();
  }
 for (int i=0; i<10; i++)
	{
	histo_aqgc[i]->Write();
	}
   h1->Write();
/*  for (int i=0; i<561; i++)
  {
    ChargedHist[i]->Write();
  }

  wjet->SetName("W1+jets");
  wjet->SetTitle("W1+jets");
  //wjet->SetLineColor(TColor::GetColor(222,90,106));
  //wjet->SetFillColor(TColor::GetColor(222,90,106));
  wjet->SetLineColor(TColor::GetColor(248,206,104));
  wjet->SetFillColor(TColor::GetColor(248,206,104));	
  wjet->SetLineWidth(0);
  wjet->Write();
  wjetup->SetName("shape_W+jetsUp");
  wjetup->SetTitle("shape_W+jetsUp");
  wjetup->Write();
  wjetdown->SetName("shape_W+jetsDown");
  wjetdown->SetTitle("shape_W+jetsDown");
  wjetdown->Write();
  wjetup1->SetName("shape2_W+jetsUp");
  wjetup1->SetTitle("shape2_W+jetsUp");
  wjetup1->Write();
  wjetdown1->SetName("shape2_W+jetsDown");
  wjetdown1->SetTitle("shape2_W+jetsDown");
  wjetdown1->Write();
  wjetup2->SetName("shape3_W+jetsUp");
  wjetup2->SetTitle("shape3_W+jetsUp");
  wjetup2->Write();
  wjetdown2->SetName("shape3_W+jetsDown");
  wjetdown2->SetTitle("shape3_W+jetsDown");
  wjetdown2->Write();
  wjetup3->SetName("shape4_W+jetsUp");
  wjetup3->SetTitle("shape4_W+jetsUp");
  wjetup3->Write();
  wjetdown3->SetName("shape4_W+jetsDown");
  wjetdown3->SetTitle("shape4_W+jetsDown");
  wjetdown3->Write();
  wjetup4->SetName("shape5_W+jetsUp");
  wjetup4->SetTitle("shape5_W+jetsUp");
  wjetup4->Write();
  wjetdown4->SetName("shape5_W+jetsDown");
  wjetdown4->SetTitle("shape5_W+jetsDown");
  wjetdown4->Write();


histo_VVjjQCD_EWK_CMS_QCDScaleBounding_Up->Write();
histo_VVjjQCD_EWK_CMS_QCDScaleBounding_Down->Write();
histo_VVjjQCD_EWK_CMS_PDFScaleBounding_Up->Write();
histo_VVjjQCD_EWK_CMS_PDFScaleBounding_Down->Write();
histo_diboson_EWK_CMS_QCDScaleBounding_Up->Write();
histo_diboson_EWK_CMS_QCDScaleBounding_Down->Write();
histo_diboson_EWK_CMS_PDFScaleBounding_Up->Write();
histo_diboson_EWK_CMS_PDFScaleBounding_Down->Write();
*/

    //-------------------------------------------------------------------------------------
    //
    //		Create Data card
    //
    //-------------------------------------------------------------------------------------
    #if 0
    char outputLimitsShape[200];
    sprintf(outputLimitsShape,"histo_limits_WV.txt");
    ofstream newcardShape;
    newcardShape.open(outputLimitsShape);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");

    newcardShape << Form("shapes * * WVchannel_datacard.root $PROCESS $PROCESS_$SYSTEMATIC\n");
    newcardShape << Form("shapes data_obs * WVchannel_datacard.root histo_Data\n");
    newcardShape << Form("shapes Higgs * WVchannel_datacard.root histo_Higgs_M$MASS histo_Higgs_M$MASS_$SYSTEMATIC\n");
    newcardShape << Form("Observation %d\n", -1/*(int)histo_Data->GetBinContent(nb)*/);
    //newcardShape << Form("bin wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
    newcardShape << Form("process aQGC Wjet WV top Zjet\n");
    newcardShape << Form("process 0 1 2 3 4\n");
    newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",-1.,-1.,-1.,-1.,-1.,-1.) ;
    //-------------------------------------------------------------------------------------

    newcardShape.close();
    #endif

  f.Write();
  f.Close();

  //-----------------------------------------------------
  //
  //	Add bin-by-bin uncertanities
  //
  //-----------------------------------------------------
/*  TString command1 = "./add_stat_shapes.py --filter Vjets --prefix Vjets_bbb " + OutPutRootFileName + OutRootFileSuffix + ".root WVchannel_datacard_BBB2.root";
  system(command1);

  char command2[3000];
  sprintf(command2,"./add_stat_shapes.py --filter diboson --prefix diboson_bbb WVchannel_datacard_BBB2.root WVchannel_datacard_BBB3.root");
  system(command2);

  char command3[3000];
  sprintf(command3,"./add_stat_shapes.py --filter VVjjQCD --prefix VVjjQCD_bbb WVchannel_datacard_BBB3.root WVchannel_datacard_BBB4.root");
  system(command3);

  char command4[3000];
  sprintf(command4,"./add_stat_shapes.py --filter top --prefix top_bbb WVchannel_datacard_BBB4.root WVchannel_datacard_BBB5.root");
  system(command4);

  char command5[3000];
  sprintf(command5,"./add_stat_shapes.py --filter W1+jets --prefix W1+jets_bbb WVchannel_datacard_BBB5.root WVchannel_datacard_BBB6.root");
  system(command5);

  TString command6 = "rm WVchannel_datacard_BBB2.root WVchannel_datacard_BBB3.root WVchannel_datacard_BBB4.root WVchannel_datacard_BBB5.root; mv WVchannel_datacard_BBB6.root " + OutPutRootFileName + ".root";
  system(command6);
*/}

void WVChannel_GetCard_WithHiggsDistributions()
{
  int start_s=clock(); 
  //model("DibosonBoostedElMuSamples13TeV_2019_03_29_03h24.txt",
  model("files2016_WV.txt",
	"ch1_splitted_TF1_WV");

  int stop_s=clock();
  cout << "time: " << double(stop_s-start_s)/(double(CLOCKS_PER_SEC)*60.0) <<" min" << endl;

}

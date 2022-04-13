#include "setTDRStyle.h"
#include <sstream>
#include <fstream>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include <vector>
#include "Rtypes.h"
#include "TColor.h"
#include "TVectorF.h"
#include <cstdlib>
#include "TGraphErrors.h"

struct shower_reco_perf {
  TGraphErrors* h_response;
  TGraphErrors* h_linearity;
  TGraphErrors* h_resolution;
  TGraphErrors* h_resolution_corr;
};


void setTDRStyle();
shower_reco_perf getPerformance(std::vector<TTree*> trees_, std::vector<unsigned int> energies_, TString var_, int color_,TString suffix, TString region, TString cut="0==0");
shower_reco_perf getPerformancePosition(std::vector<TTree*> trees_, std::vector<float> pos_, float cpEnergy_, TString var_, int color_, TString suffix,TString cut="0==0",bool isEta=false);
float calcUncorUnc(float num, float num_err, float den, float den_err);
TH1F *getVarRatio(TTree *tr, string en_str_, TString varNum, TString varDen, int nbins, float xmin, float xmax, TString name);
void addOvFlow(TH1F *h);
float getMeanOfBranch(TTree *tree, TString varName, int nbins, float xmin, float xmax);
std::string replaceFirstOccurrence(
    std::string& original,
    const std::string& toReplace,
    const std::string& replaceWith)
{
    std::string s = original;
    std::size_t pos = s.find(toReplace);
    if (pos == std::string::npos) {
      std::cerr << "There is no " << toReplace << " in " << s;
      return s;
    }
    return s.replace(pos, toReplace.length(), replaceWith);
};

void performance(std::string inputFileTemplate = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/apsallid/CalibrationStudies/ResolutionTrees_20220327/CE_E_Front_120um/hgc_singlephoton_eENERGYGeV_nopu.root", std::string region = "CE_E_Front_120um") {
  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(1);
  TH1::SetDefaultSumw2(kTRUE);

  std::vector<unsigned int> energies_;
  energies_.push_back(10.);
  energies_.push_back(100.);
  energies_.push_back(120.);
  energies_.push_back(15.);
  energies_.push_back(150.);
  energies_.push_back(180.);
  energies_.push_back(2.);
  energies_.push_back(200.);
  energies_.push_back(250.);
  energies_.push_back(3.);
  energies_.push_back(30.);
  energies_.push_back(300.);
  energies_.push_back(400.);
  energies_.push_back(500.);
  energies_.push_back(5.);
  energies_.push_back(50.);
  energies_.push_back(8.);
  energies_.push_back(80.);

  std::vector<TTree*> trees_nopu_;
  for (unsigned int itree = 0; itree<energies_.size(); ++itree) {
    std::string inputFile = replaceFirstOccurrence(inputFileTemplate, "ENERGY", to_string(energies_[itree]));
    std::cout << "Input file : " << inputFile + region << std::endl;
    TFile *f_ = TFile::Open(inputFile.c_str());
    TTree *t_ = (TTree*)f_->Get("finaltree");
    trees_nopu_.push_back(t_);
  }
  std::cout << " performance: recable \n";
  shower_reco_perf perf_recable        = getPerformance(trees_nopu_,energies_,"reconstructable_energy_" + region,4,"erecable", region, "reconstructable_energy_" + region + "<9999999");
  //std::cout << " performance: clusterized \n";
  // shower_reco_perf perf_uncalib     = getPerformance(trees_nopu_,energies_,"uncalib_energy_" + region,1,"euncalib",region,"uncalib_energy_" + region+"<9999999");
  //shower_reco_perf perf_had_clusterized = getPerformance(trees_nopu_,energies_,"cp_clusterized_e",6,"clusterized");
////  std::cout << " performance: reco clusterized \n";
////  shower_reco_perf perf_had_fromLCs_clusterized = getPerformance(trees_nopu_,energies_,"lc_fromCP_clusterized_e",7,"lc_clusterized");
////  std::cout << " performance: reco clusterized (only LC size more than 2) \n";
////  shower_reco_perf perf_had_fromLCs_clusterized_Size2 = getPerformance(trees_nopu_,energies_,"lc_fromCP_clusterizedLCSizeMoreThan2_e",52,"lc_clusterized_size2");
////  std::cout << " performance: reco clusterized (only LC size more than 3) \n";
////  shower_reco_perf perf_had_fromLCs_clusterized_Size3 = getPerformance(trees_nopu_,energies_,"lc_fromCP_clusterizedLCSizeMoreThan3_e",2,"lc_clusterized_size3");
//
//  // response w.r.t. cp_e (energy of the Caloparticle)
////  gStyle->SetOptTitle(0);
//  TCanvas *c_response = new TCanvas("c_response","c_response",500,500);
//  perf_had_nopu.h_response->GetYaxis()->SetTitle("E_{RECO}/E_{GEN}");
//  perf_had_nopu.h_response->GetXaxis()->SetTitle("E_{GEN} [GeV]");
//  perf_had_nopu.h_response->GetYaxis()->SetRangeUser(0.,1.3);
//  perf_had_nopu.h_response->Draw("AP");
//  perf_had_recable.h_response->Draw("P sames");
//  perf_had_clusterized.h_response->Draw("P sames");
////  perf_had_fromLCs_clusterized.h_response->SetMarkerStyle(22);
////  perf_had_fromLCs_clusterized.h_response->Draw("P sames");
///*
//  TCanvas *c_response_LargeLCs = new TCanvas("c_response_LargeLCs","c_response_LargeLCs",500,500);
//  perf_had_nopu.h_response->GetYaxis()->SetTitle("E_{RECO}/E_{GEN}");
//  perf_had_nopu.h_response->GetXaxis()->SetTitle("E_{GEN} [GeV]");
//  perf_had_nopu.h_response->GetYaxis()->SetRangeUser(0.,1.3);
//  perf_had_nopu.h_response->Draw("AP");
//  perf_had_recable.h_response->Draw("P sames");
//  perf_had_fromLCs_clusterized.h_response->Draw("P sames");
//  perf_had_fromLCs_clusterized_Size2.h_response->SetMarkerStyle(21);
//  perf_had_fromLCs_clusterized_Size2.h_response->Draw("P sames");
//  perf_had_fromLCs_clusterized_Size3.h_response->SetMarkerStyle(23);
//  perf_had_fromLCs_clusterized_Size3.h_response->Draw("P sames");
//*/
//
//  gStyle->SetOptStat(1);
//  TCanvas *c_linearity = new TCanvas("c_linearity","c_linearity",500,500);
//  perf_had_nopu.h_linearity->GetYaxis()->SetTitle("E_{RECO} [GeV]");
//  perf_had_nopu.h_linearity->GetXaxis()->SetTitle("E_{GEN} [GeV]");
//  perf_had_nopu.h_linearity->GetYaxis()->SetRangeUser(0.,300.);
//  perf_had_nopu.h_linearity->Draw("AP");
//  perf_had_recable.h_linearity->Draw("P sames");
//  perf_had_clusterized.h_linearity->Draw("P sames");
//
//
//  // resolution
//
  gStyle->SetOptStat(0);
  TCanvas *c_resolution = new TCanvas("c_resolution","c_resolution",500,500);
  perf_recable.h_resolution->Draw("AP");
  perf_recable.h_resolution->GetYaxis()->SetRangeUser(0.,0.5);
  // perf_recable.h_resolution->GetYaxis()->SetTitle("#sigma(E_{RECO}) / <E_{RECO}>");
  perf_recable.h_resolution->GetYaxis()->SetTitle("#sigma(E)/E");
  perf_recable.h_resolution->GetXaxis()->SetTitle("E_{GEN}");
  // perf_uncalib.h_resolution->Draw("P sames");
  // perf_had_clusterized.h_resolution->Draw("P sames");

  gStyle->SetOptStat(0);
  TCanvas *c_resolution_corr = new TCanvas("c_resolution_corr","c_resolution_corr",500,500);
  perf_recable.h_resolution_corr->Draw("AP");
  perf_recable.h_resolution_corr->GetYaxis()->SetRangeUser(0.,0.5);
  //perf_recable.h_resolution_corr->GetYaxis()->SetTitle("#sigma(E_{RECO}) / <E_{RECO}>");
  perf_recable.h_resolution_corr->GetYaxis()->SetTitle("#sigma(E)/E");
  perf_recable.h_resolution_corr->GetXaxis()->SetTitle("E_{GEN}");
  // perf_uncalib.h_resolution_corr->Draw("P sames");
  // perf_had_clusterized.h_resolution_corr->Draw("P sames");
  return;
}

void mainfunctionB() {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(1);
  TH1::SetDefaultSumw2(kTRUE);

  //  TFile *f_nopu = TFile::Open("/data/hgcal-0/user/gouskos/trees/default_version/hgc_tree_singlepi_e50GeV_pu200.root", "READONLY" );
  TFile *f_nopu = TFile::Open("/data/hgcal-0/user/gouskos/trees/default_version/hgc_tree_singlepi_e50GeV_pu200.root", "READONLY" );
  TTree *t_nopu = (TTree*)f_nopu->Get("pumitigation/tree");

  TH1F *h1_ = new TH1F("h1_","h1_",9,10.,100); h1_->SetLineColor(1);
  t_nopu->Project("h1_","smc_had_e"); 

  TH1F *h2_ = new TH1F("h2_","h2_",9,10.,100); h2_->SetLineColor(2);
  t_nopu->Project("h2_","smc_sk_had_e"); 

  TH1F *h3_ = new TH1F("h3_","h3_",9,10.,100); h3_->SetLineColor(4);
  t_nopu->Project("h3_","smc_time_had_e"); 
  
  TCanvas *c = new TCanvas("c","c",500,500);
  h1_->DrawNormalized("HIST");
  h2_->DrawNormalized("HIST sames");
  h3_->DrawNormalized("HIST sames");

}



void performanceAreaScan(float energy, bool isEta) {

  TH1::SetDefaultSumw2(kTRUE);

  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2(kTRUE);


  std::vector<float> r_;
  if (!isEta) {
    //r_.push_back(25);
    r_.push_back(30);
    r_.push_back(35);
    r_.push_back(40);
    r_.push_back(50);
    //r_.push_back(55);
    r_.push_back(60);
    r_.push_back(70);
    r_.push_back(80);
    r_.push_back(90);
    r_.push_back(100);
    r_.push_back(110);
    r_.push_back(120);
    r_.push_back(130);
    r_.push_back(140);
    //r_.push_back(145);
    //r_.push_back(150);
    //r_.push_back(155);
  }
  else {
      r_.push_back(1.5);
    //r_.push_back(1.57);
    r_.push_back(1.64);
    r_.push_back(1.69);
    r_.push_back(1.8);
    //    r_.push_back(1.89);
    r_.push_back(1.98);
    r_.push_back(2.1);
    r_.push_back(2.22);
    r_.push_back(2.36);
    r_.push_back(2.55);
    r_.push_back(2.8);
    r_.push_back(2.9);
    //r_.push_back(3.05);
  }


  std::vector<unsigned int> colors_;
  colors_.push_back(1);
  colors_.push_back(2);
  colors_.push_back(4);
  colors_.push_back(6);
  colors_.push_back(8);
  
  stringstream en_str_tmp_; en_str_tmp_ << energy;
  string en_str_; en_str_tmp_ >> en_str_;
 

  std::vector<TTree*> trees;
  for (unsigned int itree = 0; itree<r_.size(); ++itree) {

    stringstream r_str_tmp_; r_str_tmp_ << r_[itree];
    string r_str_; r_str_tmp_ >> r_str_;
    
    //TFile *f_ = TFile::Open("/data/hgcal-0/user/gouskos/samples/forCTDots2020_111X/photons_closeby_fixedenergy_scaneta/trees/hgc_perftree_singlephoton_e"+en_str_+"GeV_r"+(TString)r_str_+"_nopu.root", "READONLY" );
    //    TFile *f_ = TFile::Open("/data2/user/gouskos/samples/forCTDots2020_111X/photons_closeby_fixedenergy_scaneta/trees/hgc_perftree_singlephoton_e"+en_str_+"GeV_r"+(TString)r_str_+"_nopu.root", "READONLY" );
    //    TFile *f_ = TFile::Open("/data2/user/gouskos/samples/forCTDots2020_111X/photons_closeby_fixedenergy_scaneta/trees_stepSize2/hgc_perftree_singlephoton_e"+en_str_+"GeV_r"+(TString)r_str_+"_nopu.root", "READONLY" );
    //TFile *f_ = TFile::Open("/data2/user/gouskos/CMSSW_11_1_X_2020-02-18-1100/src/hgc_analyzer_for_ml/HGCPerformanceAnalyzer/test/hgc_perfstudiestree_singlephoton_e"+en_str_+"GeV_r"+(TString)r_str_+"_nopu.root", "READONLY" );
    //TFile *f_ = TFile::Open("/data2/user/gouskos/CMSSW_11_1_X_2020-02-28-2300/src/hgc_analyzer_for_ml/HGCPerformanceAnalyzer/test/hgc_perfstudiestree_singlephoton_e"+en_str_+"GeV_r"+(TString)r_str_+"_nopu.root", "READONLY" );
    //    std::cout << "/data2/user/gouskos/samples/forCTDots2020_111X/singlepi_fixedenergy_scaneta/trees/hgc_perfstudiestree_singlepi_e"+en_str_+"GeV_r"+(TString)r_str_+"_nopu.root" << "\n";

    //    TFile *f_ = TFile::Open("/data2/user/gouskos/samples/forCTDots2020_HGCALDPGMar18/singlepi_fixedenergy_scaneta/trees/hgc_perfstudiestree_singlepi_e"+en_str_+"GeV_r"+(TString)r_str_+"_nopu.root", "READONLY" );
    TFile *f_ = TFile::Open("");///data2/user/gouskos/samples/forCTDots2020_HGCALDPGMar18/photons_closeby_fixedenergy_scaneta/trees/hgc_perftree_singlephoton_e"+en_str_+"GeV_r"+(TString)r_str_+"_nopu.root", "READONLY" );
    //TFile *f_ = TFile::Open("/data2/user/gouskos/samples/forCTDots2020_HGCALDPGMar18/photons_closeby_fixedenergy_scaneta/trees/hgc_perftree_singlephoton_e"+en_str_+"GeV_r"+(TString)r_str_+"_nopu.root", "READONLY" );
    TTree *t_ = (TTree*)f_->Get("performancestudies/tree");
    trees.push_back(t_);
  }

 
  //  shower_reco_perf getPerformancePosition(std::vector<TTree*> trees_, std::vector<unsigned int> pos_, float cpEnergy_, TString var_, int color_, TString suffix)
  shower_reco_perf perf_had_nopu        = getPerformancePosition(trees,r_,energy,"ts_energy",4,"nopu","0==0",isEta);                    std::cout << " performance: nopu \n";
  shower_reco_perf perf_had_recable     = getPerformancePosition(trees,r_,energy,"cp_recable_e",1,"recable","0==0",isEta);              std::cout << " performance: recable \n";
  shower_reco_perf perf_had_clusterized = getPerformancePosition(trees,r_,energy,"cp_clusterized_e",2,"clusterized","0==0",isEta);      std::cout << " performance: clusterized \n";
  //shower_reco_perf perf_had_tickled     = getPerformancePosition(trees,r_,energy,"ts_sumenergy",6,"tickled","0==0",isEta);              std::cout << " performance: ticled \n";
  //  shower_reco_perf perf_had_nopunoint   = getPerformancePosition(trees,r_,energy,"ts_energy",8,"nopunoint","cp_nsimclusters==1",isEta); std::cout << " performance: nopunoint \n";

  // response
  TCanvas *c_response = new TCanvas("c_response","c_response",500,500);
  perf_had_nopu.h_response->GetYaxis()->SetTitle("E_{RECO}/E_{GEN}");
  //perf_had_nopu.h_response->GetXaxis()->SetTitle("R [cm]");
  perf_had_nopu.h_response->GetXaxis()->SetTitle("#eta_{GEN}");
  perf_had_nopu.h_response->GetYaxis()->SetRangeUser(0.,1.1);
  perf_had_nopu.h_response->Draw("AP");
  perf_had_recable.h_response->Draw("P sames");
  perf_had_clusterized.h_response->Draw("P sames");
  //perf_had_tickled.h_response->Draw("P sames");
  //  perf_had_nopunoint.h_response->Draw("P sames");

  // linearity
  gStyle->SetOptStat(1);
  TCanvas *c_linearity = new TCanvas("c_linearity","c_linearity",500,500);
  perf_had_nopu.h_linearity->GetYaxis()->SetTitle("E_{RECO} [GeV]");
  //perf_had_nopu.h_linearity->GetXaxis()->SetTitle("R [cm]");
  perf_had_nopu.h_linearity->GetXaxis()->SetTitle("#eta_{GEN}");
  perf_had_nopu.h_linearity->GetYaxis()->SetRangeUser(20.,160.);
  perf_had_nopu.h_linearity->Draw("AP");
  perf_had_recable.h_linearity->Draw("P sames");
  perf_had_clusterized.h_linearity->Draw("P sames");
  //perf_had_tickled.h_linearity->Draw("P sames");
  //  perf_had_nopunoint.h_linearity->Draw("P sames");

  // fit linearity
  //perf_had_nopu.h_linearity->Fit("pol1");
  //TF1 *fit_nopu = perf_had_nopu.h_linearity->GetFunction("pol1");
  //fit_nopu->SetLineColor(1); fit_nopu->SetLineStyle(2);
  //fit_nopu->Draw("same");

  // resolution
  gStyle->SetOptStat(0);
  TCanvas *c_resolution = new TCanvas("c_resolution","c_resolution",500,500);
  perf_had_nopu.h_resolution->Draw("AP");
  perf_had_nopu.h_resolution->GetYaxis()->SetRangeUser(0.,0.5);
  perf_had_nopu.h_resolution->GetYaxis()->SetTitle("#sigma(E_{RECO}/E_{GEN}) / <E_{RECO}/E_{GEN}>");
  //perf_had_nopu.h_resolution->GetXaxis()->SetTitle("R [cm]");
  perf_had_nopu.h_resolution->GetXaxis()->SetTitle("#eta_{GEN}");
  perf_had_recable.h_resolution->Draw("P sames");
  perf_had_clusterized.h_resolution->Draw("P sames");
  //perf_had_tickled.h_resolution->Draw("P sames");
  //perf_had_nopunoint.h_resolution->Draw("P sames");
  

  /*
  const int n = trees.size();
  float gr_deta_0p025_[n]; float gr_deta_0p05_[n]; float gr_deta_0p075_[n];
  float gr_dphi_0p025_[n]; float gr_dphi_0p05_[n]; float gr_dphi_0p075_[n];
  float gr_r_[n]; 
  int bins = 100; float xmin = 0; float xmax = 20.;
  for (unsigned int itree = 0; itree<trees.size(); ++itree) {
    
    gr_r_[itree] = r_[itree];

    gr_deta_0p025_[itree] = getMeanOfBranch(trees[itree],"cl_e_deta0p025/cp_e",bins,xmin,xmax);
    gr_deta_0p05_[itree]  = getMeanOfBranch(trees[itree],"cl_e_deta0p05/cp_e",bins,xmin,xmax);
    gr_deta_0p075_[itree] = getMeanOfBranch(trees[itree],"cl_e_deta0p075/cp_e",bins,xmin,xmax);

    gr_dphi_0p025_[itree] = getMeanOfBranch(trees[itree],"cl_e_dphi0p025/cp_e",bins,xmin,xmax);
    gr_dphi_0p05_[itree]  = getMeanOfBranch(trees[itree],"cl_e_dphi0p05/cp_e",bins,xmin,xmax);
    gr_dphi_0p075_[itree] = getMeanOfBranch(trees[itree],"cl_e_dphi0p075/cp_e",bins,xmin,xmax);

  }

  TGraph* gr_deta_0p025 = new TGraph(n,gr_r_,gr_deta_0p025_);
  TGraph* gr_deta_0p05  = new TGraph(n,gr_r_,gr_deta_0p05_);
  TGraph* gr_deta_0p075 = new TGraph(n,gr_r_,gr_deta_0p075_);

  TGraph* gr_dphi_0p025 = new TGraph(n,gr_r_,gr_dphi_0p025_);
  TGraph* gr_dphi_0p05  = new TGraph(n,gr_r_,gr_dphi_0p05_);
  TGraph* gr_dphi_0p075 = new TGraph(n,gr_r_,gr_dphi_0p075_);

  TCanvas* c_gr = new TCanvas("c_gr","c_gr",500,500);
  gr_deta_0p025->GetXaxis()->SetTitle("CaloParticle R [cms]");
  gr_deta_0p025->GetYaxis()->SetTitle("Energy fraction for #DeltaX<Y");
  gr_deta_0p025->GetYaxis()->SetRangeUser(0.7,1.05);
  gr_deta_0p025->SetMarkerStyle(20); gr_deta_0p025->SetMarkerColor(1);
  gr_deta_0p05->SetMarkerStyle(21);  gr_deta_0p05->SetMarkerColor(1);
  gr_deta_0p075->SetMarkerStyle(22); gr_deta_0p075->SetMarkerColor(1);
  gr_deta_0p025->Draw("AP");
  gr_deta_0p05->Draw("P sames");
  gr_deta_0p075->Draw("P sames");

  gr_dphi_0p025->SetMarkerStyle(24); gr_dphi_0p025->SetMarkerColor(4);
  gr_dphi_0p05->SetMarkerStyle(25);  gr_dphi_0p05->SetMarkerColor(4);
  gr_dphi_0p075->SetMarkerStyle(26); gr_dphi_0p075->SetMarkerColor(4);
  gr_dphi_0p025->Draw("P sames");
  gr_dphi_0p05->Draw("P sames");
  gr_dphi_0p075->Draw("P sames");
  */

}



// void performanceStd(){

//   TH1::SetDefaultSumw2(kTRUE);

//   setTDRStyle();
//   gROOT->SetBatch(false);
//   gStyle->SetOptStat(0);
//   gStyle->SetOptFit(0);
//   gStyle->SetPalette(1);
//   gStyle->SetOptTitle(0);
//   TH1::SetDefaultSumw2(kTRUE);

//   std::vector<unsigned int> energies_;
//   energies_.push_back(10.);
//   energies_.push_back(20.);
//   energies_.push_back(50.);
//   energies_.push_back(100.);
//   energies_.push_back(200.);
//   energies_.push_back(300.);

//   std::vector<unsigned int> colors_;
//   colors_.push_back(1);
//   colors_.push_back(2);
//   colors_.push_back(4);
//   colors_.push_back(6);
//   colors_.push_back(8);

//   // h/e+h ratios
//   std::vector<TH1F*> h_en_sci_rechit_ov_cp; h_en_sci_rechit_ov_cp.clear();
//   std::vector<TH1F*> h_en_sci_rechit_ov_recable; h_en_sci_rechit_ov_recable.clear();
//   std::vector<TH1F*> h_en_sci_rechitEmMC_ov_emMC; h_en_sci_rechitEmMC_ov_emMC.clear();
//   std::vector<TH1F*> h_en_sci_rechitHadMC_ov_hadMC; h_en_sci_rechitHadMC_ov_hadMC.clear();

//   // get trees
//   std::vector<TTree*> trees_nopu_;
//   std::vector<TTree*> trees_pu_;
//   for (unsigned int itree = 0; itree<energies_.size(); ++itree) {

//     stringstream en_str_tmp_; en_str_tmp_ << energies_[itree];
//     string en_str_; en_str_tmp_ >> en_str_;
 
//     //TFile *f_nopu = TFile::Open("/afs/cern.ch/work/g/gouskos/private/hgcal/CMSSW_11_0_X_2019-09-02-2300/src/hgc_tree_singlepi_e"+(TString)en_str_+"GeV.root", "READONLY" );
//     //TFile *f_nopu = TFile::Open("/data/hgcal-0/user/gouskos/CMSSW_11_0_0_pre9/src/trees/hgc_tree_singlepi_e"+(TString)en_str_+"GeV.root", "READONLY" );
//     //TFile *f_pu = TFile::Open("/data/hgcal-0/user/gouskos/CMSSW_11_0_0_pre9/src/trees/hgc_tree_singlepi_e"+(TString)en_str_+"GeV_pu200.root", "READONLY" );
//     //TFile *f_nopu = TFile::Open("/data/hgcal-0/user/gouskos/trees/default_version/hgc_tree_singlepi_e"+(TString)en_str_+"GeV.root", "READONLY" );
//     //TFile *f_nopu = TFile::Open("/data/hgcal-0/user/gouskos/trees/default_version/patatrack_v1/hgc_tree_singlepi_e"+(TString)en_str_+"GeV.root", "READONLY" );
//     //TFile *f_pu = TFile::Open("/data/hgcal-0/user/gouskos/trees/default_version/patatrack_v1/hgc_tree_singlepi_e"+(TString)en_str_+"GeV_pu200.root", "READONLY" );
//     //TFile *f_pu = TFile::Open("/data/hgcal-0/user/gouskos/trees/default_version/patatrack_v0/hgc_tree_singlepi_e"+(TString)en_str_+"GeV_pu200.root", "READONLY" );
//     //TFile *f_pu = TFile::Open("/data/hgcal-0/user/gouskos/trees/default_version/hgc_tree_singlepi_eta2p4_e"+(TString)en_str_+"GeV_pu200.root", "READONLY" );
//     //    TFile *f_pu = TFile::Open("/data/hgcal-0/user/gouskos/trees/default_version/hgc_tree_singlepi_e"+(TString)en_str_+"GeV_pu200.root", "READONLY" );

//     //TFile *f_nopu = TFile::Open("/data/hgcal-0/user/gouskos/samples/forCTDots2020_111X/trees/photons_closeby_hgcalcenter/hgc_perftree_singlephoton_e"+(TString)en_str_+"GeV_nopu.root", "READONLY" );
    
//     //TFile *f_nopu = TFile::Open("/data2/user/gouskos/samples/forCTDots2020_111X/singlepi_hgcalcenter/trees/hgc_perftree_singlepi_e"+(TString)en_str_+"GeV_nopu.root", "READONLY" );
//     //TFile *f_nopu = TFile::Open("/data2/user/gouskos/samples/forCTDots2020_HGCALDPGMar18/photons_closeby_hgcalcenter/trees/hgc_perftree_singlephoton_e"+(TString)en_str_+"GeV_nopu.root", "READONLY" );
//     TFile *f_nopu = TFile::Open("/data2/user/gouskos/samples/forCTDots2020_HGCALDPGMar18/singlepi_hgcalcenter/trees/hgc_perftree_singlepi_e"+(TString)en_str_+"GeV_nopu.root", "READONLY" );

//     TTree *t_nopu = (TTree*)f_nopu->Get("performance/tree");
//     //    TTree *t_nopu = (TTree*)f_nopu->Get("hgcanalyzer/tree");
//     //    TTree *t_pu = (TTree*)f_pu->Get("pumitigation/tree");
//     trees_nopu_.push_back(t_nopu);
//     //    trees_pu_.push_back(t_pu);

//   }

//   //  shower_reco_perf perf_had_nopu = getPerformance(trees_nopu_,energies_,"mc_em_e",1);
//   std::cout << " performance: nopu \n";
//   shower_reco_perf perf_had_nopu        = getPerformance(trees_nopu_,energies_,"ts_energy",4,"nopu");
//   std::cout << " performance: recable \n";
//   shower_reco_perf perf_had_recable     = getPerformance(trees_nopu_,energies_,"cp_recable_e",1,"recable");
//   std::cout << " performance: clusterized \n";
//   shower_reco_perf perf_had_clusterized = getPerformance(trees_nopu_,energies_,"cp_clusterized_e",2,"clusterized");
//   //  std::cout << " performance: ticled \n";  
//   //shower_reco_perf perf_had_tickled     = getPerformance(trees_nopu_,energies_,"ts_sumenergy",6,"tickled");
//   std::cout << " performance: nopunoint \n";
//   shower_reco_perf perf_had_nopunoint   = getPerformance(trees_nopu_,energies_,"ts_energy",8,"nopunoint","cp_nsimclusters==1");
//   //shower_reco_perf perf_had_punosk  = getPerformance(trees_pu_,energies_,"smc_had_match_e",8,"punosk"); std::cout << " performance: pu (no sk) \n";  
//   //shower_reco_perf perf_had_pu      = getPerformance(trees_pu_,energies_,"smc_had_match_sk_e",2,"pu");  std::cout << " performance: pu (sk) \n";
  

//   TCanvas *c_response = new TCanvas("c_response","c_response",500,500);
//   perf_had_nopu.h_response->GetYaxis()->SetTitle("E_{RECO}/E_{GEN}");
//   perf_had_nopu.h_response->GetXaxis()->SetTitle("E_{GEN} [GeV]");
//   perf_had_nopu.h_response->GetYaxis()->SetRangeUser(0.,1.3);
//   perf_had_nopu.h_response->Draw("AP");
//   //perf_had_pu.h_response->Draw("P sames");
//   perf_had_recable.h_response->Draw("P sames");
//   perf_had_clusterized.h_response->Draw("P sames");
//   //perf_had_tickled.h_response->Draw("P sames");
//   perf_had_nopunoint.h_response->Draw("P sames");
//   //perf_had_punosk.h_response->Draw("P sames");

//   gStyle->SetOptStat(1);
//   TCanvas *c_linearity = new TCanvas("c_linearity","c_linearity",500,500);
//   perf_had_nopu.h_linearity->GetYaxis()->SetTitle("E_{RECO} [GeV]");
//   perf_had_nopu.h_linearity->GetXaxis()->SetTitle("E_{GEN} [GeV]");
//   perf_had_nopu.h_linearity->GetYaxis()->SetRangeUser(0.,300.);
//   perf_had_nopu.h_linearity->Draw("AP");
//   //perf_had_pu.h_linearity->Draw("P sames");
//   perf_had_recable.h_linearity->Draw("P sames");
//   perf_had_clusterized.h_linearity->Draw("P sames");
//   //perf_had_tickled.h_linearity->Draw("P sames");
//   perf_had_nopunoint.h_linearity->Draw("P sames");
//   //perf_had_punosk.h_linearity->Draw("P sames");

//   // fit linearity
//   perf_had_recable.h_linearity->Fit("pol1");
//   TF1 *fit_recable = perf_had_recable.h_linearity->GetFunction("pol1"); 
//   fit_recable->SetLineColor(1); fit_recable->SetLineStyle(2);
//   fit_recable->Draw("same");

//   perf_had_nopu.h_linearity->Fit("pol1");
//   TF1 *fit_nopu = perf_had_nopu.h_linearity->GetFunction("pol1"); 
//   fit_nopu->SetLineColor(4); fit_nopu->SetLineStyle(2);
//   fit_nopu->Draw("same");

//   perf_had_nopunoint.h_linearity->Fit("pol1");
//   TF1 *fit_nopunoint = perf_had_nopunoint.h_linearity->GetFunction("pol1"); 
//   fit_nopunoint->SetLineColor(8); fit_nopunoint->SetLineStyle(2);
//   fit_nopunoint->Draw("same");
  
//   /*
//   perf_had_punosk.h_linearity->Fit("pol1");
//   TF1 *fit_punosk = perf_had_punosk.h_linearity->GetFunction("pol1"); 
//   fit_punosk->SetLineColor(8); fit_punosk->SetLineStyle(2);
//   fit_punosk->Draw("same");
//   */

//   gStyle->SetOptStat(0);
//   TCanvas *c_resolution = new TCanvas("c_resolution","c_resolution",500,500);
//   perf_had_nopu.h_resolution->Draw("AP");
//   perf_had_nopu.h_resolution->GetYaxis()->SetRangeUser(0.,0.3);
//   perf_had_nopu.h_resolution->GetYaxis()->SetTitle("#sigma(E_{RECO}/E_{GEN}) / <E_{RECO}/E_{GEN}>");
//   perf_had_nopu.h_resolution->GetXaxis()->SetTitle("E_{GEN} [GeV]");
//   //  perf_had_pu.h_resolution->Draw("P sames");
//   perf_had_recable.h_resolution->Draw("P sames");
//   perf_had_clusterized.h_resolution->Draw("P sames");
//   //perf_had_tickled.h_resolution->Draw("P sames");
//   perf_had_nopunoint.h_resolution->Draw("P sames");
//   //perf_had_punosk.h_resolution->Draw("P sames");


//   TF1 *f_res_recable = new TF1("f_res_recable","sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)))");
//   f_res_recable->SetLineColor(1);
//   perf_had_recable.h_resolution->Fit(f_res_recable);
//   f_res_recable->Draw("same");

//   TF1 *f_res_ts = new TF1("f_res_ts","sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)))");
//   f_res_ts->SetLineColor(4);
//   perf_had_nopu.h_resolution->Fit(f_res_ts);
//   f_res_ts->Draw("same");

//   TF1 *f_res_nopunoint = new TF1("f_res_nopunoint","sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)))");
//   f_res_nopunoint->SetLineColor(8);
//   perf_had_nopunoint.h_resolution->Fit(f_res_nopunoint);
//   f_res_nopunoint->Draw("same");

//   /*
//   TF1 *f_res_ts_noint = new TF1("f_res_ts_noint","sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)))");
//   f_res_ts_noint->SetLineColor(8);
//   perf_had_nopunoint.h_resolution->Fit(f_res_ts_noint);
//   f_res_ts_noint->Draw("same");
//   */
  
//   /*  TCanvas *c_resolution_uncor = new TCanvas("c_resolution_uncor","c_resolution_uncor",500,500);
//   perf_had_nopu.h_resolution_cor->Draw("AP");
//   perf_had_pu.h_resolution_cor->Draw("P sames");
//   perf_had_nopu.h_resolution_cor->Draw("P sames");*/

// }


shower_reco_perf getPerformance(std::vector<TTree*> trees_, std::vector<unsigned int> energies_, TString var_, int color_, TString suffix, TString region, TString cut="0==0") {

  const int nbins = energies_.size();
  float e[nbins];
  float mean[nbins];               float mean_err[nbins];
  float sigma[nbins];              float sigma_err[nbins];
  float resp[nbins];               float resp_err[nbins];
  float resp_sigma[nbins];         float resp_sigma_err[nbins];
  float reso[nbins];               float reso_err[nbins];
  float en_reco_resolution[nbins]; float en_reco_resolution_err[nbins];

  for (unsigned int itree = 0; itree<trees_.size(); ++itree) {

    e[itree] = float(energies_[itree]);
    TString en_str_ = std::to_string(energies_[itree]);

    // Plot only reco energy
    TString nameEn = "en_"+en_str_+"_"+var_+"_"+suffix;
    TH1F *h_en;
    if (var_.Contains("reconstructable_energy")){
      h_en = new TH1F("h_"+nameEn,"h_"+nameEn,100,e[itree]-e[itree]*0.5, e[itree]+e[itree]*0.5);
    } else if(var_.Contains("uncalib_energy_")){
      h_en = new TH1F("h_"+nameEn,"h_"+nameEn,100,0., 100000.);
    }   
    h_en->SetName("h_"+nameEn);
    std::cout << " var_ = " << var_ << "\n";
    trees_[itree]->Project("h_"+nameEn,var_,cut);

//    TCanvas *c_en = new TCanvas("c_"+nameEn,"c_"+nameEn,500,500);
//    h_en->Draw("HIST E0");

    // Fit with gaussian
    h_en->Fit("gaus");
    TF1 *fit_en = h_en->GetFunction("gaus");
    fit_en->Draw("same");

    mean[itree] = fit_en->GetParameter(1);
    mean_err[itree] = fit_en->GetParError(1);
    sigma[itree] = fit_en->GetParameter(2);
    sigma_err[itree] = fit_en->GetParError(2);

    reso[itree]      = sigma[itree]/mean[itree];
    reso_err[itree]  = calcUncorUnc(sigma[itree],sigma_err[itree],mean[itree],mean_err[itree]);
    
    std::cout << "Energy: " << e[itree] << "\n";
    std::cout << "    Mean: "<< mean[itree] << ", Sigma: " << sigma[itree] << std::endl;
    std::cout << "    Resolution: " << reso[itree] << std::endl;

    // Plot response ( reco energy / generated energy )
    TString name = "resp_"+en_str_+"_"+var_+"_"+suffix;
    TH1F *h_resp = new TH1F("h_"+name+"_resp","h_"+name+"_resp",100,0.,2.);
    h_resp->SetName("h_"+name);
    //TString var_resp = "("+var_+"/cp_e)"; std::cout << " var_resp = " << var_resp << "\n";
    TString var_resp = "SumEoverEgen_" + region; 
    std::cout << " var_resp = " << var_resp << "\n";
    trees_[itree]->Project("h_"+name,var_resp,cut);

//    TCanvas *c_resp = new TCanvas("c_"+name,"c_"+name,500,500);
//    h_resp->Draw("HIST E0");

    // Fit with gaussian
    h_resp->Fit("gaus");
    TF1 *fit_resp = h_resp->GetFunction("gaus");
    fit_resp->Draw("same");

    resp[itree]              = fit_resp->GetParameter(1); 
    resp_err[itree]          = fit_resp->GetParError(1);
    resp_sigma[itree]        = fit_resp->GetParameter(2);
    resp_sigma_err[itree]    = fit_resp->GetParError(2);

    std::cout << "Energy response for E=: " << e[itree] << "\n";
    std::cout << "    Mean: "<< resp[itree] << ", Error: " << resp_err[itree] << std::endl;
    std::cout << "    Sigma: " << resp_sigma[itree] << ", Error: "<< resp_sigma_err[itree] << std::endl;


    //DA RISCRIVERE, fa cagare



//    // for response
//    int nbins = 20; float xmin = -1.*(float)energies_[itree]; float xmax = 2.*(float)energies_[itree]; //0.7
//    TH1F *h_resp = new TH1F("h_"+name,"h_"+name,nbins,xmin,xmax);
//    trees_[itree]->Project("h_"+name,var_,cut);
//    std::cout << " First print: " << name << " " << var_ << "\n";    
//
//    double error_ov=0.; double integral_ov = h_resp->IntegralAndError(nbins,nbins+1,error_ov);
//    double error_un=0.; double integral_un = h_resp->IntegralAndError(0,1,error_un);
//    h_resp->SetBinContent(nbins,integral_ov); h_resp->SetBinError(nbins,error_ov);
//    h_resp->SetBinContent(1,integral_un);  h_resp->SetBinError(1,error_un);
//
//    // for response correction
//    /*    int bins_cor = 20; 
//    TH1F *h_resp_cor = new TH1F("h_"+name+"_cor","h_"+name++"_cor",bins_cor,0.,1.5);
//    trees_[itree]->Project("h_"+name+"_cor","("+var_+"/cp_e)");
//
//    double error_ov_cor=0.; double integral_ov_cor = h_resp_cor->IntegralAndError(nbins_cor,nbins+1,error_ov_cor);
//    double error_un_cor=0.; double integral_un_cor = h_resp_cor->IntegralAndError(0,1,error_un_cor);
//    h_resp_cor->SetBinContent(nbins_cor,integral_ov_cor); h_resp_cor->SetBinError(nbins_cor,error_ov_cor);
//    h_resp_cor->SetBinContent(1,integral_un_cor);  h_resp_cor->SetBinError(1,error_un_cor);*/
//
//
//    h_resp->Fit("gaus");
//    TF1 *fit = h_resp->GetFunction("gaus");  
//
//    mean[itree] = fit->GetParameter(1); 
//    //mean[itree] = h_resp->GetMean();
//    mean_err[itree] = fit->GetParError(1);
//    
//
//    response[itree]             = mean[itree]/e[itree];
//    response_err[itree]         = mean_err[itree]/e[itree]; 
//    response_sigma[itree]       = sigma[itree]/e[itree]; 
//    resolution[itree]           = sigma[itree]/((1./response[itree])*mean_cor[itree]);
//    resolution_err[itree]       = calcUncorUnc(sigma[itree],sigma_err[itree],mean_cor[itree],mean_cor_err[itree]);
//    //resolution[itree]           = sigma[itree]/((1./response[itree]));
//    //resolution_err[itree]       = calcUncorUnc(sigma[itree],sigma_err[itree],(1./response[itree])*mean[itree],(1./response[itree])*mean_err[itree]);
//    resolution_uncor[itree]     = sigma[itree]/(mean[itree]);
//    resolution_uncor_err[itree] = calcUncorUnc(sigma[itree],sigma_err[itree],mean[itree],mean_err[itree]);

  }
  
  TGraphErrors* h_response_         = new TGraphErrors(nbins,e,resp,0,resp_err);
  TGraphErrors* h_linearity_        = new TGraphErrors(nbins,e,mean,0,mean_err);
  TGraphErrors* h_resolution_       = new TGraphErrors(nbins,e,reso,0,reso_err);

  h_response_->SetLineColor(color_);         h_response_->SetMarkerColor(color_);
  h_linearity_->SetLineColor(color_);        h_linearity_->SetMarkerColor(color_);
  h_resolution_->SetLineColor(color_);       h_resolution_->SetMarkerColor(color_);

  // Fit to extract param to rescale
  TCanvas *c_lin = new TCanvas("c_lin_"+var_, "c_lin_"+var_, 500, 500);
  c_lin->cd();
  h_linearity_->Draw("AP");

  TFitResultPtr frp_lin(0);
  auto f_lin = new TF1("myLinFunc", "[0]+[1]*x",energies_[0], energies_[trees_.size()]);
  frp_lin = h_linearity_->Fit(f_lin, "S");

  std::cout << "Linearity fit: " << "\n";
  std::cout << "    Param[0]: "<< frp_lin->Parameter(0) << ", Param[1]: " << frp_lin->Parameter(1) << std::endl;

  // Rescale
  float mean_corr[nbins];          float mean_corr_err[nbins];
  float sigma_corr[nbins];         float sigma_corr_err[nbins];
  float reso_corr[nbins];          float reso_corr_err[nbins];

  float energy_tree;
  for (unsigned int itree = 0; itree<trees_.size(); ++itree) {
    trees_[itree]->SetBranchAddress(var_,&energy_tree);
    Long64_t nentries = trees_[itree]->GetEntries();
//    for (Long64_t i=0;i<nentries;i++) {
    for (Long64_t i=0;i<10;i++) {
      trees_[itree]->GetEntry(i);
//      std::cout << energy_tree << std::endl;
    }
  }


  for (unsigned int itree = 0; itree<trees_.size(); ++itree) {

    e[itree] = float(energies_[itree]);
    TString en_str_ = std::to_string(energies_[itree]);

    // Plot only reco energy
    TString nameEn_corr = "en_"+en_str_+"_"+var_+"_"+suffix+"_corr";
    TH1F *h_en_corr = new TH1F("h_"+nameEn_corr,"h_"+nameEn_corr,100,e[itree]-e[itree]*0.5, e[itree]+e[itree]*0.5);
    h_en_corr->SetName("h_"+nameEn_corr);
    TString par0 = TString(std::to_string(frp_lin->Parameter(0)));
    TString par1 = TString(std::to_string(frp_lin->Parameter(1)));
    TString varCorr_ = var_+"/(1.*"+par1+")-1.*"+par0;
    std::cout << " varCorr_ = " << varCorr_ << "\n";
    trees_[itree]->Project("h_"+nameEn_corr,varCorr_,cut);

    TCanvas *c_en_corr = new TCanvas("c_"+nameEn_corr,"c_"+nameEn_corr,500,500);
    h_en_corr->SetMarkerColor(kRed+2);
    h_en_corr->SetLineColor(kRed+2);
    h_en_corr->Draw("HIST E0");

    // Fit with gaussian
    h_en_corr->Fit("gaus");
    TF1 *fit_en_corr = h_en_corr->GetFunction("gaus");
    fit_en_corr->Draw("same");

    mean_corr[itree]      = fit_en_corr->GetParameter(1);
    mean_corr_err[itree]  = fit_en_corr->GetParError(1);
    sigma_corr[itree]     = fit_en_corr->GetParameter(2);
    sigma_corr_err[itree] = fit_en_corr->GetParError(2);

    reso_corr[itree]      = sigma_corr[itree]/mean_corr[itree];
    reso_corr_err[itree]  = calcUncorUnc(sigma_corr[itree],sigma_corr_err[itree],mean_corr[itree],mean_corr_err[itree]);
    
    std::cout << "Rescaled Energy: " << e[itree] << "\n";
    std::cout << "    Mean: "<< mean_corr[itree] << ", Sigma: " << sigma_corr[itree] << std::endl;
    std::cout << "    Resolution: " << reso_corr[itree] << std::endl;

  }

  TGraphErrors* h_resolution_corr_ = new TGraphErrors(nbins,e,reso_corr,0,reso_corr_err);

  h_resolution_corr_->SetLineColor(color_); h_resolution_corr_->SetMarkerColor(color_);

  // Save all in shower_reco_perf
  shower_reco_perf shower_reco_perf_tmp;
  shower_reco_perf_tmp.h_response          = h_response_;
  shower_reco_perf_tmp.h_linearity         = h_linearity_;
  shower_reco_perf_tmp.h_resolution        = h_resolution_;
  shower_reco_perf_tmp.h_resolution_corr   = h_resolution_corr_;


  return shower_reco_perf_tmp;
}


shower_reco_perf getPerformancePosition(std::vector<TTree*> trees_, std::vector<float> pos_, float cpEnergy_, TString var_, int color_, TString suffix, TString cut="0==0", bool isEta=false) {

  const int nbins = pos_.size();
  float pos[nbins];
  float mean[nbins];             float mean_err[nbins];
  float mean_cor[nbins];         float mean_cor_err[nbins];
  float sigma[nbins];            float sigma_err[nbins];
  float response[nbins];         float response_err[nbins];
  float resolution[nbins];       float resolution_err[nbins];
  float resolution_uncor[nbins]; float resolution_uncor_err[nbins];

  for (unsigned int itree = 0; itree<trees_.size(); ++itree) {

    stringstream pos_str_tmp_; pos_str_tmp_ << pos_[itree];
    string pos_str_; pos_str_tmp_ >> pos_str_;
  
    //    pos_str_.replace(pos_str_.find("."), sizeof(".") -1, "p");

    TString name = "resp_"+(TString)pos_str_+"_"+var_+"_"+suffix;
    std::cout << name << "\n";
    //(string)name.replace((string)name.find("."), sizeof(".") -1, "p");
    //s.replace(s.find("$name"), sizeof("$name") - 1, "Somename");
    std::cout << name << "\n";

    // for response
    int nbins = 200; float xmin = 0.5*(float)cpEnergy_; float xmax = 1.5*(float)cpEnergy_; //0.7
    TH1F *h_resp = new TH1F("h_"+name,"h_"+name,nbins,xmin,xmax);
    std::cout << trees_[itree] << "\n";
    trees_[itree]->Project("h_"+name,var_,cut);
    std::cout << " First print: " << name << " " << var_ << " " << pos_[itree] << "\n";    

    double error_ov=0.; double integral_ov = h_resp->IntegralAndError(nbins,nbins+1,error_ov);
    double error_un=0.; double integral_un = h_resp->IntegralAndError(0,1,error_un);
    h_resp->SetBinContent(nbins,integral_ov); h_resp->SetBinError(nbins,error_ov);
    h_resp->SetBinContent(1,integral_un);  h_resp->SetBinError(1,error_un);

    // for response correction
    /*    int bins_cor = 20; 
    TH1F *h_resp_cor = new TH1F("h_"+name+"_cor","h_"+name++"_cor",bins_cor,0.,1.5);
    trees_[itree]->Project("h_"+name+"_cor","("+var_+"/cp_e)");

    double error_ov_cor=0.; double integral_ov_cor = h_resp_cor->IntegralAndError(nbins_cor,nbins+1,error_ov_cor);
    double error_un_cor=0.; double integral_un_cor = h_resp_cor->IntegralAndError(0,1,error_un_cor);
    h_resp_cor->SetBinContent(nbins_cor,integral_ov_cor); h_resp_cor->SetBinError(nbins_cor,error_ov_cor);
    h_resp_cor->SetBinContent(1,integral_un_cor);  h_resp_cor->SetBinError(1,error_un_cor);*/


    // fit with gaussian
    h_resp->Fit("gaus");
    TF1 *fit = h_resp->GetFunction("gaus");  

    // -1.*320.*TMath::Tan(2.*(TMath::ATan(TMath::Exp(1.57))))
    /*
    if (isEta)
      pos[itree] = -1.*320.*TMath::Tan(2.*(TMath::ATan(TMath::Exp(pos_[itree]))));
    else {
      pos[itree] = (float)pos_[itree];
    }
    */
    // convert R to eta
    double thetahalf = 0.5*TMath::ASin(((float)pos_[itree])/(sqrt((320.*320.)+((float)pos_[itree]*(float)pos_[itree]))));
    pos[itree] = -1.*TMath::Log( TMath::Tan(thetahalf) );
    //pos[itree] = (float)pos_[itree];

    mean[itree] = fit->GetParameter(1); 
    //mean[itree] = h_resp->GetMean();
    mean_err[itree] = fit->GetParError(1);
    

    //    stringstream resp_cor_str_; resp_cor_str_ << 1./(mean[itree]/e[itree]);
    stringstream resp_cor_str_; resp_cor_str_ << 1./(mean[itree]);
    string resp_cor_str; resp_cor_str_ >> resp_cor_str;

    //    TH1F *h_resp_cor = new TH1F("h_"+name+"_cor","h_"+name+"_cor",100,0.3,1.5);
    TH1F *h_resp_cor = new TH1F("h_"+name+"_cor","h_"+name+"_cor",200,0.05,2.);
    //    TString resp_cor = "("+var_+"/cp_e)*"+resp_cor_str; std::cout << " resp_cor = " << resp_cor << "\n";
    TString resp_cor = "("+var_+"/cp_e)"; std::cout << " resp_cor = " << resp_cor << "\n";
    trees_[itree]->Project("h_"+name+"_cor",resp_cor,cut);
    std::cout << " mean = " << h_resp->GetMean() << " " << h_resp->GetRMS()<< "\n";
    std::cout << " mean cor = " << h_resp_cor->GetMean() << " " << h_resp_cor->GetRMS() << "\n";

    h_resp_cor->Fit("gaus");
    TF1 *fit_cor = h_resp_cor->GetFunction("gaus");
    h_resp_cor->SetName("h_resp_cor_"+name);

//    TCanvas *c_resp_cor = new TCanvas("c_resp_cor_"+name,"c_resp_cor_"+name,500,500);
//    h_resp_cor->Draw("HIST E0");

    
    mean_cor[itree] = fit_cor->GetParameter(1); 
    //mean_cor[itree] = h_resp_cor->GetMean(); 
    mean_cor_err[itree] = fit_cor->GetParError(1);
    sigma[itree] = fit_cor->GetParameter(2);
    //sigma[itree] = h_resp_cor->GetRMS();
    sigma_err[itree] = fit_cor->GetParError(2);

    response[itree]             = mean[itree]/cpEnergy_;
    response_err[itree]         = mean_err[itree]/cpEnergy_; 
    resolution_err[itree]       = calcUncorUnc(sigma[itree],sigma_err[itree],mean_cor[itree],mean_cor_err[itree]);
    resolution[itree]           = sigma[itree]/((1./response[itree])*mean_cor[itree]);
    //resolution[itree]           = sigma[itree]/((1./response[itree]));
    //resolution_err[itree]       = calcUncorUnc(sigma[itree],sigma_err[itree],(1./response[itree])*mean[itree],(1./response[itree])*mean_err[itree]);
    resolution_uncor[itree]     = sigma[itree]/(mean[itree]);
    resolution_uncor_err[itree] = calcUncorUnc(sigma[itree],sigma_err[itree],mean[itree],mean_err[itree]);

  }
  
  TGraphErrors* h_response_         = new TGraphErrors(nbins,pos,response,0,response_err);
  TGraphErrors* h_linearity_        = new TGraphErrors(nbins,pos,mean,0,mean_err);
  TGraphErrors* h_resolution_       = new TGraphErrors(nbins,pos,resolution,0,resolution_err);
  TGraphErrors* h_resolution_uncor_ = new TGraphErrors(nbins,pos,resolution_uncor,0,resolution_uncor_err);;

  h_response_->SetLineColor(color_);         h_response_->SetMarkerColor(color_);
  h_linearity_->SetLineColor(color_);        h_linearity_->SetMarkerColor(color_);
  h_resolution_->SetLineColor(color_);       h_resolution_->SetMarkerColor(color_);
  h_resolution_uncor_->SetLineColor(color_); h_resolution_uncor_->SetMarkerColor(color_);
  /*
  f_lin_em->SetLineColor(1);
  f_lin_had->SetLineColor(2);
  f_lin_recable->SetLineColor(4);
  f_res_em->SetLineColor(1);
  f_res_had->SetLineColor(2);
  f_res_recable->SetLineColor(4);
  gr_lin_recable_->Fit("f_lin_recable");
  gr_lin_em_->Fit("f_lin_em");  
  gr_lin_had_->Fit("f_lin_had");  
  gr_res_recable_->Fit("f_res_recable");
  gr_res_em_->Fit("f_res_em");  
  gr_res_had_->Fit("f_res_had");  
  */

  shower_reco_perf shower_reco_perf_tmp;
  shower_reco_perf_tmp.h_response         = h_response_;
  shower_reco_perf_tmp.h_linearity        = h_linearity_;
  shower_reco_perf_tmp.h_resolution       = h_resolution_;
  shower_reco_perf_tmp.h_resolution_corr  = h_resolution_uncor_;


  return shower_reco_perf_tmp;
}


float getMeanOfBranch(TTree *tree, TString varName, int nbins, float xmin, float xmax) {

  TH1F *h_ = new TH1F("h_","h_",nbins,xmin,xmax);
  tree->Project("h_",varName);
  float mean = h_->GetMean();
  h_->Delete();
  return mean;
}


/*
std::vector<TCanvas*> plotRecoPerformance(std::vector<had_shower_reco> mc_response) {

  const int nbins = mc_response.size();
  float e[nbins];
  float mean_recable[nbins]; float mean_recable_err[nbins];
  float res_recable[nbins];  float res_recable_err[nbins];
  float mean_em[nbins]; float mean_em_err[nbins];
  float res_em[nbins];  float res_em_err[nbins];
  float mean_had[nbins]; float mean_had_err[nbins];
  float res_had[nbins];  float res_had_err[nbins];


  for (unsigned int i0=0; i0<nbins; ++i0) {
    e[i0] = (float)mc_response[i0].cp_e;
    mean_recable[i0] = mc_response[i0].mean_recable; 
    mean_recable_err[i0] = mc_response[i0].mean_recable_err;
    mean_em[i0] = mc_response[i0].mean_em; 
    mean_em_err[i0] = mc_response[i0].mean_em_err;
    mean_had[i0] = mc_response[i0].mean_had; 
    mean_had_err[i0] = mc_response[i0].mean_had_err;
    //float resp_cor = e[i0]/mean[i0];
    float resp_cor_recable = e[i0]/mean_recable[i0];
    float resp_cor_em = e[i0]/mean_em[i0];
    float resp_cor_had = e[i0]/mean_had[i0];
 
    res_recable[i0] = mc_response[i0].sigma_recable/((resp_cor_recable)*mean_recable[i0]); 
    res_recable_err[i0] = calcUncorUnc(mc_response[i0].sigma_recable,mc_response[i0].sigma_recable_err,
				       (resp_cor_recable)*mc_response[i0].mean_recable,(resp_cor_recable)*mc_response[i0].mean_recable_err);

    res_em[i0] = mc_response[i0].sigma_em/((resp_cor_em)*mean_em[i0]); 
    res_em_err[i0] = calcUncorUnc(mc_response[i0].sigma_em,mc_response[i0].sigma_em_err,(resp_cor_em)*mc_response[i0].mean_em,(resp_cor_em)*mc_response[i0].mean_em_err);

    res_had[i0] = mc_response[i0].sigma_had/((resp_cor_had)*mean_had[i0]); 
    res_had_err[i0] = calcUncorUnc(mc_response[i0].sigma_had,mc_response[i0].sigma_had_err,(resp_cor_had)*mc_response[i0].mean_had,(resp_cor_had)*mc_response[i0].mean_had_err);


  }  

  TGraphErrors *gr_lin_em_      = new TGraphErrors(nbins,e,mean_em,0,mean_em_err);
  TGraphErrors *gr_lin_had_     = new TGraphErrors(nbins,e,mean_had,0,mean_had_err);
  TGraphErrors *gr_lin_recable_ = new TGraphErrors(nbins,e,mean_recable,0,mean_recable_err);
  TGraphErrors *gr_res_em_      = new TGraphErrors(nbins,e,res_em,0,res_em_err);
  TGraphErrors *gr_res_had_     = new TGraphErrors(nbins,e,res_had,0,res_had_err);
  TGraphErrors *gr_res_recable_ = new TGraphErrors(nbins,e,res_recable,0,res_recable_err);
  
  TF1  *f_lin_recable = new TF1("f_lin_recable","pol1");
  TF1  *f_lin_em      = new TF1("f_lin_em","pol1");
  TF1  *f_lin_had     = new TF1("f_lin_had","pol1");
  TF1  *f_res_recable = new TF1("f_res_recable","sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)))");
  TF1  *f_res_em      = new TF1("f_res_em","sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)))");
  TF1  *f_res_had     = new TF1("f_res_had","sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)))");
  //TF1  *f_res = new TF1("f_res","sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)) + ([2]/x)*([2]/x))");
  //TF1  *f_res = new TF1("f_res","sqrt(([1]/sqrt(x))*([1]/sqrt(x)) + ([2]/x)*([2]/x) + [0]*[0])");
  f_lin_em->SetLineColor(1);
  f_lin_had->SetLineColor(2);
  f_lin_recable->SetLineColor(4);
  f_res_em->SetLineColor(1);
  f_res_had->SetLineColor(2);
  f_res_recable->SetLineColor(4);
  gr_lin_recable_->Fit("f_lin_recable");
  gr_lin_em_->Fit("f_lin_em");  
  gr_lin_had_->Fit("f_lin_had");  
  gr_res_recable_->Fit("f_res_recable");
  gr_res_em_->Fit("f_res_em");  
  gr_res_had_->Fit("f_res_had");  

  TCanvas *c_lin = new TCanvas("c_lin","c_lin",500,500);
  gr_lin_em_->GetYaxis()->SetTitle("Reconstructed energy [GeV]");
  gr_lin_em_->GetXaxis()->SetTitle("Generated energy [GeV]");
  gr_lin_em_->GetYaxis()->SetRangeUser(0.,350.);
  gr_lin_em_->SetMarkerColor(1); gr_lin_em_->SetMarkerSize(1.3); gr_lin_em_->SetMarkerStyle(20); gr_lin_em_->SetLineColor(1);
  gr_lin_had_->SetMarkerColor(2); gr_lin_had_->SetMarkerSize(1.3); gr_lin_had_->SetMarkerStyle(22); gr_lin_had_->SetLineColor(2);
  gr_lin_recable_->SetMarkerColor(4); gr_lin_recable_->SetMarkerSize(1.3); gr_lin_recable_->SetMarkerStyle(21); gr_lin_recable_->SetLineColor(4);
  gr_lin_em_->Draw("AP");
  gr_lin_had_->Draw("P");
  gr_lin_recable_->Draw("P");
  //f_lin_em->Draw("sames");
  //f_lin_had->Draw("sames");
  //f_lin_recable->Draw("sames");
  

  TCanvas *c_res = new TCanvas("c_res","c_res",500,500);
  gr_res_em_->GetYaxis()->SetTitle("#sigma(E)/E [Reconstructed]");
  gr_res_em_->GetXaxis()->SetTitle("Generated energy [GeV]");
  gr_res_em_->GetYaxis()->SetRangeUser(0.,0.30);
  gr_res_em_->SetMarkerColor(1); gr_res_em_->SetMarkerSize(1.3); gr_res_em_->SetMarkerStyle(20); gr_res_em_->SetLineColor(1);
  gr_res_had_->SetMarkerColor(2); gr_res_had_->SetMarkerSize(1.3); gr_res_had_->SetMarkerStyle(22); gr_res_had_->SetLineColor(2);
  gr_res_recable_->SetMarkerColor(4); gr_res_recable_->SetMarkerSize(1.3); gr_res_recable_->SetMarkerStyle(21); gr_res_recable_->SetLineColor(4);
  gr_res_em_->Draw("AP");
  gr_res_had_->Draw("P");
  gr_res_recable_->Draw("P");

  std::vector<TCanvas*> c_;
  c_.push_back(c_lin);
  c_.push_back(c_res);

  return c_;
}
*/
float calcUncorUnc(float num, float num_err, float den, float den_err) {
  return ((num/den)*sqrt( ((num_err/num)*(num_err/num)) + ((den_err/den)*(den_err/den)) ));
}

void addOvFlow(TH1F *h) {

  double error_ov=0.; double integral_ov = h->IntegralAndError(h->GetNbinsX(),h->GetNbinsX()+1,error_ov);
  double error_un=0.; double integral_un = h->IntegralAndError(0,1,error_un);
  h->SetBinContent(h->GetNbinsX(),integral_ov); h->SetBinError(h->GetNbinsX(),error_ov);
  h->SetBinContent(1,integral_un);  h->SetBinError(1,error_un);

}


TH1F *getVarRatio(TTree *tr, string en_str_, TString varNum, TString varDen, int nbins, float xmin, float xmax, TString name) {

  TString name_ = "h_"+varNum+"_"+varDen+"_"+(TString)en_str_+"_"+name;
  TH1F *h_ = new TH1F(name_,name_,nbins,xmin,xmax);
  tr->Project(name_,varNum+"/"+varDen);

  return h_;
}

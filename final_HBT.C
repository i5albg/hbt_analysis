#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TVector3.h"
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <vector>
#include <TLorentzVector.h>
#include "THnSparse.h"
#include <cstring>
#include <ctime>
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <vector>
#include <map>
#include "TFrame.h"
#include "TBenchmark.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TDatime.h"
#include "setSPRACEColors.C"
#include <stdlib.h>
#include <algorithm>	

using namespace ROOT::Math;

using namespace std;

//to reject a range in the fit -- in principle did not reject any range
Double_t reject_range_min = 0.0;
Double_t reject_range_max = 0.00001;


//Exponential function + long range
Double_t func1_exp(Double_t* x, Double_t* par){
Double_t v = 0;
if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
else{v= par[0]*(1 + par[1]*exp(-par[2]*x[0]/0.1973))*(1+par[3]*x[0]);}
return v;
}

//Gaussian function + long range
Double_t func2_gauss(Double_t* x, Double_t* par){
Double_t v = 0;
if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/0.1973,2.)))*(1+par[3]*x[0]);}
return v;
}

//Levy function 

Double_t func3_levy(Double_t* x, Double_t* par){
Double_t v = 0;
if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/0.1973,par[4])))*(1+par[3]*x[0]);}
return v;
}

void final_HBT(){

auto fileName = "/home/isadora/Documents/root-framework/HBT_analysis/HiForestAOD_DATAtest_500k.root"; //access tree here: https://www.dropbox.com/s/7yzuzy6cgpktg16/HiForestAOD_DATAtest_500k.root?dl=0
auto outFileName = "HBT_histos.root"; //output file saving the histograms, fit functions, etc...
TFile *f = TFile::Open(fileName);
TFile *output = new TFile(outFileName,"recreate");



//auto treeName= "demo/HBT";
//auto t_1 = f->Get<TTree>(treeName);
 TTree *t_1 = (TTree*)f->Get("demo/HBT");

TH1D * h_qinv_sig_SS = new TH1D("hqinvSigSS","q_{inv} SameEvt SS tracks",50,0,1);
TH1D * h_qinv_sig_OS = new TH1D("hqinvSigOS","q_{inv} SameEvt OS tracks",50,0,1);
TH1D * h_qinv_sig_SS_Ccorr = new TH1D("hqinvSigSS_Ccorr","q_{inv} SameEvt SS tracks Ccor",50,0,1); //aplying Coulomb correction
TH1D * h_qinv_sig_OS_Ccorr = new TH1D("hqinvSigOS_Ccorr","q_{inv} SameEvt OS tracks Ccor",50,0,1); //aplying Coulomb correction

std::ostringstream os; 

auto c = new TCanvas("c", "", 600, 600);
c->cd();


t_1->Draw("qinvSigSS>>hqinvSigSS","HFsumET>0 && HFsumET<900");
h_qinv_sig_SS->Write();

h_qinv_sig_SS->GetXaxis()->SetTitle("q_{inv} [GeV]");
h_qinv_sig_SS->GetYaxis()->SetTitle("Number of pairs / bin");
h_qinv_sig_SS->GetXaxis()->SetTitleOffset(1.15);
h_qinv_sig_SS->GetYaxis()->SetTitleOffset(1.15);
h_qinv_sig_SS->GetYaxis()->SetNdivisions(810);
h_qinv_sig_SS->GetXaxis()->SetRangeUser(0.0,0.5);

gPad->SetLogy(1);

c->Update();


t_1->Draw("qinvSigOS>>hqinvSigOS","HFsumET>0 && HFsumET<900");	  
h_qinv_sig_OS->Write();
c->Update();


t_1->Draw("qinvSigSS>>hqinvSigSS_Ccorr","coulombWSS*(HFsumET>0 && HFsumET<900)");
h_qinv_sig_SS_Ccorr->Write();
c->Update();

t_1->Draw("qinvSigOS>>hqinvSigOS_Ccorr","coulombWOS*(HFsumET>0 && HFsumET<900)");
h_qinv_sig_OS_Ccorr->Write();
c->Update();

TH1D * h_qinv_sig_SS_Ccorr_clone = (TH1D *)h_qinv_sig_SS_Ccorr->Clone("h_qinv_sig_SS_Ccorr_clone");
TH1D * h_qinv_sig_OS_Ccorr_clone = (TH1D *)h_qinv_sig_OS_Ccorr->Clone("h_qinv_sig_OS_Ccorr_clone");

Int_t bin_for_normInt_min = h_qinv_sig_SS_Ccorr_clone->GetXaxis()->FindBin(0.0999); //interval for normalization of qinv distributions
Int_t bin_for_normInt_max = h_qinv_sig_SS_Ccorr_clone->GetXaxis()->FindBin(0.9999);
Double_t int_num_controlRegion = h_qinv_sig_SS_Ccorr_clone->Integral(bin_for_normInt_min,bin_for_normInt_max);
Double_t int_den_controlRegion = h_qinv_sig_OS_Ccorr_clone->Integral(bin_for_normInt_min,bin_for_normInt_max);

h_qinv_sig_OS_Ccorr_clone->Scale(int_num_controlRegion/int_den_controlRegion);
h_qinv_sig_SS_Ccorr_clone->Divide(h_qinv_sig_OS_Ccorr_clone); //Single ratio: SS Divided by OS
h_qinv_sig_SS_Ccorr_clone->SetStats(0);
h_qinv_sig_SS_Ccorr_clone->SetTitle("");
h_qinv_sig_SS_Ccorr_clone->SetLineColor(1);
h_qinv_sig_SS_Ccorr_clone->SetMarkerColor(1);
h_qinv_sig_SS_Ccorr_clone->SetMarkerStyle(20);
h_qinv_sig_SS_Ccorr_clone->SetMarkerSize(0.65);
h_qinv_sig_SS_Ccorr_clone->GetXaxis()->SetRangeUser(0.0,.99);
h_qinv_sig_SS_Ccorr_clone->GetYaxis()->SetRangeUser(0.95,1.405);
h_qinv_sig_SS_Ccorr_clone->GetYaxis()->SetNdivisions(805);
h_qinv_sig_SS_Ccorr_clone->GetXaxis()->SetNdivisions(805);
h_qinv_sig_SS_Ccorr_clone->GetYaxis()->SetTitle("C_{BE}(q)");
h_qinv_sig_SS_Ccorr_clone->GetXaxis()->SetTitle("q [GeV]");
h_qinv_sig_SS_Ccorr_clone->GetXaxis()->SetTitleOffset(1.15);
h_qinv_sig_SS_Ccorr_clone->GetYaxis()->SetTitleOffset(1.2);
h_qinv_sig_SS_Ccorr_clone->Write();


//Setting function for double ratio
TF1 *f_exp = new TF1("f_exp",func1_exp,0.01, 0.99, 4); //try other fits ranges
f_exp->SetParameters(1.0,1.0,4.0,0.0); //initialize parameters
f_exp->SetParName(0,"Const");
f_exp->SetParName(1,"#lambda");
f_exp->SetParName(2,"R (fm)");
f_exp->SetParName(3,"#epsilon");
f_exp->SetLineColor(SPcolors[SPred]); 
f_exp->SetLineWidth(2);  
//h_qinv_sig_SS_Ccorr_clone->Draw("pError");
h_qinv_sig_SS_Ccorr_clone->GetYaxis()->SetTitle("C_{BE}(q_{inv})");
h_qinv_sig_SS_Ccorr_clone->GetXaxis()->SetTitle("q_{inv} [GeV]");
TFitResultPtr res_exp;
ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
res_exp = h_qinv_sig_SS_Ccorr_clone->Fit(f_exp, "S R"); 
f_exp->Draw("SameL");

TF1 *f_gauss = new TF1("f_gauss",func2_gauss,0.01, 0.99, 4); //try other fits ranges
	f_gauss->SetParameters(1.0,1.0,4.0,0.0); //initialize parameters
	f_gauss->SetParName(0,"Const");
	f_gauss->SetParName(1,"#lambda"); 
	f_gauss->SetParName(2,"R (fm)");
	f_gauss->SetParName(3,"#epsilon");
	f_gauss->SetLineColor(SPcolors[SPdgreen]); 
	f_gauss->SetLineStyle(9);
	f_gauss->SetLineWidth(3);
	TFitResultPtr res_gauss;
	//h_qinv_sig_SS_Ccorr_clone->Fit(f_gauss, "S R Q N");
	res_gauss = h_qinv_sig_SS_Ccorr_clone->Fit(f_gauss, "S R M");
	f_gauss->Draw("SameL");
	gPad->SetLogy(0); 

TF1 *f_levy = new TF1("f_levy",func3_levy,0.01, 0.99, 5); //try other fits ranges
	f_levy->SetParameters(1,1.,4,0,1.0); //initialize parameters
	f_levy->SetParName(0,"Const");
//f_levy->SetParLimits(0,0.0,2.0);
	f_levy->SetParName(1,"#lambda"); 
    f_levy->SetParLimits(1,0.0,2.0);
	f_levy->SetParName(2,"R (fm)");
	f_levy->SetParName(3,"#epsilon");
     f_levy->SetParLimits(4,1.0,2.0);
	f_levy->SetParName(4,"#aplha");
	f_levy->SetLineColor(SPcolors[SPblue]); 
	f_levy->SetLineWidth(4);
	f_levy->SetLineStyle(7);
	TFitResultPtr res_levy;
	res_levy = h_qinv_sig_SS_Ccorr_clone->Fit(f_levy, "S R");
	f_levy->Draw("SameL");
	gPad->SetLogy(0);


	f_exp->Write();
	f_gauss->Write();
	f_levy->Write();
  

	TLine *line_at_one = new TLine(1.,1.,1.,1.);
	line_at_one->SetLineColor(1);
	line_at_one->SetLineWidth(2);
	line_at_one->SetLineStyle(3); 
	line_at_one->SetX1(0);
	line_at_one->SetX2(0.99); 
	line_at_one->Draw(); 

 	TLatex *text2 = new TLatex();
   text2->SetNDC();
   text2->SetTextColor(kBlack);
   text2->SetTextFont(42);
   text2->SetTextSize(0.04);
   text2->DrawLatex(0.15, 0.82, "Centrality 35-100%");
   
   TLatex *text = new TLatex();
   text->SetNDC();
   text->SetTextColor(kBlack);
   text->SetTextFont(62);  // Define a fonte do texto (opcional)
   text->SetTextSize(0.025);  // Define o tamanho do texto (opcional)
   text->DrawLatex(0.15, 0.86, "CMS Open Data PbPb #sqrt{s_{NN}} = 2.76 TeV");
   
        ////////////INICIO//////
	auto lexp = new TLegend(0.51,0.7,0.81,0.8);
	lexp->SetBorderSize(0);
	lexp->SetTextSize(0.030);
	lexp->AddEntry(f_exp,"Exponential fit","l");
    lexp->AddEntry(f_exp, Form("R = %.4f #pm %.4f fm", f_exp->GetParameter(2), f_exp->GetParError(2)), "");
    lexp->AddEntry(f_exp, Form("#epsilon = %.4f #pm %.4f", f_exp->GetParameter(3), f_exp->GetParError(3)), "");

	
	auto lgauss = new TLegend(0.52,0.54,0.68,0.68);
	lgauss->SetBorderSize(0);
	lgauss->SetTextSize(0.030);
	lgauss->AddEntry(f_gauss,"Gaussian fit","l");
	lgauss->AddEntry(f_gauss, Form("R = %.4f #pm %.4f fm", f_gauss->GetParameter(2), f_gauss->GetParError(2)), "");
    lgauss->AddEntry(f_gauss, Form("#epsilon = %.4f #pm %.4f", f_gauss->GetParameter(3), f_gauss->GetParError(3)), "");
    lgauss->AddEntry(f_gauss, Form("#lambda = %.4f #pm %.4f", f_gauss->GetParameter(1), f_gauss->GetParError(1)), "");


  	auto llevy = new TLegend(0.15,0.6,0.38,0.8);
	llevy->SetBorderSize(0);
	llevy->SetTextSize(0.030);
	llevy->AddEntry(f_levy,"Levy fit","l");
	llevy->AddEntry(f_levy, Form("R = %.4f #pm %.4f fm", f_levy->GetParameter(2), f_levy->GetParError(2)), "");
    llevy->AddEntry(f_levy, Form("#epsilon = %.4f #pm %.4f", f_levy->GetParameter(3), f_levy->GetParError(3)), "");
    llevy->AddEntry(f_levy, Form("#lambda= %.4f #pm %.4f", f_levy->GetParameter(1), f_levy->GetParError(1)),"");
    llevy->AddEntry(f_levy, Form("#alpha = %.4f #pm %.4f", f_levy->GetParameter(4), f_levy->GetParError(4)), "");
      ////////////FIM///////////////////

	//tex_open->Draw();                                                                                                            
	//tex_ener->Draw();
	//tex->Draw();
	lexp->Draw();
	lgauss->Draw();
	llevy->Draw();
	c->Update();

}	
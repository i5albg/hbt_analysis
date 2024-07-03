/********************************************************************
 *                  	   Macro HBT All Centralitie	              *
    *
 *******************************************************************/
 ///To run it do: root -l -b -q macro_HBT_AllCentralities.C
 
 
/*******************************************************************			        
 *                               INCLUDES                          *  
 *******************************************************************/
#include <iostream>
#include <string>
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

//Intervalo para rejeitar em um ajuste - inicialmente nenhum intervalo é rejeitado
Double_t reject_range_min = 0.0;
Double_t reject_range_max = 0.00001;

//Exponential function + long range
Double_t func1_exp(Double_t* x, Double_t* par){
Double_t v = 0;
if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
else{v= par[0]*(1 + par[1]*exp(-par[2]*x[0]/0.1973))*(1+par[3]*x[0]);}
return v;
}

Double_t func2_gauss(Double_t* x, Double_t* par){
Double_t v = 0;
if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/0.1973,2.)))*(1+par[3]*x[0]);}
return v;
}

//Levy function 
//Double_t func3_levy(Double_t* x, Double_t* par){
//Double_t v = 0;
//if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
//else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/0.1973,par[4])))*(1+par[3]*x[0]);}
//return v;
//}

//Personalizar Histogramas
void func_hist_custom_qinv(TH1D *h){
   h->GetXaxis()->SetTitle("q_{inv} [GeV]");
   h->GetYaxis()->SetTitle("Number of pairs / bin");
   h->GetXaxis()->CenterTitle(); // Centraliza o título do eixo x
   h->GetYaxis()->CenterTitle(); // Centraliza o título do eixo y

}

void func_hist_custom_sr(TH1D *h){
   h->GetXaxis()->SetTitle("q_{inv} [GeV]");
   h->GetYaxis()->SetTitle("Single Ratio / (0.02 Gev)");
   h->SetLineColor(4); // Define a cor da linha para aul
   h->SetLineWidth(1); // Define a largura da linha para 2 pixels
   h->SetLineStyle(3); // Define o estilo da linha
   h->SetMarkerColor(kBlue); // Define a cor dos marcadores para aul
   h->SetMarkerSize(0.5); // Defina um valor menor que o padrão (geralmente 2)
   h->SetMarkerStyle(4); // Define o estilo dos marcadores para círculos  
   h->SetStats(0);
   h->GetXaxis()->CenterTitle(); // Centraliza o título do eixo x
   h->GetYaxis()->CenterTitle(); // Centraliza o título do eixo y

}


//main 
void all_cent_teste(){
// open file and read the tree
auto fileName = "/home/isadora/Documents/root-framework/HBT_analysis/HiForestAOD_DATA_HBT_withMixing_500k.root"; 

auto outFileName = "HBT_histos_AllCentralities.root"; //arquivo de saída salvando os histogramas, funções de ajuste, etc...

TFile *f = TFile::Open(fileName);//Abre o arquivo ROOT especificado pelo nome do arquivo e atribui o ponteiro TFile a f. Isso permite acessar o conteúdo do arquivo.

TFile *output = new TFile(outFileName,"recreate");//Cria um novo arquivo ROOT de saída com o nome especificado por outFileName. O modo "recreate" indica que se o arquivo já existir, ele será recriado, excluindo qualquer conteúdo anterior.

auto treeName = "demo/HBT";//Define o nome da árvore (tree) que será lida a partir do arquivo ROOT. Neste caso, o nome da árvore é "demo/HBT".

auto t_1 = f->Get<TTree>(treeName);//Obtém a árvore (tree) especificada pelo nome treeName a partir do arquivo ROOT f e atribui o ponteiro TTree a t_1. Isso permite acessar os dados da árvore para análise e processamento posterior.

/******************************************************/
//     Cria diretórios para salvar em formato PDF.
/******************************************************/
std::string folder_name="sig_SS_Corr";
gSystem->Exec(Form("mkdir %s",folder_name.c_str()));

std::string _name="bkg_OS_Corr";
gSystem->Exec(Form("mkdir %s",_name.c_str()));

std::string folder_ ="SS_Over_OS";
gSystem->Exec(Form("mkdir %s",folder_.c_str()));

std::string fol_ ="funcGaussName_sr_SS_Over_OS";
gSystem->Exec(Form("mkdir %s",fol_.c_str()));

std::string f_ ="funcExpName_sr_SS_Over_OS";
gSystem->Exec(Form("mkdir %s",f_.c_str()));

/********************************************************/
/*Essa é uma descrição dos intervalos de centrality e (soma da energia transversa do Forward Hadron) correspondentes, conforme apresentado na Figura 1 do artigo "https://arxiv.org/pdf/1107.4800.pdf". Aqui está a correspondência dos intervalos de centrality e HFsumET mencionados:

0-5% (centrality) -> 3380-5000GeV (HFsumET)
5-10% -> 2770-3380GeV
10-15% -> 2280-2770GeV
15-20% -> 1850-2280GeV
20-25% -> 1500-1850GeV
25-30% -> 1210-1500GeV
30-35% -> 900-1210GeV
35-60% -> 240-900GeV
60-100% -> 0.0-240GeV
Durante a análise, você pode mesclar alguns desses intervalos, por exemplo, 10-20% pode ser mesclado para 1850-2770GeV, já que a mistura foi feita seguindo essa divisão de intervalos. 

35-100% -> 0.0-900GeV

*/
/***********************************************************/

//interval for normalization of qinv distributions
Double_t qinv_min = 0.099;
Double_t qinv_max = 0.999;

//Fit range
Double_t fit_qinv_min = 0.02;
Double_t fit_qinv_max = 1.0;//Double_t fit_qinv_max = 0.50;

//Histo nBins and range
Int_t histNbins=50;// 
Float_t histLowerRange=0;
Float_t histUpperRange=1; //maximum is 1, since saved only till this value in the tree to save space

//Declaração de ponteiros para objetos do tipo TH1D e TF1	
const unsigned int nCentBins =10;
TH1D * h_qinv_sig_SS[nCentBins];
TH1D * h_qinv_sig_SS_Corr[nCentBins];
TH1D * h_qinv_bkg_OS[nCentBins];
TH1D * h_qinv_bkg_OS_Corr[nCentBins];
TH1D * h_sr_SS_Over_OS[nCentBins];
TF1  *f_exp_sr_SS_Over_OS[nCentBins];
TF1  *f_gauss_sr_SS_Over_OS[nCentBins];

std::ostringstream os;
 
//definem os intervalos de centrality e os cortes correspondentes na variável HFsumET (soma da energia na região de Forward Hadron Calorimeter).
TString cent_bins[nCentBins]={"0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-60%","60-100%","35-100%"};
TCut hFsumEtcut[nCentBins]={"HFsumET>3380 && HFsumET<5000","HFsumET>2770 && HFsumET<3380","HFsumET>2280 && HFsumET<2770","HFsumET>1850 && HFsumET<2280","HFsumET>1500 && HFsumET<1850","HFsumET>1210 && HFsumET<1500","HFsumET>900 && HFsumET<1210","HFsumET>240 && HFsumET<900","HFsumET>0 && HFsumET<240","HFsumET>0 && HFsumET<900"};


//Correção da interação Coulombiana
TCut hFsumEtcut_coulombWSS[nCentBins]={"coulombWSS*(HFsumET>3380 && HFsumET<5000)","coulombWSS*(HFsumET>2770 && HFsumET<3380)","coulombWSS*(HFsumET>2280 && HFsumET<2770)","coulombWSS*(HFsumET>1850 && HFsumET<2280)","coulombWSS*(HFsumET>1500 && HFsumET<1850)","coulombWSS*(HFsumET>1210 && HFsumET<1500)","coulombWSS*(HFsumET>900 && HFsumET<1210)","coulombWSS*(HFsumET>240 && HFsumET<900)","coulombWSS*(HFsumET>0 && HFsumET<240)","coulombWSS*(HFsumET>0 && HFsumET<900)"};
TCut hFsumEtcut_coulombWOS[nCentBins]={"coulombWOS*(HFsumET>3380 && HFsumET<5000)","coulombWOS*(HFsumET>2770 && HFsumET<3380)","coulombWOS*(HFsumET>2280 && HFsumET<2770)","coulombWOS*(HFsumET>1850 && HFsumET<2280)","coulombWOS*(HFsumET>1500 && HFsumET<1850)","coulombWOS*(HFsumET>1210 && HFsumET<1500)","coulombWOS*(HFsumET>900 && HFsumET<1210)","coulombWOS*(HFsumET>240 && HFsumET<900)","coulombWOS*(HFsumET>0 && HFsumET<240)","coulombWOS*(HFsumET>0 && HFsumET<900)"};

//Armazenam os nomes dos histogramas que serão criados para cada classe de centrality.
std::string histName_sig_SS[nCentBins]={"hqinvSigSS_0to5","hqinvSigSS_5to10","hqinvSigSS_10to15","hqinvSigSS_15to20","hqinvSigSS_20to25","hqinvSigSS_25to30","hqinvSigSS_30to35","hqinvSigSS_35to60","hqinvSigSS_60to100","hqinvSigSS_35to100"};
std::string histName_bkg_OS[nCentBins]={"hqinvBkgOS_0to5","hqinvBkgOS_5to10","hqinvBkgOS_10to15","hqinvBkgOS_15to20","hqinvBkgOS_20to25","hqinvBkgOS_25to30","hqinvBkgOS_30to35","hqinvBkgOS_35to60","hqinvBkgOS_60to100","hqinvBkgOS_35to100"};
std::string histName_sig_SS_Corr[nCentBins]={"hqinvSigSS_Corr_0to5","hqinvSigSS_Corr_5to10","hqinvSigSS_Corr_10to15","hqinvSigSS_Corr_15to20","hqinvSigSS_Corr_20to25","hqinvSigSS_Corr_25to30","hqinvSigSS_Corr_30to35","hqinvSigSS_Corr_35to60","hqinvSigSS_Corr_60to100","hqinvSigSS_Corr_35to100"};
std::string histName_bkg_OS_Corr[nCentBins]={"hqinvBkgOS_Corr_0to5","hqinvBkgOS_Corr_5to10","hqinvBkgOS_Corr_10to15","hqinvBkgOS_Corr_15to20","hqinvBkgOS_Corr_20to25","hqinvBkgOS_Corr_25to30","hqinvBkgOS_Corr_30to35","hqinvBkgOS_Corr_35to60","hqinvBkgOS_Corr_60to100","hqinvBkgOS_Corr_35to100"};


//Armazenam os títulos dos histogramas para cada classe de centrality.
TString histTitle_sig_SS[nCentBins]={"q_{inv} SameEvt SS tracks, 0-5%","q_{inv} SameEvt SS tracks, 5-10%","q_{inv} SameEvt SS tracks, 10-15%","q_{inv} SameEvt SS tracks, 15-20%","q_{inv} SameEvt SS tracks, 20-25%","q_{inv} SameEvt SS tracks, 25-30%","q_{inv} SameEvt SS tracks, 30-35%","q_{inv} SameEvt SS tracks, 35-60%","q_{inv} SameEvt SS tracks, 60-100%","q_{inv} SameEvt SS tracks, 35-100%"};

TString histTitle_bkg_OS[nCentBins]={"q_{inv} OpposiSig OS tracks, 0-5%","q_{inv} OpposiSig OS tracks, 5-10%","q_{inv} OpposiSig OS tracks, 10-15%","q_{inv} OpposiSig OS tracks, 15-20%","q_{inv} OpposiSig OS tracks, 20-25%","q_{inv} OpposiSig OS tracks, 25-30%","q_{inv} OpposiSig OS tracks, 30-35%","q_{inv} OpposiSig OS tracks, 35-60%","q_{inv} OpposiSig OS tracks, 60-100%","q_{inv} OpposiSig OS tracks, 35-100%"};

TString histTitle_sig_SS_Corr[nCentBins]={"q_{inv} SameEvt SS_Corr tracks, 0-5%","q_{inv} SameEvt SS_Corr tracks, 5-10%","q_{inv} SameEvt SS_Corr tracks, 10-15%","q_{inv} SameEvt SS_Corr tracks, 15-20%","q_{inv} SameEvt SS_Corr tracks, 20-25%","q_{inv} SameEvt SS_Corr tracks, 25-30%","q_{inv} SameEvt SS_Corr tracks, 30-35%","q_{inv} SameEvt SS_Corr tracks, 35-60%","q_{inv} SameEvt SS_Corr tracks, 60-100%","q_{inv} SameEvt SS_Corr tracks, 35-100%"};

TString histTitle_bkg_OS_Corr[nCentBins]={"q_{inv} OpposiSig OS_Corr tracks, 0-5%","q_{inv} OpposiSig OS_Corr tracks, 5-10%","q_{inv} OpposiSig OS_Corr tracks, 10-15%","q_{inv} OpposiSig OS_Corr tracks, 15-20%","q_{inv} OpposiSig OS_Corr tracks, 20-25%","q_{inv} OpposiSig OS_Corr tracks, 25-30%","q_{inv} OpposiSig OS_Corr tracks, 30-35%","q_{inv} OpposiSig OS_Corr tracks, 35-60%","q_{inv} OpposiSig OS_Corr tracks, 60-100%","q_{inv} OpposiSig OS_Corr tracks, 35-100%"};

std::string histName_sr_SS_Over_OS[nCentBins]={"hsrSSoverOS_0to5","hsrSSoverOS_5to10","hsrSSoverOS_10to15","hsrSSoverOS_15to20","hsrSSoverOS_20to25","hsrSSoverOS_25to30","hsrSSoverOS_30to35","hsrSSoverOS_35to60","hsrSSoverOS_60to100","hsrSSoverOS_35to100"};
TString histTitle_sr_SS_Over_OS[nCentBins]={"Centrality: 0-5%","Centrality: 5-10%","Centrality: 10-15%","Centrality: 15-20%","Centrality: 20-25%","Centrality: 25-30%","Centrality: 30-35%","Centrality: 35-60%","Centrality: 60-100%","Centrality: 35-100%"};


TString funcExpName_sr_SS_Over_OS[nCentBins]={"fExpSrSSoverOS_0to5","fExpSrSSoverOS_5to10","fExpSrSSoverOS_10to15","fExpSrSSoverOS_15to20","fExpSrSSoverOS_20to25","fExpSrSSoverOS_25to30","fExpSrSSoverOS_30to35","fExpSrSSoverOS_35to60","fExpSrSSoverOS_60to100","fExpSrSSoverOS_35to100"};
TString funcGaussName_sr_SS_Over_OS[nCentBins]={"fGaussSrSSoverOS_0to5","fGaussSrSSoverOS_5to10","fGaussSrSSoverOS_10to15","fGaussSrSSoverOS_15to20","fGaussSrSSoverOS_20to25","fGaussSrSSoverOS_25to30","fGaussSrSSoverOS_30to35","fGaussSrSSoverOS_35to60","fGaussSrSSoverOS_60to100","fGaussSrSSoverOS_35to100"};


///Canvas will be updated on-the-fly

auto c = new TCanvas("c", "c", 600, 600);
c->cd();

TLatex *tex_cms = new TLatex();
tex_cms->SetNDC();
tex_cms->SetTextColor(kBlack);
tex_cms->SetTextFont(62);  // Define a fonte do texto (opcional)
tex_cms->SetTextSize(0.025);  // Define o tamanho do texto (opcional)
tex_cms->DrawLatex(0.15, 0.82, "CMS Open Data PbPb #sqrt{s_{NN}} = 2.76 TeV");

//Loop in all  bins - do everything: create histos, single ratios and fits
for(unsigned int i_centbin=0; i_centbin<nCentBins; i_centbin++){

//Estão sendo criados objetos histograma e funções ajuste para cada classe de centralidade

   h_qinv_sig_SS[i_centbin] = new TH1D((TString)histName_sig_SS[i_centbin],"",histNbins,histLowerRange,histUpperRange);
   h_qinv_sig_SS_Corr[i_centbin] = new TH1D((TString)histName_sig_SS_Corr[i_centbin],"",histNbins,histLowerRange,histUpperRange);
   
   h_qinv_bkg_OS[i_centbin] = new TH1D((TString)histName_bkg_OS[i_centbin],"",histNbins,histLowerRange,histUpperRange);
   h_qinv_bkg_OS_Corr[i_centbin] = new TH1D((TString)histName_bkg_OS_Corr[i_centbin],"",histNbins,histLowerRange,histUpperRange);
   
   h_sr_SS_Over_OS[i_centbin] = new TH1D((TString)histName_sr_SS_Over_OS[i_centbin],"",histNbins,histLowerRange,histUpperRange);
   
   f_exp_sr_SS_Over_OS[i_centbin] = new TF1(funcExpName_sr_SS_Over_OS[i_centbin],func1_exp,fit_qinv_min,fit_qinv_max,4);  
   f_gauss_sr_SS_Over_OS[i_centbin] = new TF1(funcGaussName_sr_SS_Over_OS[i_centbin],func2_gauss,fit_qinv_min,fit_qinv_max,4);
   

/****************Preenchimento dos histogramas*********************/
/***************************************************************** 
 *                     q_inv distributions                       *
 *****************************************************************/ 
   t_1->Draw(Form("qinvSigSS>>%s",histName_sig_SS[i_centbin].c_str()),hFsumEtcut[i_centbin],"goff");
   func_hist_custom_qinv(h_qinv_sig_SS[i_centbin]);
   h_qinv_sig_SS[i_centbin]->Write();
   h_qinv_sig_SS[i_centbin]->Draw();
   h_qinv_sig_SS[i_centbin]->GetXaxis()->SetTitleOffset(1.15);
	h_qinv_sig_SS[i_centbin]->GetYaxis()->SetTitleOffset(1.15);
	h_qinv_sig_SS[i_centbin]->GetYaxis()->SetNdivisions(810);
	h_qinv_sig_SS[i_centbin]->GetXaxis()->SetRangeUser(0.0,0.99);
   h_qinv_sig_SS[i_centbin]->SetLineColor(kRed);
   h_qinv_sig_SS[i_centbin]->SetLineWidth(2); 
   h_qinv_sig_SS[i_centbin]->SetLineStyle(7); // Define o estilo da linha
   h_qinv_sig_SS[i_centbin]->SetStats(0);
 
   gPad->SetLogy(1);
   c->Update(); 
   

   t_1->Draw(Form("qinvSigSS>>%s",histName_sig_SS_Corr[i_centbin].c_str()),hFsumEtcut_coulombWSS[i_centbin],"goff");
   func_hist_custom_qinv(h_qinv_sig_SS_Corr[i_centbin]);
   h_qinv_sig_SS_Corr[i_centbin]->Write();
   h_qinv_sig_SS_Corr[i_centbin]->Draw("same");
   h_qinv_sig_SS_Corr[i_centbin]->SetLineColor(kBlack);
   h_qinv_sig_SS_Corr[i_centbin]->SetLineWidth(2); 
   //h_qinv_sig_SS_Corr[i_centbin]->SetMarkerColor(kBlack);
   h_qinv_sig_SS_Corr[i_centbin]->SetStats(0);
      //Marcadores dos eixos  x e y
   c->SetTickx(1); 
   c->SetTicky(1);
   c->Update();
   
   TLatex *tex_SSOS = new TLatex();
   tex_SSOS->SetNDC();
   tex_SSOS->SetTextColor(kBlack);
   tex_SSOS->SetTextFont(42);
   tex_SSOS->SetTextSize(0.02);
   tex_SSOS->DrawLatex(0.15, 0.78, histTitle_sr_SS_Over_OS[i_centbin]);

   
   //Cria texto dentro do gráfico - lado superior esquerdo


   
   //Cria texto dentro do gráfico - lado superior esquerdo
   //TLatex *text7 = new TLatex();
   //text7->SetNDC();
   //text7->SetTextColor(kBlack);
   //text7->SetTextFont(62);  // Define a fonte do texto (opcional)
   //text7->SetTextSize(0.02);  // Define o tamanho do texto (opcional)
   //text7->DrawLatex(0.52, 0.92, "q_{inv} Same Event: Same Sign(SS_{++,--})");
   //text7->DrawLatex(0.53, 0.91, "Method: Pairs of opposite charges(SS_{++,--})");
   

   //Legenda
   TLegend *leg_uCcorr = new TLegend(0.6, 0.65, 0.9, 0.85);
   //leg_uCcorr->AddEntry(h_qinv_sig_SS[i_centbin], "Without Coulomb correction", "l");
   leg_uCcorr->AddEntry(h_qinv_sig_SS[i_centbin], "Uncorrected - Coulomb", "l");
   // Definindo o tamanho da fonte da legenda  
   leg_uCcorr->SetTextSize(0.02); 
   leg_uCcorr->SetBorderSize(0);
   leg_uCcorr->SetFillColorAlpha(0, 0);
   leg_uCcorr->Draw();
   
   
   //Legenda
   TLegend *leg_Ccorr = new TLegend(0.6, 0.62, 0.9, 0.78);
   leg_Ccorr->AddEntry(h_qinv_sig_SS_Corr[i_centbin], "Corrected - Coulomb", "l");
   // Definindo o tamanho da fonte da legenda  
   leg_Ccorr->SetTextSize(0.02); 
   leg_Ccorr->SetBorderSize(0);
   leg_Ccorr->SetFillColorAlpha(0, 0);
   leg_Ccorr->Draw();
   
/****************************************************************/    
   c->SaveAs(Form("%s/%s_%s.pdf",folder_name.c_str(),histName_sig_SS[i_centbin].c_str(),histName_sig_SS_Corr[i_centbin].c_str()));
   c->Update();
   
/**************************************************************/

   t_1->Draw(Form("qinvSigOS>>%s",histName_bkg_OS[i_centbin].c_str()),hFsumEtcut[i_centbin],"goff"); //I use "qinvSigOS" but this was not a proper name saved in the tree...should be BkgOS
   func_hist_custom_qinv(h_qinv_bkg_OS[i_centbin]);
   h_qinv_bkg_OS[i_centbin]->Write();
   h_qinv_bkg_OS[i_centbin]->Draw();
   
   //Marcadores dos eixos  x e y
   c->SetTickx(); 
   c->SetTicky();
   
   h_qinv_bkg_OS[i_centbin]->SetLineColor(kRed);
   h_qinv_bkg_OS[i_centbin]->SetLineStyle(7); // Define o estilo da linha
   h_qinv_bkg_OS[i_centbin]->SetLineWidth(2); 
   h_qinv_bkg_OS[i_centbin]->SetStats(0);
   //c->Update();
  
   t_1->Draw(Form("qinvSigOS>>%s",histName_bkg_OS_Corr[i_centbin].c_str()),hFsumEtcut_coulombWOS[i_centbin],"goff"); //I use "qinvSigOS" but this was not a proper name saved in the tree...should be BkgOS
   func_hist_custom_qinv(h_qinv_bkg_OS_Corr[i_centbin]);
   h_qinv_bkg_OS_Corr[i_centbin]->Write();
   h_qinv_bkg_OS_Corr[i_centbin]->Draw("same");
   h_qinv_bkg_OS_Corr[i_centbin]->SetLineColor(kBlack);
   h_qinv_bkg_OS_Corr[i_centbin]->SetLineWidth(2); 
   //h_qinv_bkg_OS_Corr[i_centbin]->SetMarkerColor(2);
   h_qinv_bkg_OS_Corr[i_centbin]->SetStats(0);
   gPad->SetLogy(1);
   
   
  // TLatex *text6 = new TLatex();
   //text6->SetNDC();
   ///text6->SetTextColor(kBlack);
   //text6->SetTextFont(42);
  // text6->SetTextSize(0.02);
   //text6->DrawLatex(0.15, 0.78,histTitle_sr_SS_Over_OS[i_centbin]);
   
   
  
    //Cria texto dentro do gráfico - lado superior esquerdo
   
   //text5->DrawLatex(0.15, 0.92, " q_{inv} Same Event(OS_{+-}) track");
   //text5->DrawLatex(0.55, 0.91, "Method: Opposite Sign(OS_{+- 


    //Legenda
   TLegend *leg_UCcorbkg = new TLegend(0.6, 0.65, 0.9, 0.85);
   leg_UCcorbkg->AddEntry(h_qinv_bkg_OS[i_centbin], "Uncorrected - Coulomb", "l");
   // Definindo o tamanho da fonte da legenda  
   leg_UCcorbkg->SetTextSize(0.02); 
   leg_UCcorbkg->SetBorderSize(0);
   leg_UCcorbkg->SetFillColorAlpha(0, 0);
   leg_UCcorbkg->Draw();
   

    //Legenda
   TLegend *leg_Ccorbkg = new TLegend(0.6, 0.62, 0.9, 0.78);
   leg_Ccorbkg->AddEntry(h_qinv_bkg_OS_Corr[i_centbin], "Corrected - Coulomb", "l");
   // Definindo o tamanho da fonte da legenda  
   leg_Ccorbkg->SetTextSize(0.02); 
   leg_Ccorbkg->SetBorderSize(0);
   leg_Ccorbkg->SetFillColorAlpha(0, 0);
   leg_Ccorbkg->Draw();
   
   c->SaveAs(Form("%s/%s_%s.pdf",_name.c_str(),histName_bkg_OS[i_centbin].c_str(),histName_bkg_OS_Corr[i_centbin].c_str()));
   c->Update();   
    
   
   
/***************************************************************** 
 *                     SINGLE RATIO                              *
 *****************************************************************/   
   Int_t bin_for_normInt_min = h_qinv_sig_SS_Corr[i_centbin]->GetXaxis()->FindBin(qinv_min); //interval for normalization of qinv distributions
   Int_t bin_for_normInt_max = h_qinv_sig_SS_Corr[i_centbin]->GetXaxis()->FindBin(qinv_max);
   
   //Não entendi essa parte, posso apagar depois.
   Int_t bin_for_normInt_min_bkg = h_qinv_bkg_OS_Corr[i_centbin]->GetXaxis()->FindBin(qinv_min); //interval for normalization of qinv distributions
   Int_t bin_for_normInt_max_bkg = h_qinv_bkg_OS_Corr[i_centbin]->GetXaxis()->FindBin(qinv_max);
   
  
   Double_t int_num_controlRegion = h_qinv_sig_SS_Corr[i_centbin]->Integral(bin_for_normInt_min,bin_for_normInt_max); 
   Double_t int_den_controlRegion = h_qinv_bkg_OS_Corr[i_centbin]->Integral(bin_for_normInt_min,bin_for_normInt_max);
   
  //h_qinv_sig_SS_Corr[i_centbin]->Sumw2();
  //h_qinv_bkg_OS_Corr[i_centbin]->Sumw2();

   h_sr_SS_Over_OS[i_centbin]->Divide(h_qinv_sig_SS_Corr[i_centbin],h_qinv_bkg_OS_Corr[i_centbin],int_num_controlRegion,int_den_controlRegion); 
   
   func_hist_custom_sr(h_sr_SS_Over_OS[i_centbin]);
   h_sr_SS_Over_OS[i_centbin]->SetMarkerStyle(20); 
   h_sr_SS_Over_OS[i_centbin]->SetMarkerSize(0.65); 
   h_sr_SS_Over_OS[i_centbin]->SetMarkerColor(kBlack);
   h_sr_SS_Over_OS[i_centbin]->GetXaxis()->SetRangeUser(0.0,0.99);
	h_sr_SS_Over_OS[i_centbin]->GetYaxis()->SetRangeUser(0.75,1.405);
	h_sr_SS_Over_OS[i_centbin]->GetYaxis()->SetNdivisions(805);
	h_sr_SS_Over_OS[i_centbin]->GetXaxis()->SetNdivisions(805);
   h_sr_SS_Over_OS[i_centbin]->Write();
   // Desenhe o histograma no TCanvas.Q
  // h_sr_SS_Over_OS[i_centbin]->Draw();
   
   // Desativa a escala logarítmica no eixo y
   gPad->SetLogy(0);
   
/**********************************************************/   
   

//fit exp
	TF1 *f_exp = new TF1("f_exp",func1_exp,0.02, 0.99, 4); //try other fits ranges
	f_exp->SetParameters(1.0,1.0,4.0,0.0); //initialize parameters
	f_exp->SetParName(0,"Const");
	f_exp->SetParName(2,"R (fm)");
	f_exp->SetParName(3,"#delta");
	f_exp->SetLineColor(kRed); 
	f_exp->SetLineStyle(1);
	f_exp->SetLineWidth(2);  
	TFitResultPtr res_exp;;
   h_sr_SS_Over_OS[i_centbin]->Fit(f_exp, "S R Q N");
   res_exp= h_sr_SS_Over_OS[i_centbin]->Fit(f_exp, "S R M");
   f_exp->Draw("SameL");
	

   TF1 *f_gauss = new TF1("f_gauss",func2_gauss,0.02, 0.99, 4); //try other fits ranges
	f_gauss->SetParameters(1.0,1.0,4.0,0.0); //initialize parameters
	f_gauss->SetParName(0,"Const");
	f_gauss->SetParName(1,"#lambda"); 
	f_gauss->SetParName(2,"R (fm)");
	f_gauss->SetParName(3,"#delta");
	f_gauss->SetLineColor(SPcolors[SPdgreen]); 
	f_gauss->SetLineStyle(9);
	f_gauss->SetLineWidth(3);
	TFitResultPtr res_gauss;
   h_sr_SS_Over_OS[i_centbin]->Fit(f_gauss, "S R ");
   res_gauss = h_sr_SS_Over_OS[i_centbin]->Fit(f_gauss, "S R M");
   f_gauss->Draw("SameL");
	gPad->SetLogy(0);

   f_exp->Write();
   f_gauss->Write();
   
   //Legenda
   TLegend *legend = new TLegend();
   legend->AddEntry(h_sr_SS_Over_OS[i_centbin], "Data", "p");
   // Definindo o tamanho da fonte da legenda  
   legend->SetTextSize(0.030); 
   legend->SetBorderSize(0);
   legend->SetFillColorAlpha(0, 0);
   
   


   TLegend *lgauss = new TLegend(0.52,0.46087,0.68,0.60435);
	lgauss->SetBorderSize(0);
	lgauss->SetTextSize(0.030);
	lgauss->AddEntry(f_gauss,"Gaussian fit","l");
   lgauss->AddEntry(f_gauss, Form("R = %.4f #pm %.4f fm",f_gauss->GetParameter(2), f_gauss->GetParError(2)), "");
   lgauss->AddEntry(f_gauss, Form("#delta = %.4f #pm %.4f", f_gauss->GetParameter(3), f_gauss->GetParError(3)), "");
   lgauss->AddEntry(f_gauss, Form("#lambda = %.4f", f_gauss->GetParameter(1)), "");
	

	auto lexp = new TLegend(0.51,0.65,0.81,0.75);
	lexp->SetBorderSize(0);
	lexp->SetTextSize(0.030);
	lexp->AddEntry(f_exp,"Exponential fit","l");
	lexp->AddEntry(f_exp, Form("R = %.4f #pm %.4f fm", f_exp->GetParameter(2), f_exp->GetParError(2)), "");
   lexp->AddEntry(f_exp, Form("#delta = %.4f #pm %.4f", f_exp->GetParameter(3), f_exp->GetParError(3)), "");
   


   
    //Cria texto dentro do gráfico - lado superior esquerdo
   TLatex *text = new TLatex();
   text->SetNDC();
   text->SetTextColor(kBlack);
   text->SetTextFont(62);  // Define a fonte do texto (opcional)
   text->SetTextSize(0.025);  // Define o tamanho do texto (opcional)
   text->DrawLatex(0.15, 0.82, "CMS Open Data PbPb #sqrt{s_{NN}} = 2.76 TeV");


   TLatex *text2 = new TLatex();
   text2->SetNDC();
   text2->SetTextColor(kBlack);
   text2->SetTextFont(42);
   text2->SetTextSize(0.02);
   text2->DrawLatex(0.15, 0.78, histTitle_sr_SS_Over_OS[i_centbin]);
   
   TLatex *text_bin = new TLatex();
   text_bin->SetNDC();
   text_bin->SetTextColor(kBlack);
   text_bin->SetTextFont(42);
   text_bin->SetTextSize(0.02);
   text_bin->DrawLatex(0.15, 0.8, "Normalization bins: (0.9,1.7)");

   
    //Cria texto dentro do gráfico - lado superior esquerdo
   TLatex *text9 = new TLatex();
   text9->SetNDC();
   text9->SetTextColor(kBlack);
   text9->SetTextFont(62);  // Define a fonte do texto (opcional)
   text9->SetTextSize(0.02);  // Define o tamanho do texto (opcional)
   text9->DrawLatex(0.48, 0.91, "Method: Pairs of opposite charges(SS_{++,--}/OS_{+-})");


   //Cria uma linha reta em y = 1
   TLine *line = new TLine(c->GetUxmin(), 1, c->GetUxmax(), 1);
   line->SetLineColor(kBlack); // Define a cor da linha (opcional)
   line->SetLineStyle(2); // Define o estilo da linha (opcional)
   line->Draw();
/*************************************************************/   
  
  // Salve o TCanvas como um arquivo PDF
   

   // Atualize o TCanvas

  
   lgauss->Draw();
   lexp->Draw();
   c->SaveAs(Form("%s/%s.pdf",folder_.c_str(),histName_sr_SS_Over_OS[i_centbin].c_str()));
   c->Update();
   
   

}


}

/***********************************************************************************   
DEfinicao de dividir dois argumentos: https://root.cern.ch/doc/master/classTH1.html#a4ef2329285beb6091515b836e160bd9f
************************************************************************************/

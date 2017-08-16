#include <vector>
#include <string>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>

//TString dirname = "170712_smear_MMC_VBF_million";

void HistPlotMET ()
  //const std::string& dirname)
{
  // Load the libraries for all packages
  //gROOT->ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C");

//To make our life simpler, define the plot parameters at the beginning of our code:
// the name of the histogram we are trying to plot
  //std::string histName = "smear_MMC";


//The way our code will work is that we first create all the histograms, and then plot them. So we need some variables to hold the variables in between (you'll also need to add an include for <vector> in the beginning):

// these vectors contain the data and mc histograms
 // std::vector<TH1F*> data;
  //, mc;

//Now we use the ability of SampleHandler to loop over events to get all the histograms we need:
// create and initialize the sample handler

  //TFile *v = new TFile("/afs/cern.ch/user/g/gpan/public/170726_MMC_VBF_million/hist-VBF.root");
  TFile *g = new TFile("/afs/cern.ch/user/g/gpan/public/170726_ggH_1million_MMC2016MC15C/hist-ggh.root");
  //TFile *z = new TFile("/afs/cern.ch/user/g/gpan/public/170718_MMC_Ztautau_2million/hist-Ztautau.root");

  //f->ls();

/*
  TH1F *v_truth_MMC = (TH1F*) v->Get("MMC_truth");
       float v_nocut_norm_truth = v_truth_MMC->Integral();
       v_truth_MMC->Scale(1/v_nocut_norm_truth);
       std::cout<<"VBF entries no cut truth: "<<std::endl;
       std::cout<<v_nocut_norm_truth<<std::endl;
  TH1F *v_truth_MMC_40 = (TH1F*) v->Get("MMC_40_truth");
  TH1F *v_truth_MMC_30 = (TH1F*) v->Get("MMC_30_1_truth");
  TH1F *v_truth_MMC_eta = (TH1F*) v->Get("MMC_eta_truth");
  TH1F *v_truth_MMC_met = (TH1F*) v->Get("MMC_met_truth");
  TH1F *v_truth_MMC_dR = (TH1F*) v->Get("MMC_dR_truth");
  TH1F *v_truth_MMC_deta = (TH1F*) v->Get("MMC_deta_truth");
  TH1F *v_truth_MMC_cent = (TH1F*) v->Get("MMC_cent_truth");
        float v_full_norm_truth = v_truth_MMC_cent->Integral();
        v_truth_MMC_cent->Scale(1/v_full_norm_truth);
        std::cout<<"\nVBF entries full cuts truth: "<<std::endl;
        std::cout<<v_full_norm_truth<<std::endl;




  TH1F *v_smear_MMC = (TH1F*) v->Get("smear_MMC");
       float v_nocut_norm = v_smear_MMC->Integral();
       v_smear_MMC->Scale(1/v_nocut_norm);
       std::cout<<"VBF entries no cut: "<<std::endl;
       std::cout<<v_nocut_norm<<std::endl;
  TH1F *v_smear_MMC_40 = (TH1F*) v->Get("smear_MMC_40");
  TH1F *v_smear_MMC_30 = (TH1F*) v->Get("smear_MMC_30_1");
  TH1F *v_smear_MMC_eta = (TH1F*) v->Get("smear_MMC_eta");
  TH1F *v_smear_MMC_met = (TH1F*) v->Get("smear_MMC_met");
  TH1F *v_smear_MMC_dR = (TH1F*) v->Get("smear_MMC_dR");
  TH1F *v_smear_MMC_deta = (TH1F*) v->Get("smear_MMC_deta");
  TH1F *v_smear_MMC_cent = (TH1F*) v->Get("smear_MMC_cent");
        float v_full_norm = v_smear_MMC_cent->Integral();
        std::cout<<"\nVBF entries full cuts: "<<std::endl;
        std::cout<<v_full_norm<<std::endl;
        v_smear_MMC_cent->Scale(1/v_full_norm);

        */


  TH1F *g_truth_MMC = (TH1F*) g->Get("MMC_truth");
      /*  float g_nocut_norm_truth = g_truth_MMC->Integral();
        g_truth_MMC->Scale(1/g_nocut_norm_truth);
        std::cout<<"\n ggH entries no cut truth: "<<std::endl;
        std::cout<<g_nocut_norm_truth<<std::endl;*/
  TH1F *g_truth_MMC_40 = (TH1F*) g->Get("MMC_40_truth");
  TH1F *g_truth_MMC_30 = (TH1F*) g->Get("MMC_30_1_truth");
  TH1F *g_truth_MMC_eta = (TH1F*) g->Get("MMC_eta_truth");
  TH1F *g_truth_MMC_met = (TH1F*) g->Get("MMC_met_truth");
  TH1F *g_truth_MMC_dR = (TH1F*) g->Get("MMC_dR_truth");
  TH1F *g_truth_MMC_deta = (TH1F*) g->Get("MMC_deta_truth");
  TH1F *g_truth_MMC_cent = (TH1F*) g->Get("MMC_cent_truth");
      /*  float g_full_norm_truth = g_truth_MMC_cent->Integral();
        std::cout<<"\n ggH entries full cut truth: "<<std::endl;
        std::cout<<g_full_norm_truth<<std::endl;
        g_truth_MMC_cent->Scale(1/g_full_norm_truth);*/

  TH1F *g_smear_MMC = (TH1F*) g->Get("smear_MMC");
      /*  float g_nocut_norm = g_smear_MMC->Integral();
        g_smear_MMC->Scale(1/g_nocut_norm);
        std::cout<<"\n ggH entries no cut: "<<std::endl;
        std::cout<<g_nocut_norm<<std::endl;*/
  TH1F *g_smear_MMC_40 = (TH1F*) g->Get("smear_MMC_40");
  TH1F *g_smear_MMC_30 = (TH1F*) g->Get("smear_MMC_30_1");
  TH1F *g_smear_MMC_eta = (TH1F*) g->Get("smear_MMC_eta");
  TH1F *g_smear_MMC_met = (TH1F*) g->Get("smear_MMC_met");
  TH1F *g_smear_MMC_dR = (TH1F*) g->Get("smear_MMC_dR");
  TH1F *g_smear_MMC_deta = (TH1F*) g->Get("smear_MMC_deta");
  TH1F *g_smear_MMC_cent = (TH1F*) g->Get("smear_MMC_cent");
       /* float g_full_norm = g_smear_MMC_cent->Integral();
        std::cout<<"\n ggH entries full cut: "<<std::endl;
        std::cout<<g_full_norm<<std::endl;
        g_smear_MMC_cent->Scale(1/g_full_norm);*/


/*

  TH1F *z_truth_MMC = (TH1F*) z->Get("MMC_truth");
        float z_nocut_norm_truth = z_truth_MMC->Integral();
        z_truth_MMC->Scale(1/z_nocut_norm_truth);
        std::cout<<"\n Ztautau entries no cut truth: "<<std::endl;
        std::cout<<z_nocut_norm_truth<<std::endl;
  TH1F *z_truth_MMC_40 = (TH1F*) z->Get("MMC_40_truth");
  TH1F *z_truth_MMC_30 = (TH1F*) z->Get("MMC_30_1_truth");
  TH1F *z_truth_MMC_eta = (TH1F*) z->Get("MMC_eta_truth");
  TH1F *z_truth_MMC_met = (TH1F*) z->Get("MMC_met_truth");
  TH1F *z_truth_MMC_dR = (TH1F*) z->Get("MMC_dR_truth");
  TH1F *z_truth_MMC_deta = (TH1F*) z->Get("MMC_deta_truth");
  TH1F *z_truth_MMC_cent = (TH1F*) z->Get("MMC_cent_truth");
        float z_full_norm_truth = z_truth_MMC_cent->Integral();
        std::cout<<"\n Ztautau entries full cut truth: "<<std::endl;
        std::cout<<z_full_norm_truth<<std::endl;
        z_truth_MMC_cent->Scale(1/z_full_norm_truth);

  TH1F *z_smear_MMC = (TH1F*) z->Get("smear_MMC");
        float z_nocut_norm = z_smear_MMC->Integral();
        z_smear_MMC->Scale(1/z_nocut_norm);
        std::cout<<"\n Ztautau entries no cut: "<<std::endl;
        std::cout<<z_nocut_norm<<std::endl;
  TH1F *z_smear_MMC_40 = (TH1F*) z->Get("smear_MMC_40");
  TH1F *z_smear_MMC_30 = (TH1F*) z->Get("smear_MMC_30_1");
  TH1F *z_smear_MMC_eta = (TH1F*) z->Get("smear_MMC_eta");
  TH1F *z_smear_MMC_met = (TH1F*) z->Get("smear_MMC_met");
  TH1F *z_smear_MMC_dR = (TH1F*) z->Get("smear_MMC_dR");
  TH1F *z_smear_MMC_deta = (TH1F*) z->Get("smear_MMC_deta");
  TH1F *z_smear_MMC_cent = (TH1F*) z->Get("smear_MMC_cent");
        float z_full_norm = z_smear_MMC_cent->Integral();
        std::cout<<"\n Ztautau entries fullcut: "<<std::endl;
        std::cout<<z_full_norm<<std::endl;
        z_smear_MMC_cent->Scale(1/z_full_norm);

*/

  TCanvas* c1 = new TCanvas("c1","ggH MMC MMC2016MC15C",1500,500);
  c1->Divide(2,1);

  gStyle->SetOptStat(0);
  //c1->Divide(1,2);
  
  c1->cd(1);

/*
  v_smear_MMC->Draw("HIST"); v_smear_MMC->SetLineColor(kRed); v_smear_MMC->SetTitle("MMC with no cuts");
      v_smear_MMC->SetTitle("No cuts on MMC with smeared inputs");
  g_smear_MMC->Draw("same HIST"); g_smear_MMC->SetLineColor(kViolet);
  z_smear_MMC->Draw("same HIST"); z_smear_MMC->SetLineColor(kBlue);

   auto legend1 = new TLegend(0.65,0.7,0.90,0.9);
   legend1->AddEntry(v_smear_MMC,"VBF, no cuts");
   legend1->AddEntry(g_smear_MMC,"ggH, no cuts");
   legend1->AddEntry(z_smear_MMC,"Z#rightarrow#tau#tau, no cuts");
   legend1->Draw(); */
/*

  v_truth_MMC->Draw("HIST"); v_truth_MMC->SetLineColor(kRed); v_truth_MMC->SetTitle("MMC with no cuts");
      v_truth_MMC->GetXaxis()->SetRange(10,200);
      v_truth_MMC->SetTitle("No cuts on MMC with truth inputs");
  g_truth_MMC->Draw("same HIST"); g_truth_MMC->SetLineColor(kViolet);
      g_truth_MMC->GetXaxis()->SetRange(10,200); 
  z_truth_MMC->Draw("same HIST"); z_truth_MMC->SetLineColor(kBlue);
      z_truth_MMC->GetXaxis()->SetRange(10,200);

      

   auto legend1 = new TLegend(0.65,0.7,0.90,0.9);
   legend1->AddEntry(v_truth_MMC,"VBF, no cuts");
   legend1->AddEntry(g_truth_MMC,"ggH, no cuts");
   legend1->AddEntry(z_truth_MMC,"Z#rightarrow#tau#tau, no cuts");
   legend1->Draw();

  v_truth_MMC->Draw(); v_truth_MMC->GetXaxis()->SetRange(10,200); v_truth_MMC->SetTitle("VBF MMC with truth inputs");
    v_truth_MMC->GetYaxis()->SetTitle("counts");
  v_truth_MMC_40->Draw("same"); v_truth_MMC_40->GetXaxis()->SetRange(10,200);
  v_truth_MMC_30->Draw("same"); v_truth_MMC_30->GetXaxis()->SetRange(10,200);
  v_truth_MMC_eta->Draw("same"); v_truth_MMC_eta->GetXaxis()->SetRange(10,200);
  v_truth_MMC_met->Draw("same"); v_truth_MMC_met->GetXaxis()->SetRange(10,200);
  v_truth_MMC_dR->Draw("same"); v_truth_MMC_dR->GetXaxis()->SetRange(10,200);
  v_truth_MMC_deta->Draw("same"); v_truth_MMC_deta->GetXaxis()->SetRange(10,200);
  v_truth_MMC_cent->Draw("same"); v_truth_MMC_cent->GetXaxis()->SetRange(10,200);

   auto legendt = new TLegend(0.65,0.7,0.90,0.9);
   legendt->AddEntry(v_truth_MMC,"ME_{T} > 20 GeV and #eta < 2.5");
   legendt->AddEntry(v_truth_MMC_40,"+ #tau_{0} p_{T} > 40 GeV");
   legendt->AddEntry(v_truth_MMC_30,"+ #tau_{1} p_{T} > 30 GeV");
   legendt->AddEntry(v_truth_MMC_eta,"+ #eta veto");
   legendt->AddEntry(v_truth_MMC_met,"+ ME_{T} > 20 GeV");
   legendt->AddEntry(v_truth_MMC_dR,"+ 0.8 < dR < 2.4");
   legendt->AddEntry(v_truth_MMC_deta,"+ d#eta < 1.5");
   legendt->AddEntry(v_truth_MMC_cent,"+ ME_{T} centrality");
   legendt->Draw();

  v_smear_MMC->Draw("HIST"); v_smear_MMC->GetXaxis()->SetRange(10,200);
  v_smear_MMC_40->Draw("same HIST"); v_smear_MMC_40->GetXaxis()->SetRange(10,200);
  v_smear_MMC_30->Draw("same HIST"); v_smear_MMC_30->GetXaxis()->SetRange(10,200);
  v_smear_MMC_eta->Draw("same HIST"); v_smear_MMC_eta->GetXaxis()->SetRange(10,200);
  v_smear_MMC_met->Draw("same HIST"); v_smear_MMC_met->GetXaxis()->SetRange(10,200);
  v_smear_MMC_dR->Draw("same HIST"); v_smear_MMC_dR->GetXaxis()->SetRange(10,200);
  v_smear_MMC_deta->Draw("same HIST"); v_smear_MMC_deta->GetXaxis()->SetRange(10,200);
  v_smear_MMC_cent->Draw("same HIST"); v_smear_MMC_cent->GetXaxis()->SetRange(10,200);

   auto legend = new TLegend(0.65,0.7,0.90,0.9);
   legend->AddEntry(v_smear_MMC,"ME_{T} > 20 GeV and #eta < 2.5");
   legend->AddEntry(v_smear_MMC_40,"+ #tau_{0} p_{T} > 40 GeV");
   legend->AddEntry(v_smear_MMC_30,"+ #tau_{1} p_{T} > 30 GeV");
   legend->AddEntry(v_smear_MMC_eta,"+ #eta veto");
   legend->AddEntry(v_smear_MMC_met,"+ ME_{T} > 20 GeV");
   legend->AddEntry(v_smear_MMC_dR,"+ 0.8 < dR < 2.4");
   legend->AddEntry(v_smear_MMC_deta,"+ d#eta < 1.5");
   legend->AddEntry(v_smear_MMC_cent,"+ ME_{T} centrality");
   legend->Draw(); 

  */

  g_truth_MMC->Draw(); g_truth_MMC->GetXaxis()->SetRange(10,200); g_truth_MMC->SetTitle("ggH MMC with truth inputs, MMC2016MC15C");
    g_truth_MMC->GetYaxis()->SetTitle("counts");
  g_truth_MMC_40->Draw("same"); g_truth_MMC_40->GetXaxis()->SetRange(10,200);
  g_truth_MMC_30->Draw("same"); g_truth_MMC_30->GetXaxis()->SetRange(10,200);
  g_truth_MMC_eta->Draw("same"); g_truth_MMC_eta->GetXaxis()->SetRange(10,200);
  g_truth_MMC_met->Draw("same"); g_truth_MMC_met->GetXaxis()->SetRange(10,200);
  g_truth_MMC_dR->Draw("same"); g_truth_MMC_dR->GetXaxis()->SetRange(10,200);
  g_truth_MMC_deta->Draw("same"); g_truth_MMC_deta->GetXaxis()->SetRange(10,200);
  g_truth_MMC_cent->Draw("same"); g_truth_MMC_cent->GetXaxis()->SetRange(10,200);

   auto legendt = new TLegend(0.65,0.7,0.90,0.9);
   legendt->AddEntry(g_truth_MMC,"ME_{T} > 20 GeV and #eta < 2.5");
   legendt->AddEntry(g_truth_MMC_40,"+ #tau_{0} p_{T} > 40 GeV");
   legendt->AddEntry(g_truth_MMC_30,"+ #tau_{1} p_{T} > 30 GeV");
   legendt->AddEntry(g_truth_MMC_eta,"+ #eta veto");
   legendt->AddEntry(g_truth_MMC_met,"+ ME_{T} > 20 GeV");
   legendt->AddEntry(g_truth_MMC_dR,"+ 0.8 < dR < 2.4");
   legendt->AddEntry(g_truth_MMC_deta,"+ d#eta < 1.5");
   legendt->AddEntry(g_truth_MMC_cent,"+ ME_{T} centrality");
   legendt->Draw();

  c1->cd(2);

  g_smear_MMC->Draw("HIST"); g_smear_MMC->GetXaxis()->SetRange(10,200); g_smear_MMC->SetTitle("ggH MMC with smeared inputs, MMC2016MC15C");
  g_smear_MMC_40->Draw("same HIST"); g_smear_MMC_40->GetXaxis()->SetRange(10,200);
  g_smear_MMC_30->Draw("same HIST"); g_smear_MMC_30->GetXaxis()->SetRange(10,200);
  g_smear_MMC_eta->Draw("same HIST"); g_smear_MMC_eta->GetXaxis()->SetRange(10,200);
  g_smear_MMC_met->Draw("same HIST"); g_smear_MMC_met->GetXaxis()->SetRange(10,200);
  g_smear_MMC_dR->Draw("same HIST"); g_smear_MMC_dR->GetXaxis()->SetRange(10,200);
  g_smear_MMC_deta->Draw("same HIST"); g_smear_MMC_deta->GetXaxis()->SetRange(10,200);
  g_smear_MMC_cent->Draw("same HIST"); g_smear_MMC_cent->GetXaxis()->SetRange(10,200);

   auto legend = new TLegend(0.65,0.7,0.90,0.9);
   legend->AddEntry(g_smear_MMC,"ME_{T} > 20 GeV and #eta < 2.5");
   legend->AddEntry(g_smear_MMC_40,"+ #tau_{0} p_{T} > 40 GeV");
   legend->AddEntry(g_smear_MMC_30,"+ #tau_{1} p_{T} > 30 GeV");
   legend->AddEntry(g_smear_MMC_eta,"+ #eta veto");
   legend->AddEntry(g_smear_MMC_met,"+ ME_{T} > 20 GeV");
   legend->AddEntry(g_smear_MMC_dR,"+ 0.8 < dR < 2.4");
   legend->AddEntry(g_smear_MMC_deta,"+ d#eta < 1.5");
   legend->AddEntry(g_smear_MMC_cent,"+ ME_{T} centrality");
   legend->Draw();

   /*

  z_truth_MMC->Draw(); z_truth_MMC->GetXaxis()->SetRange(10,200); z_truth_MMC->SetTitle("Z#rightarrow#tau#tau MMC with truth inputs");
  z_truth_MMC->GetYaxis()->SetTitle("counts");
  z_truth_MMC_40->Draw("same"); z_truth_MMC_40->GetXaxis()->SetRange(10,200);
  z_truth_MMC_30->Draw("same"); z_truth_MMC_30->GetXaxis()->SetRange(10,200);
  z_truth_MMC_eta->Draw("same"); z_truth_MMC_eta->GetXaxis()->SetRange(10,200);
  z_truth_MMC_met->Draw("same"); z_truth_MMC_met->GetXaxis()->SetRange(10,200);
  z_truth_MMC_dR->Draw("same"); z_truth_MMC_dR->GetXaxis()->SetRange(10,200);
  z_truth_MMC_deta->Draw("same"); z_truth_MMC_deta->GetXaxis()->SetRange(10,200);
  z_truth_MMC_cent->Draw("same"); z_truth_MMC_cent->GetXaxis()->SetRange(10,200);

   auto legendt = new TLegend(0.65,0.7,0.90,0.9);
   legendt->AddEntry(z_truth_MMC,"ME_{T} > 20 GeV and #eta < 2.5");
   legendt->AddEntry(z_truth_MMC_40,"+ #tau_{0} p_{T} > 40 GeV");
   legendt->AddEntry(z_truth_MMC_30,"+ #tau_{1} p_{T} > 30 GeV");
   legendt->AddEntry(z_truth_MMC_eta,"+ #eta veto");
   legendt->AddEntry(z_truth_MMC_met,"+ ME_{T} > 20 GeV");
   legendt->AddEntry(z_truth_MMC_dR,"+ 0.8 < dR < 2.4");
   legendt->AddEntry(z_truth_MMC_deta,"+ d#eta < 1.5");
   legendt->AddEntry(z_truth_MMC_cent,"+ ME_{T} centrality");
   legendt->Draw();
  
  z_smear_MMC->Draw("HIST"); z_smear_MMC->GetXaxis()->SetRange(10,200);
  z_smear_MMC_40->Draw("same HIST"); z_smear_MMC_40->GetXaxis()->SetRange(10,200);
  z_smear_MMC_30->Draw("same HIST"); z_smear_MMC_30->GetXaxis()->SetRange(10,200);
  z_smear_MMC_eta->Draw("same HIST"); z_smear_MMC_eta->GetXaxis()->SetRange(10,200);
  z_smear_MMC_met->Draw("same HIST"); z_smear_MMC_met->GetXaxis()->SetRange(10,200);
  z_smear_MMC_dR->Draw("same HIST"); z_smear_MMC_dR->GetXaxis()->SetRange(10,200);
  z_smear_MMC_deta->Draw("same HIST"); z_smear_MMC_deta->GetXaxis()->SetRange(10,200);
  z_smear_MMC_cent->Draw("same HIST"); z_smear_MMC_cent->GetXaxis()->SetRange(10,200);

   auto legend = new TLegend(0.65,0.7,0.90,0.9);
   legend->AddEntry(z_smear_MMC,"ME_{T} > 20 GeV and #eta < 2.5");
   legend->AddEntry(z_smear_MMC_40,"+ #tau_{0} p_{T} > 40 GeV");
   legend->AddEntry(z_smear_MMC_30,"+ #tau_{1} p_{T} > 30 GeV");
   legend->AddEntry(z_smear_MMC_eta,"+ #eta veto");
   legend->AddEntry(z_smear_MMC_met,"+ ME_{T} > 20 GeV");
   legend->AddEntry(z_smear_MMC_dR,"+ 0.8 < dR < 2.4");
   legend->AddEntry(z_smear_MMC_deta,"+ d#eta < 1.5");
   legend->AddEntry(z_smear_MMC_cent,"+ ME_{T} centrality");
   legend->Draw();*/


  /*
  v_smear_MMC_cent->Draw("HIST"); v_smear_MMC_cent->SetLineColor(kRed); v_smear_MMC_cent->GetXaxis()->SetTitle("MMC Mass (GeV)");
    v_smear_MMC_cent->SetTitle("Full cuts on MMC with smeared inputs");
  g_smear_MMC_cent->Draw("same HIST"); g_smear_MMC_cent->SetLineColor(kViolet);
  z_smear_MMC_cent->Draw("same HIST"); z_smear_MMC_cent->SetLineColor(kBlue);

   auto legend = new TLegend(0.65,0.7,0.90,0.9);
   legend->AddEntry(v_smear_MMC_cent,"VBF, all cuts");
   legend->AddEntry(g_smear_MMC_cent,"ggH, all cuts");
   legend->AddEntry(z_smear_MMC_cent,"Z#rightarrow#tau#tau, all cuts");
   legend->Draw(); */

/*
  v_truth_MMC_cent->Draw("HIST"); v_truth_MMC_cent->SetLineColor(kRed); v_truth_MMC_cent->GetXaxis()->SetTitle("MMC Mass (GeV)");
    v_truth_MMC_cent->SetTitle("Full cuts on MMC with truth inputs");
    v_truth_MMC_cent->GetXaxis()->SetRange(10,200);
  g_truth_MMC_cent->Draw("same HIST"); g_truth_MMC_cent->SetLineColor(kViolet);
    g_truth_MMC_cent->GetXaxis()->SetRange(10,200);
  z_truth_MMC_cent->Draw("same HIST"); z_truth_MMC_cent->SetLineColor(kBlue);
    z_truth_MMC_cent->GetXaxis()->SetRange(10,200);

   auto legend = new TLegend(0.65,0.7,0.90,0.9);
   legend->AddEntry(v_truth_MMC_cent,"VBF, all cuts");
   legend->AddEntry(g_truth_MMC_cent,"ggH, all cuts");
   legend->AddEntry(z_truth_MMC_cent,"Z#rightarrow#tau#tau, all cuts");
   legend->Draw();

   */

  //SH::SampleHandler sh;
  //sh.load (dirname + "/hist");
  //sh.load ("/afs/cern.ch/user/g/gpan/public/170712_smear_MMC_VBF_million/hist");
  // add any extra meta-data you desire here

  //for (unsigned samplesIter = 0; samplesIter != sh.size(); ++ samplesIter)
  //{
    // the sample we are plotting
   // SH::Sample *sample = sh.at(0);

    // determine whether this is data or MC
    // bool isMC = sample->name().find ("data") == std::string::npos;

    // retrieve the histogram
    //TH1F *hist = (TH1F*) sample->readHist(histName);

    // store the histogram for later.
    /*if (isMC){
      mc.push_back (hist);
    }*/
    //else{
      //data.push_back (hist);
    //}
  //} // end for loop over samples sh.size

//Since we want to make nice data-MC comparison stack plots we now have to do a bit of addition to get the correct plots:
// add together the histograms for proper stack-plots
  /*for (unsigned iter = 1; iter < data.size(); ++ iter){
    data[0]->Add (data[iter]);
  }
  */
 /* for (unsigned iter = 1; iter < mc.size(); ++ iter){
    mc[mc.size()-1-iter]->Add (mc[mc.size()-iter]);
  }*/

//And now just plot the histograms in the right order in a decent style:
// now plot all mc histograms
  /*for (unsigned iter = 0; iter < mc.size(); ++ iter)
  {
    mc[iter]->SetLineColor (2 + iter);
    mc[iter]->SetFillColor (2 + iter);
    mc[iter]->SetFillStyle (3001);
    mc[iter]->SetTitle (0);
    mc[iter]->Draw (iter ? "SAME HIST" : "HIST");
  }*/

  // and if we have data, plot that too
  /*if (!data.empty())
  {
    data[0]->SetMarkerStyle(20);
    data[0]->SetTitle (0);
    //data[0]->Draw (mc.empty() ? "ERR" : "SAME ERR");
  }*/



   ///////////////////////////////Carolyn's Plots////////////////////////////////////

   h_tau0ptMMC = new TH1F("MMC_no_cuts","MMC_no_cuts",200,0,200);
  h_tau0ptcut0MMC = new TH1F("MMC_cut_0pt","MMC_cut_0pt",200,0,200);
  h_tau0ptcut0cut1MMC = new TH1F("MMC_cut_1pt", "MMC_cut_1pt",200,0,200);
  h_tau0ptcutetaMMC = new TH1F("MMC_cut_eta","MMC_cut_eta",200,0,200);   
  h_tau0ptmetMMC = new TH1F("MMC_cut_met","MMC_cut_met",200,0,200);
  h_tau0ptdrMMC = new TH1F("MMC_cut_dr","MMC_cut_dr",200,0,200); 
  h_tau0ptdetaMMC = new TH1F("MMC_cut_deta","MMC_cut_deta",200,0,200); 
  h_tau0ptetmissMMC = new TH1F("MMC_met_centrality","MMC_met_centrality",200,0,200); 

  h_tau0ptMMCtruth = new TH1F("MMC_no_cuts_truth","MMC_no_cuts_truth",200,0,200);
  h_tau0ptcut0MMCtruth = new TH1F("MMC_cut_0pt)truth","MMC_cut_0pt_truth",200,0,200);
  h_tau0ptcut0cut1MMCtruth = new TH1F("MMC_cut_1pt_truth", "MMC_cut_1pt_truth",200,0,200);
  h_tau0ptcutetaMMCtruth = new TH1F("MMC_cut_eta_truth","MMC_cut_eta_truth",200,0,200);   
  h_tau0ptmetMMCtruth = new TH1F("MMC_cut_met_truth","MMC_cut_met_truth",200,0,200);
  h_tau0ptdrMMCtruth = new TH1F("MMC_cut_dr_truth","MMC_cut_dr_truth",200,0,200); 
  h_tau0ptdetaMMCtruth = new TH1F("MMC_cut_deta_truth","MMC_cut_deta_truth",200,0,200); 
  h_tau0ptetmissMMCtruth = new TH1F("MMC_met_centrality_truth","MMC_met_centrality_truth",200,0,200); 

  h_beforet = new TH1F("no_cuts","no_cuts",140,0,140);
  h_aftert = new TH1F("all_cuts","all_cuts",140,0,140);
  h_beforem = new TH1F("no_cuts_smeared", "no_cuts_smeared",140,0,140);
  h_afterm = new TH1F("all_cuts_smeared","all_cuts_smeared",140,0,140);
}
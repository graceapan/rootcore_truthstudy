#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <DiTauMassTools/MissingMassCalculator.h>
#include "DiTauMassTools/MissingMassTool.h"
#include "DiTauMassTools/DiTauMassToolsDict.h"   
#include "DiTauMassTools/IMissingMassTool.h" 
#include <truth_study/truth_study.h>

#include <string>


#include <TStyle.h>
#include <TRandom3.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TVector2.h>
#include <TH2.h>
#include <TLatex.h>

/*
#ifdef ROOTCORE
// Framework include(s):
#include "PathResolver/PathResolver.h"
#endif // ROOTCORE
*/


// this is needed to distribute the algorithm to the workers
ClassImp(truth_study)


TRandom3 m_tauRandom(7);
TRandom m_METRandom(7);
TString m_layout = "gold";
TString m_layout2 = "silver";
double PI = 3.14159265358979;
float nominal = 200;
float m_avgMu = 199.0; //pileup
//TH1F  *m_SumEtH[4][6];
//TString PUfile = "/afs/cern.ch/user/g/gpan/public/RootCore_truthstudy/packages/sumetPU_mu200_ttbar_gold.root";






double DeltaPhi(const double phi1, const double phi2){
    return PI -  std::fabs(std::fabs(phi1 - phi2) - PI );
}

bool METinDitauDPhi(double METphi, TLorentzVector tau1, TLorentzVector tau2 ){

    double dphi  = DeltaPhi(tau1.Phi(), tau2.Phi()) ;
    double dphi1 = DeltaPhi( tau1.Phi(), METphi);
    double dphi2 = DeltaPhi( tau2.Phi(), METphi);

    if ((std::max(dphi1, dphi2) <= dphi) && (dphi1+dphi2 <= PI)){
      return true;
    }
    else if(dphi1 < PI/4 || dphi2 < PI/4){
      return true;
    }
    return false;
 }

//smears the TauEfficiency to represent how good we are at correctly identifying the 
//types of taus AFTER we have triggered on them
float getTauEfficiency(float etMeV, float eta, short prong) {
  // prong variable must be 1 or 3
  float Reco_Pt_cut = 20000;
  float Reco_eta_cut = 2.5;
  if (etMeV<Reco_Pt_cut || abs(eta)>Reco_eta_cut){
    //std::cout << "WARNING: No Reco SF for these taus - either pt or eta are out of the acceptance !" << std::endl;
    return 0;
  }
  if (prong==1)
    return 0.55;
  else if (prong==3) 
    return 0.50;
  else {
    //std::cout << "getTauEfficiency: no efficiency available for prong=" << prong << std::endl;
      return 0;
  }
}


//eMeV --> pt_vis, eta-->eta_vis, and short prong--> n_charged_pion
//used in getTauSmearedEnergy
float getTauEnergyResolution(float eMeV, float eta, short prong) {
  float stocB, stocE;
  // grab parameters for mu=40 and mu=140 points 
  // then interpolate mu between 40 and 140. 
  // do not extrapolate below 40 since we have no studies there
  // linear extrapolate above 140, but only up to 200
  float mu1 = 40.0;
  float mu2 = 140.0;
  float mumin = 40.0;
  float mumax = 200.0;
  //float m_avgMu = 199.0; //pileup
  float mu = m_avgMu;

  if (mu<mumin){
    //std::cout << "WARNING: mu outside range: " << mu << " (mumin: " << mumin << "), choosing " << mumin << std::endl;
    mu = mumin;
  }
  if (mu>mumax){
    //std::cout << "WARNING: mu outside range: " << mu << " (mumax: " << mumax << "), choosing " << mumax << std::endl;
    mu = mumax;
  }
  float stocB1, stocB2;
  float stocE1, stocE2;
  if (prong == 1){
    stocB1 = 0.96;
    stocB2 = 1.31;
    stocE1 = 0.04;
    stocE2 = 0.02;
  }
  else if(prong == 3){
    stocB1 = 1.08;
    stocB2 = 1.35;
    stocE1 = 0.04;
    stocE2 = 0.04;
  }
  else{
    //  std::cout << "getTauEnergyResolution: no energy resolution available for prong=" << prong << " eta=" << eta << std::endl;
    return 0;
  }
  
  //interpolate parameters
  stocB = stocB1 + (mu-mu1)/(mu2-mu1) * (stocB2-stocB1);
  stocE = stocE1 + (mu-mu1)/(mu2-mu1) * (stocE2-stocE1);
  
  double relativeUncertainty = hypot(stocB/sqrt(eMeV/(1000.)), stocE);
  //  return initE * (m_random3.Gaus(1, (stocB/sqrt(initE/1000.)))) * (m_random3.Gaus(1, stocE)) ;
  return eMeV*relativeUncertainty;
}


//smears visible transverse momentum (energy units) 
float getTauSmearedEnergy(float eMeV, float eta, short prong) {
  //int blah = plusminus.Integer(10000)%2;
  float correction = eMeV+m_tauRandom.Gaus(0.0,getTauEnergyResolution(eMeV,eta,prong));
  // cout<<"\n original: ";
  // cout<<eMeV;
  // cout<<"corrected: ";
  // cout<<correction;
  return correction;
 }


//represents how good we are at identifying/triggering on taus to begin with
//before narrowing down the taus to where they came from, what they're like, etc
//gold, silver, bronze refer to tracker performance. bronze is most conservative, gold is aggressive.
float getDiTauTriggerEfficiency(float et1MeV, float et2MeV, float eta1, float eta2, short prong1, short prong2) {

  float et1 = et1MeV;
  float et2 = et2MeV;
  if ( et1 < et2 ) {
    et1 = et2MeV;
    et2 = et1MeV;
  }
  float minPt1 = 40000.;
  float minPt2 = 30000.;
  float eff = 0.80; //single tau efficiency
  if ( m_layout2 == "silver" ){
       //|| m_layout == bronze) {
    minPt1 = 40000.;
    minPt2 = 30000.;
  }
  // if ( m_layout == bronze ) eff = 0.75;
  if ( fabs(eta1)>2.5 || fabs(eta2)>2.5 || prong1>3 || prong2>3 ) return 0.;
  if ( et1 > minPt1 && et2 > minPt2 ) 
    return eff*eff;
  return 0;
}
 


///////////////////////////////End of Smearing Functions///////////////////////////////////////////

  TString filename ="/afs/cern.ch/user/g/gpan/public/RootCore_truthstudy/packages/sumetPU_mu200_ttbar_gold.root";

  //std::cout << "Loading Missing ET histogram file " << filename << std::endl;

  //full path name of histogram
  std::string METFile = filename.Data();

  /*
  #ifdef ROOTCORE
  // Get file from data path
  METFile = PathResolverFindCalibFile(METFile);
  std::cout << "Found Missing ET histogram file: " << METFile << std::endl;
  #endif // ROOTCORE
  */

  TFile *m_infile = new TFile(METFile.c_str(),"READ");

  TH1F* m_SumEtH = static_cast<TH1F*>(m_infile->Get("h_sumetPU"));


 float getMETResolution(float sumEtMeV, float systValue) {
   // The sumEtMeV input should include any pileup contribution
   // This is done correctly in getMETSmeared (see above)
   double METPUreso=0;
   // double METPUresoNom=0;
   // double METPUresoVar=0;
   // double METPUresoTh=0;
   // double METPUDelta=0;
   float sumEtGeV = sumEtMeV/1000.;

   int PUcondition;
    // typedef enum {mu60=0, mu80=1, mu140=2, mu200=3} mu;
   if (m_avgMu <= 60.) PUcondition = 0;
   else if (m_avgMu <= 80.) PUcondition = 1;
   else if (m_avgMu <= 140.) PUcondition = 2;
   else PUcondition = 3;

   if (PUcondition!=3 || systValue !=nominal){
     std::cout<<"getMETResolution error: The only parametrization approved is for the nominal working point <mu>=200!"<<std::endl;
     std::cout<<"       The tool will return 0"<<std::endl;
     return 0.;
   }


   double fitParsGold[4] = {36.133, 0.00125496, 3.13886e-06, -3.79944e-10};
   double fitParsSilver[4] = {37.5551, -0.000627389, 6.11304e-06, -6.52893e-10};
   double fitParsBronze[4] = {39.049, -0.00202729, 8.30795e-06, -9.24871e-10};

   if (m_layout=="gold") {
     METPUreso = ((fitParsGold[3]*sumEtGeV + fitParsGold[2])*sumEtGeV + fitParsGold[1])*sumEtGeV + fitParsGold[0];
   } else if (m_layout=="silver") {
     METPUreso = ((fitParsSilver[3]*sumEtGeV + fitParsSilver[2])*sumEtGeV + fitParsSilver[1])*sumEtGeV + fitParsSilver[0];
   } else if (m_layout=="bronze") {
     METPUreso = ((fitParsBronze[3]*sumEtGeV + fitParsBronze[2])*sumEtGeV + fitParsBronze[1])*sumEtGeV + fitParsBronze[0];
   } else {
     std::cout << "getMETResolution error: layout " << m_layout << " is not implemented.  Returning 0." << std::endl;
     METPUreso = 0.;
   }

/*
  if (systValue==resoDown ){
    METPUreso=METPUresoNom-METPUDelta;
    if (METPUreso<0) METPUreso=0;//to be sure we do not have negative smearings.
  }
  if (systValue==resoUp){
    METPUreso=METPUresoNom+METPUDelta;
  }
  */

   return METPUreso*1000.; // conversion to MeV
 }


TVector2 getMETSmeared(float sumEtMeV, float METxMeV, float METyMeV, float systValue) {



   double sumETPU;

   //pile up ET plus truth level ET
   double sumET;
   double METPUreso;
  
   double METPUsmearX;
   double METPUsmearY;

   double first;
   double second;

   TVector2 tmpMET;

   int PUcondition;
   //  typedef enum {mu60=0, mu80=1, mu140=2, mu200=3} mu;
   if (m_avgMu <= 60.) PUcondition = 0;
   else if (m_avgMu <= 80.) PUcondition = 1;
   else if (m_avgMu <= 140.) PUcondition = 2;
   else PUcondition = 3;

   if (PUcondition!=3 || systValue != nominal){
     std::cout<<"getMETSmeared error: The only parametrization approved is for the nominal working point <mu>=200!"<<std::endl;
     std::cout<<"       The tool will return the truth MET"<<std::endl;
     tmpMET.SetX(METxMeV);
     tmpMET.SetY(METyMeV);
     return tmpMET;
   }
 
   sumETPU = m_SumEtH->GetRandom();

   //pileup ET plus truth level ET
   //sumEtMeV is truth level sum ET
   sumET = (sumETPU) + sumEtMeV; // in MeV
  
//   //this smears sumET
   METPUreso = getMETResolution(sumET, nominal);

   //then get gaussian 
   METPUsmearX = m_METRandom.Gaus(0,METPUreso);
   METPUsmearY = m_METRandom.Gaus(0,METPUreso);

  first = METxMeV + METPUsmearX;
  second = METyMeV + METPUsmearY;

   tmpMET.SetX(first);
   tmpMET.SetY(second);

   return tmpMET;
 }


///////////////////////////////End of MET Smearing Functions///////////////////////////////////////////



truth_study :: truth_study ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode truth_study :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  //  MissingMassCalculator fMMC;
  //  fMMC.SetCalibrationSet(MMCCalibrationSet::MMC2016MC15C);



  return EL::StatusCode::SUCCESS;
}



EL::StatusCode truth_study :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are

  invarm = new TH1F("invarm","#tau#tau Invariant Mass, VBF",100, 0, 100);
  wk()->addOutput (invarm);

  invarm_40 = new TH1F("invarm_40","#tau_{0} p_{T}>40 GeV, VBF",100,0, 100);
   wk()->addOutput (invarm_40);
  invarm_30 = new TH1F("invarm_30_1","+#tau_{1}>30 GeV, VBF",100, 0, 100);
   wk()->addOutput (invarm_30);
  invarm_eta = new TH1F("invarm_eta","+#eta veto, VBF", 100,0,100);
   wk()->addOutput (invarm_eta);
  invarm_met = new TH1F("invarm_met","+ME_{T}>20 GeV VBF",100,0,100);
   wk()->addOutput (invarm_met);
  invarm_dR = new TH1F("invarm_dR","+dRcut",100,0,100);
   wk()->addOutput (invarm_dR);
  invarm_deta = new TH1F("invarm_deta","+d#eta",100,0,100);
   wk()->addOutput (invarm_deta);
  invarm_cent = new TH1F("invarm_cent","+E_{miss} cent",100,0,100);
   wk()->addOutput (invarm_cent);

  smear_invarm = new TH1F("smear_invarm","#tau#tau Invariant Mass: #tau eff, pT, and trigger smears, VBF",100, 0, 100);
     wk()->addOutput (smear_invarm);
  smear_invarm_40 = new TH1F("smear_invarm_40","#tau_{0} p_{T}>40 GeV, VBF",100,0, 100);
      wk()->addOutput (smear_invarm_40);
  smear_invarm_30 = new TH1F("smear_invarm_30_1","+#tau_{1}>30 GeV, VBF",100, 0, 100);
      wk()->addOutput (smear_invarm_30);
  smear_invarm_eta = new TH1F("smear_invarm_eta","+#eta veto, VBF", 100,0,100);
      wk()->addOutput (smear_invarm_eta);
  smear_invarm_met = new TH1F("smear_invarm_met","+ME_{T}>20 GeV VBF",100,0,100);
      wk()->addOutput (smear_invarm_met);
  smear_invarm_dR = new TH1F("smear_invarm_dR","+dRcut",100,0,100);
      wk()->addOutput (smear_invarm_dR);
  smear_invarm_deta = new TH1F("smear_invarm_deta","+d#eta",100,0,100);
      wk()->addOutput (smear_invarm_deta);
  smear_invarm_cent = new TH1F("smear_invarm_cent","+E_{miss} cent",100,0,100);
      wk()->addOutput (smear_invarm_cent);


//////////////////////////////////// MMC with TRUTH PARAMETERS ////////////////////////////////////
//////////////////////////////////// MMC with TRUTH PARAMETERS  ////////////////////////////////////

      

  MMC_truth = new TH1F("MMC_truth","MMC truth no cuts",200, 0, 200);
   wk()->addOutput (MMC_truth);
  MMC_40_truth = new TH1F("MMC_40_truth","MMC truth #tau_{0} p_{T}>40 GeV",200,0, 200);
   wk()->addOutput (MMC_40_truth);
  MMC_30_truth = new TH1F("MMC_30_1_truth","MMC truth +#tau_{1}>30 GeV",200, 0, 200);
   wk()->addOutput (MMC_30_truth);
  MMC_eta_truth = new TH1F("MMC_eta_truth","MMC truth +#eta veto", 200,0,200);
   wk()->addOutput (MMC_eta_truth);
  MMC_met_truth = new TH1F("MMC_met_truth","MMC truth +ME_{T}>20 GeV ",200,0,200);
   wk()->addOutput (MMC_met_truth);
  MMC_dR_truth = new TH1F("MMC_dR_truth","MMC truth +dRcut",200,0,200);
   wk()->addOutput (MMC_dR_truth);
  MMC_deta_truth = new TH1F("MMC_deta_truth","MMC truth +d#eta",200,0,200);
   wk()->addOutput (MMC_deta_truth);
  MMC_cent_truth = new TH1F("MMC_cent_truth","MMC truth +E_{miss} cent",200,0,200);
   wk()->addOutput (MMC_cent_truth);

   


//////////////////////////////////// MMC with SMEARED PARAMETERS  ////////////////////////////////////
//////////////////////////////////// MMC with SMEARED PARAMETERS  ////////////////////////////////////


  smear_MMC = new TH1F("smear_MMC","Z#rightarrow#tau#tau MMC with smeared inputs",200, 0, 200);
     wk()->addOutput (smear_MMC);
  smear_MMC_40 = new TH1F("smear_MMC_40","MMC #tau_{0} p_{T}>40 GeV",200,0, 200);
      wk()->addOutput (smear_MMC_40);
  smear_MMC_30 = new TH1F("smear_MMC_30_1","MMC +#tau_{1}>30 GeV",200, 0, 200);
      wk()->addOutput (smear_MMC_30);
  smear_MMC_eta = new TH1F("smear_MMC_eta","MMC +#eta veto", 200,0,200);
      wk()->addOutput (smear_MMC_eta);
  smear_MMC_met = new TH1F("smear_MMC_met","MMC +ME_{T}>20 GeV ",200,0,200);
      wk()->addOutput (smear_MMC_met);
  smear_MMC_dR = new TH1F("smear_MMC_dR","MMC +dRcut",200,0,200);
      wk()->addOutput (smear_MMC_dR);
  smear_MMC_deta = new TH1F("smear_MMC_deta","MMC +d#eta",200,0,200);
      wk()->addOutput (smear_MMC_deta);
  smear_MMC_cent = new TH1F("smear_MMC_cent","MMC with all smeared cuts",200,0,200);
      wk()->addOutput (smear_MMC_cent);

//////////////////////////////////// SMEARED MET  ////////////////////////////////////
//////////////////////////////////// SMEARED MET  ////////////////////////////////////

  met_hist = new TH1F("met_hist","ME_{T}, unsmeared",100,-100,100);
      wk()->addOutput (met_hist);
  smear_met_hist = new TH1F("smear_met_hist","ME_{T}, smeared",100,-100,100);
      wk()->addOutput (smear_met_hist);
  metx_hist = new TH1F("metx_hist","ME_{Tx}, unsmeared",100,-100,100);
      wk()->addOutput (metx_hist);
  smear_metx_hist = new TH1F("smear_metx_hist","ME_{Tx}, smeared",100,-100,100);
      wk()->addOutput (smear_metx_hist);
  mety_hist = new TH1F("mety_hist","ME_{Ty} et, unsmeared",100,-100,100);
      wk()->addOutput (mety_hist);
  smear_mety_hist = new TH1F("smear_mety_hist","ME_{Ty}, smeared",100,-100,100);
      wk()->addOutput (smear_mety_hist);

  met_diff_hist = new TH1F("met_diff_hist","ME_{T}, truth - smeared",200,-200,200);
      wk()->addOutput (met_diff_hist);
  met_diffx_hist = new TH1F("met_diffx_hist","ME_{T}, truth - smeared (X)",200,-200,200);
      wk()->addOutput (met_diffx_hist);
  met_diffy_hist = new TH1F("met_diffy_hist","ME_{T}, truth - smeared (Y)",200,-200,200);
      wk()->addOutput (met_diffy_hist);


/*
  MMC_values_unsmear_cut = new TH1F("MMC_values_unsmear_cut","MMC values unsmeared cuts",150,0,150);
      wk()->addOutput (MMC_values_unsmear_cut);

  MMC_values_smear_cut = new TH1F("MMC_values_smear_cut","MMC values smeared cuts",150,0,150);
      wk()->addOutput (MMC_values_smear_cut);

      */

  return EL::StatusCode::SUCCESS;

/*

     TCanvas* c1 = new TCanvas("c1","VBF TrueTaupT",3000,1000);
   c1->Divide(1,2);

  

   c1->cd(1);
   gStyle->SetOptStat(0);
   invarm->Draw();
   invarm->GetXaxis()->SetTitle("Invariant Mass (GeV)");
   invarm->GetYaxis()->SetTitle("counts");
   invarm_40->Draw("same"); invarm_40->SetLineColor(kRed);
   invarm_30->Draw("same"); invarm_30->SetLineColor(kBlue);
   invarm_eta->Draw("same"); invarm_eta->SetLineColor(kGreen);
   //invarm_met->Draw("same"); invarm_met->SetLineColor(kYellow+1);
   invarm_dR->Draw("same"); invarm_dR->SetLineColor(kViolet);
   invarm_deta->Draw("same"); invarm_deta->SetLineColor(kAzure);
   invarm_cent->Draw("same"); invarm_cent->SetLineColor(kCyan);


   c1->cd(2);
   gStyle->SetOptStat(0);
   smear_invarm->Draw();
   smear_invarm->GetXaxis()->SetTitle("Invariant Mass (GeV)");
   smear_invarm->GetYaxis()->SetTitle("counts");
   smear_invarm_40->Draw("same"); smear_invarm_40->SetLineColor(kRed);
   smear_invarm_30->Draw("same"); smear_invarm_30->SetLineColor(kBlue+1);
   smear_invarm_eta->Draw("same"); smear_invarm_eta->SetLineColor(kGreen);
   // smear_invarm_met->Draw("same"); smear_invarm_met->SetLineColor(kYellow+1);
   smear_invarm_dR->Draw("same"); smear_invarm_dR->SetLineColor(kViolet);
   smear_invarm_deta->Draw("same"); smear_invarm_deta->SetLineColor(kAzure);
   smear_invarm_cent->Draw("same"); smear_invarm_cent->SetLineColor(kCyan);

   auto smearlegend = new TLegend(0.65,0.7,0.90,0.9);
   smearlegend->AddEntry(smear_invarm,"ME_{T} > 20 GeV, #eta < 2.5");
   smearlegend->AddEntry(smear_invarm_40,"+ #tau_{0} p_{T} > 40 GeV");
   smearlegend->AddEntry(smear_invarm_30,"+ #tau_{1} p_{T} > 30 GeV");
   smearlegend->AddEntry(smear_invarm_eta,"+ #eta veto");
   //smearlegend->AddEntry(smear_invarm_met,"+ ME_{T} > 20 GeV");
   smearlegend->AddEntry(smear_invarm_dR,"+ 0.8 < dR < 2.4");
   smearlegend->AddEntry(smear_invarm_deta,"+ d#eta < 1.5");
   smearlegend->AddEntry(smear_invarm_cent,"+ ME_{T} centrality");
   smearlegend->Draw();
   */
}



EL::StatusCode truth_study :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode truth_study :: changeInput (bool firstFile)
{

  TTree *fChain=wk()->tree();
  fChain->SetBranchAddress("n_avg_int", &n_avg_int, &b_n_avg_int);
  fChain->SetBranchAddress("HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo", &HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo, &b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo);
  fChain->SetBranchAddress("HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I", &HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I, &b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I);
  fChain->SetBranchAddress("HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I_J25", &HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I_J25, &b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I_J25);
  fChain->SetBranchAddress("HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I", &HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I, &b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I);
  fChain->SetBranchAddress("HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I_J25", &HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I_J25, &b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I_J25);
  fChain->SetBranchAddress("HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1TAU20IM_2TAU12IM", &HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1TAU20IM_2TAU12IM, &b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1TAU20IM_2TAU12IM);
  fChain->SetBranchAddress("NOMINAL_pileup_combined_weight", &NOMINAL_pileup_combined_weight, &b_NOMINAL_pileup_combined_weight);
  fChain->SetBranchAddress("NOMINAL_pileup_random_run_number", &NOMINAL_pileup_random_run_number, &b_NOMINAL_pileup_random_run_number);
  fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
  fChain->SetBranchAddress("jet_NOMINAL_central_jets_global_effSF_JVT", &jet_NOMINAL_central_jets_global_effSF_JVT, &b_jet_NOMINAL_central_jets_global_effSF_JVT);
  fChain->SetBranchAddress("jet_NOMINAL_central_jets_global_ineffSF_JVT", &jet_NOMINAL_central_jets_global_ineffSF_JVT, &b_jet_NOMINAL_central_jets_global_ineffSF_JVT);
  fChain->SetBranchAddress("jet_NOMINAL_forward_jets_global_effSF_JVT", &jet_NOMINAL_forward_jets_global_effSF_JVT, &b_jet_NOMINAL_forward_jets_global_effSF_JVT);
  fChain->SetBranchAddress("jet_NOMINAL_forward_jets_global_ineffSF_JVT", &jet_NOMINAL_forward_jets_global_ineffSF_JVT, &b_jet_NOMINAL_forward_jets_global_ineffSF_JVT);
  fChain->SetBranchAddress("jet_NOMINAL_global_effSF_MVX", &jet_NOMINAL_global_effSF_MVX, &b_jet_NOMINAL_global_effSF_MVX);
  fChain->SetBranchAddress("jet_NOMINAL_global_ineffSF_MVX", &jet_NOMINAL_global_ineffSF_MVX, &b_jet_NOMINAL_global_ineffSF_MVX);
  fChain->SetBranchAddress("n_actual_int", &n_actual_int, &b_n_actual_int);
  fChain->SetBranchAddress("true_ditau_cosalpha", &true_ditau_cosalpha, &b_true_ditau_cosalpha);
  fChain->SetBranchAddress("true_ditau_deta", &true_ditau_deta, &b_true_ditau_deta);
  fChain->SetBranchAddress("true_ditau_dphi", &true_ditau_dphi, &b_true_ditau_dphi);
  fChain->SetBranchAddress("true_ditau_dpt", &true_ditau_dpt, &b_true_ditau_dpt);
  fChain->SetBranchAddress("true_ditau_dr", &true_ditau_dr, &b_true_ditau_dr);
  fChain->SetBranchAddress("true_ditau_mass", &true_ditau_mass, &b_true_ditau_mass);
  fChain->SetBranchAddress("true_ditau_ptx", &true_ditau_ptx, &b_true_ditau_ptx);
  fChain->SetBranchAddress("true_ditau_pty", &true_ditau_pty, &b_true_ditau_pty);
  fChain->SetBranchAddress("true_ditau_qxq", &true_ditau_qxq, &b_true_ditau_qxq);
  fChain->SetBranchAddress("true_ditau_scal_sum_pt", &true_ditau_scal_sum_pt, &b_true_ditau_scal_sum_pt);
  fChain->SetBranchAddress("true_ditau_vect_sum_pt", &true_ditau_vect_sum_pt, &b_true_ditau_vect_sum_pt);
  fChain->SetBranchAddress("true_ditau_vis_cosalpha", &true_ditau_vis_cosalpha, &b_true_ditau_vis_cosalpha);
  fChain->SetBranchAddress("true_ditau_vis_deta", &true_ditau_vis_deta, &b_true_ditau_vis_deta);
  fChain->SetBranchAddress("true_ditau_vis_dphi", &true_ditau_vis_dphi, &b_true_ditau_vis_dphi);
  fChain->SetBranchAddress("true_ditau_vis_dr", &true_ditau_vis_dr, &b_true_ditau_vis_dr);
  fChain->SetBranchAddress("true_ditau_vis_mass", &true_ditau_vis_mass, &b_true_ditau_vis_mass);
  fChain->SetBranchAddress("true_ditau_vis_scal_sum_pt", &true_ditau_vis_scal_sum_pt, &b_true_ditau_vis_scal_sum_pt);
  fChain->SetBranchAddress("true_ditau_vis_vect_sum_pt", &true_ditau_vis_vect_sum_pt, &b_true_ditau_vis_vect_sum_pt);
  fChain->SetBranchAddress("true_tau_0", &true_tau_0, &b_true_tau_0);
  
  fChain->SetBranchAddress("true_tau_0_charged_eta", &true_tau_0_charged_eta, &b_true_tau_0_charged_eta);
  fChain->SetBranchAddress("true_tau_0_charged_m", &true_tau_0_charged_m, &b_true_tau_0_charged_m);
  fChain->SetBranchAddress("true_tau_0_charged_phi", &true_tau_0_charged_phi, &b_true_tau_0_charged_phi);
  fChain->SetBranchAddress("true_tau_0_charged_pt", &true_tau_0_charged_pt, &b_true_tau_0_charged_pt);
  fChain->SetBranchAddress("true_tau_0_classifierParticleOrigin", &true_tau_0_classifierParticleOrigin, &b_true_tau_0_classifierParticleOrigin);
  fChain->SetBranchAddress("true_tau_0_classifierParticleType", &true_tau_0_classifierParticleType, &b_true_tau_0_classifierParticleType);
  fChain->SetBranchAddress("true_tau_0_decay_mode", &true_tau_0_decay_mode, &b_true_tau_0_decay_mode);
  fChain->SetBranchAddress("true_tau_0_et", &true_tau_0_et, &b_true_tau_0_et);
  fChain->SetBranchAddress("true_tau_0_eta", &true_tau_0_eta, &b_true_tau_0_eta);
  fChain->SetBranchAddress("true_tau_0_eta_vis", &true_tau_0_eta_vis, &b_true_tau_0_eta_vis);
  fChain->SetBranchAddress("true_tau_0_eta_vis_charged", &true_tau_0_eta_vis_charged, &b_true_tau_0_eta_vis_charged);
  fChain->SetBranchAddress("true_tau_0_eta_vis_neutral", &true_tau_0_eta_vis_neutral, &b_true_tau_0_eta_vis_neutral);
  fChain->SetBranchAddress("true_tau_0_isEle", &true_tau_0_isEle, &b_true_tau_0_isEle);
  fChain->SetBranchAddress("true_tau_0_isHadTau", &true_tau_0_isHadTau, &b_true_tau_0_isHadTau);
  fChain->SetBranchAddress("true_tau_0_isJet", &true_tau_0_isJet, &b_true_tau_0_isJet);
  fChain->SetBranchAddress("true_tau_0_isMuon", &true_tau_0_isMuon, &b_true_tau_0_isMuon);
  fChain->SetBranchAddress("true_tau_0_isTau", &true_tau_0_isTau, &b_true_tau_0_isTau);


  fChain->SetBranchAddress("true_tau_0_m", &true_tau_0_m, &b_true_tau_0_m);
  fChain->SetBranchAddress("true_tau_0_m_vis", &true_tau_0_m_vis, &b_true_tau_0_m_vis);
  fChain->SetBranchAddress("true_tau_0_m_vis_charged", &true_tau_0_m_vis_charged, &b_true_tau_0_m_vis_charged);
  fChain->SetBranchAddress("true_tau_0_m_vis_neutral", &true_tau_0_m_vis_neutral, &b_true_tau_0_m_vis_neutral);
  fChain->SetBranchAddress("true_tau_0_mother_pdgId", &true_tau_0_mother_pdgId, &b_true_tau_0_mother_pdgId);
  fChain->SetBranchAddress("true_tau_0_mother_status", &true_tau_0_mother_status, &b_true_tau_0_mother_status);
  fChain->SetBranchAddress("true_tau_0_n_charged", &true_tau_0_n_charged, &b_true_tau_0_n_charged);
  fChain->SetBranchAddress("true_tau_0_n_charged_pion", &true_tau_0_n_charged_pion, &b_true_tau_0_n_charged_pion);
  fChain->SetBranchAddress("true_tau_0_n_neutral", &true_tau_0_n_neutral, &b_true_tau_0_n_neutral);
  fChain->SetBranchAddress("true_tau_0_n_neutral_pion", &true_tau_0_n_neutral_pion, &b_true_tau_0_n_neutral_pion);
  fChain->SetBranchAddress("true_tau_0_neutral_eta", &true_tau_0_neutral_eta, &b_true_tau_0_neutral_eta);
  fChain->SetBranchAddress("true_tau_0_neutral_m", &true_tau_0_neutral_m, &b_true_tau_0_neutral_m);
  fChain->SetBranchAddress("true_tau_0_neutral_phi", &true_tau_0_neutral_phi, &b_true_tau_0_neutral_phi);
  fChain->SetBranchAddress("true_tau_0_neutral_pt", &true_tau_0_neutral_pt, &b_true_tau_0_neutral_pt);
  fChain->SetBranchAddress("true_tau_0_neutrino_eta", &true_tau_0_neutrino_eta, &b_true_tau_0_neutrino_eta);
  fChain->SetBranchAddress("true_tau_0_neutrino_m", &true_tau_0_neutrino_m, &b_true_tau_0_neutrino_m);
  fChain->SetBranchAddress("true_tau_0_neutrino_phi", &true_tau_0_neutrino_phi, &b_true_tau_0_neutrino_phi);
  fChain->SetBranchAddress("true_tau_0_neutrino_pt", &true_tau_0_neutrino_pt, &b_true_tau_0_neutrino_pt);
  fChain->SetBranchAddress("true_tau_0_origin", &true_tau_0_origin, &b_true_tau_0_origin);
  fChain->SetBranchAddress("true_tau_0_pdgId", &true_tau_0_pdgId, &b_true_tau_0_pdgId);
  fChain->SetBranchAddress("true_tau_0_phi", &true_tau_0_phi, &b_true_tau_0_phi);
  fChain->SetBranchAddress("true_tau_0_phi_vis", &true_tau_0_phi_vis, &b_true_tau_0_phi_vis);
  fChain->SetBranchAddress("true_tau_0_phi_vis_charged", &true_tau_0_phi_vis_charged, &b_true_tau_0_phi_vis_charged);
  fChain->SetBranchAddress("true_tau_0_phi_vis_neutral", &true_tau_0_phi_vis_neutral, &b_true_tau_0_phi_vis_neutral);
  fChain->SetBranchAddress("true_tau_0_pt", &true_tau_0_pt, &b_true_tau_0_pt);
  fChain->SetBranchAddress("true_tau_0_pt_vis", &true_tau_0_pt_vis, &b_true_tau_0_pt_vis);
  fChain->SetBranchAddress("true_tau_0_pt_vis_charged", &true_tau_0_pt_vis_charged, &b_true_tau_0_pt_vis_charged);
  fChain->SetBranchAddress("true_tau_0_pt_vis_neutral", &true_tau_0_pt_vis_neutral, &b_true_tau_0_pt_vis_neutral);
  fChain->SetBranchAddress("true_tau_0_pz", &true_tau_0_pz, &b_true_tau_0_pz);
  fChain->SetBranchAddress("true_tau_0_q", &true_tau_0_q, &b_true_tau_0_q);
  fChain->SetBranchAddress("true_tau_0_status", &true_tau_0_status, &b_true_tau_0_status);
  fChain->SetBranchAddress("true_tau_0_type", &true_tau_0_type, &b_true_tau_0_type);
  fChain->SetBranchAddress("true_tau_1", &true_tau_1, &b_true_tau_1);
  fChain->SetBranchAddress("true_tau_1_charged_eta", &true_tau_1_charged_eta, &b_true_tau_1_charged_eta);
  fChain->SetBranchAddress("true_tau_1_charged_m", &true_tau_1_charged_m, &b_true_tau_1_charged_m);
  fChain->SetBranchAddress("true_tau_1_charged_phi", &true_tau_1_charged_phi, &b_true_tau_1_charged_phi);
  fChain->SetBranchAddress("true_tau_1_charged_pt", &true_tau_1_charged_pt, &b_true_tau_1_charged_pt);
  fChain->SetBranchAddress("true_tau_1_classifierParticleOrigin", &true_tau_1_classifierParticleOrigin, &b_true_tau_1_classifierParticleOrigin);
  fChain->SetBranchAddress("true_tau_1_classifierParticleType", &true_tau_1_classifierParticleType, &b_true_tau_1_classifierParticleType);
  fChain->SetBranchAddress("true_tau_1_decay_mode", &true_tau_1_decay_mode, &b_true_tau_1_decay_mode);
  fChain->SetBranchAddress("true_tau_1_et", &true_tau_1_et, &b_true_tau_1_et);
  fChain->SetBranchAddress("true_tau_1_eta", &true_tau_1_eta, &b_true_tau_1_eta);
  fChain->SetBranchAddress("true_tau_1_eta_vis", &true_tau_1_eta_vis, &b_true_tau_1_eta_vis);
  fChain->SetBranchAddress("true_tau_1_eta_vis_charged", &true_tau_1_eta_vis_charged, &b_true_tau_1_eta_vis_charged);
  fChain->SetBranchAddress("true_tau_1_eta_vis_neutral", &true_tau_1_eta_vis_neutral, &b_true_tau_1_eta_vis_neutral);
  fChain->SetBranchAddress("true_tau_1_isEle", &true_tau_1_isEle, &b_true_tau_1_isEle);
  fChain->SetBranchAddress("true_tau_1_isHadTau", &true_tau_1_isHadTau, &b_true_tau_1_isHadTau);
  fChain->SetBranchAddress("true_tau_1_isJet", &true_tau_1_isJet, &b_true_tau_1_isJet);
  fChain->SetBranchAddress("true_tau_1_isMuon", &true_tau_1_isMuon, &b_true_tau_1_isMuon);
  fChain->SetBranchAddress("true_tau_1_isTau", &true_tau_1_isTau, &b_true_tau_1_isTau);
  fChain->SetBranchAddress("true_tau_1_m", &true_tau_1_m, &b_true_tau_1_m);
  fChain->SetBranchAddress("true_tau_1_m_vis", &true_tau_1_m_vis, &b_true_tau_1_m_vis);
  fChain->SetBranchAddress("true_tau_1_m_vis_charged", &true_tau_1_m_vis_charged, &b_true_tau_1_m_vis_charged);
  fChain->SetBranchAddress("true_tau_1_m_vis_neutral", &true_tau_1_m_vis_neutral, &b_true_tau_1_m_vis_neutral);
  fChain->SetBranchAddress("true_tau_1_mother_pdgId", &true_tau_1_mother_pdgId, &b_true_tau_1_mother_pdgId);
  fChain->SetBranchAddress("true_tau_1_mother_status", &true_tau_1_mother_status, &b_true_tau_1_mother_status);
  fChain->SetBranchAddress("true_tau_1_n_charged", &true_tau_1_n_charged, &b_true_tau_1_n_charged);
  fChain->SetBranchAddress("true_tau_1_n_charged_pion", &true_tau_1_n_charged_pion, &b_true_tau_1_n_charged_pion);
  fChain->SetBranchAddress("true_tau_1_n_neutral", &true_tau_1_n_neutral, &b_true_tau_1_n_neutral);
  fChain->SetBranchAddress("true_tau_1_n_neutral_pion", &true_tau_1_n_neutral_pion, &b_true_tau_1_n_neutral_pion);
  fChain->SetBranchAddress("true_tau_1_neutral_eta", &true_tau_1_neutral_eta, &b_true_tau_1_neutral_eta);
  fChain->SetBranchAddress("true_tau_1_neutral_m", &true_tau_1_neutral_m, &b_true_tau_1_neutral_m);
  fChain->SetBranchAddress("true_tau_1_neutral_phi", &true_tau_1_neutral_phi, &b_true_tau_1_neutral_phi);
  fChain->SetBranchAddress("true_tau_1_neutral_pt", &true_tau_1_neutral_pt, &b_true_tau_1_neutral_pt);
  fChain->SetBranchAddress("true_tau_1_neutrino_eta", &true_tau_1_neutrino_eta, &b_true_tau_1_neutrino_eta);
  fChain->SetBranchAddress("true_tau_1_neutrino_m", &true_tau_1_neutrino_m, &b_true_tau_1_neutrino_m);
  fChain->SetBranchAddress("true_tau_1_neutrino_phi", &true_tau_1_neutrino_phi, &b_true_tau_1_neutrino_phi);
  fChain->SetBranchAddress("true_tau_1_neutrino_pt", &true_tau_1_neutrino_pt, &b_true_tau_1_neutrino_pt);
  fChain->SetBranchAddress("true_tau_1_origin", &true_tau_1_origin, &b_true_tau_1_origin);
  fChain->SetBranchAddress("true_tau_1_pdgId", &true_tau_1_pdgId, &b_true_tau_1_pdgId);
  fChain->SetBranchAddress("true_tau_1_phi", &true_tau_1_phi, &b_true_tau_1_phi);
  fChain->SetBranchAddress("true_tau_1_phi_vis", &true_tau_1_phi_vis, &b_true_tau_1_phi_vis);
  fChain->SetBranchAddress("true_tau_1_phi_vis_charged", &true_tau_1_phi_vis_charged, &b_true_tau_1_phi_vis_charged);
  fChain->SetBranchAddress("true_tau_1_phi_vis_neutral", &true_tau_1_phi_vis_neutral, &b_true_tau_1_phi_vis_neutral);
  fChain->SetBranchAddress("true_tau_1_pt", &true_tau_1_pt, &b_true_tau_1_pt);
  fChain->SetBranchAddress("true_tau_1_pt_vis", &true_tau_1_pt_vis, &b_true_tau_1_pt_vis);
  fChain->SetBranchAddress("true_tau_1_pt_vis_charged", &true_tau_1_pt_vis_charged, &b_true_tau_1_pt_vis_charged);
  fChain->SetBranchAddress("true_tau_1_pt_vis_neutral", &true_tau_1_pt_vis_neutral, &b_true_tau_1_pt_vis_neutral);
  fChain->SetBranchAddress("true_tau_1_pz", &true_tau_1_pz, &b_true_tau_1_pz);
  fChain->SetBranchAddress("true_tau_1_q", &true_tau_1_q, &b_true_tau_1_q);
  fChain->SetBranchAddress("true_tau_1_status", &true_tau_1_status, &b_true_tau_1_status);
  fChain->SetBranchAddress("true_tau_1_type", &true_tau_1_type, &b_true_tau_1_type);
  fChain->SetBranchAddress("truth_met_et", &truth_met_et, &b_truth_met_et);
  fChain->SetBranchAddress("truth_met_etx", &truth_met_etx, &b_truth_met_etx);
  fChain->SetBranchAddress("truth_met_ety", &truth_met_ety, &b_truth_met_ety);
  fChain->SetBranchAddress("truth_met_phi", &truth_met_phi, &b_truth_met_phi);
  fChain->SetBranchAddress("truth_met_sig", &truth_met_sig, &b_truth_met_sig);
  fChain->SetBranchAddress("truth_met_sumet", &truth_met_sumet, &b_truth_met_sumet);
  fChain->SetBranchAddress("weight_mc", &weight_mc, &b_weight_mc);
  fChain->SetBranchAddress("weight_total", &weight_total, &b_weight_total);

  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode truth_study :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  int entry_base = 0;
   int entry_40 = 0;
   int entry_30 = 0;
   int entry_eta = 0;
   int entry_met = 0;
   int entry_dR = 0;
   int entry_deta = 0;
   int entry_cent = 0;

   int counter = 0;


   double smear_entry_base = 0;
   double smear_entry_40 = 0;
   double smear_entry_30 = 0;
   double smear_entry_eta = 0;
   double smear_entry_met = 0;
   double smear_entry_dR = 0;
   double smear_entry_deta = 0;
   double smear_entry_cent = 0;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode truth_study :: execute ()
{


  counter++;
  std::cout<<"\nLoop number: "<<counter<<std::endl;

   TLorentzVector v0;
   TLorentzVector v1;
   TLorentzVector v0_smearpT;
   TLorentzVector v1_smearpT;
  
   wk()->tree()->GetEntry (wk()->treeEntry());

  //put in old cutflow here


  if ((true_tau_0_n_charged_pion == 1 || true_tau_0_n_charged_pion == 3) && (true_tau_1_n_charged_pion ==1 || true_tau_1_n_charged_pion==3) 
    &&  truth_met_et > 20 
    && abs(true_tau_0_eta_vis)<2.5 && abs(true_tau_1_eta_vis)<2.5){

     float weight = getTauEfficiency(truth_met_et*1000,true_tau_0_eta_vis,true_tau_0_n_charged_pion);
     float trig_weight = getDiTauTriggerEfficiency(true_tau_0_pt_vis*1000,true_tau_1_pt_vis*1000, true_tau_0_eta_vis, true_tau_1_eta_vis,true_tau_0_n_charged_pion,true_tau_1_n_charged_pion);
     
     float xc_ggH = 5467; //in fb at 125 GeV, N3LO
     float xc_VBF = 4278; //in fb at 125 GeV 
     float lum = 3000;

     float fullweight = weight*trig_weight;

     //smearing the pT
     float smearpT_0 = getTauSmearedEnergy(true_tau_0_pt_vis*1000,true_tau_0_eta_vis,true_tau_0_n_charged_pion)/1000.0;
     float smearpT_1 = getTauSmearedEnergy(true_tau_1_pt_vis*1000,true_tau_1_eta_vis,true_tau_1_n_charged_pion)/1000.0;

     //forming the unsmeared lorentz vectors
       v0.SetPtEtaPhiM(true_tau_0_pt_vis,true_tau_0_eta_vis,true_tau_0_phi_vis,true_tau_0_m_vis); 
       v1.SetPtEtaPhiM(true_tau_1_pt_vis,true_tau_1_eta_vis,true_tau_1_phi_vis,true_tau_1_m_vis);

      //forming the lorentz vectors with smeared pT
      v0_smearpT.SetPtEtaPhiM(smearpT_0,true_tau_0_eta_vis,true_tau_0_phi_vis,true_tau_0_m_vis); 
      v1_smearpT.SetPtEtaPhiM(smearpT_1,true_tau_1_eta_vis,true_tau_1_phi_vis,true_tau_1_m_vis);   

      //calculating invariant mass unsmeared and smeared
      //double mag = sqrt(v0.Dot(v1));
      double smearmag = sqrt(v0_smearpT.Dot(v1_smearpT));

      //invarm->Fill(mag);
      //invarm_initial->Fill(mag);
      //smear_invarm->Fill(smearmag,fullweight);
      //smear_invarm_initial->Fill(smearmag,fullweight);

      entry_base++;
      smear_entry_base+=fullweight;

////////////////////////////////////// MMC for smeared MET //////////////////////////////////////
////////////////////////////////////// MMC for smeared MET //////////////////////////////////////



      MissingMassCalculator fMMC;
      fMMC.SetUseEfficiencyRecovery(1);
      fMMC.SetUseTailCleanup(0);
      fMMC.SetCalibrationSet(MMCCalibrationSet::UPGRADE); 
      //fMMC.SetSumEt(truth_met_sumet);   
      //fMMC.Apply();    

      //***ALL UNITS NEED TO BE EXPRESSED IN GeV (not MeV as ATLAS default)***

      TVector2 smeared_met = getMETSmeared(truth_met_sumet*1000,truth_met_etx*1000,truth_met_ety*1000,nominal)/1000;
      /*std::cout<<"X: \n"<< truth_met_etx << std::endl;
      std::cout<<"X smeared: \n"<< smeared_met.X() << std::endl;
      std::cout<<"Y: \n" << truth_met_ety << std::endl;
      std::cout<<"Y smeared: \n" << smeared_met.Y() << std::endl; */

      double smear_met_et;
      smear_met_et = hypot(smeared_met.X(),smeared_met.Y());
      double met_diff = truth_met_et - smear_met_et;
      double met_diffx = truth_met_etx - smeared_met.X();
      double met_diffy = truth_met_ety - smeared_met.Y();


      met_hist->Fill(truth_met_et);
      metx_hist->Fill(truth_met_etx);
      mety_hist->Fill(truth_met_ety);
      smear_met_hist->Fill(smear_met_et);
      smear_metx_hist->Fill(smeared_met.X());
      smear_mety_hist->Fill(smeared_met.Y());
      met_diff_hist->Fill(met_diff);
      met_diffx_hist->Fill(met_diffx);
      met_diffy_hist->Fill(met_diffy);


      float METresolution = getMETResolution(truth_met_sumet*1000,nominal)/1000;

      fMMC.SetMetVec(smeared_met); // input: TVector2 for MET smeared by METSmearing package

      fMMC.SetMetScanParams(0.0,METresolution,METresolution); // METresolution should be obtained from getMETandRes() in METSmearing package
      
      //TLorentzVector smeared_tau0_vec; //can be either el, mu or tau according to the channel         
      //TLorentzVector smeared_tau1_vec; //can be either el, mu or tau according to the channel         
      fMMC.SetVisTauVec(0,v0_smearpT); // smeared tau 4-vector; please read instructions above to properly set tau mass
      fMMC.SetVisTauVec(1,v1_smearpT); // smeared tau 4-vector; please read instructions above to properly set tau mass

      int tau0Type, tau1Type; 
      if (true_tau_0_n_charged_pion=1){tau0Type=10;}
      if (true_tau_0_n_charged_pion=3){tau0Type=30;}
      if (true_tau_1_n_charged_pion=1){tau1Type=10;}
      if (true_tau_1_n_charged_pion=3){tau1Type=30;}
      fMMC.SetVisTauType(0,tau0Type); // tau*Type=10 for 1-prongs and tau*Type=30 for 3-prongs
      fMMC.SetVisTauType(1,tau1Type);     


////////////////////////////////////// MMC for truth MET //////////////////////////////////////
////////////////////////////////////// MMC for truth MET //////////////////////////////////////

      MissingMassCalculator fMMC_truth;
      fMMC_truth.SetUseEfficiencyRecovery(1);
      fMMC_truth.SetUseTailCleanup(0);
      fMMC_truth.SetCalibrationSet(MMCCalibrationSet::UPGRADE);    
      //fMMC_truth.SetSumEt(truth_met_sumet);   
      //fMMC_truth.Apply(); 

      float METresolution_truth = 0;

      TVector2 unsmeared_met;
      unsmeared_met.SetX(truth_met_etx); unsmeared_met.SetY(truth_met_ety);
      fMMC_truth.SetMetVec(unsmeared_met); // input: TVector2 for MET smeared by METSmearing package

      fMMC_truth.SetMetScanParams(0.0,METresolution_truth,METresolution_truth); // METresolution should be obtained from getMETandRes() in METSmearing package
      
      //TLorentzVector smeared_tau0_vec; //can be either el, mu or tau according to the channel         
      //TLorentzVector smeared_tau1_vec; //can be either el, mu or tau according to the channel         
      fMMC_truth.SetVisTauVec(0,v0); // smeared tau 4-vector; please read instructions above to properly set tau mass
      fMMC_truth.SetVisTauVec(1,v1); // smeared tau 4-vector; please read instructions above to properly set tau mass

      

      /*int tau0Type_truth, tau1Type_truth; 
      if (true_tau_0_n_charged_pion=1){tau0Type_truth=10;}
      if (true_tau_0_n_charged_pion=3){tau0Type_truth=30;}
      if (true_tau_1_n_charged_pion=1){tau1Type_truth=10;}
      if (true_tau_1_n_charged_pion=3){tau1Type_truth=30;} */
      //fMMC_truth.SetVisTauType(0,tau0Type); // tau*Type=10 for 1-prongs and tau*Type=30 for 3-prongs
      //fMMC_truth.SetVisTauType(1,tau1Type);

      


/////////GETTING ALL MMC (SMEARED AND TRUTH) RESULTS//////////


      int mmc_stattest = fMMC.RunMissingMassCalculator(); // run MMC
      int    mmc_stat = fMMC.GetFitStatus();
      std::cout<<mmc_stat<<std::endl;   
      
      int mmc_stattest_truth = fMMC_truth.RunMissingMassCalculator(); // run MMC
      int    mmc_stat_truth = fMMC_truth.GetFitStatus();
      std::cout<<mmc_stat_truth<<std::endl; 

      if(mmc_stat==1)
        {
          double mmc_mass1 = fMMC.GetFittedMass(MMCFitMethod::MLM); // recommended way to access MMC mas
          std::cout<<"MMC Mass: " << mmc_mass1 << std::endl;    
          double mmc_mass1_truth = fMMC_truth.GetFittedMass(MMCFitMethod::MLM); // recommended way to access MMC mas
          std::cout<<"\n MMC Mass (truth param): " << mmc_mass1_truth << std::endl;  
          smear_MMC->Fill(mmc_mass1,fullweight);
          MMC_truth->Fill(mmc_mass1_truth);

          

        if (true_tau_0_pt_vis > 40){

          //invarm_40->Fill(mag);
          entry_40++;
          MMC_40_truth->Fill(mmc_mass1_truth);
          //tau1_precut -> Fill(true_tau_1_pt_vis);

          if (true_tau_1_pt_vis > 30){

            //invarm_30_initial->Fill(mag);

            //invarm_30->Fill(mag);
            MMC_30_truth->Fill(mmc_mass1_truth);
             // eta_precut->Fill(true_tau_0_eta_vis);
              entry_30++;
             // tau1_cut->Fill(true_tau_1_pt_vis);


            if ((abs(true_tau_0_eta_vis)<1.37 || abs(true_tau_0_eta_vis)>1.52) 
                && (abs(true_tau_1_eta_vis)<1.37 || abs(true_tau_1_eta_vis)>1.52)){

              //invarm_eta->Fill(mag);
              MMC_eta_truth->Fill(mmc_mass1_truth);

                //eta_veto->Fill(true_tau_0_eta_vis);
                //met_precut->Fill(truth_met_et);
                entry_eta++;

                    //dR_precut -> Fill(v0.DeltaR(v1));

                if (truth_met_et > 20.0){

                  //invarm_met -> Fill(mag);
                  MMC_met_truth->Fill(mmc_mass1_truth);
                  entry_met++;

                  if (v0.DeltaR(v1)<2.4 && v0.DeltaR(v1) > 0.8){

                    //invarm_dR->Fill(mag);
                    MMC_dR_truth->Fill(mmc_mass1_truth);
                    //dR_cut->Fill(v0.DeltaR(v1));
                    entry_dR++;
                    //deta_precut->Fill(true_ditau_deta);

                    if (true_ditau_deta<1.5){

                      //invarm_deta->Fill(mag);
                      MMC_deta_truth->Fill(mmc_mass1_truth);
                      //deta_cut->Fill(true_ditau_deta);
                      entry_deta++;

                      if (METinDitauDPhi(truth_met_phi,v0,v1)==1){
                        //m_vis_cent->Fill(true_tau_0_m_vis);
                        //m1_vis_cent->Fill(true_tau_1_m_vis);
                        //invarm_cent->Fill(mag);
                        MMC_cent_truth->Fill(mmc_mass1_truth);
                        //invarm_final->Fill(mag);
                        //met_cent->Fill(truth_met_et);
                        entry_cent++;
                      }
                    }
                  }
                }
              }
            }
          }

          
  

        if (smearpT_0 > 40){

                 //smear_invarm_40->Fill(smearmag,fullweight);
                 smear_MMC_40->Fill(mmc_mass1,fullweight);
                 smear_entry_40+=fullweight; 
                //tau1_precut -> Fill(true_tau_1_pt_vis);

          //cout<<smearpT_1-30;
          //cout<<"\n";

          if (smearpT_1  > 30){

            //smear_invarm_30->Fill(smearmag,fullweight);
            smear_MMC_30->Fill(mmc_mass1,fullweight);
            //eta_precut->Fill(true_tau_0_eta_vis);
            smear_entry_30+=fullweight;
            //tau1_cut->Fill(true_tau_1_pt_vis);



            if ((abs(true_tau_0_eta_vis)<1.37 || abs(true_tau_0_eta_vis)>1.52) 
                && (abs(true_tau_1_eta_vis)<1.37 || abs(true_tau_1_eta_vis)>1.52)){

                //smear_invarm_eta->Fill(smearmag,fullweight);
                smear_MMC_eta->Fill(mmc_mass1,fullweight);

                //eta_veto->Fill(true_tau_0_eta_vis);
                //met_precut->Fill(truth_met_et);
                smear_entry_eta+=fullweight;

                    //dR_precut -> Fill(v0.DeltaR(v1));

                if (smear_met_et > 20.0){

                  //smear_invarm_met->Fill(smearmag,fullweight);
                  smear_MMC_met->Fill(mmc_mass1,fullweight);
                  smear_entry_met+=fullweight;


                  if (v0.DeltaR(v1)<2.4 && v0.DeltaR(v1) > 0.8){

                      //smear_invarm_dR->Fill(smearmag,fullweight);
                      smear_MMC_dR->Fill(mmc_mass1,fullweight);
                      //dR_cut->Fill(v0.DeltaR(v1));
                      smear_entry_dR+=fullweight;
                      //deta_precut->Fill(true_ditau_deta);

                    if (true_ditau_deta<1.5){

                      //smear_invarm_deta->Fill(smearmag,fullweight);
                      smear_MMC_deta->Fill(mmc_mass1,fullweight);
                      //deta_cut->Fill(true_ditau_deta);
                      smear_entry_deta+=fullweight;

                      if (METinDitauDPhi(truth_met_phi,v0,v1)==1){
                        //smear_invarm_cent->Fill(smearmag,fullweight);
                        //smear_invarm_final->Fill(smearmag,fullweight);
                        smear_MMC_cent->Fill(mmc_mass1,fullweight);
                        //met_cent->Fill(truth_met_et);
                        smear_entry_cent+=fullweight;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

  
  //std::cout<<true_tau_0_pt_vis<<std::endl;


  return EL::StatusCode::SUCCESS;


}



EL::StatusCode truth_study :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode truth_study :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.


  std::cout << "\n ////////////// Start Unsmeared ////////////// \n";

  std::cout << "Total events: ";
  std::cout << entry_base;

  std::cout << "\nafter #tau_{0} > 40 GeV: ";
  std::cout << entry_40;

  std::cout << "\nafter #tau_{1} > 30 GeV:";
  std::cout << entry_30;

  std::cout << "\nafter #eta veto: ";
  std::cout << entry_eta;

  std::cout << "\n after ME_{T} > 20 GeV: ";
  std::cout<<entry_met;

  std::cout << "\nafter 0.8< dR <2.4: ";
  std::cout << entry_dR;

  std::cout << "\nafter d#eta<1.5: ";
  std::cout << entry_deta;

  std::cout << "\nafter ME_{miss} cent: ";
  std::cout << entry_cent;

  std::cout << "\n ////////////// Finish Unsmeared ////////////// \n";


  std::cout << "\n ////////////// Start Smeared ////////////// \n";

  std::cout << "Total events: ";
  std::cout << smear_entry_base;

  std::cout << "\nafter #tau_{0} > 40 GeV: ";
  std::cout << smear_entry_40;

  std::cout << "\nafter #tau_{1} > 30 GeV:";
  std::cout << smear_entry_30;

  std::cout << "\nafter #eta veto: ";
  std::cout << smear_entry_eta;

  std::cout << "\n after ME_{T} > 20 GeV: ";
  std::cout<< smear_entry_met;

  std::cout << "\nafter 0.8< dR <2.4: ";
  std::cout << smear_entry_dR;

  std::cout << "\nafter d#eta<1.5: ";
  std::cout << smear_entry_deta;

  std::cout << "\nafter ME_{miss} cent: ";
  std::cout << smear_entry_cent;

  std::cout << "\n ////////////// Finish Smeared ////////////// \n";


   
   TCanvas* c1 = new TCanvas("c1","Ztautau UPGRADE sumet",2000,1000);
   c1->Divide(2,1);

   c1->cd(1);
   gStyle->SetOptStat(0);

   MMC_truth->Draw();
   MMC_truth->GetXaxis()->SetTitle("MMC Mass(GeV)"); MMC_truth->SetTitle("MMC Z#rightarrow#tau#tau truth with sumET input ");
   smear_MMC->GetYaxis()->SetTitle("counts");
   MMC_40_truth->Draw("same"); MMC_40_truth->SetLineColor(kRed);
   MMC_30_truth->Draw("same"); MMC_30_truth->SetLineColor(kBlue);
   MMC_eta_truth->Draw("same"); MMC_eta_truth->SetLineColor(kGreen);
   MMC_met_truth->Draw("same"); MMC_met_truth->SetLineColor(kYellow+1);
   MMC_dR_truth->Draw("same"); MMC_dR_truth->SetLineColor(kViolet);
   MMC_deta_truth->Draw("same"); MMC_deta_truth->SetLineColor(kAzure);
   MMC_cent_truth->Draw("same"); MMC_cent_truth->SetLineColor(kCyan);
    c1->Update();

   auto legend_truth = new TLegend(0.65,0.7,0.90,0.9);
   legend_truth->AddEntry(MMC_truth,"ME_{T} > 20 GeV and #eta < 2.5");
   legend_truth->AddEntry(MMC_40_truth,"+ #tau_{0} p_{T} > 40 GeV");
   legend_truth->AddEntry(MMC_30_truth,"+ #tau_{1} p_{T} > 30 GeV");
   legend_truth->AddEntry(MMC_eta_truth,"+ #eta veto");
   legend_truth->AddEntry(MMC_met_truth,"+ ME_{T} > 20 GeV");
   legend_truth->AddEntry(MMC_dR_truth,"+ 0.8 < dR < 2.4");
   legend_truth->AddEntry(MMC_deta_truth,"+ d#eta < 1.5");
   legend_truth->AddEntry(MMC_cent_truth,"+ ME_{T} centrality");
   legend_truth->Draw();

   c1->cd(2);

   gStyle->SetOptStat(0);
   smear_MMC->Draw();
   smear_MMC->GetXaxis()->SetTitle("MMC Mass(GeV)"); smear_MMC->SetTitle("MMC Z#rightarrow#tau#tau smeared with sumEt input");
   smear_MMC->GetYaxis()->SetTitle("counts");
   smear_MMC_40->Draw("same"); smear_MMC_40->SetLineColor(kRed);
   smear_MMC_30->Draw("same"); smear_MMC_30->SetLineColor(kBlue);
   smear_MMC_eta->Draw("same"); smear_MMC_eta->SetLineColor(kGreen);
   smear_MMC_met->Draw("same"); smear_MMC_met->SetLineColor(kYellow+1);
   smear_MMC_dR->Draw("same"); smear_MMC_dR->SetLineColor(kViolet);
   smear_MMC_deta->Draw("same"); smear_MMC_deta->SetLineColor(kAzure);
   smear_MMC_cent->Draw("same"); smear_MMC_cent->SetLineColor(kCyan);
    c1->Update();

   auto legend = new TLegend(0.65,0.7,0.90,0.9);
   legend->AddEntry(smear_MMC,"ME_{T} > 20 GeV and #eta < 2.5");
   legend->AddEntry(smear_MMC_40,"+ #tau_{0} p_{T} > 40 GeV");
   legend->AddEntry(smear_MMC_30,"+ #tau_{1} p_{T} > 30 GeV");
   legend->AddEntry(smear_MMC_eta,"+ #eta veto");
   legend->AddEntry(smear_MMC_met,"+ ME_{T} > 20 GeV");
   legend->AddEntry(smear_MMC_dR,"+ 0.8 < dR < 2.4");
   legend->AddEntry(smear_MMC_deta,"+ d#eta < 1.5");
   legend->AddEntry(smear_MMC_cent,"+ ME_{T} centrality");
   legend->Draw();
   c1->Update();


   c1->SaveAs("170727_MMC_Ztautau_2million_UPGRADE_sumet.pdf");




//////////////////////////////////////////////////////////////////////////////


  /*

  TCanvas* c1 = new TCanvas("c1","MET, half smear",2000,800);
    c1->Divide(3,1);

c1->cd(1);

   gStyle->SetOptStat(0);


  met_hist->Draw();
    met_hist->SetTitle("ME_{T}, ggH");
    met_hist->GetXaxis()->SetTitle("MET (GeV)");
    met_hist->GetYaxis()->SetTitle("counts");
  smear_met_hist->Draw("same"); smear_met_hist->SetLineColor(kRed);
  met_diff_hist->Draw("same"); met_diff_hist->SetLineColor(kBlue);

  auto legend = new TLegend(0.65,0.7,0.90,0.9);
  legend->AddEntry(met_hist,"truth MET, ggH");
  legend->AddEntry(smear_met_hist,"smeared MET, ggH");
  legend->AddEntry(met_diff_hist,"truth - smeared MET");
  legend->Draw();
  c1->Update();

c1->cd(2);

  metx_hist->Draw();
    metx_hist->SetTitle("ME_{T}x, ggH"); 
    metx_hist->GetXaxis()->SetTitle("MET (GeV)");
    metx_hist->GetYaxis()->SetTitle("counts");
  smear_metx_hist->Draw("same"); smear_metx_hist->SetLineColor(kRed);
  met_diffx_hist->Draw("same"); met_diffx_hist->SetLineColor(kBlue);

  auto legendx = new TLegend(0.65,0.7,0.90,0.9);
  legendx->AddEntry(metx_hist,"truth METx, ggH");
  legendx->AddEntry(smear_metx_hist,"smeared METx ggH");
  legendx->AddEntry(met_diffx_hist,"truth - smeared METx");
  legendx->Draw();
  c1->Update();

c1->cd(3);

  mety_hist->Draw();
    mety_hist->SetTitle("ME_{T}y, ggH");
    mety_hist->GetXaxis()->SetTitle("MET (GeV)");
    mety_hist->GetYaxis()->SetTitle("counts");
  smear_mety_hist->Draw("same"); smear_mety_hist->SetLineColor(kRed);
  met_diffy_hist->Draw("same"); met_diffy_hist->SetLineColor(kBlue);

  auto legendy = new TLegend(0.65,0.7,0.90,0.9);
  legendy->AddEntry(mety_hist,"truth METy, ggH");
  legendy->AddEntry(smear_mety_hist,"smeared METy, ggH");
  legendy->AddEntry(met_diffy_hist,"truth - smeared METy");
  legendy->Draw();
  c1->Update();



c1->SaveAs("170727_MET_ggH_sumet.pdf");

*/

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode truth_study :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.




/*

*/

  return EL::StatusCode::SUCCESS;
}

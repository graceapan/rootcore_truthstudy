#ifndef truth_study_truth_study_H
#define truth_study_truth_study_H
#include <TH1.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include "TBranch.h"
#include <vector>
#include <TVector2.h>
#include <vector>
#include <TLorentzVector.h>
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include <EventLoop/Algorithm.h>
#include <DiTauMassTools/MissingMassCalculator.h>
#include <DiTauMassTools/MissingMassCalculator.h>
#include "DiTauMassTools/MissingMassTool.h"
#include "DiTauMassTools/DiTauMassToolsDict.h"
#include "DiTauMassTools/IMissingMassTool.h"


class truth_study : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  TH1F *invarm; //!
  TH1F *invarm_40; //!
  TH1F *invarm_30; //!
  TH1F *invarm_eta; //!
  TH1F *invarm_met; //!
  TH1F *invarm_dR; //!
  TH1F *invarm_deta; //!
  TH1F *invarm_cent; //!

  TH1F *smear_invarm; //!
  TH1F *smear_invarm_40; //!
  TH1F *smear_invarm_30; //!
  TH1F *smear_invarm_eta; //!
  TH1F *smear_invarm_met; //!
  TH1F *smear_invarm_dR; //!
  TH1F *smear_invarm_deta; //!
  TH1F *smear_invarm_cent; //!

  TH1F *met_hist; //!
  TH1F *smear_met_hist; //!
  TH1F *metx_hist; //!
  TH1F *smear_metx_hist; //!
  TH1F *mety_hist; //!
  TH1F *smear_mety_hist; //!
  TH1F *met_diff_hist; //!
  TH1F *met_diffx_hist; //!
  TH1F *met_diffy_hist; //!


  

  TH1F *MMC_truth; //!
  TH1F *MMC_40_truth; //!
  TH1F *MMC_30_truth; //!
  TH1F *MMC_eta_truth; //!
  TH1F *MMC_met_truth; //!
  TH1F *MMC_dR_truth; //!
  TH1F *MMC_deta_truth; //!
  TH1F *MMC_cent_truth; //!

  

  TH1F *smear_MMC; //!
  TH1F *smear_MMC_40; //!
  TH1F *smear_MMC_30; //!
  TH1F *smear_MMC_eta; //!
  TH1F *smear_MMC_met; //!
  TH1F *smear_MMC_dR; //!
  TH1F *smear_MMC_deta; //!
  TH1F *smear_MMC_cent; //!

   int entry_base = 0; //!
   int entry_40 = 0; //!
   int entry_30 = 0; //!
   int entry_eta = 0; //!
   int entry_met = 0; //!
   int entry_dR = 0; //!
   int entry_deta = 0; //!
   int entry_cent = 0; //!

   int counter=0; //!

   double smear_entry_base = 0; //!
   double smear_entry_40 = 0; //!
   double smear_entry_30 = 0; //!
   double smear_entry_eta = 0; //!
   double smear_entry_met = 0; //!
   double smear_entry_dR = 0; //!
   double smear_entry_deta = 0; //!
   double smear_entry_cent = 0; //!

   TCanvas *c1; //!
   TLegend *legend; //!


  //  MissingMassCalculator fMMC;
  TTree          *fChain;   //!   
  Int_t           fCurrent;  //!  

  Float_t         true_tau_0_neutral_eta;   //!         
  Float_t         true_tau_0_neutral_m;  //!         
  Float_t         true_tau_0_neutral_phi;  //!         
  Float_t         true_tau_0_neutral_pt;  //!         
  UInt_t          HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo;   //!         
  UInt_t          HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I;   //!         
  UInt_t          HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I_J25;  //!         
  UInt_t          HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I;  //!         
  UInt_t          HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I_J25;  //!         
  UInt_t          HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1TAU20IM_2TAU12IM;  //!         
  Float_t         NOMINAL_pileup_combined_weight;  //!         
  UInt_t          NOMINAL_pileup_random_run_number;  //!         
  ULong64_t       event_number;  //!         
  Float_t         jet_NOMINAL_central_jets_global_effSF_JVT;   //!         
  Float_t         jet_NOMINAL_central_jets_global_ineffSF_JVT;  //!         
  Float_t         jet_NOMINAL_forward_jets_global_effSF_JVT;  //!         
  Float_t         jet_NOMINAL_forward_jets_global_ineffSF_JVT;  //!         
  Float_t         jet_NOMINAL_global_effSF_MVX;  //!         
  Float_t         jet_NOMINAL_global_ineffSF_MVX;  //!         
  Float_t         n_actual_int;  //!         

  Float_t         n_avg_int;//!

  Float_t         n_avg_int_cor;  //!         
  Int_t           n_bjets;  //!         
  Int_t           n_electrons;  //!         
  Int_t           n_jets;  //!         
  Int_t           n_jets_bad;  //!         
  Int_t           n_muons;  //!         
  Int_t           n_photons;  //!         
  Int_t           n_pvx;  //!         
  Int_t           n_taus;  //!         
  Int_t           n_taus_loose;  //!         
  Int_t           n_taus_medium;  //!         
  Int_t           n_taus_tight;  //!         
  UInt_t          n_truth_gluon_jets;  //!         
  UInt_t          n_truth_jets;  //!         
  UInt_t          n_truth_jets_pt20_eta45;  //!         
  UInt_t          n_truth_quark_jets;  //!         
  Int_t           n_vx;  //!         
  UInt_t          run_number;  //!         
  UInt_t          sherpa_n_sherpa_truth_jets_pt20_eta45;   //!         
  Float_t         sherpa_weight;  //!         
  Int_t           true_ditau;  //!         
  Float_t         true_ditau_CP_alphaminus_ip;  //!         
  Float_t         true_ditau_CP_alphaminus_ip_rho;  //!         
  Float_t         true_ditau_CP_alphaminus_rho_rho;  //!         
  Float_t         true_ditau_CP_ip_tau0_mag;  //!         
  Float_t         true_ditau_CP_ip_tau0_x_ip;  //!         
  Float_t         true_ditau_CP_ip_tau0_y_ip;  //!         
  Float_t         true_ditau_CP_ip_tau0_z_ip;  //!         
  Float_t         true_ditau_CP_ip_tau1_mag;  //!         
  Float_t         true_ditau_CP_ip_tau1_x_ip;  //!         
  Float_t         true_ditau_CP_ip_tau1_y_ip;  //!         
  Float_t         true_ditau_CP_ip_tau1_z_ip;  //!         
  Float_t         true_ditau_CP_phi_star_cp_a1_rho;  //!         
  Float_t         true_ditau_CP_phi_star_cp_ip;  //!         
  Float_t         true_ditau_CP_phi_star_cp_ip_rho;  //!         
  Float_t         true_ditau_CP_phi_star_cp_ip_rho_opt;  //!         
  Float_t         true_ditau_CP_phi_star_cp_rho_ip;  //!         
  Float_t         true_ditau_CP_phi_star_cp_rho_rho;  //!         
  Float_t         true_ditau_CP_tau0_upsilon;  //!         
  Float_t         true_ditau_CP_tau1_upsilon;  //!         
  Float_t         true_ditau_cosalpha;  //!         
  Float_t         true_ditau_deta;  //!         
  Float_t         true_ditau_dphi;  //!         
  Float_t         true_ditau_dpt;  //!         
  Float_t         true_ditau_dr;  //!         
  Float_t         true_ditau_mass;  //!         
  Float_t         true_ditau_ptx;  //!         
  Float_t         true_ditau_pty;  //!         
  Float_t         true_ditau_qxq;  //!         
  Float_t         true_ditau_scal_sum_pt;  //!         
  Float_t         true_ditau_vect_sum_pt;  //!         
  Float_t         true_ditau_vis_cosalpha;  //!         
  Float_t         true_ditau_vis_deta;  //!         
  Float_t         true_ditau_vis_dphi;  //!         
  Float_t         true_ditau_vis_dr;  //!         
  Float_t         true_ditau_vis_mass;  //!         
  Float_t         true_ditau_vis_scal_sum_pt;  //!         
  Float_t         true_ditau_vis_vect_sum_pt;  //!         
  UInt_t          true_tau_0;  //!         
  Float_t         true_tau_0_charged_eta;  //!         
  Float_t         true_tau_0_charged_m;  //!         
  Float_t         true_tau_0_charged_phi;  //!         
  Float_t         true_tau_0_charged_pt;  //!         
  UInt_t          true_tau_0_classifierParticleOrigin;  //!         
  UInt_t          true_tau_0_classifierParticleType;  //!         
  Int_t           true_tau_0_decay_mode;  //!         
  Float_t         true_tau_0_et;  //!         
  Float_t         true_tau_0_eta;  //!         
  Float_t         true_tau_0_eta_vis;  //!         
  Double_t        true_tau_0_eta_vis_charged;  //!         
  Double_t        true_tau_0_eta_vis_neutral;  //!         
  Int_t           true_tau_0_n_neutral_pion;  //!         
  Int_t           true_tau_0_isEle;  //!         
  Int_t           true_tau_0_isHadTau;  //!         
  Int_t           true_tau_0_isJet;  //!         
  Int_t           true_tau_0_isMuon;  //!         
  Int_t           true_tau_0_isTau;  //!         
  Float_t         true_tau_0_m;  //!         
  Float_t         true_tau_0_m_vis;  //!         
  Double_t        true_tau_0_m_vis_charged;  //!         
  Double_t        true_tau_0_m_vis_neutral;  //!         
  Int_t           true_tau_0_mother_pdgId;  //!         
  Int_t           true_tau_0_mother_status;  //!         
  Int_t           true_tau_0_n_charged;  //!         
  Int_t           true_tau_0_n_charged_pion;  //!         
  Int_t           true_tau_0_n_neutral;  //!         
  Float_t         true_tau_0_neutrino_eta;  //!         
  Float_t         true_tau_0_neutrino_m;  //!         
  Float_t         true_tau_0_neutrino_phi;  //!         
  Float_t         true_tau_0_neutrino_pt;  //!         
  Int_t           true_tau_0_origin;  //!         
  Int_t           true_tau_0_pdgId;  //!         
  Float_t         true_tau_0_phi;  //!         
  Float_t         true_tau_0_phi_vis;  //!         
  Double_t        true_tau_0_phi_vis_charged;  //!         
  Double_t        true_tau_0_phi_vis_neutral;  //!         
  Float_t         true_tau_0_pt;  //!         
  Float_t         true_tau_0_pt_vis;  //!         
  Double_t        true_tau_0_pt_vis_charged;  //!         
  Double_t        true_tau_0_pt_vis_neutral;  //!         
  Float_t         true_tau_0_pz;  //!         
  Float_t         true_tau_0_q;  //!         
  Int_t           true_tau_0_status;  //!         
  Int_t           true_tau_0_type;  //!         
  UInt_t          true_tau_1;  //!         
  Float_t         true_tau_1_charged_eta;  //!         
  Float_t         true_tau_1_charged_m;  //!         
  Float_t         true_tau_1_charged_phi;  //!         
  Float_t         true_tau_1_charged_pt;  //!         
  UInt_t          true_tau_1_classifierParticleOrigin;  //!         
  UInt_t          true_tau_1_classifierParticleType;  //!         
  Int_t           true_tau_1_decay_mode;  //!         
  Float_t         true_tau_1_et;  //!         
  Float_t         true_tau_1_eta;  //!         
  Float_t         true_tau_1_eta_vis;  //!         
  Double_t        true_tau_1_eta_vis_charged;  //!         
  Double_t        true_tau_1_eta_vis_neutral;  //!         
  Int_t           true_tau_1_isEle;  //!         
  Int_t           true_tau_1_isHadTau;  //!         
  Int_t           true_tau_1_isJet;  //!         
  Int_t           true_tau_1_isMuon;  //!         
  Int_t           true_tau_1_isTau;  //!         
  Float_t         true_tau_1_m;  //!         
  Float_t         true_tau_1_m_vis;  //!         
  Double_t        true_tau_1_m_vis_charged;  //!         
  Double_t        true_tau_1_m_vis_neutral;  //!         
  Int_t           true_tau_1_mother_pdgId;  //!         
  Int_t           true_tau_1_mother_status;  //!         
  Int_t           true_tau_1_n_charged;  //!         
  Int_t           true_tau_1_n_charged_pion;  //!         
  Int_t           true_tau_1_n_neutral;  //!         
  Int_t           true_tau_1_n_neutral_pion;  //!         
  Float_t         true_tau_1_neutral_eta;  //!         
  Float_t         true_tau_1_neutral_m;  //!         
  Float_t         true_tau_1_neutral_phi;  //!         
  Float_t         true_tau_1_neutral_pt;  //!         
  Float_t         true_tau_1_neutrino_eta;  //!         
  Float_t         true_tau_1_neutrino_m;  //!           
  Float_t         true_tau_1_neutrino_phi;//!         
  Float_t         true_tau_1_neutrino_pt;//!         
  Int_t           true_tau_1_origin;//!         
  Int_t           true_tau_1_pdgId;//!         
  Float_t         true_tau_1_phi;//!         
  Float_t         true_tau_1_phi_vis;//!         
  Double_t        true_tau_1_phi_vis_charged;//!         
  Double_t        true_tau_1_phi_vis_neutral;//!         
  Float_t         true_tau_1_pt;//!         
  Float_t         true_tau_1_pt_vis;//!         
  Double_t        true_tau_1_pt_vis_charged;//!         
  Double_t        true_tau_1_pt_vis_neutral;//!         
  Float_t         true_tau_1_pz;//!         
  Float_t         true_tau_1_q;//!         
  Int_t           true_tau_1_status;//!         
  Int_t           true_tau_1_type;//!         
  Float_t         truth_met_et;//!         
  Float_t         truth_met_etx;//!         
  Float_t         truth_met_ety;//!         
  Float_t         truth_met_phi;//!         
  Float_t         truth_met_sig;//!         
  Float_t         truth_met_sumet;//!         
  Double_t        weight_mc;//!         
  Double_t        weight_total;//!         

  // List of branches                                                                                                                                                                              
  TBranch        *b_n_avg_int;   //!                                          
  TBranch        *b_true_tau_1_classifierParticleOrigin;   //!                                                                                                                                        
  TBranch        *b_true_tau_1_classifierParticleType;   //!                                                                                                                                          
  TBranch        *b_true_tau_1_decay_mode;   //!                                                                                                                                                      
  TBranch        *b_true_tau_1_et;   //!                                                                                                                                                              
  TBranch        *b_true_tau_1_eta;   //!                                                                                                                                                             
  TBranch        *b_true_tau_1_eta_vis;   //!                                                                                                                                                         
  TBranch        *b_true_tau_1_eta_vis_charged;   //!                                                                                                                                                 
  TBranch        *b_true_tau_1_eta_vis_neutral;   //!                                                                                                                                                 
  TBranch        *b_true_tau_1_isEle;   //!   
  TBranch        *b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo;   //!                                                                                                                       
  TBranch        *b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I;   //!                                                                                              
  TBranch        *b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_03dR30_L1DR_TAU20ITAU12I_J25;   //!                                                                                          
  TBranch        *b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I;   //!                                                                                                     
  TBranch        *b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1DR_TAU20ITAU12I_J25;   //!                                                                                                 
  TBranch        *b_HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1TAU20IM_2TAU12IM;   //!                                                                                                    
  TBranch        *b_NOMINAL_pileup_combined_weight;   //!                                                                                                                                          
  TBranch        *b_NOMINAL_pileup_random_run_number;   //!                                                                                                                                        
  TBranch        *b_event_number;   //!                                                                                                                                                            
  TBranch        *b_jet_NOMINAL_central_jets_global_effSF_JVT;   //!                                                                                                                               
  TBranch        *b_jet_NOMINAL_central_jets_global_ineffSF_JVT;   //!                                                                                                                             
  TBranch        *b_jet_NOMINAL_forward_jets_global_effSF_JVT;   //!                                                                                                                               
  TBranch        *b_jet_NOMINAL_forward_jets_global_ineffSF_JVT;   //!                                                                                                                             
  TBranch        *b_jet_NOMINAL_global_effSF_MVX;   //!                                                                                                                                            
  TBranch        *b_jet_NOMINAL_global_ineffSF_MVX;   //!                                                                                                                                          
  TBranch        *b_n_actual_int;   //!                                                                                                     
  /*
  TBranch        *b_n_avg_int_cor;   //!                                                                                                                                                              
  TBranch        *b_n_bjets;                                                          
  TBranch        *b_n_electrons;   //!                                                                                                                                                             
  TBranch        *b_n_jets;   //!                                                                                                                                                                  
  TBranch        *b_n_jets_bad;   //!                                                                                                                                                              
  TBranch        *b_n_muons;   //!                                                                                                                                                                 
  TBranch        *b_n_photons;   //!                                                                                                                                                               
  TBranch        *b_n_pvx;   //!                                                                                                                                                                   
  TBranch        *b_n_taus;   //!                                                                                                                                                                  
  TBranch        *b_n_taus_loose;   //!                                                                                                                                                            
  TBranch        *b_n_taus_medium;   //!                                                                                                                                                           
  TBranch        *b_n_taus_tight;   //!                                                                                                                                                            
  TBranch        *b_n_truth_gluon_jets;   //!                                                                                                                                                      
  TBranch        *b_n_truth_jets;   //!                                                                                                                                                            
  TBranch        *b_n_truth_jets_pt20_eta45;   //!                                                                                                                                                 
  TBranch        *b_n_truth_quark_jets;   //!                                                                                                                                                      
  TBranch        *b_n_vx;   //!                                                                                                                                                                    
  TBranch        *b_run_number;   //!                                                                                                                                                              
  TBranch        *b_sherpa_n_sherpa_truth_jets_pt20_eta45;   //!                                                                                                                                   
  TBranch        *b_sherpa_weight;   //!                                                                                                                                                           


  

*/
  TBranch        *b_true_ditau;   //!                                                                                                                                                              
  TBranch        *b_true_ditau_CP_ip_tau0_x_ip;   //!                                                                                                                                              
  TBranch        *b_true_ditau_CP_ip_tau0_y_ip;   //!                                                                                                                                              
  TBranch        *b_true_ditau_CP_ip_tau0_z_ip;   //!                                                                                                                                              
  TBranch        *b_true_ditau_CP_ip_tau1_mag;   //!                                                                                                                                               
  TBranch        *b_true_ditau_CP_ip_tau1_x_ip;   //!                                                                                                                                              
  TBranch        *b_true_ditau_CP_ip_tau1_y_ip;   //!                                                                                                                                              
  TBranch        *b_true_ditau_CP_ip_tau1_z_ip;   //!                                                                                                                                              
  TBranch        *b_true_ditau_CP_phi_star_cp_a1_rho;   //!                                                                                                                                        
  TBranch        *b_true_ditau_CP_phi_star_cp_ip;   //!                                                                                                                                            
  TBranch        *b_true_ditau_CP_phi_star_cp_ip_rho;   //!                                                                                                                                        
  TBranch        *b_true_ditau_CP_phi_star_cp_ip_rho_opt;   //!                                                                                                                                    
  TBranch        *b_true_ditau_CP_phi_star_cp_rho_ip;   //!                                                                                                                                        
  TBranch        *b_true_ditau_CP_phi_star_cp_rho_rho;   //!                                                                                                                                       
  TBranch        *b_true_ditau_CP_tau0_upsilon;   //!                                                                                                                                              
  TBranch        *b_true_ditau_CP_tau1_upsilon;   //!                                                                                                                                              
  TBranch        *b_true_ditau_cosalpha;   //!                                                                                                                                                     
  TBranch        *b_true_ditau_deta;   //!                                                                                                                                                         
  TBranch        *b_true_ditau_dphi;   //!                                                                                                                                                         
  TBranch        *b_true_ditau_dpt;   //!                                                                                                                                                          
  TBranch        *b_true_ditau_dr;   //!                                                                                                                                                           
  TBranch        *b_true_ditau_mass;   //!                                                                                                                                                         
  TBranch        *b_true_ditau_ptx;   //!                                                                                                                                           
  TBranch        *b_true_ditau_pty;   //!                                                                                                                                                          
  TBranch        *b_true_ditau_qxq;   //!                                                                                                                                                          
  TBranch        *b_true_ditau_scal_sum_pt;   //!                                                                                                                                                  
  TBranch        *b_true_ditau_vect_sum_pt;   //!                                                                                                                                                  
  TBranch        *b_true_ditau_vis_cosalpha;   //!                                                                                                                                                 
  TBranch        *b_true_ditau_vis_deta;   //!                                                                                                                                                     
  TBranch        *b_true_ditau_vis_dphi;   //!                                                                                                                                                     
  TBranch        *b_true_ditau_vis_dr;   //!                                                                                                                                                       
  TBranch        *b_true_ditau_vis_mass;   //!                                                                                                                                                     
  TBranch        *b_true_ditau_vis_scal_sum_pt;   //!                                                                                                                                              
  TBranch        *b_true_ditau_vis_vect_sum_pt;   //!                                                                                                                                              
  TBranch        *b_true_tau_0;   //!                                                                                                                                                              
  TBranch        *b_true_tau_0_charged_eta;   //!                                                                                                                                                  
  TBranch        *b_true_tau_0_charged_m;   //!                                                                                                                                                    
  TBranch        *b_true_tau_0_charged_phi;   //!                                                                                                                                                  
  TBranch        *b_true_tau_0_charged_pt;   //!                                                                                                                                                   
  TBranch        *b_true_tau_0_classifierParticleOrigin;   //!                                                                                                                                     
  TBranch        *b_true_tau_0_classifierParticleType;   //!                                                                                                                                       
  TBranch        *b_true_tau_0_decay_mode;   //!                                                                                                                                                   
  TBranch        *b_true_tau_0_et;   //!                                                                                                                                                           
  TBranch        *b_true_tau_0_eta;   //!                                                                                                                                                          
  TBranch        *b_true_tau_0_eta_vis;   //!                                                                                                                                  
  TBranch        *b_true_tau_0_eta_vis_charged;   //!                                                                                                                                              
  TBranch        *b_true_tau_0_eta_vis_neutral;   //!                                                                                                                                              
  TBranch        *b_true_tau_0_isEle;   //!                                                                                                                                                        
  TBranch        *b_true_tau_0_isHadTau;   //!                                                                                                                                                     
  TBranch        *b_true_tau_0_isJet;   //!                                                                                                                                                        
  TBranch        *b_true_tau_0_isMuon;   //!                                                                                                                                                       
  TBranch        *b_true_tau_0_isTau;   //!                                                                                                                                                        
  TBranch        *b_true_tau_0_m;   //!                                                                                                                                                            
  TBranch        *b_true_tau_0_m_vis;   //!                                                                                                                                                        
  TBranch        *b_true_tau_0_m_vis_charged;   //!                                                                                                                                                
  TBranch        *b_true_tau_0_m_vis_neutral;   //!                                                                                                                                                
  TBranch        *b_true_tau_0_mother_pdgId;   //!                                                                                                                                                 
  TBranch        *b_true_tau_0_mother_status;   //!                                                                                                                                                
  TBranch        *b_true_tau_0_n_charged;   //!                                                                                                                                                    
  TBranch        *b_true_tau_0_n_charged_pion;   //!                                                                                                                                               
  TBranch        *b_true_tau_0_n_neutral;   //!                                                                                                                                                    
  TBranch        *b_true_tau_0_n_neutral_pion;   //!                                                                                                                                               
  TBranch        *b_true_tau_0_neutral_eta;   //!                                                                                                                                                  
  TBranch        *b_true_tau_0_neutral_m;   //!                                                                                                                                                    
  TBranch        *b_true_tau_0_neutral_phi;   //!                                                                                                                                                  
  TBranch        *b_true_tau_0_neutral_pt;   //!                                                                                                                                                   
  TBranch        *b_true_tau_0_neutrino_eta;   //!                                                                       
  TBranch        *b_true_tau_0_neutrino_m;   //!                                                                                                                                                   
  TBranch        *b_true_tau_0_neutrino_phi;   //!                                                                                                                                                 
  TBranch        *b_true_tau_0_neutrino_pt;   //!                                                                                                                                                  
  TBranch        *b_true_tau_0_origin;   //!                                                                                                                                                       
  TBranch        *b_true_tau_0_pdgId;   //!                                                                                                                                                        
  TBranch        *b_true_tau_0_phi;   //!                                                                                                                                                          
  TBranch        *b_true_tau_0_phi_vis;   //!                                                                                                                                                      
  TBranch        *b_true_tau_0_phi_vis_charged;   //!                                                                                                                                              
  TBranch        *b_true_tau_0_phi_vis_neutral;   //!                                                                                                                                              
  TBranch        *b_true_tau_0_pt;   //!                                                                                                                                                           
  TBranch        *b_true_tau_0_pt_vis;   //!                                                                                                                                                       
  TBranch        *b_true_tau_0_pt_vis_charged;   //!                                                                                                                                               
  TBranch        *b_true_tau_0_pt_vis_neutral;   //!                                                                                                                                               
  TBranch        *b_true_tau_0_pz;   //!                                                                                                                                                           
  TBranch        *b_true_tau_0_q;   //!                                                                                                                                                            
  TBranch        *b_true_tau_0_status;   //!                                                                                                                                                       
  TBranch        *b_true_tau_0_type;   //!                                                                                                                                                         
  TBranch        *b_true_tau_1;   //!                                                                                                                                                              
  TBranch        *b_true_tau_1_charged_eta;   //!                                                                                                                                                  
  TBranch        *b_true_tau_1_charged_m;   //!                                                                                                                                                    
  TBranch        *b_true_tau_1_charged_phi;   //!                                                                                                                                                  
  TBranch        *b_true_tau_1_charged_pt;   //!                                                                         
  TBranch        *b_true_tau_1_isHadTau;   //!                             
  TBranch        *b_true_tau_1_isJet;   //!                                                                                                                                                        
  TBranch        *b_true_tau_1_isMuon;   //!                                                                                                                                                       
  TBranch        *b_true_tau_1_isTau;   //!                                                                                                                                                        
  TBranch        *b_true_tau_1_m;   //!                                                                                                                                                            
  TBranch        *b_true_tau_1_m_vis;   //!                                                                                                                                                        
  TBranch        *b_true_tau_1_m_vis_charged;   //!                                                                                                                                                
  TBranch        *b_true_tau_1_m_vis_neutral;   //!                                                                                                                                                
  TBranch        *b_true_tau_1_mother_pdgId;   //!                                                                                                                                                 
  TBranch        *b_true_tau_1_mother_status;   //!                                                                                                                                                
  TBranch        *b_true_tau_1_n_charged;   //!                                                                                                                                                    
  TBranch        *b_true_tau_1_n_charged_pion;   //!                                                                                                                                               
  TBranch        *b_true_tau_1_n_neutral;   //!                                                                                                                                                    
  TBranch        *b_true_tau_1_n_neutral_pion;   //!                                                                                                                                               
  TBranch        *b_true_tau_1_neutral_eta;   //!                                                                                                                                                  
  TBranch        *b_true_tau_1_neutral_m;   //!                                                                                                                                                    
  TBranch        *b_true_tau_1_neutral_phi;   //!                                                                                                                                                  
  TBranch        *b_true_tau_1_neutral_pt;   //!                                                                                                                                                   
  TBranch        *b_true_tau_1_neutrino_eta;   //!                                                                                                                                                 
  TBranch        *b_true_tau_1_neutrino_m;   //!                                                                                                                                                   
  TBranch        *b_true_tau_1_neutrino_phi;   //!                                                                                                                                                 
  TBranch        *b_true_tau_1_neutrino_pt;   //!                                                                                                                                                  
  TBranch        *b_true_tau_1_origin;   //!                                                 
  TBranch        *b_true_tau_tau_1_pz;   //!                                                                                            
  TBranch        *b_true_tau_1_pdgId;   //!                                                                                                                                                           
  TBranch        *b_true_tau_1_phi;   //!                                                                                                                                                             
  TBranch        *b_true_tau_1_phi_vis;   //!                                                                                                                                                         
  TBranch        *b_true_tau_1_phi_vis_charged;   //!                                                                                                                                                 
  TBranch        *b_true_tau_1_phi_vis_neutral;   //!                                                                                                                                                 
  TBranch        *b_true_tau_1_pt;   //!                                                                                                                                                              
  TBranch        *b_true_tau_1_pt_vis;   //!                                                                                                                                                          
  TBranch        *b_true_tau_1_pt_vis_charged;   //!                                                                                                                                                  
  TBranch        *b_true_tau_1_pt_vis_neutral;   //!                                                                                                                                                  
  TBranch        *b_true_tau_1_pz;   //!                                                                                                                                                              
  TBranch        *b_true_tau_1_q;   //!                                                                                                                                                               
  TBranch        *b_true_tau_1_status;   //!                                                                                                                                                          
  TBranch        *b_true_tau_1_type;   //!                                                                                                                                                            
  TBranch        *b_truth_met_et;   //!                                                                                                                                                               
  TBranch        *b_truth_met_etx;   //!                                                                                                                                                              
  TBranch        *b_truth_met_ety;   //!                                                                                                                                                              
  TBranch        *b_truth_met_phi;   //!                                                                                                                                                              
  TBranch        *b_truth_met_sig;   //!                                                                                                                                                              
  TBranch        *b_truth_met_sumet;   //!                                                                                                                                                            
  TBranch        *b_weight_mc;   //!                                                                                                                                                                  
  TBranch        *b_weight_total;   //!                                                              


  // this is a standard constructor
  truth_study ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(truth_study, 1);
};

#endif

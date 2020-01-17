
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include <vector>
#include "math.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include <numeric>     
#include "coordinateTools.h"

using TMath::ATan;
using TMath::Exp;
//---------------------------------------------------------------------Global Vars

const int           asize=35;
const long int      bsize=700;
const int           MixingNumber=9;//+1
const double        PI = 3.14159265359;

//---------------------------------------------------------------------CUTS
const float PtTrihi     = 4.0; 
const float PtTrilow    = 1.0;
const float PtAsshi     = 4.0; 
const float PtAsslow    = 1.0;
const float Thetahigh   = 0.35; 
const float Thetalow    = 0.0;
const float EtaCut      = 2.4;
const float jetEtaCut   = 1.0;
const float Vlim        = 20.0;
const float VzSepLim    = 1.5;
const float jetPtCut    = 100;

const int     DeltaPBin =30; 
const int     DeltaEBin =30; 
const float   etabound  = 2;

//---------------------------------------------------------------------BOOL FUNCTIONS
bool F_eventpass(
    int     HLT_J100E5,
    float   jetVz,
    float   jetPt[asize],
    int     ngen
    ){
    if(HLT_J100E5 == 0)             return false;                          
    if(fabs(jetVz) > Vlim)          return false;
    if(ngen > 18)                   return false;
    float max_element =0;
    for(int i=0; i<asize; i++){
        float this_element = jetPt[i];
        if(this_element > max_element){
            max_element=jetPt[i];
        }
    }
    if(max_element < jetPtCut)        return false;
    return true;
}

bool F_jetpass(
    float   jetEta[asize], 
    float   jetPt[asize],   
    int     ijet
    ){
    if(fabs( jetEta[ijet] ) >jetEtaCut)   return false;
    if(fabs( jetPt[ijet]  ) <jetPtCut)    return false;
    return true;
}

bool F_trigpassROTATED(
    float               jetPt[asize],
    float               jetEta[asize],
    float               jetPhi[asize],
    std::vector<int>    *partChg, 
    std::vector<float>  *partEta,
    std::vector<float>  *partPt,
    std::vector<float>  *partPhi,
    long int            XXtrk,
    int                 ijet
    ){
    if( 2*ATan(Exp(-(etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk])))) > Thetahigh)        return false;
    if( 2*ATan(Exp(-(etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk])))) < Thetalow)       return false;
    if(ptWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk]) >PtTrihi)                           return false;
    if(ptWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk]) >PtTrihi)                           return false;
    if((*partChg)[XXtrk] == 0)                                                                                                                  return false;
    if((etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk])) >EtaCut)                           return false;
    return true;
}

bool F_asspassROTATED(
    float               jetPt[asize],
    float               jetEta[asize],
    float               jetPhi[asize],
    std::vector<int>    *partChg, 
    std::vector<float>  *partEta,
    std::vector<float>  *partPt,
    std::vector<float>  *partPhi,
    long int            XXtrk,
    int                 ijet
    ){
    if( 2*ATan(Exp(-(etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk])))) > Thetahigh)        return false;
    if( 2*ATan(Exp(-(etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk])))) < Thetalow)       return false;
    if(ptWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk]) >PtAsshi)                           return false;
    if(ptWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk]) >PtAsshi)                           return false;
    if((*partChg)[XXtrk] == 0)                                                                                                                  return false;
    if((etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk])) >EtaCut)                           return false;
    return true;
}

bool F_trigpassBEAM(
    std::vector<int>    *partChg, 
    std::vector<float>  *partEta,
    std::vector<float>  *partPt,
    long int            XXtrk 
    ){
    //if( 2*ATan(Exp(-(trkEta[XXtrk]))) >Thetahi)         return false;
    //if( 2*ATan(Exp(-(trkEta[XXtrk]))) <Thetalow)        return false;
    if((*partPt)[XXtrk] >PtTrihi)                          return false;
    if((*partPt)[XXtrk] <PtTrilow)                         return false;
    if((*partChg)[XXtrk] == 0)                             return false;
    if((*partEta)[XXtrk] > EtaCut)                         return false;
    return true;
}

bool F_asspassBEAM(
    std::vector<int>    *partChg, 
    std::vector<float>  *partEta,
    std::vector<float>  *partPt,
    long int            YYtrk 
    ){
    //if( 2*ATan(Exp(-(trkEta[YYtrk]))) >Thetahi)         return false;
    //if( 2*ATan(Exp(-(trkEta[YYtrk]))) <Thetalow)        return false;
    if((*partPt)[YYtrk] >PtAsshi)                          return false;
    if((*partPt)[YYtrk] <PtAsslow)                         return false;
    if((*partChg)[YYtrk] == 0)                             return false;
    if((*partEta)[YYtrk] > EtaCut)                         return false;
    return true;
}
//-------------------------------------------------------------------MAIN_CODE-------------------------------------------------------------------------------
void PythiaRegCorr(){ 

    TH1::SetDefaultSumw2(kTRUE);


    TH2D* hSigBEAM = new TH2D("2Signal BEAM","2Signal BEAM",  DeltaPBin, -etabound, etabound, DeltaEBin, -.5*PI, 1.5*PI);
    TH2D* hMixBEAM = new TH2D("2Mix BEAM","2Mix BEAM",        DeltaPBin, -etabound, etabound, DeltaEBin, -.5*PI, 1.5*PI);
    TH2D* hRatBEAM = new TH2D("2R BEAM","2R BEAM",            DeltaPBin, -etabound, etabound, DeltaEBin, -.5*PI, 1.5*PI);

    //TH2D* hSigROTATED = new TH2D("2Signal BEAM","2Signal BEAM",  DeltaPBin, -etabound, etabound, DeltaEBin, -.5*PI, 1.5*PI);
    //TH2D* hMixROTATED = new TH2D("2Mix BEAM","2Mix BEAM",        DeltaPBin, -etabound, etabound, DeltaEBin, -.5*PI, 1.5*PI);
    //TH2D* hRatROTATED = new TH2D("2R BEAM","2R BEAM",            DeltaPBin, -etabound, etabound, DeltaEBin, -.5*PI, 1.5*PI);

    int HLT_J100E5=0;                       int BHLT_J100E5=0;
    int ngen =     0;                       int Bngen =     0;
    float jetVz =  0;                       float BjetVz =  0;

    float jetEta[asize] = {0};              float BjetEta[asize] = {0};
    float jetPt[asize] =  {0};              float BjetPt[asize] =  {0};
    float jetPhi[asize] = {0};              float BjetPhi[asize] = {0};

    std::vector<int>    * partChg=0;        std::vector<int>    * BpartChg=0;
    std::vector<float>  * partPt=0;         std::vector<float>  * BpartPt=0;
    std::vector<float>  * partEta=0;        std::vector<float>  * BpartEta=0;
    std::vector<float>  * partPhi=0;        std::vector<float>  * BpartPhi=0;
    

    TFile* f1 = new TFile("/storage1/users/ptg2/PythiaMini.root");
    TFile* f2 = new TFile("/storage1/users/ptg2/PythiaMini.root");

    TTree* jetGEN =      (TTree*) f1->Get("ak4PFJetAnalyzer/t" );                   TTree* BjetGEN =      (TTree*) f2->Get("ak4PFJetAnalyzer/t" );
    TTree* partGEN =     (TTree*) f1->Get("HiGenParticleAna/hi");                   TTree* BpartGEN =     (TTree*) f2->Get("HiGenParticleAna/hi");
    TTree* hltanalysis = (TTree*) f1->Get("hltanalysis/HltTree");                   TTree* Bhltanalysis = (TTree*) f2->Get("hltanalysis/HltTree");

//-------------------------------------------------------------------JET_TREE
    jetGEN->SetBranchStatus("*",0);                                                 BjetGEN->SetBranchStatus("*",0);            

    jetGEN->SetBranchStatus("geneta",1);                                            BjetGEN->SetBranchStatus("geneta",1);       
    jetGEN->SetBranchStatus("genpt",1);                                             BjetGEN->SetBranchStatus("genpt",1);        
    jetGEN->SetBranchStatus("genphi",1);                                            BjetGEN->SetBranchStatus("genphi",1);       
    jetGEN->SetBranchStatus("vz"    ,1);                                            BjetGEN->SetBranchStatus("vz"    ,1);
    jetGEN->SetBranchStatus("ngen"  ,1);                                            BjetGEN->SetBranchStatus("ngen"  ,1);

    jetGEN->SetBranchAddress("geneta",jetEta);                                      BjetGEN->SetBranchAddress("geneta", BjetEta);       
    jetGEN->SetBranchAddress("genpt",jetPt);                                        BjetGEN->SetBranchAddress("genpt",  BjetPt);        
    jetGEN->SetBranchAddress("genphi",jetPhi);                                      BjetGEN->SetBranchAddress("genphi", BjetPhi);       
    jetGEN->SetBranchAddress("vz",&jetVz);                                          BjetGEN->SetBranchAddress("vz",    &BjetVz);
    jetGEN->SetBranchAddress("ngen",&ngen);                                         BjetGEN->SetBranchAddress("ngen",  &Bngen);

//-------------------------------------------------------------------PARTICLE_TREE
    partGEN->SetBranchStatus("*",0);                                                BpartGEN->SetBranchStatus("*",0);

    partGEN->SetBranchStatus("chg",1);                                              BpartGEN->SetBranchStatus("chg",1);
    partGEN->SetBranchStatus("phi",1);                                              BpartGEN->SetBranchStatus("phi",1);
    partGEN->SetBranchStatus("eta",1);                                              BpartGEN->SetBranchStatus("eta",1);
    partGEN->SetBranchStatus("pt",1);                                               BpartGEN->SetBranchStatus("pt",1);

    partGEN->SetBranchAddress("chg",&partChg);                                      BpartGEN->SetBranchAddress("chg", &BpartChg);
    partGEN->SetBranchAddress("phi",&partPhi);                                      BpartGEN->SetBranchAddress("phi", &BpartPhi);
    partGEN->SetBranchAddress("eta",&partEta);                                      BpartGEN->SetBranchAddress("eta", &BpartEta);
    partGEN->SetBranchAddress("pt",&partPt);                                        BpartGEN->SetBranchAddress("pt",  &BpartPt);

//-------------------------------------------------------------------HLT_TREE
    hltanalysis->SetBranchStatus("*",0);                                            Bhltanalysis->SetBranchStatus("*",0);
    hltanalysis->SetBranchStatus("HLT_HIAK4CaloJet100_v1",1);                       Bhltanalysis->SetBranchStatus("HLT_HIAK4CaloJet100_v1",1);
    hltanalysis->SetBranchAddress("HLT_HIAK4CaloJet100_v1",&HLT_J100E5);            Bhltanalysis->SetBranchAddress("HLT_HIAK4CaloJet100_v1",&BHLT_J100E5);



    //*************** CHANGE BEFORE COMPILING ***************
    //*************** CHANGE BEFORE COMPILING ***************
    //
    const long int ALLEVENTS = jetGEN->GetEntries();
    long int nevents =5000;
    long int progressmeter1 = floor(nevents/50);
    long int progressmeter2 = floor(nevents/50);
    //
    //***************
    //***************


//-------------------------------------------------------------------MAIN_LOOPS
            for(long int  ievent=0; ievent<nevents; ievent++){
                long int percentdone =floor(100*ievent/nevents);
                if(ievent%progressmeter1==0) cout<< percentdone << "percent done for Signal " <<endl;
                
                jetGEN->GetEntry(ievent);
                hltanalysis->GetEntry(ievent);
                partGEN->GetEntry( ievent );

                if(!F_eventpass(HLT_J100E5,jetVz,jetPt, ngen)) continue;

                long int NNtrk = partPt->size();
                long int Ntrig  =0;

for(int ijet=0; ijet<ngen; ijet++){
if( !F_jetpass(jetEta, jetPt, ijet) ) continue;
            
                for(long int  XXtrk=0; XXtrk < NNtrk; XXtrk++ ){
                    if(!F_trigpassROTATED(jetPt, jetEta, jetPhi, partChg, partEta, partPt, partPhi, XXtrk, ijet)) continue;

                    Ntrig=Ntrig+1;
                }

                for(long int  XXtrk=0;XXtrk < NNtrk; XXtrk++ ){
                    if(!F_trigpassROTATED(jetPt, jetEta, jetPhi, partChg, partEta, partPt, partPhi, XXtrk, ijet)) continue;

float THISeventEta = etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk]);
float THISeventPhi = phiWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk]);
//                    float THISeventEta = (*partEta)[XXtrk];
//                    float THISeventPhi = (*partPhi)[XXtrk];

                    for(long int  YYtrk=0; YYtrk< NNtrk; YYtrk++ ){
                        if(!F_asspassROTATED(jetPt, jetEta, jetPhi, partChg, partEta, partPt, partPhi, YYtrk, ijet)) continue;

float YYTHISeventEta = etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[YYtrk],(*partEta)[YYtrk],(*partPhi)[YYtrk]);
float YYTHISeventPhi = phiWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[YYtrk],(*partEta)[YYtrk],(*partPhi)[YYtrk]);

                        float deltaEta = THISeventEta - YYTHISeventEta;
                        float deltaPhi = THISeventPhi - YYTHISeventPhi;

                        if(deltaPhi<0) deltaPhi+=2*PI;
                        if(deltaPhi>1.5*PI) deltaPhi -= 2*PI;

                        hSigBEAM->Fill(deltaEta, deltaPhi, 1./Ntrig);

                        //float negDeltaEta = -1* deltaEta;
                        //float negDeltaPhi = -1* deltaPhi;
                        //if(negDeltaPhi<-0.5*PI) negDeltaPhi += 2*PI;
                        //hSigBEAM->Fill(negDeltaEta, negDeltaPhi, 1./Ntrig);
                    }
                }
            }
            }

long int RegPass=0;
long int ConePass=0;

            for(long int  ievent=0; ievent<nevents; ievent++){
                int percentdone =floor(100 * ievent/nevents);
                if(ievent%progressmeter2==0) cout<< percentdone << "percent done for Background " <<endl;

                jetGEN->GetEntry(ievent);
                hltanalysis->GetEntry(ievent);
                partGEN->GetEntry(ievent);      
                
                if(!F_eventpass(HLT_J100E5,jetVz,jetPt, ngen)) continue;
   
                long int NNtrk = partPt->size();

                std::vector<long int> numbers;

                for(long int imix=ievent+1; imix < ALLEVENTS; imix++){ 
                    if(imix == ievent)                          continue;
                    if((imix-ievent)%1000==0 && imix > 10) cout<< "IMIX is deep: " << imix  <<endl;
                    if(numbers.size() > MixingNumber)                     break;

                    BpartGEN->GetEntry(imix);
                    Bhltanalysis->GetEntry(imix);
                    BjetGEN->GetEntry(imix);
                    if(!F_eventpass(BHLT_J100E5,BjetVz,BjetPt, Bngen))    continue;                    
                    if(fabs(jetVz-BjetVz) > VzSepLim )                    continue;
                    numbers.push_back (imix);
                }
                for(int qmix = 0; qmix <  numbers.size(); qmix++){
                    long int jmix = numbers[qmix];
                    long int Ntrig = 0;

for(int ijet=0; ijet<ngen; ijet++){
if( !F_jetpass(jetEta, jetPt, ijet) ) continue;

                    for(long int  XXtrk=0;  XXtrk < NNtrk; XXtrk++){
                        if(!F_trigpassROTATED(jetPt, jetEta, jetPhi, partChg, partEta, partPt, partPhi, XXtrk, ijet)) continue;
                        Ntrig=Ntrig+1;
                    }

                    BjetGEN->GetEntry(jmix);
                    BpartGEN->GetEntry(jmix);

                    for(long int  XXtrk=0;  XXtrk < NNtrk; XXtrk++){
                        if(!F_trigpassROTATED(jetPt, jetEta, jetPhi, partChg, partEta, partPt, partPhi, XXtrk, ijet)) continue;

float THISeventEta = etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk]);
float THISeventPhi = phiWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*partPt)[XXtrk],(*partEta)[XXtrk],(*partPhi)[XXtrk]);

//                        float THISeventEta = (*partEta)[XXtrk];
//                        float THISeventPhi = (*partPhi)[XXtrk];

                        long int BNNtrk = BpartPt->size();

//for(int KKjet=0; KKjet<Bngen; KKjet++){
if( !F_jetpass(BjetEta, BjetPt, KKjet) ) continue;

                        for(long int  YYtrk=0;  YYtrk < BNNtrk; YYtrk++){                      
                            if(!F_asspassROTATED(jetPt, jetEta, jetPhi, BpartChg, BpartEta, BpartPt, BpartPhi, YYtrk, ijet)) continue;
                            //RegPass +=1; 
                            //if( 2*ATan(Exp(-(etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*BpartPt)[YYtrk],(*BpartEta)[YYtrk],(*BpartPhi)[YYtrk])))) > Thetahigh ) continue;
                            //ConePass +=1;

float YYMIXeventEta = etaWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*BpartPt)[YYtrk],(*BpartEta)[YYtrk],(*BpartPhi)[YYtrk]);
float YYMIXeventPhi = phiWRTJet(jetPt[ijet],jetEta[ijet],jetPhi[ijet],(*BpartPt)[YYtrk],(*BpartEta)[YYtrk],(*BpartPhi)[YYtrk]);

                            float deltaEta = THISeventEta - YYMIXeventEta;
                            float deltaPhi = THISeventPhi - YYMIXeventPhi;

                            if(deltaPhi<0) deltaPhi+=2*PI;
                            if(deltaPhi>1.5*PI) deltaPhi -= 2*PI;

                            hMixBEAM->Fill(deltaEta, deltaPhi, (1./Ntrig));
                        }


//                    }
                    }
                }
                }
            }


cout << RegPass  << endl;
cout << ConePass << endl;

            long int PSIG = 0;
            long int PMIX = 0;
            PMIX=hMixBEAM->GetEntries();
            PSIG=hSigBEAM->GetEntries();

            cout << " ********** " << endl;
            cout << " ********** " << endl;

            cout << " Outputs: " << endl;

            cout << " ********** " << endl;
            cout << " ********** " << endl;

            cout << (PSIG) << " Signal fills " << endl;
            cout << PMIX << " Background fills " << endl;

            TCanvas* c1 = new TCanvas("c1", "", 1000, 1000);
            c1->cd();
            hSigBEAM->GetYaxis()->SetTitle("Delta Phi");
            hSigBEAM->GetXaxis()->SetTitle("Delta Eta");
            hSigBEAM->Draw("SURF2");

            TCanvas* c2 = new TCanvas("c2", "", 1000, 1000);
            c2->cd();
            hMixBEAM->GetYaxis()->SetTitle("Delta Phi");
            hMixBEAM->GetXaxis()->SetTitle("Delta Eta");
            hMixBEAM->Draw("SURF2");

            TCanvas* c3 = new TCanvas("c3", "", 1000, 1000);
            hRatBEAM = (TH2D*)hSigBEAM->Clone();
            hRatBEAM->SetName("hratio");
            hRatBEAM->SetTitle("hratio");
            hRatBEAM->Divide(hMixBEAM);
            c3->cd();
            hRatBEAM->GetYaxis()->SetTitle("Delta Phi");
            hRatBEAM->GetXaxis()->SetTitle("Delta Eta");
            hRatBEAM->GetZaxis()->SetRangeUser(0, .5);
            hRatBEAM->Draw("SURF2");

}

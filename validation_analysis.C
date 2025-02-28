#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include <string>
#include <iostream>
#include <map>
#include <set>

std::map<std::string, int> runMap = {
        {"thr_scan", 560692},
        {"dig_scan", 560721},
        {"vresetd_1d", 560724},
        {"vresetd_2d", 560725},
        {"ps_1d", 560726},
        {"ps_2d", 560727},
        {"tot_full", 560728},
        {"vcasn", 560729},
        {"ithr", 560730}
};

// Tree variables
short int chipid[1024], row[1024], thr[1024], nhits[1024], vresetd[1024], strobedel[1024];
float  noise[1024];
bool success[1024];
unsigned char charge[1024];
const int nVar = 6; // max number of variables to plot
const int nScan = 9; // number of scans 
int nRandom = 1000; // number of random entries to be picked up from each tree

//histograms
std::map<std::string, std::vector<TH1F*>> hHisto = {
    {"thr_scan", {}},
    {"dig_scan", {}},
    {"vresetd_1d", {}},
    {"vresetd_2d", {}},
    {"ps_1d", {}},
    {"ps_2d", {}},
    {"tot_full", {}},
    {"vcasn", {}},
    {"ithr", {}}
};

std::map<std::string, std::vector<TH1F*>> hHistoRef = {
    {"thr_scan", {}},
    {"dig_scan", {}},
    {"vresetd_1d", {}},
    {"vresetd_2d", {}},
    {"ps_1d", {}},
    {"ps_2d", {}},
    {"tot_full", {}},
    {"vcasn", {}},
    {"ithr", {}}
};

std::map<std::string, std::vector<string>> histoX = {
    {"thr_scan", {"chipid","row","thr","noise","success","empty1"}},
    {"dig_scan", {"chipid","row","n_hits","empty1","empty2","empty3"}},
    {"vresetd_1d", {"chipid","row","n_hits","vresetd","empty1","empty2"}},
    {"vresetd_2d", {"chipid","row","thr","noise","vresetd","success"}},
    {"ps_1d", {"chipid","row","n_hits","strobedel","empty1","empty2"}},
    {"ps_2d", {"chipid","row","n_hits","strobedel","charge","empty1"}},
    {"tot_full", {"chipid","row","n_hits","strobedel","charge","empty1"}},
    {"vcasn", {"chipid","row","vcasn","noise","success","empty1"}},
    {"ithr", {"chipid","row","ithr","noise","success","empty1"}}
};
std::map<std::string, std::vector<int>> histoBins = {
    {"thr_scan", {6700, 512, 100, 50, 2, 1}},
    {"dig_scan", {6700, 512, 100, 1, 1, 1}},
    {"vresetd_1d", {6700, 512, 100, 255, 1, 1}},
    {"vresetd_2d", {6700, 512, 100, 50, 255, 2}},
    {"ps_1d", {6700, 512, 100, 45, 1, 1}},
    {"ps_2d", {6700, 512, 100, 500, 170,1}},
    {"tot_full", {6700, 512, 100, 250, 170, 1}},
    {"vcasn",{6700, 512, 30, 50, 2, 1}},
    {"ithr", {6700, 512, 50, 50, 2, 1}}
};
std::map<std::string, std::vector<double>> histoMin = {
    {"thr_scan", {0, 0, 10, 0, 0, 0}},
    {"dig_scan", {0, 0, 0, 0, 0, 0}},
    {"vresetd_1d", {0, 0, 0, 0, 0, 0}},
    {"vresetd_2d", {0, 0, 0, 0, 0, 0}},
    {"ps_1d", {0, 0, 0, 0, 0, 0}},
    {"ps_2d", {0, 0, 0, 0, 0, 0}},
    {"tot_full", {0, 0, 0, 0, 0, 0}},
    {"vcasn", {0, 0, 20, 0, 0, 0}},
    {"ithr", {0, 0, 10, 0, 0, 0}}
};
std::map<std::string, std::vector<double>> histoMax = {
    {"thr_scan", {6700, 512, 500, 100, 2, 1}},
    {"dig_scan", {6700, 512, 300,1, 1, 1}},
    {"vresetd_1d", {6700, 512, 300, 255, 1, 1}},
    {"vresetd_2d", {6700, 512, 500, 100, 255, 2}},
    {"ps_1d", {6700, 512, 300, 450, 1, 1}},
    {"ps_2d", {6700, 512, 300, 2100, 170, 1}},
    {"tot_full", {6700, 512, 300, 1200, 170, 1}},
    {"vcasn", {6700, 512, 80, 100, 2, 1}},
    {"ithr", {6700, 512, 120, 100, 2, 1}}
};

// random tree entries
std::set<int> entries;


void initHisto();
void fillHistograms(std::string runType, TTree* myTree, bool isRef);
TTree* getTree(std::string runType, bool isRef, std::string pathRef, int fileNumber);
std::set<int> generateRandomNum(int n, int a, int b);

// MAIN
void validation_analysis(std::string pathRef){

    // pathRef is the base path to know where the reference Trees are 
    // Trees to compare to reference tree are supposed to be in the local folder (where O2 was run)

    // init histo 
    initHisto();

    for (const auto &m: runMap){ // loop on runs or runtypes
        std::cout<<"\n\n==> Processing "<<m.first<<std::endl;
        // Get data and fill histograms
        for (int fileNumber = 0; fileNumber < 3; fileNumber++){ //loop on root files 
            // Get reference data 
            std::cout<<"... Getting tree data of ref run - file "<<fileNumber<<std::endl;
            TTree *dataTreeRef = getTree(m.first, true, pathRef, fileNumber); //fileNumber can be 0,1,2 (3 files per run)
            if(dataTreeRef->GetEntries()==0){
                std::cout<<"WARN: tree is empty"<<std::endl;
            }
            std::cout<<"... Filling histograms"<<std::endl;
            if (dataTreeRef->GetEntries()-1 < nRandom){
                std::cout<<"WARN: too many random entries. Setting nRandom = (tree entries - 1) = "<<(dataTreeRef->GetEntries()-1)<<std::endl;
                nRandom = dataTreeRef->GetEntries()-1;
            }
            entries = generateRandomNum(nRandom, 0, dataTreeRef->GetEntries()-1); //random entries
            fillHistograms(m.first, dataTreeRef, 1);
            
            std::cout<<"... Getting tree data of current run - file "<<fileNumber<<std::endl;
            TTree *dataTree = getTree(m.first, false, pathRef, fileNumber); //fileNumber can be 0,1,2 (3 files per run)
            if(dataTree->GetEntries()==0){
                std::cout<<"WARN: tree is empty"<<std::endl;
            }
            std::cout<<"... Filling histograms"<<std::endl;
            fillHistograms(m.first, dataTree, 0);
        }

        // Make all ratios with Ref runs and print all in PDF file 
        std::cout<<"... Making ratios and storing plots"<<std::endl;
        for(size_t ihist = 0; ihist<hHisto[m.first].size(); ihist++){
            TH1F *hRatio = (TH1F*)hHisto[m.first][ihist]->Clone("hRatio");
            hRatio->SetStats(0);
            hRatio->SetTitle(Form("Ratio THIS run / Ref run -> %s -> %s; %s; THIS / REF", m.first.c_str(), histoX[m.first][ihist].c_str(), histoX[m.first][ihist].c_str()));
            // ratio 
            hRatio->Divide(hHistoRef[m.first][ihist]);
            // style 
            hRatio->SetLineColor(kRed);
            hRatio->SetMarkerColor(kRed);
            hRatio->SetLineWidth(2);
            // draw 
            TCanvas cnv("cnv","cnv",1000,200); 
            cnv.SetGridx();
            cnv.SetGridy();
            hRatio->Draw("HIST");
            hRatio->GetYaxis()->SetRangeUser(0.996,1.004);
            if(m.first=="dig_scan" && ihist==0) {
                cnv.SaveAs("validation.pdf[");
            }
            cnv.SaveAs("validation.pdf");
            if(m.first=="vresetd_2d" && ihist==hHisto[m.first].size()-1){
                cnv.SaveAs("validation.pdf]");
            }
            delete hRatio;
        }
        std::cout<<"... DONE"<<std::endl;
    }

}

std::set<int> generateRandomNum(int n, int a, int b){
    TRandom3 randGen(0);
    std::set<int> unique_numbers;

    if (n==b){
        for(int i=a; i<=n; i++){
            unique_numbers.insert(i);
        }
        return unique_numbers;
    }
    
    while ((int)unique_numbers.size() < n) {
        int rnd = randGen.Integer(b - a + 1) + a; // Generate number in [a, b]
        unique_numbers.insert(rnd);
    }
    return unique_numbers;
} 

void fillHistograms(std::string runType, TTree* myTree, bool isRef) {

    for(auto i : entries){ // take random entries from the tree
        myTree->GetEntry(i);
        for(size_t iv = 0; iv<histoX[runType].size(); iv++){
            for (int ic=0; ic<1024; ic++){
                if(histoX[runType][iv]=="chipid"){
                    isRef ? hHistoRef[runType][iv]->Fill(chipid[ic]) : hHisto[runType][iv]->Fill(chipid[ic]);
                } else if (histoX[runType][iv]=="row"){
                    isRef ? hHistoRef[runType][iv]->Fill(row[ic]) : hHisto[runType][iv]->Fill(row[ic]);
                } else if (histoX[runType][iv]=="n_hits"){
                    isRef ? hHistoRef[runType][iv]->Fill(nhits[ic]) : hHisto[runType][iv]->Fill(nhits[ic]);
                } else if (histoX[runType][iv]=="thr"){
                    isRef ? hHistoRef[runType][iv]->Fill(thr[ic]) : hHisto[runType][iv]->Fill(thr[ic]);
                } else if (histoX[runType][iv]=="noise"){
                    isRef ? hHistoRef[runType][iv]->Fill(noise[ic]) : hHisto[runType][iv]->Fill(noise[ic]);
                } else if (histoX[runType][iv]=="vresetd"){
                    isRef ? hHistoRef[runType][iv]->Fill(vresetd[ic]) : hHisto[runType][iv]->Fill(vresetd[ic]);
                } else if (histoX[runType][iv]=="success"){
                    isRef ? hHistoRef[runType][iv]->Fill((int)success[ic]) : hHisto[runType][iv]->Fill((int)success[ic]);
                } else if (histoX[runType][iv]=="strobedel"){
                    isRef ? hHistoRef[runType][iv]->Fill(strobedel[ic]) : hHisto[runType][iv]->Fill(strobedel[ic]);
                } else if (histoX[runType][iv]=="charge"){
                    isRef ? hHistoRef[runType][iv]->Fill((int)charge[ic]) : hHisto[runType][iv]->Fill((int)charge[ic]);
                } else { //empty
                    continue;
                }
            }
        }
    }
}



TTree* getTree(std::string runType, bool isRef, std::string pathRef, int fileNumber){

    // run type can be: thr_scan, dig_scan, vresetd_1d, vresetd_2d, ps_1d, ps_2d, tot_full, vcasn, ithr 

    TFile *inFile;
    if(isRef) {
        inFile = TFile::Open(Form("%s/%d_%s/unknown_%d/%d_0_flpits2_modSel%d.root", pathRef.c_str(), runMap[runType],runType.c_str(),runMap[runType],runMap[runType],fileNumber));
    } else {
        inFile = TFile::Open(Form("unknown_%d/%d_0_flpits2_modSel%d.root",runMap[runType],runMap[runType],fileNumber));
    }
    TTree *myTree = (TTree*)inFile->Get("ITS_calib_tree");
    myTree->SetBranchAddress("chipid",&chipid[0]);
    myTree->SetBranchAddress("row",&row[0]);
    if (runType == "thr_scan"){
        myTree->SetBranchAddress("thr",&thr[0]);
        myTree->SetBranchAddress("noise",&noise[0]);
        myTree->SetBranchAddress("success",&success[0]);
    } else if (runType == "dig_scan"){
        myTree->SetBranchAddress("n_hits",&nhits[0]);
    } else if (runType == "vresetd_1d") {
        myTree->SetBranchAddress("n_hits",&nhits[0]);
        myTree->SetBranchAddress("vresetd",&vresetd[0]);
    } else if (runType == "vresetd_2d") {
        myTree->SetBranchAddress("thr",&thr[0]);
        myTree->SetBranchAddress("noise",&noise[0]);
        myTree->SetBranchAddress("success",&success[0]);
        myTree->SetBranchAddress("vresetd",&vresetd[0]);
    } else if (runType == "ps_1d") {
        myTree->SetBranchAddress("n_hits",&nhits[0]);
        myTree->SetBranchAddress("strobedel",&strobedel[0]);
    } else if (runType == "ps_2d" || runType == "tot_full") {
        myTree->SetBranchAddress("n_hits",&nhits[0]);
        myTree->SetBranchAddress("strobedel",&strobedel[0]);
        myTree->SetBranchAddress("charge",&charge[0]);
    }

    return myTree;
}

void initHisto(){
    for (auto m: runMap){
        for (int var=0; var<nVar; var++){
            std::string histoname = Form("%s_%s",m.first.c_str(),histoX[m.first][var].c_str());
            std::string histonameRef = histoname + "_ref";
            std::string histoaxes = Form("%s -> %s; %s; Counts", m.first.c_str(),histoX[m.first][var].c_str(),histoX[m.first][var].c_str());
            TH1F* myHist = new TH1F(histoname.c_str(), histoaxes.c_str(), histoBins[m.first][var], histoMin[m.first][var], histoMax[m.first][var]);
            TH1F* myHistRef = new TH1F(histonameRef.c_str(), histoaxes.c_str(), histoBins[m.first][var], histoMin[m.first][var], histoMax[m.first][var]);
            //fill all bins with 1 entry (help interpreting the ratios)
            for (int ibin=1; ibin<=histoBins[m.first][var]; ibin++){
                myHist->SetBinContent(ibin,1);
                myHistRef->SetBinContent(ibin,1);
            }
            hHisto[m.first].push_back(myHist);
            hHistoRef[m.first].push_back(myHistRef);
        }
    }
}
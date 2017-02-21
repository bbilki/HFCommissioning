#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TRandom3.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TColor.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TClonesArray.h"

#include <sys/types.h>
#include <dirent.h>

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TVector3.h>
#include <TMultiGraph.h>
#include <TH2Poly.h>
#include <TLine.h>

using namespace std;using namespace ROOT::Math;

struct edata
{
	vector <Int_t> *ieta;
	vector <Int_t> *iphi;
	vector <Int_t> *depth;
	vector <vector <Int_t>> *pulse;
	vector <vector <Int_t>> *tdc;
	vector <vector <Int_t>> *capid;
	vector <vector <Int_t>> *histo;
};
edata ed;

struct semapin
{
	int crate;
	int slot;
	int fiber;
	int channel;
	string PMTname;
	string Channelname;
	string winchester;
	int qiecard;
	string qiecrate;
	int qieslot;
	int PMTid;
	int BBid;
	int ieta;
	int iphi;
	int depth;
	int box;
	string boxname;
	string boxbarcode;
	int VA;
	int VB;
	int VC;
	
	float SX5gain;
	float SX5gainerr;
	float B904pedestal;
	float B904pedestalsigma;
};
semapin semin;

struct semap
{
	int crate;
	int slot;
	int fiber;
	int channel;
	string PMTname;
	string Channelname;
	string winchester;
	int qiecard;
	int ieta;
	int iphi;
	int depth;
	int box;
	string boxname;
	
	float SX5gain;
	float SX5gainerr;
	float B904pedestal;
};
semap sem;vector <semap> SEM;

struct pedlist
{
	int mapind;
	float pedmean;
	float pedsigma;
	TH1F* h;
};

struct ledlist
{
	int mapind;
	float Qmean;
	float Qsigma;
	float npemean;
	float npeerr;
	float gain;
	float gainerr;
	float nbadcapid;
	TH1F* p;//pulse
	TH1F* pn;//pulse norm
	TH1F* t;//tdc
	TH1F* q;//charge
};

//from JM 070516

float adc2fC_QIE10[256]={

  // =========== RANGE 0 ===========

  // --------- subrange 1 ---------
  -14.45,-11.35,-8.25,-5.15,-2.05,1.05,4.15,7.25,10.35,13.45,16.55,19.65,22.75,25.85,28.95,32.05,
  // --------- subrange 2 ---------
  36.7,42.9,49.1,55.3,61.5,67.7,73.9,80.1,86.3,92.5,98.7,104.9,111.1,117.3,123.5,129.7,135.9,142.1,148.3,154.5,
  // --------- subrange 3 ---------
  163.8,176.2,188.6,201.0,213.4,225.8,238.2,250.6,263.0,275.4,287.8,300.2,312.6,325.0,337.4,349.8,362.2,374.6,387.0,399.4,411.8,
  // --------- subrange 4 ---------
  430.4,455.2,480.0,504.8,529.6,554.4,579.2,
  // =========== RANGE 1 ===========

  // --------- subrange 1 ---------
  529.4,554.2,579.0,603.8,628.6,653.4,678.2,703.0,727.8,752.6,777.4,802.2,827.0,851.8,876.6,901.4,
  // --------- subrange 2 ---------
  938.6,988.2,1037.8,1087.4,1137.0,1186.6,1236.2,1285.8,1335.4,1385.0,1434.6,1484.2,1533.8,1583.4,1633.0,1682.6,1732.2,1781.8,1831.4,1881.0,
  // --------- subrange 3 ---------
  1955.4,2054.6,2153.8,2253.0,2352.2,2451.4,2550.6,2649.8,2749.0,2848.2,2947.4,3046.6,3145.8,3245.0,3344.2,3443.4,3542.6,3641.8,3741.0,3840.2,3939.4,
  // --------- subrange 4 ---------
  4088.2,4286.6,4485.0,4683.4,4881.8,5080.2,5278.6,
  // =========== RANGE 2 ===========

  // --------- subrange 1 ---------
  4879.2,5077.6,5276.0,5474.4,5672.8,5871.2,6069.6,6268.0,6466.4,6664.8,6863.2,7061.6,7260.0,7458.4,7656.8,7855.2,
  // --------- subrange 2 ---------
  8152.8,8549.6,8946.4,9343.2,9740.0,10136.8,10533.6,10930.4,11327.2,11724.0,12120.8,12517.6,12914.4,13311.2,13708.0,14104.8,14501.6,14898.4,15295.2,15692.0,
  // --------- subrange 3 ---------
  16287.2,17080.8,17874.4,18668.0,19461.6,20255.2,21048.8,21842.4,22636.0,23429.6,24223.2,25016.8,25810.4,26604.0,27397.6,28191.2,28984.8,29778.4,30572.0,31365.6,32159.2,
  // --------- subrange 4 ---------
  33349.6,34936.8,36524.0,38111.2,39698.4,41285.6,42872.8,
  // =========== RANGE 3 ===========

  // --------- subrange 1 ---------
  39693.5,41280.5,42867.5,44454.5,46041.5,47628.5,49215.5,50802.5,52389.5,53976.5,55563.5,57150.5,58737.5,60324.5,61911.5,63498.5,
  // --------- subrange 2 ---------
  65879.0,69053.0,72227.0,75401.0,78575.0,81749.0,84923.0,88097.0,91271.0,94445.0,97619.0,100793.0,103967.0,107141.0,110315.0,113489.0,116663.0,119837.0,123011.0,126185.0,
  // --------- subrange 3 ---------
  130946.0,137294.0,143642.0,149990.0,156338.0,162686.0,169034.0,175382.0,181730.0,188078.0,194426.0,200774.0,207122.0,213470.0,219818.0,226166.0,232514.0,238862.0,245210.0,251558.0,257906.0,
  // --------- subrange 4 ---------
  267428.0,280124.0,292820.0,305516.0,318212.0,330908.0,343604.0

};

int HFMBoxMap[37]={0,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19};

int RunNo=0;int RunType=0;

int getmap()
{
//	ifstream semapinf("semapex.txt");
	ifstream semapinf("semapex_v2.txt");
	while(!semapinf.eof())
	{
		semapinf>>semin.crate>>semin.slot>>semin.fiber>>semin.channel>>semin.PMTname>>semin.Channelname>>semin.winchester>>semin.qiecard>>semin.qiecrate>>semin.qieslot>>semin.PMTid>>semin.BBid>>semin.ieta>>semin.iphi>>semin.depth>>semin.box>>semin.boxname>>semin.boxbarcode>>semin.VA>>semin.VB>>semin.VC>>semin.SX5gain>>semin.SX5gainerr>>semin.B904pedestal>>semin.B904pedestalsigma;
		
		sem.crate=semin.crate;
		sem.slot=semin.slot;
		sem.fiber=semin.fiber;
		sem.channel=semin.channel;
		sem.PMTname=semin.PMTname;
		sem.Channelname=semin.Channelname;
		sem.winchester=semin.winchester;
		sem.qiecard=semin.qiecard;
		sem.ieta=semin.ieta;
		sem.iphi=semin.iphi;
		sem.depth=semin.depth;
		sem.box=semin.box;
		sem.boxname=semin.boxname;
		sem.SX5gain=semin.SX5gain;
		sem.SX5gainerr=semin.SX5gainerr;
		sem.B904pedestal=semin.B904pedestal;
		SEM.push_back(sem);
	}
	semapinf.close();
}

int findQ(int b)
{
	int retb=-1;
	if(b>=1 && b<=9) return 1;
	else if(b>=10 && b<=18) return 2;
	else if(b>=19 && b<=27) return 3;
	else if(b>=28 && b<=36) return 4;
}

// int print()
// {
// 	char hname[500];
// 	sprintf(hname,"../NTuples/N_%d.root",RunNo);
// 	TFile* inroot=new TFile(hname);
// 	TTree *tree = (TTree*)inroot->Get("Events");
// 	if(RunType==1 || RunType==2)
// 	{
// 		tree->SetBranchAddress("ieta",&ed.ieta);
// 		tree->SetBranchAddress("iphi",&ed.iphi);
// 		tree->SetBranchAddress("depth",&ed.depth);
// 		tree->SetBranchAddress("pulse",&ed.pulse);
// 		tree->SetBranchAddress("tdc",&ed.tdc);
// 		tree->SetBranchAddress("capid",&ed.capid);
// 	}
// 	for(int i=0;i<tree->GetEntries();i++)
// 	{
// 		tree->GetEntry(i);
// 		for(int i1=0;i1<ed.ieta->size();i1++)
// 		{
// 			cout<<ed.ieta->at(i1)<<" "<<ed.iphi->at(i1)<<" "<<ed.depth->at(i1)<<endl;
// 			cout<<"                   "<<ed.pulse->at(i1)[0]<<" "<<ed.pulse->at(i1)[1]<<" "<<ed.pulse->at(i1)[2]<<" "<<ed.pulse->at(i1)[3]<<" "<<ed.pulse->at(i1)[4]<<" "<<ed.pulse->at(i1)[5]<<" "<<ed.pulse->at(i1)[6]<<" "<<ed.pulse->at(i1)[7]<<" "<<ed.pulse->at(i1)[8]<<" "<<ed.pulse->at(i1)[9]<<endl;
// 			cout<<"                   "<<ed.capid->at(i1)[0]<<" "<<ed.capid->at(i1)[1]<<" "<<ed.capid->at(i1)[2]<<" "<<ed.capid->at(i1)[3]<<" "<<ed.capid->at(i1)[4]<<" "<<ed.capid->at(i1)[5]<<" "<<ed.capid->at(i1)[6]<<" "<<ed.capid->at(i1)[7]<<" "<<ed.capid->at(i1)[8]<<" "<<ed.capid->at(i1)[9]<<endl;
// 			cout<<"                   "<<ed.tdc->at(i1)[0]<<" "<<ed.tdc->at(i1)[1]<<" "<<ed.tdc->at(i1)[2]<<" "<<ed.tdc->at(i1)[3]<<" "<<ed.tdc->at(i1)[4]<<" "<<ed.tdc->at(i1)[5]<<" "<<ed.tdc->at(i1)[6]<<" "<<ed.tdc->at(i1)[7]<<" "<<ed.tdc->at(i1)[8]<<" "<<ed.tdc->at(i1)[9]<<endl;
// 		}
// 	}
// 	inroot->Close();
// }

int plotpeds()
{
	char hname[500];
	sprintf(hname,"../NTuples/N_%d.root",RunNo);
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	tree->SetBranchAddress("ieta",&ed.ieta);
	tree->SetBranchAddress("iphi",&ed.iphi);
	tree->SetBranchAddress("depth",&ed.depth);
	tree->SetBranchAddress("pulse",&ed.pulse);
	tree->SetBranchAddress("tdc",&ed.tdc);
	tree->SetBranchAddress("capid",&ed.capid);
	
	pedlist pl;vector <pedlist> PL;int plind=-1;
	tree->GetEntry(0);
	for(int i1=0;i1<ed.ieta->size();i1++)
	{
		plind=-1;
		for(int i2=0;i2<PL.size();i2++)
		{
			if(SEM[PL[i2].mapind].ieta==ed.ieta->at(i1) && SEM[PL[i2].mapind].iphi==ed.iphi->at(i1) && SEM[PL[i2].mapind].depth==ed.depth->at(i1))
			{
				plind=i2;break;
			}
		}
		if(plind==-1)
		{
			pl.mapind=-1;
			for(int i2=0;i2<SEM.size();i2++)
			{
				if(SEM[i2].ieta==ed.ieta->at(i1) && SEM[i2].iphi==ed.iphi->at(i1) && SEM[i2].depth==ed.depth->at(i1))
				{
					pl.mapind=i2;break;
				}
			}
// 			if(pl.mapind==-1) cout<<ed.ieta->at(i1)<<" "<<ed.iphi->at(i1)<<" "<<ed.depth->at(i1)<<endl;
			pl.pedmean=0.;
			pl.pedsigma=0.;
			sprintf(hname,"%d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
			pl.h=new TH1F(hname,hname,10,-0.5,9.5);
			pl.h->SetLineWidth(2);pl.h->SetLineColor(1);
			PL.push_back(pl);
		}
	}
	for(int i=0;i<tree->GetEntries();i++)
// 	for(int i=0;i<10;i++)
	{
		tree->GetEntry(i);
		if(i%200==0) cout<<i<<" / "<<tree->GetEntries()<<endl;
		for(int i1=0;i1<ed.ieta->size();i1++)
		{
			plind=-1;
			for(int i2=0;i2<PL.size();i2++)
			{
				if(SEM[PL[i2].mapind].ieta==ed.ieta->at(i1) && SEM[PL[i2].mapind].iphi==ed.iphi->at(i1) && SEM[PL[i2].mapind].depth==ed.depth->at(i1))
				{
					plind=i2;break;
				}
			}
			for(int i2=0;i2<10;i2++)
			{
				PL[plind].h->Fill(ed.pulse->at(i1)[i2]);
			}
		}
	}
	double HMAX[2][4]={{0.}};
	int iside=0;int iquad=0;
	TF1* tf=new TF1("gaus","gaus",0.,10.);
	for(int i2=0;i2<PL.size();i2++)
	{
		iside=(SEM[PL[i2].mapind].box>0?0:1);
		iquad=(abs(SEM[PL[i2].mapind].box)-1)/9;
		if(PL[i2].h->GetMean()>0)
		{
			PL[i2].h->Fit(tf,"q","q",0.,10.);
			PL[i2].pedmean=tf->GetParameter(1);
			PL[i2].pedsigma=tf->GetParameter(2);
			tf->SetLineStyle(2);
			if(PL[i2].h->GetBinContent(PL[i2].h->GetMaximumBin())>HMAX[iside][iquad]){HMAX[iside][iquad]=PL[i2].h->GetBinContent(PL[i2].h->GetMaximumBin());}
		}
		else
		{
			PL[i2].pedmean=0.;
			PL[i2].pedsigma=0.;
		}
	}
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<4;i2++)
		{
			HMAX[i1][i2]+=20.;
		}
	}
// 	for(int i1=0;i1<SEM.size();i1++){cout<<SEM[i1].boxname<<" "<<SEM[i1].PMTname<<" "<<SEM[i1].Channelname<<endl;}
// 	for(int i1=0;i1<PL.size();i1++){cout<<SEM[PL[i1].mapind].boxname<<" "<<SEM[PL[i1].mapind].PMTname<<" "<<SEM[PL[i1].mapind].Channelname<<endl;}
	{
		TFile* outPeds=new TFile("Pedestals.root","recreate");
		char cname[400];
		char cname1[400];
		char bname[100];
		string pmtnames[24]={"A2","A4","A6","A8","A1","A3","A5","A7","B2","B4","B6","B8","B1","B3","B5","B7","C2","C4","C6","C8","C1","C3","C5","C7"};
		string chnames[2]={"1+2","3+4"};
		string side[2]={"P","M"};
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<4;i2++)
			{
				sprintf(cname,"HF%s_Q%d_Ped_details.pdf",side[i1].c_str(),(i2+1));
				for(int i3=0;i3<9;i3++)
				{
					TCanvas* cc1=new TCanvas("cc1","cc1",4500,6000);
					gStyle->SetOptStat(0);
					gStyle->SetTitleFontSize(0.1);
					cc1->Divide(8,6,0,0);
					int cc1i=1;
					if((i2*9+i3+1)<10){sprintf(bname,"HF%s_0%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
					else{sprintf(bname,"HF%s_%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
					string bname1(bname);
					for(int is1=0;is1<24;is1++)
					{
						for(int is2=0;is2<2;is2++)
						{
							bool found=false;plind=-1;
							for(int ik1=0;ik1<PL.size();ik1++)
							{
								if(SEM[PL[ik1].mapind].PMTname==pmtnames[is1] && SEM[PL[ik1].mapind].Channelname==chnames[is2] && SEM[PL[ik1].mapind].boxname==bname1)
								{
									plind=ik1;break;
								}
							}
							iside=(SEM[PL[plind].mapind].box>0?0:1);
							iquad=(abs(SEM[PL[plind].mapind].box)-1)/9;
							cc1->cd(cc1i);
							PL[plind].h->Draw();
							PL[plind].h->GetYaxis()->SetRangeUser(0.,HMAX[iside][iquad]);//PL[plind].h->SetLineWidth(2);PL[plind].h->SetLineColor(1);
							cc1i++;
							PL[plind].h->Write();
						}
					}
					if(i3==0) sprintf(cname1,"%s(",cname);
					else if(i3<8) sprintf(cname1,"%s",cname);
					else sprintf(cname1,"%s)",cname);
					cc1->SaveAs(cname1);
					delete cc1;
				}
				sprintf(hname,"mv %s ../Plots/%d",cname,RunNo);system(hname);
			}
		}
		outPeds->Close();
		for(int ik1=0;ik1<PL.size();ik1++)
		{
			PL[ik1].h=0;
		}
		sprintf(hname,"mv Pedestals.root ../Histos/%d",RunNo);system(hname);
	}
	
	ofstream outfile("log.txt");
	outfile<<"Run : "<<RunNo<<endl;
	{
		TFile* outGraphs=new TFile("PedestalSummaries.root","recreate");
		char cname[400];
		char bname[100];
		TGraphErrors* tg[2][4][9][2];
		string pmtnames[24]={"A2","A4","A6","A8","A1","A3","A5","A7","B2","B4","B6","B8","B1","B3","B5","B7","C2","C4","C6","C8","C1","C3","C5","C7"};
		string chnames[2]={"1+2","3+4"};
		string side[2]={"P","M"};
		int np=0;
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<4;i2++)
			{
				sprintf(cname,"HF%s_Q%d_Pedestals.pdf",side[i1].c_str(),(i2+1));
				TCanvas* cc1=new TCanvas("cc1","cc1",900,900);
				gStyle->SetOptStat(0);
				cc1->Divide(3,3,0,0);
				int cc1i=1;
				double ymin=1000.;double ymax=0.;
				for(int i3=0;i3<9;i3++)
				{
					for(int i4=0;i4<2;i4++)
					{
						sprintf(hname,"HF%s Q%d %d %s",side[i1].c_str(),i2+1,(i2*9+i3+1),chnames[i4].c_str());
						tg[i1][i2][i3][i4]=new TGraphErrors();
						tg[i1][i2][i3][i4]->SetName(hname);tg[i1][i2][i3][i4]->SetTitle(hname);
						tg[i1][i2][i3][i4]->SetMarkerColor(i4+1);tg[i1][i2][i3][i4]->SetMarkerStyle(i4+24);tg[i1][i2][i3][i4]->SetLineColor(i4+1);
						tg[i1][i2][i3][i4]->GetXaxis()->SetTitle("PMT ID");tg[i1][i2][i3][i4]->GetXaxis()->CenterTitle();
						tg[i1][i2][i3][i4]->GetYaxis()->SetTitle("Pedestal Mean (Sigma)");tg[i1][i2][i3][i4]->GetYaxis()->CenterTitle();
					
						if((i2*9+i3+1)<10){sprintf(bname,"HF%s_0%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
						else{sprintf(bname,"HF%s_%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
						string bname1(bname);
						
						for(int is1=0;is1<24;is1++)
						{
							bool found=false;plind=-1;
							for(int ik1=0;ik1<PL.size();ik1++)
							{
								if(SEM[PL[ik1].mapind].PMTname==pmtnames[is1] && SEM[PL[ik1].mapind].Channelname==chnames[i4] && SEM[PL[ik1].mapind].boxname==bname1)
								{
									plind=ik1;break;
								}
							}
							tg[i1][i2][i3][i4]->SetPoint(is1+1,((double)is1+1),PL[plind].pedmean);
							tg[i1][i2][i3][i4]->SetPointError(is1+1,0.,PL[plind].pedsigma);
							if(PL[plind].pedmean>ymax) ymax=PL[plind].pedmean;
							if(PL[plind].pedmean<ymin) ymin=PL[plind].pedmean;
							
							if(PL[plind].pedmean>5.5)
							{
								outfile<<"Pedestal too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
							if(PL[plind].pedmean<3.)
							{
								outfile<<"Pedestal too low: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
							if(PL[plind].pedsigma>1.)
							{
								outfile<<"Pedestal sigma too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
							if(PL[plind].pedsigma<0.6)
							{
								outfile<<"Pedestal sigma too low: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
						}
					}
				}
				ymax+=0.5;ymin-=0.5;
				if(ymax<10.) ymax=10.;
				if(ymin>0.) ymin=0.;
				for(int i3=0;i3<9;i3++)
				{
					cc1->cd(cc1i);
					sprintf(hname,"HF%s Q%d %d",side[i1].c_str(),i2+1,(i2*9+i3+1));
					tg[i1][i2][i3][0]->SetTitle(hname);
					tg[i1][i2][i3][0]->Draw("AP");
					tg[i1][i2][i3][0]->GetYaxis()->SetTitle("Pedestal ADC");tg[i1][i2][i3][0]->GetYaxis()->CenterTitle();
					tg[i1][i2][i3][0]->GetXaxis()->SetTitle("PMT ID");tg[i1][i2][i3][0]->GetXaxis()->CenterTitle();
					gPad->SetGridy(1);gPad->SetGridx(1);
					tg[i1][i2][i3][0]->GetYaxis()->SetRangeUser(ymin,ymax);tg[i1][i2][i3][0]->GetXaxis()->SetRangeUser(0.5,24.5);
					tg[i1][i2][i3][1]->Draw("same P");
					cc1i++;
					tg[i1][i2][i3][0]->Write();tg[i1][i2][i3][1]->Write();
				}
				cc1->SaveAs(cname);
				delete cc1;
				sprintf(hname,"mv %s ../Plots/%d",cname,RunNo);system(hname);
			}
		}
		outGraphs->Close();
		sprintf(hname,"mv PedestalSummaries.root ../Histos/%d",RunNo);system(hname);
	}
	{
		string pnames[12][4]={{"A2_3+4","A4_3+4","A6_3+4","A8_3+4"},{"A2_1+2","A4_1+2","A6_1+2","A8_1+2"},{"A1_3+4","A3_3+4","A5_3+4","A7_3+4"},{"A1_1+2","A3_1+2","A5_1+2","A7_1+2"},{"B2_3+4","B4_3+4","B6_3+4","B8_3+4"},{"B2_1+2","B4_1+2","B6_1+2","B8_1+2"},{"B1_3+4","B3_3+4","B5_3+4","B7_3+4"},{"B1_1+2","B3_1+2","B5_1+2","B7_1+2"},{"C2_3+4","C4_3+4","C6_3+4","C8_3+4"},{"C2_1+2","C4_1+2","C6_1+2","C8_1+2"},{"C1_3+4","C3_3+4","C5_3+4","C7_3+4"},{"C1_1+2","C3_1+2","C5_1+2","C7_1+2"}};
		TFile *inrootP=new TFile("HFMG.root");
		TMultiGraph *mg;
		TFile* outPoly=new TFile("Poly.root","recreate");
		TCanvas* cc1=new TCanvas("cc1","cc1",600,600);
		int bin=0;
		gStyle->SetOptStat(0);
		gStyle->SetTitleFontSize(0.06);
		gStyle->SetPalette(kRainBow);
		TH2Poly *HFPMean= new TH2Poly("HFP pedestal mean","HFP pedestal mean",-500,500,-500,500);
		TH2Poly *HFPSigma= new TH2Poly("HFP pedestal sigma","HFP pedestal sigma",-500,500,-500,500);
		TH2Poly *HFPMeanDiff= new TH2Poly("HFP (Ped_{P5}-Ped_{B904})/Ped_{B904}","HFP (Ped_{P5}-Ped_{B904})/Ped_{B904}",-500,500,-500,500);
		TH2Poly *HFMMean= new TH2Poly("HFM pedestal mean","HFM pedestal mean",-500,500,-500,500);
		TH2Poly *HFMSigma= new TH2Poly("HFM pedestal sigma","HFM pedestal sigma",-500,500,-500,500);
		TH2Poly *HFMMeanDiff= new TH2Poly("HFM (Ped_{P5}-Ped_{B904})/Ped_{B904}","HFM (Ped_{P5}-Ped_{B904})/Ped_{B904}",-500,500,-500,500);
		inrootP->GetObject("HF",mg);
		TLine *line1 = new TLine(-500.,0.,500.,0.);line1->SetLineColor(1);line1->SetLineWidth(2);
		TLine *line2 = new TLine(0.,-500.,0.,500.);line2->SetLineColor(1);line2->SetLineWidth(2);
		double pdiff=0.;
		for(int i1=0;i1<36;i1++)
		{
			for(int i2=0;i2<12;i2++)
			{
				for(int i3=0;i3<4;i3++)
				{
					sprintf(hname,"B%d_%s",(i1+1),pnames[i2][i3].c_str());
					bin=HFPMean->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFPSigma->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFPMeanDiff->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFMMean->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFMSigma->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFMMeanDiff->AddBin(mg->GetListOfGraphs()->FindObject(hname));
				}
			}
		}
		for(int i2=0;i2<PL.size();i2++)
		{
			pdiff=(PL[i2].pedmean-SEM[PL[i2].mapind].B904pedestal)/SEM[PL[i2].mapind].B904pedestal;
			if(SEM[PL[i2].mapind].box>0)
			{
				sprintf(hname,"B%d_%s_%s",SEM[PL[i2].mapind].box,SEM[PL[i2].mapind].PMTname.c_str(),SEM[PL[i2].mapind].Channelname.c_str());
				HFPMean->Fill(hname,PL[i2].pedmean);
				HFPSigma->Fill(hname,PL[i2].pedsigma);
				HFPMeanDiff->Fill(hname,pdiff);
			}
			else
			{
				sprintf(hname,"B%d_%s_%s",HFMBoxMap[abs(SEM[PL[i2].mapind].box)],SEM[PL[i2].mapind].PMTname.c_str(),SEM[PL[i2].mapind].Channelname.c_str());
				HFMMean->Fill(hname,PL[i2].pedmean);
				HFMSigma->Fill(hname,PL[i2].pedsigma);
				HFMMeanDiff->Fill(hname,pdiff);
			}
			if(fabs(pdiff)>0.1)
			{
				outfile<<"Pedestal difference too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
			}
		}
		HFPMean->SetMinimum(0.);HFPMean->SetMaximum(8.);HFPMean->Draw("colz");HFPMean->GetXaxis()->SetLabelSize(0);HFPMean->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFPMean.png");
		HFPSigma->SetMinimum(0.4);HFPSigma->SetMaximum(1.2);HFPSigma->Draw("colz");HFPSigma->GetXaxis()->SetLabelSize(0);HFPSigma->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFPSigma.png");
		HFPMeanDiff->SetMinimum(-0.2);HFPMeanDiff->SetMaximum(0.2);HFPMeanDiff->Draw("colz");HFPMeanDiff->GetXaxis()->SetLabelSize(0);HFPMeanDiff->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFPMeanDiff.png");
		HFMMean->SetMinimum(0.);HFMMean->SetMaximum(8.);HFMMean->Draw("colz");HFMMean->GetXaxis()->SetLabelSize(0);HFMMean->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFMMean.png");
		HFMSigma->SetMinimum(0.4);HFMSigma->SetMaximum(1.2);HFMSigma->Draw("colz");HFMSigma->GetXaxis()->SetLabelSize(0);HFMSigma->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFMSigma.png");
		HFMMeanDiff->SetMinimum(-0.2);HFMMeanDiff->SetMaximum(0.2);HFMMeanDiff->Draw("colz");HFMMeanDiff->GetXaxis()->SetLabelSize(0);HFMMeanDiff->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFMMeanDiff.png");
		
		HFPMean->Write();
		HFPSigma->Write();
		HFPMeanDiff->Write();
		HFMMean->Write();
		HFMSigma->Write();
		HFMMeanDiff->Write();
		outPoly->Close();
		inrootP->Close();
		
		sprintf(hname,"mv Poly.root ../Histos/%d",RunNo);system(hname);
		sprintf(hname,"mv *.png ../Plots/%d",RunNo);system(hname);
	}
	outfile.close();
	sprintf(hname,"mv log.txt ../Plots/%d",RunNo);system(hname);
	inroot->Close();
}

int plotleds()
{
	char hname[500];
	sprintf(hname,"../NTuples/N_%d.root",RunNo);
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	tree->SetBranchAddress("ieta",&ed.ieta);
	tree->SetBranchAddress("iphi",&ed.iphi);
	tree->SetBranchAddress("depth",&ed.depth);
	tree->SetBranchAddress("pulse",&ed.pulse);
	tree->SetBranchAddress("tdc",&ed.tdc);
	tree->SetBranchAddress("capid",&ed.capid);
	
	ledlist pl;vector <ledlist> PL;int plind=-1;
	tree->GetEntry(0);
	for(int i1=0;i1<ed.ieta->size();i1++)
	{
		plind=-1;
		for(int i2=0;i2<PL.size();i2++)
		{
			if(SEM[PL[i2].mapind].ieta==ed.ieta->at(i1) && SEM[PL[i2].mapind].iphi==ed.iphi->at(i1) && SEM[PL[i2].mapind].depth==ed.depth->at(i1))
			{
				plind=i2;break;
			}
		}
		if(plind==-1)
		{
			pl.mapind=-1;
			for(int i2=0;i2<SEM.size();i2++)
			{
				if(SEM[i2].ieta==ed.ieta->at(i1) && SEM[i2].iphi==ed.iphi->at(i1) && SEM[i2].depth==ed.depth->at(i1))
				{
					pl.mapind=i2;break;
				}
			}
			pl.Qmean=0.;
			pl.Qsigma=0.;
			pl.npemean=0.;
			pl.npeerr=0.;
			pl.gain=0.;
			pl.gainerr=0.;
			pl.nbadcapid=0.;
			sprintf(hname,"Pulse %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
			pl.p=new TH1F(hname,hname,10,-0.5,9.5);
			pl.p->GetXaxis()->SetTitle("TS (x25 ns)");pl.p->GetXaxis()->CenterTitle();
			pl.p->GetYaxis()->SetTitle("Mean Charge per TS (fC)");pl.p->GetYaxis()->CenterTitle();
			pl.p->SetFillColor(4);pl.p->SetLineColor(4);
			sprintf(hname,"PulseNorm %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
			pl.pn=new TH1F(hname,hname,10,-0.5,9.5);
			
			sprintf(hname,"TDC %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
			pl.t=new TH1F(hname,hname,10,-0.5,9.5);
			pl.t->GetXaxis()->SetTitle("TS (x25 ns)");pl.t->GetXaxis()->CenterTitle();
			pl.t->GetYaxis()->SetTitle("Mean TDC per TS");pl.t->GetYaxis()->CenterTitle();
			pl.t->SetLineColor(2);
			
			sprintf(hname,"Charge %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
			pl.q=new TH1F(hname,hname,10000,0.,250000.);
			pl.q->GetXaxis()->SetTitle("Charge (fC)");pl.q->GetXaxis()->CenterTitle();
			pl.q->GetYaxis()->SetTitle("Entries / 25 fC");pl.q->GetYaxis()->CenterTitle();
			pl.q->SetLineColor(1);pl.q->SetLineWidth(2);
			PL.push_back(pl);
		}
	}
	double ped=0.;double sig=0.;bool capidOK=true;int cap0=0;
	for(int i=0;i<tree->GetEntries();i++)
// 	for(int i=0;i<10;i++)
	{
		tree->GetEntry(i);
		if(i%200==0) cout<<i<<" / "<<tree->GetEntries()<<endl;
		for(int i1=0;i1<ed.ieta->size();i1++)
		{
			plind=-1;
			for(int i2=0;i2<PL.size();i2++)
			{
				if(SEM[PL[i2].mapind].ieta==ed.ieta->at(i1) && SEM[PL[i2].mapind].iphi==ed.iphi->at(i1) && SEM[PL[i2].mapind].depth==ed.depth->at(i1))
				{
					plind=i2;break;
				}
			}
			ped=0.;sig=0.;capidOK=true;
			for(int i2=0;i2<10;i2++)
			{
// 				cout<<i2<<" "<<ed.pulse->at(i1)[i2]<<" "<<adc2fC_QIE10[ed.pulse->at(i1)[i2]]<<endl;
				PL[plind].p->Fill(i2,adc2fC_QIE10[ed.pulse->at(i1)[i2]]);
				PL[plind].t->Fill(i2,ed.tdc->at(i1)[i2]);
				PL[plind].pn->Fill(i2);
				if(i2<=2) ped+=adc2fC_QIE10[ed.pulse->at(i1)[i2]];
				if(i2>=3 && i2<=7) sig+=adc2fC_QIE10[ed.pulse->at(i1)[i2]];
				if(i2>0)
				{
					if(!((ed.capid->at(i1)[i2]-ed.capid->at(i1)[i2-1])==1 || (ed.capid->at(i1)[i2]-ed.capid->at(i1)[i2-1])==-3))
					{
						capidOK=false;
					}
				}
			}
			ped/=3.;
			sig-=(5.*ped);
			PL[plind].q->Fill(sig);
			if(!capidOK) PL[plind].nbadcapid+=1.;
		}
	}
	double Qxmax[2][4]={{0.}};double Qymax[2][4]={{0.}};
	double Pmax[2][4]={{0.}};
	int iside=0;int iquad=0;float fnevt=((float)tree->GetEntries());
	TF1* tf=new TF1("gaus","gaus",0.,250000.);
	for(int i2=0;i2<PL.size();i2++)
	{
		iside=(SEM[PL[i2].mapind].box>0?0:1);
		iquad=(abs(SEM[PL[i2].mapind].box)-1)/9;
		PL[i2].p->Divide(PL[i2].pn);
		PL[i2].t->Divide(PL[i2].pn);
		PL[i2].nbadcapid/=fnevt;
		if(PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin())>Pmax[iside][iquad]){Pmax[iside][iquad]=PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin());}
		if(PL[i2].q->GetMean()>5)
		{
			PL[i2].q->Fit(tf,"q","q",0.,250000.);
			PL[i2].q->Fit(tf,"q","q",tf->GetParameter(1)-1.5*tf->GetParameter(2),tf->GetParameter(1)+1.5*tf->GetParameter(2));
			PL[i2].Qmean=tf->GetParameter(1);
			PL[i2].Qsigma=tf->GetParameter(2);
			PL[i2].npemean=1.15*pow(tf->GetParameter(1)/tf->GetParameter(2),2.);
			PL[i2].npeerr=2.*PL[i2].npemean*sqrt(pow(tf->GetParError(1)/tf->GetParameter(1),2.)+pow(tf->GetParError(2)/tf->GetParameter(2),2.));
			PL[i2].gain=((pow(tf->GetParameter(2),2.)/(tf->GetParameter(1)*1.15)));
			PL[i2].gainerr=sqrt(pow(tf->GetParError(1),2.)+pow(sqrt(2.)*tf->GetParError(2),2.));
			tf->SetLineStyle(2);
			if(PL[i2].q->GetBinContent(PL[i2].q->GetMaximumBin())>Qymax[iside][iquad]){Qymax[iside][iquad]=PL[i2].q->GetBinContent(PL[i2].q->GetMaximumBin());}
			if(PL[i2].q->GetBinCenter(PL[i2].q->FindLastBinAbove(0))>Qxmax[iside][iquad]){Qxmax[iside][iquad]=PL[i2].q->GetBinCenter(PL[i2].q->FindLastBinAbove(0));}
// 			if(PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin())>1000) cout<<SEM[PL[i2].mapind].boxname<<" "<<SEM[PL[i2].mapind].PMTname<<" "<<PL[i2].p->GetBinContent(PL[i2].p->GetMaximumBin())<<endl;
		}
	}
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<4;i2++)
		{
			Pmax[i1][i2]+=100.;Qxmax[i1][i2]+=50.;Qymax[i1][i2]+=20.;
		}
	}
	{
		TFile* outLeds=new TFile("LEDDetails.root","recreate");
		char cname[400];
		char cname1[400];
		char dname[400];
		char dname1[400];
		char ename[400];
		char ename1[400];
		char bname[100];
		string pmtnames[24]={"A2","A4","A6","A8","A1","A3","A5","A7","B2","B4","B6","B8","B1","B3","B5","B7","C2","C4","C6","C8","C1","C3","C5","C7"};
		string chnames[2]={"1+2","3+4"};
		string side[2]={"P","M"};
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<4;i2++)
			{
				sprintf(cname,"HF%s_Q%d_PulseShapes.pdf",side[i1].c_str(),(i2+1));
				sprintf(dname,"HF%s_Q%d_Charges.pdf",side[i1].c_str(),(i2+1));
				sprintf(ename,"HF%s_Q%d_TDCShapes.pdf",side[i1].c_str(),(i2+1));
				for(int i3=0;i3<9;i3++)
				{
					TCanvas* cc1=new TCanvas("cc1","cc1",4500,6000);
					TCanvas* cc2=new TCanvas("cc2","cc2",4500,6000);
					TCanvas* cc3=new TCanvas("cc3","cc3",4500,6000);
					gStyle->SetOptStat(0);
					gStyle->SetTitleFontSize(0.1);
					cc1->Divide(8,6,0,0);
					cc2->Divide(8,6,0,0);
					cc3->Divide(8,6,0,0);
					int cc1i=1;int cc2i=1;int cc3i=1;
					if((i2*9+i3+1)<10){sprintf(bname,"HF%s_0%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
					else{sprintf(bname,"HF%s_%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
					string bname1(bname);
					for(int is1=0;is1<24;is1++)
					{
						for(int is2=0;is2<2;is2++)
						{
							bool found=false;plind=-1;
							for(int ik1=0;ik1<PL.size();ik1++)
							{
								if(SEM[PL[ik1].mapind].PMTname==pmtnames[is1] && SEM[PL[ik1].mapind].Channelname==chnames[is2] && SEM[PL[ik1].mapind].boxname==bname1)
								{
									plind=ik1;break;
								}
							}
							
							iside=(SEM[PL[plind].mapind].box>0?0:1);
							iquad=(abs(SEM[PL[plind].mapind].box)-1)/9;
							
							cc1->cd(cc1i);
// 							PL[plind].p->Divide(PL[plind].pn);
							PL[plind].p->Draw("hist");
							PL[plind].p->GetYaxis()->SetRangeUser(-10.,Pmax[iside][iquad]);
							cc1i++;
							PL[plind].p->Write();
							
							cc2->cd(cc2i);
							PL[plind].q->Draw();
							PL[plind].q->GetYaxis()->SetRangeUser(0.,Qymax[iside][iquad]);PL[plind].q->GetXaxis()->SetRangeUser(0.,Qxmax[iside][iquad]);
							cc2i++;
							PL[plind].q->Write();
							
							cc3->cd(cc3i);
							PL[plind].t->Draw();
							PL[plind].t->GetYaxis()->SetRangeUser(0.,75.);
							cc3i++;
							PL[plind].t->Write();
						}
					}
					if(i3==0) sprintf(cname1,"%s(",cname);
					else if(i3<8) sprintf(cname1,"%s",cname);
					else sprintf(cname1,"%s)",cname);
					cc1->SaveAs(cname1);
					delete cc1;
					
					if(i3==0) sprintf(dname1,"%s(",dname);
					else if(i3<8) sprintf(dname1,"%s",dname);
					else sprintf(dname1,"%s)",dname);
					cc2->SaveAs(dname1);
					delete cc2;
					
					if(i3==0) sprintf(ename1,"%s(",ename);
					else if(i3<8) sprintf(ename1,"%s",ename);
					else sprintf(ename1,"%s)",ename);
					cc3->SaveAs(ename1);
					delete cc3;
				}
				sprintf(hname,"mv %s ../Plots/%d",cname,RunNo);system(hname);
				sprintf(hname,"mv %s ../Plots/%d",dname,RunNo);system(hname);
				sprintf(hname,"mv %s ../Plots/%d",ename,RunNo);system(hname);
			}
		}
		outLeds->Close();
		for(int ik1=0;ik1<PL.size();ik1++)
		{
			PL[ik1].p=0;PL[ik1].pn=0;PL[ik1].q=0;PL[ik1].t=0;
		}
		sprintf(hname,"mv LEDDetails.root ../Histos/%d",RunNo);system(hname);
	}
	ofstream outfile("log.txt");
	outfile<<"Run : "<<RunNo<<endl;
	{
		TFile* outGraphs=new TFile("LEDSummaries.root","recreate");
		char cname[400];
		char cname1[400];
		char bname[100];
		TGraphErrors* tg1[2][4][9][2];
		TGraphErrors* tg2[2][4][9][2];
		TGraphErrors* tg3[2][4][9][2];
		string pmtnames[24]={"A2","A4","A6","A8","A1","A3","A5","A7","B2","B4","B6","B8","B1","B3","B5","B7","C2","C4","C6","C8","C1","C3","C5","C7"};
		string chnames[2]={"1+2","3+4"};
		string side[2]={"P","M"};
		int np=0;
		double gdiff=0.;double gdifferr=0.;
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<4;i2++)
			{
				sprintf(cname,"HF%s_Q%d_LEDs.pdf",side[i1].c_str(),(i2+1));
				TCanvas* cc1=new TCanvas("cc1","cc1",900,900);
				TCanvas* cc2=new TCanvas("cc2","cc2",900,900);
				TCanvas* cc3=new TCanvas("cc3","cc3",900,900);
				gStyle->SetOptStat(0);
				cc1->Divide(3,3,0,0);
				cc2->Divide(3,3,0,0);
				cc3->Divide(3,3,0,0);
				int cc1i=1;int cc2i=1;int cc3i=1;
				double ymin1=1000.;double ymax1=0.;
				double ymin2=1000.;double ymax2=0.;
				double ymin3=1000.;double ymax3=0.;
				for(int i3=0;i3<9;i3++)
				{
					for(int i4=0;i4<2;i4++)
					{
						sprintf(hname,"HF%s Q%d %d %s NPE",side[i1].c_str(),i2+1,(i2*9+i3+1),chnames[i4].c_str());
						tg1[i1][i2][i3][i4]=new TGraphErrors();
						tg1[i1][i2][i3][i4]->SetName(hname);tg1[i1][i2][i3][i4]->SetTitle(hname);
						tg1[i1][i2][i3][i4]->SetMarkerColor(i4+1);tg1[i1][i2][i3][i4]->SetMarkerStyle(i4+24);tg1[i1][i2][i3][i4]->SetLineColor(i4+1);
						tg1[i1][i2][i3][i4]->GetXaxis()->SetTitle("PMT ID");tg1[i1][i2][i3][i4]->GetXaxis()->CenterTitle();
						tg1[i1][i2][i3][i4]->GetYaxis()->SetTitle("NPE");tg1[i1][i2][i3][i4]->GetYaxis()->CenterTitle();
						
						sprintf(hname,"HF%s Q%d %d %s Gains",side[i1].c_str(),i2+1,(i2*9+i3+1),chnames[i4].c_str());
						tg2[i1][i2][i3][i4]=new TGraphErrors();
						tg2[i1][i2][i3][i4]->SetName(hname);tg2[i1][i2][i3][i4]->SetTitle(hname);
						tg2[i1][i2][i3][i4]->SetMarkerColor(i4+1);tg2[i1][i2][i3][i4]->SetMarkerStyle(i4+24);tg2[i1][i2][i3][i4]->SetLineColor(i4+1);
						tg2[i1][i2][i3][i4]->GetXaxis()->SetTitle("PMT ID");tg2[i1][i2][i3][i4]->GetXaxis()->CenterTitle();
						tg2[i1][i2][i3][i4]->GetYaxis()->SetTitle("Gain (fC)");tg2[i1][i2][i3][i4]->GetYaxis()->CenterTitle();
						
						sprintf(hname,"HF%s Q%d %d %s GainDiff",side[i1].c_str(),i2+1,(i2*9+i3+1),chnames[i4].c_str());
						tg3[i1][i2][i3][i4]=new TGraphErrors();
						tg3[i1][i2][i3][i4]->SetName(hname);tg3[i1][i2][i3][i4]->SetTitle(hname);
						tg3[i1][i2][i3][i4]->SetMarkerColor(i4+1);tg3[i1][i2][i3][i4]->SetMarkerStyle(i4+24);tg3[i1][i2][i3][i4]->SetLineColor(i4+1);
						tg3[i1][i2][i3][i4]->GetXaxis()->SetTitle("PMT ID");tg3[i1][i2][i3][i4]->GetXaxis()->CenterTitle();
						tg3[i1][i2][i3][i4]->GetYaxis()->SetTitle("(G_{P5}-G_{SX5})/G_{SX5}");tg3[i1][i2][i3][i4]->GetYaxis()->CenterTitle();
					
						if((i2*9+i3+1)<10){sprintf(bname,"HF%s_0%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
						else{sprintf(bname,"HF%s_%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
						string bname1(bname);
						
						for(int is1=0;is1<24;is1++)
						{
							bool found=false;plind=-1;
							for(int ik1=0;ik1<PL.size();ik1++)
							{
								if(SEM[PL[ik1].mapind].PMTname==pmtnames[is1] && SEM[PL[ik1].mapind].Channelname==chnames[i4] && SEM[PL[ik1].mapind].boxname==bname1)
								{
									plind=ik1;break;
								}
							}
							tg1[i1][i2][i3][i4]->SetPoint(is1+1,((double)is1+1),PL[plind].npemean);
							tg1[i1][i2][i3][i4]->SetPointError(is1+1,0.,PL[plind].npeerr);
							if(PL[plind].npemean>ymax1) ymax1=PL[plind].npemean;
							if(PL[plind].npemean<ymin1) ymin1=PL[plind].npemean;
							
							tg2[i1][i2][i3][i4]->SetPoint(is1+1,((double)is1+1),PL[plind].gain);
							tg2[i1][i2][i3][i4]->SetPointError(is1+1,0.,PL[plind].gainerr);
							if(PL[plind].gain>ymax2) ymax2=PL[plind].gain;
							if(PL[plind].gain<ymin2) ymin2=PL[plind].gain;
							
							gdiff=(PL[plind].gain-SEM[PL[plind].mapind].SX5gain)/SEM[PL[plind].mapind].SX5gain;
							gdifferr=gdiff*sqrt(pow(PL[plind].gainerr/PL[plind].gain,2.)+pow(SEM[PL[plind].mapind].SX5gainerr/SEM[PL[plind].mapind].SX5gain,2.));
							tg3[i1][i2][i3][i4]->SetPoint(is1+1,((double)is1+1),gdiff);
							tg3[i1][i2][i3][i4]->SetPointError(is1+1,0.,gdifferr);
							if(gdiff>ymax3) ymax3=gdiff;
							if(gdiff<ymin3) ymin3=gdiff;
							
							if(PL[plind].npemean>35.)
							{
								outfile<<"NPE too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
							if(PL[plind].npemean<5.)
							{
								outfile<<"NPE too low: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
							if(PL[plind].gain>35.)
							{
								outfile<<"Gain too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
							if(PL[plind].gain<15.)
							{
								outfile<<"Gain too low: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
							if(fabs(gdiff)>0.2)
							{
								outfile<<"Gain difference too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
						}
					}
				}
				ymax1+=5.;ymin1-=5.;
				ymax2+=5.;ymin2-=5.;
				for(int i3=0;i3<9;i3++)
				{
					cc1->cd(cc1i);
					sprintf(hname,"HF%s Q%d %d NPE",side[i1].c_str(),i2+1,(i2*9+i3+1));
					tg1[i1][i2][i3][0]->SetTitle(hname);
					tg1[i1][i2][i3][0]->Draw("AP");
					tg1[i1][i2][i3][0]->GetYaxis()->SetTitle("NPE");tg1[i1][i2][i3][0]->GetYaxis()->CenterTitle();
					tg1[i1][i2][i3][0]->GetXaxis()->SetTitle("PMT ID");tg1[i1][i2][i3][0]->GetXaxis()->CenterTitle();
					gPad->SetGridy(1);gPad->SetGridx(1);
					tg1[i1][i2][i3][0]->GetYaxis()->SetRangeUser(ymin1,ymax1);tg1[i1][i2][i3][0]->GetXaxis()->SetRangeUser(0.5,24.5);
					tg1[i1][i2][i3][1]->Draw("same P");
					cc1i++;
					tg1[i1][i2][i3][0]->Write();tg1[i1][i2][i3][1]->Write();
					
					cc2->cd(cc2i);
					sprintf(hname,"HF%s Q%d %d Gain",side[i1].c_str(),i2+1,(i2*9+i3+1));
					tg2[i1][i2][i3][0]->SetTitle(hname);
					tg2[i1][i2][i3][0]->Draw("AP");
					tg2[i1][i2][i3][0]->GetYaxis()->SetTitle("Gain (fC)");tg2[i1][i2][i3][0]->GetYaxis()->CenterTitle();
					tg2[i1][i2][i3][0]->GetXaxis()->SetTitle("PMT ID");tg2[i1][i2][i3][0]->GetXaxis()->CenterTitle();
					gPad->SetGridy(1);gPad->SetGridx(1);
					tg2[i1][i2][i3][0]->GetYaxis()->SetRangeUser(ymin2,ymax2);tg2[i1][i2][i3][0]->GetXaxis()->SetRangeUser(0.5,24.5);
					tg2[i1][i2][i3][1]->Draw("same P");
					cc2i++;
					tg2[i1][i2][i3][0]->Write();tg2[i1][i2][i3][1]->Write();
					
					cc3->cd(cc3i);
					sprintf(hname,"HF%s Q%d %d GainDiff",side[i1].c_str(),i2+1,(i2*9+i3+1));
					tg3[i1][i2][i3][0]->SetTitle(hname);
					tg3[i1][i2][i3][0]->Draw("AP");
					tg3[i1][i2][i3][0]->GetYaxis()->SetTitle("(G_{P5}-G_{SX5})/G_{SX5}");tg3[i1][i2][i3][0]->GetYaxis()->CenterTitle();
					tg3[i1][i2][i3][0]->GetXaxis()->SetTitle("PMT ID");tg3[i1][i2][i3][0]->GetXaxis()->CenterTitle();
					gPad->SetGridy(1);gPad->SetGridx(1);
					tg3[i1][i2][i3][0]->GetYaxis()->SetRangeUser(-1.,1.);tg3[i1][i2][i3][0]->GetXaxis()->SetRangeUser(0.5,24.5);
					tg3[i1][i2][i3][1]->Draw("same P");
					cc3i++;
					tg3[i1][i2][i3][0]->Write();tg3[i1][i2][i3][1]->Write();
				}
				sprintf(cname1,"%s(",cname);
				cc1->SaveAs(cname1);
				sprintf(cname1,"%s",cname);
				cc2->SaveAs(cname1);
				sprintf(cname1,"%s)",cname);
				cc3->SaveAs(cname1);
				delete cc1;delete cc2;delete cc3;
				sprintf(hname,"mv %s ../Plots/%d",cname,RunNo);system(hname);
			}
		}
		outGraphs->Close();
		sprintf(hname,"mv LEDSummaries.root ../Histos/%d",RunNo);system(hname);
	}
	{
		TFile* outGraphs=new TFile("CapIDSummaries.root","recreate");
		char cname[400];
		char cname1[400];
		char bname[100];
		TGraphErrors* tg1[2][4][9][2];
		string pmtnames[24]={"A2","A4","A6","A8","A1","A3","A5","A7","B2","B4","B6","B8","B1","B3","B5","B7","C2","C4","C6","C8","C1","C3","C5","C7"};
		string chnames[2]={"1+2","3+4"};
		string side[2]={"P","M"};
		int np=0;
		double gdiff=0.;double gdifferr=0.;
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<4;i2++)
			{
				sprintf(cname,"HF%s_Q%d_CapIDs.pdf",side[i1].c_str(),(i2+1));
				TCanvas* cc1=new TCanvas("cc1","cc1",900,900);
				gStyle->SetOptStat(0);
				cc1->Divide(3,3,0,0);
				int cc1i=1;
				for(int i3=0;i3<9;i3++)
				{
					for(int i4=0;i4<2;i4++)
					{
						sprintf(hname,"HF%s Q%d %d %s Bad CapID Fraction",side[i1].c_str(),i2+1,(i2*9+i3+1),chnames[i4].c_str());
						tg1[i1][i2][i3][i4]=new TGraphErrors();
						tg1[i1][i2][i3][i4]->SetName(hname);tg1[i1][i2][i3][i4]->SetTitle(hname);
						tg1[i1][i2][i3][i4]->SetMarkerColor(i4+1);tg1[i1][i2][i3][i4]->SetMarkerStyle(i4+24);tg1[i1][i2][i3][i4]->SetLineColor(i4+1);
						tg1[i1][i2][i3][i4]->GetXaxis()->SetTitle("PMT ID");tg1[i1][i2][i3][i4]->GetXaxis()->CenterTitle();
						tg1[i1][i2][i3][i4]->GetYaxis()->SetTitle("Bad CapID Fraction");tg1[i1][i2][i3][i4]->GetYaxis()->CenterTitle();
					
						if((i2*9+i3+1)<10){sprintf(bname,"HF%s_0%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
						else{sprintf(bname,"HF%s_%d_Q%d",side[i1].c_str(),(i2*9+i3+1),(i2+1));}
						string bname1(bname);
						
						for(int is1=0;is1<24;is1++)
						{
							bool found=false;plind=-1;
							for(int ik1=0;ik1<PL.size();ik1++)
							{
								if(SEM[PL[ik1].mapind].PMTname==pmtnames[is1] && SEM[PL[ik1].mapind].Channelname==chnames[i4] && SEM[PL[ik1].mapind].boxname==bname1)
								{
									plind=ik1;break;
								}
							}
							tg1[i1][i2][i3][i4]->SetPoint(is1+1,((double)is1+1),PL[plind].nbadcapid);
							tg1[i1][i2][i3][i4]->SetPointError(is1+1,0.,0.);
							
							if(PL[plind].nbadcapid>0.1)
							{
								outfile<<"CapID rotation failure too high: "<<SEM[PL[plind].mapind].boxname<<" "<<SEM[PL[plind].mapind].PMTname<<" "<<SEM[PL[plind].mapind].Channelname<<" QIE10 card: "<<SEM[PL[plind].mapind].qiecard<<" Winchester: "<<SEM[PL[plind].mapind].winchester<<endl;
							}
						}
					}
				}
				for(int i3=0;i3<9;i3++)
				{
					cc1->cd(cc1i);
					sprintf(hname,"HF%s Q%d %d Bad CI Frac",side[i1].c_str(),i2+1,(i2*9+i3+1));
					tg1[i1][i2][i3][0]->SetTitle(hname);
					tg1[i1][i2][i3][0]->Draw("AP");
					tg1[i1][i2][i3][0]->GetYaxis()->SetTitle("Bad CapID Fraction");tg1[i1][i2][i3][0]->GetYaxis()->CenterTitle();
					tg1[i1][i2][i3][0]->GetXaxis()->SetTitle("PMT ID");tg1[i1][i2][i3][0]->GetXaxis()->CenterTitle();
					gPad->SetGridy(1);gPad->SetGridx(1);
					tg1[i1][i2][i3][0]->GetYaxis()->SetRangeUser(-0.1,1.1);tg1[i1][i2][i3][0]->GetXaxis()->SetRangeUser(0.5,24.5);
					tg1[i1][i2][i3][1]->Draw("same P");
					cc1i++;
					tg1[i1][i2][i3][0]->Write();tg1[i1][i2][i3][1]->Write();
				}
				sprintf(cname1,"%s",cname);
				cc1->SaveAs(cname1);
				delete cc1;
				sprintf(hname,"mv %s ../Plots/%d",cname,RunNo);system(hname);
			}
		}
		outGraphs->Close();
		sprintf(hname,"mv CapIDSummaries.root ../Histos/%d",RunNo);system(hname);
	}
	{
		string pnames[12][4]={{"A2_3+4","A4_3+4","A6_3+4","A8_3+4"},{"A2_1+2","A4_1+2","A6_1+2","A8_1+2"},{"A1_3+4","A3_3+4","A5_3+4","A7_3+4"},{"A1_1+2","A3_1+2","A5_1+2","A7_1+2"},{"B2_3+4","B4_3+4","B6_3+4","B8_3+4"},{"B2_1+2","B4_1+2","B6_1+2","B8_1+2"},{"B1_3+4","B3_3+4","B5_3+4","B7_3+4"},{"B1_1+2","B3_1+2","B5_1+2","B7_1+2"},{"C2_3+4","C4_3+4","C6_3+4","C8_3+4"},{"C2_1+2","C4_1+2","C6_1+2","C8_1+2"},{"C1_3+4","C3_3+4","C5_3+4","C7_3+4"},{"C1_1+2","C3_1+2","C5_1+2","C7_1+2"}};
		TFile *inrootP=new TFile("HFMG.root");
		TMultiGraph *mg;
		TFile* outPoly=new TFile("Poly.root","recreate");
		TCanvas* cc1=new TCanvas("cc1","cc1",600,600);
		int bin=0;
		gStyle->SetOptStat(0);
		gStyle->SetTitleFontSize(0.06);
		gStyle->SetPalette(kRainBow);
		TH2Poly *HFPNPE= new TH2Poly("HFP NPE","HFP NPE",-500,500,-500,500);
		TH2Poly *HFPGain= new TH2Poly("HFP Gains","HFP Gains",-500,500,-500,500);
		TH2Poly *HFPGainDiff= new TH2Poly("HFP (G_{P5}-G_{SX5})/G_{SX5}","HFP (G_{P5}-G_{SX5})/G_{SX5}",-500,500,-500,500);
		TH2Poly *HFMNPE= new TH2Poly("HFM NPE","HFM NPE",-500,500,-500,500);
		TH2Poly *HFMGain= new TH2Poly("HFM Gains","HFM Gains",-500,500,-500,500);
		TH2Poly *HFMGainDiff= new TH2Poly("HFM (G_{P5}-G_{SX5})/G_{SX5}","HFM (G_{P5}-G_{SX5})/G_{SX5}",-500,500,-500,500);
		inrootP->GetObject("HF",mg);
		TLine *line1 = new TLine(-500.,0.,500.,0.);line1->SetLineColor(1);line1->SetLineWidth(2);
		TLine *line2 = new TLine(0.,-500.,0.,500.);line2->SetLineColor(1);line2->SetLineWidth(2);
		for(int i1=0;i1<36;i1++)
		{
			for(int i2=0;i2<12;i2++)
			{
				for(int i3=0;i3<4;i3++)
				{
					sprintf(hname,"B%d_%s",(i1+1),pnames[i2][i3].c_str());
					bin=HFPNPE->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFPGain->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFPGainDiff->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFMNPE->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFMGain->AddBin(mg->GetListOfGraphs()->FindObject(hname));
					bin=HFMGainDiff->AddBin(mg->GetListOfGraphs()->FindObject(hname));
				}
			}
		}
		for(int i2=0;i2<PL.size();i2++)
		{
			if(SEM[PL[i2].mapind].box>0)
			{
				sprintf(hname,"B%d_%s_%s",SEM[PL[i2].mapind].box,SEM[PL[i2].mapind].PMTname.c_str(),SEM[PL[i2].mapind].Channelname.c_str());
				HFPNPE->Fill(hname,PL[i2].npemean);
				HFPGain->Fill(hname,PL[i2].gain);
				HFPGainDiff->Fill(hname,(PL[i2].gain-SEM[PL[i2].mapind].SX5gain)/SEM[PL[i2].mapind].SX5gain);
			}
			else
			{
				sprintf(hname,"B%d_%s_%s",HFMBoxMap[abs(SEM[PL[i2].mapind].box)],SEM[PL[i2].mapind].PMTname.c_str(),SEM[PL[i2].mapind].Channelname.c_str());
				HFMNPE->Fill(hname,PL[i2].npemean);
				HFMGain->Fill(hname,PL[i2].gain);
				HFMGainDiff->Fill(hname,(PL[i2].gain-SEM[PL[i2].mapind].SX5gain)/SEM[PL[i2].mapind].SX5gain);
			}
		}
		HFPNPE->SetMinimum(0.);HFPNPE->SetMaximum(30.);HFPNPE->Draw("colz");HFPNPE->GetXaxis()->SetLabelSize(0);HFPNPE->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFPNPE.png");
		HFPGain->SetMinimum(0.);HFPGain->SetMaximum(50.);HFPGain->Draw("colz");HFPGain->GetXaxis()->SetLabelSize(0);HFPGain->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFPGain.png");
		HFPGainDiff->SetMinimum(-1.);HFPGainDiff->SetMaximum(1.);HFPGainDiff->Draw("colz");HFPGainDiff->GetXaxis()->SetLabelSize(0);HFPGainDiff->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFPGainDiff.png");
		HFMNPE->SetMinimum(0.);HFMNPE->SetMaximum(30.);HFMNPE->Draw("colz");HFMNPE->GetXaxis()->SetLabelSize(0);HFMNPE->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFMNPE.png");
		HFMGain->SetMinimum(0.);HFMGain->SetMaximum(50.);HFMGain->Draw("colz");HFMGain->GetXaxis()->SetLabelSize(0);HFMGain->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFMGain.png");
		HFMGainDiff->SetMinimum(-1.);HFMGainDiff->SetMaximum(1.);HFMGainDiff->Draw("colz");HFMGainDiff->GetXaxis()->SetLabelSize(0);HFMGainDiff->GetYaxis()->SetLabelSize(0);
		line1->Draw("same");line2->Draw("same");cc1->SaveAs("HFMGainDiff.png");
		
		HFPNPE->Write();
		HFPGain->Write();
		HFPGainDiff->Write();
		HFMNPE->Write();
		HFMGain->Write();
		HFMGainDiff->Write();
		outPoly->Close();
		inrootP->Close();
		
		sprintf(hname,"mv Poly.root ../Histos/%d",RunNo);system(hname);
		sprintf(hname,"mv *.png ../Plots/%d",RunNo);system(hname);
	}
	outfile.close();
	sprintf(hname,"mv log.txt ../Plots/%d",RunNo);system(hname);
	inroot->Close();
}

int main(int argc, char *argv[])
{
	int opt=atoi(argv[1]);
	RunNo=atoi(argv[2]);
	RunType=opt;
// 	RunType=atoi(argv[3]);
	
	getmap();
	
// 	if(opt==0)
// 	{
// 		print();
// 	}
	if(opt==1)
	{
		plotpeds();
	}
	else if(opt==2)
	{
		plotleds();
	}
}


















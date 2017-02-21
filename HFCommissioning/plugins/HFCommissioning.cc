// -*- C++ -*-
//
// Package:    HCALCommissioning2017/HFCommissioning
// Class:      HFCommissioning
// 
/**\class HFCommissioning HFCommissioning.cc HCALCommissioning2017/HFCommissioning/plugins/HFCommissioning.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Burak Bilki
//         Created:  Sun, 12 Feb 2017 23:06:35 GMT
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
#include "DataFormats/HcalDetId/interface/HcalOtherDetId.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "DataFormats/HcalDigi/interface/QIE10DataFrame.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalCalibDetId.h"
#include "EventFilter/HcalRawToDigi/interface/AMC13Header.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUHTRData.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/// Per Event Header Structure
struct eventHeader
{
	uint32_t cdf0;
	uint32_t cdf1;
	uint32_t cdf2;
	uint32_t cdf3;
	uint32_t h0;
	uint32_t h1;
	uint32_t h2;
	uint32_t h3;
};

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

class HFCommissioning : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
	public:
		explicit HFCommissioning(const edm::ParameterSet&);
		~HFCommissioning();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		TFile *_file;
		TTree *_tree;
		int histoFED;
		int EID;
		int numChannels;
		string outFileName;
		int runType;
		string digiCollection;
		
		edm::EDGetTokenT<HcalDataFrameContainer<QIE10DataFrame> > tok_QIE10DigiCollection_;
		edm::EDGetTokenT<HFDigiCollection> hf_token;
		edm::EDGetTokenT<FEDRawDataCollection> raw_token;  
		edm::Handle<QIE10DigiCollection> qie10DigiCollection;
		edm::Handle<FEDRawDataCollection> raw_collection;  
};

HFCommissioning::HFCommissioning(const edm::ParameterSet& iConfig)
{
	tok_QIE10DigiCollection_ = consumes<HcalDataFrameContainer<QIE10DataFrame> >(edm::InputTag("hcalDigis"));
	hf_token = consumes<HFDigiCollection>(edm::InputTag("hcalDigis"));
	raw_token = consumes<FEDRawDataCollection>(edm::InputTag("source"));
	runType = iConfig.getParameter<int>("RunType");//1:pedestal 2:LED 3:histogram
	outFileName=iConfig.getUntrackedParameter<string>("OutFileName");
	histoFED = iConfig.getParameter<int>("histoFED");
	_file = new TFile(outFileName.c_str(), "recreate");
	_tree = new TTree("Events", "Events");
	_tree->Branch("ieta", &ed.ieta);
	_tree->Branch("iphi", &ed.iphi);
	_tree->Branch("depth", &ed.depth);
	if(runType==1 || runType==2)//pedestal or LED
	{
		_tree->Branch("pulse", &ed.pulse);
		_tree->Branch("tdc", &ed.tdc);
		_tree->Branch("capid", &ed.capid);
	}
	else if(runType==3)//histogram
	{
		_tree->Branch("histo", &ed.histo);
	}
	numChannels=0;
	EID=0;
}

HFCommissioning::~HFCommissioning()
{
	_file->cd();
	_file->Write();
	_file->Close();
}

void HFCommissioning::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace hcal;
	
	iEvent.getByToken(tok_QIE10DigiCollection_,qie10DigiCollection);
	const QIE10DigiCollection& qie10dc=*(qie10DigiCollection);
	iEvent.getByToken(raw_token,raw_collection);  
	
	edm::ESHandle<HcalElectronicsMap> item;
	edm::ESHandle<HcalDbService> pSetup;
	iSetup.get<HcalDbRecord>().get(pSetup);
	iSetup.get<HcalElectronicsMapRcd>().get(item);
	const FEDRawData& raw = raw_collection->FEDData(histoFED);
	
	vector <Int_t> vpulse;
	vector <Int_t> vtdc;
	vector <Int_t> vcapid;
	if(runType==1 || runType==2)//pedestal or LED
	{
		for (unsigned int j=0; j < qie10dc.size(); j++)
		{
			QIE10DataFrame qie10df = static_cast<QIE10DataFrame>(qie10dc[j]);
			DetId detid = qie10df.detid();
			HcalDetId hcaldetid = HcalDetId(detid);
			ed.ieta->push_back(hcaldetid.ieta());
			ed.iphi->push_back(hcaldetid.iphi());
			ed.depth->push_back(hcaldetid.depth());
			int nTS = qie10df.samples();
			vpulse.clear();vtdc.clear();vcapid.clear();
			for(int i=0; i<nTS; ++i)
			{
// 				eda.pulse[i].push_back(qie10df[i].adc());
// 				eda.tdc[i].push_back(qie10df[i].le_tdc());
// 				eda.capid[i].push_back(qie10df[i].capid());
				
				vpulse.push_back(qie10df[i].adc());
				vtdc.push_back(qie10df[i].le_tdc());
				vcapid.push_back(qie10df[i].capid());
			}
			ed.pulse->push_back(vpulse);
			ed.tdc->push_back(vtdc);
			ed.capid->push_back(vcapid);
		}
	}
	else if(runType==3)//histogram -- need to edit for P5
	{
// 		//the histos
// 		const struct eventHeader* eh =(const struct eventHeader*)(raw.data());
// 		const uint32_t* pData = (const uint32_t*) raw.data(); 
// 		uint32_t numHistos  = ((eh->h3)>>16)&0xFFFF;
// 		uint32_t numBins    = ((eh->h3)>>1)&0x0000FFFE; //includes overflow and header word
// 		uint32_t fiber   = 0;
// 		uint32_t channel = 0;
// 		pData+=8;
// 		int iz=0;
// 		int nH=-1;
// 		
// 		for (unsigned int iHist = 0; iHist<numHistos; iHist++)
// 		{
// 			if(iHist==96) iz++;
// 			fiber   = (*pData>>7)&0x1F;
// 			channel = (*pData>>2)&0x1F;
// 			
// 			bool found=false;
// 			for(int is1=0;is1<3;is1++)
// 			{
// 				for(int is2=0;is2<24;is2++)
// 				{
// 					for(int is3=0;is3<2;is3++)
// 					{
// 						if(MAP2Ch[is1][is2][is3][0]==((int)fiber) && MAP2Ch[is1][is2][is3][1]==((int)channel) && MAP2Ch[is1][is2][is3][2]==(iz+1)) {found=true;break;}
// 					}
// 					if(found){break;}
// 				}
// 				if(found){break;}
// 			}
// 			if(found)
// 			{
// 				nH++;
// 				eda.ieta.push_back((int)fiber);
// 				eda.iphi.push_back((int)channel);
// 				eda.depth.push_back(iz+1);
// 			}
// 			for(unsigned int iBin = 0; iBin<numBins+1; iBin++)
// 			{
// 				if(iBin<60 && found)
// 				{
// 					eda.histo[iBin].push_back(pData[iBin+1]);
// 				}
// 			}
// 			if(iHist<(numHistos-1))
// 			{
// 				pData+=(numBins+2);
// 			}
// 		}
	}
	_tree->Fill();
	if(runType==1 || runType==2)//pedestal or LED
	{
		ed.ieta->clear();
		ed.iphi->clear();
		ed.depth->clear();
		ed.pulse->clear();
		ed.tdc->clear();
		ed.capid->clear();
	}
	else if(runType==3)//histogram -- need to edit for P5
	{
		
	}
	EID++;
}

void HFCommissioning::beginJob(){}

void HFCommissioning::endJob(){}

void HFCommissioning::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFCommissioning);

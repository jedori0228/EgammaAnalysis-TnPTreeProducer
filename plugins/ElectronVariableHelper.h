#ifndef _ELECTRONVARIABLEHELPER_H
#define _ELECTRONVARIABLEHELPER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
//#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

typedef edm::View<reco::Candidate> CandView;

template <class T>
class ElectronVariableHelper : public edm::EDProducer {
 public:
  explicit ElectronVariableHelper(const edm::ParameterSet & iConfig);
  virtual ~ElectronVariableHelper() ;

  virtual float getEffArea(float scEta);
  virtual float PassTrigMVAHNTightv4(float thismva, float scEta);
  virtual float PassTrigMVAHNLoose(float thismva, float scEta);

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
  
private:
  edm::EDGetTokenT<std::vector<T> > probesToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<BXVector<l1t::EGamma> > l1EGTkn;
  edm::EDGetTokenT<CandView> pfCandToken_;
  edm::EDGetTokenT<double> rhoLabel_;
  edm::EDGetToken electronsMiniAODToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
};

template<class T>
ElectronVariableHelper<T>::ElectronVariableHelper(const edm::ParameterSet & iConfig) :
  probesToken_(consumes<std::vector<T> >(iConfig.getParameter<edm::InputTag>("probes"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  l1EGTkn(consumes<BXVector<l1t::EGamma> >(iConfig.getParameter<edm::InputTag>("l1EGColl"))),
  rhoLabel_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"))),
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap")))
 {

  produces<edm::ValueMap<float> >("chi2");
  produces<edm::ValueMap<float> >("dz");
  produces<edm::ValueMap<float> >("dxy");
  produces<edm::ValueMap<float> >("dxysig");
  produces<edm::ValueMap<float> >("mva");
  produces<edm::ValueMap<float> >("HNMVATight");
  produces<edm::ValueMap<float> >("HNMVALoose");
  produces<edm::ValueMap<float> >("reliso03");
  produces<edm::ValueMap<float> >("passConversionVeto");
  produces<edm::ValueMap<float> >("missinghits");
  produces<edm::ValueMap<float> >("l1e");
  produces<edm::ValueMap<float> >("l1et");
  produces<edm::ValueMap<float> >("l1eta");
  produces<edm::ValueMap<float> >("l1phi");
  produces<edm::ValueMap<float> >("pfPt");

  if( iConfig.existsAs<edm::InputTag>("pfCandColl") ) {
    pfCandToken_ = consumes<CandView>(iConfig.getParameter<edm::InputTag>("pfCandColl"));
  }

  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electronsMiniAOD"));

}

template<class T>
ElectronVariableHelper<T>::~ElectronVariableHelper()
{}

template<class T>
void ElectronVariableHelper<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // read input
  edm::Handle<std::vector<T> > probes;
  edm::Handle<reco::VertexCollection> vtxH;
  
  iEvent.getByToken(probesToken_, probes);
  iEvent.getByToken(vtxToken_, vtxH);
  const reco::VertexRef vtx(vtxH, 0);

  reco::Vertex pv;
  if (vtxH->size()) pv = vtxH->at(0);
  GlobalPoint pVertex(pv.position().x(),pv.position().y(),pv.position().z());

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  double rhoIso = std::max(*(rhoHandle.product()), 0.0);

  edm::Handle<BXVector<l1t::EGamma> > l1Cands;
  iEvent.getByToken(l1EGTkn, l1Cands);

  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronsMiniAODToken_,electrons);
  edm::Handle<edm::ValueMap<float> > mvaValues;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  
  edm::Handle<CandView> pfCands;
  if( !pfCandToken_.isUninitialized() ) iEvent.getByToken(pfCandToken_,pfCands);

  // prepare vector for output
  std::vector<float> chi2Vals;
  std::vector<float> dzVals;
  std::vector<float> dxyVals;
  std::vector<float> dxysigVals;
  std::vector<float> mvaVals;
  std::vector<float> HNMVATightVals;
  std::vector<float> HNMVALooseVals;
  std::vector<float> reliso03Vals;
  std::vector<float> passConversionVetoVals;
  std::vector<float> mhVals;
  std::vector<float> l1EVals;
  std::vector<float> l1EtVals;
  std::vector<float> l1EtaVals;
  std::vector<float> l1PhiVals;
  std::vector<float> pfPtVals;

  typename std::vector<T>::const_iterator probe, endprobes = probes->end();

  int j=0;
  for (probe = probes->begin(); probe != endprobes; ++probe) {
    
    const auto el = electrons->ptrAt(j);

    chi2Vals.push_back(probe->gsfTrack()->normalizedChi2());
    dzVals.push_back(probe->gsfTrack()->dz(vtx->position()));
    dxyVals.push_back(probe->gsfTrack()->dxy(vtx->position()));

    TrajectoryStateOnSurface eleTSOS;
    reco::TransientTrack eletranstrack = trackBuilder->build( probe->gsfTrack() );
    eleTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), pVertex, eletranstrack.field());
    float this_sip = 0.;
    if(eleTSOS.isValid()) {
      std::pair<bool, Measurement1D> eleIPpair = IPTools::signedTransverseImpactParameter(eletranstrack, eleTSOS.globalDirection(), pv);
      float eleSignificanceIP = eleIPpair.second.significance();
      this_sip = eleSignificanceIP;
    }
    dxysigVals.push_back(this_sip);

    passConversionVetoVals.push_back( probe->passConversionVeto() );

    float scEta = probe->superCluster()->eta();
    float ecalpt = probe->ecalDrivenMomentum().pt();
    float elEffArea = getEffArea(scEta);
    float chIso = probe->pfIsolationVariables().sumChargedHadronPt;
    float nhIso = probe->pfIsolationVariables().sumNeutralHadronEt;
    float phIso = probe->pfIsolationVariables().sumPhotonEt;

    float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*elEffArea) )/ ecalpt;
    reliso03Vals.push_back( relIso );

    float mvaval = (*mvaValues)[el];
    mvaVals.push_back( mvaval );

    HNMVATightVals.push_back( PassTrigMVAHNTightv4(mvaval, scEta) );
    HNMVALooseVals.push_back( PassTrigMVAHNLoose(mvaval, scEta) );

    mhVals.push_back(float(probe->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)));

    float l1e = 999999.;    
    float l1et = 999999.;
    float l1eta = 999999.;
    float l1phi = 999999.;
    float pfpt = 999999.;
    float dRmin = 0.3;
    for (std::vector<l1t::EGamma>::const_iterator l1Cand = l1Cands->begin(0); l1Cand != l1Cands->end(0); ++l1Cand) {

      float dR = deltaR(l1Cand->eta(), l1Cand->phi() , probe->superCluster()->eta(), probe->superCluster()->phi());
      if (dR < dRmin) {
	dRmin = dR;
	l1e = l1Cand->energy();
	l1et = l1Cand->et();
        l1eta = l1Cand->eta();
        l1phi = l1Cand->phi();
      }
    }
    if( pfCands.isValid() )
    for( size_t ipf = 0; ipf < pfCands->size(); ++ipf ) {
        auto pfcand = pfCands->ptrAt(ipf);
	if( abs( pfcand->pdgId() ) != 11 ) continue;
	float dR = deltaR(pfcand->eta(), pfcand->phi() , probe->eta(), probe->phi());
	if( dR < 0.0001 ) pfpt = pfcand->pt();
    }

    l1EVals.push_back(l1e);
    l1EtVals.push_back(l1et);
    l1EtaVals.push_back(l1eta);
    l1PhiVals.push_back(l1phi);
    pfPtVals.push_back(pfpt);

    ++j;  
  }

  
  // convert into ValueMap and store
  std::unique_ptr<edm::ValueMap<float> > chi2ValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler chi2Filler(*chi2ValMap);
  chi2Filler.insert(probes, chi2Vals.begin(), chi2Vals.end());
  chi2Filler.fill();
  iEvent.put(std::move(chi2ValMap), "chi2");

  std::unique_ptr<edm::ValueMap<float> > dzValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler dzFiller(*dzValMap);
  dzFiller.insert(probes, dzVals.begin(), dzVals.end());
  dzFiller.fill();
  iEvent.put(std::move(dzValMap), "dz");

  std::unique_ptr<edm::ValueMap<float> > dxyValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler dxyFiller(*dxyValMap);
  dxyFiller.insert(probes, dxyVals.begin(), dxyVals.end());
  dxyFiller.fill();
  iEvent.put(std::move(dxyValMap), "dxy");

  std::unique_ptr<edm::ValueMap<float> > dxysigValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler dxysigFiller(*dxysigValMap);
  dxysigFiller.insert(probes, dxysigVals.begin(), dxysigVals.end());
  dxysigFiller.fill();
  iEvent.put(std::move(dxysigValMap), "dxysig");

  std::unique_ptr<edm::ValueMap<float> > mvaValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler mvaFiller(*mvaValMap);
  mvaFiller.insert(probes, mvaVals.begin(), mvaVals.end());
  mvaFiller.fill();
  iEvent.put(std::move(mvaValMap), "mva");

  std::unique_ptr<edm::ValueMap<float> > HNMVATightValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler HNMVATightFiller(*HNMVATightValMap);
  HNMVATightFiller.insert(probes, HNMVATightVals.begin(), HNMVATightVals.end());
  HNMVATightFiller.fill();
  iEvent.put(std::move(HNMVATightValMap), "HNMVATight");

  std::unique_ptr<edm::ValueMap<float> > HNMVALooseValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler HNMVALooseFiller(*HNMVALooseValMap);
  HNMVALooseFiller.insert(probes, HNMVALooseVals.begin(), HNMVALooseVals.end());
  HNMVALooseFiller.fill();
  iEvent.put(std::move(HNMVALooseValMap), "HNMVALoose");

  std::unique_ptr<edm::ValueMap<float> > reliso03ValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler reliso03Filler(*reliso03ValMap);
  reliso03Filler.insert(probes, reliso03Vals.begin(), reliso03Vals.end());
  reliso03Filler.fill();
  iEvent.put(std::move(reliso03ValMap), "reliso03");

  std::unique_ptr<edm::ValueMap<float> > passConversionVetoValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler passConversionVetoFiller(*passConversionVetoValMap);
  passConversionVetoFiller.insert(probes, passConversionVetoVals.begin(), passConversionVetoVals.end());
  passConversionVetoFiller.fill();
  iEvent.put(std::move(passConversionVetoValMap), "passConversionVeto");

  std::unique_ptr<edm::ValueMap<float> > mhValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler mhFiller(*mhValMap);
  mhFiller.insert(probes, mhVals.begin(), mhVals.end());
  mhFiller.fill();
  iEvent.put(std::move(mhValMap), "missinghits");

  std::unique_ptr<edm::ValueMap<float> > l1EValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler l1EFill(*l1EValMap);
  l1EFill.insert(probes, l1EVals.begin(), l1EVals.end());
  l1EFill.fill();
  iEvent.put(std::move(l1EValMap), "l1e");

  std::unique_ptr<edm::ValueMap<float> > l1EtValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler l1EtFill(*l1EtValMap);
  l1EtFill.insert(probes, l1EtVals.begin(), l1EtVals.end());
  l1EtFill.fill();
  iEvent.put(std::move(l1EtValMap), "l1et");

  std::unique_ptr<edm::ValueMap<float> > l1EtaValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler l1EtaFill(*l1EtaValMap);
  l1EtaFill.insert(probes, l1EtaVals.begin(), l1EtaVals.end());
  l1EtaFill.fill();
  iEvent.put(std::move(l1EtaValMap), "l1eta");

  std::unique_ptr<edm::ValueMap<float> > l1PhiValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler l1PhiFill(*l1PhiValMap);
  l1PhiFill.insert(probes, l1PhiVals.begin(), l1PhiVals.end());
  l1PhiFill.fill();
  iEvent.put(std::move(l1PhiValMap), "l1phi");

  std::unique_ptr<edm::ValueMap<float> > pfPtValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler pfPtFill(*pfPtValMap);
  pfPtFill.insert(probes, pfPtVals.begin(), pfPtVals.end());
  pfPtFill.fill();
  iEvent.put(std::move(pfPtValMap), "pfPt");

  
}

template<class T>
float ElectronVariableHelper<T>::getEffArea(float scEta){
  // ElectronEffectiveArea::ElectronEffectiveAreaTarget electronEATarget;
  // if ( runOnMC_ ) electronEATarget = ElectronEffectiveArea::kEleEAFall11MC;
  // else electronEATarget = ElectronEffectiveArea::kEleEAData2012;
  // if( dR < 0.35)
  //   return ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, scEta, electronEATarget);
  // else
  //   return ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, scEta, electronEATarget);

  // new effArea  https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt
  float absEta = std::abs(scEta);
  if ( 0.0000 >= absEta && absEta < 1.0000 ) return 0.1703;
  if ( 1.0000 >= absEta && absEta < 1.4790 ) return 0.1715;
  if ( 1.4790 >= absEta && absEta < 2.0000 ) return 0.1213;
  if ( 2.0000 >= absEta && absEta < 2.2000 ) return 0.1230;
  if ( 2.2000 >= absEta && absEta < 2.3000 ) return 0.1635;
  if ( 2.3000 >= absEta && absEta < 2.4000 ) return 0.1937;
  if ( 2.4000 >= absEta && absEta < 5.0000 ) return 0.2393;
  return 0;
}

template<class T>
float ElectronVariableHelper<T>::PassTrigMVAHNTightv4(float thismva, float scEta){

  float mva_cut=0.93;
  if(fabs(scEta) > 1.479) mva_cut=0.93;
  else if(fabs(scEta) > 0.8) mva_cut=0.825;
  else mva_cut=0.9;

  if(thismva > mva_cut) return 1.;
  return 0.;

}

template<class T>
float ElectronVariableHelper<T>::PassTrigMVAHNLoose(float thismva, float scEta){

  if( (fabs(scEta) < 0.8)                            && thismva > -0.02) return 1.;
  if( (fabs(scEta) > 0.8) && (fabs(scEta)  < 1.479)  && thismva > -0.52) return 1.;
  if( (fabs(scEta) < 2.5) && (fabs(scEta)  > 1.479)  && thismva > -0.52) return 1.;

  return 0.;

}










#endif

./pileupCalc.py \
-i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt \
--inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt \
--calcMode true --minBiasXsec 69200 \
--maxPileupBin 75 --numPileupBins 75 \
pileup_2016_35fb.root

./pileupCalc.py \
-i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt \
--inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt \
--calcMode true --minBiasXsec 69200 \
--maxPileupBin 99 --numPileupBins 99 \
pileup_2017_41fb.root


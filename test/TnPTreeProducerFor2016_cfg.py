import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

process = cms.Process("tnpEGM")

###################################################################
## argument line options
###################################################################
varOptions = VarParsing('analysis')

varOptions.register(
    "isMC", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute MC efficiencies"
    )

varOptions.register(
    "doEleID", True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Include tree for photon ID SF"
    )

varOptions.register(
    "doPhoID", True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Include tree for photon ID SF"
    )

varOptions.register(
    "doTrigger", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Include tree for Trigger SF"
    )

varOptions.register(
    "doRECO", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Include tree for Reco SF"
    )

varOptions.register(
    "calibEn", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "use EGM smearer to calibrate photon and electron energy"
    )

#### HLTname is HLT2 in reHLT samples
varOptions.register(
    "HLTname", "HLT",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "HLT process name (default HLT)"
    )

varOptions.register(
    #"GT","auto",
    #"GT","101X_dataRun2_Prompt_v9",
    #"GT","94X_dataRun2_ReReco_EOY17_v6",
    #"GT","80X_dataRun2_2016LegacyRepro_v4",
    "GT","80X_mcRun2_asymptotic_2016_TrancheIV_v8",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global Tag to be used"
    )

varOptions.register(
    "includeSUSY", False,
#    "isAOD", True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "add also the variables used by SUSY"
    )

varOptions.register(
    "isAOD", False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "use AOD"
    )

varOptions.parseArguments()


###################################################################
## Define TnP inputs 
###################################################################

options = dict()
options['useAOD']               = cms.bool(varOptions.isAOD)

options['HLTProcessName']       = varOptions.HLTname

### set input collections
options['ELECTRON_COLL']        = "slimmedElectrons"
options['PHOTON_COLL']          = "slimmedPhotons"
options['SUPERCLUSTER_COLL']    = "reducedEgamma:reducedSuperClusters" ### not used in AOD
if options['useAOD']:
    options['ELECTRON_COLL']    = "gedGsfElectrons"
    options['PHOTON_COLL'  ]    = "gedPhotons"


options['ELECTRON_CUTS']        = "ecalEnergy*sin(superClusterPosition.theta)>5.0 &&  (abs(-log(tan(superClusterPosition.theta/2)))<2.5)"
options['SUPERCLUSTER_CUTS']    = "abs(eta)<2.5 &&  et>5.0"
options['PHOTON_CUTS']          = "(abs(-log(tan(superCluster.position.theta/2)))<=2.5) && pt> 10"
options['ELECTRON_TAG_CUTS']    = "(abs(-log(tan(superCluster.position.theta/2)))<=2.1) && !(1.4442<=abs(-log(tan(superClusterPosition.theta/2)))<=1.566) && pt >= 30.0"

options['MAXEVENTS']            = cms.untracked.int32(varOptions.maxEvents) 
#options['MAXEVENTS']            = 2000
options['DoTrigger']            = cms.bool( varOptions.doTrigger )
options['DoRECO']               = cms.bool( varOptions.doRECO    )
options['DoEleID']              = cms.bool( varOptions.doEleID   )
options['DoPhoID']              = cms.bool( varOptions.doPhoID   )

options['OUTPUTEDMFILENAME']    = 'edmFile.root'
options['DEBUG']                = cms.bool(False)
options['isMC']                 = cms.bool(False)
options['UseCalibEn']           = varOptions.calibEn

options['addSUSY']               = varOptions.includeSUSY
if options['useAOD']: 
    options['addSUSY']               = cms.bool(False)

if (varOptions.isMC):
    options['isMC']                = cms.bool(True)
    options['OUTPUT_FILE_NAME']    = "TnPTree_mc.root"
    if varOptions.isAOD :  options['OUTPUT_FILE_NAME']    = "TnPTree_mc_aod.root"
    options['TnPPATHS']            = cms.vstring("HLT_Ele27_eta2p1_WPTight_Gsf_v*") #FOR 2016
    options['TnPHLTTagFilters']    = cms.vstring("hltEle27erWPTightGsfTrackIsoFilter") #FOR 2016
    options['HLTFILTERTOMEASURE']  = cms.vstring("hltEle27WPTightGsfTrackIsoFilter","hltEG175HEFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['GLOBALTAG']           = '94X_mcRun2_asymptotic_v3'
else:
    options['OUTPUT_FILE_NAME']    = "TnPTree_data.root"
    options['TnPPATHS']            = cms.vstring("HLT_Ele27_eta2p1_WPTight_Gsf_v*") #FOR 2016
    options['TnPHLTTagFilters']    = cms.vstring("hltEle27erWPTightGsfTrackIsoFilter") #FOR 2016
    options['HLTFILTERTOMEASURE']  = cms.vstring("hltEle27WPTightGsfTrackIsoFilter","hltEG175HEFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['GLOBALTAG']           = '94X_dataRun2_v10'

if varOptions.GT != "auto" :
    options['GLOBALTAG'] = varOptions.GT


###################################################################
## Define input files for test local run
###################################################################
#from EgammaAnalysis.TnPTreeProducer.etc.tnpInputTestFiles_cff import filesMiniAOD_Preliminary2018 as inputs
#if options['useAOD'] : from EgammaAnalysis.TnPTreeProducer.etc.tnpInputTestFiles_cff import filesAOD_23Sep2016 as inputs #switch to 2017 samples if want to cmsRun on AOD
#if options['useAOD'] : 
#from EgammaAnalysis.TnPTreeProducer.etc.tnpInputTestFiles_cff import filesAOD_Preliminary2017 as inputs #switch to 2017 samples if want to cmsRun on AOD
#from EgammaAnalysis.TnPTreeProducer.etc.tnpInputTestFiles_cff import filesMiniAOD_2016 as inputs #switch to 2017 samples if want to cmsRun on AOD
#from EgammaAnalysis.TnPTreeProducer.etc.tnpInputTestFiles_cff import filesMiniAOD_2016 as inputs #switch to 2017 samples if want to cmsRun on AOD
#if options['useAOD'] : from EgammaAnalysis.TnPTreeProducer.etc.tnpInputTestFiles_cff import filesAOD_empty as inputs #switch to 2017 samples if want to cmsRun on AOD

options['INPUT_FILE_NAME'] = cms.untracked.vstring(
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/A653C523-898B-E811-AE8D-0025905B85DC.root',
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/FE6AFFE3-798B-E811-86D1-0025905A613C.root',
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/326B5EDD-798B-E811-948C-0CC47A7C35D2.root',
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/944E4220-598B-E811-8FDC-0025905B85B8.root',
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/64E54575-5D8B-E811-8507-0025905B859A.root',
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/42D5EB6C-018B-E811-84E5-0CC47A78A478.root',
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/B81D0682-5F8A-E811-B3C1-0025905A6066.root',
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/7438C576-018B-E811-B4B7-0025905A60F2.root',
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/DAF9304C-F88A-E811-BAEB-0025905A60F2.root',
'/store/data/Run2016H/SingleElectron/MINIAOD/17Jul2018-v1/00000/E4193B73-018B-E811-B737-0CC47A7C35F8.root',
)
if varOptions.isMC:
  options['INPUT_FILE_NAME'] = cms.untracked.vstring(
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/ACDA5D95-3EDF-E811-AC6F-842B2B6AEE8B.root',
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/EE95C7CD-3ADF-E811-AEEF-842B2B6F8647.root',
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/C68F0587-3ADF-E811-9850-842B2B6AE7AB.root',
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/B0A6C7CD-3ADF-E811-9FD0-842B2B6F8647.root',
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/74DB5D95-3EDF-E811-AC93-842B2B6AEE8B.root',
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/C40DE067-44DF-E811-8A35-D4AE5280690B.root',
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/34DF9F61-44DF-E811-BBB5-D4AE5280690B.root',
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/F47E26F1-46DF-E811-B835-D4AE529017EB.root',
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/567020F4-46DF-E811-BBA7-D4AE529017EB.root',
'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/6ACED6F6-4DDF-E811-A191-782BCB3BCE9F.root',
)
#for file in open("file.list").readlines():
#    inputs['data'].append(file.strip())

###################################################################
## import TnP tree maker pythons and configure for AODs
###################################################################
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff') 
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options['GLOBALTAG'] , '')

import EgammaAnalysis.TnPTreeProducer.egmTreesSetup_cff as tnpSetup
tnpSetup.setupTreeMaker(process,options)



###################################################################
## Init and Load
###################################################################
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = options['INPUT_FILE_NAME'],
                            )
process.maxEvents = cms.untracked.PSet( input = options['MAXEVENTS'])
#process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(options['MAXEVENTS'] )
#)

if options['addSUSY']    : print "  -- Including variables for SUSY       -- "
if options['DoTrigger'] : print "  -- Producing HLT (trigger ele) efficiency tree -- "
if options['DoRECO']    : print "  -- Producing RECO SF tree        -- "
if options['DoEleID']   : print "  -- Producing electron SF tree    -- "
if options['DoPhoID']   : print "  -- Producing photon SF tree      -- "

    
###################################################################
## Define sequences and TnP pairs
###################################################################
process.cand_sequence = cms.Sequence( process.init_sequence + process.tag_sequence )
if options['addSUSY']                         : process.cand_sequence += process.susy_ele_sequence
if options['DoEleID'] or options['DoTrigger'] : process.cand_sequence += process.ele_sequence
if options['DoPhoID']                         : process.cand_sequence += process.pho_sequence
if options['DoTrigger']                       : process.cand_sequence += process.hlt_sequence
if options['DoRECO']                          : process.cand_sequence += process.sc_sequence

process.tnpPairs_sequence = cms.Sequence()
if options['DoTrigger'] : process.tnpPairs_sequence *= process.tnpPairingEleHLT
if options['DoRECO']    : process.tnpPairs_sequence *= process.tnpPairingEleRec
if options['DoEleID']   : process.tnpPairs_sequence *= process.tnpPairingEleIDs
if options['DoPhoID']   : process.tnpPairs_sequence *= process.tnpPairingPhoIDs

##########################################################################
## TnP Trees
##########################################################################
import EgammaAnalysis.TnPTreeProducer.egmTreesContent_cff as tnpVars
if options['useAOD']: tnpVars.setupTnPVariablesForAOD()
tnpVars.mcTruthCommonStuff.isMC = cms.bool(varOptions.isMC)

process.tnpEleTrig = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                    tnpVars.CommonStuffForGsfElectronProbe, tnpVars.mcTruthCommonStuff,
                                    tagProbePairs = cms.InputTag("tnpPairingEleHLT"),
                                    probeMatches  = cms.InputTag("genProbeEle"),
                                    allProbes     = cms.InputTag("probeEle"),
                                    flags = cms.PSet(
                                        passingHLT        = cms.InputTag("probeElePassHLT"),
                                        passingLoose94X   = cms.InputTag("probeEleCutBasedLoose94X" ),
                                        passingMedium94X  = cms.InputTag("probeEleCutBasedMedium94X"),
                                        passingTight94X   = cms.InputTag("probeEleCutBasedTight94X" ),
                                        passingHEEP       = cms.InputTag("probeEleHEEP"),
                                        ),
                                    )

process.tnpEleReco = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                    tnpVars.mcTruthCommonStuff, tnpVars.CommonStuffForSuperClusterProbe, 
                                    tagProbePairs = cms.InputTag("tnpPairingEleRec"),
                                    probeMatches  = cms.InputTag("genProbeSC"),
                                    allProbes     = cms.InputTag("probeSC"),
                                    flags         = cms.PSet(
        passingRECO   = cms.InputTag("probeSCEle", "superclusters"),
        passingRECOEcalDriven   = cms.InputTag("probeSCEle", "superclustersEcalDriven"),
        passingRECOTrackDriven   = cms.InputTag("probeSCEle", "superclustersTrackDriven")
        ),

                                    )

process.tnpEleIDs = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                    tnpVars.mcTruthCommonStuff, tnpVars.CommonStuffForGsfElectronProbe,
                                    tagProbePairs = cms.InputTag("tnpPairingEleIDs"),
                                    probeMatches  = cms.InputTag("genProbeEle"),
                                    allProbes     = cms.InputTag("probeEle"),
                                    flags         = cms.PSet(

                                        passingVeto80X    = cms.InputTag("probeEleCutBasedVeto80X"  ),
                                        passingLoose80X   = cms.InputTag("probeEleCutBasedLoose80X" ),
                                        passingMedium80X  = cms.InputTag("probeEleCutBasedMedium80X"),
                                        passingTight80X   = cms.InputTag("probeEleCutBasedTight80X" ),
                                        passingMVA80Xwp90 = cms.InputTag("probeEleMVA80Xwp90" ),
                                        passingMVA80Xwp80 = cms.InputTag("probeEleMVA80Xwp80" ),

                                        passingVeto94X    = cms.InputTag("probeEleCutBasedVeto94X"  ),
                                        passingLoose94X   = cms.InputTag("probeEleCutBasedLoose94X" ),
                                        passingMedium94X  = cms.InputTag("probeEleCutBasedMedium94X"),
                                        passingTight94X   = cms.InputTag("probeEleCutBasedTight94X" ),

                                        passingVeto94XV2    = cms.InputTag("probeEleCutBasedVeto94XV2"  ),
                                        passingLoose94XV2   = cms.InputTag("probeEleCutBasedLoose94XV2" ),
                                        passingMedium94XV2  = cms.InputTag("probeEleCutBasedMedium94XV2"),
                                        passingTight94XV2   = cms.InputTag("probeEleCutBasedTight94XV2" ),

                                        passingMVA94XwpLnoiso = cms.InputTag("probeEleMVA94XwpLnoiso" ),
                                        passingMVA94Xwp90noiso = cms.InputTag("probeEleMVA94Xwp90noiso" ),
                                        passingMVA94Xwp80noiso = cms.InputTag("probeEleMVA94Xwp80noiso" ),
                                        passingMVA94XwpLiso = cms.InputTag("probeEleMVA94XwpLiso" ),
                                        passingMVA94Xwp90iso = cms.InputTag("probeEleMVA94Xwp90iso" ),
                                        passingMVA94Xwp80iso = cms.InputTag("probeEleMVA94Xwp80iso" ),
                                        passingMVA94XwpLnoisoV2 = cms.InputTag("probeEleMVA94XwpLnoisoV2" ),
                                        passingMVA94Xwp90noisoV2 = cms.InputTag("probeEleMVA94Xwp90noisoV2" ),
                                        passingMVA94Xwp80noisoV2 = cms.InputTag("probeEleMVA94Xwp80noisoV2" ),
                                        passingMVA94XwpLisoV2 = cms.InputTag("probeEleMVA94XwpLisoV2" ),
                                        passingMVA94Xwp90isoV2 = cms.InputTag("probeEleMVA94Xwp90isoV2" ),
                                        passingMVA94Xwp80isoV2 = cms.InputTag("probeEleMVA94Xwp80isoV2" ),

                                        passingMVA94XwpHZZisoV2 = cms.InputTag("probeEleMVA94XwpHZZisoV2" ),

                                        passingHEEP       = cms.InputTag("probeEleHEEP"),

                                        )
                                    )

process.tnpPhoIDs = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                    tnpVars.mcTruthCommonStuff, tnpVars.CommonStuffForPhotonProbe,
                                    tagProbePairs = cms.InputTag("tnpPairingPhoIDs"),                                                                                         
                                    probeMatches  = cms.InputTag("genProbePho"),
                                    allProbes     = cms.InputTag("probePho"),
                                    flags         = cms.PSet(
#                                        passingLoose      = cms.InputTag("probePhoCutBasedLoose"),
#                                        passingMedium     = cms.InputTag("probePhoCutBasedMedium"),
#                                        passingTight      = cms.InputTag("probePhoCutBasedTight"),
#                                        passingMVA        = cms.InputTag("probePhoMVA"),
                                         passingLoose80X   = cms.InputTag("probePhoCutBasedLoose80X"),
                                         passingMedium80X  = cms.InputTag("probePhoCutBasedMedium80X"),
                                         passingTight80X   = cms.InputTag("probePhoCutBasedTight80X"),
                                         passingMVA80Xwp90 = cms.InputTag("probePhoMVA80Xwp90"),
                                         passingMVA80Xwp80 = cms.InputTag("probePhoMVA80Xwp80"),
                                         
                                         passingLoose94X   = cms.InputTag("probePhoCutBasedLoose94X"),
                                         passingMedium94X  = cms.InputTag("probePhoCutBasedMedium94X"),
                                         passingTight94X   = cms.InputTag("probePhoCutBasedTight94X"),

                                         #passingLoose100XV2   = cms.InputTag("probePhoCutBasedLoose100XV2"),
                                         #passingMedium100XV2  = cms.InputTag("probePhoCutBasedMedium100XV2"),
                                         #passingTight100XV2   = cms.InputTag("probePhoCutBasedTight100XV2"),

                                         passingMVA94Xwp90 = cms.InputTag("probePhoMVA94Xwp90"),
                                         passingMVA94Xwp80 = cms.InputTag("probePhoMVA94Xwp80"),

                                         passingMVA94XV2wp90 = cms.InputTag("probePhoMVA94XV2wp90"),
                                         passingMVA94XV2wp80 = cms.InputTag("probePhoMVA94XV2wp80"),
                                        )
                                    )

## add pass HLT-safe flag, available for miniAOD only
if not options['useAOD'] :
    setattr( process.tnpEleTrig.flags, 'passingHLTsafe', cms.InputTag("probeEleHLTsafe" ) )
    setattr( process.tnpEleIDs.flags , 'passingHLTsafe', cms.InputTag("probeEleHLTsafe" ) )

# Add SUSY variables to the "variables", add SUSY IDs to the "flags"
if options['addSUSY'] :
    setattr( process.tnpEleIDs.variables , 'el_miniIsoChg', cms.string("userFloat('miniIsoChg')") )
    setattr( process.tnpEleIDs.variables , 'el_miniIsoAll', cms.string("userFloat('miniIsoAll')") )
    setattr( process.tnpEleIDs.variables , 'el_ptRatio', cms.string("userFloat('ptRatio')") )
    setattr( process.tnpEleIDs.variables , 'el_ptRatioUncorr', cms.string("userFloat('ptRatioUncorr')") )
    setattr( process.tnpEleIDs.variables , 'el_ptRel', cms.string("userFloat('ptRel')") )
    setattr( process.tnpEleIDs.variables , 'el_MVATTH', cms.InputTag("susyEleVarHelper:electronMVATTH") )   
    setattr( process.tnpEleIDs.variables , 'el_sip3d', cms.InputTag("susyEleVarHelper:sip3d") )
    def addFlag(name):
        setattr( process.tnpEleIDs.flags, 'passing'+name, cms.InputTag('probes'+name ) )
    from EgammaAnalysis.TnPTreeProducer.electronsExtrasSUSY_cff  import workingPoints
    for wp in workingPoints: addFlag(wp)


tnpSetup.customize( process.tnpEleTrig , options )
tnpSetup.customize( process.tnpEleIDs  , options )
tnpSetup.customize( process.tnpPhoIDs  , options )
tnpSetup.customize( process.tnpEleReco , options )


process.tree_sequence = cms.Sequence()
if (options['DoTrigger']): process.tree_sequence *= process.tnpEleTrig
if (options['DoRECO'])   : process.tree_sequence *= process.tnpEleReco
if (options['DoEleID'])  : process.tree_sequence *= process.tnpEleIDs
if (options['DoPhoID'])  : process.tree_sequence *= process.tnpPhoIDs


##########################################################################
## PATHS
##########################################################################
process.out = cms.OutputModule("PoolOutputModule", 
                               fileName = cms.untracked.string(options['OUTPUTEDMFILENAME']),
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
                               )
process.outpath = cms.EndPath(process.out)
if (not options['DEBUG']):
    process.outpath.remove(process.out)

process.evtCounter = cms.EDAnalyzer('SimpleEventCounter')

process.p = cms.Path(
        process.evtCounter        +
        process.hltFilter         +
        process.cand_sequence     + 
        process.tnpPairs_sequence +
        process.mc_sequence       +
        process.tree_sequence 
        )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(options['OUTPUT_FILE_NAME']),
    closeFileFast = cms.untracked.bool(True)
    )

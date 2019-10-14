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
    "GT","",
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
options['DataYear']             = 2018
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
#options['ELECTRON_TAG_CUTS']    = "(abs(-log(tan(superCluster.position.theta/2)))<=2.1) && !(1.4442<=abs(-log(tan(superClusterPosition.theta/2)))<=1.566) && pt >= 30.0"
#### I think we don't need sceta<2.1 for 2017
#### If not, we can apply it later
options['ELECTRON_TAG_CUTS']    = "!(1.4442<=abs(-log(tan(superClusterPosition.theta/2)))<=1.566) && pt >= 30.0"

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
    options['TnPPATHS']            = cms.vstring("HLT_Ele32_WPTight_Gsf_v*")
    options['TnPHLTTagFilters']    = cms.vstring("hltEle32WPTightGsfTrackIsoFilter")
    options['HLTFILTERTOMEASURE']  = cms.vstring("hltEle32WPTightGsfTrackIsoFilter","hltEG200HEFilter","hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['GLOBALTAG']           = '102X_upgrade2018_realistic_v18'
else:
    options['OUTPUT_FILE_NAME']    = "TnPTree_data.root"
    options['TnPPATHS']            = cms.vstring("HLT_Ele32_WPTight_Gsf_v*") 
    options['TnPHLTTagFilters']    = cms.vstring("hltEle32WPTightGsfTrackIsoFilter")
    options['HLTFILTERTOMEASURE']  = cms.vstring("hltEle32WPTightGsfTrackIsoFilter","hltEG200HEFilter","hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter")
    options['TnPHLTProbeFilters']  = cms.vstring()
    options['GLOBALTAG']           = '102X_dataRun2_Sep2018ABC_v2' ## for periodD prompt, '102X_dataRun2_Prompt_v13'

if varOptions.GT != "" :
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
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/11697BCC-C4AB-204B-91A9-87F952F9F2C6.root',
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/B4E2FD3B-1B1A-5A43-A4A5-09F81DBC6467.root',
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/B95B0723-63AF-5C4B-BE74-8042562F684A.root',
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/7E268FF6-C75F-034E-B3FE-547E23A1F96E.root',
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/B30A6292-160E-5C44-9961-67C774C35A10.root',
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/8939432A-9601-6A45-8D03-68FCC77DEB24.root',
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/AD4C6C13-03D5-9B4C-82F5-F9903C59E0F1.root',
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/73C922DB-1246-754E-84FB-DA709C321F07.root',
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/63BDE40E-08E1-6A46-926D-A2EC232CC3D6.root',
'/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/1A2556C2-5EC8-9049-B8E7-65836677D338.root',
)
if varOptions.isMC:
  options['INPUT_FILE_NAME'] = cms.untracked.vstring(
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/042C8EE9-9431-5443-88C8-77F1D910B3A5.root',
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/08540B6C-AA39-3F49-8FFE-8771AD2A8885.root',
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/8F103C41-A7BA-754F-923E-B5C102366249.root',
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/973AE986-070D-ED40-9A99-393E4E212670.root',
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/872CAA08-8946-AA43-A812-F6E5963D917B.root',
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/604CA4EB-3C6D-D849-A471-72BDB283F4CE.root',
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/E9AD25CE-9067-374E-B1E4-490979612D55.root',
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/92B45BD0-8EFF-D443-B015-2C328C0E70A8.root',
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/B4A14390-2C2E-9847-B1F5-84BE85DEEC54.root',
'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/07E16CC8-CEDE-B94C-A68F-976E312C37D9.root',
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
print '@@@@ GT = '+str(options['GLOBALTAG'])
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

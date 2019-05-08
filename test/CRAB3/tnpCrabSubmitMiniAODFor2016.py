from CRABClient.UserUtilities import config, getUsernameFromSiteDB

# this will use CRAB client API
from CRABAPI.RawCommand import crabCommand

# this now standard CRAB configuration

from WMCore.Configuration import Configuration

config = config()

submitVersion = "2016"
config.General.transferLogs = False
config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = '/afs/cern.ch/work/j/jskim/EGammaTnP/For2016_CMSSW_10_2_5/src/EgammaAnalysis/TnPTreeProducer/test/TnPTreeProducerFor2016_cfg.py'
config.JobType.sendExternalFolder     = True

config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/EGammaTnPNtuple_Ele115Added/2016/' % (getUsernameFromSiteDB())
config.Site.storageSite = 'T2_KR_KNU'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_%s' % submitVersion

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    ##### submit MC

    config.JobType.pyCfgParams  = ['doTrigger=True', 'doEleID=True', 'isMC=True', 'doPhoID=False']
    config.Data.splitting = 'FileBased'
    config.General.requestName  = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM'
    config.Data.unitsPerJob = 1
    submit(config)


    config.General.requestName  = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM'
    submit(config)

    # /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSI
    # /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM
    # /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM


    ##### now submit DATA

    config.Data.splitting     = 'LumiBased'
    config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
    config.Data.unitsPerJob   = 100
    config.JobType.pyCfgParams  = ['doTrigger=True', 'doEleID=True', 'isMC=False', 'doPhoID=False']
 
    config.General.requestName  = 'Legacy2016_RunB'
    config.Data.inputDataset    = '/SingleElectron/Run2016B-17Jul2018_ver2-v1/MINIAOD'
    submit(config)    

    config.General.requestName  = 'Legacy2016_RunC'
    config.Data.inputDataset    = '/SingleElectron/Run2016C-17Jul2018-v1/MINIAOD'
    submit(config)  

    config.General.requestName  = 'Legacy2016_RunD'
    config.Data.inputDataset    = '/SingleElectron/Run2016D-17Jul2018-v1/MINIAOD'
    submit(config)  

    config.General.requestName  = 'Legacy2016_RunE'
    config.Data.inputDataset    = '/SingleElectron/Run2016E-17Jul2018-v1/MINIAOD'
    submit(config)  

    config.General.requestName  = 'Legacy2016_RunF'
    config.Data.inputDataset    = '/SingleElectron/Run2016F-17Jul2018-v1/MINIAOD'
    submit(config)  

    config.General.requestName  = 'Legacy2016_RunG'
    config.Data.inputDataset    = '/SingleElectron/Run2016G-17Jul2018-v1/MINIAOD'
    submit(config)  

    config.General.requestName  = 'Legacy2016_RunH'
    config.Data.inputDataset    = '/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD'
    submit(config)  

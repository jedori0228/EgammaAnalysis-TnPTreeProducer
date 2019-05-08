from CRABClient.UserUtilities import config, getUsernameFromSiteDB

# this will use CRAB client API
from CRABAPI.RawCommand import crabCommand

# this now standard CRAB configuration

from WMCore.Configuration import Configuration

config = config()

submitVersion = "2017"
config.General.transferLogs = False
config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = '/afs/cern.ch/work/j/jskim/EGammaTnP/For2016_CMSSW_10_2_5/src/EgammaAnalysis/TnPTreeProducer/test/TnPTreeProducerFor2017_cfg.py'
config.JobType.sendExternalFolder     = True

config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/EGammaTnPNtuple_Ele115Added/2017/' % (getUsernameFromSiteDB())
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
    config.General.requestName  = 'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
    config.Data.unitsPerJob = 1
    #submit(config)

    config.General.requestName  = 'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext1'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
    #submit(config)

    config.General.requestName  = 'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8'
    config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
    #submit(config)

    ##### now submit DATA

    config.Data.splitting     = 'LumiBased'
    config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
    config.Data.unitsPerJob   = 100
    config.JobType.pyCfgParams  = ['doTrigger=True', 'doEleID=True', 'isMC=False', 'doPhoID=False']

    periods = ["B", "C", "D", "E", "F"]
    for period in periods:
      config.General.requestName  = 'Run2017_period'+period
      config.Data.inputDataset    = '/SingleElectron/Run2017'+period+'-31Mar2018-v1/MINIAOD'
      submit(config)




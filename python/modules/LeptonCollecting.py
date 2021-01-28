import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class LeptonCollecting(Module):

    def __init__(
        self,
        tightMuonsCollection = lambda event: Collection(event, "Muon"),
        tightElectronsCollection = lambda event: Collection(event, "Electron"),
        looseMuonCollection = lambda event: Collection(event, "Muon"),
        looseElectronCollection = lambda event: Collection(event, "Electron"),
        outputName = "Leptons",
        globalOptions={"isData": False, "year": 2016},
        storeKinematics=["pt", "eta", "phi", "charge", "isMuon", "isElectron", "relIso", "dxy", "dz", 'dxysig', 'dzsig']
    ):

        self.globalOptions = globalOptions
        self.tightMuonsCollection = tightMuonsCollection
        self.tightElectronsCollection = tightElectronsCollection
        self.looseMuonCollection = looseMuonCollection
        self.looseElectronCollection = looseElectronCollection
        self.outputName = outputName
        self.storeKinematics = storeKinematics

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("nleading"+self.outputName, "I")
        self.out.branch("nsubleading"+self.outputName, "I")

        self.out.branch(self.outputName+"_muonmuon", "I")
        self.out.branch(self.outputName+"_electronelectron", "I")
        self.out.branch(self.outputName+"_muonelectron", "I")
        self.out.branch(self.outputName+"_electronmuon", "I")
        self.out.branch(self.outputName+"_muonjets", "I")
        self.out.branch(self.outputName+"_electronjets", "I")

        for variable in self.storeKinematics:
            self.out.branch("leading"+self.outputName+"_"+variable, "F", lenVar="nleading"+self.outputName)
            self.out.branch("subleading"+self.outputName+"_"+variable, "F", lenVar="nsubleading"+self.outputName)

        #for variable in self.storeKinematics:
            #self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        tightMuons = self.tightMuonsCollection(event)
        tightElectrons = self.tightElectronsCollection(event)

        looseMuons = self.looseMuonCollection(event)
        looseElectrons = self.looseElectronCollection(event)


        for lepton in tightMuons+looseMuons:
            lepton.isMuon = 1
            lepton.isElectron = 0
            lepton.relIso = lepton.pfRelIso04_all

        for lepton in tightElectrons+looseElectrons:
            lepton.isMuon = 0
            lepton.isElectron = 1
            lepton.relIso = lepton.pfRelIso03_all

        tightLeptons = []
        looseLeptons = []

        tightLeptons = tightMuons+tightElectrons
        looseLeptons = looseMuons+looseElectrons

        tightLeptons = sorted(tightLeptons, key=lambda x: x.pt, reverse=True)
        # select leading only, move subleading to "loose"
        looseLeptons.extend(tightLeptons[1:])

        if len(tightLeptons) > 0:
            tightLeptons = [tightLeptons[0]]
        else:
            tightLeptons = []
        looseLeptons = sorted(looseLeptons, key=lambda x: x.pt, reverse=True)

        muonmuon = 0
        electronelectron = 0
        muonelectron = 0
        electronmuon = 0
        muonjets = 0
        electronjets = 0


        ## flavour categorisation :
        if len(tightLeptons) > 0 and len(looseLeptons) > 0:
            # Ensure pt(l1) > pt(l2)
            if tightLeptons[0].pt < looseLeptons[0].pt:
                return False
                
            if tightLeptons[0].isMuon and looseLeptons[0].isMuon:
                muonmuon = 1

            elif tightLeptons[0].isElectron and looseLeptons[0].isElectron:
                electronelectron= 1

            elif tightLeptons[0].isMuon and looseLeptons[0].isElectron:
                muonelectron = 1
               
            elif tightLeptons[0].isElectron and looseLeptons[0].isMuon:
                electronmuon = 1

        elif len(tightLeptons) > 0:
            if tightLeptons[0].isMuon:
                muonjets = 1
            elif tightLeptons[0].isElectron:
                electronjets = 1

        if muonmuon or muonelectron or muonjets:
            if event.IsoMuTrigger_flag:
                setattr(event, "isTriggered", 1)
            else:
                setattr(event, "isTriggered", 0)
        elif electronelectron or electronmuon or electronjets:
            if event.IsoElectronTrigger_flag:
                setattr(event, "isTriggered", 1)
            else:
                setattr(event, "isTriggered", 0)
        else:
            setattr(event, "isTriggered", 0)

        for lepton in tightLeptons+looseLeptons:
            if lepton.dxyErr < 1e-6:
                 setattr(lepton, "dxysig", -1.)
            else:
                 setattr(lepton, "dxysig", math.fabs(lepton.dxy)/math.fabs(lepton.dxyErr))

            if lepton.dzErr < 1e-6:
                 setattr(lepton, "dzsig", -1.)
            else:
                 setattr(lepton, "dzsig", math.fabs(lepton.dz)/math.fabs(lepton.dzErr))

        self.out.fillBranch("nleading"+self.outputName, len(tightLeptons))
        self.out.fillBranch("nsubleading"+self.outputName, len(looseLeptons))

        for variable in self.storeKinematics:
            self.out.fillBranch("leading"+self.outputName+"_"+variable,map(lambda lepton: getattr(lepton,variable),tightLeptons))
            self.out.fillBranch("subleading"+self.outputName+"_"+variable,map(lambda lepton: getattr(lepton,variable),looseLeptons))

        self.out.fillBranch(self.outputName+"_muonmuon", muonmuon)
        self.out.fillBranch(self.outputName+"_electronelectron", electronelectron)
        self.out.fillBranch(self.outputName+"_muonelectron", muonelectron)
        self.out.fillBranch(self.outputName+"_electronmuon", electronmuon)
        self.out.fillBranch(self.outputName+"_muonjets", muonjets)
        self.out.fillBranch(self.outputName+"_electronjets", electronjets)

        setattr(event, "leading"+self.outputName, tightLeptons)
        setattr(event, "subleading"+self.outputName, looseLeptons)

        return True

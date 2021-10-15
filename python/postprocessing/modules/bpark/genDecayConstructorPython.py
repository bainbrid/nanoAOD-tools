import ROOT
import numpy as np
import itertools
import os
from collections import Counter
from time import sleep


ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
_rootLeafType2rootBranchType = { 'UChar_t':'b', 'Char_t':'B', 'UInt_t':'i', 'Int_t':'I', 'Float_t':'F', 'Double_t':'D', 'ULong64_t':'l', 'Long64_t':'L', 'Bool_t':'O' }

class genDecayConstructorPython(Module):

    def __init__(self,momPdgId,daughtersPdgId,outputMomColl,outputDaughterColls,intermediateDecay=[],trgMuonPtEtaThresholds=[],selectTrgMuon=False,excludeTrgMuon=False):
        #super(genDecayConstructor,self).__init__()
        self.momPdgId= momPdgId
        self.daughtersPdgId = daughtersPdgId
        self.branches = ["pdgId","pt","eta","phi","mass","charge"]                
        self.outputDaughterColls=outputDaughterColls
        self.outputMomColl = outputMomColl
        self.interDecay=intermediateDecay
        self.trgMuonPtEta = trgMuonPtEtaThresholds
        self.interMom=[]
        self.selectTrgMuon=selectTrgMuon
        self.excludeTrgMuon=excludeTrgMuon
        if len(self.trgMuonPtEta):
          print "WARNING! evts without a gen mu filling trg requirments WILL be skipped."
          sleep(2)
        if self.selectTrgMuon and self.excludeTrgMuon:
           raise Exception ("cannot select and exclude trg mu simultaneously. Result propbably invalid")
        pass

    def beginJob(self):
        print self.interDecay 
        for dpdgId in self.daughtersPdgId:
          momChain = 0
          print "daughter",dpdgId
          for decay in self.interDecay:
             print "decay",decay
             decay = decay.split("->")
             mom = int( decay[0] ) 
             daughters = map(int,decay[1].split(","))
             for daughter in daughters:
               if int(daughter) != dpdgId: continue
               momChain = int(mom)
          self.interMom.append(momChain)
        print "moms",self.interMom
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for br in self.branches:
          self.out.branch("%s_%s"%(self.outputMomColl,br),'F')
          for daughter in self.outputDaughterColls:
            self.out.branch("%s_%s"%(daughter,br), 'F')
        pass 


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


     
    def filterBranchNames(self,branches,collection):
        out = []
        for br in branches:
            name = br.GetName()
            if not name.startswith(collection+'_'): continue
            out.append(name.replace(collection+'_',''))
            self.branchType[out[-1]] = br.FindLeaf(br.GetName()).GetTypeName()
        return out
   
    

    def analyze(self, event):
      """process event, return True (go to next module) or False (fail, go to next event)"""
        
      genparts = Collection(event,"GenPart")
      
      PdgIdDaughters =[]
      BidxDaughters = []
      
      trgMuons=[]
      
      #trigger threshold trgMuonPtEta Pt and Eta. User provides that
      if len(self.trgMuonPtEta)==1:
        trgMuons = [genObj for genObj in genparts if abs(getattr(genObj,"pdgId")) ==13  and getattr(genObj,"pt") > self.trgMuonPtEta[0]  ]
        if len(trgMuons)==0: 
            return False
      if len(self.trgMuonPtEta)==2:
        trgMuons = [genObj for genObj in genparts if abs(getattr(genObj,"pdgId")) ==13  and getattr(genObj,"pt") > self.trgMuonPtEta[0] and abs(getattr(genObj,"eta")) < self.trgMuonPtEta[1]  ]
        if len(trgMuons)==0: 
            return False
      for genObj in genparts:
        objPdgId =  getattr(genObj,"pdgId")
        match = False
        matchConj = False
        chainMom = 0;
        for idau,dauPdgId in enumerate(self.daughtersPdgId):
          if objPdgId == dauPdgId: 
             match = True
             chainMom = self.interMom[idau]
          if objPdgId == -1*dauPdgId: 
             matchConj = True
             chainMom = self.interMom[idau]
        
         
        if not match and not matchConj:
           continue    
        
        momIdx =  getattr(genObj,"genPartIdxMother")
        # check if the mom actually exists
        if momIdx<0 or momIdx > len(genparts):
           continue

        # get gen particlle mom pdg
        tempMomPdgId = getattr(genparts[momIdx],"pdgId")
        if chainMom==0:
          if ( match and tempMomPdgId == self.momPdgId ) or ( matchConj and tempMomPdgId == -1*self.momPdgId ):
            PdgIdDaughters.append(genObj)
            BidxDaughters.append(momIdx)
        else: 
            if abs(chainMom) !=  abs( tempMomPdgId): continue;
            grandMomIdx = getattr(genparts[momIdx],"genPartIdxMother")
            if grandMomIdx<0 or grandMomIdx > len(genparts):
               continue
            grandMomPdgId = getattr(genparts[grandMomIdx],"pdgId")
            if abs(grandMomPdgId) == abs(self.momPdgId):
              PdgIdDaughters.append(genObj)
              BidxDaughters.append(grandMomIdx)

      if len(BidxDaughters)<len(self.daughtersPdgId):
         return False
      realBidx = -1 
      uniques= Counter(BidxDaughters) 
      for idx,occur in  uniques.items(): 
         if occur == len(self.daughtersPdgId): 
            realBidx = int(idx) 
      if realBidx<0: 
        return False;

      # get final particles
      finalStateParts = [ daughter for Bidx,daughter in zip(BidxDaughters, PdgIdDaughters) if realBidx == Bidx ]
      
      realB = genparts[realBidx]
     
      #remove or select triggering muon     
      if len(trgMuons)>0:        
        finalmuons = [ part for part in finalStateParts if abs(getattr(part,"pdgId"))==13]
        if self.excludeTrgMuon and all(elem in finalmuons for elem in trgMuons):
           return False
        if self.selectTrgMuon and not( any(elem in finalmuons for elem in trgMuons) ):
           return False    


      for br in self.branches: 
        if br != "charge":   
          self.out.fillBranch("%s_%s"%(self.outputMomColl,br), getattr(realB,br) )
      else: 
          self.out.fillBranch("%s_%s"%(self.outputMomColl,br),getattr(realB,"pdgId")/abs(getattr(realB,"pdgId"))  )
        

      signFactor = getattr(realB,"pdgId") / self.momPdgId

      for part in finalStateParts:
        try:
          bridx = self.daughtersPdgId.index(signFactor*getattr(part,"pdgId"))
        except ValueError:
          reutrn False

        for br in self.branches:
           if br != "charge":
             self.out.fillBranch("%s_%s"%(self.outputDaughterColls[bridx],br),getattr(part,br) )
           else :
              if abs(part.pdgId)!=11 and abs(part.pdgId)!=13 and abs(part.pdgId)!=15:
                 self.out.fillBranch("%s_%s"%(self.outputDaughterColls[bridx],br),-1*part.pdgId/abs(part.pdgId) )
              else:
                self.out.fillBranch("%s_%s"%(self.outputDaughterColls[bridx],br),part.pdgId/abs(part.pdgId) )

      return True
      


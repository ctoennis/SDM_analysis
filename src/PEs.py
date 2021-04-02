import numpy as np
import ingredients as ING
import Likelihoods as LH

class PseudoExperiment():

    'This class will contain a single pseudo experiment. Containiing functions generate signal and background events according to histograms provided to them'

    def __init__(self):

        self.events=[]
        self.sigs=[]
        self.bgs=[]
        self.nstrue=0
        self.nbgtrue=0
        self.ntot=0

    def GenBG(self,ingredient,energy,n_bg,test = ING.H1D.Empty(-2.0,3.0,50),test_ene = ING.H1D.Empty(-2.0,3.0,50),tio=False):

        'This function generates background events. The first two arguments provide the angular and energy distribution of background events. The third argument should contain the signal energy distribution to avoid errors in TS calculation'

        for i in range(int(n_bg)):

            ran = ingredient.GetRandom()
            ene = energy.GetRandom()

            if tio:

                while test.GetContent(ran) == 0.0:

                    ran = ingredient.GetRandom()

                while test_ene.GetContent(ene) == 0.0:

                    ene = energy.GetRandom()

#            print ran

            self.events.append([ran,ene])
            self.bgs.append(self.ntot+i)

        self.nbgtrue+=int(n_bg)
        self.ntot+=int(n_bg)        

    def GenSG(self,ingredient,energy,n_sg,test = ING.H1D.Empty(-2.0,3.0,50),test_ene = ING.H1D.Empty(-2.0,3.0,50),tio=False):

        'This function generates signal events. The first two arguments provide the angular and energy distribution of signal events. The third argument should contain the background energy distribution to avoid errors in TS calculation'

        for i in range(int(n_sg)):

            ran = ingredient.GetRandom()
            ene = energy.GetRandom()

            if tio:

                while test.GetContent(ran) == 0.0:

                    ran = ingredient.GetRandom()

                while test_ene.GetContent(ene) == 0.0:

                    ene = energy.GetRandom()

#            print ene

            self.events.append([ran,ene])
            self.sigs.append(self.ntot+i)
            
        self.nstrue+=n_sg
        self.ntot+=n_sg


    def GenEvt(self,event):

        'This function adds an event from a given event set. Used fro unblinding'

        self.events.append(event)
        self.bgs.append(len(self.bgs))
        self.nbgtrue += 1
        self.ntot    += 1

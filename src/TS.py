import numpy as np
import ingredients as ING
import Likelihoods as LH
import PEs as PE
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar

class TSdist():

    def __init__(self,s_signal, s_ene_signal, s_background, s_ene_background, signal, ene_signal, background, ene_background, N = 100, out = 0.0):

        self.TS             = ING.H1D.Empty(-10.0,5.0,300)
        self.NS             = ING.H1D.Empty(0.0,50.0,100)
        self.signal         = signal
        self.background     = background
        self.ene_signal     = ene_signal
        self.ene_background = ene_background
        self.PSF            = s_signal
        self.BPSF           = s_background
        self.ENE_signal     = s_ene_signal
        self.ENE_background = s_ene_background
        self.N              = N
        self.out            = out
        self.unblind        = 0.0

    def AddPE(self,n):

        it = PE.PseudoExperiment()
        it.GenBG(self.background,self.ene_background,self.N)
        it.GenSG(self.signal,self.ene_signal,n,tio=True,test=self.background,test_ene=self.ene_background)
        Lik = LH.StandardLH(it,self.PSF,self.BPSF,self.ENE_signal,self.ENE_background,self.out)
        res = minimize_scalar(Lik,bounds=(0.0,100.0), method = "bounded")
#        print(res.x)
        print n, res.x
#        print((Lik(10.0)))
        fin = Lik(res.x)

        print res.x , (Lik(0)-fin)

        if fin < Lik(0):

            self.TS.Fill(np.log10(Lik(0)-fin))

#            print np.log10(Lik(0)-fin), res.x

        else: 
            
            self.TS.Fill(-10.0)
            
#            print -5.0, Lik(0), fin

        self.NS.Fill(res.x)

    def Unblind(self,dist):

        'get the unblinded TS value'

        it = PE.PseudoExperiment()

        for event in dist:

            it.GenEvt(event)

        Lik = LH.StandardLH(it,self.PSF,self.BPSF,self.ENE_signal,self.ENE_background,self.out)
        res = minimize_scalar(Lik,bounds=(0.0,100.0), method = "bounded")

        fin = Lik(res.x)

        print res.x , (Lik(0)-fin)

        if fin < Lik(0):

            self.unblind = (np.log10(Lik(0)-fin))

#            print np.log10(Lik(0)-fin), res.x

        else:

            self.unblind = -10 
        
class TSdist_twosample():

    def __init__(self, s_signal1, s_signal2, ENE_signal1, ENE_signal2, s_background1, s_background2, ENE_background1, ENE_background2, signal1, signal2, ene_signal1, ene_signal2, background1, background2, ene_background1, ene_background2, w1, w2, N1, N2, out1 = 0.0, out2 = 0.0):

        self.TS1             = ING.H1D.Empty(-5.0,5.0,100)
        self.TS2             = ING.H1D.Empty(-5.0,5.0,100)
        self.NS1             = ING.H1D.Empty(0.0,50.0,100)
        self.NS2             = ING.H1D.Empty(0.0,50.0,100)
        self.signal1         = signal1
        self.signal2         = signal2
        self.ene_signal1     = ene_signal1
        self.ene_signal2     = ene_signal2
        self.background1     = background1
        self.background2     = background2
        self.ene_background1 = ene_background1
        self.ene_background2 = ene_background2
        self.PSF1            = s_signal1
        self.ENE_signal1     = ENE_signal1
        self.BPSF1           = s_background1
        self.ENE_background1 = ENE_background1
        self.PSF2            = s_signal2
        self.ENE_signal2     = ENE_signal2
        self.BPSF2           = s_background2
        self.ENE_background2 = ENE_background2
        self.N1              = N1
        self.N2              = N2
        self.out1            = out1
        self.out2            = out2
        self.w1              = w1
        self.w2              = w2

    def AddPE(self,n,pois=False):

        it1 = PE.PseudoExperiment()
        it1.GenBG(self.background1,self.ene_background1,self.N1)

        if pois:

            it1.GenSG(self.signal1,self.ene_signal1,np.random.poisson(n*self.w1/(self.w1+self.w2)),test=self.background1,tio=True)

        else:

            #it1.GenSG(self.signal1,n*self.w1/(self.w1+self.w2),test=self.background1,tio=True)
            it1.GenSG(self.signal1,self.ene_signal1,n,test=self.background1,tio=True)

#        print it1.nbgtrue 
#        print it1.nstrue

        it2 = PE.PseudoExperiment()
        it2.GenBG(self.background2,self.ene_background2,self.N2)

        if pois:

            it2.GenSG(self.signal2,self.ene_signal2,np.random.poisson(self.w2/(self.w1+self.w2) * n),test=self.background2,tio=True)

        else:

            #it2.GenSG(self.signal2,self.w2/(self.w1+self.w2) * n,test=self.background2,tio=True)
            it2.GenSG(self.signal2,self.ene_signal2,n,test=self.background2,tio=True)
            

#        print 'now bg'

#        print it2.nbgtrue 
#        print it2.nstrue

        Lik1 = LH.StandardLH(it1,self.PSF1,self.BPSF1,self.ENE_signal1,self.ENE_background1,cone=0.9)
        Lik2 = LH.StandardLH(it2,self.PSF2,self.BPSF2,self.ENE_signal2,self.ENE_background2,cone=0.0)
        res1 = minimize_scalar(Lik1,bounds=(0.0,100),method='bounded')
        res2 = minimize_scalar(Lik2,bounds=(0.0,100),method='bounded')
        print res1.x
        print res2.x
        #print(Lik(res.x))
#        print(Lik([0.0,0.0]))
#        print((Lik(0.0)-Lik(res.x)))
        self.TS1.Fill((Lik1(0.0)-Lik1(res1.x)))
        self.TS2.Fill((Lik2(0.0)-Lik2(res2.x)))
#        self.NS.Fill(res.x)

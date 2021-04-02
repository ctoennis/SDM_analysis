import sys
import os
import numpy as np
import ingredients as ing
import PEs as pe

def StandardLH(data,signal,backg,ene_signal,ene_backg,out,cone=-1.0):

    'Standard likelihood function'
    
    def fx(n_s):

        lk = 0.0

        ntot = len(data.events)
        
        for event in data.events:

            if event[0] > cone:

                lk -= np.log10(n_s/(float(ntot)+n_s) * signal.GetContent(event[0]) * ene_signal.GetContent(event[1]) + float(ntot)/(float(ntot)+n_s) * backg.GetContent(event[0]) * ene_backg.GetContent(event[1]))

        return lk

    return fx

def ExtendedLH(data,signal,backg,ene_signal,ene_backg,out,cone=-1.0):

    'the extended likelihood'

    def fx(n_s):

        lk = 0.

        ntot = len(data.events)

        for event in data.events:

            lk -= np.log10(n_s * signal.GetContent(event[0]) * ene_signal.GetContent(event[1]) + float(ntot) * backg.GetContent(event[0]) * ene_backg.GetContent(event[1])) + n_s 

        return lk

    return fx

def Twosample_Standard(data1,data2,w1,w2,signal1,signal2,backg1,backg2,out1=0.0,out2=0.0):

    'likelihood to combine two samples'

    def fx(n_s):

        lk = 0.

        ntot1 = len(data1.events)
        ntot2 = len(data2.events) 
        N_s   = float(n_s)
        N_s1  = N_s#*(w1/(w1+w2))
        N_s2  = N_s#*(w2/(w1+w2))

        print n_s

        for event in data2.events:

            if event > -1.0:

                lk -= np.log10(N_s2/(float(ntot2)+N_s2) * 1.0 * signal2.GetContent(event) + float(ntot2)/(float(ntot2)+N_s2) * backg2.GetContent(event))

            else:

                lk -= (np.log10(float(ntot2)/(float(ntot2)+N_s2) * backg2.GetContent(event)))

        return lk

    return fx

                

            

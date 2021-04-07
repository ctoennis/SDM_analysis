#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00


import os
import src/ING as ING
import matplotlib as mlp
import numpy as np
mlp.use('Agg')
import matplotlib.pyplot as plt

mass = [100,250,350,500,750,1000,2500,5000,7500,10000,25000,50000,75000]
channel = [5,8,11,13]
life = [1, 10, 100, 1000, 10000]
medmass = [1, 10, 100, 1000, 10000]

masslist     = [100, 250, 350, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000]
medmasslist  = [1, 10, 100, 1000, 10000]
lifetimelist = [1, 10, 100, 1000, 10000]
channellist  = [5, 8, 11, 13]

ACC     = []
NS      = []
FLUX    = []
SI      = []
SD      = []
nlim    = []
dNdE_e   = []
dNdE_ae  = []
dNdE_mu  = []
dNdE_amu = []
N_e      = []
N_ae     = []
N_mu     = []
N_amu    = []

fluxconv2   = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/conversion_darksusy/SDCS.txt")

for Mass in mass: #[100,200,350,500,750,1000,2500,5000,7500,10000,25000,50000,75000,100000]:

    a1= []
    a2= []
    a3= []
    a4= []
    aa1= []
    aa2= []
    aa3= []
    aa4= []

    for Channel in channel: #[5,8,11,13]:

        b1= []
        b2= []
        b3= []
        b4= []
        bb1= []
        bb2= []
        bb3= []
        bb4= []

        for Life in life: # [0.00042,0.042,0.42,4.2]:

            c1= []
            c2= []
            c3= []
            c4= []
            cc1= []
            cc2= []
            cc3= []
            cc4= []

            for Medmass in medmass: #[10,100,1000,10000]:

                if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(Mass) + "-ch" + str(Channel) + "-x" + str(Medmass) + "-l" + str(Life) + "-009080-000000-read-earth_mu.dat"):
		    
                    c1.append( ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(Mass) + "-ch" + str(Channel) + "-x" + str(Medmass) + "-l" + str(Life) + "-009080-000000-read-earth_e.dat"))

		    c2.append( ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(Mass) + "-ch" + str(Channel) + "-x" + str(Medmass) + "-l" + str(Life) + "-009080-000000-read-earth_ae.dat"))

		    c3.append( ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(Mass) + "-ch" + str(Channel) + "-x" + str(Medmass) + "-l" + str(Life) + "-009080-000000-read-earth_mu.dat"))

		    c4.append( ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/spectra_wimpsim/wa-sunsum-m" + str(Mass) + "-ch" + str(Channel) + "-x" + str(Medmass) + "-l" + str(Life) + "-009080-000000-read-earth_amu.dat"))
                    cc1.append(c3[-1].Integral(c3[-1].points[0][0],c3[-1].points[-1][0]))
                    cc2.append(c4[-1].Integral(c4[-1].points[0][0],c4[-1].points[-1][0]))
                    cc3.append(sum([x[1] for x in c3[-1].points])/200.)
                    cc4.append(sum([x[1] for x in c4[-1].points])/200.)

                else:

                    c1.append("none")
                    c2.append("none")
                    c3.append("none")
                    c4.append("none")
                    cc1.append("none")
                    cc2.append("none")
                    cc3.append("none")
                    cc4.append("none")

            b1.append(c1)
            b2.append(c2)
            b3.append(c3)
            b4.append(c4)
            bb1.append(cc1)
            bb2.append(cc2)
            bb3.append(cc3)
            bb4.append(cc4)

        a1.append(b1)
        a2.append(b2)
        a3.append(b3)
        a4.append(b4)
        aa1.append(bb1)
        aa2.append(bb2)
        aa3.append(bb3)
        aa4.append(bb4)

    dNdE_e.append(a1)
    dNdE_ae.append(a2)
    dNdE_mu.append(a3)
    dNdE_amu.append(a4)
    N_e.append(aa1)
    N_ae.append(aa2)
    N_mu.append(aa3)
    N_amu.append(aa4)


for Channel in range(len(channel)):

    nlom = []
    CCC  = []

    for Life in range(len(life)):

        nlum = []
	BCC  = []

        for Medmass in range(len(medmass)):

	    if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/new_sens/ns_limit-med" + str(Medmass) + "-gl" + str(Life) + "-ch" + str(Channel) + ".txt"):

		nlum.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/new_sens/ns_limit-med" + str(Medmass) + "-gl" + str(Life) + "-ch" + str(Channel) + ".txt"))

	    else:

		nlum.append("none")

	    if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Acceptance/rev_med" + str(medmass[Medmass]) + "_lifetime" + str(life[Life]) + "_channel" + str(channel[Channel]) + "_track.txt"):
		
		BCC.append(ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Acceptance/rev_med" + str(medmass[Medmass]) + "_lifetime" + str(life[Life]) + "_channel" + str(channel[Channel]) + "_track.txt"))

		

	    else:

		BCC.append("none")
	      
	nlom.append(nlum)
	CCC.append(BCC)

    nlim.append(nlom)
    ACC.append(CCC)


ardid      = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/mardid.txt")

antares    = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/ANTARES.txt")

SDCS_ardid = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/mardid_SDCS.txt")

SICS_ardid = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/mardid_SICS.txt")

SDCS_HAWC  = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/HAWC.txt")

SDCS_HAWC_tau  = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/HAWC_tau.txt")

SDCS_IC_bb = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/SDCS_IC_bb_new.txt")
SDCS_IC_WW = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/SDCS_IC_WW_new.txt")
SDCS_IC_tt = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/SDCS_IC_tt_new.txt")

AU = 149597870.71

plt.rc('font', size=24)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rc('text', usetex=True)



for Channel in range(4):

    FLUX_a = []
    SD_a   = []

    for Life in range(5):

        FLUX_b = []
        SD_b   = []

        for Medmass in range(4):

	    FLUX_p = []
	    SD_p = []
  
	    if ACC[Channel][Life][Medmass] == "none":

		
		for Mass in range(len(mass)):

		    FLUX_p.append([mass[Mass],1])

                    SD_p.append([mass[Mass],1])

		FLUX_i = ING.Graph(FLUX_p)
		SD_i   = ING.Graph(SD_p)

		FLUX_b.append(FLUX_i)
		SD_b.append(SD_i)

		continue

            if nlim[Channel][Life][Medmass] == "none":

		for Mass in range(len(mass)):

                    FLUX_p.append([mass[Mass],1])

                    SD_p.append([mass[Mass],1])

                FLUX_i = ING.Graph(FLUX_p)
                SD_i   = ING.Graph(SD_p)

                FLUX_b.append(FLUX_i)
                SD_b.append(SD_i)

		continue
		

	    for Mass in range(len(mass)):



	        if N_mu[Mass][Channel][Life][Medmass] == "none" or (N_mu[Mass][Channel][Life][Medmass]+N_amu[Mass][Channel][Life][Medmass]) == 0.0:


		    continue


		if ACC[Channel][Life][Medmass].Eval(mass[Mass]) == 0.0:

		    continue

		flux = nlim[Channel][Life][Medmass].Eval(mass[Mass])*2*365.0*86400*1000000.0/(ACC[Channel][Life][Medmass].Eval(mass[Mass])) #the factor of 2 accounts for the treatment of neutrinos and antineutrinos 

		capt = flux*4.0*np.pi*np.power(AU,2)/((N_mu[Mass][Channel][Life][Medmass]+N_amu[Mass][Channel][Life][Medmass])*365.0*86400) #capture rate 

	        FLUX_p.append([mass[Mass],flux])
		
		SD_p.append([mass[Mass],capt/fluxconv2.Eval(mass[Mass])])

	    FLUX_i = ING.Graph(FLUX_p)
	    SD_i   = ING.Graph(SD_p)

	    FLUX_b.append(FLUX_i)
	    SD_b.append(SD_i)

	FLUX_a.append(FLUX_b)
	SD_a.append(SD_b)

    FLUX.append(FLUX_a)
    SD.append(SD_a)

fig, ax = plt.subplots()


plt.loglog([x[0] for x in FLUX[0][3][1].points],[x[1] for x in FLUX[0][3][1].points],"purple",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3, bb',linewidth=1.5)
plt.loglog([x[0] for x in FLUX[1][4][2].points],[x[1] for x in FLUX[1][4][2].points],"k",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3 WW',linewidth=1.5)
plt.loglog([x[0] for x in FLUX[2][4][1].points],[x[1] for x in FLUX[2][4][1].points],"r",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3 tt',linewidth=1.5)
plt.loglog([x[0] for x in FLUX[3][4][1].points],[x[1] for x in FLUX[3][4][1].points],"darkgreen",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3 nunu',linewidth=1.5)

ax.tick_params(axis='both', length=12, width=1)

ax.tick_params(axis='both', which='minor', length=9)
plt.title(r"Flux sensitivity")
plt.ylim((1e5,1e10))
plt.xlim((1e2,1e5))
legend = plt.legend(loc='upper right')#,bbox_to_anchor=(1.1, 1.0))

for label in legend.get_lines():
    label.set_linewidth(1.5)

#plt.text(230.0,5e5,"IceCube Preliminary",color="r")

plt.ylabel(r"Flux [Yr$^{-1}$km$^{-2}]$")
plt.xlabel(r"M$_{\rm DM}$ [GeV]")
plt.savefig("/home/ctoennis/analyses/standard_analysis_framework/flux_new.png",format="png")

plt.clf()

fig, ax = plt.subplots()

#print(SD[0])


plt.loglog([x[0] for x in SDCS_ardid.points],[x[1] for x in SDCS_ardid.points],"purple",label=r'Antares',linewidth=1.5,linestyle='--')
plt.loglog([x[0] for x in SD[0][4][2].points],[x[1] for x in SD[0][4][2].points],"darkgreen",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3, $V \to b \bar b$',linewidth=1.5)
plt.loglog([x[0] for x in SD[1][4][2].points],[x[1] for x in SD[1][4][2].points],"b",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3, $V \to W^+ W^-$',linewidth=1.5)
plt.loglog([x[0] for x in SD[2][4][1].points[1:]],[x[1] for x in SD[2][4][1].points[1:]],"r",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3, $V \to \tau \bar \tau$',linewidth=1.5)
plt.loglog([x[0] for x in SD[3][4][1].points[1:]],[x[1] for x in SD[3][4][1].points[1:]],"k",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3, $V \to \nu \bar \nu$',linewidth=1.5)
plt.loglog([x[0] for x in SDCS_HAWC_tau.points],[x[1] for x in SDCS_HAWC_tau.points],"r",label=r'HAWC 4 year, $V \to 2 \tau$ ',linewidth=1.5,linestyle='-.')
plt.xlim((1e2,1e5))
plt.ylim((1e-9,1e0))

plt.title(r"Sensitivity")

ax.tick_params(axis='both', length=12, width=1)

ax.tick_params(axis='both', which='minor', length=9)
plt.gca().yaxis.set_minor_locator(mlp.ticker.LogLocator(base=10,subs=range(10)))
legend = plt.legend(loc='upper right')

for label in legend.get_lines():
    label.set_linewidth(1.5)

#plt.text(430.0,5e-8,"IceCube Preliminary",color="r")

plt.ylabel(r"$\sigma_{\rm SD}$ [pb]")
plt.xlabel(r"M$_{\rm DM}$ [GeV]")
plt.savefig("/home/ctoennis/analyses/standard_analysis_framework/SDCS_new.png",format="png")

plt.clf()

fig, ax = plt.subplots()

plt.loglog([x[0] for x in SD[0][3][2].points],[x[1] for x in SD[0][3][2].points],"darkgreen",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.43, $V \to b \bar b$',linewidth=1.5)
plt.loglog([x[0] for x in SD[1][3][2].points],[x[1] for x in SD[1][3][2].points],"b",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.43, $V \to W^+ W^-$',linewidth=1.5)
plt.loglog([x[0] for x in SD[2][3][1].points[1:]],[x[1] for x in SD[2][3][1].points[1:]],"r",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.43, $V \to \tau \bar \tau$',linewidth=1.5)
plt.loglog([x[0] for x in SD[3][3][1].points[1:]],[x[1] for x in SD[3][3][1].points[1:]],"k",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.43, $V \to \nu \bar \nu$',linewidth=1.5)

plt.xlim((1e2,1e5))
plt.ylim((1e-9,1e0))

plt.title(r"Sensitivity")

ax.tick_params(axis='both', length=12, width=1)

ax.tick_params(axis='both', which='minor', length=9)
plt.gca().yaxis.set_minor_locator(mlp.ticker.LogLocator(base=10,subs=range(10)))
legend = plt.legend(loc='upper right')

for label in legend.get_lines():
    label.set_linewidth(1.5)

#plt.text(430.0,5e-8,"IceCube Preliminary",color="r")

plt.ylabel(r"$\sigma_{\rm SD}$ [pb]")
plt.xlabel(r"M$_{\rm DM}$ [GeV]")
plt.savefig("/home/ctoennis/analyses/standard_analysis_framework/SDCS_1s.png",format="png")

plt.clf()


fig, ax = plt.subplots()
plt.loglog([x[0] for x in SDCS_IC_bb.points],[x[1] for x in SDCS_IC_bb.points],"darkgreen",label=r'Jeff Regular DM $b \bar b$',linewidth=1.5,linestyle='-.')
plt.loglog([x[0] for x in SDCS_IC_WW.points],[x[1] for x in SDCS_IC_WW.points],"b",label=r'Jeff Regular DM $W^+ W^- $ ',linewidth=1.5,linestyle='-.')
plt.loglog([x[0] for x in SDCS_IC_tt.points],[x[1] for x in SDCS_IC_tt.points],"r",label=r'Jeff Regular DM $\tau \bar \tau $ ',linewidth=1.5,linestyle='-.')
plt.loglog([x[0] for x in SD[0][0][2].points],[x[1] for x in SD[0][0][2].points],"darkgreen",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.0043, $V \to b \bar b$',linewidth=1.5)
plt.loglog([x[0] for x in SD[1][0][2].points],[x[1] for x in SD[1][0][2].points],"b",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.0043, $V \to W^+ W^-$',linewidth=1.5)
plt.loglog([x[0] for x in SD[2][0][2].points[1:]],[x[1] for x in SD[2][0][2].points[1:]],"r",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.0043, $V \to \tau \bar \tau$',linewidth=1.5)
plt.xlim((2e2,1e4))
plt.ylim((1e-6,1e2))

plt.title(r"Sensitivity")

ax.tick_params(axis='both', length=12, width=1)

ax.tick_params(axis='both', which='minor', length=9)

legend = plt.legend(loc='upper right')
plt.gca().yaxis.set_minor_locator(mlp.ticker.LogLocator(base=10,subs=range(10)))
for label in legend.get_lines():
    label.set_linewidth(1.5)

#plt.text(430.0,5e-8,"IceCube Preliminary",color="r")

plt.ylabel(r"$\sigma_{\rm SD}$ [pb]")
plt.xlabel(r"M$_{\rm DM}$ [GeV]")
plt.savefig("/home/ctoennis/analyses/standard_analysis_framework/SDCS_new_short.png",format="png")



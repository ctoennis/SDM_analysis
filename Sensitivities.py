#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00

import ingredients as ING
import matplotlib as mlp
import numpy as np
mlp.use('Agg')
import matplotlib.pyplot as plt

ACC     = []
NS      = []
FLUX    = []
SI      = []
SD      = []
SDQ     = []

N_anu   = {
   
    1     : [0.6454527052109089, 0.4771622587189924, 0.3364361461264027, 0.07341974563399863],
    100   : [0.7196091449693658, 0.586229625649324, 0.43900865987863497, 0.10950276425027326],
    1000  : [0.8074306202923366, 0.7850870866423472, 0.7276274092481274, 0.48843091103503694],
    10000 : [0.9029815603283379, 0.9252520371364517, 0.929108033616832, 0.91286407856481]

}

N_nu    = {

    1     : [0.7321528225532506, 0.5681874254664947, 0.4112306182996759, 0.08227194311962774],
    100   : [0.7870679281591981, 0.6715309452153471, 0.5183512827299065, 0.12346864950749073],
    1000  : [0.8379634432110211, 0.8022419454366084, 0.7522059132783657, 0.5143513808528943],
    10000 : [0.9029815603283379, 0.9252520371364517, 0.929108033616832, 0.91286407856481]
   
}

fluxconv   = ING.Graph.FromFile("/home/ctoennis/analyses/SDM/rawfactor.txt")
fluxconv2   = ING.Graph.FromFile("/home/ctoennis/analyses/SDM/rawfactor_2.txt")

ardid      = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/mardid.txt")

antares    = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/ANTARES.txt")

SDCS_ardid = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/mardid_SDCS.txt")

SICS_ardid = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/mardid_SICS.txt")

SDCS_HAWC  = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/HAWC.txt")

SDCS_HAWC_tau  = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/HAWC_tau.txt")

#print ardid.points[1:]

AU = 149597870.7

conversion = 1e-36

plt.rc('font', size=24)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize
plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rc('text', usetex=True)

plt.xticks(fontsize=14)

for LIFE in [1, 100, 1000, 10000]:

    ACC_i = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/"+str(LIFE)+"_track.txt")
    ACC.append(ACC_i)
    NS_i = ING.Graph.FromFile("/home/ctoennis/analyses/standard_analysis_framework/"+str(LIFE)+"_track_flux.txt")
    NS.append(NS_i)
    
    FLUX_p = []
    SD_p   = []
    SI_p   = []
    SD_q   = []

    for a,b,c,d in zip(ACC_i.points,NS_i.points,N_nu[LIFE],N_anu[LIFE]):

        FLUX_p.append([a[0],b[1]*365.0*86400*1000000.0/a[1]]) #*2.0*np.pi)])
	
	print(FLUX_p[-1])

#	SD_p.append([a[0],b[1]*365.0*86400*1000000.0*fluxconv.Eval(a[0])/(a[1]*(c+d))])
        SD_p.append([a[0],b[1]*1000000.0*8.0*np.pi*np.power(AU,2)/(fluxconv2.Eval(a[0])*a[1]*(c+d))])
 
        SD_q.append([a[0],b[1]*1000000.0*np.power(a[0]/50.0,2)*8.0*np.pi*np.power(AU,2)*np.power(10.0,-26)/(9.0*a[1]*(c+d))])
#	SD_q.append([a[0],b[1]*1000000.0*np.power(a[0]/1000.0,2)*8.0*np.pi*np.power(AU,2)*np.power(10.0,-24)/(2.77*a[1]*(c+d))])
	SI_p.append([a[0],b[1]*1000000.0*np.power(a[0]/1000.0,2)*8.0*np.pi*np.power(AU,2)*np.power(10.0,-24)/(4270.0*a[1]*(c+d))])
#	print([a[0],b[1]*1000000.0*np.power(a[0]/1000.0,2)*4.0*np.pi*np.power(AU,2)*np.power(10.0,-24)/(4270.0*2.0*a[1]*(c+d))])
#	print([a[0],b[1]*365.0*86400*1000000.0/a[1],a[1]],b[1])
#	print(b[1])

    FLUX_i = ING.Graph(FLUX_p)
    SD_i   = ING.Graph(SD_p)
    SI_i   = ING.Graph(SI_p)
    SD_I   = ING.Graph(SD_q)

    FLUX.append(FLUX_i)
    SD.append(SD_i)
    SI.append(SI_i)
    SDQ.append(SD_I)

fig, ax = plt.subplots()



plt.loglog([x[0] for x in FLUX[0].points],[x[1] for x in FLUX[0].points],"purple",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.00043',linewidth=1.5)
plt.loglog([x[0] for x in FLUX[1].points],[x[1] for x in FLUX[1].points],"k",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.043',linewidth=1.5)
plt.loglog([x[0] for x in FLUX[2].points],[x[1] for x in FLUX[2].points],"r",label=r'$\gamma c \tau / R_{\rm sun}$ = 0.43',linewidth=1.5)
plt.loglog([x[0] for x in FLUX[3].points],[x[1] for x in FLUX[3].points],"darkgreen",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3',linewidth=1.5)
#plt.loglog([x[0] for x in ardid.points],[x[1] for x in ardid.points],"r",label=r'External IceCube public data analysis',linewidth=1.5,linestyle='-.')
#plt.loglog([x[0] for x in antares.points],[x[1] for x in antares.points],"blue",label=r'ANTARES',linewidth=1.5,linestyle='-.')

ax.tick_params(axis='both', length=12, width=1)

ax.tick_params(axis='both', which='minor', length=9)
plt.title(r"Flux sensitivity")
plt.ylim((1e5,1e10))
plt.xlim((2e2,1e4))
legend = plt.legend(loc='upper right')#,bbox_to_anchor=(1.1, 1.0))

for label in legend.get_lines():
    label.set_linewidth(1.5)

plt.text(230.0,5e5,"IceCube work in progress",color="r")

plt.ylabel(r"Flux [Yr$^{-1}$km$^{-2}]$")
plt.xlabel(r"M$_{\rm DM}$ [GeV]")
plt.savefig("/home/ctoennis/analyses/standard_analysis_framework/flux_ps.pdf",format="pdf")
plt.savefig("/home/ctoennis/analyses/standard_analysis_framework/flux_ps.png",format="png")


plt.clf()

fig, ax = plt.subplots()


print([x[1]*conversion for x in SD[3].points])
#plt.loglog([x[0] for x in SD[0].points],[x[1] for x in SD[0].points],"purple",label=r'$c \gamma \tau / R_{\rm sun}$ = 0.00043',linewidth=1.5)
#plt.loglog([x[0] for x in SD[1].points],[x[1] for x in SD[1].points],"k",label=r'$c \gamma \tau / R_{\rm sun}$ = 0.043',linewidth=1.5)
#plt.loglog([x[0] for x in SD[2].points],[x[1] for x in SD[2].points],"r",label=r'$c \gamma \tau / R_{\rm sun}$ = 0.43',linewidth=1.5)
#plt.loglog([x[0] for x in SDQ[3].points],[x[1]*conversion for x in SDQ[3].points],"g",label=r'old conversion',linewidth=1.5)
plt.loglog([x[0] for x in SDCS_ardid.points],[x[1]*conversion for x in SDCS_ardid.points],"b",label=r'ANTARES 2017, $V \to \nu \bar \nu$',linewidth=1.5,linestyle='-.')
plt.loglog([x[0] for x in SDCS_HAWC_tau.points],[x[1]*conversion for x in SDCS_HAWC_tau.points],"purple",label=r'HAWC 4 year, $V \to 2 \tau$ ',linewidth=1.5,linestyle='-.')
plt.loglog([x[0] for x in SD[3].points],[x[1]*conversion for x in SD[3].points],"k",label=r'This analysis, $\gamma c \tau / R_{\rm sun}$ = 4.3, $V \to \nu \bar \nu$',linewidth=1.5)
plt.loglog([x[0] for x in SDCS_HAWC.points],[x[1]*conversion for x in SDCS_HAWC.points],"darkgreen",label=r'HAWC 4 year, $V \to 2 \gamma$ ',linewidth=1.5,linestyle='-.')
#plt.loglog([x[0] for x in SDCS_HAWC_tau.points],[x[1] for x in SDCS_HAWC_tau.points],"r",label=r'HAWC 4 year, $V \to 2 \tau$ ',linewidth=1.5,linestyle='-.')
plt.xlim((2e2,1e4))
plt.ylim((1e-9*conversion,1e-3*conversion))

plt.title(r"Sensitivity")

ax.tick_params(axis='both', length=12, width=1)

ax.tick_params(axis='both', which='minor', length=9)

legend = plt.legend(loc='upper left',bbox_to_anchor=(0.1, 1.025))

for label in legend.get_lines():
    label.set_linewidth(1.5)

plt.text(430.0,5e-35,"IceCube work in progress",color="r")

plt.ylabel(r"$\sigma_{\rm SD}$ [cm$^2$]",labelpad=0)
plt.xlabel(r"M$_{\rm DM}$ [GeV]")
plt.savefig("/home/ctoennis/analyses/standard_analysis_framework/SDCS_ps.pdf",format="pdf")
plt.savefig("/home/ctoennis/analyses/standard_analysis_framework/SDCS_ps.png",format="png")

plt.clf()

fig, ax = plt.subplots()

#plt.loglog([x[0] for x in SI[0].points],[x[1] for x in SI[0].points],"purple",label=r'$c \gamma \tau / R_{\rm sun}$ = 0.00043',linewidth=1.5)
#plt.loglog([x[0] for x in SI[1].points],[x[1] for x in SI[1].points],"k",label=r'$c \gamma \tau / R_{\rm sun}$ = 0.043',linewidth=1.5)
#plt.loglog([x[0] for x in SI[2].points],[x[1] for x in SI[2].points],"r",label=r'$c \gamma \tau / R_{\rm sun}$ = 0.43',linewidth=1.5)
#plt.loglog([x[0] for x in SI[3].points],[x[1] for x in SI[3].points],"K",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3, $V \to \nu \bar \nu$',linewidth=1.5)
plt.loglog([x[0] for x in SICS_ardid.points],[x[1]*conversion for x in SICS_ardid.points],"b",label=r'Ardid et.al. JCAP04(2017)010',linewidth=1.5,linestyle='-.')
plt.loglog([x[0] for x in SI[3].points],[x[1]*conversion for x in SI[3].points],"K",label=r'$\gamma c \tau / R_{\rm sun}$ = 4.3, $V \to \nu \bar \nu$',linewidth=1.5)
plt.ylim((1e-12*conversion,1e-8*conversion))
plt.xlim((2e2,1e4))

plt.title(r"Sensitivity")

ax.tick_params(axis='both', length=12, width=1)

ax.tick_params(axis='both', which='minor', length=9)

legend = plt.legend(loc='upper left',bbox_to_anchor=(0.25, 0.55))

#legend = plt.legend(loc='upper left',bbox_to_anchor=(0.4, 1.05))

for label in legend.get_lines():
    label.set_linewidth(1.5)

plt.text(230.0,2e-12,"IceCube Preliminary",color="r")

plt.ylabel(r"$\sigma_{\rm SI}$ [cm$^2$]",labelpad=0)
plt.xlabel(r"M$_{\rm DM}$ [GeV]")
plt.savefig("/home/ctoennis/analyses/standard_analysis_framework/SICS_ps.pdf",format="pdf")

plt.clf()



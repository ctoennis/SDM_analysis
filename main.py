#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT icerec/V05-01-00

#this is the main code for PE generation and to produce TS distributions

import sys
import os
import argparse

import PEs as PE #in this file the PE generation is handled
import ingredients as ING #histogram implementation
import Likelihoods as LLH #contains different likelihood functions, mainly used by TestStatistics.py
import TestStatistic as TS #contains a TS distribution class that can generate PEs and add the TS to a histogram

#first specify the WIMP mass and lifetime case you want to simulate. Aslo get the number of fake signal events to be generated

Life = 3
Mass = 13
n_s  = 0
Channel = 0
Medmass = 0
loop = 0
feature =  0

parser = argparse.ArgumentParser(description='Process some parameters.')
parser.add_argument('-ns', dest='n_s', default=0, type=int, choices = range(101), metavar='n_s', help='average number of signal events')
parser.add_argument('-m', dest='m', default=0, type=int, choices = range(14), metavar='mass', help='mass number')
parser.add_argument('-l', dest='l', default=0, type=int, choices = range(5), metavar='lifetime', help='lifetime number')
parser.add_argument('-med', dest='med', default=0, type=int, choices = range(5), metavar='medmass', help='mediator mass number')
parser.add_argument('-ch', dest='ch', default=0, type=int, choices = range(5), metavar='channel', help='channel number')
parser.add_argument('-loop', dest='loop', default=0, type=int, choices = range(2), metavar='loop', help='loopmode switch')

args    = parser.parse_args()
Life    = args.l
Mass    = args.m
n_s     = args.n_s
Medmass = args.med
Channel = args.ch
loop    = args.loop

#Get the likelihood components

masslist     = [100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000]
medmasslist  = [1, 10, 100, 1000, 10000]
lifetimelist = [1, 10, 100, 1000, 10000]
channellist  = [5, 8, 11, 13]

PSF_track = []
ENE_track = []

for mass in masslist: #[100,200,350,500,750,1000,2500,5000,7500,10000,25000,50000,75000,100000]:

    a3 = []
    b3 = []

    for channel in channellist: #[5,8,11,13]:

        a2 = []
	b2 = []

        for life in lifetimelist: #[0.00042,0.042,0.42,4.2]:

	    a1 = []
	    b1 = []

	    for medmass in medmasslist: #[10,100,1000,10000]:

	    #	        print "/home/ctoennis/analyses/standard_analysis_framework/northern_ing/WIMP/ENE_m" + str(mass) + "_c" + str(channel) + "_l" + str(life) + "_med" + str(medmass) + ".txt"

	        if os.path.isfile("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/wimp_newspectra/ENERGY_track_m" + str(mass) +  "_med" + str(medmass) + "_l" + str(life) + "_ch" + str(channel) + ".txt"):
	  	    
		    PSF_temp = ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/wimp_newspectra/PSF_WIMPSIM_m" + str(mass) + "_med" + str(medmass) + "_l" + str(life) + "_ch" + str(channel) + ".txt")
		    PSF_temp.Scale(len(PSF_temp.content)*1.0/PSF_temp.integral)
	            a1.append(PSF_temp)
		    b1.append(ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/wimp_newspectra/ENERGY_track_m" + str(mass) +  "_med" + str(medmass) + "_l" + str(life) + "_ch" + str(channel) + ".txt"))
			
		else:

		    a1.append("none")
		    b1.append("none")

	    a2.append(a1)
	    b2.append(b1)

	a3.append(a2)
	b3.append(b2)

    PSF_track.append(a3)
    ENE_track.append(b3)

#print PSF_track[Mass][Channel][Life][Medmass]

BG_track      = ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Background_WIMP.txt")
ENE_BG_track  = ING.H1D.FromFile("/home/ctoennis/analyses/standard_analysis_framework/northern_ing/Background_WIMP_ene.txt")

#normalising

ENE_BG_track.Scale(ENE_BG_track.nbin/ENE_BG_track.integral)

n_bg_track   = sum(BG_track.content)/2.0

BG_track.Scale(BG_track.nbin/BG_track.integral)

#print BG_track.integral
#print ENE_BG_track.integral
#print sum(PSF_track[Mass][Channel][Life][Medmass].content)
#print sum(ENE_track[Mass][Channel][Life][Medmass].content)

#see if you want to do one single job or loop over all mass cases

if loop == 0:

    dist = TS.TSdist(PSF_track[Mass][Channel][Life][Medmass],ENE_track[Mass][Channel][Life][Medmass],BG_track,ENE_BG_track,PSF_track[Mass][Channel][Life][Medmass],ENE_track[Mass][Channel][Life][Medmass],BG_track,ENE_BG_track,n_bg_track,0.0) #This is a TS distribution class

    if PSF_track[Mass][Channel][Life][Medmass] == "none" or ENE_track[Mass][Channel][Life][Medmass] == "none":

	exit()

    while dist.TS.integral < 4000:

        dist.AddPE(n_s) #Generate PE and add entry to TS distribution

    dist.TS.Scale(1.0/dist.TS.integral) #normalize
    dist.TS.Write("/home/ctoennis/analyses/standard_analysis_framework/new_TS/TS_m" + str(masslist[Mass]) + "-med" + str(medmasslist[Medmass]) + "-gl" + str(lifetimelist[Life]) + "-ch" + str(channellist[Channel]) + "-ns"+ str(n_s)+".txt")
    dist.NS.Write("/home/ctoennis/analyses/standard_analysis_framework/new_TS/NS_m" + str(masslist[Mass]) + "-med" + str(medmasslist[Medmass]) + "-gl" + str(lifetimelist[Life]) + "-ch" + str(channellist[Channel]) + "-ns"+ str(n_s)+".txt")


else:

    for Mass in range(len(masslist)):

        if PSF_track[Mass][Channel][Life][Medmass] == "none" or ENE_track[Mass][Channel][Life][Medmass] == "none":

	    continue
 
        dist = TS.TSdist(PSF_track[Mass][Channel][Life][Medmass],ENE_track[Mass][Channel][Life][Medmass],BG_track,ENE_BG_track,PSF_track[Mass][Channel][Life][Medmass],ENE_track[Mass][Channel][Life][Medmass],BG_track,ENE_BG_track,n_bg_track,0.0)

        while dist.TS.integral < 10000:

            dist.AddPE(n_s)

	dist.TS.Scale(1.0/dist.TS.integral)
	dist.TS.Write("/home/ctoennis/analyses/standard_analysis_framework/new_TS/TS_m" + str(masslist[Mass]) + "-med" + str(medmasslist[Medmass]) + "-gl" + str(lifetimelist[Life]) + "-ch" + str(channellist[Channel]) + "-ns"+ str(n_s)+".txt")
	dist.NS.Write("/home/ctoennis/analyses/standard_analysis_framework/new_TS/NS_m" + str(masslist[Mass]) + "-med" + str(medmasslist[Medmass]) + "-gl" + str(lifetimelist[Life]) + "-ch" + str(channellist[Channel]) + "-ns"+ str(n_s)+".txt")

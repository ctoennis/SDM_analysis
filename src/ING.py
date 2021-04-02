from itertools import izip
import numpy as np

class H1D:

    'This is a histogram class for one  dimensional histograms. H1D.content gives you the bin contents. low and high are the ayis borders. nbin gives the number of bins. A histograms is initialised with the low and the high value and the bin contents in an array or list in that order.'


    def __init__(self,low,high,give=[0]):

        self.content    = give
        self.low        = float(low)
        self.high       = float(high)
        self.nbin       = len(give)
        self.integral   = float(sum(give))
        self.ran        = float(high) - float(low)
        self.acumen     = np.cumsum(give)
        self.edges      = map(lambda x: low + x*(high-low)/self.nbin,range(0,self.nbin+1))
        self.center     = map(lambda x: low + (high-low)/(2*self.nbin)+x*(high-low)/(self.nbin),range(0,self.nbin))

    @classmethod
    def Empty(cls,low,high,nbins):

        'make empty histograms with nbins bins'
        
        return cls(low,high,[0]*nbins)

    @classmethod
    def FromFile(cls,filename):

        'load a histogram from a file'

        infile = open(filename)
        header = 0
        content = []
        low=0.0
        high=1.0

        for line in infile:

            if header == 0:

                if not line:

                    continue

                elif line == "Ran an die Bouletten! \n":

                    header = 1

                elif line.strip()[0] == "#":

                    continue

                else:

                    split = [x for x in line.split("=")]

                    if len(split) == 2:

                        if split[0].strip() == "low":

                            low = float(split[1])

                        if split[0].strip() == "high":

                            high = float(split[1])


                    else:

                        continue                    

            else:

                content.append(float(line))

        return cls(low,high,content)

    def GetBin(self,value):

        'GetBin(value) returns the number of the bin at value.'

        corner = 0

        for step in range(1,int(np.log2(self.nbin))+1):

            if value >= self.edges[corner+int(self.nbin/np.power(2,step))]:              #look at which half, then quarter, then 1/8......

                corner += int(self.nbin/np.power(2,step))

        for step in range(0,int(np.log2(self.nbin))):                

            if corner < self.nbin - 1:

                if value >= self.edges[corner+1]:                                               #For uncomfortable bin numbers

                    corner+=1

            else:

                break

        if corner >= self.nbin:

            corner = self.nbin - 1

        return corner


    def GetRandom(self):

        'This function returns a random value according to the histogram contents statistics'

        ranme = np.random.ranf()    

        corner = 0

        for step in range(1,int(np.log2(self.nbin))+1):                                 #minimum number of steps

            if ranme >= self.acumen[corner-1+int(self.nbin/np.power(2,step))]/self.integral:    #look at which half, then quarter, then 1/8......

                corner += int(self.nbin/np.power(2,step))

        for step in range(0,int(np.log2(self.nbin))):                

            if ranme >= self.acumen[corner]/self.integral:

                corner += 1

        return self.edges[corner] + np.random.ranf()*self.ran/self.nbin

    def Fill(self,value,weight=1.):

        'Fill(value,peso) fills the histogram at value with the weight peso. The default weight is 1.'

        din                 = self.GetBin(value)
        self.content[din]  += weight
        self.integral      += weight
        self.acumen         = np.cumsum(self.content)

    def SetContent(self,value,index):

        'SetContent(Value,Index) sets the content of bin with Index to Value'

        self.content[index] = value
        self.integral       = sum(self.content)
        self.acumen         = np.cumsum(self.content)

    def GetContent(self,value):

        'GetContent(Value) gets the content of the bin at value.'

        if value < self.edges[0] or value > self.edges[-1]:

            return(0.0)

        din = self.GetBin(value)

        return float(self.content[din])

    def Integral(self,first,last,Mode="bin"):

        'Integral(first,last,mode) gives the integral fromfirst to last. if mode is bin first and last are taken as bin number, value takes them as values'

        integral = 0.0

        if Mode == "bin":

            integral = sum(self.content[first:last+1])

        elif Mode == "value":

            integral = sum(self.content[int(self.GetBin(first)):int(self.GetBin(last)+1)])

        return integral

    def AddH1(self,hist,weight=1.0):

        'add a histogram with the same binning to the histogram on which this method is called using a weight.'

        newcont         = [x+y*weight for x,y in izip(self.content, hist.content)]
#        print [x + y for x,y in izip(self.content, hist.content)]

        self.content    = newcont
        self.integral  += sum(newcont)
        self.acumen     = np.cumsum(self.content)

    def Scale(self,weight):

        newcont         = map(lambda x: x*weight,self.content)
        self.content    = newcont
        self.integral   = float(sum(self.content))
        self.acumen     = np.cumsum(self.content)

    def GetQuantile(self,quant):

        'get a quantile of the histogram'

        corner = self.nbin//2

        while self.acumen[corner]/self.integral <= quant and self.acumen[corner+1]/self.integral >= quant:

            if self.acumen[corner]/self.integral <= quant:

                corner += 1

            else:

                corner -= 1

        return [self.edges[corner+1],corner]

    def GetFCInterval(self,CL):

        'get feltman cousins interval for confidence level CL'

        sume = 0.0
        low  = self.nbin
        high = 0

        CL *= sum(self.content)

        indexes = sorted(zip(range(self.nbin),self.content), key = lambda x: -x[1])

        for i in indexes:

            sume += i[1]
            
            if low  > i[0]:

                low = i[0]

            if high < i[0]:

                high = i[0]

            if sume > CL:

                break

        return (low,high)

    def Multiply(self,factor):

        'multiply content by factor'

        newcont = np.multiply(self.content,factor)
        self.content    = newcont
        self.integral   = float(sum(self.content))
        self.acumen     = np.cumsum(self.content)

    def Write(self,filename):

        'write histogram to file'

        infile = open(filename,'w')

        infile.write("low   =   " + str(self.low) + "\n")
        infile.write("high  =   " + str(self.high) + "\n")
        infile.write("Ran an die Bouletten! \n")
        
        for line in self.content:
            
            infile.write(str(line) + "\n")


        
class H2D:

    """This is a histogram class for two dimensional histograms. H2D.content gives you the bin contents. lowx/y and highx/y are the ayis borders. nbinx/y gives the number of bins.
    A histograms is initialised with the low and the high value and the bin contents in an array or list in that order."""

    def __init__(self,lowx,highx,lowy,highy,give=[[0]]):

        self.content    = give
        self.lowx       = float(lowx)
        self.highx      = float(highx)
        self.lowy       = float(lowy)
        self.highy      = float(highy)
        self.nbinx      = len(give)
        self.nbiny      = len(give[0])
        self.transpose  = map(list,zip(*self.content))
        self.integralx  = map(sum,give)
        self.integraly  = map(sum,self.transpose)
        self.integral   = sum(self.integralx)
        self.ranx       = float(highx) - float(lowx)
        self.rany       = float(highy) - float(lowy)
        self.acumenx    = np.cumsum(self.integralx)
        self.acumeny    = map(np.cumsum,self.content)
        self.edgesx     = map(lambda x: lowx + x*(highx-lowx)/self.nbinx,range(0,self.nbinx+1))
        self.edgesy     = map(lambda x: lowy + x*(highy-lowy)/self.nbiny,range(0,self.nbiny+1))

    @classmethod
    def Empty(cls,lowx,highx,lowy,highy,nbinsx,nbinsy):

        'make empty histograms with nbins bins'
        content     = [[0]*nbinsy]*nbinsx
        lowx        = float(lowx)
        highx       = float(highx)
        lowy        = float(lowy)
        highy       = float(highy)

        return cls(lowx,highx,lowy,highy,content)

    @classmethod
    def FromFile(cls,filename):

        'load histogram from file'

        infile = open(filename)
        header = 0
        content = []
        lowx = 0.0
        highx = 1.0
        lowy = 0.0
        highy = 1.0

        for line in infile:

            if header == 0:

                if not line:

                    continue

                elif line == "Ran an die Bouletten!\n":

                    header = 1

                elif line[0] == "#":

                    continue

                else:

                    split = [x.strip() for x in line.split("=")]

                    if len(split) == 2:
                      
                        if split[0] == "lowx":

                            lowx = float(split[1])

                        if split[0] == "highx":

                            highx = float(split[1])

                        if split[0] == "lowy":

                            lowy = float(split[1])

                        if split[0] == "highy":

                            highy = float(split[1])
                    else:

                        continue                    

            else:

                content.append([float(x) for x in line])

        return cls(lowx,highx,lowy,highy,content)

    def GetBin(self,valuex,valuey):

        'GetBin(valuex,valuey) returns the number of the bin at valuex/valuey.'

        cornerx = 0
        cornery = 0

        for step in range(1,int(np.log2(self.nbinx))+1):

            if valuex >= self.edgesx[cornerx+int(self.nbinx/np.power(2,step))]:              #look at which half, then quarter, then 1/8......

                cornerx += int(self.nbinx/np.power(2,step))

        for step in range(0,int(np.log2(self.nbinx))):                

            if valuex >= self.edgesx[cornerx+1]:                                               #For uncomfortable bin numbers

                cornerx+=1

        if cornerx >= self.nbinx:

            cornerx = self.nbinx - 1

        for step in range(1,int(np.log2(self.nbiny))+1):

            if valuey >= self.edgesy[cornery+int(self.nbiny/np.power(2,step))]:              #look at which half, then quarter, then 1/8......

                cornery += int(self.nbiny/np.power(2,step))

        for step in range(0,int(np.log2(self.nbiny))):                

            if cornery < self.nbiny:

                if valuey >= self.edgesy[cornery+1]:                                               #For uncomfortable bin numbers

                    cornery+=1

        if cornery >= self.nbiny:

            cornery = self.nbiny - 1

        return [cornerx,cornery]

    def GetRandom(self):

        'This function returns a random value according to the histogram contents statistics'

        ranmex = np.random.ranf()
        ranmey = np.random.ranf()        

        cornerx = 0

        for step in range(1,int(np.log2(self.nbinx))+1):                                 #minimum number of steps

            if ranmex >= self.acumenx[cornerx-1+int(self.nbinx/np.power(2,step))]/self.integral:    #look at which half, then quarter, then 1/8......

                cornerx += int(self.nbinx/np.power(2,step))

        for step in range(0,int(np.log2(self.nbinx))):                

            if ranmex >= self.acumenx[cornerx]/self.integral:

                cornerx += 1

        cornery = 0

        for step in range(1,int(np.log2(self.nbiny))+1):                                 #minimum number of steps

            if ranmey >= self.acumeny[cornerx][cornery-1+int(self.nbiny/np.power(2,step))]/self.integraly[cornerx]:    #look at which half, then quarter, then 1/8......

                cornery += int(self.nbiny/np.power(2,step))

        for step in range(0,int(np.log2(self.nbiny))):                

            if ranmey >= self.acumeny[cornerx][cornery]/self.integraly[cornerx]:

                cornery += 1

        return [self.edgesy[cornerx] + np.random.ranf()*self.ranx/self.nbinx,self.edgesy[cornery] + np.random.ranf()*self.rany/self.nbiny]

    def Fill(self,valuex,valuey,weight=1.):

        'Fill(value,peso) fills the histogram at value with the weight peso. The default weight is 1.'

        din=self.GetBin(valuex,valuey)

        self.content[din[0]][din[1]] += weight

        self.integralx  = map(sum,self.content)
        self.integraly  = map(sum,self.transpose)
        self.integral   = sum(self.integralx)

        self.acumenx = np.cumsum(self.integralx)

    def GetContent(self,valuex,valuey):

        'GetContent(Value) gets the content of the bin at value.'

        din = self.GetBin(valuex,valuey)

        return self.content[din[0]][din[1]]

    def AddH2(self,hist,weight=1.0):

        'add a histogram with the same binning to the histogram on which this method is called using a weight.'

        if  self.nbinsx != hist.nbinsx:

            raise IndexError

        elif self.nbinsy != hist.nbinsy:

            raise IndexError

        newcont=[]

        for line1,line2 in zip(hist.content,self.content):
            newcont.append([sum(x)*weight for x in zip(line1, line2)])
        self.content = newcont

        #Update peripherical parameters
        
        self.transpose  = map(list,zip(*self.content))
        self.integralx  = map(sum,self.content)
        self.integraly  = map(sum,self.transpose)
        self.integral   = sum(self.integralx)
        self.acumenx    = np.cumsum(self.integralx)
        self.acumeny    = map(np.cumsum,self.content)

    def Scale(self,weight):

        'Scale the histogram by factor'

        newcont         = [map(lambda x: x*weight,p) for p in self.content]
        self.content    = newcont
        self.transpose  = map(list,zip(*self.content))
        self.integralx  = map(sum,self.content)
        self.integraly  = map(sum,self.transpose)
        self.integral   = sum(self.integralx)
        self.acumenx    = np.cumsum(self.integralx)
        self.acumeny    = map(np.cumsum,self.content)

    def Write(self,filename):

        'write the histogram to a file'

        infile = open(filename,"w")

        infile.write("lowx    =" + str(self.lowx) + "\n")
        infile.write("highx   =" + str(self.highx) + "\n")
        infile.write("lowy    =" + str(self.lowy) + "\n")
        infile.write("highy   =" + str(self.highy) + "\n")

        infile.write("Ran an die Bouletten \n")

        for line in self.content:
                
            infile.write(reduce(lambda x,y: str(x) + "   " + str(y), line) + "\n")

class H3D:

    def __init__(self,low,high,give=[H2D(0.,1.,0.,1.,[[0]])]):

        'This is a 3-D histogram method derived from 2-D histogrems. Supply the z low end and high end value and a list of H2D to initialise. All methods are derived from H2D.'

        self.content        = give
        self.nbinz          = len(give)
        self.nbinx          = give[0].nbinx
        self.nbiny          = give[0].nbiny
        self.ranz           = float(high) - float(low)
        self.integralz      = [float(x.integral) for x in give]
        self.integral       = sum(self.integralz)
        self.acumen         = np.cumsum(self.integralz)
        self.edgesz         = map(lambda x: low + x*(high-low)/self.nbinz,range(0,self.nbinz+1))
        self.edgesx         = give[0].edgesx
        self.edgesy         = give[0].edgesy

    @classmethod
    def Empty(cls,lowx,highx,lowy,highy,lowz,highz,nbinsx,nbinsy,nbinsz):

        'make empty histograms with nbins bins'
        return cls(lowz,highz,[H2D.Empty(lowx,highx,lowy,highy,nbinsx,nbinsy)]*nbinsz)

    @classmethod
    def FromFile(cls,filename):

        'load histogram from file'

        infile = open(filename)
        header = 0

        for line in infile:

            if header == 0:

                if not line:

                    continue

                elif line.strip()[0] == "#":

                    continue

                else:

                    split = [x.strip() for x in line.split("=")]

                    if len(split) == 2:

                        if split[0] == "lowx":

                            lowx = split[1]
                            
                        if split[0] == "lowy":

                            lowy = split[1]

                        if split[0] == "lowz":

                            lowz = split[1]

                        if split[0] == "highx":

                            highx = split[1]
                            
                        if split[0] == "highy":

                            highy = split[1]

                        if split[0] == "highz":

                            highz = split[1]


    def GetRandom(self):

        'This function returns a random value according to the histogram contents statistics'

        ranme = np.random.ranf()    

        corner = 0

        for step in range(1,int(np.log2(self.nbinz))+1):                                 #minimum number of steps

            if ranme >= self.acumen[corner+int(self.nbinz/np.power(2,step))]/self.integral:    #look at which half, then quarter, then 1/8......

                corner += int(self.nbinz/np.power(2,step))

        for step in range(0,int(np.log2(self.nbinz))):                

            if ranme >= self.acumen[corner+1]/self.integral:

                corner += 1

        return self.content[corner].GetRandom() + [self.edgesz[corner] + np.random.ranf()*self.ranz/self.nbin]

    def GetBin(self,valuex,valuey,valuez):

        'GetBin(value) returns the number of the bin at value.'

        corner = 0

        for step in range(1,int(np.log2(self.nbinz))+1):

            if value >= self.edgesz[corner+int(self.nbin/np.power(2,step))]:              #look at which half, then quarter, then 1/8......

                corner += int(self.nbinz/np.power(2,step))

        for step in range(0,int(np.log2(self.nbinz))):                

            if value >= self.edgesz[corner+1]:                                               #For uncomfortable bin numbers

                corner+=1

        return self.content[0].GetBin(valuex,valuey) + [corner]

    def GetBinZ(self,valuez):

        'GetBin(value) returns the number of the z-bin at value.'

        corner = 0

        for step in range(1,int(np.log2(self.nbinz))+1):

            if valuez >= self.edgesz[corner+int(self.nbinz/np.power(2,step))]:              #look at which half, then quarter, then 1/8......

                corner += int(self.nbinz/np.power(2,step))

        for step in range(0,int(np.log2(self.nbinz))):                

            if valuez >= self.edgesz[corner+1]:                                               #For uncomfortable bin numbers

                corner+=1

        if corner > self.nbinz:

            corner = nbinz

        return corner

    def Fill(self,valuex,valuey,valuez,weight=1.):

        'Fill(value,peso) fills the histogram at value with the weight peso. The default weight is 1.'

        din=self.GetBinZ(valuez)

        self.content[din].Fill(valuex,valuey,weight)

        self.integralz  = [float(x.integral) for x in self.content]
        self.integral   = sum(self.integralz)
        self.acumen     = np.cumsum(self.integralz)

    def GetContent(self,valuex,valuey,valuez):

        'GetContent(Value) gets the content of the bin at value.'

        din = self.GetBinz(valuez)

        return self.content[din[0]].GetContent(valuex,valuey)

    def AddH3(self,hist,weight=1.0):

        'add a histogram with the same binning to the histogram on which this method is called using a weight.'

        newcont=[]

        for line1,line2 in izip(hist.content,self.content):
            line3 = line2
            line3.AddH2(line1,weight)
            newcont.append(line3)

        self.content = newcont

        #Update peripherical parameters
        
        self.integralz  = [float(x.integral) for x in self.content]
        self.integral   = sum(self.integral)
        self.acumen     = np.cumsum(self.integralz)

    def Scale(self,weight):

        newcont=[x.Scale(weight) for x in self.content]
        self.content=newcont
        
        self.integralz  = [float(x.integral) for x in self.content]
        self.integral   = sum(self.integral)
        self.acumen     = np.cumsum(self.integralz)

class Graph:

    'This is a graph class that allows evaluation at a point and integration as well as scaling'

    def __init__(self,points):

        self.points=sorted(points,key=lambda x: x[0])
        self.N=len(self.points)

    @classmethod
    def FromFile(cls,filename):

        'load graph from file'

        infile      = open(filename)
        points = []
        
        for line in infile:

            split2 = [float(x.strip()) for x in line.split()]

            if len(split2) == 2:
            
                points.append(split2)

        return cls(points)

    def Eval(self,value):

        'get graph value at specific point'

        corner = 0

	if value < self.points[0][0]:

	    return 0.0

	if value > self.points[self.N-1][0]:

	    return 0.0

        if value == self.points[self.N-1][0]:

            return self.points[self.N-1][1]

        for step in range(1,int(np.log2(self.N))+1):                                 #minimum number of steps

            if value >= self.points[corner+int(self.N/np.power(2,step))][0]:    #look at which half, then quarter, then 1/8......

                corner += int(self.N/np.power(2,step))

        for step in range(0,int(np.log2(self.N))):                

            if value >= self.points[corner+1][0]:

                corner += 1

        rest    = value - self.points[corner][0]
        out     = self.points[corner][1] + rest * (self.points[corner+1][1] - self.points[corner][1])/(self.points[corner+1][0] - self.points[corner][0])

        return out

    def Integral(self,low,high):

        'get integral of graph in limits'

        cornerlow = 0

        if low <= self.points[0][0]:
            
            cornerlow = 0

        else:

            for step in range(1,int(np.log2(self.N))+1):                                 #minimum number of steps

                if low >= self.points[cornerlow+int(self.N/np.power(2,step))][0]:    #look at which half, then quarter, then 1/8......

                    cornerlow += int(self.N/np.power(2,step))

            for step in range(0,int(np.log2(self.N))):

                if low >= self.points[cornerlow+1][0]:

                    cornerlow += 1
    
        cornerhigh = 0

        if high >= self.points[-1][0]:

            cornerhight = self.N - 1

        else:

            for step in range(1,int(np.log2(self.N))+1):                                 #minimum number of steps

                if high >= self.points[cornerhigh+int(self.N/np.power(2,step))][0]:    #look at which half, then quarter, then 1/8......

                    cornerhigh += int(self.N/np.power(2,step))

            for step in range(0,int(np.log2(self.N))):

                if high >= self.points[cornerhigh+1][0]:

                    cornerhigh += 1

        integral = 0.0

        for point,point2 in zip(self.points[cornerlow+1:cornerhigh],self.points[cornerlow+2:cornerhigh+1]):

            integral += ((point[1]+point2[1])/2.0)*(point2[0]-point[0])

        integral += ((self.points[cornerlow][1]+self.points[cornerlow+1][1])/2.0)*(self.points[cornerlow+1][0]-low)

        integral += ((self.points[cornerhigh][1]+self.points[cornerhigh+1][1])/2.0)*(high-self.points[cornerhigh][0])

        return integral
            

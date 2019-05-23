

def ReadFrequencyMatrix(Infile):
    import csv
    f = open(Infile)
    reader = csv.reader(f)
    keys = []
    values = []
    next(reader)
    for line in reader:
        keys.append(int(line[0]))
        values.append(int(line[1])) 
    f.close()
    FrequencyMatrix =  dict(zip(keys, values))
    return FrequencyMatrix

def AccumulationCurve(FrequencyMatrix, TotalSamples=1473):
    #Generates a species accumulation curve using hypergeometric distribution
    from scipy.stats import hypergeom 
    import numpy
    #hypergeom.pmf(Occurences,TotalSamples,SetSize,Sampled)
    SetSize = FrequencyMatrix.keys()
    N_j = FrequencyMatrix.values()
    #Singletons = numpy.zeros(TotalSamples)
    #Doubletons = numpy.zeros(TotalSamples)
    OccurenceMatrix = numpy.zeros((int(FrequencyMatrix.keys()[-1]),TotalSamples))
    for Sampled in range(TotalSamples):
        print Sampled
        for Ntons in range((int(FrequencyMatrix.keys()[-1]))):
            OccurenceMatrix[Ntons,Sampled] = numpy.dot(N_j,hypergeom.pmf(Ntons+1,TotalSamples,SetSize,Sampled+1))
    return OccurenceMatrix
    
def func(x, a, b, c): #quadratic function, for fitting purposes
    return a*x**2 + b*x + c

def func2(x, a, b, c): #log quadratic function, for fitting purposes
    import numpy as np
    return np.e**(a*np.log(x)*np.log(x)) * x**b * c

def Rsquared(y_real, y_est): #Calculates the coefficient of correlation  
    SST = 0
    SSReg = 0
    y_real_mean = y_real.mean()
    for i in range(len(y_real)):
       SST = SST + (y_real[i] - y_real_mean)**2
       SSReg = SSReg + (y_est[i] - y_real[i])**2
    RR = 1- SSReg/SST
    return RR

def FindSamplesRequired(OccurenceMatrix, popt, percentagearray, fitting='log'):
    #Finds the number of samples required in order for the percentage of new genes per sample to be lower than a certain threshold
    #percentagearray = an array of thresholds, and popt are the fitting parameters
    import numpy as np
    from scipy.optimize import root
    
    NumberOfSamples = np.zeros(len(percentagearray))
    for i in range(len(percentagearray)):
        if fitting == 'log':
            NumberOfSamples[i] = root(lambda x: np.e**func(np.log(x),*popt)/OccurenceMatrix[0,0]*100-percentagearray[i], 1000).x
        else:
            NumberOfSamples[i] = root(lambda x:func2(x,*popt)/OccurenceMatrix[0,0]*100-percentagearray[i], 1000).x
            
    return NumberOfSamples

def PlotDiscoveryCurve(OccurenceMatrix, title, path, filetag, fitting='log'):
    #Least square error criteria can either be in logarithmic or in linear space
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    
    #Calculates the number of new genes per sample based on OccurenceMatrix
    Influx = [OccurenceMatrix[0,0]]
    Influx.extend(OccurenceMatrix[0,1:] - OccurenceMatrix[0,:-1])
    Influx = np.array(Influx)
    Efflux = [0]
    Efflux.extend(OccurenceMatrix[1:,1:].sum(0) - OccurenceMatrix[1:,:-1].sum(0))
    Efflux = np.array(Efflux)
    ydata = Influx + Efflux
    xdata = np.arange(1.0,float(OccurenceMatrix.shape[1]+1))
    xdata2 = np.logspace(0,8,1000)
    
    fig = plt.figure(dpi=300, figsize=(8, 6))
    plt.xlabel('Samples')
    plt.ylabel('Percentage (%)')
    plt.title('Discovery curve in percentage of new genes, '+ title)
    #Plot Gene discovery curve
    plt.loglog(xdata, ydata/OccurenceMatrix[0,0]*100, color = 'green', label = 'Percentage of new genes discovered in sample')
    
    if fitting=='log': #fitting to the log(data) to a quadratic function (right biased)
        popt, pcov = curve_fit(func, np.log(xdata), np.log(ydata))
        plt.loglog(xdata2, np.e**func(np.log(xdata2), *popt)/OccurenceMatrix[0,0]*100, color = 'red', label = 'Fitting f(x) ='+ format(-popt[1], '.0f')+'*x^' + format(-popt[0], '.3f'))
        RR= Rsquared(np.array(ydata),np.e**func(np.log(xdata), *popt))
    else: #fitting data to a Log(quadratic) function (left biased)    
        popt2, pcov2 = curve_fit(func2, xdata, ydata)
        plt.loglog(xdata2, func2(xdata2, *popt2)/OccurenceMatrix[0,0]*100, label = 'Fitting f(x) ='+ format(-popt[1], '.0f')+'*x^' + format(-popt[0], '.3f'))
        RR= Rsquared(np.array(ydata),np.array(func2(xdata, *popt2)))
   
    plt.show()    
    fig.savefig(path+filetag+'GeneDiscovery.png')
    
    return popt, RR
 

def PlotExtrapolator2(OccurenceMatrix, power10, title, path, filetag):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    import scipy.integrate as integrate

    Influx = [OccurenceMatrix[0,0]]
    Influx.extend(OccurenceMatrix[0,1:] - OccurenceMatrix[0,:-1])
    Influx = np.array(Influx)
    Efflux = [0]
    Efflux.extend(OccurenceMatrix[1:,1:].sum(0) - OccurenceMatrix[1:,:-1].sum(0))
    Efflux = np.array(Efflux)
    ydata = Influx + Efflux
    xdata = np.arange(1.0,float(OccurenceMatrix.shape[1]+1))
    popt, pcov = curve_fit(func, np.log(xdata), np.log(ydata))
    
    
    xdata2 = np.logspace(0,power10,1000)
    AccumulationCurve = np.zeros(len(xdata2))
   
    Asymptote = integrate.quad(lambda x: np.e**func(np.log(x),*popt), 0, np.inf)[0]
    for i in range(len(xdata2)):
        if i>0:
            AccumulationCurve[i] = max(AccumulationCurve[i-1], integrate.quad(lambda x: np.e**func(np.log(x),*popt), 0, xdata2[i])[0])
        else:
            AccumulationCurve[i] = integrate.quad(lambda x: np.e**func(np.log(x),*popt), 0, xdata2[i])[0]
    poptPercentage, pcov = curve_fit(func, np.log(xdata), np.log(OccurenceMatrix[0,]/OccurenceMatrix.sum(0)))    
    poptSingletons,pcov = curve_fit(func, np.log(xdata), np.log(OccurenceMatrix[0,:]))   
 
    fig = plt.figure(dpi=300, figsize = (16,6))
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    #plt.ylim(0.1,1E9)
    plt.title('Extrapolated rarefaction curve using quadratic log estimator, '+ title)
    plt.semilogx(xdata2, AccumulationCurve, color = (1, 127/256.0,127/256.0), label = 'Extrapolated accumulation curve')
    plt.semilogx(xdata2, np.e**func(np.log(xdata2), *poptPercentage)*AccumulationCurve, color = 'red', label = 'Extrapolated Singletons')
    plt.semilogx(xdata2,  AccumulationCurve -  np.e**func(np.log(xdata2), *poptPercentage)*AccumulationCurve, color = (1, 191/256.0,191/256.0), label = 'Extrapolated Non-singletons')
    plt.semilogx(xdata, OccurenceMatrix[1:,].sum(0), color = (117/256.0, 117/256.0, 117/256.0), label = 'Non-singletons' )
    plt.semilogx(xdata, OccurenceMatrix.sum(0), color = (65/256.0,113/256.0, 156/256.0), label = 'Rarefaction curve')
    plt.semilogx(xdata, OccurenceMatrix[0,], color = 'black', label='Singletons')
    plt.legend()
    plt.show()  
    fig.savefig(path+filetag+'Extrapolator2.png')
    
    
    fig2 = plt.figure(dpi=300, figsize = (12,9))
    plt.title('Fraction of singletons, extrapolated')
    plt.xlabel('Samples')
    plt.ylabel('Fraction')
    #plt.loglog(xdata2, np.e**func(np.log(xdata2), *poptSingletons)/AccumulationCurve, label='extrapolated fraction of singletons')
    plt.loglog(xdata,OccurenceMatrix[0,]/OccurenceMatrix.sum(0), label = 'Fraction of singletons')  
    plt.loglog(xdata2, np.e**func(np.log(xdata2), *poptPercentage), label='Extrapolated fraction of singletons')
    plt.legend()
    fig2.savefig(path+filetag+'SingletonFractions.png')
    
    return Asymptote

def PlotEstimators(OccurenceMatrix, title, path, filetag):
    #Unfortunately, the code in this function requires some editing... I have not automated it. 
    
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    NumSamples = OccurenceMatrix.shape[1]
    fig1= plt.figure(figsize=(12,9), dpi=240)
    ax = plt.subplot(111)
    plt.title('Rarefaction curve + Estimation, '+title)
    plt.ylim((0,1.7E8))
    plt.xlabel('Number of samples')
    plt.ylabel('Number of genes')
    
    ax.plot(range(NumSamples),OccurenceMatrix.sum(0), color='#000000', label='Rarefaction curve')
    ax.plot(range(NumSamples),OccurenceMatrix[0,], color='#404040', label='Rarefaction curve - Singletons')
    ax.plot(range(NumSamples),OccurenceMatrix[1:,].sum(0), color='#808080', label = 'Rarefaction curve - Non-singletons')
             
    Chao2Data = pd.read_csv('./95Chao2.csv')
    Chao2Data = Chao2Data.as_matrix()
    ax.plot(range(NumSamples),Chao2Data[:,1], color='#ff0000', label = 'Chao2 with 95%CI')
    ax.fill_between(range(NumSamples), Chao2Data[:,3], Chao2Data[:,4], alpha=0.2, edgecolor='#ff4040', facecolor='#ff6060', linewidth=1, antialiased=True)

    ChaoLeeACEData = pd.read_csv('./95ChaoLeeACE.csv')
    ChaoLeeACEData = ChaoLeeACEData.as_matrix()
    ax.plot(range(NumSamples),ChaoLeeACEData[:,1], color='#ff4040', label= 'ChaoLee ACE with 95%CI')
    ax.fill_between(range(NumSamples), ChaoLeeACEData[:,3], ChaoLeeACEData[:,4], alpha=0.2, edgecolor='#ff8080', facecolor='#ffa0a0', linewidth=1, antialiased=True)

    ChaoLeeACE2Data = pd.read_csv('./95ChaoLeeACE2.csv')
    ChaoLeeACE2Data = ChaoLeeACE2Data.as_matrix()
    ax.plot(range(NumSamples),ChaoLeeACE2Data[:,1], color='#ff8080', label = 'ChaoLee ACE2 with 95%CI')
    ax.fill_between(range(NumSamples), ChaoLeeACE2Data[:,3], ChaoLeeACE2Data[:,4], alpha=0.2, edgecolor='#ffc0c0', facecolor='#ffe0e0', linewidth=1, antialiased=True)

    Jackknife3Data = pd.read_csv('./95Jackknife3.csv')
    Jackknife3Data = Jackknife3Data.as_matrix()
    ax.plot(range(NumSamples),Jackknife3Data[:,1], color='#0000ff', label = 'Jackknife 3rd order 95%CI')
    ax.fill_between(range(NumSamples), Jackknife3Data[:,3], Jackknife3Data[:,4], alpha=0.2, edgecolor='#4040ff', facecolor='#6060ff', linewidth=1, antialiased=True)                     

    Jackknife6Data = pd.read_csv('./95Jackknife6.csv')
    Jackknife6Data = Jackknife6Data.as_matrix()
    ax.plot(range(NumSamples),Jackknife6Data[:,1], color='#4040ff', label = 'Jackknife 6th order 95%CI')
    ax.fill_between(range(NumSamples), Jackknife6Data[:,3], Jackknife6Data[:,4], alpha=0.2, edgecolor='#8080ff', facecolor='#a0a0ff', linewidth=1, antialiased=True)

    Jackknife9Data = pd.read_csv('./95Jackknife9.csv')
    Jackknife9Data = Jackknife9Data.as_matrix()
    ax.plot(range(NumSamples),Jackknife9Data[:,1], color='#8080ff', label = 'Jackknife 9th order 95%CI')
    ax.fill_between(range(NumSamples), Jackknife9Data[:,3], Jackknife9Data[:,4], alpha=0.2, edgecolor='#c0c0ff', facecolor='#e0e0ff', linewidth=1, antialiased=True)                     
    
    ChaoBunge3Data = pd.read_csv('./95ChaoBunge3.csv')
    ChaoBunge3Data = ChaoBunge3Data.as_matrix()
    ChaoBunge3Data[ChaoBunge3Data>1E10] = np.inf
    ChaoBunge3Data[ChaoBunge3Data<-1E10] = -np.inf
    ax.plot(range(NumSamples),ChaoBunge3Data[:,1], color='#00ff00', label = 'ChaoBunge 3rd order 95%CI')
    ax.fill_between(range(NumSamples), ChaoBunge3Data[:,3], ChaoBunge3Data[:,4], alpha=0.2, edgecolor='#40ff40', facecolor='#60ff60', linewidth=1, antialiased=True)                     

    ChaoBunge6Data = pd.read_csv('./95ChaoBunge6.csv')
    ChaoBunge6Data = ChaoBunge6Data.as_matrix()
    ChaoBunge6Data[ChaoBunge6Data>1E10] = np.inf
    ChaoBunge6Data[ChaoBunge6Data<-1E10] = -np.inf
    ax.plot(range(NumSamples),ChaoBunge6Data[:,1], color='#40ff40', label = 'ChaoBunge 6th order 95%CI')
    ax.fill_between(range(NumSamples), ChaoBunge6Data[:,3], ChaoBunge6Data[:,4], alpha=0.2, edgecolor='#80ff80', facecolor='#a0ffa0', linewidth=1, antialiased=True)                     

    ChaoBunge9Data = pd.read_csv('./95ChaoBunge9.csv')
    ChaoBunge9Data = ChaoBunge9Data.as_matrix()
    ChaoBunge9Data[ChaoBunge9Data>1E10] = np.inf
    ChaoBunge9Data[ChaoBunge9Data<-1E10] = -np.inf
    ax.plot(range(NumSamples),ChaoBunge9Data[:,1], color='#80ff80', label = 'ChaoBunge 9th order 95%CI')
    ax.fill_between(range(NumSamples), ChaoBunge9Data[:,3], ChaoBunge9Data[:,4], alpha=0.2, edgecolor='#c0ffc0', facecolor='#e0ffe0', linewidth=1, antialiased=True)                     

    ax.legend() #loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':10}
    #plt.gca().legend((,'Rarefaction curve - Singletons', 'Rarefaction curve - Non-singletons', 'Chao2 with 95%CI', 'ChaoLee ACE with 95%CI', 'ChaoLee ACE2 with 95%CI', 'Jackknife 2nd order 95%CI', 'Jackknife 5th order 95%CI', 'Jackknife 10th order 95%CI',  'ChaoBunge 2nd order 95%CI', 'ChaoBunge 5th order 95%CI', 'ChaoBunge 9th order 95%CI'))
    #matplotlib.pyplot.errorbar(range(1476),Chao2Data2[:,1],Chao2Data2[:,2], capsize = 10)
    plt.show()
    fig1.savefig(path+filetag+'Estimators.png')

    x = {'Chao2':Chao2Data[-1,1], 'ChaoLeeAce':ChaoLeeACEData[-1,1], 'ChaoLeeAce2':ChaoLeeACE2Data[-1,1], 'Jackknife3':Jackknife3Data[-1,1],'Jackknife6':Jackknife6Data[-1,1], 'Jackknife9':Jackknife9Data[-1,1], 'ChaoBunge3':ChaoBunge3Data[-1,1], 'ChaoBunge6':ChaoBunge6Data[-1,1], 'ChaoBunge9':ChaoBunge9Data[-1,1]}    
    return x 

def PlotExtrapolator(OccurenceMatrix, title, path, filetag, estimate=0):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    NumSamples = OccurenceMatrix.shape[1]
    Species = OccurenceMatrix.sum(0)
    xdata = np.arange(1.0,float(len(Species)+1))
    ydata = Species
    popt, pcov = curve_fit(func, xdata, ydata)
    SST = 0
    SSReg = 0
    for i in range(len(Species)):
        SST = SST + (Species[i] - Species.mean())**2
        SSReg = SSReg + (func(xdata[i], *popt) - Species[i])**2
    Rsquared = 1- SSReg/SST

    fig = plt.figure(figsize=(12,9), dpi=240)

    plt.xlabel('Samples')
    plt.ylabel('Genes')    
    ax = plt.subplot(111)
    ax.plot(range(NumSamples),OccurenceMatrix.sum(0), label='Rarefaction curve', color = (65/256.0,113/256.0, 156/256.0))
    ax.plot(range(NumSamples),OccurenceMatrix[0,], label='Rarefaction curve - Singletons', color = (0/256.0, 0/256.0, 0/256.0) )
    ax.plot(range(NumSamples),OccurenceMatrix[1:,].sum(0), label = 'Rarefaction curve - Non-singletons', color = (117.0/256.0, 117.0/256.0, 117.0/256.0) )
    
    if estimate == 0:
        plt.title('Rarefaction curve, '+title)
        ax.legend()
        plt.show()
        fig.savefig(path+filetag+'Rarefaction.png')
    else:
        plt.title('Rarefaction curve + Curve fitting, '+title)
        ax.plot(xdata, func(xdata, *popt), color = '#ff8000', label = 'Fitting f(x) ='+ format(popt[1], '.0f')+'*x^' + format(-popt[0], '.3f'))
        ax.text(0.85, 0.98,r'$R^2$ = '+format(Rsquared, '.3f'),
                horizontalalignment='right',
                verticalalignment='top',
                transform = ax.transAxes)
        ax.legend() #loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':10})
        plt.show()
        fig.savefig(path+filetag+'Extrapolator.png')
    #plt.ylim((0,1E8))

    #This is the estimate for the human population
    return func(7.6E9, *popt)


def main(x, NumSamples = 1473):
    '''
    #x is a list of lists.     
    #x[i][0]: STRING Read from here to create frequency matrices
    #x[i][1]: STRING Write occurence matrix into this file for R to read
    #x[i][2]: STRING Additional plot title tag
    #x[i][3]: STRING Path to write files
    #x[i][4]: STRING File name tag (without extension)
    #x[i][5]: STRING '11111', binary decisions (1 = yes, 0 = 0) to plot...
        x[i][5][0] Estimators (ChaoBunge, etc.), 
        x[i][5][1] Heap's Law extrapolator, 
        x[i][5][2] simple rarefaction curve, 
        x[i][5][3] Gene discovery curve (in percentage)
        x[i][5][4] asymptotic number of genes extrapolation. 
    '''
    #Uses data from a file to generate Occurence Matrix and plots it
    import numpy as np
    import os
    
    power10 = 15 #Extrapolates up to 10**15 samples using quadratic log estimation
    percentagearray = [1,0.1,0.01,0.001,0.0001,0.00001,0.000001] #Finds number of samples required for 1%, 0.1%, 0.01%, 0.001%, and 0.0001% of samples to be new
       
    OccurenceMatrices = {}
    for i in range(len(x)):
        if x[i][5] != '00000':
            print "Reading Frequency Matrix"
            FrequencyMatrix = ReadFrequencyMatrix(x[i][0])
            print "Generating Accumulation Curves"
            OccurenceMatrix = AccumulationCurve(FrequencyMatrix, NumSamples)
            OccurenceMatrices[i] = OccurenceMatrix
            #OccurenceMatrix is an MxN matrix that describes the expected number of times we'd pick up M-1tons after N samples
            #Check
            print "max k-ton is: "+str(len(OccurenceMatrix))
            print "Plotting Disaggregated Accumulation Curves" 
        
            if x[i][5][0] == '1':
                np.savetxt(x[i][1], OccurenceMatrix, delimiter = ',')
                print x[i]
                print "Plotting estimators"
                #This saves the OccurenceMatrix into a file so that R can read it (you need to run the R code separately)
                #The R script RollingSpeciesEstimators creates rolling estimates, and puts results into various files.
                #raw_input("Press Enter to continue... (RUN RollingSpeciesEstimators.R SCRIPT FIRST)")
                os.system('Rscript ./RollingSpeciesEstimators.r')
            #    os.system('C:/Program^ Files/R/R-3.5.0/bin/R.exe CMD BATCH --vanilla --slave RollingSpeciesEstimators.r')
                PlotEstimators(OccurenceMatrix, x[i][2], x[i][3], x[i][4])
            if x[i][5][1] == '1': #Plots rarefaction curve + estimator using Heap's (power) Law 
                PlotExtrapolator(OccurenceMatrix, x[i][2], x[i][3], x[i][4],1)
            if x[i][5][2] == '1': #Plots rarefaction curve
                PlotExtrapolator(OccurenceMatrix, x[i][2], x[i][3], x[i][4],0) 
            if x[i][5][3] == '1': 
                #Plots the gene discovery curve + quadratic estimator
                popt,RR = PlotDiscoveryCurve(OccurenceMatrix, x[i][2], x[i][3], x[i][4])  
                print popt, RR
                #Finds samples required to have 1%, 0.1%, 0.01%, 0.001%, and 0.0001% new genes per sample
                NumberOfSamples = FindSamplesRequired(OccurenceMatrix, popt, percentagearray)
                print NumberOfSamples
            if x[i][5][4] == '1':
                #Estimates total number of genes, and plots up to 10**power10 samples
                Asymptote = PlotExtrapolator2(OccurenceMatrix, power10, x[i][2], x[i][3], x[i][4])
                print Asymptote
            PercentagePlot(OccurenceMatrix,'Amino Acids, 50%','./','Amino50')
    return x


def PercentagePlot(OccurenceMatrix, title, path, filetag):
    #For testing purposes only
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    import numpy as np
    NumSamples = OccurenceMatrix.shape[1]
    xdata = np.arange(1.0,float(OccurenceMatrix.shape[1]+1))
    xdata2 = np.logspace(0,10,1000)
    poptPercentage, pcov = curve_fit(func, np.log(xdata), np.log(OccurenceMatrix[0,]/OccurenceMatrix.sum(0)))
    fig = plt.figure(figsize=(12,9), dpi=240)
    ax = plt.subplot(111)
    plt.title('Extrapolation of fraction of singleton genes, ' + title)
    #plt.ylim((0,1E8))
    plt.xlabel('Samples')
    plt.ylabel('Percentage (%)')    
    plt.ylim((0,110))
    #plt.xlim((0,NumSamples))
    ax.semilogx(xdata,OccurenceMatrix[0,]/OccurenceMatrix.sum(0)*100, color = 'black', label='Rarefaction curve - Singletons (fraction)', )
    ax.semilogx(xdata2, np.e**func(np.log(xdata2), *poptPercentage)*100, color = 'red')
    #ax.legend()
    plt.show()  
    fig.savefig(path+filetag+'RarefactionPercentagePlot.png')

    SST = 0
    SSReg = 0
    meanPercentageSingletons = OccurenceMatrix[0,]/OccurenceMatrix.sum(0)
    meanPercentageSingletons = meanPercentageSingletons.mean()
    for i in range(OccurenceMatrix.shape[1]):
       SST = SST + (OccurenceMatrix[0,i]/OccurenceMatrix.sum(0)[i] - meanPercentageSingletons)**2
       SSReg = SSReg + (np.e**func(np.log(xdata[i]), *poptPercentage) - OccurenceMatrix[0,i]/OccurenceMatrix.sum(0)[i])**2
    print SST, SSReg
    Rsquared = 1- SSReg/SST
    print poptPercentage
    print Rsquared


def PlotNTons(OccurenceMatrix, title,  path, filetag, indices = [0]):
#Deprecated... for testing purposes only
    import matplotlib.pyplot as plt
    
    NumSamples = OccurenceMatrix.shape[1]
    fig = plt.figure(figsize=(12,9), dpi=240)
    ax = plt.subplot(111)
    plt.title('Logplot of rarefaction curve, '+title)
    plt.ylim((1E0,1E8))
    plt.xlabel('Number of samples')
    plt.ylabel('Number of genes')    
    for j in indices:           
        ax.semilogy(range(1,NumSamples),OccurenceMatrix[j,1:], label='%s-tons' %(j+1))
    ax.legend()
    plt.show()  
    fig.savefig(path+filetag+'Rarefaction.png')

if __name__ == '__main__':
    x = [["50perc_cluster_sizes.csv", 'OccurenceMatrix.csv', 'Amino acid gene catalog at 50% identity', './', 'Amino50','11111']]
    main(x,1473)
    


    


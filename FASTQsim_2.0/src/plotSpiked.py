#code to plot background dataset statistics, linear regression curves, and spiked read statistics. 
#stores data to be plotted, while code is executing. 
#When code execution is complete, generates plot of background data, 
#curve fit, and in silico-generated spiked reads. 
import matplotlib.pyplot as plt; 
class plotSpiked(object):
    plotName=""; 
    pos=None; 
    #plot for background data
    readLengths=None
    curveFitReadLengths=None 
    insertProb=None
    curveFitInsertProb=None 
    insertSize=None
    curveFitInsertSize=None 
    insertRepeat=None
    curveFitInsertRepeat=None
    delProb=None 
    curveFitDelProb=None
    delSize=None
    curveFitDelSize=None 
    delRepeat=None 
    curveFitDelRepeat=None
    mutationProb=None
    curveFitMutationProb=None 
    
    #plot for spiked reads. 
    actualPosCount=None
    actualReadLengths=None 
    actualInsertPos=None 
    actualDelPos=None
    actualInsertSize=None 
    actualRepeatSize=None
    actualDelSize=None
    actualMutPos=None 

    def plotReadLengthDist(self):
        fig=plt.figure(); 
        plt.subplots_adjust(hspace=0.8)
        ax=fig.add_subplot(211); 
        ax.bar(self.readLengths.keys(),self.readLengths.values(),color='b',label='Empirical');
        ax.plot(self.curveFitReadLengths.keys(), self.curveFitReadLengths.values(),color='r',label='Regression',linewidth=3); 
        ax.set_title(self.plotName+' Read Length\nBackground Data',fontsize=18)
        ax.set_xlabel('Read Length (in bases)',fontsize=16)
        ax.set_ylabel('Number of Reads',fontsize=16)
        ax.legend(loc=4)

        plt.setp(ax.get_xticklabels(),fontsize=16) 
        plt.setp(ax.get_yticklabels(),fontsize=16)
        ax.grid=True; 
        bx=fig.add_subplot(212); 
        bx.bar(self.actualReadLengths.keys(),self.actualReadLengths.values(),color='b',label='Spiked'); 
        bx.set_title(self.plotName+' Read Length\nSpiked Data',fontsize=18)
        bx.set_xlabel('Read Length (in bases)',fontsize=16) 
        bx.set_ylabel('Number of Reads',fontsize=16)
        plt.setp(bx.get_xticklabels(),fontsize=16)
        plt.setp(bx.get_yticklabels(),fontsize=16) 
        bx.grid=True
        plt.savefig(self.plotName+'SpikingReadLength.png',bbox_inches=0)

    def plotInsertionRate(self):
        fig=plt.figure()
        plt.subplots_adjust(hspace=0.8)
        ax=fig.add_subplot(211); 
        for key in self.insertProb: 
            self.insertProb[key]=self.insertProb[key]/float(self.pos[key]); 
        curveFitUpdated=dict(); 
        for key in self.curveFitInsertProb: 
            if key in self.pos: 
                curveFitUpdated[key]=self.curveFitInsertProb[key]/float(self.pos[key]) 
        ax.bar(self.insertProb.keys(),self.insertProb.values(),color='b',log=True,label='Empirical'); 
        ax.plot(curveFitUpdated.keys(),curveFitUpdated.values(),color='r',label='Regression',linewidth=3); 
        ax.set_ylim([10e-5,1])
        ax.set_title(self.plotName+' Insertion Probability\nBackground Data',fontsize=18) 
        ax.set_xlabel('Base Position Along Read',fontsize=16)
        ax.set_ylabel('Probability of Insertion',fontsize=16)
        plt.setp(ax.get_xticklabels(),fontsize=16) 
        plt.setp(ax.get_yticklabels(),fontsize=16)
        ax.legend(loc=4); 
        plt.setp(ax.get_yticklabels(),fontsize=16)
        ax.grid=True; 

        bx=fig.add_subplot(212);
        for key in self.actualInsertPos: 
            self.actualInsertPos[key]=self.actualInsertPos[key]/float(self.actualPosCount[key]); 
        bx.bar(self.actualInsertPos.keys(),self.actualInsertPos.values(),color='b',log=True,label='Spiked')
        bx.set_ylim([10e-5,1]) 
        bx.set_title(self.plotName+' Insertion Probability\nSpiked Data',fontsize=18) 
        bx.set_xlabel('Base Position Along Read',fontsize=16) 
        bx.set_ylabel('Probability of Insertion',fontsize=16)
        plt.setp(bx.get_xticklabels(),fontsize=16) 
        plt.setp(bx.get_yticklabels(),fontsize=16)
        bx.grid=True 
        plt.savefig(self.plotName+'SpikingInsertionRate.png',bbox_inches=0) 

    def plotInsertionSize(self):
        fig=plt.figure()
        plt.subplots_adjust(hspace=0.8)
        #graph background insertion size statistics 
        barwidth=0.2
        sizeKeys=[i-barwidth for i in self.insertSize.keys()]
        sizeValues=[i for i in self.insertSize.values()]
        ax=fig.add_subplot(211); 
        ax.bar(sizeKeys,sizeValues,width=barwidth,color='b',align='center',log=True,label='Back. Insert.')
        ax.bar(self.insertRepeat.keys(),self.insertRepeat.values(),width=barwidth,color='r',align='center',log=True,label='Back. Repeat Insertion')

        ax.grid=True
        ax.set_title(self.plotName+' Insertion Size\nBackground Data',fontsize=18)
        ax.set_xlabel('Insertion Length (in bases)',fontsize=16)
        ax.set_ylabel('Number of Insertions',fontsize=16)
        plt.setp(ax.get_xticklabels(),fontsize=16) 
        plt.setp(ax.get_yticklabels(),fontsize=16)
        ax.legend(loc=1)
        ax.grid=True; 
        #plot spiked insertion size statistics 
        bx=fig.add_subplot(212); 
        totalsizevals=[i[0] for i in self.actualInsertSize.values()]; 
        repeatsizevals=[i[1] for i in self.actualInsertSize.values()]; 
        sizeKeys=[i-barwidth for i in self.actualInsertSize.keys()]; 
        sizeValues=[i for i in totalsizevals]; 
        bx.bar(sizeKeys,sizeValues,width=barwidth,color='b',log=True,label='Spiked Insert.')
        bx.bar(self.actualInsertSize.keys(),repeatsizevals,width=barwidth,color='r',log=True,label='Spiked Repeat Insertion')
        bx.set_title(self.plotName+' Insertion Size\nSpiked Data',fontsize=18) 
        bx.set_xlabel('Insertion Length (in bases)',fontsize=16) 
        bx.set_ylabel('Number of Insertions',fontsize=16)
        plt.setp(bx.get_xticklabels(),fontsize=16) 
        plt.setp(bx.get_yticklabels(),fontsize=16) 
        bx.legend(loc=1); 
        bx.grid=True; 
        plt.savefig(self.plotName+'SpikingInsertionSize.png',bbox_inches=0)     

    def plotDeletionRate(self): 
        fig=plt.figure()
        plt.subplots_adjust(hspace=0.8)
        ax=fig.add_subplot(211); 
        for key in self.delProb: 
            self.delProb[key]=self.delProb[key]/float(self.pos[key])
        ax.bar(self.delProb.keys(),self.delProb.values(),color='b',label='Empirical',log=True);
        curveFitUpdated=dict(); 
        for key in self.curveFitDelProb: 
            if key in self.pos: 
                curveFitUpdated[key]=self.curveFitDelProb[key]/float(self.pos[key]); 

        ax.plot(curveFitUpdated.keys(),curveFitUpdated.values(),color='r',label='Regression',linewidth=3); 
        ax.set_ylim([10e-5,1])
        ax.set_title(self.plotName+' Deletion Probability\nBackground Data',fontsize=18) 
        ax.set_xlabel('Base Position Along Read',fontsize=16)
        ax.set_ylabel('Probability of Deletion',fontsize=16)
        plt.setp(ax.get_xticklabels(),fontsize=16) 
        plt.setp(ax.get_yticklabels(),fontsize=16)
        ax.legend(loc=4); 
        ax.grid=True; 

        bx=fig.add_subplot(212); 
        for key in self.actualDelPos: 
            self.actualDelPos[key]=self.actualDelPos[key]/float(self.actualPosCount[key]); 
        bx.bar(self.actualDelPos.keys(),self.actualDelPos.values(),color='b',log=True,label='Spiked')
        bx.set_ylim([10e-5,1]) 
        bx.set_title(self.plotName+' Deletion Probability\nSpiked Data',fontsize=18) 
        bx.set_xlabel('Base Position Along Read',fontsize=16) 
        bx.set_ylabel('Probability of Deletion',fontsize=16)
        plt.setp(bx.get_xticklabels(),fontsize=16) 
        plt.setp(bx.get_yticklabels(),fontsize=16) 
        bx.grid=True 
        plt.savefig(self.plotName+'SpikingDeletionRate.png',bbox_inches=0) 


    def plotDeletionSize(self): 
        fig=plt.figure()
        plt.subplots_adjust(hspace=0.8)
        #graph background insertion size statistics 
        barwidth=0.2
        sizeKeys=[i-barwidth for i in self.delSize.keys()]
        sizeValues=[i for i in self.delSize.values()] 
        ax=fig.add_subplot(211); 
        ax.bar(sizeKeys,sizeValues,width=barwidth,color='b',align='center',log=True,label='Back. Deletions.')
        ax.bar(self.delRepeat.keys(),self.delRepeat.values(),width=barwidth,color='r',align='center',log=True,label='Back. Repeat Deletions') 
        ax.grid=True
        
        ax.set_title(self.plotName+' Deletion Size\nBackground Data',fontsize=18)
        ax.set_xlabel('Deletion Length (in bases)',fontsize=16)
        ax.set_ylabel('Number of Deletions',fontsize=16)
        plt.setp(ax.get_xticklabels(),fontsize=16) 
        plt.setp(ax.get_yticklabels(),fontsize=16)
        ax.legend(loc=1)
        ax.grid=True; 
        #plot spiked insertion size statistics 
        bx=fig.add_subplot(212); 
        totalsizevals=[i[0] for i in self.actualDelSize.values()]; 
        repeatsizevals=[i[1] for i in self.actualDelSize.values()]; 
        sizeKeys=[i-barwidth for i in self.actualDelSize.keys()]; 
        sizeValues=[i for i in totalsizevals]; 
        bx.bar(sizeKeys,sizeValues,width=barwidth,color='b',log=True,label='Spiked Deletions')
        bx.bar(self.actualDelSize.keys(),repeatsizevals,width=barwidth,color='r',log=True,label='Spiked Repeat Deletions')
        bx.set_title(self.plotName+' Deletions Size\nSpiked Data',fontsize=18) 
        bx.set_xlabel('Deletion Length (in bases)',fontsize=16) 
        bx.set_ylabel('Number of Deletions',fontsize=16) 
        plt.setp(ax.get_xticklabels(),fontsize=16) 
        plt.setp(ax.get_yticklabels(),fontsize=16)
        bx.legend(loc=1);
        bx.grid=True; 
        plt.savefig(self.plotName+'SpikingDeletionSize.png',bbox_inches=0)     


    def plotMutationRate(self): 
        fig=plt.figure()
        plt.subplots_adjust(hspace=0.8)
        ax=fig.add_subplot(211); 
        for key in self.mutationProb: 
            self.mutationProb[key]=self.mutationProb[key]/float(self.pos[key])
        ax.bar(self.mutationProb.keys(),self.mutationProb.values(),color='b',log=True,label='Empirical'); 
        curveFitUpdated=dict(); 
        for key in self.curveFitMutationProb: 
            if key in self.pos: 
                curveFitUpdated[key]=self.curveFitMutationProb[key]/float(self.pos[key])
        ax.plot(curveFitUpdated.keys(),curveFitUpdated.values(),color='r',label='Regression',linewidth=3); 
        ax.set_ylim([10e-5,1])
        ax.set_title(self.plotName+' Mutation Probability\nBackground Data',fontsize=18) 
        ax.set_xlabel('Base Position Along Read',fontsize=16)
        ax.set_ylabel('Probability of\nSingle Base Mutation',fontsize=16)
        plt.setp(ax.get_xticklabels(),fontsize=16) 
        plt.setp(ax.get_yticklabels(),fontsize=16) 
        ax.legend(loc=4) 
        ax.grid=True; 

        bx=fig.add_subplot(212); 
        for key in self.actualMutPos: 
            self.actualMutPos[key]=self.actualMutPos[key]/float(self.actualPosCount[key]); 
        bx.bar(self.actualMutPos.keys(),self.actualMutPos.values(),color='b',log=True,label='Spiked')
        bx.set_ylim([10e-5,1]) 
        bx.set_title(self.plotName+' Mutation Probability\nSpiked Data',fontsize=18) 
        bx.set_xlabel('Base Position Along Read',fontsize=16) 
        bx.set_ylabel('Probability of\nSingle Base Mutation',fontsize=16) 
        plt.setp(bx.get_xticklabels(),fontsize=16) 
        plt.setp(bx.get_yticklabels(),fontsize=16)
        bx.grid=True 
        plt.savefig(self.plotName+'SpikingMutationRate.png',bbox_inches=0) 


    def plotAll(self): 
        self.plotReadLengthDist()
        self.plotInsertionRate()
        self.plotInsertionSize()
        self.plotDeletionRate()
        self.plotDeletionSize()
        self.plotMutationRate()

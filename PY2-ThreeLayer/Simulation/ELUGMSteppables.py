from __future__ import division
from PySteppables import *       
import math                      
import numpy                   
import CompuCell
import sys
import random
RNG=random.SystemRandom()  #draw true random sequence, overkill but why not?
#Motility Variables
CtoM=52  #cell to adhesion value, self-consistent across simulations as 26 but x2 due to double interface formation
BASAL=100 #baseline motility for L929, due to their extremely motile behaviour
SCF=0.5 #self-attenuator weighing basal motility vs loss of motility due to adhesion
#Self-Cutoff
ENDMCS=20000 #call runtime here directly
#Mitosis Variables
RADAVG=3 #average radius of the gaussian distribution to choose random radius
RADDEV=.5 #standard deviation of target radius, too low and division couples, too high and you'll lose cells at the start
MTFORCEMIN=-3*10**(-3.88) #negative mitosis driving force fluctuation, usually only need to change the exponential part
MTFORCEMAX=4*10**(-3.88)  #positive mitosis driving force fluctuation, usually change only the exponential part
#Signaling Variables
CONEXPSCF=10000 #Steady state expression of ligand expressed on a sender cell. This ligand is unaffected by signaling.
THETA=0 #time lag for expression of your constitutive, non-signaling affected ligand, start at 0 for simplicity, but can be adjusted depending on experiment results if known for generalizability
XI=1000 #controls how fast the sender cells reaches steady state for your constitutive, non-signaling affected ligand
FASTAPPROX=5000 #force approx for function of above variables at the time step, saves calling the mcs and doing the caluclation, purely computational speed effeciency

ALPHAYG=1 #controls how much your reporter synthesis magnitude due to signal S; can be set to 1 if decay is set properly
BETAYG=1750 #threshold of signal required to generate a response in your cell due to signaling
EPSILONYG=1000 #modulates how sharp the response is due to signaling, can turn synthesis to linear or heavi-side theta like if desired
KAPPAYG=25000 #general decay constant
THRESHOLDUPYG=5263 #activation threshold to change state  
THRESHOLDDOYG=0 #deactivation threshold to revert state, as no clear deactivation/unsorting occured in reference experiments

ALPHABR=1 #controls how much your reporter synthesis magnitude due to signal S; can be set to 1 if decay is set properly
BETABR=921.181 #threshold of signal required to generate a response in your cell due to signaling
EPSILONBR=526.389 #modulates how sharp the response is due to signaling, can turn synthesis to linear or heavi-side theta like if desired
KAPPABR=25000 #general decay constant
THRESHOLDUPBR=5301 #activation threshold to change state, the choice of paramters in this paragraph render it similar to that of the previous paragraph
THRESHOLDDOBR=0 #deactivation threshold to revert state, as no clear deactivation/unsorting occured in reference experiments

#Single Cell Trace Variables
MARKEDCELLS=[22,233,51,228] # ID of cells to track if you desire single cell points tracked, change to fit setup

#Sampling and Comp Speed
RESOL=100 #Data sampling rate, choose to satisfy nyquist theorem if necessary
USEDNODES=32 #Choose a power of 2, otherwise the grids overlap and your simulation will eventually randomly crash, follow the recommendations given in the manual by developers


class ELUGMSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):

        self.pW1 = self.addNewPlotWindow(
            _title='Calibrate',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Psi',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW1.addPlot('IBR', _style='Dots', _color='green', _size=3) 
        self.pW1.addPlot('IYG', _style='Dots', _color='red', _size=3)
        self.pW1.addPlot('SBR', _style='Dots', _color='yellow', _size=5)        
        self.pW1.addPlot('SYG', _style='Dots', _color='white', _size=5)   #plots for morphospace/hetergenity measure
        
        self.pW2 = self.addNewPlotWindow(
            _title='Psi',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Psi',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW2.addPlot('Hamiltonian', _style='Dots', _color='white', _size=3) #measure the system total energy
 
        self.pW3 = self.addNewPlotWindow(
            _title='Types',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Count',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW3.addPlot('Y', _style='Dots', _color='gray', _size=3)
        self.pW3.addPlot('G', _style='Dots', _color='green', _size=3)
        self.pW3.addPlot('B', _style='Dots', _color='blue', _size=3)   
        self.pW3.addPlot('R', _style='Dots', _color='red', _size=3)
        self.pW3.addPlot('AG', _style='Dots', _color='green', _size=5)
        self.pW3.addPlot('AR', _style='Dots', _color='red', _size=5) #count the types of cells in the simulation, along with how mny activate due to signaling
        
        self.pW4 = self.addNewPlotWindow(
            _title='Point System',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Count',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW4.addPlot('YG', _style='Dots', _color='green', _size=3)
        self.pW4.addPlot('BR', _style='Dots', _color='blue', _size=3) #measure average points per cell type, can be tied to average flurosence intensity
 
        self.pW5 = self.addNewPlotWindow(
            _title='Single Point System',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Count',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW5.addPlot('1', _style='Dots', _color='blue', _size=3)
        self.pW5.addPlot('2', _style='Dots', _color='gray', _size=3)
        self.pW5.addPlot('3', _style='Dots', _color='green', _size=3)   
        self.pW5.addPlot('4', _style='Dots', _color='red', _size=3)
        self.pW5.addPlot('1T', _style='Dots', _color='blue', _size=3)
        self.pW5.addPlot('2T', _style='Dots', _color='gray', _size=3)
        self.pW5.addPlot('3T', _style='Dots', _color='green', _size=3)   
        self.pW5.addPlot('4T', _style='Dots', _color='red', _size=3)  #single cell point traces 

        self.pW6 = self.addNewPlotWindow(
            _title='Sphericity',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Count',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW6.addPlot('BF', _style='Dots', _color='blue', _size=3)
        self.pW6.addPlot('MR', _style='Dots', _color='green', _size=3)  #bright field and merged field sphericity      
                
        global YtoY,YtoG,GtoY,YtoB,BtoY,YtoR,RtoY,GtoG,GtoB,BtoG,GtoR,RtoG,BtoB,BtoR,RtoB,RtoR #the adhesion matrix, call these values and store them for motility code
        YtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','Y']))
        YtoG=GtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','G']))
        YtoB=BtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','B']))
        YtoR=RtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','R']))
        GtoG=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','G','Type2','G']))
        GtoB=BtoG=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','G','Type2','B']))
        GtoR=RtoG=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','G','Type2','R']))
        BtoB=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','B','Type2','B']))
        BtoR=RtoB=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','B','Type2','R']))
        RtoR=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','R','Type2','R']))                        

        for cell in self.cellList:
            cell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV) #assign the cells a random target radius
            cell.lambdaSurface=2.5                    #temporary value, will be changed in later code
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2  #spherical surface area
            cell.lambdaVolume=2.5                     #temporary value, will be changed in later code
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
            cell.dict["PTS"]=[0]                      #initial points for cell, feel free to randomize if desired
            cell.dict["P"]=[0,0]                      #activation counter, counts how many cells are active due to signaling at any given time

    def step(self,mcs):                 
        NUMTY=0 #number of type Y
        NUMTG=0 #number of type G
        NUMTB=0 #number of type B
        NUMTR=0 #number of type R
        
        YGPTS=0 #number of Y+G points
        BRPTS=0 #number of B+R points
        
        SYPSI=0 #system hamiltonian over all interaction over configuration

        CSAYGBR=0 #common surface area of YG to BR
        CSABRYG=0 #common surface area of BR to YG
        CSAYGYG=0 #common surface area of YG to YG
        CSABRBR=0 #common surface area of BR to BR
        YGCBR=0 #how many Y or G cells are in contact with B or R cell?
        BRCYG=0 #how many B or R cells are in contact with Y or G cells?
        YGCYG=0 #how many Y or G cells are in contact with Y or G cells?
        BRCBR=0 #how many B or R cells are in contact with B or R cells?
        
        SUMBFSF=0 #total bright field surface area
        SUMBFVL=0 #total bright field volume
        SUMMRSF=0 #total color field surface area
        SUMMRVL=0 #total color field volume
        
        NAR=0 #number of activated red cells due to signaling
        NAG=0 #number of activated green cells due to signaling
        
        if mcs==1:
            self.changeNumberOfWorkNodes(USEDNODES) #set to necessary computational nodes                  

        for cell in self.cellList: #iterate over cell list
            CSAY=0 #each cell detect how much sirface area it shares with Y cells
            CSAG=0 #each cell detect how much sirface area it shares with G cells
            CSAB=0 #each cell detect how much sirface area it shares with B cells
            CSAR=0 #each cell detect how much sirface area it shares with R cells
            CSAM=0 #each cell detect how much sirface area it shares with medium
            
            PTSY=0 #each cell gains points from neighbor type Y
            PTSG=0 #each cell gains points from neighbor type G
            PTSB=0 #each cell gains points from neighbor type B
            PTSR=0 #each cell gains points from neighbor type R
            DTRES=0 #change in reporter due to signal S

            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell): #iterate for each cell its neighbors
                if neighbor is None: #none type refers to medium 
                    continue
                if neighbor.type==1: #gray cells
                    CSAY+=commonSurfaceArea #total common surface area with gray cells
                    PTSY+=0                 #points gray cells send to receiver
                if neighbor.type==2: #green cells
                    CSAG+=commonSurfaceArea #total common surface area with green cells
                    PTSG+=commonSurfaceArea*neighbor.dict["PTS"][0]/(neighbor.surface) #PHI and L as in the text
                if neighbor.type==3: #blue cells
                    CSAB+=commonSurfaceArea #total common surface area with blue cells                     
                    if mcs<FASTAPPROX: #constitutive ligand function times PHI
                        PTSB+=commonSurfaceArea*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface
                    else: #fast computational approximation
                        PTSB+=commonSurfaceArea*CONEXPSCF/neighbor.surface
                if neighbor.type==4: #red cells
                    CSAR+=commonSurfaceArea #total common surface area with red cells 
                    if mcs<FASTAPPROX: #constitutive ligand function times PHI
                        PTSR+=commonSurfaceArea*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface
                    else: #fast computational approximation
                        PTSR+=commonSurfaceArea*CONEXPSCF/neighbor.surface
            CSAM=cell.surface-(CSAY+CSAG+CSAB+CSAR) #alternative method to calculate common surface area with medium                 

# VETTING CODE                                        
#             if cell.id==3:
#                 for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell): #iterate for each cell its neighbors
#                     if neighbor:
#                         print "NID", neighbor.id, "TY", neighbor.type, "NCSA", commonSurfaceArea, "DICT", neighbor.dict["PTS"][0], "NS", neighbor.surface
#                 print "CID", cell.id, "CT",cell.type, "CS", cell.surface
#                 print CSAY, CSAG, PTSG, CSAB, PTSB, CSAR, PTSR, CSAM
           
            if (cell.type==1 or cell.type==2): #which cells receive what type of signal
                DTRES=(1/(ALPHAYG+math.exp(-((PTSR+PTSB)-BETAYG)/EPSILONYG)))-(1/KAPPAYG)*cell.dict["PTS"][0] #del reporter
                cell.dict["PTS"][0]+=DTRES #update total reporter
            if (cell.type==3 or cell.type==4): #which cells receive what type of signal
                DTRES=(1/(ALPHABR+math.exp(-((PTSY+PTSG)-BETABR)/EPSILONBR)))-(1/KAPPABR)*cell.dict["PTS"][0] #del reporter
                cell.dict["PTS"][0]+=DTRES #update total reporter
                
            if cell.type==1: #change Y cell state
                if cell.dict["PTS"][0]>=THRESHOLDUPYG:
                    cell.type=2
                    cell.dict["P"][0]=1 #iterate one for activating
            if cell.type==2: #change G cell state
                if cell.dict["PTS"][0]<THRESHOLDDOYG:
                    cell.type=1
                    cell.dict["P"][0]=0 #set to 0 for deactivating
            if cell.type==3: #change B cell state
                if cell.dict["PTS"][0]>=THRESHOLDUPBR:
                    cell.type=4
                    cell.dict["P"][1]=1 #iterate one for activating
            if cell.type==4: #change R cell state
                if cell.dict["PTS"][0]<THRESHOLDDOBR:
                    cell.type=3
                    cell.dict["P"][1]=0 #set to 0 for deactivating

            if cell.type==1: #gray cells
                SUMBFSF+=CSAM #grays cell surface area count under bright field
                SUMBFVL+=cell.volume #gray cell volume count under bright field              
                CSAYGBR+=(CSAB+CSAR)/cell.surface #total common surface area of YG cells to BR cells normailzed to YG cell surface
                if (CSAB+CSAR)>0:                 # count the number of YG cells in contact with BR cells
                    YGCBR+=1
                CSAYGYG+=(CSAY+CSAG)/cell.surface # count the number of YG cells in contact with YG cells normailzed to YG cell surface Think of as proportions
                if (CSAY+CSAG)>0:                 # count the number of Yg cells in contact with YG cells
                    YGCYG+=1                    
                cell.lambdaSurface=2.2            #change depending on cell adhesitivity
                cell.lambdaVolume=2.2             #change depending on cell adhesitivity  
                NUMTY+=1                          #count the number if gray cells
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+YtoY*CSAY+YtoG*CSAG+YtoB*CSAB+YtoR*CSAR)/cell.surface #corrected cell motility, tune based on adhesive neighbors, vetted
                YGPTS+=cell.dict["PTS"][0] #count number of points YG cells have

            if cell.type==2: #green cells
                NAG+=cell.dict["P"][0] #number of activated green cells due to activation at any given time step
                SUMBFSF+=CSAM          #green cells surface area under bright field
                SUMBFVL+=cell.volume   #green cell volume under bright field
                SUMMRSF+=(CSAM+CSAY+CSAB) #green cells count under color field, medium blue and gray are invisible
                SUMMRVL+=cell.volume      #green cells contribute to color field volume         
                CSAYGBR+=(CSAB+CSAR)/cell.surface #total common surface area of YG cells to BR cells normailzed to YG cell surface
                if (CSAB+CSAR)>0:                 # count the number of YG cells in contact with BR cells
                    YGCBR+=1
                CSAYGYG+=(CSAY+CSAG)/cell.surface # count the number of YG cells in contact with YG cells normailzed to YG cell surface Think of as proportions
                if (CSAY+CSAG)>0:                 # count the number of Yg cells in contact with YG cells
                    YGCYG+=1
                cell.lambdaSurface=1.0            #change depending on cell adhesitivity
                cell.lambdaVolume=1.0             #change depending on cell adhesitivity    
                NUMTG+=1                          #count the number of green cells
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+GtoY*CSAY+GtoG*CSAG+GtoB*CSAB+GtoR*CSAR)/cell.surface #corrected cell motility, tune based on adhesive neighbors, vetted
                YGPTS+=cell.dict["PTS"][0] #count number of points YG cells have
             
            if cell.type==3: #blue cells            
                SUMBFSF+=CSAM # blue surface area contributes to bright field surface area
                SUMBFVL+=cell.volume #blue volume contributes to bright field volume                
                CSABRYG+=(CSAY+CSAG)/cell.surface #total common surface area of BR cells to YG cells normalized
                if (CSAY+CSAG)>0:                 # count number of BR cells in contact with YG cels
                    BRCYG+=1
                CSABRBR+=(CSAB+CSAR)/cell.surface #total common surface area of BR cells to BR cells normalized
                if (CSAB+CSAR)>0:                 # count number of BR cells in contact with BR cels
                    BRCBR+=1                    
                cell.lambdaSurface=2.2           #change depending on cell adhesitivity
                cell.lambdaVolume=2.2            #change depending on cell adhesitivity      
                NUMTB+=1                         #count number of blue cells
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+BtoY*CSAY+BtoG*CSAG+BtoB*CSAB+BtoR*CSAR)/cell.surface # corrected cell motility, vetted
                BRPTS+=cell.dict["PTS"][0] #count number of points BR cells have
                
            if cell.type==4: #red cells
                NAR+=cell.dict["P"][1] #number of activated red cells due to activation at any given time step
                SUMBFSF+=CSAM          #red cells are visible under bright field and thus conribute surface area
                SUMBFVL+=cell.volume   #red cell volume contributes to bright field                
                SUMMRSF+=(CSAM+CSAY+CSAB) #red cells contribute to color field, medium blue and gray are invisible
                SUMMRVL+=cell.volume      #red cells contribute to color field volume             
                CSABRYG+=(CSAY+CSAG)/cell.surface #common surface area of blue cells to gray cells
                if (CSAY+CSAG)>0:                 #count number of blue red cells that see y g cells
                    BRCYG+=1
                CSABRBR+=(CSAB+CSAR)/cell.surface #common surface area of BR cells that see BR cells
                if (CSAB+CSAR)>0:                 #count the number of BR cells that are in contact with BRBR cells
                    BRCBR+=1                      
                cell.lambdaSurface=1.0            #change depending on cell adhesitivity
                cell.lambdaVolume=1.0             #change depending on cell adhesitivity      
                NUMTR+=1                          #count number of red cells
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+RtoY*CSAY+RtoG*CSAG+RtoB*CSAB+RtoR*CSAR)/cell.surface #corrected cell motility, vetted
                BRPTS+=cell.dict["PTS"][0] #count number of points of BR cells
                    
            SYPSI+=cell.fluctAmpl #sensitive measure to sorting events
            
#vetting code
            #print cell.id, cell.type, cell.fluctAmpl

            if mcs%RESOL==0: #record points for signle cells traces        
                if cell.id==MARKEDCELLS[0]:
                    self.pW5.addDataPoint("1", mcs, cell.dict["PTS"][0])        
                    self.pW5.addDataPoint("1T", mcs, cell.type) 
                if cell.id==MARKEDCELLS[1]:
                    self.pW5.addDataPoint("2", mcs, cell.dict["PTS"][0])        
                    self.pW5.addDataPoint("2T", mcs, cell.type)
                if cell.id==MARKEDCELLS[2]:
                    self.pW5.addDataPoint("3", mcs, cell.dict["PTS"][0])        
                    self.pW5.addDataPoint("3T", mcs, cell.type)                           
                if cell.id==MARKEDCELLS[3]:
                    self.pW5.addDataPoint("4", mcs, cell.dict["PTS"][0])        
                    self.pW5.addDataPoint("4T", mcs, cell.type)  
        
        if mcs%RESOL==0:
            
            self.pW2.addDataPoint("Hamiltonian", mcs, SYPSI)
            
            self.pW3.addDataPoint("Y", mcs, NUMTY)
            self.pW3.addDataPoint("G", mcs, NUMTG) 
            self.pW3.addDataPoint("B", mcs, NUMTB)   
            self.pW3.addDataPoint("R", mcs, NUMTR)
            self.pW3.addDataPoint("AR", mcs, NAR)
            self.pW3.addDataPoint("AG", mcs, NAG)

            self.pW4.addDataPoint("YG", mcs, YGPTS/(NUMTY+NUMTG))
            self.pW4.addDataPoint("BR", mcs, BRPTS/(NUMTB+NUMTR))
            
            if YGCYG==0:
                self.pW1.addDataPoint("IYG", mcs, 0)
            if YGCYG>0:
                self.pW1.addDataPoint("IYG", mcs, CSAYGYG/(YGCYG))
            if BRCBR==0:
                self.pW1.addDataPoint("IBR", mcs, 0)
            if BRCBR>0:
                self.pW1.addDataPoint("IBR", mcs, CSABRBR/(BRCBR))
                
            if YGCBR==0:
                self.pW1.addDataPoint("SYG", mcs, 0)
            if YGCBR>0:
                self.pW1.addDataPoint("SYG", mcs, CSAYGBR/(YGCBR))
            if BRCYG==0:
                self.pW1.addDataPoint("SBR", mcs, 0)
            if BRCYG>0:
                self.pW1.addDataPoint("SBR", mcs, CSABRYG/(BRCYG))
            
            if SUMBFSF==0:
                self.pW6.addDataPoint("BF", mcs, 0)                
            if SUMBFSF>0: 
                self.pW6.addDataPoint("BF", mcs, ((math.pi**(1/3))*(6*SUMBFVL)**(2/3))/SUMBFSF)
            if SUMMRSF==0:
                self.pW6.addDataPoint("MR", mcs, 0)
            if SUMMRSF>0:
                self.pW6.addDataPoint("MR", mcs, ((math.pi**(1/3))*(6*SUMMRVL)**(2/3))/SUMMRSF)
             
            if mcs==ENDMCS:
                fileName = "Calibrate" + str(mcs) + ".txt"
                self.pW1.savePlotAsData(fileName)                
                fileName = "PSI" + str(mcs) + ".txt"
                self.pW2.savePlotAsData(fileName)
                fileName = "FOU" + str(mcs) + ".txt"
                self.pW3.savePlotAsData(fileName)
                fileName = "SIG" + str(mcs) + ".txt"
                self.pW4.savePlotAsData(fileName)
                fileName = "SCSIG" + str(mcs) + ".txt"
                self.pW5.savePlotAsData(fileName)                 
                fileName = "Sphericity" + str(mcs) + ".txt"
                self.pW6.savePlotAsData(fileName)                 
                self.stopSimulation()                 
    def finish(self):
        pass

from PySteppablesExamples import MitosisSteppableBase

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        self.setParentChildPositionFlag(0) #randomize child cell position, see developer manual
    def step(self,mcs):        
        cells_to_divide=[]          #gen cells to divide list
        for cell in self.cellList:
            cell.dict["RDM"]+=RNG.uniform(MTFORCEMIN,MTFORCEMAX) #make cells grow in target radius by this much
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2 #spherical surface area
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
            if cell.volume>2*(4/3)*math.pi*RADAVG**3: #divide at two times the mean radius initialized with               
                cells_to_divide.append(cell)           #add these cells to divide list
                
        for cell in cells_to_divide:
            self.divideCellRandomOrientation(cell)  #divide the cells

    def updateAttributes(self):
        self.parentCell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV) #reassign new target radius
        self.parentCell.targetVolume=(4/3)*math.pi*self.parentCell.dict["RDM"]**3 #new target volume
        self.parentCell.targetSurface=4*math.pi*self.parentCell.dict["RDM"]**2 #new target surface area
        self.cloneParent2Child()  #copy characterstics to child cell, indlucig signaling
        self.childCell.dict["P"][0]=0 #reset the activation counter, we dont care about cells from activated parent
        self.childCell.dict["P"][1]=0 #reset the activation counter, we dont care about cells from activated parent
        
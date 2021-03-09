#from __future__ import division
from cc3d.core.PySteppables import *
import math                      
import numpy                   
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

    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
    def start(self):

        global YtoY,YtoG,GtoY,YtoB,BtoY,YtoR,RtoY,GtoG,GtoB,BtoG,GtoR,RtoG,BtoB,BtoR,RtoB,RtoR #the adhesion matrix, call these values and store them for motility code
        YtoYC1=self.get_xml_element('YtoY')
        YtoGC1=GtoYC1=self.get_xml_element('YtoG')
        YtoBC1=BtoYC1=self.get_xml_element('YtoB')
        YtoRC1=RtoYC1=self.get_xml_element('YtoR')
        GtoGC1=self.get_xml_element('GtoG')
        GtoBC1=BtoGC1=self.get_xml_element('GtoB')
        GtoRC1=RtoGC1=self.get_xml_element('GtoR')
        BtoBC1=self.get_xml_element('BtoB')
        BtoRC1=RtoBC1=self.get_xml_element('BtoR')
        RtoRC1=self.get_xml_element('RtoR')
        
        YtoY=float(YtoYC1.cdata)
        YtoG=GtoY=float(YtoGC1.cdata)
        YtoB=BtoY=float(YtoBC1.cdata)
        YtoR=RtoY=float(YtoRC1.cdata)
        GtoG=float(GtoGC1.cdata)
        GtoB=BtoG=float(GtoBC1.cdata)
        GtoR=RtoG=float(GtoRC1.cdata)
        BtoB=float(BtoBC1.cdata)
        BtoR=RtoB=float(BtoRC1.cdata)
        RtoR=float(RtoRC1.cdata)

        for cell in self.cell_list:
            cell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV) #assign the cells a random target radius
            cell.lambdaSurface=2.5                    #temporary value, will be changed in later code
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2  #spherical surface area
            cell.lambdaVolume=2.5                     #temporary value, will be changed in later code
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
            cell.dict["PTS"]=[0]                      #initial points for cell, feel free to randomize if desired
            cell.dict["P"]=[0,0]                      #activation counter, counts how many cells are active due to signaling at any given time

    def step(self,mcs):                         
        if mcs==1:
            self.change_number_of_work_nodes(USEDNODES) #set to necessary computational nodes                  

        for cell in self.cell_list: #iterate over cell list
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

            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell): #iterate for each cell its neighbors
                if neighbor is None: #none type refers to medium 
                    continue
                if neighbor.type==1: #gray cells
                    CSAY+=common_surface_area #total common surface area with gray cells
                    PTSY+=0                 #points gray cells send to receiver
                if neighbor.type==2: #green cells
                    CSAG+=common_surface_area #total common surface area with green cells
                    PTSG+=common_surface_area*neighbor.dict["PTS"][0]/(neighbor.surface) #PHI and L as in the text
                if neighbor.type==3: #blue cells
                    CSAB+=common_surface_area #total common surface area with blue cells                     
                    if mcs<FASTAPPROX: #constitutive ligand function times PHI
                        PTSB+=common_surface_area*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface
                    else: #fast computational approximation
                        PTSB+=common_surface_area*CONEXPSCF/neighbor.surface
                if neighbor.type==4: #red cells
                    CSAR+=common_surface_area #total common surface area with red cells 
                    if mcs<FASTAPPROX: #constitutive ligand function times PHI
                        PTSR+=common_surface_area*(CONEXPSCF/(1+math.exp(-(mcs-THETA)/XI)))/neighbor.surface
                    else: #fast computational approximation
                        PTSR+=common_surface_area*CONEXPSCF/neighbor.surface
            CSAM=cell.surface-(CSAY+CSAG+CSAB+CSAR) #alternative method to calculate common surface area with medium                 

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
                cell.lambdaSurface=2.2            #change depending on cell adhesitivity
                cell.lambdaVolume=2.2             #change depending on cell adhesitivity                          #count the number if gray cells
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+YtoY*CSAY+YtoG*CSAG+YtoB*CSAB+YtoR*CSAR)/cell.surface #corrected cell motility, tune based on adhesive neighbors, vetted

            if cell.type==2: #green cells
                cell.lambdaSurface=1.0            #change depending on cell adhesitivity
                cell.lambdaVolume=1.0             #change depending on cell adhesitivity    
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+GtoY*CSAY+GtoG*CSAG+GtoB*CSAB+GtoR*CSAR)/cell.surface #corrected cell motility, tune based on adhesive neighbors, vetted
            
            if cell.type==3: #blue cells                              
                cell.lambdaSurface=2.2           #change depending on cell adhesitivity
                cell.lambdaVolume=2.2            #change depending on cell adhesitivity      
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+BtoY*CSAY+BtoG*CSAG+BtoB*CSAB+BtoR*CSAR)/cell.surface # corrected cell motility, vetted
               
            if cell.type==4: #red cells                    
                cell.lambdaSurface=1.0            #change depending on cell adhesitivity
                cell.lambdaVolume=1.0             #change depending on cell adhesitivity      
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+RtoY*CSAY+RtoG*CSAG+RtoB*CSAB+RtoR*CSAR)/cell.surface #corrected cell motility, vetted
              
    def finish(self):
        pass



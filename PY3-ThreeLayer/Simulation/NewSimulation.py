from cc3d import CompuCellSetup
        
from .NewSimulationSteppables import ELUGMSteppable
CompuCellSetup.register_steppable(steppable=ELUGMSteppable(frequency=1))

CompuCellSetup.run()

from cc3d import CompuCellSetup
from BaseSteppables import BaseSteppable
CompuCellSetup.register_steppable(steppable=BaseSteppable(frequency=1))
from BaseSteppables import MitosisSteppable
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))
CompuCellSetup.run()

# from PrirodaOptimizer import Conformer, PrirodaOptimizer
from PrirodaOptimiser.conformer import Conformer
from PrirodaOptimiser.transformer import PrirodaOptimizer


c = Conformer.from_xyz('./PrirodaOptimiser/data/Alanine.xyz')
p = PrirodaOptimizer()

print(p._prepare_input(c, tmp_dir="/tmp"))

from PrirodaOptimizer import Conformer, Priroda_Optimizer


c = Conformer.from_xyz('Alanine.xyz')
p = Priroda_Optimizer()

print(p._prepare_input(c))
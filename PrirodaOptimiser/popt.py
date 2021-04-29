from itertools import islice

# class PrirodaOptimizer():
    
#     def __init__(self, relative=False, basis="L1", memory=200, mp_dir="/tmp", tmp_ram=0):
#         self.basis = basis
#         self.memory = memory
#         self.mp_dir = mp_dir
#         self.tmp_ram = tmp_ram

#     def fit(self, x, y=None):
#         return self

#     def transform(self, x, y=None):
#         """Точка входа, на вход список молекул с начальным приближением, выход оптимизированные молекулы"""
#         ...


def validator(f):
    def inner(self, atoms, coords, charge, multiplicity):
        if not (-8 < charge < 8):
            raise ValueError('Not valid charge')
        if not multiplicity in (1,2,3,4,5):
            raise ValueError('Not valid multiplicity')
        if len(atoms)!=len(coords):
            raise BaseException('Number of atoms is not equals number of coords')
        valid_atoms = {"C", "N", "O", "F", "Cl", "Br", "I", "H", "S", "P", "B", "Si"}
        for atom in atoms:
            if atom not in valid_atoms:
                raise ValueError('Not valid type of atom')

        return f(self, atoms, coords, charge, multiplicity)
    return inner


class Conformer:
    @validator
    def __init__(self, atoms, coords, charge, multiplicity):
        self.atoms = atoms
        self.coords = coords
        self.charge = charge
        self.multiplicity = multiplicity
    
    
    @classmethod
    def from_xyz(cls, file, charge=0, multiplicity=1):  
        
        atoms, coords= [], []
        
        if isinstance(file, str):
            inp = open(file, 'r')
        else:
            inp = file
        
        numAtoms = int(inp.readline().rstrip('\n'))
        for line in islice(inp, 1, numAtoms+1):
            atom, *coord = line.rstrip('\n').split()
            atoms.append(atom)
            coords.append(tuple(map(float, coord)))
            
        if inp:
            inp.close()
            
        return cls(tuple(atoms), tuple(coords), charge, multiplicity)


if __name__ == '__main__':
    with open('data/Alanine.xyz', 'r') as inp:
        test = Conformer.from_xyz(inp)
    print(test.__dict__)
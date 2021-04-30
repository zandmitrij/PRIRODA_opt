from itertools import islice


class Validator():
    def __get__(self, obj, objtype=None):
        return getattr(obj, self.name)

    def __set_name__(self, owner, name):
        self.name = '_' + name


class Atom_Validator(Validator):
    atoms = {'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I'}

    def __set__(self, obj, value):
        for i in value:
            if i not in self.atoms:
                raise ValueError('Not valid type of atom')
        setattr(obj, self.name, value)


class Coords_Validator(Validator):  # проверка на все что угодно: тюплы, флоуты, длина атомов и координат
    def __init__(self, atoms):
        self.atoms = atoms
        # сначала атомы, потом координаты

    def __set__(self, obj, value):
        if not isinstance(value, tuple):
            ValueError('Not valid type of value')
        for i in value:
            if not (len(i)) != 3:
                raise ValueError
            for j in i:
                if not isinstance(j, float):
                    raise TypeError
        # if len(getattr(obj)) != len(self.atoms):  # ???
        #      raise BaseException('Number of atoms is not equals number of coords')
        # getatt(obj) ??? len(self.atoms) ???

        setattr(obj, self.name, value)


class Charge_Validator(Validator):
    def __set__(self, obj, value):
        if not isinstance(value, int):
            raise ValueError('Not valid type of charge')
        if not -8 < value < 8:
            raise ValueError('Not valid charge')
        setattr(obj, self.name, value)


class Multiplicity_Validator(Validator):
    def __set__(self, obj, value):
        if not isinstance(value, int):
            raise ValueError('Not valid type of multiplicity')
        if value < 1:
            raise ValueError('Not valid multiplicity')
        setattr(obj, self.name, value)


class Conformer:
    atoms = Atom_Validator()
    charge = Charge_Validator()
    multiplicity = Multiplicity_Validator()
    coords = Coords_Validator('atoms')

    def __init__(self, atoms: list, coords: tuple, charge: int, multiplicity: int):
        self.atoms = atoms
        self.coords = coords
        self.charge = charge
        self.multiplicity = multiplicity

    @classmethod
    def from_xyz(cls, file, charge=0, multiplicity=1):
        # file - pathlib / путь до файла / открытый файл
        # parser xyz
        # itertools islice(generator, n)
        # line.split() = sfff symbol,x,y,z

        atoms, coords = [], []
        file_open = False
        if isinstance(file, str):
            inp = open(file, 'r')
            file_open = True
        else:
            inp = file

        numAtoms = int(inp.readline().rstrip('\n'))
        for line in islice(inp, 1, numAtoms + 1):
            atom, *coord = line.rstrip('\n').split()
            atoms.append(atom)
            coords.append(tuple(map(float, coord)))

        if file_open:
            inp.close()

        return cls(tuple(atoms), tuple(coords), charge, multiplicity)


# with open('../Alanine.xyz', 'r') as input:
#     test = Conformer.from_xyz(input)
# print(test.__dict__)
# print(test.atoms)
# print(test.coords)

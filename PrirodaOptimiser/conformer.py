from itertools import islice


class _Validator():
    def __get__(self, obj, objtype=None):
        return getattr(obj, self.name)

    def __set_name__(self, owner, name):
        self.name = '_' + name


class AtomValidator(_Validator):
    atoms = {'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I'}

    def __set__(self, obj, value):
        for i in value:
            if i not in self.atoms:
                raise ValueError('Not valid type of atom')
        setattr(obj, self.name, value)


class CoordsValidator(_Validator):
    def __set__(self, instance, value):
        if not isinstance(value, tuple):
            raise TypeError(f'{value}, не tuple')
        for i in value:

            if not isinstance(i, tuple):
                raise TypeError(f'{i}, не tuple')
            if len(i) != 3:
                raise ValueError(f'{i}, не 3 значания')

            for temp in i:
                if not isinstance(temp, float):
                    raise TypeError
        setattr(instance, self.name, value)


class ChargeValidator(_Validator):
    def __set__(self, obj, value):
        if not isinstance(value, int):
            raise ValueError('Not valid type of charge')
        if not -8 < value < 8:
            raise ValueError('Not valid charge')
        setattr(obj, self.name, value)


class MultiplicityValidator(_Validator):
    def __set__(self, obj, value):
        if not isinstance(value, int):
            raise ValueError('Not valid type of multiplicity')
        if value < 1:
            raise ValueError('Not valid multiplicity')
        setattr(obj, self.name, value)


class Conformer:
    atoms = AtomValidator()
    charge = ChargeValidator()
    multiplicity = MultiplicityValidator()
    coords = CoordsValidator()

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

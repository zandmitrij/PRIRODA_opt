from itertools import islice
from .atom_map import atoms, atom_map, reverse_atom_map
from pathlib import Path
from typing import Tuple, Union, TextIO, Optional
from CGRtools import XYZRead
# import rdkit

class _Validator:
    def __get__(self, obj, objtype=None):
        return getattr(obj, self.name)

    def __set_name__(self, owner, name):
        self.name = '_' + name


class AtomValidator(_Validator):
    def __init__(self):
        self.atoms = atoms

    def __set__(self, obj, value):
        for i in value:
            if i not in self.atoms:
                raise ValueError('Not valid type of atom!')
        setattr(obj, self.name, value)


class CoordsValidator(_Validator):
    def __set__(self, instance, value):
        if not isinstance(value, tuple):
            raise TypeError("Not a tuple!")
        for i in value:

            if not isinstance(i, tuple):
                raise TypeError("Not a tuple!")
            if len(i) != 3:
                raise ValueError("Not three values received!")

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


class EnergyValidator(_Validator):
    def __set__(self, obj, value):
        if not isinstance(value, float):
            raise ValueError('Not valid type of energy')
        setattr(obj, self.name, value)


class HessianValidator(_Validator):
    def __set__(self, obj, value):
        if not isinstance(value, bool):
            raise ValueError('Not valid type of Hessian')
        setattr(obj, self.name, value)


class Conformer:
    atoms = AtomValidator()
    charge = ChargeValidator()
    multiplicity = MultiplicityValidator()
    coords = CoordsValidator()
    energy = EnergyValidator()
    hessian = HessianValidator()

    def __init__(self, atoms: Tuple[str, ...], coords: Tuple[Tuple[float, float, float], ...],
                 charge: int, multiplicity: int, hessian: bool = False, energy: float = 0.):
        self.atoms = atoms
        self.coords = coords
        self.charge = charge
        self.multiplicity = multiplicity
        self.hessian = hessian
        self.energy = energy

    @classmethod
    def from_xyz(cls, file: Union[str, Path, TextIO], charge: int = 0, multiplicity: int = 1):
        atoms, coords = [], []
        if isinstance(file, str):
            inp = open(file)
            file_open = True
        elif isinstance(file, Path):
            inp = file.open()
            file_open = True
        else:
            inp = file
            file_open = False

        num_atoms = int(inp.readline().strip())
        for line in islice(inp, 1, num_atoms + 1):
            atom, *coord = line.split()
            atoms.append(atom)
            coords.append(tuple(map(float, coord)))

        if file_open:
            inp.close()

        return cls(tuple(atoms), tuple(coords), charge, multiplicity)


    @classmethod
    def from_rdkit(cls, mol, multiplicity: int = 1, conformer: int = 0):

        # atoms: Iterator[atom] = mol.GetAtoms()
        # atom_num: int = atom.GetAtomicNum()
        # atom_charge: int = atom.GetFormalCharge()
        # conformers: Tuple[conformer, ...] = mol.GetConformers()
        # xyz: Iterable[Iterable[float]] = conformers[conformer].GetPositions()

        atoms, coords = [], []
        charge = 0
        for i in mol.GetAtoms():
            atoms.append(i.GetSymbol())
            charge += i.GetFormalCharge()  # заряд считается как сумма всех зарядов

        conformers = mol.GetConformers()
        coords = tuple(map(lambda x: tuple(x), conformers[conformer].GetPositions()))
        # берем конфермер по айди, по позиции он вытаскивает координаты всех атомов.
        # Получаем массив, преобразовываем в tuple
        # map принимает lambda функцию ко всем координатам всех атомов

        return cls(tuple(atoms), coords, charge, multiplicity)

    def to_xyz(self, file: Union[str, Path, TextIO]):
        if isinstance(file, str):
            out = open(file, 'w')
            file_open = True
        elif isinstance(file, Path):
            out = file.open()
            file_open = True
        else:
            out = file
            file_open = False

        atoms = self.atoms
        coords = self.coords
        out.write(f'{len(atoms)}\n\n')

        for atom, (x, y, z) in zip(atoms, coords):
            out.write(f'{atom} {x} {y} {z}\n')

        if file_open:
            out.close()

    def to_cgrtools(self, radius_multiplier=1.25, store_log=False, multiplicity: Optional[int] = None):
        parser = XYZRead.create_parser(radius_multiplier=radius_multiplier, store_log=store_log)
        # вызывает метод create_parser класса XYZRead. Создается функция, которая записывается в переменную parser
        matrix = []
        for a, c in zip(self.atoms, self.coords):
            atom_num = atom_map[a]
            x, y, z = c
            matrix.append((atom_num, x, y, z))
        charge = self.charge
        if multiplicity is None:
            multiplicity = self.multiplicity
        if multiplicity == 1:
            radical = 0
        elif multiplicity == 2:
            radical = 1
        else:
            radical = multiplicity

        return parser(matrix, charge, radical)
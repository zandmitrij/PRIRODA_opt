from sklearn.base import BaseEstimator, TransformerMixin
from .conformer import Conformer
from os.path import dirname, join
from multiprocessing.pool import ThreadPool
from subprocess import check_call
from tempfile import mkdtemp
from itertools import tee, dropwhile, takewhile, islice
from shutil import rmtree


atom_map = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12,
            'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23,
            'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33,
            'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44,
            'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
            'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65,
            'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76,
            'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
            'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
            'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107,
            'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Uub': 112, 'Uuq': 114}
reverse_atom_map = {str(v): k for k, v in atom_map.items()}


class PrirodaOptimizer(BaseEstimator, TransformerMixin):
    def __init__(self, relative=False, basis="L1", memory=200, tmp_dir="/tmp", tmp_ram=10, n_jobs=1, n_process=1,
                 steps=100):
        # будем считать гессиан по умолчанию False
        # basis - базисный набор,
        # relative - релятивисткий базис
        # memory  - количество оперативной памяти для расчетов 200мБ
        # tmp_ram - количество памяти для храния временных файлов расчета, если не является 0
        # n_jobs
        # n_process
        # steps
        self.relative = relative  # точно символ в символ сохранить, особенность работы с сайкитлерн
        self.basis = basis
        self.memory = memory
        self.tmp_dir = tmp_dir
        self.tmp_ram = tmp_ram
        self.n_jobs = n_jobs
        self.n_process = n_process
        self.steps = steps

    def fit(self, x, y=None):  # обучающие данные
        return self

    def transform(self, x, y=None):
        # подаем список молекул c начальным приближением, возвращаем список молекул с оптимизированной геометриией
        # храним атомы, коорд, заряд, мультиплетность, для гессиана True/False
        # x - описывает свойства, экземпляр класса Conformer
        # валидатор х, iterable: x - list, tuple
        if not isinstance(x, (list, tuple)):
            raise TypeError
        for i in x:
            if not isinstance(i, Conformer):
                raise TypeError

        jobs_1, jobs_2 = tee(self._prepare_calc(x))
        results = []
        with ThreadPool(self.n_jobs) as p:  # создали н_джобс тредов, число потоков
            for _, (dir_, _, log), conf in zip(p.imap(check_call,
                                                      (('priexec', '-np', str(self.n_process), 'pribin', inp, out) for
                                                       _, inp, out in jobs_1)), jobs_2, x):
                a, c, h, e = self._parse_output(open(log))
                results.append(Conformer(a, c, conf.charge, conf.multiplicity, h, e))
                rmtree(dir_)
        return results

    def _prepare_input(self, conf: Conformer, tmp_dir):
        basis = join(dirname(__file__), f'basis{int(self.relative)}')
        out = [f'$system memory={self.memory} disk={self.tmp_ram} path={tmp_dir} $end',
               f'$control task=optimize+hessian theory=dft four={int(self.relative)} basis={basis} $end',
               '$dft functional=pbe $end', f'$optimize steps={self.steps} $end',
               f'$molecule charge={conf.charge} mult={conf.multiplicity}',
               ' cartesian',
               f' set={self.basis}']
        for atom, (x, y, z) in zip(conf.atoms, conf.coords):
            out.append(f' {atom_map[atom]} {x:.4f} {y:.4f} {z:.4f}')
        out.append('$end')
        return '\n'.join(out)

    def _prepare_calc(self, confs):  # generator
        for conf in confs:
            tmp_dir = mkdtemp(dir=self.tmp_dir)
            task = self._prepare_input(conf, tmp_dir)
            inp = join(tmp_dir, 'input')
            out = join(tmp_dir, 'output')
            with open(inp, 'w') as f:
                f.write(task)
            yield tmp_dir, inp, out

    def _parse_output(self, log):
        atoms = []
        coords = []

        for line in takewhile(lambda x: not x.startswith('MOL>$end'),
                              islice(dropwhile(lambda x: not x.startswith('MOL> set'), log), 1, None)):
            _, atom, x, y, z = line.split()
            atoms.append(reverse_atom_map[atom])
            coords.append((float(x), float(y), float(z)))

        freq = next(islice(dropwhile(lambda x: not x.startswith(' | Mode |'), log), 2, None)).split()[3]
        energy = float(next(dropwhile(lambda x: not x.startswith(' E = '), log)).split()[2])
        gibbs = float(next(islice(dropwhile(lambda x: not x.startswith(' translational'), log), 3, None)).split()[-1])
        return tuple(atoms), tuple(coords), 'i' not in freq, energy * 627.51 + gibbs

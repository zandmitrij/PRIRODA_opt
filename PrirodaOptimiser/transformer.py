from sklearn.base import BaseEstimator, TransformerMixin
from .conformer import Conformer
from os.path import dirname, join
from multiprocessing.pool import ThreadPool
from subprocess import check_call
from tempfile import mkdtemp
from itertools import tee, dropwhile, takewhile, islice
from shutil import rmtree


atom_map = {'H': 1, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Br': 35, 'I': 53}
reverse_atom_map= {str(v) : k for k, v in atom_map.items()}


class Priroda_Optimizer(BaseEstimator, TransformerMixin):
    def __init__(self, relative=False, basis="L1", memory=200, tmp_dir="/tmp", tmp_ram=10,
                 n_jobs=1, n_process=1, steps=100):
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
        # tee - на вход принимает поток, потом записывает в файл, но в ауте тоже его показывает.
        # Взяли итератор и из него сделали 2 генератора. Если он бесконечный, то он зависнет
        with ThreadPool(self.n_jobs) as p:  # создали н_джобс тредов, число потоков
            results = []
            for _, (dir, _, log) in zip(p.imap(check_call,
                                     (('mpiexec', '-np', str(self.n_process), 'priroda', inp, out)
                                      for dir, inp, out in jobs_1)), jobs_2):
        # imap - итератор, check_call - запуск
        # imap запускает паралелльно процессы, когда первый процесс завершится, то он его вернет, переходит к сл.
        # сделали генератор задач (('mpiexec', '-np', str(self.n_process), 'priroda', inp, out)
        #                                       for dir, inp, out in jobs_1))
        # самое главное зазиповать джобс_2
                results.append(self._parse_output(log))  # parse_output - вовзвращает конформер
                rmtree(dir)  # удаляет всю директорию, несмотря на то, что в ней что-то есть

    def _prepare_input(self, conf: Conformer, tmp_dir):
        # генерирует файл для природы
        # cartesian - Декарт
        basis = join(dirname(__file__), f'basis{int:self.relative}')
        out = [f'$system memory={self.memory} disk={self.tmp_ram} path={tmp_dir} $end']
        out.append(f'$control task=optimize+hessian theory=dft four={int(self.relative)}'
                   f'basis={basis} $end')
        out.append('$dft functional=pbe $end')
        out.append(f'$optimize steps={self.steps} $end')
        out.append(f'$molecule charge={conf.charge} mult={conf.multiplicity}')
        out.append(' cartesian')
        out.append(f' set={self.basis}')

        for atom, (x, y, z) in zip(conf.atoms, conf.coords):
            out.append(f' {atom_map[atom]} {x:.4f} {y:.4f} {z:.4f}')
            out.append('$end')
            return '\n'.join(out)

    def _prepare_calc(self, confs):  # generator
        for conf in confs:
            tmp_dir = mkdtemp(dir=self.tmp_dir)  # создали путь до временной директории
            task = self._prepare_input(conf, tmp_dir)  # создали содержимое файла
            inp = join(tmp_dir, 'input')
            out = join(tmp_dir, 'output')
            with open(inp, 'w') as f:
                f.write(task)
            yield tmp_dir, inp, out

    def _parse_output(self, log):  # log - open file
        atoms, coords = [], []
        for line in takewhile(lambda x: x.startswith('MOL>  '),
                              dropwhile(lambda x: not x.startswith('MOL>  '), log)):
            # dropwhile - generator
            # takewhile - create new generator
            _, atom, x, y, z = line.split()
            atoms.append(reverse_atom_map[atom])
            coords.append((float(x), float(y), float(z)))
        next(dropwhile(lambda x: not x.startswith('| Mode'), log))
        next(log)
        hessian = 'i' not in next(log).split()[3]
        # нужен 4-ый столбец. Если в этом столбце есть буква i, то гессиан False
        energy = float(next(islice(dropwhile
                    (lambda x: not x.startswith('T = '), log), 6, 7)).split()[-1])  # попросили с 6 по 7 в генераторе
        # создать новый объект конформер, и вернуть его 
from setuptools import setup

setup(
    name='PrirodaOptimizer',
    version='0.1.0',
    packages=['PrirodaOptimizer'],
    license='CIMM.07-013',
    author="Dmitry x2, Zhenya x2, Daria, Polina, Il'nura, Roma, Anya",
    url='https://cimm.site/',
    data_files=[('bin', ['bin/pribin', 'bin/priexec'])],
    package_data={'PrirodaOptimizer': ['basis0', 'basis1']},
    install_requires=['scikit-learn'],
    zip_safe=False)

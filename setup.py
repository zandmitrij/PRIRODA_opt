from setuptools import setup

setup(
    name='PrirodaOptimizer',
    version='0.1.0',
    packages=['PrirodaOptimizer'],
    license='CIMM.07-013',
    author="Dmitry x2, Zhenya x2, Darya, Polina, Il'nura, Roma, Anya",
    url='https://cimm.site/',
    package_data={"PrirodaOptimizer": ['conformer', 'transtormer', 'data/*']},
    install_requires=['scikit-learn'],
    zip_safe=False)

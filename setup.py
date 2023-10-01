from setuptools import setup

setup(
    name='pysasa',
    version='0.0.4',
    description='A module for calculating the solvent accessible surface area of molecules',
    url='https://github.com/avanteijlingen/pySASA',
    author='Alexander van Teijlingen',
    author_email='a.vant@linuxmail.org',
    license='BSD 2-clause',
    packages=['pysasa'],
    install_requires=['ase',
                      'pandas',
                      'numpy',
                      ],

    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
    ],
)

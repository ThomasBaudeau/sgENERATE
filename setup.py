from setuptools import setup
from sgENERATE import __version__

setup(
    name='sgENERATE',
    version=__version__,
    packages=['sgENERATE'],
    package_dir={'sgENERATE': 'sgENERATE'},
    package_data={'sgENERATE':['recource']},
    url='',
    license='',
    author='Thomas.Baudeau',
    author_email='thomas.baudeau@univ-lille.fr',
    description='sgENERATE create read sample mimicking ARTIC ONT data and benchmark tool for and quantifies sgRNAs in SARS-CoV-2',
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'sgENERATE = sgENERATE.sgENERATE:main',
        ],
    },
    install_requires=['snakemake'],
)

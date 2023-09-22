from setuptools import setup
from periscope import __version__, _program

setup(
    name='sgENERATE',
    version=__version__,
    packages=['sgENERATE'],
    packages=find_packages(),
    package_dir={'periscope': 'periscope'},
    package_data={'periscope':['resources/*']},
    url='',
    license='',
    author='Thomas.Baudeau',
    author_email='thomas.baudeau@univ-lille.fr',
    description='sgENERATE create read sample mimicking ARTIC ONT data and benchmark tool for and quantifies sgRNAs in SARS-CoV-2',
    entry_points={
        'console_scripts': [
            'sgENERATE = sgENERATE.__main__:main',
        ],
    },
    install_requires=['snakemake']
)
    include_package_data=True,
)

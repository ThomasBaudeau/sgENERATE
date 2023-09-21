# sgENERATE

Snakemake pipeline to compare sgRNA finding tools 

*<h2>Authors </h2>* 

PhD Student, Lille University

Thomas BAUDEAU – thomas.baudeau@univ-lille.fr

<h2>Install </h2>*

**Installation**

Clone the repository:

<pre>
git clone https://github.com/ThomasBaudeau/BenchMapL <br> 
cd BenchMapL
</pre>


**Requirement**

* [Snakemake workflows](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) <br> 
* [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) <br> 


## How to start with BenchMapL :

The structue of the workflow is build as following:

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── workflow
    │   ├── benchmarks
    |   │   ├── benchresult.txt
    |   │   └── ...
    |   ├── data
    |   |   ├── sample
    |   |   │   ├── species_length_error-rates.fasta
    |   |   │   └── ...
    |   |   └── ref_species.fasta
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── helps_tool
    |   │   ├── tools_helpdoc.txt
    |   │   └── ...
    │   ├── plots
    |   │   ├── plot1.pdf
    |   │   └── plot2.pdf
    |   ├── resu
    |   │   ├── resu.bam
    |   │   └── ...
    |   └── Snakefile
    ├── workflow2
    |   ├── data
    |   |   ├── supl
    |   |   │   ├── model.model
    |   |   │   └── ...
    |   |   └── ref_species.fasta
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── result
    |   │   ├── nano
    |   │   └── pbsim
    |   └── Snakefile
    ├── config
        ├── config.yaml
        └── config2.yaml

The different mappeurs implemented are :

  * [Minimap2](https://github.com/lh3/minimap2)
  * [Graphmap2](https://github.com/lbcb-sci/graphmap2)
  * [Graphmap](https://github.com/isovic/graphmap)
  * [Blasr](https://github.com/PacificBiosciences/blasr)
  * [Winnowmap2](https://github.com/marbl/Winnowmap)
  * [MagicBlast](https://ncbi.github.io/magicblast/)
  * [lra](https://github.com/ChaissonLab/LRA)


 ### BenchMapL Usage:

 #### Workflow2 : Data Generation Part

  1. Add the species file of the reference in fasta format to the directory *data*. 
  > a read files of the species must be add to use Nanosim. 
  2. Rename the file with the following format : __ref__\_ __species-name__.fasta 
  > Nanosim read file : __read__\_ __species-name__.fasta
  3. Open the config2 files and add the *species-name* in the species fields.
  4. Modifies the __size__, __error_rate__ and __number__ fields to change to the desired length, error rate and coverage for each generated read file
  5. run with :
        * "__snakemake -R 1 all__"  : for Pbsim2 reads generation stage
        * "__snakemake -R 1 nano__" : for Nanosim reads generation stage
  
 #### Workflow : Mappeur Comparaison Parts

   1. Add the different datasets generated in the directory *data/sample* and the different reference files in the directory *data*.
   2. Open the config files and add the *species-name* in the species fields.
   3. Modifies the __size__ and __error_rate__ fields to change to the corresponding length, error rate and coverage of each generated read file
   4. run with :
        * "__snakemake -R 1 run-default__" : for comparing each tools in there default configuration
        * "__snakemake -R 1 run-params__" : for comparing each tools with multiple configuration
        * "__snakemake -R 1 run-perfect__"  : for comparing each tools with there best configuration


 ### BenchMapL Configuration:

#### Workflow : Adding a new configuration for a tool
1. Open the config files and add a new command name in the __name__ field of the choosen tool. The name can not contains * character
2. Add the command in the __command__ fields the name and the command position must be the same. 
> the space character must be replace by "#" and the __\___ by "§" 

#### Workflow : Change the best configuration for a tool
1. Open the config files and find the "*" character in __name__ field.
1. Delete "*" character and add it to the name of the choosen command.

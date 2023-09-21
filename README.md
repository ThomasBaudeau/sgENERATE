# sgENERATE

Snakemake pipeline to generate ONT read with sgRNA mimicking ARTIC protocol and to compare sgRNA finding tools 

<h2>Install </h2>*

**Installation**

Clone the repository:

<pre>
git clone https://github.com/ThomasBaudeau/sgENERATE <br> 
cd BenchMapL
</pre>


**Requirement**

* [Snakemake workflows](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) <br> 
* [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) <br> 


## How to start with BenchMapL :

The different Tool implemented are :

  * [Periscope](https://github.com/sheffield-bioinformatics-core/periscope)
  * [Periscope_multi](https://github.com/ThomasBaudeau/periscope_multifasta)



 ### BenchMapL Usage:

 #### Part1 : Data Generation Part

 #### Workflow : sgRNA Tool Comparaison Parts



 ### BenchMapL Configuration:

#### Workflow : Adding a new configuration for a tool
1. Open the config files and add a new command name in the __name__ field of the choosen tool. The name can not contains * character
2. Add the command in the __command__ fields the name and the command position must be the same. 

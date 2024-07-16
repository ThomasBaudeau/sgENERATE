# sgENERATE

A Snakemake pipeline to generate ONT reads with sgRNA mimicking the ARTIC protocol and to compare sgRNA finding tools.

More information on the [sgENERATE documentation](https://thomasbaudeau.github.io/sgENERATE/)


<h2>Installation</h2>*


**Requirement**

* [Snakemake workflows](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) <br> 
* [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) <br> 


**Installation Steps**

Clone the repository and set up the environment:

<pre>
git clone https://github.com/ThomasBaudeau/sgENERATE <br> 
cd sgENERATE
conda env create -f environment.yml
conda activate sgENERATE
pip install .
</pre>




## Getting Started with sgENERATE


The tools implemented are:

  * [Periscope](https://github.com/sheffield-bioinformatics-core/periscope)
  * [Periscope_multi](https://github.com/ThomasBaudeau/periscope_multifasta)



 ### sgENERATE Usage:
 <pre>
--coverage (default 5000)
--error mean error rate of the generated samples (default 0.95, i.e., 5% error)
--compare if true, don't run the benchmark and only generate the data (default false)
</pre>

# sgENERATE

A Snakemake pipeline to generate ONT reads with sgRNA mimicking the ARTIC protocol and to compare sgRNA finding tools.

More information on installation, running instructions, parameters etc is found on the [sgENERATE documentation](https://thomasbaudeau.github.io/sgENERATE/) website.


<h2>Installation</h2>


**Requirements**

* [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) <br> 
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



 ### sgENERATE Example usage:

Simply run `sgENERATE` from the command line to both simulate an ARTIC sgRNA dataset and execute a benchmark between currently available tools (periscope and periscope_multi).

To only simulate data, execute `sgENERATE --compare False`

 ### sgENERATE parameters

 Type `sgENERATE --help` for more information on available parameters.


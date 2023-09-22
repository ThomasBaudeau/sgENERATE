# sgENERATE

Snakemake pipeline to generate ONT read with sgRNA mimicking ARTIC protocol and to compare sgRNA finding tools 


<h2>Install </h2>*


**Requirement**

* [Snakemake workflows](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) <br> 
* [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) <br> 


**Installation**

Clone the repository:

<pre>
git clone https://github.com/ThomasBaudeau/sgENERATE <br> 
cd sgENERATE
pip install .
</pre>




## How to start with sgENERATE :

The different Tool implemented are :

  * [Periscope](https://github.com/sheffield-bioinformatics-core/periscope)
  * [Periscope_multi](https://github.com/ThomasBaudeau/periscope_multifasta)



 ### sgENERATE Usage:
 <pre>
--coverage (default 5000)
--error mean error rate of the generated samples (default 0.95 i.e 5% error)
--compare  if true don't run the benchmark and only make the data generation (default false) 
</pre>

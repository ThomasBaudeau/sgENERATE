from sgENERATE.recource.script.make_stats import main

try:            
    main(snakemake.input['a'],snakemake.input['b'],snakemake.input['peri'],snakemake.input['d'],snakemake.input['a2'],snakemake.input['peri2'],snakemake.input['d2'],snakemake.output[0],snakemake.input['nbread'],file3=snakemake.input['c']) 
except AttributeError:
    main(snakemake.input['a'],snakemake.input['b'],snakemake.input['peri'],snakemake.input['d'],snakemake.input['a2'],snakemake.input['peri2'],snakemake.input['d2'],snakemake.output[0],snakemake.input['nbread'])

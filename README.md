# How to 

The pipeline is designed to be executed both on a HPC cluster or a standalone server. To start, edit the file `config_main.yaml` adding the absolute PATH to the required resources. 

## Sources for the required files

TODO

### Using the same annotation in every files

## Run on Cluster HPC

The pipeline was developed and tested on a HPC cluster which uses the Slurm scheduler. If your HPC is using another scheduler, these options may not work as expected.

Inside the `cluster_config.json` add your **account_name** and desired **partition** where the job will be executed inside the `__default__` definition. 

```json
"__default__":
    {
        "account": "John Smith",
        "time": "01:30:00",
        "ncpus": 2,
        "mem": "8G",
        "partition": "somewhere",
        "out": "{logpath}/{rule}-{jobid}.out"
    }
```

This ensure that for each rule not specified inside the `cluster_config.json`, the Slurm scheduler will submit a job requesting the same walltime, ncpus and RAM. Indeed this is a generalistic setup which may results in allocating higher resources than the ones really needed, but having the pipeline blocked for walltime exaustion or not enough resources is annoying.

To launch the pipeline use 

```bash
snakemake -s Snakefile --cluster-config cluster_config.json -j <max_job_number> --use-conda --rerun-incomplete --cluster 'sbatch -p {cluster.partition} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out} -t {cluster.time}'
```

Where `<max_job_number>` is an integer representing the maximum number of jobs that could be submitted on your cluster. Even if you declare a number higher than the one permitted, it may be executed scaling automatically to the maximum possible. Note that to guarantee lower execution time, Mutect2 works on intervals, so a single job will be spawned for each interval, resulting in an high number of parallel jobs. So, a number higher than 100 is recommended. 




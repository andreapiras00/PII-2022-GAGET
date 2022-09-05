# GAGET - Genome Assembly Graph Evaluation Toolkit

This is a tool used to evaluate the alignment between an assembly graph, in the `.gfa` format, and its reference, in the `.fasta` or `.fa` format.

## Setup the environment

In order to execute the application:

1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Clone this repo by either downloading the `.zip` and extracting it, or using the `git clone` command.
3. Open a terminal, navigate to the project folder and create a new conda environment:
   ```bash
   cd /path/to/project/home/directory
   conda env create -f dependencies.yml
   ```
   This command should automatically install all required dependencies on the environment
4. Activate the newly created environment:
   ```bash
   conda activate gaget
   ```

## Excecute the software

To launch the application run the following command inside the environment:

```bash
python GAGET.py -g path/to/graph/ -r path/to/reference
```

Additional options can be found using the `--help` option:

```bash
python GAGET.py --help
```

or

```bash
python GAGET.py -h
```

### Example execution

```bash
python GAGET.py -g Test_Data/ebola_graph.gfa -r Test_Data/KM034562v1.fa
```

## System requirements

- Operating system: Linux, including Ubuntu, RedHat, CentOS 7+, and others.

**Known limitations**
This tool can currently only be excecuted on a `Linux` OS, due to a bug of `graph-tool` package with macOS.

<p align="center">
  <img width="2000" src="/figures/logo.png">
</p>

# ClusterBOSS

This is the README document for ClusterBOSS (Cluster Based On Sequence Similarity), a Python pipeline for clustering sequencing data from _in vitro_ selection experiments into families of sequence similarity. This tool clusters sequences using edit distance as a measure of sequence similarity. ClusterBOSS can be used to process nucleotides or amino acids sequencing data.

---------------------------------------

# Usage

The clustering tool is written in Python. It can be run using a Python interpreter, like the Command Line Interface (aka the Terminal) or any other specific software (e.g. PythonWin) and can be run using any version of Python 3. To run the script from the Terminal, type:

`python ClusterBOSS.py input_file d_cutoff n_min a_min c_min rec(y/n) keep_not_clustered(y/n)`

* input_file: name of input file (must include the full path to the directory where it's located).
* d_cutoff.txt: cutoff distance used to cluster.
* n_min: minimum number of sequences per peak.
* a_min: minimum abundance of sequences included in peaks.
* c_min.txt: minimum abundance of center sequences.
* rec(y/n): whether sequences should be used in more than one peak.
* keep_not_clustered(y/n): whether an output file with unclustered sequences should be generated.


# Environment setup

We recommend [using Anaconda to create a virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). Although this is not necessary, using a virtual environment can prevent version conflicts and undesired upgrades/downgrades on already existing packages. 

Once the virtual environment is active, Python and the other required dependencies can be installed there. The only dependency that is needed is [python-Levenshtein](https://pypi.org/project/python-Levenshtein/). This dependency can either be installed manually, or using the `requirements.txt` file.

To create a virtual environment and install all dependencies:

```
# create environment
conda create -n myenv python=3

# activate environment
source activate myenv

# install additional dependencies
conda install --file requirements.txt 
```

Alternatively (although not recommended), ClusterBOSS can be run from a local environment, where Python must be installed. In this case, we recommend using the Anaconda distribution of Python. See the [Anaconda documentation](https://docs.anaconda.com/anaconda/install/) for installation. 

# Input

The script requires one input files: sequencing reads. Sequencing reads are assumed to be 'galaxy-type', that is, to have 3 head lines: number of unique sequences, total number of molecules and an empty line.

# Output

The pipeline will generate an output directory, called `e+d_cutoff`, where the peak files will be stored.

# Reporting bugs

Please report any bugs to Celia Blanco (celiablanco@ucla.edu). 

When reporting bugs, please include the full output printed in the terminal when running the pipeline. 

# Citation

Evan Janzen, Yuning Shen, Ziwei Liu, Celia Blanco, Irene A. Chen. Error minimization and specificity could emerge in a genetic code as by-products of prebiotic evolution. *Submitted.*


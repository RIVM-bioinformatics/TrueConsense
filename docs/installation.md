# Installation instructions

TrueConsense is only available on Linux (or Linux-based) operating systems. MacOS may also work but is not tested. TrueConsense will *not work* on Windows.

TrueConsense will be made available for installation through Conda and Pip. However, this is currently not yet available.
We will update these docs when installation through Conda and/or pip is available.

## Prerequisistes

TrueConsense requires Python 3.7 or later to be installed on your system (or in an environment).

Other dependencies will be installed during the installation, your don't have to install them manually. These extra dependencies are as follows:

* pysam<0.16
* pandas>=1.2.3
* gffpandas>=1.2.0
* parmap>=1.5.2
* tqdm>=4.59.0
* biopython>=1.78

We strongly advise you to use a conda environment (or similar) to make sure there won't be any conflicts in package dependencies.

## Download and install from source

First start by cloning the repository and checkout out the latest released version of TrueConsense:
```bash
git clone https://github.com/RIVM-bioinformatics/TrueConsense.git; cd TrueConsense; git checkout tags/$(git tag --sort=committerdate | tail -1) >> /dev/null
```

You're now in the newly create "TrueConsense" directory.

!!! tip "Make a new Conda environment before continuing"
    If you have Conda installed on your system, please create and activate a new environment before continuing.

    Use the following command to create and activate a new Conda environment named "TrueConsense" based on the environment-recipe we provide in the github-repository

    ```bash
    conda env create -f env.yml; conda activate TrueConsense
    ```

    The "TrueConsense" conda-environment should now be active.

You can now install TrueConsense via the following command:
```bash
pip install .
```

TrueConsense should now be installed!
You can verify if installation was successful by typing `trueconsense --version` on the command-line, this should show the installed TrueConsense version.

## Pipeline/workflow inegration

You can easily integrate TrueConsense in your snakemake bioinformatics workflow if you use Conda environments in your workflow.

To do this, simply add the following structure to your conda-environment recipe, replace `{VERSION}` with the TrueConsense version you wish to use:

```yaml
dependencies:
    - pip
    - pip:
        - git+https://github.com/RIVM-bioinformatics/TrueConsense.git@{VERSION}
```

Conda will now install TrueConsense and its dependencies in the specified snakemake conda-environment.

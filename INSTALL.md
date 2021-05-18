NetSquid Snippet Installation
=============================

These are the general installation instructions for [NetSquid snippet](https://netsquid.org/snippets) packages.

Installation
------------

For Python to be able to find the NetSquid snippet package it needs to be installed as a package or added to the `PYTHONPATH` environment variable.

### Requirements

The calculation of the parameters of the satellite downlink channel depends on two external packages, *orekit*, for orbit determination, and *lowtran* for atmospheric absorption. The packages must be installed before the installation of the snippet.

#### Orekit

[Orekit](https://www.orekit.org/) is an open-source space dynamics library written in Java, with a [wrapper](https://gitlab.orekit.org/orekit-labs/python-wrapper) which allows to access its functionalities in Python.

The recommended way to install it is through [Anaconda](https://docs.continuum.io/anaconda/install/):

```shell
conda install -c conda-forge orekit
```

More information can be found in the [installation tutorial](https://gitlab.orekit.org/orekit-labs/python-wrapper/-/wikis/installation) of the orekit python-wrapper wiki.

#### Lowtran

[Lowtran](https://pypi.org/project/lowtran/) is a Python package providing direct access to the LOWTRAN7 atmospheric absorption extinction model developed by the Air Force Geophysics Laboratory (AFGL). Since the model is implemented in Fortran, it is necessary to install a suitable compiler before the package.

To install Gfortran:
* Linux: `apt install gfortran`
* Max: `brew install gcc`
* [Windows](https://www.scivision.dev/windows-gcc-gfortran-cmake-make-install/)

The Python lowtran package can be installed using pip. More information can be found in the [lowtran page](https://pypi.org/project/lowtran/).

### Install using pip

To install the package and its requirements using pip run the following command in the repository root directory:

```shell
make install
```

Note: To be able to install NetSquid and possibly other snippets on the netsquid server you need to provide your user name and password for the pypi server (*pypi.netsquid.org*); these match your forum credentials. You can store the user name and password in the environment variables NETSQUIDPYPI_USER and NETSQUIDPYPI_PWD, respectively, to prevent having to type them in manually during installation.

### Install without using pip

To install without pip run the following command in the repository root directory:

```shell
python3 setup.py install --user
```

### Install by adding to _PYTHONPATH_

Add this repository to your `PYTHONPATH` environment variable.
For example, in a bash shell do:

```shell
export PYTHONPATH=$PYTHONPATH:/path/to/this/repository
```

or add this line to your `$HOME/.bashrc`.

Running the tests
-----------------

To run all of the available unit tests including the linter, run the following command in the repository root directory:

```shell
make verify
```

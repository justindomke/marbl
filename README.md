marbl
=====

Marginal based learning of Conditional Random Field parameters.

# Overview

This implementation roughly corresponds to the algorithm described in the paper:
 * Justin Domke, [Learning Graphical Model Parameters with Approximate Marginal Inference](http://users.cecs.anu.edu.au/~jdomke/papers/2013pami.pdf), IEEE Transactions on Pattern Analysis, 2013.

# Getting started

* Make sure that you have a compiler that supports C++11.  Recent versions of g++ and clang++ are known to work.
* Install openMPIand/or openMP, if you want parallelism.  (Optional; see below)
* Download the [code](https://github.com/justindomke/marbl/archive/master.zip).
* Go to the main code directory, and edit the `make.sh` script in three ways
 * Change the `compiler` variable to the one for your system
 * Set the `use_openmp` variable to true/false, depending on if you want openMP support
* Run the `make.sh` script.  This will compile and install libLBFGS to a local directory, and then build the `infer_MRF`, `learn_CRF` and (if you have MPI installed) `learn_CRF_mpi` executables.
* Go through a few of the [examples](examples) to get some idea of what you are doing.

Note that Marbl has been compiled under Mac OS and Linux, but hasn’t been tested under Windows.  If you are able to compile under Windows, please send any information about how you did so.

# Features

* Marbl can handle quite graphs, with factors of any size, variables taking any number of values, and whatever entropy approximation you desire.
* Marbl can use openMP for parallelism on a single computer with multiple cores.
* Marbl uses [openMPI](http://www.open-mpi.org/) for parallelism, either on a single computer with multiple cores or on a cluster of machines.  This is trivial to use on a single machine to get your problem running.  (Essentially, you just install openMPI, and then compile using the provided scripts.)  Since almost all clusters provide MPI, you can also run Marbl on clusters of computers.  We have observed essentially linear speedups using up to several hundred cores.
* Marbl is a command-line tool, with inputs in very simple specified text formats.  Thus, you can use whatever (presumably high-level) language to create your problem, by writing a simple routine to produce the data and graph specifications in the appropriate format.

# Requirements

Strictly speaking, Marbl requires no external libraries.  However, learning can be parallelized, using either openMP or openMPI.  These are both optional.

### openMP

[openMP](http://openmp.org/wp/) is a platform for shared-memory computing.  (i.e. using multiple threads on a single machine).  The executable `learn_CRF` can use openMP if your compiler supports it.  Generally, it would be easier to just install a compiler that supports openMP than trying to install openMP separately.

[Most recent compilers](http://openmp.org/wp/openmp-compilers/) support openMP out of the box.  The major exception is recent versions of xcode/clang on Mac OS.  (i.e. the default gcc compiler).

### openMPI

Despite the similar name, and despite similar goals of parallel computing, openMPI is completely disjoint from openMP.  [openMPI](http://www.open-mpi.org/) is a tool for parallel computing. While openMPI is typically used on clusters, it can manage communication between processes on a single computer.  Thus, Marbl’s `learn_CRF_mpi` uses this single mechanism for parallelism, both on single machines and clusters.
  * On linux or Mac OS, openMPI is trivially installed using your favorite package management system.  For example, using [homebrew](http://brew.sh/) on Mac OS, it can be installed with `brew install openmpi`.
  * On Windows, pre-compiled binaries for Cygwin are available [here](http://www.open-mpi.org/software/ompi/v1.8/). 


# Comparison with JGML

Marbl is similar to the [JGMT](http://users.cecs.anu.edu.au/~jdomke/JGMT/) toolbox (implementing essentially the same algorithm), with the following differences:

### Language
* JGMT is implemented in a mix of Matlab and C++, and must be used through Matlab.
* Marbl is a command-line tool written in pure C++.  Interaction is via the command-line and text files specifying the data and structure of the graph.
* However, one would typically create the inputs to Marbl using a high-level language, and an interface doing this is provided in Matlab.

### Limitations
* JGMT can only handle pairwise graphs.  Marbl can handle arbitrarily large factors.
* JGMT assumes that all variables have the same possible number of values.  Marbl does not assume this.
* JGMT is reasonably efficient.  Marbl is still more reasonable.

### Parallelism
* Marbl uses openMPI for parallelism.  This can be used either on a single machine, or on a cluster of computers.
* JGMT is parallel to some degree, using Matlab’s parfor mechanism.  However, this requires the (expensive) parallel computing toolbox to use up to 12 cores, and the (yet more expensive) distributed computing toolbox to use more than that.

# Acknowledgements

Marbl includes the code for two other packages:

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), a template library for linear algebra.
* [libLBFG](http://www.chokkan.org/software/liblbfgs/) is a C library for the L-BFGS optimization algorithm.

Note that while Marbl and libLBFGS are under the [MIT license](http://opensource.org/licenses/MIT) (under which you are, roughly speaking, allowed to do whatever you want) Eigen is under the [Mozilla public license](http://www.mozilla.org/MPL/2.0), which places some (very minor) restrictions on how it (and thus Marbl) can be used.
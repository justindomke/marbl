marbl
=====

Marginal based learning of Conditional Random Field parameters.

# Overview

This implementation roughly corresponds to the algorithm described in the paper:
 * Justin Domke, [Learning Graphical Model Parameters with Approximate Marginal Inference](http://users.cecs.anu.edu.au/~jdomke/papers/2013pami.pdf), IEEE Transactions on Pattern Analysis, 2013.

# Getting started



# Requirements

Compiling Marbl requires a C++ compiler that supports C++11, and two easily-installed external projects.

For example, using [homebrew](http://brew.sh/) on Mac OS, these can be installed as `brew install liblbfgs` and `brew install openmpi`.

* [openMPI](http://www.open-mpi.org/) is a tool for parallel computing. While openMPI is typically used on clusters, it can manage communication between processes on a single computer.  Thus, Marbl uses this single mechanism for parallelism, both on single machines and clusters.
  * On linux or Mac OS, openMPI is trivially installed using your favorite package management system, as discussed above.
  * On Windows, pre-compiled binaries for Cygwin are available [here](http://www.open-mpi.org/software/ompi/v1.8/). 

* [libLBFG](http://www.chokkan.org/software/liblbfgs/) is a C library for the L-BFGS optimization algorithm.
 * On Linux or Mac OS, liblbfgs is trivially installed using your favorite package management system.

# Features

* Marbl can handle quite arbitrary problems, with factors of any size, and variables taking any number of values.
* Marbl uses [openMPI](http://www.open-mpi.org/) for parallelism.  This is trivial to use on a single machine to get your problem running.  (Essentially, you just install openMPI, and then compile using the provided scripts.)  Since almost all clusters provide MPI, you can also run Marbl on clusters of computers.  We have observed essentially linear speedups using up to several hundred cores.
* Marbl is a command-line tool, with inputs in very simple specified text formats.  Thus, you can use whatever (presumably high-level) language to create your problem, by writing a simple routine to produce the data and graph specifications in the appropriate format.

# Comparison with JGML

Marbl is similar to the [JGMT](http://users.cecs.anu.edu.au/~jdomke/JGMT/) toolbox, with the following differences:

* Language
 * JGMT is implemented in a mix of Matlab and C++, and must be used through Matlab.
 * Marbl is a command-line tool written in pure C++.  Interaction is via the command-line and text files specifying the data and structure of the graph.
 * However, one would typically create the inputs to Marbl using a high-level language, and some examples of this in Matlab are provided.

* Limitations
 * JGMT can only handle pairwise graphs.  Marbl can handle arbitrarily large factors.
 * JGMT assumes that all variables have the same possible number of values.  Marbl does not assume this.
 * JGMT is reasonably efficient.  Marbl is better still.

* Parallelism
 * Marbl uses openMPI for parallelism.  This can be used either on a single machine, or on a cluster of computers.
 * JGMT is parallel to some degree, using Matlab’s parfor mechanism.  However, this requires the (expensive) parallel computing toolbox to use up to 12 cores, and the (yet more expensive) distributed computing toolbox to use more than that.

* Acknowledgements

Marbl includes the code for Eigen.

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), a template library for linear algebra.

Note that Eigen is under the [Mozilla public license](http://www.mozilla.org/MPL/2.0), which places some (very minor) restrictions on how it (and thus Marbl) can be used.
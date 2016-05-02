# An Implementation of Multi-Input Functional Encryption #

## Prerequisites ##

We rely on

 * [FLINT](http://flintlib.org/) for pretty much everything
 * [DGS](https://bitbucket.org/malb/dgs) for sampling from discrete Gausisans
   over the Integers

We ship these libraries as submodules.

## Installation ##

    git clone --recursive https://github.com/GaloisInc/gghlite-flint.git
    cd gghlite-flint
    mkdir dgs/m4
    autoreconf -i
    export CFLAGS="-O3"
    ./configure
    make -j

## Usage ##

### GGH Params ###

Prints parameter choices to stdout.

### GGH Instance ###

Instantiates a GGHLite instance

### NIKE ###

A non-interactive key exchange (but we skip running the hash function in the end). For example,
calling

    ./applications/nike -l 80 -N 7

will run 7-partite NIKE with security level λ=80.

## Subdirectories ##

 * **applications** NIKE, parameter estimation, benchmarks
 * **dgs** discrete Gaussians over the integers
 * **dgsl** discrete Gaussians over arbitrary regular lattices and ideal lattices modulo $x^{2^k}+1$
 * **flint** a current snapshot of the FLINT library
 * **gghlite** the actual GGHLite implementation
 * **oz** efficient functions for computing in $\\mathbb{Z}[x]/(x^{2^k}+1)$ and $\\mathbb{Z}_q[x]/(x^{2^k}+1)$
 * **tests** correctness tests


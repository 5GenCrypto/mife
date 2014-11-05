# An Implementation of GGHLite #

## Prerequisits ##

We rely on

 * [FLINT](http://flintlib.org/) for pretty much everything
 * [DGS](https://bitbucket.org/malb/dgs) for sampling from discrete Gausisans
   over the Integers

We ship these libraries as submodules.

## Installation ##

    git clone --recursive https://bitbucket.org/malb/gghlite-flint
    cd gghlite-flint
    mkdir m4
    mkdir dgs/m4
    autoreconf -i
    export CFLAGS="-O3"
    ./configure
    make -j2

## Usage ##

### GGH Params ###

Prints parameter choices to stdout.

### GGH Instance ###

Instantiates a GGHLite instance

### NIKE ###

A non-interactive key exchange (but we skip running the hash function in the end). For example,
calling

    ./applications/nike -l 80 -N 7 - f

will run 7-partite NIKE with security level λ=80. The parameter -f specifies that we are not
checking if our ideal is a prime ideal.

## Subdirectories ##

 * **applications** NIKE, parameter estimation, benchmarks
 * **dgs** discrete Gaussians over the integers
 * **dgsl** discrete Gaussians over arbitrary regular lattices and ideal lattices modulo $x^{2^k}+1$
 * **flint** a current snapshot of the FLINT library
 * **gghlite** the actual GGHLite implementation
 * **oz** efficient functions for computing in $Z[x]/(x^{2^k}+1)$ and $Z_q[x]/(x^{2^k}+1)$
 * **tests** correctness tests

## Bugs ##

This code has bugs, help us fix them.


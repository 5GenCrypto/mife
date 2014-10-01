# GGHLite Implementation #

## Prerequisits ##

We rely on

 * (FLINT)[http://flintlib.org/] for pretty much everything
 * (DGS)[https://bitbucket.org/malb/dgs] for sampling from discrete Gausisans
   over the Integers

We ship these libraries as submodules.

## Installation ##

    git clone --recursive https://bitbucket.org/malb/ggh-flint
    cd ggh-flint
    mkdir m4
    autoreconf -i
    ./configure
    make -j2

## Usage ##

### GGH Params ###

Prints parameter choices to stdout.

### GGH Instance ###

Instantiates a GGHLite instance

### NIKE ###

Most steps of a non-interactive key exchange

## Bugs ##

This code has bugs, don't use it yet.


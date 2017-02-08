Gravitational Lensing with Image Mapping
========================================

This repository contains a number of tools to extract information about 
gravitational lenses from the mapping between multiple images.


Contents
--------

The following tools are available.

-   [`ptmap`](#ptmap) - compute transformation matrices and lens quantities by
    mapping points between multiple images


Manual
------

### ptmap

    usage: ptmap [-vq] [-o OUTFILE] [-m MATFILE] FILE

The `ptmap` tool reads a list of observed points from `FILE`, where each line
corresponds to one multiple image, and the first image is the reference image.
The tool then finds affine transformations that cause the least total distance
squared between the transformed points of the reference image and the observed
points in the other multiple images.

The input file format is

    pt1; pt2; pt3 [ ; pt4 ... ]

where at least three points for three images must be given. The format for each
point is

    x, y [ , dx [ , dy [ , rho ] ] ]

where at least the position `x, y` must be given. The uncertainty of the point
is given by the optional covariance `dx, dy, rho`, where `dx = 1px`, `dy = dx`
and `rho = 0` is assumed by default. Uncertainties for the reference image are
ignored. An example definition is available in the `example/points.txt` file.
The [`reg2pts`](#re2pts) tool can be used to generate point definitions and
uncertainties from SAOImage DS9 region files.

From these transformations, the lens quantities `f` and `g` are computed and
printed in the form `f g1 g2`, where each line corresponds to the images as
given in the input file.

If the `-o` flag is given, results for `f` and `g` are written to `OUTFILE`
using the same format.

If the `-m` flag is given, the transformation matrices are written to `MATFILE`
in the form `T_11 T_12 T_21 T_22`, where `T_ij` are the matrix entries. The
total number of lines will be one less than the number of images, because the
transformations are relative to the reference image.

The `-v` and `-q` flags can be used to make the output more verbose and quiet,
respectively.

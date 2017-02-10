Gravitational Lensing with Image Mapping
========================================

This repository contains a number of tools to extract information about 
gravitational lenses from the mapping between multiple images.


Contents
--------

The following tools are available.

-   [**`immap`**](#immap) -
    map reference image using relative magnification matrices
-   [**`lens2mat`**](#lens2mat) -
    convert from lens quantities to relative magnification matrices
-   [**`mat2lens`**](#mat2lens) -
    convert from relative magnification matrices to lens quantities
-   [**`ptmatch`**](#ptmatch) -
    relative magnification matrices and lens quantities from point matching
-   [**`reg2pts`**](#reg2pts) -
    generate point definitions from SAOImage DS9 region files
-   [**`regcrop`**](#regcrop) -
    crop FITS files using SAOImage DS9 region files
-   [**`srcim`**](#srcim) -
    map observed images to source plane using lens quantities


Manual
------

### immap

map reference image using relative magnification matrices

    usage: immap MATFILE IMFITS OUTFITS

The `immap` tool reads a list of relative magnification matrices from `MATFILE`
and applies the transformation from each row to the image from `IMFITS`. The
mapped images are saved as a multi-extension FITS file in `OUTFITS`.

The matrix file must contain rows of the form `T_11 T_12 T_21 T_22`, where the
`T_ij` are the matrix entries. Each row corresponds to the magnification matrix
of one additional multiple image relative to a reference image. The output will
contain the reference image as the primary HDU (index 0), and one image-mapped
version of the reference image for each row of the matrix file, in the same
order. Hence the output file will contain one more multiple image than there
are given matrices.


### lens2mat

convert from lens quantities to relative magnification matrices

    usage: lens2mat [-q] [-o OUTFILE] FILE

The `lens2mat` tool reads a list of convergence ratios and shears from `FILE`
and computes the corresponding relative magnification matrices.

The input file must contain rows of the form `f g1 g2`, where each row
corresponds to one multiple image. The first row is the reference image for
which the relative magnification matrices are computed. All values of `f` are
normalised by the value given for the reference image, which can differ from
unity.

The relative magnification matrices are printed as `T_11 T_12 T_21 T_22`, where
`T_ij` are the matrix entries. The total number of rows will be one less than
the number of images, because the transformations are relative to the reference
image.

If the `-o` flag is given, the results are written to `OUTFILE` using the same
format.

The `-q` flag can be used to suppress output.


### mat2lens

convert from relative magnification matrices to lens quantities

    usage: mat2lens [-q] [-o OUTFILE] FILE

The `mat2lens` tool reads a list of relative magnification matrices from `FILE`
and computes the corresponding convergence ratios `f` and reduced shears `g`.

The input file must contain rows of the form `T_11 T_12 T_21 T_22`, where the
`T_ij` are the matrix entries. Each row corresponds to the magnification matrix
of one additional multiple image relative to a reference image. Hence there is
one less row in the input file than there are observed multiple images.

The convergence ratios `f` and reduced shears `g` are computed and printed in
the form `f g1 g2`, where the first row is for the reference image (with unit
convergence ratio) and each subsequent row corresponds to the additional images
as given in the input file.

If the `-o` flag is given, the results are written to `OUTFILE` using the same
format.

The `-q` flag can be used to suppress output.


### ptmatch

relative magnification matrices and lens quantities from point matching

    usage: ptmatch [-vq] [-o OUTFILE] [-m MATFILE] FILE

The `ptmatch` tool reads a list of observed points from `FILE`, where each row
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
ignored. An example point definition file [is available](example/points.txt).
The [`reg2pts`](#re2pts) tool can be used to generate point definitions and
uncertainties from SAOImage DS9 region files.

From these transformations, the convergence ratios `f` and reduced shears `g`
are computed and printed in the form `f g1 g2`, where each row corresponds to
the images as given in the input file.

If the `-o` flag is given, results for `f` and `g` are written to `OUTFILE`
using the same format.

If the `-m` flag is given, the relative magnification matrices are written to
`MATFILE` in the form `T_11 T_12 T_21 T_22`, where `T_ij` are the matrix
entries. The total number of rows will be one less than the number of images,
because the transformations are relative to the reference image.

The `-v` and `-q` flags can be used to make the output more verbose and quiet,
respectively.


### srcim

map observed images to source plane using lens quantities

    usage: srcim OUTFITS INFILE IMFITS0 [IMFITS1 ...]

The `srcim` tool reads a list of convergence ratios and shears from `INFILE`
and applies the magnification matrix from each row to the images `IMFITSn`. The
number of images supplied must match the number of rows in the input file, but
individual rows can skipped by giving `--` as the image name. The source images
are saved as a multi-extension FITS file in `OUTFITS`, where each extension is
called `SRCn`, where `n` corresponds to the index of the input image.

The input file must contain rows of the form `f g1 g2`, where each row
corresponds to one multiple image. The magnification matrices are normalised in
a way that the transformation of the first image has unit determinant. The
first value for `f` can differ from unity.

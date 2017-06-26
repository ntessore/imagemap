#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, NullLocator
from matplotlib.cm import get_cmap
from scipy.ndimage import gaussian_filter

# compute quantiles from CDF
def quantile(q, x, w):
    ord = np.argsort(x)
    cdf = np.cumsum(w[ord])
    cdf = (cdf - cdf[0])/(cdf[-1] - cdf[0])
    return np.interp(q, cdf, x[ord])

# parse arguments
parser = argparse.ArgumentParser(description='distribution plot from samples')
parser.add_argument('-a', action='store_true', help='plot anchor points')
parser.add_argument('-f', action='store_true', help='plot convergence ratios')
parser.add_argument('-g', action='store_true', help='plot reduced shears')
parser.add_argument('-x', action='store_true', help='print parameter values')
parser.add_argument('-b', metavar='BINS', type=int, default=50,
                    help='number of bins')
parser.add_argument('-s', metavar='SMOOTH', type=int, default=2.5,
                    help='smoothing scale in bins')
parser.add_argument('-t', metavar='TRUFILE', type=argparse.FileType('r'),
                    help='file containing truth values')
parser.add_argument('samfile', metavar='SAMFILE', type=argparse.FileType('r'),
                    help='file containing samples')
parser.add_argument('outfile', metavar='OUTFILE', type=str,
                    help='output file for plot')
args = parser.parse_args()

# make sure that at least something is shown
if not (args.a or args.f or args.g):
    parser.error('nothing to plot, use one or more of -afg')

# load samples
samples = np.loadtxt(args.samfile)

# number of samples
nsam = len(samples)

# number of parameter columns
npar = len(samples[0]) - 2

# number of images
nimg = int(npar/5);

# effective number of samples
neff = int(np.sum(samples[:,0])**2/np.sum(samples[:,0]**2))

# output some information
print('using {:} samples from {:} images'.format(nsam, nimg))
print('effective number of samples: {:}'.format(neff))

# read truth values if given
if args.t:
    truth = np.loadtxt(args.t)

# quantiles
ql = [ 0.005, 0.025, 0.160, 0.840, 0.975, 0.995 ]

# fill colours
cm = get_cmap('bone_r')
ec = 'w'
zc = cm(0.1)
qc = cm([0.2, 0.4, 0.6, 0.4, 0.2])

# first image in plot
fimg = 0 if args.g else 1

# number of rows in plot
nrow = nimg - fimg

# number of columns in plot
ncol = 0
if args.a: ncol += 2
if args.f: ncol += 1
if args.g: ncol += 2

# reshape truth values for plot
if args.t:
    truth = truth.reshape(nrow, ncol)

# size of figure
xdim = ncol*2.0
ydim = nrow*2.4

# create figure
fig = plt.figure(figsize=(xdim, ydim))

# create histograms, skip zero image if not plotting g
for i, j in zip(range(fimg, nimg), range(0, nrow)):
    cols = []
    labs = []
    
    if args.a:
        if i > 0:
            cols += [ 5*i + 2, 5*i + 3 ]
            labs += '$x_{0:}$ $y_{0:}$'.format(i).split()
        else:
            cols += [ None, None ]
            labs += [ None, None ]
    if args.f:
        if i > 0:
            cols += [ 5*i + 4 ]
            labs += '$f_{0:}$'.format(i).split()
        else:
            cols += [ None ]
            labs += [ None ]
    if args.g:
        cols += [ 5*i + 5, 5*i + 6 ]
        labs += '$g_{{{0:},1}}$ $g_{{{0:},2}}$'.format(i).split()
    
    for k in range(ncol):
        c = cols[k]
        l = labs[k]
        
        if not c:
            continue;
        
        ax = fig.add_subplot(nrow, ncol, j*ncol + k + 1)
        
        wht = samples[:,0]
        col = samples[:,c]
        
        m = quantile(0.5, col, wht)
        s = quantile([0.16, 0.84], col, wht) - m
        r = m + 5*s
        
        n, b = np.histogram(col, bins=args.b, weights=wht, range=r)
        
        if args.s:
            x = 0.5*(b[:-1] + b[1:])
            y = gaussian_filter(n, args.s)
        else:
            x = np.array(list(zip(b[:-1], b[1:]))).flatten()
            y = np.array(list(zip(n, n))).flatten()
        
        y /= np.max(y)
        
        ax.fill_between(x, y, facecolor=zc, edgecolor=ec, lw=0.1)
        
        if args.s:
            x2 = np.linspace(x[0], x[-1], 10*args.b)
            y2 = np.interp(x2, x, y)
            q = quantile(ql, col, wht)
            for a, b, c in zip(q[:-1], q[1:], qc):
                ax.fill_between(x2, y2, where=((a <= x2) & (x2 < b)),
                                facecolor=c, edgecolor=ec, lw=0.1)
        
        ax.plot(x, y, color='k', lw=0.5)
        
        if args.t:
            ax.axvline(truth[j,k], color='k', lw=1.5, alpha=0.5, ls='dashed')
        
        if args.x:
            l += ' = ${:.3f}_{{{:+.3f}}}^{{{:+.3f}}}$'.format(m, s[0], s[1])
        
        ax.set_xlabel(l)
        ax.set_xlim(r)
        ax.set_ylim(-0.05, 1.05)
        ax.set_yticks([], [])
        ax.xaxis.label.set_fontsize(16)
        ax.xaxis.set_label_coords(0.5, -0.2)
        try:
            ax.xaxis.set_major_locator(MaxNLocator('auto', prune='both'))
        except ValueError:
            ax.xaxis.set_major_locator(MaxNLocator(4, prune='both'))

# fix the layout
plt.tight_layout(h_pad=1.0, w_pad=0.0)

# save the plot
fig.savefig(args.outfile, bbox_inches='tight', figsize=(xdim, ydim))

import numpy as np
import argparse
import os.path
import matplotlib.pyplot as plt
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits


# Skype message
# [28/03/2013 22:48:33] Family: Jane Arthur:
# /fs/nas11/other0/jane/results/30032012_y/FITS/
# 10042012_y_0030vc-Halpha.fits
# 10042012_y_0030vlist-Halpha.fits
# 10042012_y_0030m0-Halpha.fits

def main(cubename, slice_mode, display,  iwindow, vlimits):

    # Agnostic as to whether .fits extensiom given or not
    if cubename.endswith(".fits"):
        prefix, _ = os.path.splitext(cubename)
    else:
        prefix = cubename
        cubename = prefix + ".fits"

    hdu = pyfits.open(cubename)[0]
    nv, ny, nx = hdu.data.shape

    fig = plt.figure(figsize=(6, 6))
    ax = plt.subplot(111)

    xpts = np.linspace(-2.0, 2.0, nx, endpoint=False)
    ypts = np.linspace(-2.0, 2.0, ny, endpoint=False)
    vpts = np.linspace(vlimits[0], vlimits[1], nv, endpoint=False)
    if slice_mode == "isovel":
        if iwindow is None:
            iwindow = [nv/2, nv/2]
        k1, k2 = iwindow
        imslice = hdu.data[k1:k2+1, :, :].sum(axis=0)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xlim(-2.0, 2.0)
        ax.set_ylim(-2.0, 2.0)
        A, B = np.meshgrid(xpts, ypts)
        ax.set_title("v = {:.1f} to {:.1f}".format(
            vpts[k1], vpts[k2]))
    elif slice_mode.endswith("slit"):
        if slice_mode.startswith("x"):
            if iwindow is None:
                iwindow = [ny/2, ny/2]
            j1, j2 = iwindow
            imslice = hdu.data[:, j1:j2+1, :].sum(axis=1).T
            ax.set_xlabel("v")
            ax.set_xlim(*vlimits)
            ax.set_ylabel("x")
            ax.set_ylim(-2.0, 2.0)
            A, B = np.meshgrid(vpts, xpts)
            ax.set_title("y = {:.1f} to {:.1f}".format(
                ypts[j1], ypts[j2]))
        else:
            if iwindow is None:
                iwindow = [nx/2, nx/2]
            i1, i2 = iwindow
            imslice = hdu.data[:, :, i1:i2+1].sum(axis=2).T
            ax.set_xlabel("v")
            ax.set_xlim(*vlimits)
            ax.set_ylabel("y")
            ax.set_ylim(-2.0, 2.0)
            A, B = np.meshgrid(vpts, ypts)
            ax.set_title("x = {:.1f} to {:.1f}".format(
                xpts[i1], xpts[i2]))
    else:
        raise NotImplementedError("Unknown mode: " + slice_mode)

    print "Minimum/Maximum: ", imslice.max(), imslice.min()
    if display == "contour":
        ax.contour(A, B, imslice)
    elif display == "grayscale":
        ax.imshow(imslice, origin="low", cmap=plt.cm.gray_r,
                  interpolation="nearest")
    ax.grid()
    fig.savefig("{}-{}-{}-{}.pdf".format(
        prefix, slice_mode, iwindow[0], iwindow[1]))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Plot isovel and PV images from
        PPV emission cubes"""
    )
    parser.add_argument(
        "cubename", type=str,
        help="""Name of FITS file containing PPV emission cube"""
    )

    parser.add_argument(
        "--slice-mode", type=str,
        choices=("x-slit", "y-slit", "isovel"),
        default="isovel",
        help="""Mode of operation - which way to slice"""
    )

    parser.add_argument(
        "--display", type=str,
        choices=("contour", "grayscale"),
        default="contour",
        help="""How to display the image"""
    )

    parser.add_argument(
        "--iwindow", type=int, nargs=2, default=None,
        help="""Range of positions or velocities to sum over
        (in pixel units)"""
    )

    parser.add_argument(
        "--vlimits", type=float, nargs=2,
        default=[-0.78107E+02, 0.73896E+02],
        help="""Minimum and maximum velocities in cube"""
    )

    cmd_args = parser.parse_args()
    main(**vars(cmd_args))


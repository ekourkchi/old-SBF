import os, sys, string, math, uuid, time, json
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import scipy.linalg as sla
import pandas as pd
import numpy as np
from matplotlib import patches
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from shapely.geometry.polygon import LinearRing
from astropy.io import fits
from astropy import wcs
from matplotlib import cm
from matplotlib.colors import LogNorm
from datetime import datetime
import pickle
import ipywidgets as widgets
from os.path import exists
import matplotlib.gridspec as gridspec
import numpy
from astropy.table import Table
import requests
from PIL import Image
from io import BytesIO
import pylab, os, sys

import copy

from .utils import *
from .sbfTools import *

##############################################################


def eval_plots(config):

    root = config["DIR"]
    psf = config["PSF"]
    prf = config["PRF"]
    rsd = config["RSD"]
    ptm = config["PTM"]

    cwd = os.getcwd()
    os.chdir(root)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(15, 5))

    tv(psf, ax=ax1, options="sqrt")
    ax1.set_title("PSF")
    tv(prf, ax=ax2, options="sqrt")
    ax2.set_title("Profile")
    tv(rsd, ax=ax3)
    ax3.set_title("Residual")
    tv(ptm, ax=ax4)
    ax4.set_title("Point Source Mask")

    os.chdir(cwd)

    return ax1, ax2, ax3, ax4


##############################################################


def get_sbf_signal(
    DIR,
    PSF,
    PRF,
    RSD,
    PTM,
    power_file="power.c0",
    tv_file="tv.jpg",
    K_fit_range=[25, 60],
    psf_order: int = 4,
    psf_k=[0, 15],
    radius=[32, 64],
    angle=[0, 360],
    center=[256, 256],
    nrows=512,
    ncols=512,
    sky_adjustment=0,
):

    KS0 = K_fit_range[0]
    KS1 = K_fit_range[1]
    psfk0 = psf_k[0]
    psfk1 = psf_k[1]
    a0 = angle[0]
    a1 = angle[1]
    r0 = radius[0]
    r1 = radius[1]
    X0 = center[0]
    Y0 = center[1]

    cwd = os.getcwd()
    os.chdir(DIR)

    monsta_script = (
        """

        rd 1 '"""
        + str(PSF)
        + """'                   ! reading psf from library
        rd 2 '"""
        + str(PRF)
        + """'                   ! eliprof final model
        rd 3 '"""
        + str(RSD)
        + """'                   ! galaxy subtracted residual
        rd 5 '"""
        + str(PTM)
        + """' bitmap            ! hybrid point source mask  (dophot+sextractor)
        ac 3 """
        + str(sky_adjustment)
        + """
        mi 3 5
        cop 4 1
        open 6 nr="""
        + str(nrows)
        + """ nc="""
        + str(ncols)
        + """
        fluc 6 5 mask x0="""
        + str(X0)
        + """ y0="""
        + str(Y0)
        + """ r0="""
        + str(r0)
        + """ r1="""
        + str(r1)
        + """ a0="""
        + str(a0)
        + """ a1="""
        + str(a1)
        + """ ! c0     taking the donut shape and putting the mask
        cop 7 6
        mi 7 3
        !tv 7 
        fluc 6 2 window         ! taking the fft of the mask multiply by the sqrt of the galaxy profile 
        ! fluc 6 4 expect plot ... 
        fluc 6 4 expect order="""
        + str(psf_order)
        + """ k0="""
        + str(psfk0)
        + """ k1="""
        + str(psfk1)
        + """   ! expectation power spectrum for psf  0->20
        fft 8 7
        power 4 8
        fluc 6 4 fit ks0="""
        + str(KS0)
        + """ ks1="""
        + str(KS1)
        + """ plot   ! 6:expectation opower spectrum   4:the power spectrum of data   ks: the power range ->>
        print fluc file="""
        + power_file + """
        """
    )

    if tv_file is not None:
        monsta_script += """     
        tv 7 JPEG="""+ tv_file + """
        """

    run_monsta(monsta_script, "monsta.pro", "monsta.log")

    with open(power_file, "r") as f:
        lines = f.readlines()

    P0 = np.float(lines[4].rsplit(" ")[2])
    dP0 = np.float(lines[4].rsplit(" ")[5])

    os.chdir(cwd)

    if tv_file is None:
        tv_file = ""

    return (P0, dP0), os.path.join(DIR, power_file), os.path.join(DIR, tv_file)


##############################################################
def extract_sbf_result(fname):

    with open(fname, "r") as f:
        lines = f.readlines()

    P0_h = np.float(lines[4].rsplit(" ")[2])

    K = []
    P0 = []
    D = []
    Power = []
    X = []

    for i in range(7, len(lines)):

        line = lines[i]

        l = [w.strip() for w in line.split(" ") if w != ""]

        #     print(l)

        Power_ = l[7]
        P0_ = l[8]

        if len(l[7]) > 12:
            Power_ = Power_[:8]
            P0_ = l[7][8:]

        if len(P0_) > 12:
            P0_ = P0_[:8]

        try:
            K.append(np.float(l[0]))
            Power.append(np.float(Power_))
            P0.append(np.float(P0_))
            D.append(np.float(l[4]))
            X.append(np.float(l[1]))
        except Exception as e:
            print(f"An exception occurred: {str(e)}")
            print("line #{} in file: {}".format(i, fname))
            print(l[7], l[8])
            print("P0: ", P0_)

    K = np.asarray(K)  # K - Wavenumber
    P0 = np.asarray(P0)  # P0
    Power = np.asarray(Power)  # P0*E + P1
    D = np.asarray(D)  # Data/E0
    X = np.asarray(X)  # Expect(k)

    return K, P0, P0_h, Power, D, X


def plot_power(fname="test.c0", axes=None):

    K, P0, P0_h, Power, D, X = extract_sbf_result(fname)

    if axes is None:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7))
    else:
        ax1, ax2 = axes

    ax1.plot(K, Power, "k:")
    ax1.plot(K, D, ".")
    ax1.set_yscale("log")
    ax1.set_xlim(-5, 170)
    ax1.set_ylim(2, 500)
    ax1.set_yticks([3, 10, 30, 100, 300])
    ax1.set_yticklabels(["3", "10", "30", "100", "300"])

    ax2.axhline(P0_h, color="k", ls=":")
    ax2.plot(K, P0, ".")
    ax2.set_xlim(-5, 170)
    ax2.set_ylim(P0_h - 40, P0_h + 40)

    ax1.set_ylabel("Power", fontsize=16)
    ax2.set_ylabel("P0", fontsize=16)
    ax2.set_xlabel("Wavenumber", fontsize=16)

    x_ax1, y_ax1 = set_axes(ax1, fontsize=12)
    y_ax1.set_yscale("log")
    y_ax1.set_yticklabels([])
    y_ax1.set_yticks([3, 10, 30, 100, 300])

    x_ax2, y_ax2 = set_axes(ax2, fontsize=12)
    x_ax2.set_ylim(P0_h - 40, P0_h + 40)
    y_ax2.set_yticks(ax2.get_yticks())

    plt.subplots_adjust(wspace=0, hspace=0)

    return ax1, ax2, (K, X)


##############################################################
def plot_power_results(
    DIR,
    PSF,
    PRF,
    RSD,
    PTM,
    power_file="power.c0",
    tv_file="tv.jpg",
    K_fit_range=[25, 60],
):

    with open(power_file, "r") as f:
        lines = f.readlines()

    XY = [l.strip("(,)") for l in lines[1].split(" ") if l != ""]
    X0 = np.int(XY[4])
    Y0 = np.int(XY[5])

    KS0 = K_fit_range[0]
    KS1 = K_fit_range[1]

    fig = plt.figure(figsize=(13, 7))

    gs1 = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    gs1.update(left=0.05, right=0.70, wspace=0.05)
    ax1 = plt.subplot(gs1[0])
    ax2 = plt.subplot(gs1[1])

    _, _, KX = plot_power(fname=power_file, axes=(ax1, ax2))

    # gs2 = gridspec.GridSpec(2,1, height_ratios=[0.7, 0.7])
    # gs2.update(left=0.70, right=0.98, hspace=0.05)
    # ax3 = plt.subplot(gs2[0])
    # ax4 = plt.subplot(gs2[1])

    # definitions for the axes
    left, width = 0.75, 0.20
    bottom, height = 0.55, 0.37

    ax3 = plt.axes([left, bottom, width, height])
    ax4 = plt.axes([left, 0.1, width, height])

    cwd = os.getcwd()
    os.chdir(DIR)

    im_1 = mpimg.imread(tv_file)
    ax3.imshow(np.flipud(im_1))
    dx = im_1.shape[0] // 2
    dy = im_1.shape[1] // 2
    r0 = np.int(XY[9])
    r1 = np.int(XY[10])

    ax3.set_xlim(dx - r1, dx + r1)
    ax3.set_ylim(dy - r1, dy + r1)

    monsta_script = (
        """

        rd 1 '"""
        + str(PSF)
        + """'                   ! reading psf from library
        rd 2 '"""
        + str(PRF)
        + """'                   ! eliprof final model
        rd 3 '"""
        + str(RSD)
        + """'                   ! galaxy subtracted residual
        rd 5 '"""
        + str(PTM)
        + """' bitmap            ! hybrid point source mask  (dophot+sextractor)
        mi 3 5
        
        tv 3 JPEG=tmp.jpg

    """
    )

    run_monsta(monsta_script, "monsta.pro", "monsta.log")

    im = mpimg.imread("tmp.jpg")
    X = im.shape[0] - X0
    im = im[X - dx : X + dx, Y0 - dy : Y0 + dy]
    ax4.imshow(im)

    ax3.axis("off")
    ax4.axis("off")

    # drawing concentric circles on panel 4
    ax4.add_patch(plt.Circle((dx, dy), r0, color="cyan", fill=False))
    ax4.add_patch(plt.Circle((dx, dy), r1, color="cyan", fill=False))

    # highlighting the fitting range
    ax2.axvspan(KS0, KS1, color="yellow", alpha=0.2)

    os.chdir(cwd)

    return (ax1, ax2, ax3, ax4), KX


##############################################################
class SBF_widgets:
    def __init__(self, XX, YY, X0, Y0, **config) -> None:

        self.X0 = X0
        self.Y0 = Y0
        self.XX = XX
        self.YY = YY

        self._k_fit = None
        self._radius = None
        self._angle = None
        self._psf_k = None

        self.DIR = config.get("DIR")
        self.PSF = config.get("PSF")
        self.PRF = config.get("PRF")
        self.RSD = config.get("RSD")
        self.PTM = config.get("PTM")

    def k_fit(self):

        if self._k_fit is None:
            self._k_fit = widgets.IntRangeSlider(
                value=[25, 60],
                min=1,
                max=200,
                step=1,
                description="K Fit Range",
                disabled=False,
                continuous_update=False,
                orientation="horizontal",
                readout=True,
                readout_format="d",
            )
        return self._k_fit

    def radius(self):

        if self._radius is None:
            options = generate_radius_options(self.X0, self.Y0, self.XX, self.YY)
            if "64 - 128" in options:
                value = "64 - 128"
            else:
                value = "32 - 64"
            self._radius = widgets.ToggleButtons(
                options=options,
                description="Radius [pixel]",
                button_style="primary",
                value=value,
                style={"button_width": "150px"},
            )
        return self._radius

    def angle(self):

        if self._angle is None:

            self._angle = widgets.IntRangeSlider(
                value=[0, 360],
                min=0,
                max=360,
                step=1,
                description="Angle [deg]",
                disabled=False,
                continuous_update=False,
                orientation="horizontal",
                readout=True,
                readout_format="d",
            )

            self._angle_b = widgets.ToggleButtons(
                options=[
                    "0 - 360",
                    "0 - 180",
                    "180 - 360",
                    "0 - 90",
                    "90 - 180",
                    "180 - 270",
                    "270 - 360",
                ],
                description="Angle",
                disabled=False,
                button_style="success",
                style={"button_width": "120px"},
            )

            def value_func(change):
                self._angle.value = [int(x) for x in self._angle_b.value.split("-")]

            self._angle_b.observe(value_func, names="value")

        return self._angle

    def psf_k(self):

        if self._psf_k is None:
            self._psf_k = widgets.IntRangeSlider(
                value=[0, 15],
                min=0,
                max=30,
                step=1,
                description="PSF K range:",
                disabled=False,
                continuous_update=False,
                orientation="horizontal",
                readout=True,
                readout_format="d",
            )

        return self._psf_k

    def settings(self):

        self.k_fit()
        self.radius()
        self.angle()
        self.psf_k()

        all = widgets.VBox(
            [self._k_fit, self._radius, self._angle, self._angle_b, self._psf_k]
        )
        display(all)

    def get_sbf(self, psf_order=4, tv_file="tv.jpg", power_file="power.c0"):

        KS0, KS1 = self._k_fit.value

        radius_range = [int(x) for x in self._radius.value.split("-")]

        if radius_range[1] < 512:
            ndim = 512
        else:
            ndim = 2 * radius_range[1]

        P0, power_file, tv_file = get_sbf_signal(
            self.DIR,
            self.PSF,
            self.PRF,
            self.RSD,
            self.PTM,
            power_file=power_file,
            tv_file=tv_file,
            K_fit_range=[KS0, KS1],
            psf_order=psf_order,
            psf_k=self._psf_k.value,
            radius=radius_range,
            angle=self._angle.value,
            center=[self.X0, self.Y0],
            nrows=ndim,
            ncols=ndim,
        )

        _ = Logtext(power_file, "Results: Power File")
        _ = Logtext(os.path.join(self.DIR, "monsta.pro"), "Monsta Code")
        _ = Logtext(os.path.join(self.DIR, "monsta.log"), "Monsta LOG")

        (ax1, ax2, ax3, ax4), KX = plot_power_results(
            self.DIR,
            self.PSF,
            self.PRF,
            self.RSD,
            self.PTM,
            power_file=power_file,
            tv_file=tv_file,
            K_fit_range=[KS0, KS1],
        )

        return (ax1, ax2, ax3, ax4), (KX, self._psf_k.value), P0

    def get_sbf_all(self, psf_order=4, angle_range=None):

        output = {}

        KS0, KS1 = self._k_fit.value

        if angle_range is None:
            angle_range = self._angle.value

        for i, rad_range in enumerate(
            generate_radius_options(self.X0, self.Y0, self.XX, self.YY)
        ):

            radius_range = [int(x) for x in rad_range.split("-")]

            print(i, rad_range)

            if radius_range[1] < 512:
                ndim = 512
            else:
                ndim = 2 * radius_range[1]

            P0, power_file, tv_file = get_sbf_signal(
                self.DIR,
                self.PSF,
                self.PRF,
                self.RSD,
                self.PTM,
                power_file="power_c%d.txt" % i,
                tv_file="tv_c%d.jpg" % i,
                K_fit_range=[KS0, KS1],
                psf_order=psf_order,
                psf_k=self._psf_k.value,
                radius=radius_range,
                angle=angle_range,
                center=[self.X0, self.Y0],
                nrows=ndim,
                ncols=ndim,
            )

            output["c%d" % i] = P0

        return output

    def get_sbf_iter(
        self, psf_order=4, angle_range=None, radius_range=None, sky_sigma=100, n_iter=5, psf_list=None
    ):

        output = {}

        KS0, KS1 = self._k_fit.value

        if angle_range is None:
            angle_range = self._angle.value

        if radius_range is None:
            radius_range = [int(x) for x in self._radius.value.split("-")]


        for i in range(n_iter):
            sky_adjustment = sky_sigma * np.random.normal()

            
            j = 0 
            if psf_list is not None:
                N = len(psf_list)
                j = np.random.randint(N)
                psf = psf_list[j]
            else:
                psf = self.PSF

            #print(i, j, angle_range, radius_range, sky_adjustment)

            if radius_range[1] < 512:
                ndim = 512
            else:
                ndim = 2 * radius_range[1]

            P0, power_file, tv_file = get_sbf_signal(
                self.DIR,
                psf,
                self.PRF,
                self.RSD,
                self.PTM,
                power_file="power_iter.txt",
                tv_file=None,
                K_fit_range=[KS0, KS1],
                psf_order=psf_order,
                psf_k=self._psf_k.value,
                radius=radius_range,
                angle=angle_range,
                center=[self.X0, self.Y0],
                nrows=ndim,
                ncols=ndim,
                sky_adjustment=sky_adjustment,
            )

            K, P, _, _, _, _ = extract_sbf_result(power_file)
            out_dict = {}

            out_dict["P0"] = P0
            out_dict["K"] = K
            out_dict["P"] = P
            
            output["iter%03d" % i] = out_dict

        return output


##############################################################
def is_valid_radius(r, X0, Y0, XX, YY):

    if X0 - r < 0 or Y0 - r < 0:
        return False

    if X0 + r > XX or Y0 + r > YY:
        return False

    return True


##############################################################
def generate_radius_options(X0, Y0, XX, YY):

    options = []
    r0 = 32
    while is_valid_radius(r0, X0, Y0, XX, YY) and is_valid_radius(
        r0 * 2, X0, Y0, XX, YY
    ):
        options.append("{} - {}".format(r0, 2 * r0))
        r0 *= 2
    return options


##############################################################
def plot_psf_power(psf_info, K_upper_limit=70, highlight=True):

    K, Power = psf_info[0]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(K, Power)

    ax.set_title("PSF", fontsize=14)
    ax.set_xlabel("Wavenumber", fontsize=12)
    ax.set_ylabel("Power", fontsize=12)

    ax.set_xlim(-2, K_upper_limit)
    Power_ = Power[K < K_upper_limit]
    mn = np.min(Power_)
    mx = np.max(Power_)
    r = (mx - mn) * 0.15
    ax.set_ylim(mn - r, mx + r)

    if highlight:
        ax.axvspan(psf_info[1][0], psf_info[1][1], color="blue", alpha=0.2)

    return ax


##############################################################
##############################################################
##############################################################

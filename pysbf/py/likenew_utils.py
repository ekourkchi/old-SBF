import os, sys, string, math, uuid, time, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from matplotlib.ticker import NullFormatter

from .utils import *
from .likenew import *


def eq(lw):
    if "=" in lw:
        for i, w in enumerate(lw):
            if w == "=":
                break
        return lw[i + 1]
    else:
        return None


def beq(lw):
    if "=" in lw:
        for i, w in enumerate(lw):
            if w == "=":
                break
        return "_".join(lw[:i])
    else:
        return None


########################################################################################


def read_lkn6(fname):

    with open(fname, "r") as f:
        lines = f.readlines()
    ########################################################################################
    ########################################################################################
    header = {}
    ac = 0

    for i, l in enumerate(lines):
        for w in [",", "(", ")", "/", "#"]:
            l = l.replace(w, "")
        words = l.split()
        if words != []:
            if words[1] == "lines":
                header["_".join([w for w in words[1:] if w != "of"])] = int(words[0])
            if words[1] == "points":
                header["_".join([w for w in words[1:] if w != "rm"])] = int(words[0])

            if "=" in words:
                if words[0] == "Distance":
                    if "assumed" in words[1:]:
                        header["distance_assumed"] = float(eq(words))
                    if "derived" in words[1:]:
                        header["distance_derived"] = float(eq(words))
                elif words[0] == "Annulus":
                    ac += 1
                    header[beq(words) + str(ac)] = float(eq(words))
                else:
                    header[beq(words)] = float(eq(words))

        else:
            break

    ########################################################################################
    ########################################################################################
    table1 = lines[i + 1 :]

    E = []
    for e in [x.strip().split() for x in table1[2].strip().split("<r<")]:
        E = E + e
    R_dict = {}
    jj = 1
    for ii in range(0, len(E), 2):
        R_dict["r" + str(jj)] = [E[ii], E[ii + 1]]
        jj += 1

    for i1, l in enumerate(table1):
        words = l.split()
        if words[0] == "r":
            break

    table1 = table1[i1 + 1 :]
    T1 = []
    for i2, row in enumerate(table1):
        if row.split() == []:
            break
        r = [float(x) for x in row.split()]
        T1.append(r)

    T1 = np.asarray(T1)

    df1 = pd.DataFrame()

    df1["r"] = T1[:, 0]
    df1["frac"] = T1[:, 1]

    for idx in [2]:
        df1["N"] = T1[:, idx]
        df1["n_m"] = T1[:, idx + 1]
        df1["m"] = T1[:, idx + 2]

    jj = 1
    for idx in range(5, T1.shape[1], 2):
        df1["N" + str(jj)] = T1[:, idx]
        df1["n_m" + str(jj)] = T1[:, idx + 1]
        jj += 1

    ########################################################################################
    ########################################################################################
    table2 = table1[i2 + 1 :]

    columns = table2[2].strip().split()

    jj = 1
    for i in range(2, len(columns), 2):
        columns[i] = columns[i] + "_" + str(jj)
        columns[i + 1] = columns[i + 1] + "_" + str(jj)
        jj += 1

    T2 = []
    for i3, row in enumerate(table2[3:]):
        if row.split() == []:
            break
        r = [float(x) for x in row.split()]
        T2.append(r)

    T2 = np.asarray(T2)

    df2 = pd.DataFrame()

    for i, c in enumerate(columns):
        df2[c] = T2[:, i]

    ########################################################################################
    ########################################################################################
    table3 = table2[i3 + 6 :]

    columns = table3[0].strip().split()
    columns[1] = "m_bright"
    columns[2] = "m_faint"

    T3 = []
    for i4, row in enumerate(table3[1:]):
        if row.split() == []:
            break
        r = [float(x) for x in row.split()]
        T3.append(r)

    T3 = np.asarray(T3)

    df3 = pd.DataFrame()

    for i, c in enumerate(columns):
        df3[c] = T3[:, i]
    ########################################################################################
    return header, (df1, df2, df3), R_dict
    ########################################################################################


def likenew_plot_old(lkn6_name, header, tables, radii):

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))

    df1, df2, df3 = tables

    for i in range(6):

        ii = i % 3
        jj = (i - ii) // 3
        ax = axes[jj][ii]

        idx = (i + 1) % 6 + 1
        ax.plot(df1.m, df1["n_m" + str(idx)], "D")
        ax.plot(df2.m, df2["GC_" + str(idx)], "-")
        ax.plot(df2.m, df2["both_" + str(idx)], "-")
        ax.plot(df2.m, df2.gxy, "k:")

        r1 = radii["r" + str(idx)][0]
        r2 = radii["r" + str(idx)][1]
        ax.text(
            0.1, 0.85, r1 + "<r<" + r2, fontsize=16, color="k", transform=ax.transAxes
        )

        ax.set_yscale("log")
        x_ax, y_ax = set_axes(ax, xlim=(18.5, 27.5), ylim=(0.5, 2000), fontsize=14)
        y_ax.set_yscale("log")

        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        x_ax.xaxis.set_minor_locator(MultipleLocator(0.5))

        ax.yaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        plt.yticks([1, 10, 100, 1000], ("1", "10", "100", "1000"))

        plt.setp(y_ax.get_yticklabels(), visible=False)
        if ii > 0:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel(r"$n \/\/ (mag^{-1} \/ arcmin^{-2})$", fontsize=14)

        if jj != 1:
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            ax.set_xlabel(r"magnitude", fontsize=14)

    plt.subplots_adjust(hspace=0, wspace=0)

    s = "dist = %.2f        GCLF = %.4f       S/N = %.1f" % (
        header["distance_derived"],
        header["Cmax"],
        header["SN_limit"],
    )
    _ = plt.suptitle(
        s + "          " + lkn6_name, fontsize=12, y=0.91, x=0.4, color="maroon"
    )

    return axes


############################################################
def likenew_plot(lkn6_name, header, tables, radii):

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))

    df1, df2, df3 = tables

    for i in range(6):

        ii = i % 3
        jj = (i - ii) // 3
        ax = axes[jj][ii]

        idx = (i + 1) % 6 + 1
        ax.plot(df1.m, df1["n_m" + str(idx)], "D")
        ax.plot(df2.m, df2["GC_" + str(idx)], "-")
        ax.plot(df2.m, df2["both_" + str(idx)], "-")
        ax.plot(df2.m, df2.gxy, "k:")

        r1 = radii["r" + str(idx)][0]
        r2 = radii["r" + str(idx)][1]
        ax.text(
            0.1, 0.85, r1 + "<r<" + r2, fontsize=16, color="k", transform=ax.transAxes
        )

        ax.set_yscale("log")
        #         x_ax, y_ax = set_axes(ax, xlim=(18.5,27.5), ylim=(0.5,2000), fontsize=14)
        #         y_ax.set_yscale("log")
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        ax.tick_params(which="major", length=8, width=1.0, direction="in")
        #         if minor:
        ax.tick_params(
            which="minor", length=4, color="#000033", width=1.2, direction="in"
        )

        ax.set_xlim(18.5, 27.5)
        ax.set_ylim(0.5, 2000)

        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        #         x_ax.xaxis.set_minor_locator(MultipleLocator(0.5))

        #         ax.yaxis.set_major_formatter(NullFormatter())
        #         ax.yaxis.set_minor_formatter(NullFormatter())
        plt.yticks([1, 10, 100, 1000], ("1", "10", "100", "1000"))

        #         plt.setp(y_ax.get_yticklabels(), visible=False)
        if ii > 0:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel(r"$n \/\/ (mag^{-1} \/ arcmin^{-2})$", fontsize=14)

        if jj != 1:
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            ax.set_xlabel(r"magnitude", fontsize=14)

    plt.subplots_adjust(hspace=0, wspace=0)

    s = "dist = %.2f        GCLF = %.4f       S/N = %.1f" % (
        header["distance_derived"],
        header["Cmax"],
        header["SN_limit"],
    )
    _ = plt.suptitle(
        s + "          " + lkn6_name, fontsize=12, y=0.91, x=0.4, color="maroon"
    )

    return axes


## END

########################################################################################


def pylike(
    like_pack,
    gclf_width=1.4,
    distance=60,
    bright_cutoff=20.5,
    kscale=1.2,
    plot=False,
    **Config
):

    name = Config["name"]
    inFolder = Config["inFolder"]
    objRoot = Config["objRoot"]
    lkn, lkn6, ptm6, resid = like_pack

    cwd = os.getcwd()
    xcmd("cp {}/{}/{}j.dmask {}.".format(inFolder, name, name, objRoot), verbose=True)
    os.chdir(objRoot)
    likenew = LikeNew(in_folder="./")
    likenew.run(
        gclf_width=gclf_width,
        distance=distance,
        bright_cutoff=bright_cutoff,
        kscale=kscale,
        fname=lkn,
        mname="{}j.dmask".format(name),  # initial mask
        verbose=True,
    )
    os.chdir(cwd)

    if plot:
        lkn6_name = os.path.join(objRoot, lkn6)
        header, tables, radii = read_lkn6(lkn6_name)
        _ = likenew_plot(lkn6_name, header, tables, radii)


########################################################################################

########################################################################################
########################################################################################
########################################################################################

########################################################################################
########################################################################################

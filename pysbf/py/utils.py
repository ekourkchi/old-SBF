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

import numpy
from astropy.table import Table
import requests
from PIL import Image
from io import BytesIO
import pylab, os, sys

import copy


##############################################################
class MissingGalaxyCenter(Exception):
    pass


##############################################################
def get_sex_catal_ColName(catalName):

        with open(catalName, 'r') as f:

            lines = f.readlines()

        cols = {}
        i = 0 
        while lines[i].split()[0]=="#":
            j = int(lines[i].split()[1])
            name = lines[i].split()[2]
            
            n = 1
            cols[i] = [j, name, n]
            
            if i>0:
                n = cols[i][0]-cols[i-1][0]
                cols[i-1][2] = n
          
            i+=1
        
        col_names = []
        n_head = i
        for m in range(i):
            s = cols[m][2]
            name = cols[m][1]
            if s==1:
                col_names.append(name)
            else:
                for p in range(s):
                    col_names.append(name+'_'+str(p+1))
                            
        return n_head, col_names

##############################################################
def get_sextract_catal_df(catalName):

    n_head, col_names = get_sex_catal_ColName(catalName)
    catal_df = pd.read_csv(catalName, delimiter=r"\s+", skiprows=n_head, 
                            header = None, names = col_names)

    return catal_df


##############################################################
def tv(fits_name, xlim=[0,1], ylim=[0,1], ax=None, options="", zoom=1., **kwarg):

    size = 1./zoom

    img = get_img(fits_name, options=options)
    
    x_max, y_max, _ = img.shape

    xlim = [int(x_max*l) for l in xlim]
    ylim = [int(y_max*l) for l in ylim]

    if "XY" in kwarg:

        try:
            XY = kwarg["XY"]
            X0 = int(XY[0])
            Y0 = int(XY[1])

            x0 = np.max([X0 - int(size * x_max) , 0])
            x1 = np.min([X0 + int(size * x_max) , x_max])
            y0 = np.max([Y0 - int(size * y_max) , 0])
            y1 = np.min([Y0 + int(size * y_max) , y_max])

            xlim = [x0, x1]
            ylim = [y0, y1]
        except:
            print("Enter the galaxy center !")
            pass        
        
    if ax is None:
        plt.figure(figsize=(10, 10))
        plt.subplot(111)
        ax = plt.gca()

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    imgplot = ax.imshow(np.flipud(img))
    
    return ax
##############################################################

def get_img(fits_file, ax=None, options=""):

    jpg_name = 'tv.jpg'

    ## Monsta script
    script = """
    rd 1 '"""+fits_file+"""'
    tv 1 """+options+""" JPEG="""+jpg_name+"""
    q

    """

    run_monsta(script, 'tv.pro', 'tv.log')
    
    xcmd("rm tv.pro & rm tv.log &", verbose=False)
    
    img = mpimg.imread(jpg_name)

    xcmd("rm tv.jpg &", verbose=False)
    
    return img
##############################################################

def run_monsta(script, Monsta_pro, Monsta_log, monsta="monsta", silent=False):

    with open(Monsta_pro, "w") as f:
        f.write(script)

    cmd = monsta + " " + Monsta_pro
    xcmd(cmd + " > " + Monsta_log, verbose=False)

    with open(Monsta_log) as f:
        text = f.read()

    with open(Monsta_log) as f:
        text = f.read()
    Tsplit = text.upper().split()
    if "ERROR" in Tsplit or "SEGMENTATION FAULT" in Tsplit:
        if not silent:
            print("===== MONSTA =====")
            for i, line in enumerate(script.split("\n")):
                print(f"{i+1:>4}{line}")
            print("===================")
            print(text)
            return text
        else:
            return "[Warning] MONSTA NOT OK !"
    else:
        return "OK"



##############################################################


def get_extinction(ra, dec):
    URL = (
        """http://ned.ipac.caltech.edu/ffs/sticky/CmdSrv?cmd=tableSave&request=%7B%22startIdx%22%3A0%2C%22pageSize%22%3A2147483647%2C%22filters%22%3A%22%22%2C%22source%22%3A%22http%3A%2F%2Fned.ipac.caltech.edu%2Fcgi-bin%2Fcalc%3Fin_csys%3DEquatorial%26out_csys%3DEquatorial%26in_equinox%3DJ2000.0%26out_equinox%3DJ2000.0%26obs_epoch%3D2000.0%26of%3Dxml_main%26ext%3D1%26lon%3D"""
        + str(ra)
        + """d%26lat%3D"""
        + str(dec)
        + """d%22%2C%22alt_source%22%3A%22http%3A%2F%2Fned.ipac.caltech.edu%2Fcgi-bin%2Fcalc%3Fin_csys%3DEquatorial%26out_csys%3DEquatorial%26in_equinox%3DJ2000.0%26out_equinox%3DJ2000.0%26obs_epoch%3D2000.0%26of%3Dxml_main%26ext%3D1%26lon%3D27.43241667d%26lat%3D35.78563889d%22%2C%22META_INFO%22%3A%7B%22title%22%3A%22extinctions%22%2C%22tbl_id%22%3A%22extinctions_table%22%2C%22col.Refcode.PrefWidth%22%3A%2220%22%2C%22col.Refcode%20of%20the%20publications.PrefWidth%22%3A%2228%22%2C%22col.Central%20Wavelength.FmtDisp%22%3A%22%25.2f%22%2C%22col.The%20Galactic%20extinction.FmtDisp%22%3A%22%25.3f%22%2C%22resultSetRequest%22%3A%22%7B%5C%22RequestClass%5C%22%3A%5C%22ServerRequest%5C%22%2C%5C%22META_INFO%5C%22%3A%7B%5C%22col.Refcode%20of%20the%20publications.PrefWidth%5C%22%3A%5C%2228%5C%22%2C%5C%22tbl_id%5C%22%3A%5C%22extinctions_table%5C%22%2C%5C%22col.Central%20Wavelength.FmtDisp%5C%22%3A%5C%22%25.2f%5C%22%2C%5C%22col.The%20Galactic%20extinction.FmtDisp%5C%22%3A%5C%22%25.3f%5C%22%2C%5C%22title%5C%22%3A%5C%22extinctions%5C%22%2C%5C%22col.Refcode.PrefWidth%5C%22%3A%5C%2220%5C%22%7D%2C%5C%22tbl_id%5C%22%3A%5C%22extinctions_table%5C%22%2C%5C%22startIdx%5C%22%3A0%2C%5C%22pageSize%5C%22%3A1000%2C%5C%22alt_source%5C%22%3A%5C%22http%3A%5C%5C%2F%5C%5C%2Fned.ipac.caltech.edu%5C%5C%2Fcgi-bin%5C%5C%2Fcalc%3Fin_csys%3DEquatorial%26out_csys%3DEquatorial%26in_equinox%3DJ2000.0%26out_equinox%3DJ2000.0%26obs_epoch%3D2000.0%26of%3Dxml_main%26ext%3D1%26lon%3D27.43241667d%26lat%3D35.78563889d%5C%22%2C%5C%22id%5C%22%3A%5C%22IpacTableFromSource%5C%22%2C%5C%22source%5C%22%3A%5C%22http%3A%5C%5C%2F%5C%5C%2Fned.ipac.caltech.edu%5C%5C%2Fcgi-bin%5C%5C%2Fcalc%3Fin_csys%3DEquatorial%26out_csys%3DEquatorial%26in_equinox%3DJ2000.0%26out_equinox%3DJ2000.0%26obs_epoch%3D2000.0%26of%3Dxml_main%26ext%3D1%26lon%3D27.43241667d%26lat%3D35.78563889d%5C%22%2C%5C%22ffSessionId%5C%22%3A%5C%22FF-Session-1647588578055%5C%22%7D%22%2C%22resultSetID%22%3A%22DATA%22%7D%2C%22id%22%3A%22IpacTableFromSource%22%2C%22RequestClass%22%3A%22ServerRequest%22%2C%22ffSessionId%22%3A%22FF-Session-1647588578055%22%2C%22inclCols%22%3A%22%5C%22Bandpass%5C%22%2C%5C%22Central%20Wavelength%5C%22%2C%5C%22The%20Galactic%20extinction%5C%22%2C%5C%22Refcode%20of%20the%20publications%5C%22%22%7D&file_name=extinctions&file_format=csv"""
    )
    cmd = 'curl -X "POST" "' + URL + '" > tmp.csv'
    xcmd(cmd, verbose=False)
    df = pd.read_csv("tmp.csv")
    df[df.columns[:10]].head()

    return df


##############################################################
def set_axes(
    ax,
    xlim=None,
    ylim=None,
    fontsize=16,
    twinx=True,
    twiny=True,
    minor=True,
    inout="in",
):

    if not ylim is None:
        ax.set_ylim(ylim)
    else:
        ylim = ax.get_ylim()

    if not xlim is None:
        ax.set_xlim(xlim)
    else:
        xlim = ax.get_xlim()

    ax.tick_params(which="major", length=6, width=1.0, direction=inout)
    #         if minor:
    ax.tick_params(which="minor", length=3, color="#000033", width=1.0, direction=inout)

    if twiny:
        y_ax = ax.twinx()
        y_ax.set_ylim(ylim)
        y_ax.set_yticklabels([])
        y_ax.minorticks_on()
        y_ax.tick_params(which="major", length=6, width=1.0, direction=inout)
        if minor:
            y_ax.tick_params(
                which="minor", length=3, color="#000033", width=1.0, direction=inout
            )

    if twinx:
        x_ax = ax.twiny()
        x_ax.set_xlim(xlim)
        x_ax.set_xticklabels([])
        x_ax.minorticks_on()
        x_ax.tick_params(which="major", length=6, width=1.0, direction=inout)
        if minor:
            x_ax.tick_params(
                which="minor", length=3, color="#000033", width=1.0, direction=inout
            )

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    return x_ax, y_ax


##############################################################


def ellipse_polyline(ellipses, n=10000):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    st = np.sin(t)
    ct = np.cos(t)
    result = []
    for x0, y0, a, b, angle in ellipses:
        angle = np.deg2rad(angle)
        sa = np.sin(angle)
        ca = np.cos(angle)
        p = np.empty((n, 2))
        p[:, 0] = x0 + a * ca * ct - b * sa * st
        p[:, 1] = y0 + a * sa * ct + b * ca * st
        result.append(p)
    return result


## return True if circumference of two ellipses intersect
def ell_itersect(Ell_A, Ell_B):

    ell_a = (Ell_A[0][0], Ell_A[0][1], Ell_A[1], Ell_A[2], Ell_A[3])
    ell_b = (Ell_B[0][0], Ell_B[0][1], Ell_B[1], Ell_B[2], Ell_B[3])

    if not np.isnan(np.sum(ell_a)) and not np.isnan(np.sum(ell_b)):

        ellipses = [ell_a, ell_b]
        a, b = ellipse_polyline(ellipses)
        ea = LinearRing(a)
        eb = LinearRing(b)
        try:
            mp = ea.intersection(eb)
            return mp.type == "MultiPoint"
        except:
            print((ell_a, ell_b))
            return True  # likely invalid
    else:
        return False


##############################################################

## returns True if the surface of two ellipoids cross each other
def ellipsoid_intersection_test(Ell_A, Ell_B, tau=1):

    Sigma_A, mu_A = Ell_A[4], Ell_A[0]
    Sigma_B, mu_B = Ell_B[4], Ell_B[0]

    lambdas, Phi, v_squared = ellipsoid_intersection_test_helper(
        Sigma_A, Sigma_B, mu_A, mu_B
    )
    res = minimize_scalar(
        ellipsoid_K_function, bracket=[0.0, 0.5, 1.0], args=(lambdas, v_squared, tau)
    )
    return res.fun[0] >= 0


def ellipsoid_intersection_test_helper(Sigma_A, Sigma_B, mu_A, mu_B):
    lambdas, Phi = eigh(Sigma_A, b=Sigma_B)
    v_squared = np.dot(Phi.T, mu_A - mu_B) ** 2
    return lambdas, Phi, v_squared


def ellipsoid_K_function(ss, lambdas, v_squared, tau):
    ss = np.array(ss).reshape((-1, 1))
    lambdas = np.array(lambdas).reshape((1, -1))
    v_squared = np.array(v_squared).reshape((1, -1))
    return 1.0 - (1.0 / tau ** 2) * np.sum(
        v_squared * ((ss * (1.0 - ss)) / (1.0 + ss * (lambdas - 1.0))), axis=1
    )


##############################################################

## This function allows to execute the OS commands
def xcmd(cmd, verbose=True):
    """Runs an OS command
    :param cmd: terminal command
    :type cmd: ``str``
    :param verbose: printing the details, default True
    :type verbose: ``boolean``
    :return: OS outputs
    :rtype: ``str``
    """

    if verbose:
        print("\n" + cmd)

    tmp = os.popen(cmd)
    output = ""
    for x in tmp:
        output += x
    if "abort" in output:
        failure = True
    else:
        failure = tmp.close()
    if False:
        print("execution of %s failed" % cmd)
        print("error is as follows", output)
        sys.exit()
    else:
        return output


##############################################################


def getimages(ra, dec, size=240, filters="grizy"):

    """Query ps1filenames.py service to get a list of images

    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """

    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = (
        "{service}?ra={ra}&dec={dec}&size={size}&format=fits" "&filters={filters}"
    ).format(**locals())
    table = Table.read(url, format="ascii")
    return table


def geturl(
    ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False
):

    """Get URL for images in the table

    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """

    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg", "png", "fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra, dec, size=size, filters=filters)
    url = (
        "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
        "ra={ra}&dec={dec}&size={size}&format={format}"
    ).format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table["filter"]]
    table = table[numpy.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0, len(table) // 2, len(table) - 1]]
        for i, param in enumerate(["red", "green", "blue"]):
            url = url + "&{}={}".format(param, table["filename"][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table["filename"]:
            url.append(urlbase + filename)
    return url


def getcolorim(ra, dec, size=240, output_size=None, filters="grizy", format="jpg"):

    """Get color image at a sky position

    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """

    if format not in ("jpg", "png"):
        raise ValueError("format must be jpg or png")
    url = geturl(
        ra,
        dec,
        size=size,
        filters=filters,
        output_size=output_size,
        format=format,
        color=True,
    )
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


def getgrayim(ra, dec, size=240, output_size=None, filter="g", format="jpg"):

    """Get grayscale image at a sky position

    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filter = string with filter to extract (one of grizy)
    format = data format (options are "jpg", "png")
    Returns the image
    """

    if format not in ("jpg", "png"):
        raise ValueError("format must be jpg or png")
    if filter not in list("grizy"):
        raise ValueError("filter must be one of grizy")

    #     print(ra,dec, size, filter)

    url = geturl(
        ra, dec, size=size, filters=filter, output_size=output_size, format=format
    )

    #     print(url)
    r = requests.get(url[0])
    im = Image.open(BytesIO(r.content))

    return im


##############################################################


def createDir(folderPath):
    """generating a directory/folder if it doesn't exist
    :param folderPath: path to the desired folder
    :type folderPath: ``str``
    :return: True is created, False if the folder already exists
    :rtype: ``str``
    """

    if not os.path.exists(folderPath):
        os.makedirs(folderPath)
        return True
    else:
        False


##############################################################


def unit_vecort(angle):

    theta = angle * np.pi / 180.0
    ux = np.cos(theta)
    uy = np.sin(theta)

    return (ux, uy)


def get_Sigma(Smajor, Sminor, angle):

    u1 = unit_vecort(angle)
    u2 = unit_vecort(angle - 90)
    R = np.vstack((u1, u2))
    L = R.T
    X = np.diag((1.0 / Smajor ** 2, 1.0 / Sminor ** 2))
    A = np.matmul(np.matmul(L, X), R)

    return np.linalg.inv(A)


##############################################################


def make_Ellipse(center, Smajor, Sminor, angle):

    Sigma = get_Sigma(Smajor, Sminor, angle)
    return [np.asarray(center), Smajor, Sminor, angle, Sigma]


def plot_Ellipse(Ell, ax=None, **kwargs):

    if ax is None:
        ax = plt.gca()

    e = patches.Ellipse(
        Ell[0], width=2 * Ell[1], height=2 * Ell[2], angle=Ell[3], **kwargs
    )
    ax.add_patch(e)

    return ax


##############################################################


def list_Ell(data_frame):

    Ell_list = []

    for i in range(len(data_frame)):

        df = data_frame.iloc[i]
        Smajor = df.Rmaj
        Sminor = Smajor * (1.0 - df.ellip)
        angle = df.alpha + 90
        center = (df.x0, df.y0)
        Sigma = get_Sigma(Smajor, Sminor, angle)
        Ell = make_Ellipse(center, Smajor, Sminor, angle)

        Ell_list.append(Ell)

    return Ell_list


##############################################################

# df is a Pandas series that holds one ellipse
def plot_E(df, ax=None, **kwargs):

    Smajor = df.Rmaj
    Sminor = Smajor * (1.0 - df.ellip)

    angle = df.alpha + 90
    center = (df.x0, df.y0)
    Sigma = get_Sigma(Smajor, Sminor, angle)
    Ell = make_Ellipse(center, Smajor, Sminor, angle)

    if ax is None:
        ax = plt.gca()

    plot_Ellipse(Ell, ax=ax, **kwargs)

    return ax


##############################################################

##############################################################

def imOpen(inFits):

    hdu_list = fits.open(inFits)
    imarray = hdu_list[0].data
    header = hdu_list[0].header

    return imarray, header


##############################################################


def seg2mask(inFits, outMask, overwrite=True, seg_num=0, good_segments=[0], object_mask=False):

    imarray, header = imOpen(inFits)

    for good_pix in good_segments:
        imarray[((imarray == good_pix))] = -1  # good pixel
    if seg_num != 0:
        imarray[((imarray == seg_num))] = -1  # good pixel

    if object_mask:
        imarray[((imarray == 0))] = -1  # backgrounds --> good pixel
    imarray[(imarray != -1)] = 0  # masked
    imarray[imarray == -1] = 1  # good pixel

    fits.writeto(outMask, np.float32(imarray), header, overwrite=overwrite)

    return imarray


##############################################################

# def seg2mask(inFits, outMask, overwrite=True, seg_num=0):

#     imarray, header = imOpen(inFits)
#     imarray[imarray==0] = -1
#     imarray[imarray>seg_num] = 0
#     imarray[imarray==-1] = 1

#     fits.writeto(outMask, np.float32(imarray), header, overwrite=overwrite)

#     return imarray


##############################################################


def plot_2darray(array, ax=None, shape=10):

    x_max, y_max = array.shape

    if ax is None:
        plt.figure(figsize=(shape, shape))
        plt.subplot(111)
        ax = plt.gca()

    ax.set_xlim([0, x_max])
    ax.set_ylim([0, y_max])

    imgplot = ax.imshow(array)

    return ax


##############################################################


def Xellipses(ells):
    N = len(ells)
    n_cross = 0
    for n in range(0, N - 1):
        ell_a = ells[n]
        for m in range(n + 1, N):
            ell_b = ells[m]
            cross = ell_itersect(ell_a, ell_b)
            if cross:
                n_cross += 1
    return n_cross


##############################################################
def open_log_df(logName):

    with open(logName, "r") as f:
        lines = f.readlines()

    r = logName.rsplit("/", 1)[0]
    if r==logName:
        R = "./"
    else:
        R = r + '/'

    if exists(R+"log.tmp"):
        # print("backing up log.temp")
        xcmd("cp "+R+"log.tmp "+R+"log.tmp.back", verbose=False)
        xcmd("rm "+R+"log.tmp", verbose=False)

    with open(R+"log.tmp", "w") as f:

        First = True
        for l in lines:
            if "######" in l:
                if First:
                    First = False
                    continue
                else:
                    break
            else:
                f.write(l)

    df = pd.read_csv(R+"log.tmp")
    for col in df.columns:
        df = df.rename(columns={col: col.strip()})

    for col in df.columns:
        df[col] = df[col].apply(lambda x: x.strip())

    df.set_index("index", inplace=True)
    df.index = [x.strip() for x in df.index]

    return df


def get_obj_params(df):

    params = {}

    main_key = "backSextract"
    params[main_key] = {}
    params[main_key]["threshold"] = np.float(df.loc["BXT"].value.strip())

    main_key = "naiveSextract"
    params[main_key] = {}
    params[main_key]["minarea"] = np.float(df.loc["NXM"].value.strip())
    params[main_key]["threshold"] = np.float(df.loc["NXT"].value.strip())
    params[main_key]["smooth"] = np.float(df.loc["NXS"].value.strip())

    main_key = "basic_elliprof"
    params[main_key] = {}
    params[main_key]["r0"] = np.float(df.loc["BER"].value.strip())
    params[main_key]["c_kron"] = np.float(df.loc["BEC"].value.strip())
    params[main_key]["sky_factor"] = np.float(df.loc["BES"].value.strip())
    params[main_key]["k_ellipse"] = np.float(df.loc["BEK"].value.strip())
    params[main_key]["option"] = df.loc["BEO"].value.strip()

    main_key = "second_elliprof"
    params[main_key] = {}
    params[main_key]["r0"] = np.float(df.loc["SER"].value.strip())
    params[main_key]["c_kron"] = np.float(df.loc["SEC"].value.strip())
    params[main_key]["sky_factor"] = np.float(df.loc["SEF"].value.strip())
    params[main_key]["k_ellipse"] = np.float(df.loc["SEK"].value.strip())
    params[main_key]["option"] = df.loc["SEO"].value.strip()
    params[main_key]["minarea"] = np.float(df.loc["SEM"].value.strip())
    params[main_key]["threshold"] = np.float(df.loc["SET"].value.strip())
    params[main_key]["smooth"] = np.float(df.loc["SES"].value.strip())
    params[main_key]["renuc"] = np.float(df.loc["SEN"].value.strip())

    return params

#########################################################
## `get_RMS`
# Defining the `rms` of the flux deviations from the `r^1/4` profile, extrapolated in the outer regions of the target galaxy
#########################################################
def get_RMS(obj, r0, r1, nr, sky_factor, options=""):
    """

    Returns:
        - rms: the rms of deviations
        - n_cross: number of ellipses crossing each other

    """

    sky = int(sky_factor * obj.sky_med)
    n_cross = 0

    if (
        obj.elliprof(
            r0,
            r1,
            nr=nr,
            sky=sky,
            niter=10,
            mask=1,
            model_mask=0,
            model=1000,
            options=options,
        )
        != "OK"
    ):
        n_cross += 1

    model = 1000
    n_cross += Xellipses(obj.list_ellipses(model=1000))
    root = obj.objRoot
    suffix = ".%03d" % model

    ellipseFile = root + "/elliprof" + suffix
    df = pd.read_csv(ellipseFile, delimiter=r"\s+", skiprows=7)
    df = df.apply(pd.to_numeric, errors="coerce")
    x = df.Rmaj ** 0.25
    y = 2.5 * np.log10(df.I0)

    maxX = np.max(x)
    minX = np.min(x)
    dx = maxX - minX
    x1 = 0.70 * dx + minX
    x2 = maxX - 0.10 * dx
    x3 = maxX - 0.10 * dx
    x0 = x[((x < x2) & (x > x1))]
    y0 = y[((x < x2) & (x > x1))]

    m, b = np.polyfit(x0, y0, 1)

    x_data = x[((x >= x3))]
    y_data = y[((x >= x3))]
    y_model = m * x_data + b

    rms = np.sqrt(np.mean((y_data.values - y_model.values) ** 2))

    return rms, n_cross


#########################################################

def get_f(obj, r0, r1, nr, options=""):
    """


    Returns:
        - func: a function that gets the sky_factor and returns the rms of deviations
        This function in its heart uses function `get_RMS`

    """

    def func(sky_factor):

        rms, n_cross = get_RMS(obj, r0, r1, nr, sky_factor, options=options)

        sig = rms

        if sig > 10 or np.isnan(sig) or n_cross > 0:
            sig = 10

        return -sig

    return func


#########################################################

def populate_dict(objDict, inDict):

    if objDict is None:
        objDict = {}

    for key in inDict:
        objDict[key] = inDict[key]

    return objDict

####################################### Set Axes
def set_axes(ax, xlim=None, ylim=None, fontsize=16, twinx=True, twiny=True, minor=True, inout='in'):
        
        if not ylim is None:
            ax.set_ylim(ylim)
        else:
            ylim = ax.get_ylim() 
            
        if not xlim is None:    
            ax.set_xlim(xlim) 
        else:
            xlim = ax.get_xlim()
            
        ax.tick_params(which='major', length=8, width=1., direction=inout)
#         if minor:
        ax.tick_params(which='minor', length=4, color='#000033', width=1.2, direction=inout)  
        
        if twiny:
            y_ax = ax.twinx()
            y_ax.set_ylim(ylim)
            y_ax.set_yticklabels([])
            y_ax.minorticks_on()
            y_ax.tick_params(which='major', length=8, width=1., direction=inout)
            if minor:
                y_ax.tick_params(which='minor', length=4, color='#000033', width=1.2, direction=inout) 
        
        if twinx:
            x_ax = ax.twiny()
            x_ax.set_xlim(xlim)
            x_ax.set_xticklabels([])
            x_ax.minorticks_on()
            x_ax.tick_params(which='major', length=8, width=1.2, direction=inout)
            if minor:
                x_ax.tick_params(which='minor', length=4, color='#000033', width=1.2, direction=inout)  
                
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)


        return x_ax, y_ax
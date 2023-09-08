from .utils import *


class ReidualMaksedImageNotFound(Exception):
    pass


##############################################################
def DoPhot(
    model=0,
    renuc=1,
    residual=None,
    image_in=None,
    image_out=None,
    object_out=None,
    **Config
):

    root = Config["objRoot"]
    name = Config["name"]
    inFolder = Config["inFolder"]
    configFolder = Config["configFolder"]
    X0 = Config["X0"]
    Y0 = Config["Y0"]
    SKY = Config["SKY"]

    suffix = ".%03d" % model

    if residual is None:
        residName = root + "/resid" + suffix
    else:
        residName = root + "/" + residual

    modelName = root + "/model" + suffix

    objFits = inFolder + "{}/{}j.fits".format(name, name)
    hdu_list = fits.open(objFits)
    image_data = hdu_list[0].data
    TOP = str(int(np.ceil(np.max(image_data))))  #
    BOTTOM = str(int(np.floor(np.min(image_data))))  # ~1000, -2000

    THRESHMAX = TOP  # 1850000   # around saturation level
    THRESHMIN = 1000
    SCALEAPRADIUS = 3.0
    RDNOISE = 200
    FWHM = 1.4  # same most of the time

    if image_out is None:
        IMAGE_OUT = "dpr" + suffix  ## temp file to use in ds9
    else:
        IMAGE_OUT = image_out

    if object_out is None:
        OBJECTS_OUT = "dpo" + suffix  ## the most important
    else:
        OBJECTS_OUT = object_out

    FINISHFILE = "dfin" + suffix  ## remove this (for running dophot mulitple times)

    if image_in is None:
        IMAGE_IN = (
            "masked_se" + suffix
        )  # "resid" + suffix   ## can be resid + (extended objects + Dmask --> dophot)
    else:
        IMAGE_IN = image_in

    if not os.path.exists(root + IMAGE_IN):
        print("input image not found: ", IMAGE_IN)
        raise ReidualMaksedImageNotFound

    IMAGE_VAR = "model" + suffix
    # IMAGE_VAR   =  "dvar" + suffix  ## renuc factor 2x the one used for sextractor

    variance = modelName

    if renuc is not None and renuc != 1:

        IMAGE_VAR = "model" + suffix + "_do_renuc_" + str(renuc)
        variance = root + "/" + IMAGE_VAR

        script = (
            """
        rd 1 """
            + modelName  # prf file
            + """
        rd 2 """
            + objFits  #    xxxj.fits raw
            + """
        mc 1 """
            + str(renuc)
            + """

        ai 1 2
        wd 1 """
            + variance
            + """
        q
        
        """
        )

        run_monsta(script, root + "obj.pro", root + "obj.log")

    ## making dpar files
    dophot_params = (
        """
= """
        + name
        + """ GO-16262
SKY = """
        + SKY
        + """     Approximate mean sky value in data numbers.
FWHM = """
        + str(FWHM)
        + """             Approx FWHM of objects (pixels) along major axis.
AXIS_RATIO = 0.9       For star objects.  AR=b/a; b=minor axis.
TILT = 90.0            Angle of major axis in degrees; +x=0; +y=90.
EPERDN = 1.0           Electrons per data number.
RDNOISE = """
        + str(RDNOISE)
        + """         Readout noise in electrons.
TOP = """
        + TOP
        + """          Maximum allowed data value in data numbers.
=
AUTOSCALE = 'YES'      Auto-scaling of fitting and masking radii by FWHM.
AUTOSCALEAP = 'YES'    Auto-scaling of star and sky aperture radii by FWHM.
AUTOTHRESH = 'NO'      Auto-scaling of thresholds and bottom.
AUTOBLIT = 'YES'       Auto-scaling of obliteration parameters.
AUTOMEDSCALE = 'YES'   Auto-scaling for median filter box size.
= If autothresh=no take these parameters.
BOTTOM = """
        + BOTTOM
        + """       Lowest allowed data value in data numbers.
THRESHMAX = """
        + str(THRESHMAX)
        + """    Value of maximum threshold.
THRESHMIN = """
        + str(THRESHMIN)
        + """     Value of lowest threshold.
=
PARAMS_DEFAULT = 'paramdefault'   Default parameters file name.
IMAGE_IN = '"""
        + IMAGE_IN
        + """'            Input image name. 
=  If desired a file of extra variance can be specified, which is
=  added to the normal variance of (DATA/eADU)+(RN/eADU)^2.  Note
=  that it must have proper ADU^2 units of variance 
IMAGE_VAR = '"""
        + IMAGE_VAR
        + """'            Input extra variance image
IMAGE_OUT = '"""
        + IMAGE_OUT
        + """'             Output image name.
OBJECTS_OUT = '"""
        + OBJECTS_OUT
        + """'           Output object list file name.
CENTINTMAX = """
        + TOP
        + """              Obliterate if central intensity exceeds this.
TOPSAT = """
        + TOP
        + """                  Ignore pixels with higher counts in saturated stars.
=
OBJTYPE_IN = 'INTERNAL'     Input format: (COMPLETE, INTERNAL, BINARY)
OBJTYPE_OUT = 'COMPLETE'    Output format: (COMPLETE, INTERNAL,INCOMPLETE, BINARY,DAOPHOT)
=   If autoscaleap=yes determine parameters using
SCALEAPRADIUS = """
        + str(SCALEAPRADIUS)
        + """    Size of aperture phot radius in units of FWHM.
END
"""
    )

    IMAGE_OUT = root + IMAGE_OUT
    OBJECTS_OUT = root + OBJECTS_OUT
    DPAR_FILE = root + "dpar" + suffix
    with open(DPAR_FILE, "w") as f:
        f.write(dophot_params)

    ## run dophot
    cwd = os.getcwd()
    os.chdir(root)
    xcmd("cp " + configFolder + "/dophot/paramdefault .", verbose=False)
    xcmd("dophot dpar.000 2>&1 | tee dophot.log", verbose=False)
    xcmd("rm ./paramdefault", verbose=False)
    os.chdir(cwd)

    log = root + "/dophot.log"

    with open(log, "r") as f:
        lines = f.readlines()
        if "error" in "".join(lines).lower():
            print("Execution not completed. Please check the log file ... ")
            print(
                "\n Make sure that SExtractor_foreground has been already executed \n and the foreground masked residual is available."
            )
            print("Check the dpar file for some clues ...")
            return log, residName, DPAR_FILE, root + IMAGE_IN, None, None

    return log, residName, DPAR_FILE, root + IMAGE_IN, IMAGE_OUT, OBJECTS_OUT


##############################################################
def read_dpo_file(dpo_file, Xctr, Yctr, ZP=0):

    dpo_data = np.genfromtxt(dpo_file)

    df = pd.DataFrame(
        dpo_data[:, :10],
        columns=[
            "id",
            "star",
            "xpos",
            "ypos",
            "mag",
            "magerr",
            "xxx",
            "AA",
            "BB",
            "PA",
        ],
    )

    df["radius"] = ((Xctr - df.xpos) ** 2 + (Yctr - df.ypos) ** 2) ** 0.5
    df["magauto"] = df.mag + ZP
    df["tilt"] = np.radians(df.PA)
    df["fwhm"] = df.AA

    df["kron"] = df.fwhm

    df["cxx"] = (np.cos(df.tilt) ** 2 / df.AA ** 2) + (
        np.sin(df.tilt) ** 2 / df.BB ** 2
    )
    df["cyy"] = (np.sin(df.tilt) ** 2 / df.AA ** 2) + (
        np.cos(df.tilt) ** 2 / df.BB ** 2
    )
    df["cxy"] = (
        2.0 * np.cos(df.tilt) * np.sin(df.tilt) * (1.0 / df.AA ** 2 - 1.0 / df.BB ** 2)
    )

    # Additional correction at the faint end based on Cho et al. (2016) comparison
    df["magautocor"] = df.apply(
        lambda row: row.magauto - 0.038 * (row.magauto - 24)
        if row.magauto > 24
        else row.magauto,
        axis=1,
    )

    df["mtot"] = df.magautocor - 0.539
    df["mtoterr"] = df.magerr

    df["sestar"] = 1
    df["sestar"][((df.star == 2) & (df.magerr > -0.01) & (df.magerr < 0.3))] = 0
    df["sestar"][((df.star == 5) & (df.magerr > 0.0) & (df.magerr < 0.3))] = 0

    df["isgal"] = -1
    df["isgal"][((df.star == 3) & (df.magerr > 0.0) & (df.magerr < 0.3))] = 0
    df["isgal"][((df.star == 2) & (df.magerr > -0.01) & (df.magerr < 0.3))] = 1
    df["isgal"][((df.star == 7) & (df.magerr > 0.0) & (df.magerr < 0.3))] = 0
    df["isgal"][((df.star == 4) & (df.magerr > 0.0) & (df.magerr < 0.3))] = 0
    df["isgal"][((df.star == 4) & (df.magerr > 0.0) & (df.magerr < 0.3))] = 1

    return df


##############################################################
def make_do_lkn(
    dpo_file,
    model=0,
    r_aperture=0,
    exclude457=False,
    lkn_file_name=None,
    ZP=0,
    **Config
):

    name = Config["name"]
    objRoot = Config["objRoot"]
    Xctr = Config["X0"]
    Yctr = Config["Y0"]

    header = """
    # Region file format: DS9 version 4.1
    image
    """

    suffix = ""

    if not model is None:
        suffix = ".%03d" % model

    if lkn_file_name is None:
        lkn_file_name = objRoot + name + "_do_lkn" + suffix

    DF = read_dpo_file(dpo_file, Xctr, Yctr, ZP)

    df = DF.copy()

    if exclude457:
        df = df[((df.star != 4) & (df.star != 5) & (df.star != 7))]

    gc = len(df[df.isgal == 0])
    gal = len(df[df.isgal == 1])

    nfile = 0

    with open(lkn_file_name, "w") as alk:

        for i in range(len(df)):
            R = df.iloc[i]
            alk.write(
                "%d %8.2f %7.2f %8.2f  1 %9.5f %8.5f %9.5f  %5s %7.4f %d %8.4f %6.4f %d %d\n"
                % (
                    R.id,
                    R.xpos,
                    R.ypos,
                    R.radius,
                    R.cxx,
                    R.cyy,
                    R.cxy,
                    R.kron,
                    R.AA,
                    R.id,
                    R.mtot,
                    R.mtoterr,
                    R.sestar,
                    R.star,
                )
            )
            nfile += 1

    print("wrote: ", lkn_file_name)
    print("number of lines: ", nfile)
    print("GCs: ", gc)
    print("galaxies: ", gal)

    if exclude457:
        print("Excluding dophot object types 4, 5, and 7")
    else:
        print("Including dophot object types 4, 5, and 7")

    with open("ds9.reg", "w") as ds9_file:
        ds9_file.write(header)
        ds9_file.write(
            "circle(%.4f, %.4f, %4f) # color=yellow \n" % (Xctr, Yctr, r_aperture)
        )
        for i in range(len(df)):

            R = df.iloc[i]

            if R.radius > 20:
                if R.isgal == 1:
                    ds9_file.write(
                        "ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) # color=green \n"
                        % (
                            R.xpos + 0.5,
                            R.ypos + 0.5,
                            R.kron * R.AA,
                            R.kron * R.BB,
                            R.PA,
                        )
                    )
                elif R.isgal == 0:
                    ds9_file.write(
                        "circle(%.4f, %.4f, 3) # color=cyan \n"
                        % (R.xpos + 0.5, R.ypos + 0.5)
                    )

                else:
                    ds9_file.write(
                        "circle(%.4f, %.4f, 3) # color=red \n"
                        % (R.xpos + 0.5, R.ypos + 0.5)
                    )

    return DF, lkn_file_name

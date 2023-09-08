from .utils import *
import copy


##############################################################
def SE_foreground(
    model=0,
    smooth=None,
    minArea=10,
    mask=None,
    thresh=2,
    good_segments=[0],
    r_aperture=100,
    renuc=1,
    **Config
):

    root = Config["objRoot"]
    name = Config["name"]
    inFolder = Config["inFolder"]
    configFolder = Config["configFolder"]
    X0 = Config["X0"]
    Y0 = Config["Y0"]

    suffix = ".%03d" % model

    if mask is None:
        suffix_mask = ".%03d" % model
    else:
        suffix_mask = ".%03d" % mask

    residName = root + "/resid" + suffix
    modelName = root + "/model" + suffix
    segment = root + "/objCheck" + suffix + ".segment"
    objName = root + "/objCheck" + suffix
    objCatal = root + "/objCatal" + suffix
    maskName = root + "/mask" + suffix_mask
    sex_obj_maskName = root + "/mask_se" + suffix
    model_mask = root + "/mask" + ".%03d" % model
    tmp = root + "/smooth"

    objFits = inFolder + "{}/{}j.fits".format(name, name)
    Dmask = inFolder + "{}/{}j.dmask".format(name, name)

    maskName = Dmask

    script = (
        """
    rd 1 '"""
        + objFits  # original fits image
        + """'
    rd 2 '"""
        + modelName  # prf file
        + """'
    rd 3 """
        + maskName
        + """
    si 1 2
    mi 1 3
    wd 1 '"""
        + sex_obj_maskName
        + """'
    q
    """
    )
    run_monsta(script, root + "obj.pro", root + "obj.log")

    if smooth is not None:
        script = (
            """
        rd 1 """
            + sex_obj_maskName
            + """
        smooth 1 fw="""
            + str(smooth)
            + """
        wd 1 """
            + tmp
            + """
        q
        
        """
        )

        run_monsta(script, root + "obj.pro", root + "obj.log")
        sex_obj_maskName = tmp

    sex_configFolder = configFolder + "sextractor/wfc3j.inpar"
    PARAMETERS_NAME = configFolder + "sextractor/sbf.param"
    FILTER_NAME = configFolder + "sextractor/gauss_2.0_5x5.conv"
    STARNNW_NAME = configFolder + "sextractor/default.nnw"

    sex_cmd = (
        """sex """
        + sex_obj_maskName  # residual image
        + """ -c """
        + sex_configFolder
        + """ -CHECKIMAGE_NAME """
        + segment  ### "./back.fits"
    )
    sex_cmd += " -CATALOG_NAME  " + objCatal
    sex_cmd += " -DETECT_MINAREA " + str(minArea)
    sex_cmd += " -DETECT_THRESH " + str(thresh)
    sex_cmd += " -ANALYSIS_THRESH " + str(thresh)
    sex_cmd += " -BACK_SIZE 32"
    sex_cmd += " -CHECKIMAGE_TYPE SEGMENTATION "
    sex_cmd += (
        " -PARAMETERS_NAME "
        + PARAMETERS_NAME
        + " -FILTER_NAME "
        + FILTER_NAME
        + " -STARNNW_NAME "
        + STARNNW_NAME
    )

    variance = modelName

    if True:  # renuc is not None and renuc != 1:
        variance = root + "/model" + suffix + "_se_renuc_" + str(renuc)
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

    sex_cmd += " -WEIGHT_IMAGE  " + variance
    sex_cmd += " -WEIGHT_TYPE  MAP_VAR"

    # runnint Sextractor
    xcmd(sex_cmd + " > " + root + "sextractor.log", verbose=False)

    df = get_sextract_catal_df(objCatal)
    df["rc"] = np.sqrt((df.X_IMAGE - X0) ** 2 + (df.Y_IMAGE - Y0) ** 2)

    df_ignore = df[df.rc <= r_aperture]
    good_segments = df_ignore["NUMBER"].values
    df = df[df.rc > r_aperture]

    imarray = seg2mask(segment, objName, good_segments=good_segments, object_mask=True)

    ## Monsta script
    sex_obj_masked = root + "/masked_se" + suffix
    script = (
        """
    rd 1 """
        + residName
        + """
    rd 2 """
        + objName
        + """
    rd 5 '"""
        + configFolder
        + """/common.mask'
    mi 1 2
    mi 1 5
    wd 1 """
        + sex_obj_masked
        + """ 
    tv 1 JPEG="""
        + sex_obj_masked
        + """.jpg
    q
    
    """
    )
    # print(root+'obj'+suffix+'.pro', root+'obj'+suffix+'.log')

    run_monsta(script, root + "obj" + suffix + ".pro", root + "obj" + suffix + ".log")

    print(root + "obj" + suffix + ".pro")

    sex_obj_maskName = root + "/mask_se" + suffix

    return objCatal, df, objName, sex_obj_maskName, sex_obj_masked, residName, segment


##############################################################


def make_se_lkn(
    catal_df,
    model=None,
    star_f=0.7,
    r_aperture=0,
    filter="j",
    zp=0,
    Extended=False,
    sfx="se",
    **Config
):

    name = Config["name"]
    objRoot = Config["objRoot"]
    X0 = Config["X0"]
    Y0 = Config["Y0"]

    header = """
    # Region file format: DS9 version 4.1
    image
    """

    #### ---> Outputs_n0439/n0439_se.lknj
    gal = 0
    gc = 0
    nline = 0

    lkn_file_name = objRoot + name + "_" + sfx + "_lkn" + filter

    if not model is None:
        suffix = ".%03d" % model
        lkn_file_name = objRoot + name + "_" + sfx + "_lkn" + filter + suffix

    lkn_file = open(lkn_file_name, "w")

    Xctr = X0
    Yctr = Y0
    # I already put m1star in the sextractor .inpar files
    #    m1zpt = names[name]['m1star814']
    m1zpt = zp

    #    test.write('#  0  1         2         3     4      5          6       7       8       9      10      11       12\n')
    #    test.write('#  N  magaut   +/-       m4    +/-   maut_cor   m4cor   star     m5     +/-    m5cor  aut_cor  isoerr\n')
    lkn_file.write(
        "%4s %8.2f %7.2f %8.2f  0 %9.5f %9.5f %9.5f  %5s %7.4f %5s %8.4f %6.4f\n"
        % ("0", Xctr, Yctr, 0, 0, 0, 0, "0.00", 0, "0", 0, 0)
    )

    ingore_id_list = []

    id = catal_df.NUMBER.values
    xpos = catal_df.X_IMAGE.values - 0.5
    ypos = catal_df.Y_IMAGE.values - 0.5
    radius = ((Xctr - xpos) ** 2 + (Yctr - ypos) ** 2) ** 0.5
    magauto = catal_df.MAG_AUTO.values + m1zpt
    magautoerr = catal_df.MAGERR_AUTO.values
    AA = catal_df.A_IMAGE.values  # A_IMAGE
    BB = catal_df.B_IMAGE.values  # B_IMAGE
    fwhm = catal_df.FWHM_IMAGE.values
    star = catal_df.CLASS_STAR.values  # Class_star S/G (1=galaxy, 0=star)
    pa = catal_df.THETA_IMAGE.values

    cxx = catal_df.CXX_IMAGE.values
    cyy = catal_df.CYY_IMAGE.values
    cxy = catal_df.CXY_IMAGE.values
    kron = catal_df.KRON_RADIUS.values
    magiso = catal_df.MAG_ISO.values + m1zpt
    magisoerr = catal_df.MAGERR_ISO.values

    m4 = catal_df.MAG_APER_4.values + m1zpt
    m5 = catal_df.MAG_APER_5.values + m1zpt
    m4err = catal_df.MAGERR_APER_4.values
    m5err = catal_df.MAGERR_APER_5.values

    with open("ds9.reg", "w") as ds9_file:
        ds9_file.write(header)

        if r_aperture > 0:
            ds9_file.write(
                "circle(%.4f, %.4f, %4f) # color=red \n" % (Xctr, Yctr, r_aperture)
            )

        for i in range(len(id)):

            ignore = False

            #            aut_cor = max(0.05, 0.1 + 0.01*(magauto-19.5))
            # We haven't determined a correction for auto mags for WFC3/IR yet
            aut_cor = 0

            # From Cho et al. 2016
            magautocor = magauto[i] - aut_cor
            m4cor = (
                m4[i] - 0.259
            )  # aperture correction particular filter F110W, particular instrument WFC3-IR
            m5cor = m5[i] - 0.259

            isStar = False
            if star[i] >= star_f and m5err[i] < 50.0 and magautoerr[i] < 50.0:
                mtot = m5cor
                mtoterr = m5err[i]
                gc = gc + 1
                isStar = True
            else:
                mtot = magautocor
                mtoterr = magautoerr[i]
                gal = gal + 1

            merr = min(mtoterr, magisoerr[i])

            # already extinction corrected
            # add check to make sure it's a reasonable & useful object:
            # changed merr>50 to merr>0.3 to produce a cleaner sample.
            # Ignore objects within 10 pix of the galaxy center.

            if merr > 0.3 or radius[i] < 10:
                ignore = True

            # Add a check to remove the bad objects with FWHM=0 that SE sometimes finds
            # These single pixel object might be bad pixels or cosmic rays
            # in masked areas of WFC3/IR data;
            # reduced fwhm limit to 20 pix to avoid huge residual "galaxies" (maybe only needed for Coma)

            if fwhm[i] < 0.01 or fwhm[i] > 10:
                ignore = True

            if not ignore:
                if isStar:
                    ds9_file.write(
                        "circle(%.4f, %.4f, 3) # color=cyan \n"
                        % (xpos[i] + 0.5, ypos[i] + 0.5)
                    )
                else:
                    ds9_file.write(
                        "ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) # color=green \n"
                        % (
                            xpos[i] + 0.5,
                            ypos[i] + 0.5,
                            kron[i] * AA[i],
                            kron[i] * BB[i],
                            pa[i],
                        )
                    )

                if not isStar or not Extended:
                    if radius[i] != 0:
                        lkn_file.write(
                            "%5s %8.2f %7.2f %8.2f  1 %9.5f %8.5f %9.5f  %5s %7.4f %5s %8.4f %6.4f %6.3f %6.4f\n"
                            % (
                                id[i],
                                xpos[i],
                                ypos[i],
                                radius[i],
                                cxx[i],
                                cyy[i],
                                cxy[i],
                                kron[i],
                                AA[i],
                                id[i],
                                mtot,
                                mtoterr,
                                star[i],
                                merr,
                            )
                        )
                        nline += 1

            else:
                ds9_file.write(
                    "ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) # color=yellow \n"
                    % (
                        xpos[i] + 0.5,
                        ypos[i] + 0.5,
                        kron[i] * AA[i],
                        kron[i] * BB[i],
                        pa[i],
                    )
                )
                ingore_id_list.append(id[i])

    lkn_file.close()

    print("wrote: ", lkn_file_name)
    print("number of lines: ", nline)
    print("# of GCs: ", gc)
    print("# of galaxies: ", gal)
    gal = 0
    gc = 0

    return lkn_file_name, ingore_id_list


##############################################################


def edit_ptm6_mask(like_pack, rad=5, self=False, **Config):

    objRoot = Config["objRoot"]
    lkn, lkn6, ptm6, resid = like_pack

    if self:
        infile = "{}b".format(ptm6)
    else:
        infile = ptm6

    cwd = os.getcwd()
    os.chdir(objRoot)

    if not os.path.exists(infile):
        xcmd("cp {} {}".format(ptm6, infile), verbose=True)

    monsta_script = (
        """
        rd 1 """
        + resid
        + """
        rd 2 '"""
        + infile
        + """
        mi 2 1 
        tv 2 -1000 1000
        tidy 2 defer rad="""
        + str(rad)
        + """
        di 2 2 
        wd 2 """
        + ptm6
        + """b bitmap
    """
    )
    run_monsta(monsta_script, "monsta.pro", "monsta.log")
    now = datetime.now().strftime("%Y%m%d%H%M%S")
    xcmd("cp {}b {}b.{}".format(ptm6, ptm6, now), verbose=False)
    print("modified mask name: {}b".format(ptm6))
    print("backup: {}b.{}".format(ptm6, now))
    os.chdir(cwd)

    return "{}b".format(ptm6)


##############################################################
def resid_ptm6(resid, all_masks, outfile=None, **Config):

    objRoot = Config["objRoot"]
    name = Config["name"]

    if outfile is None:
        outfile = "{}_resid_masked".format(name)

    cwd = os.getcwd()
    os.chdir(objRoot)
    monsta_script = (
        """
        rd 1 """
        + resid
    )

    for mask in all_masks:
        monsta_script += (
            """
        rd 2 """
            + mask
            + """
        mi 1 2
        """
        )

    monsta_script += """wd 1 """ + outfile

    run_monsta(monsta_script, "monsta.pro", "monsta.log")
    print("residual file ready for dophot: ", outfile)
    os.chdir(cwd)

    return outfile


##############################################################
def display_ptm6_mask(like_pack, **Config):

    objRoot = Config["objRoot"]
    lkn, lkn6, ptm6, resid = like_pack

    cwd = os.getcwd()
    os.chdir(objRoot)
    monsta_script = (
        """
        rd 1 """
        + resid
        + """
        rd 2 '"""
        + ptm6
        + """
        mi 2 1 
        rd 3 '"""
        + ptm6
        + """b
        mi 3 1
        tv 1 -1000 1000
        tv 2 -1000 1000
        tv 3 -1000 1000
        tidy 1 defer rad=5
    """
    )
    run_monsta(monsta_script, "monsta.pro", "monsta.log")
    os.chdir(cwd)


##############################################################
def ds9_scale_image(reg_file="./ds9.reg"):

    os.system("xpaset -p ds9 regions delete all &")
    os.system("xpaset -p ds9 regions " + reg_file + " &")

    # http://ds9.si.edu/doc/ref/xpa.html
    os.system("xpaset -p ds9 regions " + reg_file + " &")
    os.system("xpaset -p ds9 cmap value 7.9 0.47 &")

    # set image scale
    os.system("xpaset -p ds9 scale log &")
    os.system("xpaset -p ds9 scale limits -4100 200000 &")


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################

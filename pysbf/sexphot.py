from .utils import *



def unnamed(catal_df, **Config):

    name = Config["name"]
    objRoot = Config["objRoot"]
    X0 = Config["X0"]
    Y0 = Config["Y0"]

    header = '''
    # Region file format: DS9 version 4.1
    image
    '''

    #### ---> Outputs_n0439/n0439_se.lknj

    filter = 'j'
    gal = 0 
    gc = 0


    alk = open(objRoot+name+'_se.lkn'+filter, 'w')

    Xctr = X0
    Yctr = Y0
    # I already put m1star in the sextractor .inpar files
    #    m1zpt = names[name]['m1star814']
    m1zpt = 0.0

    #    test.write('#  0  1         2         3     4      5          6       7       8       9      10      11       12\n')
    #    test.write('#  N  magaut   +/-       m4    +/-   maut_cor   m4cor   star     m5     +/-    m5cor  aut_cor  isoerr\n')
    alk.write('%4s %8.2f %7.2f %8.2f  0 %9.5f %9.5f %9.5f  %5s %7.4f %5s %8.4f %6.4f\n' \
            %('0', Xctr, Yctr, 0, 0,0,0, '0.00', 0, '0', 0, 0))

    xpos_lst = []
    ypos_lst = []

    id = catal_df.NUMBER.values
    xpos = catal_df.X_IMAGE.values - 0.5
    ypos = catal_df.Y_IMAGE.values - 0.5
    radius = ((Xctr-xpos)**2 + (Yctr-ypos)**2)**0.5 
    magauto = catal_df.MAG_AUTO.values + m1zpt
    magautoerr = catal_df.MAGERR_AUTO.values
    AA = catal_df.A_IMAGE.values       # A_IMAGE  
    BB = catal_df.B_IMAGE.values       # B_IMAGE 
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

    with open('ds9.reg', 'w') as f:
        f.write(header)
        for i in range(len(catal_df)):       

        #            aut_cor = max(0.05, 0.1 + 0.01*(magauto-19.5))
        # We haven't determined a correction for auto mags for WFC3/IR yet
                aut_cor = 0

        # From Cho et al. 2016
                magautocor = magauto[i] - aut_cor
                m4cor = m4[i] - 0.259  # aperture correction particular filter F110W, particular instrument WFC3-IR
                m5cor = m5[i] - 0.259
                
                isStar = False
                if star[i] >= 0.7 and m5err[i] < 50. and magautoerr[i] < 50. :
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
                
                # if merr>0.3 or radius[i]<10:
                #     continue

                # Add a check to remove the bad objects with FWHM=0 that SE sometimes finds 
                # These single pixel object might be bad pixels or cosmic rays
                # in masked areas of WFC3/IR data;
                # reduced fwhm limit to 20 pix to avoid huge residual "galaxies" (maybe only needed for Coma)
                
        #         if fwhm[i]<0.01 or fwhm[i]>10:
        # #                gal = gal - 1
        #             continue
            
                if isStar:
                    f.write("circle(%.4f, %.4f, 3) # color=red \n"%(xpos[i]+0.5,ypos[i]+0.5))
                else:
                    f.write("ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) # color=green \n"%(xpos[i]+0.5,ypos[i]+0.5,kron[i]*AA[i], kron[i]*BB[i], pa[i]))
                    

                alk.write('%5s %8.2f %7.2f %8.2f  1 %9.5f %8.5f %9.5f  %5s %7.4f %5s %8.4f %6.4f %6.3f %6.4f\n' \
                            %(id[i], xpos[i], ypos[i], radius[i], cxx[i], cyy[i], cxy[i], kron[i], AA[i], id[i], mtot, mtoterr, star[i], merr))

                xpos_lst.append(xpos[i])
                ypos_lst.append(ypos[i])
                
    alk.close()

    print('wrote: ', objRoot+name+'_se.lkn'+filter)
    print('GCs: ',gc)
    print('galaxies: ',gal)
    gal=0
    gc=0

    string name '/media/Data/Home/PanStarrs/Jan/HI/augment/SBF/codes/n0679/ps/n0679'
    rd 1 '{name}g.fits'
    rd 2 '{name}z.fits'
    clip 1 nan=0
    clip 2 nan=0
    !tv 1
    set nr=100/10$nint
    string nr '%i2.0' nr
    cop 3 1
    cop 4 2
    tv 3 zoom=2
    elliprof 3 tv model rmstar x0=1500 y0=1500 r0=7 r1=100 nr=nr niter=5
    tv 4 zoom=2
    elliprof 4 tv model rmstar x0=1500 y0=1500 r0=7 r1=100 nr=nr niter=5
    cop 5 1
    si 5 3
    cop 6 2
    si 6 4
    wd 5 '{name}g.resid'
    wd 6 '{name}z.resid'

    ! run SExtractor
    % sex {name}g.resid -c config/sextractor/psg.inpar -CATALOG_NAME Outputs_n0679//psg.cat -CHECKIMAGE_NAME Outputs_n0679//psg.obj
    % sex {name}z.resid -c config/sextractor/psz.inpar -CATALOG_NAME Outputs_n0679//psz.cat -CHECKIMAGE_NAME Outputs_n0679//psz.obj

    rd 7 Outputs_n0679//psg.obj
    rd 8 Outputs_n0679//psz.obj
    !tv 7
    !tv 8
    di 7 7
    di 8 8
    mi 8 7
    wd 8 Outputs_n0679/n0679.PS.mask bitmap

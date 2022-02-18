
        string name 'n0679'
        rd 1 '/media/Data/Home/PanStarrs/Jan/HI/augment/SBF/codes//n0679/n0679j.fits'
        sc 1 3098.085864257813                            ! sky subtraction
        rd 2 Outputs_n0679//mask.001
        mi 1 2
        tv 1 sqrt JPEG=Outputs_n0679//n0679.000.jpg

        
        cop 6 2
        mc 6 0
        ac 6 1
        si 6 2
        rd 7 Outputs_n0679//model.000
        sc 7 3098.085864257813   ! aky subtraction
        mi 7 6
        cop 3 1 
        ai 3 7 
        elliprof 3  model rmstar x0=564.689 y0=562.985 r0=9.0 r1=260 nr=17 niter=10  COS3X=0
        print elliprof file=Outputs_n0679//elliprof.000
        cop 4 1                               ! object
        si 4 3                                ! object - model
        ac 3 3098.085864257813                  
        !mi 3 2 
        mi 4 2                                ! multiply by mask
        wd 3 Outputs_n0679//model.000
        wd 4 Outputs_n0679//resid.000
        tv 4 JPEG=Outputs_n0679//resid.000.jpg
        tv 3 JPEG=Outputs_n0679//model.000.jpg
        q
        
        
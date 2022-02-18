
        string name 'n3268'
        rd 1 '/media/Data/Home/PanStarrs/Jan/HI/augment/SBF/wfc3-16262//n3268/n3268j.fits'
        sc 1 4144.910888671875                            ! sky subtraction
        rd 2 Outputs_n3268//mask.001
        mi 1 2
        tv 1 sqrt JPEG=Outputs_n3268//n3268.000.jpg

        
        cop 6 2
        mc 6 0
        ac 6 1
        si 6 2
        rd 7 Outputs_n3268//model.000
        sc 7 4144.910888671875   ! aky subtraction
        mi 7 6
        cop 3 1 
        ai 3 7 
        elliprof 3  model rmstar x0=564.826 y0=566.165 r0=9.0 r1=279 nr=19 niter=10  COS3X=0
        print elliprof file=Outputs_n3268//elliprof.000
        cop 4 1                               ! object
        si 4 3                                ! object - model
        ac 3 4144.910888671875                  
        !mi 3 2 
        mi 4 2                                ! multiply by mask
        wd 3 Outputs_n3268//model.000
        wd 4 Outputs_n3268//resid.000
        tv 4 JPEG=Outputs_n3268//resid.000.jpg
        tv 3 JPEG=Outputs_n3268//model.000.jpg
        q
        
        
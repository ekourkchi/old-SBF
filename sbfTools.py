import os, sys, string, math
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

import copy



##############################################################
def set_axes(ax, xlim=None, ylim=None, fontsize=16, twinx=True, twiny=True, minor=True, inout='in'):
        
        if not ylim is None:
            ax.set_ylim(ylim)
        else:
            ylim = ax.get_ylim() 
            
        if not xlim is None:    
            ax.set_xlim(xlim) 
        else:
            xlim = ax.get_xlim()
            
        ax.tick_params(which='major', length=6, width=1., direction=inout)
#         if minor:
        ax.tick_params(which='minor', length=3, color='#000033', width=1.0, direction=inout)  
        
        if twiny:
            y_ax = ax.twinx()
            y_ax.set_ylim(ylim)
            y_ax.set_yticklabels([])
            y_ax.minorticks_on()
            y_ax.tick_params(which='major', length=6, width=1., direction=inout)
            if minor:
                y_ax.tick_params(which='minor', length=3, color='#000033', width=1.0, direction=inout) 
        
        if twinx:
            x_ax = ax.twiny()
            x_ax.set_xlim(xlim)
            x_ax.set_xticklabels([])
            x_ax.minorticks_on()
            x_ax.tick_params(which='major', length=6, width=1.0, direction=inout)
            if minor:
                x_ax.tick_params(which='minor', length=3, color='#000033', width=1.0, direction=inout)     

        for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(fontsize) 
        for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(fontsize) 
        
        return x_ax, y_ax

##############################################################

def ellipse_polyline(ellipses, n=10000):
    t = np.linspace(0, 2*np.pi, n, endpoint=False)
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
            return (mp.type=="MultiPoint")
        except:
            print((ell_a, ell_b))
            return True   # likely invalid
    else:
        return False

##############################################################

## returns True if the surface of two ellipoids cross each other
def ellipsoid_intersection_test(Ell_A, Ell_B, tau=1):
    
    Sigma_A, mu_A = Ell_A[4], Ell_A[0]
    Sigma_B, mu_B = Ell_B[4], Ell_B[0]
    
    lambdas, Phi, v_squared = ellipsoid_intersection_test_helper(Sigma_A, Sigma_B, mu_A, mu_B)
    res = minimize_scalar(ellipsoid_K_function,
                          bracket=[0.0, 0.5, 1.0],
                          args=(lambdas, v_squared, tau))
    return (res.fun[0] >= 0)


def ellipsoid_intersection_test_helper(Sigma_A, Sigma_B, mu_A, mu_B):
    lambdas, Phi = eigh(Sigma_A, b=Sigma_B)
    v_squared = np.dot(Phi.T, mu_A - mu_B) ** 2
    return lambdas, Phi, v_squared


def ellipsoid_K_function(ss, lambdas, v_squared, tau):
    ss = np.array(ss).reshape((-1,1))
    lambdas = np.array(lambdas).reshape((1,-1))
    v_squared = np.array(v_squared).reshape((1,-1))
    return 1.-(1./tau**2)*np.sum(v_squared*((ss*(1.-ss))/(1.+ss*(lambdas-1.))), axis=1)

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

    if verbose: print('\n'+cmd)

    tmp=os.popen(cmd)
    output=''
    for x in tmp: output+=x
    if 'abort' in output:
        failure=True
    else:
        failure=tmp.close()
    if False:
        print('execution of %s failed' % cmd)
        print('error is as follows', output)
        sys.exit()
    else:
        return output

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
    
    theta = angle * np.pi / 180.
    ux = np.cos(theta)
    uy = np.sin(theta)
    
    return (ux, uy)

def get_Sigma(Smajor, Sminor, angle):

    u1 = unit_vecort(angle)
    u2 = unit_vecort(angle-90)
    R = np.vstack((u1,u2))
    L = R.T
    X = np.diag((1./Smajor**2, 1./Sminor**2))
    A = np.matmul(np.matmul(L, X), R)

    return np.linalg.inv(A)

##############################################################

def make_Ellipse(center, Smajor, Sminor, angle):
    
    Sigma = get_Sigma(Smajor, Sminor, angle)
    return [np.asarray(center), Smajor, Sminor, angle, Sigma]

def plot_Ellipse(Ell, ax=None, **kwargs):
    
    if ax is None:
        ax = plt.gca()
    
    e = patches.Ellipse(Ell[0], width=2*Ell[1], height=2*Ell[2], angle=Ell[3], **kwargs)
    ax.add_patch(e)
    
    return ax

##############################################################

def list_Ell(data_frame):
    
    Ell_list = []
    
    for i in range(len(data_frame)):
        
        df = data_frame.iloc[i]
        Smajor = df.Rmaj
        Sminor = Smajor*(1.-df.ellip)
        angle = df.alpha+90
        center = (df.x0, df.y0)
        Sigma = get_Sigma(Smajor, Sminor, angle)
        Ell = make_Ellipse(center, Smajor, Sminor, angle)
        
        Ell_list.append(Ell)
    
    return Ell_list

##############################################################

# df is a Pandas series that holds one ellipse
def plot_E(df, ax=None, **kwargs):
    
    Smajor = df.Rmaj
    Sminor = Smajor*(1.-df.ellip)
    
    angle = df.alpha+90
    center = (df.x0, df.y0)
    Sigma = get_Sigma(Smajor, Sminor, angle)
    Ell = make_Ellipse(center, Smajor, Sminor, angle)

    if ax is None:
        ax = plt.gca()
        
    plot_Ellipse(Ell, ax=ax, **kwargs)
    
    return ax

##############################################################

def tv(fits_file, ax=None, options=""):

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
    x_max, y_max, _ = img.shape
        
    if ax is None:
        plt.figure(figsize=(10,10))
        plt.subplot(111)
        ax = plt.gca()

    ax.set_xlim([0, x_max])
    ax.set_ylim([0, y_max])

    imgplot = ax.imshow(np.flipud(img))

    return ax

##############################################################

def run_monsta(script, Monsta_pro, Monsta_log, monsta="/home/ehsan/Home/Monsta/bin/monsta"):
    
    with open(Monsta_pro, 'w') as f:
        f.write(script)

    cmd = monsta+' '+Monsta_pro
    xcmd(cmd + ' > '+Monsta_log, verbose=False)

    with open(Monsta_log) as f:
        text = f.read()

    if "ERROR" in text.upper().split():
        print(text)
        return Monsta_log
    return 'OK'
            
##############################################################

def imOpen(inFits):
    
    hdu_list = fits.open(inFits)
    imarray = hdu_list[0].data
    header = hdu_list[0].header   
    
    return imarray, header

##############################################################

def seg2mask(inFits, outMask, overwrite=True, seg_num=0):
    
    imarray, header = imOpen(inFits)
    imarray[((imarray==0))] = -1    # good pixel
    if seg_num!=0:
         imarray[((imarray==seg_num))] = -1    # good pixel
    imarray[(imarray!=-1)] = 0   # masked
    imarray[imarray==-1] = 1   # good pixel
    
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
    for n in range(0, N-1):
        ell_a = ells[n]
        for m in range(n+1, N):
            ell_b = ells[m]
            cross = ell_itersect(ell_a, ell_b)
            if cross:
                n_cross+=1
    return n_cross    

##############################################################


##############################################################
class ellOBJ:
    
    x0 = 0
    y0 = 0
    a  = 1
    b  = 1
    name = ""
    R = 1        # Kron Radius
    angle = 0    # Position angle (CCW/x)
    
    catalName = "catal.dat"
    objRoot = './'
    inFolder = './'
    monsta = 'monsta'
    config = './'

    
    def __init__(self, name, outFolder=None, inFolder=None, configFolder='./', automatic=True):
        
        if outFolder is None:
            outFolder = "Outputs_"+name

        
        createDir(outFolder)
        self.objRoot = outFolder+'/'
        

        if inFolder is not None:
            self.inFolder = inFolder+'/'

        self.config = config
        
        self.name = name
        
        if automatic:
            try:
                self.SExtract()
            except:
                print("Error: Could not run SExtractor on the file")
                fits_file = self.inFolder+'{}/{}j.fits'.format(name,name)
                if not os.path.exists(fits_file):
                    print("Couldn't find "+fits_file)
                return
            
            hdu_list = fits.open(self.inFolder+'{}/{}j.fits'.format(name,name))
            image_data = hdu_list[0].data
            # w = wcs.WCS(hdu_list[0].header)
            self.x_max, self.y_max = image_data.shape
                    
            self.backSextract()

            self.r_max = int(min([self.x0, self.x_max-self.x0, self.y0, self.y_max-self.y0]))


            
        
    def tv_resid(self, model=0, ax=None, options="", additions=""):
        root = self.objRoot
        suffix = '.%03d'%model
        fits_file = root+'/resid'+suffix
        return self.tv(fits_file=fits_file, ax=ax, options=options, additions=additions)

    def tv_model(self, model=0, ax=None, options="", additions=""):
        root = self.objRoot
        suffix = '.%03d'%model
        fits_file = root+'/model'+suffix
        return self.tv(fits_file=fits_file, ax=ax, options=options, additions=additions)
    
    def tv_mask(self, mask=None, ax=None, options="", additions=""):
        root = self.objRoot

        if mask is None:
            mask_name='./common.mask'
        else:
            suffix = '.%03d'%mask
            mask_name = root+'/mask'+suffix

        return self.tv(fits_file=mask_name, ax=ax, options=options, additions=additions)
    
    def tv(self, fits_file=None, ax=None, options="", additions=""):
        
        if fits_file is None:
            name = self.name
            fits_file = self.inFolder+"{}/{}j.fits".format(name, name)
        
        root = self.objRoot
        jpg_name = root+'tv.jpg'
        
        ## Monsta script
        script = """
        rd 1 '"""+fits_file+"""'
        """+additions+"""
        tv 1 """+options+""" JPEG="""+jpg_name+"""
        q
        
        """

        self.run_monsta(script, root+'tv.pro', root+'tv.log')
        
        return self.plot_jpg(jpg_name, ax=ax)
        
        
        
    
    def SExtract(self):
        
        name = self.name
        root = self.objRoot
        config = self.config
        
        catalName = root + self.catalName
        segmentation = root + 'segmentation.fits'

        sex_config = config+'sextractor/wfc3j_sex.config'
        PARAMETERS_NAME = config+'sextractor/sbf.param'
        FILTER_NAME = config+'sextractor/gauss_2.0_5x5.conv'
        STARNNW_NAME = config+'sextractor/default.nnw'

       
        cmd = 'sex -c ' + sex_config + ' '+self.inFolder+'{}/{}j.fits'.format(name,name)+' -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME '+segmentation+' -CATALOG_NAME ' + catalName + ' -PARAMETERS_NAME ' + PARAMETERS_NAME+ ' -FILTER_NAME ' + FILTER_NAME + ' -STARNNW_NAME ' + STARNNW_NAME
        
        # print(cmd)
        
        xcmd(cmd + ' > '+root+'sextractor.log', verbose=False)
        
        
        col_names = self.getColName(catalName)

        df = pd.read_csv(catalName, delimiter=r"\s+", skiprows=len(col_names), 
                         header = None, names = col_names[:18], usecols = range(18))
        
        self.x0 = df.loc[0].X_IMAGE
        self.y0 = df.loc[0].Y_IMAGE
        A  = df.loc[0].A_IMAGE
        B  = df.loc[0].B_IMAGE
        R  = df.loc[0].KRON_RADIUS
        self.angle  = df.loc[0].THETA_IMAGE
        self.a = A*R
        self.b = B*R
        self.R = R
        self.name = name

    def inputMaks(self, maskName, mask=None):

        root = self.objRoot
        if mask is None:
            outMask = root+'/composit.mask'
        else:
            outMask = root+'/mask'+'.%03d'%mask

        script = """
        rd 1 ./common.mask
        """

        if maskName is not None:
            if os.path.isfile(maskName):
                script += """
                rd 2 """+maskName+"""
                mi 1 2
                """
        script += """
        wd 1 """+outMask+""" 
        q
        
        """       
        # print(root+'obj'+suffix+'.pro', root+'obj'+suffix+'.log')

        self.run_monsta(script, root+'obj_inMask.pro', root+'obj_inMask.log')

        return outMask       
        


    def addMasks(self, maskList=None, mask=None):
        
        root = self.objRoot
        if mask is None:
            outMask = root+'/composit.mask'
        else:
            outMask = root+'/mask'+'.%03d'%mask

        
        script = """
        rd 1 ./common.mask
        """

        if maskList is not None:

            for m in maskList:
                suffix = '.%03d'%m
                inMask = root+'/mask'+suffix
                if os.path.isfile(inMask):
                    script += """
                    rd 2 """+inMask+"""
                    mi 1 2
                    """
        script += """
        wd 1 """+outMask+""" 
        q
        
        """       
        # print(root+'obj'+suffix+'.pro', root+'obj'+suffix+'.log')

        self.run_monsta(script, root+'obj'+suffix+'.pro', root+'obj'+suffix+'.log')


        return outMask


    def outerR(self, c_kron):
        return min(int(c_kron*np.sqrt(self.a*self.b)), self.r_max)

        
    def maskOpen(self, mask=None, maskFile=None):
        root = self.objRoot

        if maskFile is None:
            if mask is None:
                inMask='./common.mask'
            else:
                suffix = '.%03d'%mask
                inMask = root+'/mask'+suffix
        else:
            inMask = maskFile
        
        ## Monsta script
        script = """
        rd 1 """+inMask+"""
        wd 1 """+inMask+'.fits'+""" 
        q

        """       

        run_monsta(script, 'monsta.pro', 'monsta.log')
        xcmd("rm monsta.log; rm monsta.pro", False)

        image, header = imOpen(inMask+'.fits')
        xcmd("rm "+inMask+'.fits', False)

        return image, header


    def elliprof(self, inner_r=5, outer_r=200, sky=None, 
                 cosnx="", k=None, 
                 nr=40, niter=10,
                 model = 0, mask=None, options="", use_common=False, model_mask=None
                ):
        
        root = self.objRoot
        suffix = '.%03d'%model
        name = self.name
        
        if sky is None:
            sky = self.sky_med
        
        if mask is None:
            maskName = './common.mask'
            monsta_masking =  """
        cop 3 1 
        """
        elif model_mask is None:
            maskName = root+'/mask'+'.%03d'%mask
            monsta_masking =  """
        cop 3 1 
        """
        else:
            maskName = root+'/mask'+'.%03d'%mask
            model_mask = root+'/model'+'.%03d'%model_mask
            monsta_masking =  """
        cop 6 2
        mc 6 0
        ac 6 1
        si 6 2
        rd 7 """+model_mask+"""
        sc 7 """+str(sky)+"""   ! aky subtraction
        mi 7 6
        cop 3 1 
        ai 3 7 
        """
              
        if cosnx == "":
            kosnx = ""
        else:
            if k is None:
                kosnx = cosnx+"=0"
            else:
                kosnx = cosnx+"="+str(k)
        
        
        residName = root+'/resid'+suffix
        modelName = root+'/model'+suffix
        ellipseFile = root+'/elliprof'+suffix
        objName = root+'/'+self.name+suffix
                                                                    
        elliprof_cmd = "elliprof 3  model rmstar x0="+str(self.x0)+" y0="+str(self.y0)
        elliprof_cmd += " r0="+str(inner_r)+" r1="+str(outer_r)+" nr="+str(nr)+" niter="+str(niter)+" "+kosnx
        elliprof_cmd += " "+options

       
        objFits = self.inFolder+'{}/{}j.fits'.format(name,name)
        
        ## Monsta script
        script = """
        string name '"""+self.name+"""'
        rd 1 '"""+objFits+"""'
        sc 1 """+str(sky)+"""                            ! sky subtraction
        rd 2 """+maskName+"""
        mi 1 2
        tv 1 sqrt JPEG="""+objName+""".jpg

        """

        script += monsta_masking

        script += elliprof_cmd
        script += """
        print elliprof file="""+ellipseFile+"""
        cop 4 1                               ! object
        si 4 3                                ! object - model
        ac 3 """+str(sky)+"""                  
        !mi 3 2 
        mi 4 2                                ! multiply by mask
        wd 3 """+modelName+"""
        wd 4 """+residName+"""
        tv 4 JPEG="""+residName+""".jpg
        tv 3 JPEG="""+modelName+""".jpg
        q
        
        """
       
        Monsta_pro = root+'monsta'+suffix+'.pro'
        Monsta_log = root+'monsta'+suffix+'.log'

        # print(Monsta_pro, Monsta_log)
        
        return self.run_monsta(script, Monsta_pro, Monsta_log)
        
    def objSExtract(self, model=0, smooth=None, minArea=10, thresh=2, mask=None, renuc=1):
        
        root = self.objRoot
        suffix = '.%03d'%model
        
        if mask is None:
            suffix_mask = '.%03d'%model
        else:
            suffix_mask = '.%03d'%mask
        
        residName = root+'/resid'+suffix
        modelName = root+'/model'+suffix
        segment = root+'/objCheck'+suffix+'.segment'
        objName = root+'/objCheck'+suffix
        objCatal = root+'/objCatal'+suffix
        maskName = root+'/mask'+suffix_mask
        tmp = root+'/tmp'

        name = self.name
        objFits = self.inFolder+'{}/{}j.fits'.format(name,name)
        script = """
        rd 1 '"""+objFits+"""'
        rd 2 """+maskName+"""
        rd 3 """+modelName+"""
        si 1 3
        mi 1 2
        wd 1 '"""+residName+"""'
        q
        """
        self.run_monsta(script, root+'obj.pro', root+'obj.log')

        print(root+'obj.pro')

        
        if smooth is not None:
            script = """
            rd 1 """+residName+"""
            smooth 1 fw="""+str(smooth)+"""
            wd 1 """+tmp+"""
            q
            
            """
            residName = tmp
            self.run_monsta(script, root+'obj.pro', root+'obj.log')

        variance = modelName
        if renuc is not None and renuc!=1:
            variance = root+'/model'+suffix+'_renuc_'+str(renuc)
            script = """
            rd 1 """+modelName+"""
            mc 1 """+str(renuc)+"""
            wd 1 """+variance+"""
            q
            
            """
            
            self.run_monsta(script, root+'obj.pro', root+'obj.log')    

        config = self.config

        sex_config = config+'sextractor/wfc3j.inpar'
        PARAMETERS_NAME = config+'sextractor/sbf.param'
        FILTER_NAME = config+'sextractor/gauss_2.0_5x5.conv'
        STARNNW_NAME = config+'sextractor/default.nnw'        
        
        sex_cmd = """sex """+residName+""" -c w""" + sex_config + """ -CHECKIMAGE_NAME """+segment
        sex_cmd += " -CATALOG_NAME  "+objCatal
        sex_cmd += " -DETECT_MINAREA " +str(minArea)
        sex_cmd += " -DETECT_THRESH "+str(thresh)
        sex_cmd += " -ANALYSIS_THRESH "+str(thresh)
        sex_cmd += " -CHECKIMAGE_TYPE SEGMENTATION "
        sex_cmd += ' -PARAMETERS_NAME ' + PARAMETERS_NAME+ ' -FILTER_NAME ' + FILTER_NAME + FILTER_NAME + ' -STARNNW_NAME ' + STARNNW_NAME

        if renuc is not None:
            sex_cmd += " -WEIGHT_IMAGE  "+variance
            sex_cmd += " -WEIGHT_TYPE  MAP_VAR"
            

        print(sex_cmd)

        xcmd(sex_cmd + ' > '+root+'sextractor.log', verbose=False)
        seg2mask(segment, objName) 

        print(segment, objName, maskName)
       
        
        # ## Monsta script
        # script = """
        # rd 1 """+objName+"""
        # rd 5 './common.mask'
        # mi 1 5
        # wd 1 """+maskName+""" bitmap
        # tv 1 JPEG="""+maskName+""".jpg
        # q
        
        # """       
        # # print(root+'obj'+suffix+'.pro', root+'obj'+suffix+'.log')

        # self.run_monsta(script, root+'obj'+suffix+'.pro', root+'obj'+suffix+'.log')

    def naive_Sextract(self, minArea=10, thresh=2, mask=None, smooth=None):
        
        name = self.name
        root = self.objRoot
        fits_name = self.inFolder+'{}/{}j.fits'.format(name,name)
        odj_common = root+'/tmp'

        if smooth is not None:
            script = """
            rd 1 """+fits_name+"""
            smooth 1 fw="""+str(smooth)+"""
            wd 1 """+odj_common+"""
            q
            
            """
            fits_name = odj_common
            self.run_monsta(script, root+'obj.pro', root+'obj.log')


        if mask is None:
            suffix_mask = '.%03d'%model
        else:
            suffix_mask = '.%03d'%mask

        maskName = root+'/mask'+suffix_mask
        
        script = """
        rd 1 """+fits_name+"""
        rd 2 ./common.mask
        mi 1 2
        wd 1 """+odj_common+"""
        q

        """
        self.run_monsta(script, root+'obj.pro', root+'obj.log')  
              
        cmd = 'sex -c wfc3j_sex.config '+odj_common
        cmd += ' -BACK_SIZE 500 -DETECT_MINAREA '+str(minArea)+' -DETECT_THRESH '+str(thresh)+' -CHECKIMAGE_TYPE "SEGMENTATION" -CHECKIMAGE_NAME '
        cmd += maskName
              
        xcmd(cmd + ' > '+root+'sextractor.log', verbose=False)

        im, header = imOpen(maskName)
        x0 = int(np.round(self.x0))
        y0 = int(np.round(self.y0))

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))

        im, header = imOpen(maskName)
        ax1.imshow(np.flipud(im))
        
        mask2 = seg2mask(maskName, maskName, seg_num=im[x0,y0])      
        im, _ = imOpen(odj_common)

        return ax1, ax2
        
    def backSextract(self, thresh=0.03):
        
        name = self.name
        root = self.objRoot
        config = self.config
        fits_name = self.inFolder+'{}/{}j.fits'.format(name,name)
        odj_common = root+'/tmp'
        segmentation = root + 'back_mask.fits'

        sex_config = config+'sextractor/wfc3j_sex.config'
        PARAMETERS_NAME = config+'sextractor/sbf.param'
        FILTER_NAMEf = config+'sextractor/gauss_2.0_5x5.conv'
        STARNNW_NAME = config+'sextractor/default.nnw'

        script = """
        rd 1 """+fits_name+"""
        rd 2 ./common.mask
        mi 1 2
        wd 1 """+odj_common+"""
        q

        """
        self.run_monsta(script, root+'obj.pro', root+'obj.log')  
            
        cmd = 'sex -c ' + sex_config + ' '+odj_common
        cmd += ' -BACK_SIZE 500 -DETECT_MINAREA 4 -DETECT_THRESH ' + str(thresh) + ' -CHECKIMAGE_TYPE "SEGMENTATION" -CHECKIMAGE_NAME ' + catalName + ' -PARAMETERS_NAME '+ PARAMETERS_NAME+ ' -FILTER_NAME ' + FILTER_NAME + ' -STARNNW_NAME ' + STARNNW_NAME

        cmd += segmentation
              
        xcmd(cmd + ' > '+root+'sextractor.log', verbose=False)
        
        mask2 = seg2mask(segmentation, segmentation)      
        im, _ = imOpen(odj_common)
        
        masked_image = im*mask2
        a = masked_image
        a = a[(a!=0)]
        std = np.std(a)
        median = np.median(a)
        a = a[((a>median-3.*std)&(a<median+3.*std))]

        self.sky_med = np.median(a)
        self.sky_ave = np.mean(a)
        self.sky_std = np.std(a)

        self.masked_image = masked_image
        self.back_pixels = a

        return im, mask2
        
            
    def run_monsta(self, script, Monsta_pro, Monsta_log):

        with open(Monsta_pro, 'w') as f:
            f.write(script)

        cmd = self.monsta+' '+Monsta_pro
        xcmd(cmd + ' > '+Monsta_log, verbose=False)

        with open(Monsta_log) as f:
            text = f.read()
        Tsplit = text.upper().split()
        if "ERROR" in Tsplit or 'SEGMENTATION FAULT' in Tsplit:
            print(text)
            return text
        else:
            return "OK"
            
          
    def getColName(self, catalName):
        with open(catalName, 'r') as f:

            lines = f.readlines()

        col_names = []
        i = 0 
        while lines[i].split()[0]=="#":
            col_names.append(lines[i].split()[2])   
            i+=1
            
        return col_names
    
    def set_center(self, x0, y0):
        self.x0 = x0
        self.y0 = y0
    
    def get_center(self):
        return self.x0, self.y0

    
    def plot_jpg(self, jpg_name, ax=None):
        
        img = mpimg.imread(jpg_name)
        x_max, y_max, _ = img.shape
        
        if ax is None:
            plt.figure(figsize=(10,10))
            plt.subplot(111)
            ax = plt.gca()

        ax.set_xlim([0, self.x_max])
        ax.set_ylim([0, self.y_max])

        imgplot = ax.imshow(np.flipud(img))
        
        return ax        
        
        
    def plot_mask(self, mask=0, ax=None):   
        
        root = self.objRoot
        suffix = '.%03d'%mask
        
        jpg_name = root+'/mask'+suffix+'.jpg'
        
        return self.plot_jpg(jpg_name, ax=ax)        
    
    
    def plot_resid(self, model=0, ax=None):   
        
        root = self.objRoot
        suffix = '.%03d'%model
        
        jpg_name = root+'/resid'+suffix+'.jpg'
        
        return self.plot_jpg(jpg_name, ax=ax)

    def plot_object(self, model=0, ax=None):   
        
        root = self.objRoot
        suffix = '.%03d'%model
        
        jpg_name = root+'/'+self.name+suffix+'.jpg'    
        
        return self.plot_jpg(jpg_name, ax=ax)
    
    
    def list_ellipses(self, model=0):   
        
        root = self.objRoot
        suffix = '.%03d'%model
        
        ellipseFile = root+'/elliprof'+suffix
        df = pd.read_csv(ellipseFile, delimiter=r"\s+", skiprows=7)
        df = df.apply(pd.to_numeric, errors='coerce')
        
        return list_Ell(df)
            
    

    def plot_ellipse(self, model=0, ax=None, **kwargs):   
        
        root = self.objRoot
        suffix = '.%03d'%model
        
        ellipseFile = root+'/elliprof'+suffix
        df = pd.read_csv(ellipseFile, delimiter=r"\s+", skiprows=7)
        df = df.apply(pd.to_numeric, errors='coerce')
       
        if ax is None:
            plt.figure(figsize=(10,10))
            plt.subplot(111)
            ax = plt.gca()
            ax.set_xlim([0, self.x_max])
            ax.set_ylim([0, self.y_max])

        for i in range(len(df)):
            plot_E(df.iloc[i], ax=ax, **kwargs)
        
        return ax
    
    
    def plot_background(self):   
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14,4))

        ## objects are masked, display background pixels
        ax1 = plot_2darray(self.masked_image, ax=ax1)
        ax1.set_title(self.name, fontsize=14)
        ax1.set_xlabel("X [pixel]", fontsize=12)
        ax1.set_ylabel("Y [pixel]", fontsize=12)

        print("Back Median: %.2f"%self.sky_med)
        print("Back Mean: %.2f"%self.sky_ave)
        print("Back Stdev: %.2f"%self.sky_std)

        ## Histogram of the potential background pixel values
        ax2.hist(self.back_pixels, bins=np.linspace(self.sky_med-3*self.sky_std, self.sky_med+3*self.sky_std, 10), density=True, color='g', alpha=0.7)
        ax2.set_xlabel("pixel value", fontsize=12)
        ax2.set_ylabel("frequency", fontsize=12)
        ax2.set_title("Background", fontsize=14)

        self.tv(options="log", ax=ax3)
        
        return ax1, ax2, ax3
        
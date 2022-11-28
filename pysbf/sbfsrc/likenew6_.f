* Program to do a maximum likelihood analysis of point source density(m,r)
***
* Jul 2, 2020 J. Jensen, make bright cutoff 22
* May 29, 2020 J. Jensen, update J, H, K GCLF peak abs mags
* May 18, 2017 J. Jensen, force user to enter GCLF width sigma
* Sep 22, 2016 J. Jensen, fix bright cutoff default limits for Vega mags
* Nov 20, 2015 J. Jensen corrected K-band Vega mags
* May 20, 2013 J. Jensen HST AB mags
* likenew6.f now includes near-IR AB mags and updated galaxy number count
* slope. Uses Source Extractor input format. Made galaxy and GCLF parameters 
* interactively changeable.
*
* This program is now (13/Feb/2009) called "likenew5.f", and is meant
* to be more general than acslikenew.f.  It can read bigger images and
* much larger catalogs.  You can also input an arbitrary kscale value
* (to scale the size of the masking area in terms of Kron aperture), but
* it selects between using the elliptical aperture (cs<0.5) or a circular
* psf one (cs>0.5) based on the last catalog column, nominally CLASS_STAR;
* this can be all 1's or all 0's (e.g.) to select one type for all.
*
* Revision #2.2.1 - a few updates to the GCLF params - May 31, 2007
*
* Revision #2.2 different masking using info from Sextractor
*               new stompout subroutine; getdata modified
*               AJ -- Feb 28, 2007
*   Input parameters in catalog are now:
*   idx, x, y, r, itype, cxx, cyy, cxy, kronfac, AA, id, mag, sig, class_star
*
* Revision #2.1, accepts command line arguments for I band recalculation:
*   expects existing .lkr file and .dmask file; 
*   write .rlk file but no masks, no stdout output
*
*   Syntax: likelier galaxy dist(Mpc),   e.g.   likelier n1399 18
*
* Revision #2, using a numerical integration of the probability function
* so as to allow non-circular spatial domains of integration.
*
* Add a correction to residual variance for soft completeness limit
*
* The output format is also changed to make it more intelligible.
* 
***
      subroutine likenew(fname,n_fname,icolor,secperpixel_,fwhm,
     $                 distance, kscale, delta_, snlim, 
     $                 mlim_, mname, n_mname, yessoft, abmags, 
     $                 beta_, cnorm_, cmax_, alpha_, tot_gc,
     $                 gamma_, gnorm, tnorm, tot_gal, galpersec_,
     $                 tysonempersec_, status) 

C      parameter (maxpt=10000, maxvar=3, maxdim=2048, pi=3.14159265)
      parameter (maxpt=75000, maxvar=3, maxdim=10000, pi=3.14159265)
      parameter (nannuli=5, npred=100)
      parameter (maxhead=1000)
      parameter (maxcolor=9)
      character*1000 fname, mname, outname
      integer n_fname, icolor, status, n_mname
      character colname(0:maxcolor-1)
      real brightcut(0:maxcolor-1), brightcutvega(0:maxcolor-1)
      real r(maxpt), m(maxpt), dm(maxpt), rlim(2), mlim(4), mlim_
      real x(maxpt), y(maxpt)
      real*8 par(20), cnorm, beta, alpha, cmax, gamma, delta
      real*8 likely, kscale
      real*8 xi(maxvar*maxvar), fret
      real*8 bsc, bz, cstar(maxpt)
      real*8 cxx(maxpt),cyy(maxpt),cxy(maxpt),kron(maxpt),aa(maxpt)
      integer*2 bitmap(maxdim*maxdim)
      integer npix(maxdim)
      character*80 header(maxhead)
      real emmarg(npred), rmarg(npred), fracmarg(npred)
      real fracann(nannuli), rann(nannuli), dens(npred,2+nannuli)
      integer count(npred,2+nannuli)
      logical yessoft, relike, abmags, interact
      equivalence (par(1),beta), (par(2),alpha), (par(3),cmax)
      equivalence (par(4),cnorm), (par(5),gamma), (par(6),delta)
      external likely
      common /data/ secperpixel, nobs, x, y, r, m, dm, rlim, mlim, npix,
     $      galpersec
      common /params/ par, nfit
      common /test/ itest, niter

      data interact /.false./
*      data colname /'b','v','r','i','j','h','k'/
* jpb, 03/July/2004, accomodate the z-band
* jj, 20/Nov/2015 add Ks band
      data colname /'b','v','r','i','z','j','h','k','s'/
* YESSOFT was a bad idea...but now it lives in bestfluc
C      data yessoft /.false./

* Use AB mags by default
    !   data abmags /.true./
*      data brightcut /21.0,20.5,20.0,19.5,19.5,20.5,20.25,19.45,
      data brightcut /21.0,20.5,20.0,19.5,19.5,22,20.25,19.45,
     $      19.45/
* Vega mags defaults
      data brightcutvega /21.0,20.5,20.0,19.5,19.5,19.0,19.0,17.5,
     $      17.5/

 6000 format(1x,a,$)
 6001 format(a)
 6007 format(a,'(',f5.2,'):',$)

      itest = 0
      relike = .false.
    !   kscale = 1.2

      write(6,*) 'LIKENEW6 -- IR version 2020'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C First of all, read the data
C         fname = "u12517se.lknj"
C         icolor = 5
C Since we are using only one color from SExtractor via mkincat...
         icol = 1
         secperpixel = secperpixel_
        !  fwhm = 1.4
        !  distance = 60

!          if(fwhm.lt.0) then
!             fwhm = -fwhm
! C 'Kron radius scale factor (def 1.1): '
!             kscale = 1.1
!          end if

C This is the bright end cutoff, or set mlikenew6.so;anual
        !   mlim(4) = 22.
C         mlim(4) = brightcut(icolor)
C         mlim(4) = brightcutvega(icolor)

        !  snlim = 4.5
C         yessoft = .true.

C Enter GC width (1.2-1.4 mag)
         delta = delta_

C Now ask for a bitmap file and calculate areas.
        !  mname = "u12517sej.ptm6"

         if(yessoft) then
            write(6,*) 'Using soft cutoff bias correction.'
         end if


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call getdata(nobs,x,y,x0,y0,r,icol,m,dm,rlim,snlim,mlim,
     $             cxx,cyy,cxy,kron,aa,cstar,fname)

      nhead = maxhead
      call rfits(mname,nhead,header,ibitpix,nx,ny,bsc,bz,bitmap)
      if(ibitpix.ne.1) then
         write(6,*) mname(:lnb(mname)), ' is not a bitmap!', ibitpix
         status = 1
         return
      end if
      call bitshort(nx*ny, bitmap, bitmap)

      
C Find out what the Vista coordinate offsets are
C disallow use of CRVAL's -- too much trouble in this WCS era!
      idx = 0
      idy = 0
      do 88 j = 1,nhead
         i = index(header(j),'CNPIX1')
         if(i.gt.0) then
            read(header(j)(i+10:),'(g20.0)') dx
            idx = nint(dx)
         end if
         i = index(header(j),'CNPIX2')
         if(i.gt.0) then
            read(header(j)(i+10:),'(g20.0)') dy
            idy = nint(dy)
         end if
 88   continue

C Count the number of pixels not exised from the image as a function of r.
      do 9 i = 1,maxdim
         npix(i) = 0
 9    continue

      ignore = 3
      npixtot = 0
      do 10 j = 1+ignore,ny-ignore
         ypix = j - 1 + idy - y0
         do 11 i = 1+ignore,nx-ignore
            xpix = i - 1 + idx - x0
            k = (j-1)*nx + i
            if(bitmap(k).eq.1) then
               ir = max(1,nint(sqrt(xpix*xpix + ypix*ypix)))
               npix(ir) = npix(ir) + 1
               npixtot = npixtot + 1
            end if
 11      continue
 10   continue

C Count in bounds points and limits
      nok = 0
      rmin = 100000
      rmax = 0
      emmin = 100000
      emmax = 0
      emenv = 0
      do 12 i = 1,nobs
         rmin = amin1(rmin,r(i))
         rmax = amax1(rmax,r(i))
         emmin = amin1(emmin,m(i))
         emmax = amax1(emmax,m(i))
         emcut = 2.5*(mlim(1) - mlim(2)*alog10(1+mlim(3)/r(i)))
         emenv = amax1(emenv,emcut)
         if(r(i).lt.rlim(1).or.r(i).gt.rlim(2)) goto 12
         if(m(i).lt.mlim(4).or.m(i).gt.emcut) goto 12
         nok = nok + 1
 12   continue

      if(.not.relike) write(6,7362) nobs, nok, nint(mlim(3))
 7362 format(i5,' sources found,',i5,' used for fits;  sky radius =',i4)

C Get initial values for the parameters
C Use distance < 0 as a signal to fit for it      
      if(distance.lt.0) then
         nfit = 3
         distance = -distance
      else
         nfit = 2
      end if

C Initial guess for gxy/total number
      beta = 0.5

C Initial guess for d lnN / d lnr of GC distribution
      alpha = 0.6

C Assume icolor = 0/1/2/3/4/5/6/7/8 for B/V/R/I/z/J/H/K/Ks
C GC maximum brightness: (from Harris Extragalactic Distance Scale)
C Calibration on MW gives  B0 = -7.0 +/- 0.15
C Calibration on M31 gives  B0 = -6.85  with (m-M)V = 24.6 assumed
C Calibration on Virgo gives B0 = -7.0  with (m-M) = 31.7 assumed
C Observed Virgo B0 = 24.7, adopting (m-M) = 30.9 (15 Mpc) gives B0 = -6.2
C Since these results will be used in conjunction with fluctuation analysis
C   which finds a Virgo distance of 15 Mpc, let us abuse Prof. Harris a
C   bit and adjust his absolute magnitude from -7.0 to -6.2.
C M_V = -7.0, B-V = 0.80, V-R = 0.45, R-I = 0.50, 
C and averages from Table 2, Frogel et al 1978 ApJ 220 75,
C where V-J/V-H/V-K = 1.62,2.14,2.23
      distmod = 25 + 5*alog10(distance)
C      vabs = -7.0
C jpb, change 21/Apr/2004
C Note: z-band assumes AB calibration, while the rest are Vega.
C jpb, 31/May/2007:
C     update z-band number following Andres's results, Mz0=-8.4
C JJ, 20/Nov/2015: add Ks band
      vabs = -7.4
      if(interact) then
*         write(6,6007) ' Enter Mv for globular clusters (Vega)',vabs
         write(6,6007) ' Enter Mv for GCs (Vega), NOT used for JHK',vabs
         read(5,*) vabs
      else
         write(6,6010) 'Using Mv = ',vabs,' for GCs (NOT for JHK)'
 6010    format(1x,a,f5.2,a)
      end if

C JJ, 29 May 2020, using Nantais et al. 2006 AJ 131, 1416  
C     M_J=-9.18, M_H=-9.73, M_K=-9.98 Vega
      if(abmags) then
         write(6,*) 'Using AB mags for JHK'
         if(icolor.eq.0) emabs = vabs + 0.8
         if(icolor.eq.1) emabs = vabs - 0.0
         if(icolor.eq.2) emabs = vabs - 0.5
         if(icolor.eq.3) emabs = vabs - 1.0
         if(icolor.eq.4) emabs = vabs - 1.0
*         if(icolor.eq.5) emabs = vabs - 1.62 + 0.7595
*         if(icolor.eq.6) emabs = vabs - 2.14 + 1.2514
*         if(icolor.eq.7) emabs = vabs - 2.23 + 1.95
*         if(icolor.eq.8) emabs = vabs - 2.23 + 1.85
         if(icolor.eq.5) emabs = -9.18 + 0.7595
         if(icolor.eq.6) emabs = -9.73 + 1.2514
         if(icolor.eq.7) emabs = -9.98 + 1.95
         if(icolor.eq.8) emabs = -9.98 + 1.85

      else
         write (6,*) 'Using Vega mags'
         if(icolor.eq.0) emabs = vabs + 0.8
         if(icolor.eq.1) emabs = vabs - 0.0
         if(icolor.eq.2) emabs = vabs - 0.5
         if(icolor.eq.3) emabs = vabs - 1.0
         if(icolor.eq.4) emabs = vabs - 1.0
*         if(icolor.eq.5) emabs = vabs - 1.62
*         if(icolor.eq.6) emabs = vabs - 2.14
*         if(icolor.eq.7) emabs = vabs - 2.23
*         if(icolor.eq.8) emabs = vabs - 2.23
         if(icolor.eq.5) emabs = -9.18
         if(icolor.eq.6) emabs = -9.73
         if(icolor.eq.7) emabs = -9.98
         if(icolor.eq.8) emabs = -9.98
      end if
      cmax = distmod + emabs

C d logN / d m for gxy, Assume icolor = 0/1/2/3 for B/V/R/I, 
C Tyson gives slope = 0.48/0.42/0.39/0.34
C Gardner ApJL 415 L9 and Cowie ApJ 434 114 find more or less that
C N = 10^4 * 10^{0.3(K-19)} gxy/deg^2/mag and gamma = 0.3 for JHK
C
C May 2013: udpated galaxy slope:
C Retzlaff et al. 2010, A&A, 511, A50: GOODs survey results for JHK
C Windhorst et al. 2011, ApJ, 193, 27: HST WFC3 HUDF measurements
C N = 0.03 * 10^(0.25*AB) gxy/deg^2/mag, gamma = 0.25 for JHK
C This undercounts at brighter mags AB<20 but is a much better fit for
C AB = 21 - 25 and fainter. Note that the mag at which you get 1 gal/asec^2
C is now more like AB = 34.5
      if(icolor.eq.0) gamma = 0.48
      if(icolor.eq.1) gamma = 0.42
      if(icolor.eq.2) gamma = 0.39
      if(icolor.eq.3) gamma = 0.34
      if(icolor.eq.4) gamma = 0.35

      if(icolor.eq.5) gamma = 0.25
      if(icolor.eq.6) gamma = 0.25
      if(icolor.eq.7) gamma = 0.25
      if(icolor.eq.8) gamma = 0.25

      if(interact) then
         write(6,6007) ' Enter galaxy slope gamma ',gamma
         read(5,*) gamma
      end if         
C Galaxy density normalization from Tyson
      if(abmags) then
         galpersec = 34.5
         if(icolor.eq.0) tysonempersec = 30.25
         if(icolor.eq.1) tysonempersec = 30.45
         if(icolor.eq.2) tysonempersec = 30.55
         if(icolor.eq.3) tysonempersec = 30.60
         if(icolor.eq.4) tysonempersec = 30.60

         if(icolor.eq.5) tysonempersec = 29.37 + 0.7595
         if(icolor.eq.6) tysonempersec = 29.37 + 1.2514
         if(icolor.eq.7) tysonempersec = 29.37 + 1.95
         if(icolor.eq.8) tysonempersec = 29.37 + 1.85
      else
         galpersec = 30.5
         if(icolor.eq.0) tysonempersec = 30.25
         if(icolor.eq.1) tysonempersec = 30.45
         if(icolor.eq.2) tysonempersec = 30.55
         if(icolor.eq.3) tysonempersec = 30.60
         if(icolor.eq.4) tysonempersec = 30.60

         if(icolor.eq.5) tysonempersec = 29.37
         if(icolor.eq.6) tysonempersec = 29.37
         if(icolor.eq.7) tysonempersec = 29.37
         if(icolor.eq.8) tysonempersec = 29.37
      end if
      if(interact) then
         write(6,6007) ' Enter mag for 1 galaxy per arcsec^2 ',
     $        galpersec
         read(5,*) galpersec
         write(6,6007) ' Enter Tyson mag for 1 gal per arcsec ',
     $        tysonempersec
         read(5,*) tysonempersec
      end if

C Width of GC distribution
C jpb, 31/May/2007:
C     update this to reflect ACS results in Virgo, fornax, etc:
C      delta = 1.2
C      delta = 1.35
C      delta = 1.4
C      if(interact) then
C      write(6,6007) ' Enter GC width (1.2-1.4 mag)',delta
C      read(5,*) delta
C      else
C         write (6,6010) 'Using width delta = ',delta,' for GCLF'
C      end if

      write(6,6010) 'Using galaxy slope gamma = ',gamma,' '
      write(6,6010) 'Using ',galpersec,' normalization for 1 gal/arcsec'

C      write(6,6000) 'Enter beta, alpha, cnorm, cmax, gamma, delta: '
C      read(5,*) (par(i),i=1,6)

C Put the magnitude limit parameters into PAR
      do 21 i = 1,4
         par(i+6) = dble(mlim(i))
 21   continue

      if(itest.gt.1) then
         do 6 i = 1,nobs
            write(2,2000) i, x(i), y(i), r(i), m(i), dm(i)
 2000       format(i5,3f8.1,5f8.2)
 6       continue
         do 7 i = 1,10
            dr = (rlim(2)-rlim(1))/9
            rr = rlim(1) + dr*float(i-1)
            emcut = 2.5*(mlim(1) - mlim(2)*alog10(1+mlim(3)/rr))
            n = 0
            n23 = 0
            do 8 j = 1,nobs
               if(r(j).ge.rr.and.r(j).lt.rr+dr) then
                  n = n + 1
                  if(m(j).lt.23) n23 = n23 + 1
               end if
 8          continue
            write(3,2000) i, rr, emcut, float(n), n/(rr+0.5*dr), 
     $           float(n23), n23/(rr+0.5*dr)
 7       continue
      end if

      IF(ITEST.GT.1) THEN
 1       WRITE(6,*)
         WRITE(6,6000) 'Enter beta and alpha: '
         READ(5,*,END=2) (PAR(I),I=1,MAXVAR)
         BOZO = LIKELY(PAR)
C     WRITE(6,4000) (PAR(I),I=1,6), BOZO
 4000    FORMAT(1P2G12.3,0P4F8.3,1PG12.3)
         GOTO 1
      END IF
 2    CONTINUE

C Derive a fit to the point source distribution
      nderiv = 0
      ierr = 0
      if(.not.relike) ierr = 1
      call mini(nfit,par,likely,dum,dum,xi,nderiv,ierr)
      fret = likely(par)
      if(ierr.eq.2) then
         write(6,*) 'SINGULAR MATRIX!'
         status = 1  ! stop here
         return
      end if

      if(.not.relike) then
         write(6,*)
         write(6,4001) niter, fret
 4001    format(i5,' iterations.  Final (-) likelihood / pt =',f8.3)
         write(6,*)
         write(6,*) 'Covariance matrix:'
         do 22 j = 1,nfit
            write(6,4002) (xi(i+nfit*(j-1)),i=1,nfit)
 4002       format(1p3g12.3)
 22      continue
      end if

C Tell us about the results
      csum = par(12)
      gsum = par(13)
      gnorm = beta * (nok) / gsum
C      tnorm = gnorm * 10**(gamma*(tysonempersec-30.5))
      tnorm = gnorm * 10**(gamma*(tysonempersec-galpersec))
      gcdist = distance * 10.**(0.2*(cmax-(distmod+emabs)))

      if(.not.relike) then
         write(6,*)
         write(6,4003) ' Beta =      ', beta, '#gxy / total sources'
         write(6,*)
         if(nfit.eq.3) then
            write(6,4007) 'Distance =   ', distance, 'assumed, Mpc'
            write(6,4007) 'Distance =   ', gcdist, 'derived, Mpc'
         end if
         write(6,4003) 'Cnorm =      ', cnorm, 'GC normalization'
         write(6,4007) ' Cmax =      ', cmax, 'GC peak magnitude'
         write(6,4007) 'Delta =      ', delta, 'GC distribution width'
         write(6,4003) 'Alpha =      ', alpha, 'GC log slope vs log r'
         write(6,4004) 'Total # GC = ', csum*cnorm
         write(6,*)
         write(6,4007) 'Gamma =      ', gamma, 'Gxy log slope vs m'
C         write(6,4007) 'Gnorm =      ', gnorm, 'Gxy count / 1/" @ 30.5'
         write(6,4008) 'Gnorm =      ', gnorm, 'Gxy count / 1/" @ ',
     $      galpersec
         write(6,4008) 'Tnorm =      ', tnorm, 'Tyson Gxy count @ ',
     $      tysonempersec
         write(6,4004) 'Total # gxy =', gsum*gnorm

         tot_gc = csum*cnorm
         tot_gal  = gsum*gnorm

         cnorm_ = cnorm
         cmax_ = cmax
         alpha_ = alpha
         gamma_ = gamma
         beta_ = beta

         galpersec_ = galpersec
         tysonempersec_ = tysonempersec

 4003    format(a,1pe12.4,5x,'(',a,')')
 4007    format(a,f12.5,5x,'(',a,')')
 4004    format(a,f8.1)
 4008    format(a,f12.5,5x,'(',a,f5.2,')')
      end if

C Compute residual variance as a function of radius and luminosity function
C We will write a file called xxx.lkr which has parameters
C observed and predicted N(m), N(r), res_var(r), etc

C NANNULI is the number of annuli in which results are computed
C NMARG  is the number of bins for the marginal distributions
C NPRED  is the number of points where predicted results are reported

      ind = index(fname,'.') - 1
      if(ind.lt.0) ind = lnb(fname)
      if(relike) then
         write(outname,'(a,a,a)') fname(:ind), colname(icolor), '.rlk'
      else
         write(outname,'(a,a,a)') fname(:ind), colname(icolor), '.lkn6'
      end if
      open(unit=2,file=outname,status='unknown')

C Initialize the counts and ranges for the marginal sums of the distribution
      embin = 0.5
      em0 = embin * int(emmin/embin)
      em1 = embin * int(emmax/embin+0.999)
      nmarg = min(nint((em1-em0)/embin),npred)
      do 31 j = 1,nmarg
         emmarg(j) = em0 + (j-0.5)*embin
         do 311 i = 1,2+nannuli
            count(j,i) = 0
 311     continue
 31   continue

C Compute the mean radii and average filling of each annulus...
C The annuli are defined between rlim(1) and rlim(2), the ranges of the fits.
      do 312 j = 1,max(nmarg,nannuli)
         if(j.le.nannuli) then
            rann(j) = 0
            fracann(j) = 0
         end if
         if(j.le.nmarg) then
            rmarg(j) = 0
            fracmarg(j) = 0
         end if
 312  continue
      do 313 j = nint(rmin),nint(rmax)
         if(j.ge.nint(rlim(1)).and.j.le.nint(rlim(2))) then
            ir = int(nannuli*(j-rlim(1))/(rlim(2)-rlim(1))) + 1
            if(ir.ge.1.and.ir.le.nannuli) then
               fracann(ir) = fracann(ir) + npix(j)
               rann(ir) = rann(ir) + j*npix(j)
            end if
         end if
         ir = int(nmarg*(j-rmin)/(rmax-rmin)) + 1
         if(ir.ge.1.and.ir.le.nmarg) then
            fracmarg(ir) = fracmarg(ir) + npix(j)
            rmarg(ir) = rmarg(ir) + j*npix(j)
         end if
 313  continue
      do 314 j = 1,max(nmarg,nannuli)
         if(j.le.nannuli) then
            rann(j) = rann(j) / amax1(1.,fracann(j))
            r0 = rlim(1) + float(j-1)/nannuli*(rlim(2)-rlim(1))
            r1 = rlim(1) + float(j)/nannuli*(rlim(2)-rlim(1))
            fracann(j) = fracann(j) / (pi*(r1*r1-r0*r0))
         end if
         if(j.le.nmarg) then
            rmarg(j) = rmarg(j) / amax1(1.,fracmarg(j))
            r0 = rmin + float(j-1)/nmarg*(rmax-rmin)
            r1 = rmin + float(j)/nmarg*(rmax-rmin)
            fracmarg(j) = fracmarg(j) / (pi*(r1*r1-r0*r0))
         end if
 314  continue

C Add up the numbers falling in each bin.  Get both marginal distributions
C in both m and r, and also marginal distribution as a function of m in
C each of the annuli
C The counts are saved as COUNT(J,K), where J is the index of the marginal
C variable, and K refers to 1/2/3-N: r_distribution/m_distribution/annuli
      do 32 i = 1,nobs
         j = int(nmarg*(r(i)-rmin)/(rmax-rmin)-0.001) + 1
         if(j.ge.1.and.j.le.npred) then
            count(j,1) = count(j,1) + 1
         end if
         k = int((m(i)-em0)/embin) + 1
         if(k.ge.1.and.k.le.npred) then
            count(k,2) = count(k,2) + 1
            j = int(nannuli*(r(i)-rlim(1))/(rlim(2)-rlim(1))-0.001) + 1
            if(j.ge.1.and.j.le.nannuli) then
               count(k,j+2) = count(k,j+2) + 1
            end if
         end if
 32   continue

C Write the header information
 2003 format(i6,5x,a)
 2004 format(13x,i12,5x,a)
      write(2,2003) 26+nannuli, 'lines of header information'
      write(2,2003) nmarg, 'lines of marginal counts'
      write(2,2003) npred, 'lines of predicted fits'
      write(2,2004) nobs, 'points originally present'
      write(2,2004) nok, 'points within (r,m) fit limits'
      write(2,4007) 'Distance =   ', distance, 'distance assumed, Mpc'
      write(2,4007) 'Distance =   ', gcdist, 'distance derived, Mpc'
      write(2,4007) 'sec/pixel =  ', secperpixel, 'pixel scale'
      write(2,4007) 'Likelihood = ', fret, '-ln prob/pt of fit'
      write(2,4007) ' Beta =      ', beta, '#gxy / total sources'
      write(2,4003) 'Cnorm =      ', cnorm, 'GC normalization'
      write(2,4007) ' Cmax =      ', cmax, 'GC peak magnitude'
      write(2,4007) 'Delta =      ', delta, 'GC distribution width'
      write(2,4003) 'Alpha =      ', alpha, 'GC slope log N vs log r'
      write(2,4004) 'Total # GC = ', csum*cnorm
      write(2,4007) 'Gamma =      ', gamma, 'Gxy slope log N vs m'
      write(2,4008) 'Gnorm =      ', gnorm, 'Gxy count / 1/sec^2 @ ',
     $     galpersec
      write(2,4008) 'Tnorm =      ', tnorm, 'Tyson Gxy count @ ',
     $     tysonempersec
      write(2,4004) 'Total # gxy =', gsum*gnorm
      write(2,4007) 'rlim(1) =    ', rlim(1), 'Innermost fit radius'
      write(2,4007) 'rlim(2) =    ', rlim(2), 'Outermost fit radius'
      write(2,4007) 'SN limit =   ', snlim, 'Completeness limit...'
      write(2,4007) 'mlim(1) =    ', mlim(1), 'A: K + log(1/SN)'
      write(2,4007) 'mlim(2) =    ', mlim(2), 'B: ...'
      write(2,4005) 'mlim(3) =    ', mlim(3), 
     $     'C: mcut=2.5*(A-Blog(1+C/r)'
      write(2,4007) 'mlim(4) =    ', mlim(4), 'bright end cutoff'
 4005 format(a,1pg12.3,5x,'(',a,')')
      do 33 i = 1,nannuli
         write(2,4003) 'Annulus r =  ', rann(i), 'mean r of annulus'
 33   continue
      write(2,*)

C Write the marginal distributions
      write(2,*) 'Observed densities in number/magnitude/arcmin^2'
      write(2,*) 'Marginal distributions in magnitude and radius'
      arcminperpix = (secperpixel/60.)**2

      rscale = (rlim(2)-rlim(1))/nannuli
      write(2,2006) nint(rmin),nint(rmax),
     $     (nint(float(i-1)*rscale+rlim(1)),
     $      nint(float(i)*rscale+rlim(1)),i=1,nannuli)
 2006 format(30x,10(i4,' <r<',i4))
      write(2,2005)
 2005 format('   r   frac     N  n(/m/'')   m     ',
     $     '  N n(/m/'')     N n(/m/'')     N n(/m/'')   ',
     $     '  N n(/m/'')     N n(/m/'')     N n(/m/'')')
      do 34 j = 1,nmarg
         r0 = rmin + float(j-1)/nmarg*(rmax-rmin)
         r1 = rmin + float(j)/nmarg*(rmax-rmin)
         emc0 = mlim(4)
         emc1 = 2.5*(mlim(1) - mlim(2)*alog10(1+mlim(3)/rmarg(j)))
         area = fracmarg(j) * pi*(r1*r1-r0*r0) * arcminperpix
         dens(j,1) = count(j,1) / (emc1-emc0) / area
         area = npixtot * arcminperpix
         dens(j,2) = count(j,2) / embin / area
         do 35 i = 1,nannuli
            r0 = rlim(1) + float(i-1)/nannuli*(rlim(2)-rlim(1))
            r1 = rlim(1) + float(i)/nannuli*(rlim(2)-rlim(1))
            area = fracann(i) * pi*(r1*r1-r0*r0) * arcminperpix
            dens(j,i+2) = count(j,i+2) / embin / area
 35      continue
         write(2,2001) rmarg(j), fracmarg(j), count(j,1), dens(j,1), 
     $        emmarg(j), (count(j,i+2), dens(j,i+2),i=0,nannuli)
 2001    format(f6.1,f6.3,i6,f8.2,f6.2,10(i6,f8.2))
 34   continue

C Compute magnitude marginal distributions for gxy, gc, gxy+gc
      write(2,*)
      write(2,*) 'Model distributions in number/magnitude/arcmin^2'
      write(2,2007) nint(rmin),nint(rmax),(nint(float(i-1)*rscale+rmin),
     $     nint(float(i)*rscale+rmin),i=1,nannuli)
 2007 format(10x,'Dens: ',6(i7,' <r<',i4,1x))
      write(2,2008)
 2008 format('   m        gxy  ',6('   GC     both  '))

C Here is the normalization for the galaxy density in #/pixel^2
C      emgalsec = 30.5
C Now updated to Windhorst et al. 2011 AB mag slope: 34.5
C Note that emgalsec is set in two places in this code.
      emgalsec = galpersec
      gscale = gnorm * secperpixel**2 * 10**(-gamma*emgalsec)
      emstep = 0.1
      emextra = 2
      do 40 i = 1,npred
         em = mlim(4) - emextra + (i-1)*emstep

C Get gxy density, save as dens(1,1)
C Get GC total and annuli averages, save as dens(1,*), * = 2, 3...
C Get both total and annuli averages, save as dens(2,*), * = 2, 3...
         dens(1,1) = gscale * 10**(gamma*em) / arcminperpix
         do 411 j = 0,nannuli
            dens(1,j+2) = 0
 411     continue
         gaussian = exp(-0.5*((em-cmax)/delta)**2)
         do 41 j = nint(rmin),nint(rmax)
C Get GC density
            gcdensity = cnorm * float(j)**(-alpha) * gaussian
            if(j.ge.nint(rlim(1)).and.j.le.nint(rlim(2))) then
               ir = int(nannuli*(j-rlim(1))/(rlim(2)-rlim(1))) + 1
               if(ir.ge.1.and.ir.le.nannuli) then
                  dens(1,ir+2) = dens(1,ir+2) + npix(j)*gcdensity
               end if
            end if
            dens(1,2) = dens(1,2) + npix(j)*gcdensity
 41      continue
         dens(1,2) = dens(1,2) / npixtot / arcminperpix
         dens(2,2) = dens(1,2) + dens(1,1)
         do 412 j = 1,nannuli
            r0 = rlim(1) + float(j-1)/nannuli*(rlim(2)-rlim(1))
            r1 = rlim(1) + float(j)/nannuli*(rlim(2)-rlim(1))
            dens(1,j+2) = dens(1,j+2) / 
     $           (fracann(j)*pi*(r1*r1-r0*r0)) / arcminperpix
            dens(2,j+2) = dens(1,j+2) + dens(1,1)
 412     continue

         write(2,2010) em, dens(1,1), ((dens(j,k+2),j=1,2),k=0,nannuli)
 2010    format(f7.2,20f8.2)
 40   continue

C Compute the residual variances
      write(2,*)
      write(2,*) 'Residual variance PER PIXEL in magnitudes, V = 1*f^2'
      write(2,2020) 
 2020 format(7x,'bright  faint   Residual variance')
      write(2,2021)
 2021 format('   r      m      m    m_gc   m_gxy  m_all  slope   bias')
      rstep = nint((1.2*rmax-rmin/2) / (npred-1))
      do 45 i = 1,npred
         rr = (i-1)*rstep + nint(rmin/2)

C Bright and faint magnitude limits
         emb = mlim(4)
         emf = 2.5*(mlim(1) - mlim(2)*alog10(1+mlim(3)/rr))

C residual mbar from GC's
         offset = 0.4 * alog(10.) * delta**2
         arg = (emf-cmax+2*offset)/(sqrt(2.)*delta)
         arglog = cnorm * rr**(-alpha) * sqrt(pi/2) * delta * erfc(arg)
         embargc = cmax - offset - 1.25 * alog10(arglog)

C residual mbar from gxy's
         arglog = gscale * 10**(gamma*emf) / (alog(10.)*(0.8-gamma))
         embargxy = emf - 1.25*alog10(arglog)

C variance enhancemant from flux error bias
         if(yessoft) then
            gxylim = gscale * 10.**(gamma*emf)
            arggc = (emf-cmax)/delta
            gclim = cnorm* rr**(-alpha) * exp(-0.5*arggc*arggc)
C This is from the erroneous derivation of 8/19/91
C            fermi = 4/sqrt(2*pi/snlim)
            fermi = 4*snlim/sqrt(2*pi)
            slopelim = 1 + 
     $           (gamma*alog(10.)*gxylim-arggc*gclim)/(gxylim+gclim)
            rat = 1+5./3.*(3-slopelim)**2/(fermi*fermi-(3-slopelim)**2)
            
C residual mbar from both
            embar = -1.25 *
     $           alog10(rat*(10**(-0.8*embargc)+10**(-0.8*embargxy)))
         else
C residual mbar from both
            embar = -1.25*alog10(10**(-0.8*embargc)+10**(-0.8*embargxy))
         end if

C Now write the residual variances
         write(2,2025) rr, emb, emf, embargc, embargxy, embar, 
     $        slopelim, rat
 2025    format(f6.1,7f7.2)
 45   continue

      close(2)

* If we are just recomputing things, we are done now...
      if(relike) goto 999

C Compute a bitmap which is masked at the locations of the point sources
C Use a psf p(r) = A { exp(-1/2*r^2/s^2) + 0.4 * exp(-r/s) }
C For this psf the FWHM = 2.14*s and Total flux = 1.4 * 2 * pi * s^2 * A
      s = fwhm / 2.14

C Excise down to this surface brightness (mag/pix^2) limit
C jpb: make fainter for ACS
C      surfcut = 30
      surfcut = 32
      rmin = amax1(2.,fwhm/2)

      write(6,2026) 'Making bitmap with kscale = ', kscale
 2026 format(1x,a,f5.2)
C Punch out the bitmap around each point source.
      do 50 i = 1,nobs
         emcut = 2.5*(mlim(1) - mlim(2)*alog10(1+mlim(3)/r(i)))
         if(m(i).gt.emcut) goto 50
         ix = x(i) - idx
         iy = y(i) - idy
C The surface brightness A of this object is
         surf = m(i) + 2.5*alog10(1.4*2*pi*s*s)
C Ignore the gaussian center, and just use the exponential skirt
C         emcut = surf - 2.5*alog10(0.4*exp(-r/s))
         rc = s * alog(0.4*10**(0.4*(surfcut-surf)))
         if(rc.lt.fwhm) then
            frac = exp(-(1-rc/fwhm))
            rcut = frac*fwhm + (1-frac)*rmin
         else
            rcut = rc
         end if

         if (cstar(i).lt.0.6) then
            call stompout(nx,ny,bitmap,ix,iy,cxx(i),cyy(i),
     $           cxy(i),kscale*kron(i),kscale*kron(i)*aa(i))
         else
            call stompout_old(nx,ny,bitmap,ix,iy,rcut)
         end if

 50   continue

C Convert the array to a bitmap and write it to Nxxxx_?.ptm
      write(outname,'(a,a,a)') fname(:ind), colname(icolor), '.ptm6'
      call shortbit(nx*ny, bitmap, bitmap)
      ibitpix = 1
      call wfits(outname,nhead,header,ibitpix,nx,ny,bitmap)

C      call i2tob1(bitmap,nx*ny,bitmap)
C      call wfbitfile(nhead,header,nx,ny,bitmap,outname)

 999  status = 10
      return 
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function likely(p)
C
C Data is stored in /data/ common:
C N = npoints
C R = radius
C M = magnitude
C DM = magnitude error
C RLIM(2) = rmin, rmax
C MLIM(4) = A,B,C, mmin, where mmax = 2.5*(A-Blog(1+C/R))
C Data is guaranteed to all lie within these limits
C
C Parameters are stored in /params/ common and referred to as equivalenced
C #gxy/mag/pixel^2(m) = GNORM * 10**GAMMA * (m - mg)
C #GC/mag/pixel^2(m,r) = CNORM * r**-ALPHA * exp(-((m)-CMAX)^2/(2*DELTA^2))
C NFIT refers to how many of these parameters are being varied by the
C minimization routine. All are initially set, however.
C
C The integrated number density is forced to the number observed
C
C Beyond the basic 6, PAR also holds magnitude limit cofficients A,B,C, (7-9)
C and the lower magnitude limit (10)
      parameter (maxpt=75000, maxdim=10000)
      implicit real*8 (a-h,o-z)
      real*8 p(6)
      real*8 par(20), beta, cnorm, alpha, cmax, gamma, delta
      real*8 oldalpha, oldcmax, oldgamma, olddelta
      real*8 r1, r2, likely, betaclip
      real*4 r(maxpt), m(maxpt), dm(maxpt), rlim(2), mlim(4)
      real*4 x(maxpt), y(maxpt), secperpixel, galpersec
      integer npix(maxdim)

      real pbuf(100,100)
      character*80 pcom(10)

      equivalence (par(1),beta), (par(2),alpha), (par(3),cmax)
      equivalence (par(4),cnorm), (par(5),gamma), (par(6),delta)
      common /data/ secperpixel, nobs, x, y, r, m, dm, rlim, mlim, npix,
     $              galpersec
      common /params/ par, nfit
      COMMON /TEST/ ITEST, ncall
      data oldalpha /0d0/, oldcmax /0d0/, oldgamma /0d0/, olddelta /0d0/
      data ncall /0/

      ncall = ncall + 1

C Update the saved parameters to the new parameters
      do 10 i = 1,nfit
         par(i) = p(i)
 10   continue

      alpha = dabs(alpha)

      r1 = dble(rlim(1))
      r2 = dble(rlim(2))

      acc = 1d-10
      pi = 4d0*datan(1d0)
      twopi = 2d0*pi

C log N (sec^-2 mag^-2) = GAMMA * (m - mg)
C Tyson's prediction of how many background galaxies we should count
C is GSUM = Int { 10**(GAMMA*m) } * SECPERPIXEL^2 * 10**(-GAMMA*MG)
C We will use JT's rule: in any color there is 1 gxy/mag/sec^2 at mg=30.5
C      emgalsec = 30.5
C JJ: I'm afraid this doesn't work with our new slope and normalization
C in AB mags. Using Windhorst 2011 slope 0.25 and
C N = 0.03 * 10^(0.25*m) per AB mag per square degree
C we get one gal per square arcsec at 34.5 AB. 
C Note that emgalsec gets set in two places in this code.
      emgalsec = galpersec
      tysonscale = secperpixel**2 * 10**(-gamma*emgalsec)

C If GAMMA has changed we need to recalculate the ML normalization
      if(gamma.ne.oldgamma) then
         par(11) = 10d0**(gamma*par(10))
         gsum = 0
         do 11 i = idnint(r1),idnint(r2)
            gsum = gsum + npix(i)*grand(par,dfloat(i))
 11      continue
         gsum = gsum * tysonscale / (gamma*dlog(10d0))
         oldgamma = gamma
      end if

C If ALPHA, DELTA, or CMAX has changed, we need to redo more ML norm.
      if(alpha.ne.oldalpha.or.delta.ne.olddelta.or.cmax.ne.oldcmax) then
         par(11) = derfc((cmax-par(10))/(dsqrt(2d0)*delta))
         csum = 0
         do 12 i = idnint(r1),idnint(r2)
            csum = csum + npix(i)*crand(par,dfloat(i))
 12      continue
         csum = dsqrt(pi/2d0)*delta * csum

         oldalpha = alpha
         olddelta = delta
         oldcmax = cmax
      end if

C Make sure that mini isn't giving us a bogus BETA
      if(beta.gt.0d0.and.beta.lt.1d0) then
         betaclip = beta
      else
         betaclip = dmin1(1d0,dmax1(0d0,beta))
         write(6,'(a,f12.4)') 'LIKELY: bad value for BETA',beta
         write(6,*) betaclip
      end if
C Now we must sum over all observations their unnormalized probability
      sum = 0
      nok = 0

      do 20 i = 1,nobs

C Skip out of bounds points
         if(r(i).lt.rlim(1).or.r(i).gt.rlim(2)) goto 20
         if(m(i).lt.mlim(4)) goto 20
         emcut = 2.5*(mlim(1) - mlim(2)*alog10(1+mlim(3)/r(i)))
         if(m(i).gt.emcut) goto 20

         nok = nok + 1

         gdens = tysonscale*twopi*dble(r(i)) * 10d0**(gamma*dble(m(i)))
         carg = -0.5d0*((dble(m(i))-cmax)/delta)**2
         cdens = twopi * dble(r(i))**(1-alpha) * dexp(carg)
         prob = dlog(betaclip*gdens/gsum+(1d0-betaclip)*cdens/csum)

         sum = sum + prob
C         write(6,3672) i, r(i), m(i), gdens, carg, cdens, prob, sum
 3672    format(i4,f6.1,f6.2,1p5g12.3)
 20   continue

C Now we can produce (-) the maximum likelihood estimator

      likely = -sum / nok

C We constrain NOK = CNORM*CSUM + GNORM*GSUM, and set 
C BETA = GNORM*GSUM/(CNORM*CSUM + GNORM*GSUM) = GNORM*GSUM / NOK
C Let us keep gnorm from dropping below 0.5
      betamin = 0.5 * gsum / nok
      betamax = 1

C While we're at it we can set CNORM
      cnorm = nok*(1-betaclip)/csum

C Adjust LIKELY up to keep it (more or less) within our BETA limits
       if(beta.gt.1d0 .or. beta.lt.0d0) then
         wall = (beta-0.5)**6
      else
         dbeta = 0.01
         wallmax = dexp((beta-betamax)/dbeta)
         wallmin = dexp((betamin-beta)/dbeta)
         wall = (wallmax + wallmin)/nok
      end if
      if(wall.gt.10000) then
         write(6,'(a,1pe12.2)') 'BETA wall hit very hard   ', wall
         wall = dlog(wall)**2
      end if
      likely = likely + wall

      IF(ITEST.GT.0) THEN
         WRITE(6,4000) NCALL,(PAR(I),I=1,2),
     $        NOK,CNORM,GSUM,CSUM,SUM,LIKELY
 4000    FORMAT(I5,F7.3,F7.2,I5,1P5G11.3)
      END IF

C Communicate some interesting parameters back to main
      par(12) = csum
      par(13) = gsum

C Mongo output of the probability contours?
      if(itest.gt.1) then
         amag0 = 19
         amag1 = 25
         r0 = 10
         r1 = 500
         do 3 k = 1,100
            amag = amag0+float(k-1)/99*(amag1-amag0)
            do 4 j = 1,100
               emf = par(7)+rr*(par(8)+rr*par(9))
               rr = r0+float(j-1)/99*(r1-r0)
               gdens = twopi*rr * 10d0**(gamma*dble(amag))
               carg = -0.5d0*((dble(amag)-cmax)/delta)**2
               cdens = twopi * rr**(1-alpha) * dexp(carg)
               pbuf(j,k) = dlog(beta*gdens/gsum+(1d0-beta)*cdens/csum)
C     if(amag.gt.emf) pbuf(j,k) = -30
 4          continue
 3       continue
         pcom(1) = 'limits 10 500 19 25'
C         call mongo(1,pcom,100,100,pbuf)
      end if

      return
      end

      function grand(par,r)
C Integrand for the galaxy part of the ML normalization
      implicit real*8 (a-h,o-z)
      real*8 par(11)
* Here is an emcut calculation from mlim...
      emf = 2.5d0*(par(7)-par(8)*dlog10(1d0+par(9)/r))
      gamma = par(5)
      grand = 10d0**(gamma*emf) - par(11)
      return
      end

      function crand(par,r)
C Integrand for the cluster part of the ML normalization
      implicit real*8 (a-h,o-z)
      real*8 par(11)
      alpha = par(2)
      cmax = par(3)
      delta = par(6)
* Here is an emcut calculation from mlim...
      emf = 2.5d0*(par(7)-par(8)*dlog10(1d0+par(9)/r))
      emgral = derfc((cmax-emf)/(dsqrt(2d0)*delta)) - par(11)
      crand = r**(-alpha) * emgral
C      write(6,4000) 'C:',r,emf,cmax,emgral,crand
 4000 format(a,5f12.3)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getdata(n,x,y,x0,y0,r,ic,m,dm,rlim,snlim,mlim,
     $                   cxx,cyy,cxy,kron,aa,cstar,fname)
      parameter (maxpt=75000)
      real m(1), r(1), x(1), y(1), mlim(4), dm(1), rlim(2)
      real*8 cxx(1),cyy(1),cxy(1),kron(1),aa(1),cstar(1)
      real abcissa(2,maxpt), ordinate(maxpt)
      real readin(20)
      real*8 coeff(2), cov(2,2)
      character*(*) fname

C Read data from FNAME, expecting MACH0 format: 
C    idx, x, y, r, type, N*(I,m,dm,col)
C Compute a fit to where the magnitude signal-to-noise is SNLIM
C as a function of r:  mcut = 2.5*(MLIM(1) - MLIM(2)*log(1+mlim(3)/r))
C Also cut out the brightest few points and return the cutoff in MLIM(4)
C Determine the inner and outer limiting radii and return as RLIM
C Finally return the data values as X, Y, R, M, and DM.

C      open(unit=8,file='testlkr.out',status='unknown')

      open(unit=10,file=fname,status='old')
      n = 0
      rlim(1) = 100000
      rlim(2) = 0
      read(10,*) i, x0, y0
C AJ -- include cxx,cyy,cxy and kron radius for masks
      do 10 i = 1,maxpt
         read(10,*,end=11) idx, x(i), y(i), r(i),
     1        itype, cxx(i), cyy(i), cxy(i),kron(i),aa(i), 
     2        (readin(j),j=1,4*ic)
         m(i) = readin(4*ic-2)
         dm(i) = readin(4*ic-1)
* jpb, also get class_star from last column
         cstar(i) = readin(4*ic)
         n = n + 1
         rlim(1) = amin1(rlim(1),r(i))
         rlim(2) = amax1(rlim(2),r(i))
C JPB, Feb/2009: get rid of this - mask really bright stars by hand
CC AJ: minimum kron radius in order to better mask bright stars
C         if (kron(i) < 5.) then
C            kron(i) = 5.
C         end if
C       WRITE(8,2533) I, X(I),cxx(I),cyy(I),cxy(I),kron(I),cstar(I)
C 2533    FORMAT(I5,3F8.1,3F8.2)
 10   continue
 11   continue
      close(10)

C Set the outer limit to miss the farthest corner(s)
      rlim(2) = rlim(2) / sqrt(2.)
C If mlim has been set, use its value; otherwise force it to 20
      if(mlim(4).lt.15.or.mlim(4).gt.25) mlim(4) = 20

C
C We expect that:
C (1) This is FAINT object photometry so the errors are independent of m
C (2) The errors are a constant flux at a given radius
C (3) The errors are the square root of a mean galaxy flux + sky flux
C (4) The galaxy falls off as a power law, 1/r or 1/r^2; pick 1/r
C Therefore we expect that:
C
C        df = (constant error flux) = K * sqrt(gxy/r + sky)
C        dm = df / f
C           = df * 10**(+0.4(m-m1))
C    log dm = log df + 0.4 m - 0.4 m1
C           = 0.4 m + C1 + C2 log (1+RSKY/r)
C
C RSKY is a non-linear parameter, so let's guess at it and then iterate.
C fit all the errors which are reasonable in terms of this
C linear function of m and r, determining C1 and C2.  Then a choice of
C limiting signal to noise gives us a function of mcut(r):
C     mlim(1) = alog10(1/snlim) - C1
C     mlim(2) = C2
C     mlim(3) = rsky
C     mcut = 2.5*(mlim(1) - mlim(2)*alog10(1+mlim(3)/r(i)))
C

      rsky = 50
      drsky = 10

      do 222 itry = 1,50
         rsky = rsky + drsky
         nfit = 0
         do 20 i = 1,n
            if(m(i).lt.mlim(4)) goto 20
* what the hell was this for?
*            if(dm(i).lt.0.015.or.dm(i).gt.0.30) goto 20
            if(dm(i).lt.0.001.or.dm(i).gt.2.5) goto 20
            if(r(i).lt.rlim(1).or.r(i).gt.rlim(2)) goto 20
            nfit = nfit + 1
            abcissa(1,nfit) = 1
            abcissa(2,nfit) = alog10(1+rsky/r(i))
            ordinate(nfit) = alog10(dm(i)) - 0.4*m(i)
 20      continue

         call linearfit(nfit,2,abcissa,ordinate,cov,coeff,rms)

         if(itry.eq.1) rmsmin = rms
C         WRITE(6,*) RSKY, RMS

         rmsmin = amin1(rmsmin,rms)
         if(rms.gt.rmsmin) then
C            WRITE(17,*) ITRY
C            DO 7761 I = 1,NFIT
C               WRITE(17,*) ABCISSA(1,I), ABCISSA(2,I), ORDINATE(I)
C 7761       CONTINUE
            goto 223
         end if

         mlim(1) = alog10(1/snlim) - coeff(1)
         mlim(2) = coeff(2)
         mlim(3) = rsky

 222  continue
 223  continue

C Now go through the data and remove all unwanted points
      k = 0
      do 30 i = 1,n
C Let's just remove points with no (zero) magnitude
         if(m(i).lt.10) goto 30
         k = k + 1
         x(k) = x(i)
         y(k) = y(i)
         r(k) = r(i)
         m(k) = m(i)
         dm(k) = dm(i)
 30   continue

      n = k

C      write(6,*) mlim
C      do 40 j = 10,1000,10
C         rr = float(j)
C         if(j.gt.200.and.mod(j,50).ne.0) goto 40
C         emcut = 2.5*(mlim(1) - mlim(2)*alog10(1+mlim(3)/rr))
C         em4 = emcut + 2.5*(alog10(snlim*0.25))
C         em10 = emcut + 2.5*(alog10(snlim*0.10))
C         write(6,*) rr, emcut, em4, em10
C 40   continue

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine stompout(nx,ny,bitmap,ix,iy,
     $                    cxx,cyy,cxy,lrad,aa)
      integer*2 bitmap(nx,ny)
      real*8 cxx,cyy,cxy,lrad,aa
      ircut = int(2.0*aa+0.99)
      do 10 j = max(0,iy-ircut), min(ny-1,iy+ircut)
         do 20 i = max(0,ix-ircut), min(nx-1,ix+ircut)
            ir2 = cxx*(real(i)-real(ix))**2 + 
     $            cyy*(real(j)-real(iy))**2
     $            + cxy*(real(i)-real(ix))*(real(j)-real(iy))
C            write (6,*) ix,iy,i,j,ir2,lrad*lrad,ircut,aa
C            radius2 = (real(i)-real(ix))**2 + (real(j)-real(iy))**2

C            if((ir2.le.lrad*lrad).or.(radius2.le.rcut*rcut)) then
            if(ir2.le.lrad*lrad) then
               bitmap(i+1,j+1) = 0
            end if
 20      continue
 10   continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine stompout_old(nx,ny,bitmap,ix,iy,rcut)
      integer*2 bitmap(nx,ny)
      ircut = int(rcut+0.99)
      do 10 j = max(0,iy-ircut), min(ny-1,iy+ircut)
         do 20 i = max(0,ix-ircut), min(nx-1,ix+ircut)
            ir2 = (i-ix)**2 + (j-iy)**2
            if(ir2.le.rcut*rcut) then
               bitmap(i+1,j+1) = 0
            end if
 20      continue
 10   continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine linearfit(npt,npar,x,y,cov,coeff,rms)
* Program to do a multi-parameter least-squares linear fit:
*    Y(i) = A1 x(1,i) + A2 x(2,i) + ... + An x(n,i)
      parameter (maxpar=20)
      implicit real*8 (a-h,o-z)
      real*4 x(npar,npt), y(npt), rms
      real*8 cov(npar,npar), vec(maxpar), coeff(npar)

      if(npar.gt.maxpar) then
         write(6,*) 'LINEARFIT: ',npar,' parameters?  Get real.'
         return
      end if

      if(npt.lt.npar) then
         write(6,*) 'LINEARFIT: ',npar,npt,' params/points?'
         return
      end if

      do 18 j = 1,npar
         vec(j) = 0
         do 19 i = 1,npar
            cov(i,j) = 0
 19      continue
 18   continue

      do 20 j = 1,npt
         do 21 k = 1,npar
            vec(k) = vec(k) + dble(y(j))*dble(x(k,j))
            do 22 i = 1,npar
               cov(i,k) = cov(i,k) + dble(x(i,j))*dble(x(k,j))
 22         continue
 21      continue
 20   continue

      call invert(npar,cov,coeff,det)
      if(det.eq.0) then
         write(6,*) 'LINEARFIT: singular matrix', npt, npar
         rms = 0
         return
      end if

      do 26 j = 1,npar
         coeff(j) = 0
         do 27 i = 1,npar
            coeff(j) = coeff(j) + vec(i) * cov(i,j)
 27      continue
         cov(j,j) = dsqrt(cov(j,j))
 26   continue

      do 28 j = 1,npar
         do 29 i = 1,npar
            if(i.eq.j) goto 29
            cov(i,j) = cov(i,j) / (cov(i,i)*cov(j,j))
 29      continue
 28   continue

      sum = 0
      ave = 0
      do 30 j = 1,npt
         fit = 0
         do 31 i = 1,npar
            fit = fit + coeff(i)*x(i,j)
 31      continue
         diff = y(j) - fit
         ave = ave + diff
         sum = sum + diff*diff
 30   continue
      ave = ave / npt
      rms = dsqrt(sum/npt - ave*ave)

 999  return
      end


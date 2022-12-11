! Program dophot
!
! This program performs PSF photometry in astronomical images. 
! For a complete and deep explanation of its main characteristics, 
! and to understand how it works, read the paper 
! "DoPHOT, A CCD Photometry Program: Description and Tests"
! (Schechter, Mateo, and Saha 1993 PASP, 105, 1345).

program dophot

! JT hackery to use openmp multi-threading...
!$  use omp_lib

  IMPLICIT REAL(8) (A-H,O-Z)
  integer, parameter :: npmax=7
  integer, parameter :: ntype=9
  integer, parameter :: nff=18
!  integer, parameter :: npap=4
!JT restore aperture magnitude errors
  integer, parameter :: npap=5
  integer, parameter :: maxfil=25000
  integer, parameter :: mthresh=50
  integer, parameter :: maxfit=10
  integer, parameter :: iadd=1,isub=-1
  real(8), parameter :: rmagic=1.0e10
  real(8) :: sig(3),acc(npmax),alim(npmax),ava(npmax)
  real(8) :: accskyp(npmax),alimskyp(npmax),accskyh(npmax),alimskyh(npmax)
  character(80) :: flags(nff),files(nff) 
  logical :: enuffvpsf
  logical :: fixpos,fixshape
  real(4), allocatable :: big(:,:),noise(:,:)
  character(80), allocatable :: headerim(:),headervar(:)
  real(8) :: threshold(mthresh), nsfrac
  real(8),allocatable :: skyplanepar(:),skyhubpar(:),skymedarray(:,:)
  real(8),allocatable :: a5(:),a6(:),a7(:)
  real(8),allocatable :: starpar(:,:), shadow(:,:),shaderr(:,:)
  integer,allocatable :: itype(:)
  real(8),allocatable :: appar(:,:)
  real(8),allocatable :: probg(:)
  real(8),allocatable :: ax(:),ay(:),axinv(:),ayinv(:)
  logical :: matched
  real(4) :: timearr(2)
  logical :: enuffap,applyap
  character(24) :: date1
  external pseud2d,pseud2d_two

  istime = time()                ! function time is not standard f90
  timediff = dtime(timearr)      ! function dtime is not standard f90

  call tuneup(npmax,nff,flags,files,skyguess,fitrad,maskrad,&
     aprad,apskymin,apskymax,eperdn,rnoise,bot,top,tmin,tmax,tfac,&
     nsmin,nsfrac,nsmax,nit,npstar,nfit0,nfit1,nfit2,nsrmsmin,fac,xpnd,&
     icrit,widobl,cmax,topsat,fsub,fobl,stograt,discrim,sig,chicrit,xtra,&
     crit7,snlim,bumpcrit,sncos,enuff4,enuff7,&
     nbadleft,nbadright,nbadtop,nbadbot,lverb,acc,alim,ava,&
     beta4,beta6,pixthresh,fixpos,fixshape,&
     napcmin,napertures,aperrmax,apmagmaxerr,&
     mapparmag,mapparxy,napsurffit,apfaintmag,&
     npskyp,npskyh,nsskypmin,nsskyhmin,iskymedstep,&
     accskyp,alimskyp,accskyh,alimskyh,itskyp,itskyh,&
     npsffit,npsfmin,ncmin,norder,&
     nsmatchmin,amagmaxm,amagminm,amagerrmaxm,xoffinit,yoffinit,&
     radius,dradius,radmin)

! For initalization
  allocate (starpar(npmax,nsmax),shadow(npmax,nsmax),shaderr(npmax,nsmax))
  allocate (itype(nsmax),appar(npap,nsmax),probg(nsmax))
  allocate (skyplanepar(npskyp),skyhubpar(npskyh))
  allocate (a5(npsffit),a6(npsffit),a7(npsffit))
  allocate (ax(maxfit),ay(maxfit),axinv(maxfit),ayinv(maxfit))
  nstot = 0
  nsperf = 0
  starpar = 0.
  shadow = 0.
  shaderr = 0.
  itype = 0
  appar = 0.
  probg = 0.
  skyplanepar = 0.
  skyhubpar = 0.
  a5 = 0.
  a6 = 0.
  a7 = 0.
  ax = 0.
  ay = 0.
  axinv = 0.
  ayinv = 0.
  enuffvpsf = .false.
  enuffap = .false.
  applyap = .false.
  matched = .false.

! JT if multithreaded, hard wire the max number of threads to use
! !$  call omp_set_num_threads( 8 )
! No, use "export OMP_NUM_THREADS=8" to set it from the shell

!  Open log file if desired.
  if((lverb.gt.0).and.(files(13)(1:4).ne.'TERM')) &
       call opena(13,files(13),2,ierr)
  lverb = lverb*10 + 1

!  Open input image and get its information
  call openimage(files(1),1,istat)
  call getheader1(1,nheaderim)
  allocate (headerim(nheaderim))
  call getheader2(1,nheaderim,headerim)
  call getdimen(1,nx,ny)
  write(*,'(a,2i6,2x,a)') ' From frame header: NX, NY = ',nx,ny,files(1)(:40)
  write(13,'(a,2i6,2x,a)') ' From frame header: NX, NY = ',nx,ny
  allocate (big(nx,ny))
  call getdatar4(1,big,nx,ny)
  call closeimage(1,istat)

! Allocate noise array
  allocate (noise(nx,ny))

! No variance image read: just fill noise array
  if(flags(18)(1:1).eq."N") then
! Zero noise array
     do i = 1,ny
        do j = 1,nx
           noise(j,i) = 0.0
        end do
     end do

! If requested read extra variance image
  else
     call openimage(files(18),18,istat)
     call getheader1(18,nheadervar)
     allocate (headervar(nheadervar))
     call getheader2(18,nheadervar,headervar)
     call getdimen(18,nvx,nvy)
     write(*,'(a,2i6,2x,a)') ' From variance header: NX, NY = ',&
          nvx,nvy,files(18)(:40)
     write(13,'(a,2i6,2x,a)') ' From variance header: NX, NY = ',nvx,nvy
     if(nx.ne.nvx .or. ny.ne.nvy) then
        write(6,*) 'Error: image and var have different sizes ',nx,ny,nvx,nvy
        call exit(1)
     end if
     call getdatar4(18,noise,nx,ny)
     call closeimage(18,istat)
  end if

  call makenoise(big,nx,ny,rnoise,eperdn,top,bot,&
       nbadtop,nbadbot,nbadleft,nbadright,rmagic,noise)

  if (flags(5)(1:1).eq.'Y') then
    call warmstart(big,noise,nx,ny,flags,files,nff,ava,npstar,nsmax,npap,&
         pseud2d,beta4,beta6,fobl,isub,fsub,xpnd,fac,rmagic,lverb,&
         nstot,starpar,itype,shadow,shaderr)
    if (fixpos.or.fixshape) nfix = nstot
  end if

  call thresholds(tmin,tmax,tfac,mthresh,threshold,nthresh)
  cmin = threshold(ncmin)

  allocate (skymedarray(nx/iskymedstep,ny/iskymedstep))
  skymedarray = 0.
 
!JT timing, JTtime>0 prints CPU time at each step of each iteration
  JTtime = 0
  if(JTtime .gt. 0) then
     call cpu_time(fin)
  end if
  do iloop=1,nthresh

     thresh = threshold(iloop)
     if(lverb.gt.10) then
        write(13,*)
        write(13,*) '  Starting loop at thresh = ',thresh
     end if

     if ((flags(17)(1:1).eq.'Y').and.(iloop.eq.1))  &
          call volcano(big,noise,nx,ny,rmagic,maxfil,fitrad,top,topsat)
!JT timing
     if(JTtime .gt. 0) then
        start = fin
        call cpu_time(fin)
        write(6,'(a40,i4,2f9.3)') 'volcano complete, iter', iloop, fin-start, fin
     end if

     call skypar(nstot,nx,ny,maxfil,nsmax,npstar,&
          big,noise,starpar,shadow,itype,ava,flags,nff,rmagic,lverb,&
          npskyp,npskyh,iskymedstep,&
          nsskypmin,itskyp,accskyp,alimskyp,nsskyhmin,itskyh,accskyh,alimskyh,&
          skyplanepar,skyhubpar,skymedarray)
!JT timing
     if(JTtime .gt. 0) then
        start = fin
        call cpu_time(fin)
        write(6,'(a40,i4,2f9.3)') 'skypar1 complete, iter', iloop, fin-start, fin
     end if

     if ((flags(5)(1:1).eq.'Y').and.(iloop.eq.1))  &
          call improve(big,noise,nx,ny,nstot,fixpos,fixshape,nfix,nfit0,nfit1,&
          nit,acc,alim,enuff4,enuff7,eperdn,rnoise,&
          apmagmaxerr,aprad,apskymin,apskymax,&
          nsperf,fsub,xpnd,fac,pseud2d,beta4,beta6,flags,nff,&
          ava,a5,a6,a7,npsffit,enuffvpsf,lverb,iadd,isub,&
          maxfil,snlim,fitrad,rmagic,cmin,crit7,nsmax,npstar,npap,&
          starpar,itype,appar,probg)

     call isearch(big,noise,nx,ny,thresh,pixthresh,rmagic,flags,nff,&
          skyplanepar,npskyp,skyhubpar,npskyh,iskymedstep,skymedarray,&
          enuffvpsf,a5,a6,a7,npsffit,ava,nsmax,npstar,&
          pseud2d,beta4,beta6,maskrad,&
          crit7,bumpcrit,lverb,&
          snlim,enuff4,fitrad,maxfil,nfit1,acc,alim,nit,&
          top,topsat,cmax,icrit,fobl,cmin,widobl,discrim,sncos,&
          fsub,xpnd,isub,fac,&
          nstot,starpar,itype,nsnew)
!JT timing
     if(JTtime .gt. 0) then
        start = fin
        call cpu_time(fin)
        write(6,'(a40,i4,2f9.3)') 'isearch complete, iter', iloop, fin-start, fin
     end if

     call shappe(big,noise,nx,ny,pseud2d,pseud2d_two,beta4,beta6,&
          fixpos,fixshape,nfix,nstot,nsnew,nsperf,nsrmsmin,nsmax,npstar,lverb,&
          flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,acc,alim,nit,&
          iadd,isub,fsub,xpnd,fac,rmagic,snlim,fitrad,maxfil,maxval,&
          enuff7,sig,xtra,chicrit,nfit2,stograt,&
          starpar,shadow,shaderr,itype)
!JT timing
     if(JTtime .gt. 0) then
        start = fin
        call cpu_time(fin)
        write(6,'(a40,i4,2f9.3)') 'shappe complete, iter', iloop, fin-start, fin
     end if

     if (iloop.gt.0) &
          call stparavg(nstot,npstar,nsmax,itype,shadow,shaderr,lverb,nsperf,ava)

     if (flags(11)(1:1).eq.'Y') &
          call vpsf (nsperf,nstot,starpar,itype,shadow,shaderr,ava,npstar,nsmax,&
          files,nff,npsffit,npsfmin,rmagic,lverb,enuffvpsf,a5,a6,a7)
!JT timing
     if(JTtime .gt. 0) then
        start = fin
        call cpu_time(fin)
        write(6,'(a40,i4,2f9.3)') 'vpsf complete, iter', iloop, fin-start, fin
     end if

     call skypar(nstot,nx,ny,maxfil,nsmax,npstar,&
          big,noise,starpar,shadow,itype,ava,flags,nff,rmagic,lverb,&
          npskyp,npskyh,iskymedstep,&
          nsskypmin,itskyp,accskyp,alimskyp,nsskyhmin,itskyh,accskyh,alimskyh,&
          skyplanepar,skyhubpar,skymedarray)
!JT timing
     if(JTtime .gt. 0) then
        start = fin
        call cpu_time(fin)
        write(6,'(a40,i4,2f9.3)') 'skypar2 complete, iter', iloop, fin-start, fin
     end if

     call improve(big,noise,nx,ny,nstot,fixpos,fixshape,nfix,nfit0,nfit1,&
          nit,acc,alim,enuff4,enuff7,eperdn,rnoise,&
          apmagmaxerr,aprad,apskymin,apskymax,&
          nsperf,fsub,xpnd,fac,pseud2d,beta4,beta6,flags,nff,&
          ava,a5,a6,a7,npsffit,enuffvpsf,lverb,iadd,isub,&
          maxfil,snlim,fitrad,rmagic,cmin,crit7,nsmax,npstar,npap,&
          starpar,itype,appar,probg)
!JT timing
     if(JTtime .gt. 0) then
        start = fin
        call cpu_time(fin)
        write(6,'(a40,i4,2f9.3)') 'improve complete, iter', iloop, fin-start, fin
     end if

!    If there was a warmstart with fixpos or fixshape, say so in the last loop
     if ((fixpos.or.fixshape).and.(iloop.eq.nthresh)) then
        do i=1,nfix
           itype(i)=itype(i)+10
        end do
     end if

!    Save the photometric information of the loop
     if(flags(4)(1:5).eq.'BINAR') then
        call openab(4,files(4),2,ierr,33)
     else
        call opena(4,files(4),2,ierr)
        if(flags(4)(1:5).eq.'DAOPH') call daoph_headinit(4,ntype)
     end if
     if (flags(8)(1:1).eq.'Y') call opena(8,files(8),2,ierr)
     do i=1,nstot
        call choose_out(4,flags(4),i,itype(i),starpar(1,i),npstar,&
             appar(1,i),npap,probg(i),applyap)
        if(flags(8)(1:1).eq.'Y') then
           if(shadow(1,i).ne.0) then
              call shadow_out(i,itype(i),shadow(1,i),npstar,8)
           else
              call shadow_out(i,itype(i),starpar(1,i),npstar,8)
           end if
        end if
     end do
     close(unit=4)
     if(flags(8)(1:1).eq.'Y') close(unit=8)

!    Save the substracted image
     if(flags(3)(1:1).eq.'Y') then
        call openimage(files(2),2,istat)
        call putheader(2,headerim,nheaderim)
        call putdatar4(2,big,nx,ny)
        close(unit=2)
     end if

!    Give the information about the loop
     if(lverb.gt.10) call notifyend(13,thresh,nsnew,nstot)
!     print*,'LOOP ',iloop,' of ',nthresh
!     call notifyend(6,thresh,nsnew,nstot)

     write(6,539) iloop, nthresh, thresh, nsnew, nstot
539  format('DoPhot iter', i3, ' of', i3, ' at thresh', f8.1, '   Nfound', i8, ' total', i8)

! JT Have we found enough stars that we can exit this loop?
     if(nstot .gt. nsmin) then
        write(6,'(a,i7,a,i7)') 'Quit: Ntotal', nstot, ' exceeds requirement', nsmin
        exit
     end if

! JT Has the increase been less than nsfrac?
     if(nsnew .lt. nsfrac*(nstot-nsnew)) then
        write(6,'(a,i7,a,f5.3,a,i7)') 'Quit: Nfound', nsnew, ' less than ', nsfrac, ' times', nstot
        exit
     end if

  end do
! Here we finish the last threshold

!  Aperture correction or aperture photometry.
  if ((flags(12)(1:1).eq.'Y').or.(flags(13)(1:1).eq.'Y')) &
       call apcorr2(big,noise,nx,ny,starpar,itype,npstar,nsmax,&
       appar,npap,pseud2d,beta4,beta6,nstot,flags,files,nff,nit,&
       mapparmag,mapparxy,napsurffit,apfaintmag,&
       napcmin,norder,nbadleft,nbadright,nbadtop,nbadbot,iadd,isub,&
       fsub,xpnd,fac,rmagic,lverb,eperdn,rnoise,&
       napertures,aprad,apskymin,apskymax,aperrmax,enuffap)

  if (flags(16)(1:1).eq.'N') &
       call apcorsat(big,noise,nx,ny,starpar,itype,npstar,nsmax,&
       appar,npap,pseud2d,beta4,beta6,nstot,flags,nff,&
       nbadleft,nbadright,nbadtop,nbadbot,iadd,isub,&
       fsub,xpnd,fac,rmagic,lverb,eperdn,rnoise,&
       aprad,apskymin,apskymax,aperrmax,nit,mapparmag)
     
   applyap = ((flags(12)(1:1).eq.'Y').and.enuffap).or.(flags(13)(1:1).eq.'Y')
   if(applyap) then
      if(flags(4)(1:5).eq.'BINAR') then
         call openab(4,files(4),2,ierr,33)
      else
         call opena(4,files(4),2,ierr)
         if(flags(4)(1:5).eq.'DAOPH') call daoph_headinit(4,ntype)
      end if
      do i=1,nstot
         call choose_out(4,flags(4),i,itype(i),starpar(1,i),npstar,&
              appar(1,i),npap,probg(i),applyap)
      end do
      close(unit=4)
   end if

! Automatch file.
  if((flags(9)(1:1).eq.'Y').and.(nsperf.gt.nsmatchmin)) then
     write(13,*) ' '
     write(13,*) ' Preparing for automatch'
     call findmatch2(starpar,itype,npstar,nsmax,nstot,npap,flags,files,nff,&
          amagminm,amagmaxm,amagerrmaxm,xoffinit,yoffinit,rmagic,lverb,&
          radius,dradius,radmin,nit,maxfit,norder,ax,ay,axinv,ayinv,matched)
     if(matched) then
        write(13,*) ' Automatch succeeded.'
     else
        write(13,*)
        write(13,*) ' WARNING:  Failed automatch; thresh = ', thresh
     end if
     if((flags(15)(1:1).eq.'Y').and.matched) then
       mfit=MGETMFIT(norder)
        write(13,*)
        write(13,*) '  WARNING! Applying transformation to output files! '
        write(13,*) 'norder,mfit = ',norder,mfit
        write(13,*) 'axinv = ',axinv
        write(13,*) 'ayinv = ',ayinv
        do i=1,nstot
           xorig = starpar(3,i)
           yorig = starpar(4,i)
           starpar(3,i) = SURFEVAL(xorig,yorig,axinv,mfit)
           starpar(4,i) = SURFEVAL(xorig,yorig,ayinv,mfit)
        end do
        if(flags(4)(1:5).eq.'BINAR') then
           call openab(4,files(4),2,ierr,33)
        else
           call opena(4,files(4),2,ierr)
           if(flags(4)(1:5).eq.'DAOPH') call daoph_headinit(4,ntype)
        end if
        do i=1,nstot
           call choose_out(4,flags(4),i,itype(i),starpar(1,i),npstar,&
                appar(1,i),npap,probg(i),applyap)
        end do
        close(unit=4)
     end if
  end if

  if(flags(4)(1:5).eq.'DAOPH') call daoph_final(4,ntype,files(4),nstot,npstar,ava)

! Save the information about the program run
  itimenow = time()    ! function time is an extension to f90
  elaptime = float(itimenow - istime)
  timediff = dtime(timearr)  ! function dtime is an extension to f90
  call fdate(date1) ! subroutine fdate is an extension to f90
  call finish(files,flags,nff,ava,npstar,axinv,ayinv,maxfit,itype,nstot,&
     norder,timediff,elaptime,date1)

  call exit(0)        ! subroutine exit is an extension to f90

end program dophot

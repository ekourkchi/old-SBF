! Subroutine isearch(big,noise,nfast,nslow,thresh,pixthresh,rmagic,&
!     flags,nff,skyplanepar,npskyp,skyhubpar,npskyh,iskymedstep,skymedarray,&
!     enuffvpsf,a5,a6,a7,npsffit,ava,nsmax,npstar,&
!     onestar,par1,par2,maskrad,&
!     crit7,bumpcrit,lverb,&
!     snlim,enuff4,fitrad,maxfil,nfit1,acc,alim,nit,&
!     top,topsat,cmax,icrit,fobl,cmin,widobl,discrim,sncos,&
!     fsub,xpnd,isub,fac,&
!     nstot,starpar,itype,nsnew)
! This subroutine looks for new objects in the image in a way explained in 
! paragraph 4 of section 3 of the DoPHOT paper (Schechter, Mateo, and Saha 1993
! PASP, 105, 1345).

subroutine isearch(big,noise,nfast,nslow,thresh,pixthresh,rmagic,&
     flags,nff,skyplanepar,npskyp,skyhubpar,npskyh,iskymedstep,skymedarray,&
     enuffvpsf,a5,a6,a7,npsffit,ava,nsmax,npstar,&
     onestar,par1,par2,maskrad,&
     crit7,bumpcrit,lverb,&
     snlim,enuff4,fitrad,maxfil,nfit1,acc,alim,nit,&
     top,topsat,cmax,icrit,fobl,cmin,widobl,discrim,sncos,&
     fsub,xpnd,isub,fac,&
     nstot,starpar,itype,nsnew)

! JT hackery to use openmp multi-threading...
!$  use omp_lib

  IMPLICIT REAL(8) (A-H,O-Z)
  real(4) :: big(nfast,nslow),noise(nfast,nslow)
  real(8) :: skymedarray(nfast/iskymedstep,nslow/iskymedstep)
  real(8) :: ava(npstar),skyplanepar(npskyp),skyhubpar(npskyh)
  real(8) :: starmask(-maskrad:maskrad,-maskrad:maskrad)
  real(8) :: a5(npsffit),a6(npsffit),a7(npsffit)
  real(8) :: a(npstar),err(npstar),acc(npstar),alim(npstar)
  real(8) :: xy(2,maxfil),z(maxfil),ze(maxfil)
  character(*) :: flags(nff)
  logical:: enuffvpsf,transmask,enuffpts,offpic
  logical :: cosm,toobright,toob,wipe,toofaint
  real(8) :: starpar(npstar,nsmax)
  integer :: itype(nsmax)

  integer :: it, ymin, ymax, nthread, stripass, stripwidth

  external onestar

! Number of threads, default is just one
  nthread = 1
! Number of threads if -fopenmp is enabled
!$ nthread = omp_get_max_threads()
  stripwidth = nslow / (2*nthread)

  nsprev = nstot

! Two passes with disjoint stripes in y prevents multiple threads from interfering
  do stripass=0,1

!$OMP PARALLEL DO&
!$OMP&   PRIVATE(IT, YMIN, YMAX, IX, IY, K, X, Y, BESTSKY, HIGHTHRESH, STARMASK, A,&
!$OMP&           NEWSTAR, XY,Z,ZE,NPT, AREAFILL, AREATOT, ENUFFPTS, COSM, ERR, CHI2)

! Work on each of the different stripes
     do it=0,nthread-1
        ymin = (2*it+stripass) * stripwidth
        ymax = ymin + stripwidth - 1
        if(it .eq. 0 .and. stripass .eq. 0) ymin = 1
        if(it .eq. nthread-1 .and. stripass .eq. 1) ymax = nslow

! Process all the stars; skip stars that are not in this stripe
        
        do iy = ymin,ymax
           y = float(iy)
           do ix = 1,nfast
              x = float(ix)
              if (noise(ix,iy).ge.amin1(rmagic,(thresh/pixthresh)**2)) cycle
              bestsky=STPAR1(x,y,flags,nff,skyplanepar,npskyp,skyhubpar,npskyh,&
                   skymedarray,nfast,nslow,iskymedstep)
              highthresh = bestsky+thresh
              if (big(ix,iy).lt.highthresh) cycle 
              a = 0.
              call makemask(x,y,flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,npstar,&
                   onestar,par1,par2,maskrad,starmask)
              if (.not.(TRANSMASK(x,y,bestsky,starmask,maskrad,&
                   big,noise,nfast,nslow,bumpcrit,rmagic,lverb))) cycle
              if(lverb.gt.20) &
                   write(13,*) ' Triggering on new object:',ix,iy,big(ix,iy)        
              call fillerup (x,y,big,noise,nfast,nslow,snlim,fitrad,rmagic,&
                   maxfil,npstar,xy,z,ze,npt,areafill,areatot,a)
              enuffpts = areafill.ge.enuff4*areatot
              if (.not.enuffpts) then
                 if(lverb.gt.20) write(13,*) 'Skipping: npt = ', npt
                 cycle
              end if
              if(.not.TRANSMASK(x,y,a(1),starmask,maskrad,&
                   big,noise,nfast,nslow,bumpcrit,rmagic,lverb)) then
                 if(lverb.gt.20) write(13,*) 'Failed transmask on average sky'
                 cycle
              end if
              call stpar567(x,y,flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,npstar,a)
              err=0.
              call chisq(onestar,par1,par2,xy,z,ze,npt,a,npstar,&
                   nfit1,acc,alim,nit,rmagic,lverb,err,chi2)
              if(chi2.ge.rmagic) then
                 if(lverb.gt.20) write(13,*) 'Failed to converge. chi2= ', chi2
                 cycle
              end if

!$OMP CRITICAL
              nstot = nstot+1
              newstar = nstot
              if(nstot.gt.nsmax) then
                 print *
                 print *,'The next star would overflow star arrays.'
                 print *,'Reductions complete to threshold = ',thresh
                 print *,'To continue, resize arrays and do warmstart.'
                 call exit(7)     ! subroutine exit is an ext to f90
              end if
              if(lverb.gt.20) write(13,*) 'This is star no. ',nstot
!$OMP END CRITICAL

              do k = 1,npstar
                 starpar(k,newstar) = a(k)
              end do

              if (OFFPIC(a(3),a(4),nfast,nslow)) itype(newstar) = 9

              call cosmic(big,noise,nfast,nslow,onestar,par1,par2,npstar,&
                   widobl,discrim,sncos,rmagic,lverb,&
                   starpar(1,newstar),cosm)
              toob = TOOBRIGHT(starpar(1,newstar),npstar,big,noise,nfast,nslow,&
                   top,fitrad,cmax,icrit,rmagic,lverb)
              if (cosm.or.(toob.and.(flags(16)(1:1).eq.'Y'))) then
                 itype(newstar) = 8
                 call oblit(starpar(1,newstar),npstar,'c',fobl,cosm,lverb,rmagic,&
                      nfast,nslow,big,noise)
                 cycle
              else if (toob.and.(flags(16)(1:1).eq.'N')) then
                 itype(newstar) = 9
                 call oblit2(starpar(1,newstar),npstar,topsat,fitrad,lverb,rmagic,&
                      nfast,nslow,big,noise)
              end if

              if (TOOFAINT(starpar(1,newstar),err,npstar,cmin,crit7,lverb)) &
                   itype(newstar) = 7   

              call addstar(isub,starpar(1,newstar),npstar,onestar,par1,par2,&
                   fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)

           end do
        end do

     end do    ! Stripe loop

!$OMP END PARALLEL DO

  end do    ! Alternate stripe loop

  
  nsnew = nstot - nsprev

end subroutine isearch







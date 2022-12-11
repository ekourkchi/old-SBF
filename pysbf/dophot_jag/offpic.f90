! Function offpic(a,npstar,ix,iy,nfast,nslow)
! This subroutine checks if center of the star is beyond edge of picture

function offpic(xcen,ycen,nfast,nslow)

  IMPLICIT REAL(8) (A-H,O-Z)
  logical :: offpic

  offpic = .false.
  if ((xcen.lt.1).or.(xcen.gt.nfast)) offpic = .true.
  if ((ycen.lt.1).or.(ycen.gt.nslow)) offpic = .true.

end function offpic


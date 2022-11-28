* A few character-oriented fortran utilities.
* John Tonry, 1980

* Length Non Blank of a string
      function lnb(c)
      character*(*) c
      lnb = len(c)
 10   if(lnb.gt.0.and.(c(lnb:lnb).eq.' '.or.c(lnb:lnb).eq.char(0))) then
         lnb = lnb - 1
         goto 10
      end if
      return
      end

* Write a floating number X into C with NSIG digits of precision
      subroutine fwrite(x,nsig,c)
      character*(*) c
      character*20 form, result
      l = len(c)
      write(form,1000) nsig
 1000 format('(e20.',i1,')')
      write(result,form) x
C Finish this someday when I need it more...
      return
      end

      subroutine parser(input,maxword,nword,length,word)
C Parse a string into words and return the lengths
      character*(*) input, word(maxword)
      integer length(maxword)
      logical blank, prev
      data itab /9/

      prev = .true.
      nword = 0
      do 10 i = 1,len(input)
         blank = input(i:i).eq.' ' .or. ichar(input(i:i)).eq.0 .or.
     $        ichar(input(i:i)).eq.itab
         if(.not.blank) then
            if(prev) then
               if(nword.ge.maxword) then
                  write(6,*) 'PARSE: Cannot parse so many words'
                  return
               end if
               nword = nword + 1
               length(nword) = 1
               word(nword) = input(i:i)
            else
               if(length(nword).ge.len(word(nword))) then
                  write(6,*) 'PARSE: Cannot parse such a long word'
               else
                  length(nword) = length(nword) + 1
                  word(nword)(length(nword):length(nword)) = input(i:i)
               end if
            end if
         end if
         prev = blank
 10   continue
      return
      end

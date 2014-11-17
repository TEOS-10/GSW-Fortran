!==========================================================================
subroutine indx(x,n,z,k)
!==========================================================================

!  Finds the index of the value in a monotonically increasing array
!
!  x	 :  array of monotonically increasing values
!  n     :  length of the array
!  z     :  value to be indexed
!
!  k      : index K :- if x(k) <= z < x(k+1), or
!               n-1 :- if z = x(n)
!

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n, k, ku, kl, km

real (r14), dimension(n) :: x
real (r14) :: z

if(z.gt.x(1).and.z.lt.x(n)) then
   kl=1
   ku=n
   do while (ku-kl.gt.1)
      km=0.5d0*(ku+kl)
      if(z.gt.x(km)) then
         kl=km
      else
         ku=km
      endif
   end do
   k=kl
   if(z.eq.x(k+1)) then 
     k = k+1
   end if
elseif (z.le.x(1)) then
      k = 1
elseif (z.ge.x(n)) then
      k = n-1
else
      write(*,*) 'ERROR in subroutine indx : out of range'
      write(*,*) 'z = ', z, 'n = ', n, 'x = ',x
end if

return
end subroutine

!--------------------------------------------------------------------------


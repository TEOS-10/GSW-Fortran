pure function gsw_util_sort_real (rarray) result(iarray)

! Sorts a real array and returns a list of sorted indices. The sort
! algorithm needs to be stable (should preserve the order of items of
! equal value).

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: rarray(:)    ! Values to be sorted

integer :: iarray(size(rarray))       ! Sorted ids

integer :: i, nx, depth_limit

nx = size(rarray)
depth_limit = 2*floor(log(nx*1.0_r8)/log(2.0_r8))

iarray = (/ (i, i=1,nx) /)
call introsort_index(iarray,1,nx,depth_limit)
return

contains

    !--------------------------------------------------------------------------

    pure recursive subroutine introsort_index (indx, lo, hi, depth_limit)

    ! An implementation of Introsort (also called introspective sort) - a
    ! sorting algorithm designed by David Musser in 1997. It is based on
    ! quicksort but introduces two optimizations to speed up the typical
    ! quicksort algorithm:
    ! - if the quicksort algorithm on any sub-array runs longer than expected
    !   (recursion is deeper than log(array size)) then switch to heapsort
    !   which is guaranteed to finish in O(n*log(n)) time;
    ! - if the number of elements is small enough then switch to an insertion
    !   sort which is faster on small data sets.
    ! These two optimizations allow introsort to outperform the “dumb”
    ! quicksort implementation on most arrays (even ~200 times on artificial
    ! “evil” arrays consisting of 100,000 elements).
    !
    ! Note 1: this version returns an array of sorted indices. So the sorted
    ! data is real_array(indx) - Glenn Hyland 20/8/2015
    !
    ! Note 2: this version of Introsort is stable (preserves the order of
    ! items of equal value) - Glenn Hyland 28/8/2015

    use gsw_mod_kinds

    implicit none

    integer, intent(inout) :: indx(:)
    integer, intent(in) :: lo, hi, depth_limit

    integer :: mid, p
    real (r8) :: pivot

    integer, parameter :: size_threshold = 16

    if (hi-lo .gt. size_threshold) then

        if (depth_limit .eq. 0) then
            call heapsort_index(indx(lo:hi))
            return
        end if

	mid = lo + (hi-lo)/2		   ! avoids potential overflow
        pivot = medianof3(indx,lo,mid,hi)
        call partition(indx,lo,hi,pivot,p)

        call introsort_index(indx,lo,p-1,depth_limit-1)
        call introsort_index(indx,p,hi,depth_limit-1)

    else

        call insertionsort_index(indx,lo,hi)

    end if
    return

    end subroutine introsort_index

    !--------------------------------------------------------------------------

    pure subroutine partition (indx, lo, hi, pivot, p)

    implicit none

    integer, intent(inout) :: indx(:)
    integer, intent(in) :: lo, hi
    real (r8), intent(in) :: pivot
    integer, intent(out) :: p

    integer :: left, right, t

    left = lo - 1
    right = hi + 1
    do while (left .lt. right)
        right = right - 1
        do while (pivot .lt. rarray(indx(right)))
            right = right - 1
        end do
        left = left + 1
        do while (rarray(indx(left)) .lt. pivot)
            left = left + 1
        end do
        if (left < right) then
            t = indx(left); indx(left) = indx(right); indx(right) = t
        end if
    end do
    if (left .eq. right) then
        p = left + 1
    else 
        p = left
    end if
    return
 
    end subroutine partition

    !--------------------------------------------------------------------------

    pure function medianof3 (indx, lo, mid, hi)

    implicit none

    integer, intent(in) :: indx(:), lo, mid, hi

    real (r8) :: medianof3

    if (rarray(indx(mid)) .lt. rarray(indx(lo))) then
        if (rarray(indx(hi)) .lt. rarray(indx(mid))) then
	      medianof3 = rarray(indx(mid))
        else
            if (rarray(indx(hi)) .lt. rarray(indx(lo))) then
                medianof3 = rarray(indx(hi))
            else
                medianof3 = rarray(indx(lo))
            end if
        end if
    else
        if (rarray(indx(hi)) .lt. rarray(indx(mid))) then
            if (rarray(indx(hi)) .lt. rarray(indx(lo))) then
                medianof3 = rarray(indx(lo))
            else
                medianof3 = rarray(indx(hi))
            end if
        else
            medianof3 = rarray(indx(mid))
        end if
    end if
    return

    end function medianof3

    !--------------------------------------------------------------------------

    pure subroutine insertionsort_index (indx, lo, hi)

    ! Insertion sort algorithm
    !
    ! Note: this version returns an array of sorted indices. So the sorted
    ! data is real_array(indx) - Glenn Hyland 20/8/2015

    implicit none

    integer, intent(inout) :: indx(:)
    integer, intent(in) :: lo, hi

    integer :: i, j, t

    do i = lo+1, hi
        j = i - 1
        t = indx(i)
        do while (j .ge. lo .and. rarray(indx(j)) .gt. rarray(t))
            indx(j+1) = indx(j)
            j = j - 1
        end do
        indx(j+1) = t
    end do
    return

    end subroutine insertionsort_index

    !--------------------------------------------------------------------------

    pure subroutine heapsort_index (indx)
 
    ! Heapsort algorithm
    !
    ! Note 1: this version returns an array of sorted indices. So the sorted
    ! data is real_array(indx) - Glenn Hyland 20/8/2015
    !
    ! Note 2: this version of Heapsort is stable (preserves the order of items
    ! of equal value) through the use of the lexicographic key comparison
    ! function "less_than" - Glenn Hyland 28/8/2015

    integer, intent(inout) :: indx(0:)

    integer :: start, n, bottom, t
 
    n = size(indx)
    do start = (n-2)/2, 0, -1
       call siftdown(indx,start,n);
    end do
 
    do bottom = n-1, 1, -1
       t = indx(0); indx(0) = indx(bottom); indx(bottom) = t;
       call siftdown(indx,0,bottom)
    end do
    return
 
    end subroutine heapsort_index
 
    !--------------------------------------------------------------------------

    pure subroutine siftdown (indx, start, bottom)
 
    integer, intent(inout) :: indx(0:)
    integer, intent(in) :: start, bottom

    integer :: child, root, t
 
    root = start
    do while(root*2 + 1 < bottom)

       child = root*2 + 1
       if (child+1 < bottom) then
          if (less_than(indx(child),indx(child+1))) child = child + 1
       end if
 
       if (less_than(indx(root),indx(child))) then
          t = indx(child); indx(child) = indx(root); indx(root) = t
          root = child
       else
          return
       end if  
    end do      
    return
 
    end subroutine siftdown
 
    !--------------------------------------------------------------------------

    pure logical function less_than (indx1, indx2)

    ! Lexicographic key comparison function to ensure sorting stability so
    ! the order of items of equal value is preserved.

    integer, intent(in) :: indx1, indx2

    if (rarray(indx1) .lt. rarray(indx2)) then
        less_than = .true.
    else if (rarray(indx1) .eq. rarray(indx2)) then
        if (indx1 .lt. indx2) then
            less_than = .true.
        else
            less_than = .false.
        end if
    else
        less_than = .false.
    end if
    return

    end function less_than

end function

!--------------------------------------------------------------------------

module gsw_mod_netcdf

use netcdf

use gsw_mod_kinds

integer, private :: ncid

public :: ncdf_open
public :: ncdf_close
public :: ncdf_get_dim
public :: ncdf_get_var
public :: ncdf_get_var_att

private :: ncdf_handle_err

contains

    !--------------------------------------------------------------------------

    subroutine ncdf_open (file_name)
    implicit none

    character (*), intent(in) :: file_name

    integer :: istat

    character (*), parameter :: fname = 'ncdf_open'

    istat = nf90_open(file_name, nf90_nowrite, ncid)
    if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                  'nf90_open',file_name)

    return
    end subroutine ncdf_open

    !--------------------------------------------------------------------------

    subroutine ncdf_close ()
    implicit none

    integer :: istat

    character (*), parameter :: fname = 'ncdf_close'

    istat = nf90_close(ncid)
    if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname,'nf90_close')

    return
    end subroutine ncdf_close

    !--------------------------------------------------------------------------

    integer function ncdf_get_dim (dim_name)
    implicit none

    character (*), intent(in) :: dim_name

    integer :: istat, dimid

    character (*), parameter :: fname = 'ncdf_get_dim'

    istat = nf90_inq_dimid(ncid, dim_name, dimid)
    if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                  'nf90_inq_dimid',dim_name)

    istat = nf90_inquire_dimension(ncid, dimid, len=ncdf_get_dim)
    if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                             'nf90_inquire_dimension',dim_name)
    return
    end function ncdf_get_dim

    !--------------------------------------------------------------------------

    subroutine ncdf_get_var (var_name, var0, var1, var2, var3)
    implicit none

    character (*), intent(in) :: var_name
    real (r8), intent(out), optional :: var0, var1(:), var2(:,:), var3(:,:,:)

    integer :: istat, varid

    character (*), parameter :: fname = 'ncdf_get_var'

    istat = nf90_inq_varid(ncid, var_name, varid)
    if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                  'nf90_inq_varid',var_name)

    if (present(var0)) then

       istat = nf90_get_var(ncid, varid, var0)
       if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                     'nf90_get_var',var_name)

    else if (present(var1)) then

       istat = nf90_get_var(ncid, varid, var1)
       if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                     'nf90_get_var',var_name)

    else if (present(var2)) then

       istat = nf90_get_var(ncid, varid, var2)
       if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                     'nf90_get_var',var_name)

    else if (present(var3)) then

       istat = nf90_get_var(ncid, varid, var3)
       if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                     'nf90_get_var',var_name)

    end if

    return
    end subroutine ncdf_get_var

    !--------------------------------------------------------------------------

    subroutine ncdf_get_var_att (var_name, var, att_name, att)
    implicit none

    character (*), intent(in) :: var_name, att_name
    real (r8), intent(out) :: var(:,:), att

    integer :: istat, varid

    character (*), parameter :: fname = 'ncdf_get_var_att'

    istat = nf90_inq_varid(ncid, var_name, varid)
    if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                  'nf90_inq_varid',var_name)

    istat = nf90_get_var(ncid, varid, var)
    if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                  'nf90_get_var',var_name)

    istat = nf90_get_att(ncid, varid, att_name, att)
    if (istat /= nf90_noerr) call ncdf_handle_err(istat,fname, &
                                                  'nf90_get_att',att_name)

    return
    end subroutine ncdf_get_var_att

    !--------------------------------------------------------------------------

    subroutine ncdf_handle_err (status, fname, nf90_name, var_name)
    implicit none

    integer, intent(in) :: status
    character (*), intent(in) :: fname, nf90_name
    character (*), intent(in), optional :: var_name

    if (status == nf90_noerr) return

    if (present(var_name)) then
        print*, '**Error return from ', nf90_name, ' (', var_name, ') in ', &
	        fname, ' ...'
    else
        print*, '**Error return from ', nf90_name, ' in ', fname, ' ...'
    end if
    print '(3x,a," (status=",i0,")")', trim(nf90_strerror(status)), status

    stop
    end subroutine ncdf_handle_err

end module gsw_mod_netcdf

!--------------------------------------------------------------------------

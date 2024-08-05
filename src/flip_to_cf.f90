module flip
    implicit none
    private
    public :: flip_to_cf
contains
    subroutine flip_to_cf(xin,xout,nlon,nlat,nz)

        ! Copied from flip_to_cf_and_box_smooth Feb 16, 2014

        implicit none

        integer, parameter :: RKIND = selected_real_kind(12)
        integer, intent(in) :: nlon,nlat,nz
        real (kind=RKIND), intent(in), dimension(nz,nlat,nlon) :: xin
        real (kind=RKIND), intent(out), dimension(nlon,nlat,nz) :: xout
 
        integer :: i,j,k

        do i=1,nlon
            do j=1,nlat
                do k=1,nz
                    xout(i,j,k) = xin(k,j,i)
                end do
            end do
        end do
 
    end subroutine flip_to_cf
end module flip

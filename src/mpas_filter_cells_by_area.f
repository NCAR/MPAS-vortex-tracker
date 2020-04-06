C NCLFORTSTART
      subroutine mpas_filter_cells_by_area( nEdgesOnCell,
     &                              cellsOnCell,areaCell,
     &                              maxEdges, nVertLevels, nCells,
     &                              field, radius_km )

      implicit none

      ! input

      integer, INTENT(IN) :: maxEdges,nCells, nVertLevels
      integer, INTENT(IN) :: nEdgesOnCell(nCells)
      integer, INTENT(IN) :: cellsOnCell(maxEdges,nCells)
      double precision, INTENT(IN) :: areaCell(nCells)
      real*8, INTENT(INOUT) ::  field(nVertLevels,nCells)
      real, INTENT(IN)      :: radius_km
C NCLEND

      !  local
      
      integer fpass,i,j, k
      integer fpasses(nCells)
      real pi
      parameter (pi=3.141592)
      integer nnbr   ! neighbor count
      real*8 nbr_tot ! sum of neighbor values
      ! tmp = holds new cell values for each smoothing pass
      real*8  tmp(nVertLevels,nCells)

      ! The number of smoothing passes varies for each cell.
      ! It is inversely proportional to the radius of each cell
      ! where the radius of the cell is sqrt(area/pi).
      ! Smaller cells undergo more smoothing passes while
      ! larger cells undergo fewer.
      ! This is not exactly the same as a circular filter with
      ! uniform or inverse-distance weighting, but it is similar.
      ! Using neighbors and smoothing passes may be harder to 
      ! explain but it is quicker than calculating distances from 
      ! every cell to every cell which a circular filter would
      ! require.
      !write(0,*)'Entered mpas_filter_by_area'


      ! From Apr 1, 2020 email from Michael Duda (via Wei Wang)
      ! The valid cellsOnCell indices are still 1-based; however,
      ! invalid cellsOnCell entries (i.e., those entries beyond
      ! nEdgeOnCell) were recently changed to 0. In other words, when
      ! looping through cellsOnCell for a particular cell, iCell, in the
      ! mesh, it's important to only loop from 1 to nEdgesOnCell(iCell)
      ! and not from 1 to maxEdges.

      ! Used to make sure cellsOnCell had values that ranged from 1 to
      ! nCells. But as of MPAS v7, zero is sometimes used for cell
      ! numbers. I don't know what problem this sanity check was trying
      ! to avoid, but I can't use it anymore. 
      if(maxval(cellsOnCell).gt.nCells
     c              .or.minval(cellsOnCell).lt.0) then
        write(0,*)'cellsOnCell(j,i)=',cellsOnCell(j,i)
        write(0,*)'but must be between 1 and ',nCells
        write(0,*)'minval(cellsOnCell)=',minval(cellsOnCell)
        write(0,*)'maxval(cellsOnCell)=',maxval(cellsOnCell)
        write(0,*)'stopping'
        stop
      endif

      do i=1,nCells
        fpasses(i) = radius_km/sqrt(1e-6*areaCell(i)/pi)
      end do
      !write(0,*)'mpas_filter_by_area: calculated fpasses'

      do fpass = 1, maxval(fpasses)
        ! Each smoothing pass is conservative.
        do k=1,nVertLevels

           do i=1,nCells
             tmp(k,i) = field(k,i)
             ! Smooth cell if fpass <= fpasses(i)
             ! Average it with neighbors that are also
             ! under fpass threshold.
             if(fpass <= fpasses(i))then
               nnbr = 0
               nbr_tot = 0.
               ! Average neighbors
               do j=1,nEdgesOnCell(i)
                ! is neighbor <= fpass threshold?
                if(fpass <= fpasses(cellsOnCell(j,i)))then
                  nnbr = nnbr + 1 ! increment neighbor count
                  nbr_tot = nbr_tot + field(k,cellsOnCell(j,i))
                endif
               end do
               ! Smoothed value is average of old value and 
               ! average of neighbors under fpass threshold.
               if(nnbr.gt.0) tmp(k,i)=0.5*(tmp(k,i) + nbr_tot/nnbr)
             endif
           end do !nCells

           do i=1,nCells
             field(k,i) = tmp(k,i)
           end do !nCells
        end do ! nVertLevels
      end do ! fpass
      return
      end

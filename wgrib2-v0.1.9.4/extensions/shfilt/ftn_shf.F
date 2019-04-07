#include "shfilt.h"

*.................................................................

      subroutine SPHFIL ( gfil, grid, lons, lats, nlon, nlat, 
     .                    n1, n2, grid_type, ier )
!
!     Spectral filter.
!
      implicit NONE

      integer   nlon, nlat
      integer   n1, n2

      real      grid(nlon,nlat)       ! input  grid
      real      gfil(nlon,nlat)       ! output grid
      real      lons(nlon)
      real      lats(nlat)
      integer   grid_type             ! LATLON or GAUSSIAN      

      integer   ier

!                                ---

      integer mtrunc, ntrunc, lwsha, lwshs, lwork, l1, l2, maxnt, nlon1
      logical  xwrap, xglob

      call shwrap_ ( lons, nlon, xwrap, xglob )
      if ( .not. xglob ) then
           ier = 98
           print *, ">>> spherepak requires global grids, but got"
           print *, ">>> lon(1) = ",lons(1),  
     .              "  ---  lon(nlon) = ",lons(nlon)
           return
      end if

      if ( xwrap ) then
           nlon1 = nlon - 1
      else         
           nlon1 = nlon
      end if

*     Local parameters defined from basic dimensions
*     -----------------------------------------------
      mtrunc = (nlon1+2)/2
      ntrunc = nlat
      l1     = nlat
      l2     = (nlat+1)/2
      maxnt  = 1
      if ( grid_type .eq. LATLON ) then
         lwsha  = 2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
         lwshs  = 2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
         lwork  = nlat*(maxnt*nlon+max0(3*l2,nlon))
      else if ( grid_type .eq. GAUSSIAN ) then
         lwsha  = nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
         lwshs  = nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
         lwork  = nlat*(nlon*maxnt+max0(3*l2,nlon))
      else
         ier = 97
         return
      end if

      call SPHFIL_ ( gfil, grid, nlon, nlat, n1, n2, ier,
     .               mtrunc, ntrunc, lwsha, lwshs, lwork,
     .               xwrap, grid_type )

      end

*.................................................................

      subroutine SPHFIL_ ( gfil, grid, nlon, nlat, n1, n2, ier,
     .                     mtrunc, ntrunc, lwsha, lwshs, lwork,
     .                     xwrap, grid_type )

      implicit NONE

      integer   nlon, nlat
      integer   n1, n2
      integer   mtrunc, ntrunc, lwsha, lwshs, lwork 
      logical   xwrap

      real      grid(nlon,nlat)       ! input  grid
      real      gfil(nlon,nlat)       ! output grid
      integer   grid_type             ! LATLON or GAUSSIAN
      
      integer   ier

c
c     This routine applies a spectral filter retaining total wavenumbers
c     in the range [n1,n2].
c
C  ON INPUT:
c
C     grid(nlon,nlat):       array(s) to be filtered
C
c
C  ON OUTPUT:
c
C     gfil(nlon,nlat):       filtered array (may share storage with grid)
C   
c  WORK SPACE:
c
c                           Allocated internally
c
c
c  NOTE 1:
c  ------
c
C       For the triangular truncation that SPHEREPACK uses,
C  the even and odd spectral coefficients are arranged as follows:
C
C         **********
C         *********0
C         ********00
C         *******000
C         ******0000    * => NONZERO VALUE
C      N  *****00000
C         ****000000
C         ***0000000
C         **00000000
C         *000000000
C
C             M
C
c   m=0 to mtrunc (the 2-delta-x wave)
c   n=m to ntrunc-1
c   the latitudinal wavenumber (number of zeroes) is n-m=0 to ntrunc-m-1
c
c   The even and odd spectral coefficients associated with a particular
c m and n are stored in even(m+1,n+1) and odd(m+1,n+1). IMASK should
c be set appropriately.
c
c   Select resolutions by editing file 'spec.h'. It is assumed a GEOS
c kind of grid: longitudes from -180 to 180 and latitudes from -90
c to 90 (poles are grid points). Proper reordering to comply with
c SPHEREPACK conventions is performed in this routine.
c
c NOTE 2:
c ------
c
c   The original version of this routine allowed to perform analysis
c of several levels at a time. This can be easily restored.
c
c.......................................................................


*     Local work space
*     ----------------
      real      even(mtrunc,ntrunc) ! even spectral coefficients
      real      odd(mtrunc,ntrunc)    ! odd  spectral coefficients
      real      g(nlat,nlon)

      real   wsha(lwsha)
      real   wshs(lwshs)
      real*8 work(lwork)

      integer isym, nt, i, j, m, n, nlon1

*     Costant to determine whether is necessary to initialize
*     -------------------------------------------------------
      integer lastlon, lastlat
      save lastlon, lastlat
      data lastlon, lastlat / -1, -1 /

!ams      print *, 'inside SHPFIL'
!ams      print *, 'nlon, nlat = ', nlon, nlat
!ams      print *, 'mtrunc, ntrunc = ', mtrunc, ntrunc
!ams      print *, 'lwsha, lwork = ', lwsha, lwork 

!     Is it x-wrapped?
!     ----------------
      if ( xwrap ) then
           nlon1 = nlon - 1
      else         
           nlon1 = nlon
      end if

*     Initialize spherical harmonics package
*     --------------------------------------
!!!      if ( nlon .ne. lastlon .or. nlat .ne. lastlat ) then
      if ( grid_type .eq. LATLON ) then
         call SHAECI ( nlat, nlon1, wsha, lwsha, work, lwork, ier )
         if ( ier .ne. 0 ) return
         call SHSECI ( nlat, nlon1, wshs, lwshs, work, lwork, ier)
         if ( ier .ne. 0 ) return
      else if ( grid_type .eq. GAUSSIAN ) then
         call SHAGCI ( nlat, nlon1, wsha, lwsha, work, lwork, ier )
         if ( ier .ne. 0 ) return
         call SHSGCI ( nlat, nlon1, wshs, lwshs, work, lwork, ier)
         if ( ier .ne. 0 ) return
      else
         ier = 97
         return
      end if
!!!         lastlon = nlon
!!!         lastlat = nlat
!!!      endif

*     Reorder array from GrADS to spherepak ordering
*     ----------------------------------------------
      do i=1,nlat
         do j=1,nlon1/2
            g(i,j)=grid(j+nlon1/2,nlat-i+1)
         end do
         do j=nlon1/2+1,nlon1
            g(i,j)=grid(j-nlon1/2,nlat-i+1)
         end do
      end do

*     Calculate spectral transform
*     ----------------------------
      isym = 0
      nt = 1
      if ( grid_type .eq. LATLON ) then
         call SHAEC ( nlat, nlon1, isym, nt, 
     .                g, nlat, nlon, 
     .                even, odd, mtrunc, ntrunc, 
     .                wsha, lwsha, work, lwork,ier )
      else
         call SHAGC ( nlat, nlon1, isym, nt, 
     .                g, nlat, nlon, 
     .                even, odd, mtrunc, ntrunc, 
     .                wsha, lwsha, work, lwork,ier )
      end if
      if ( ier .ne. 0 ) return


!     Chop spectral coefficients
!     --------------------------
      do m = 1, mtrunc
         do n = m, ntrunc
            if ( n .lt. n1 .OR. n .gt. n2 ) then
                 even(m,n) = 0.0
                 odd(m,n)  = 0.0
             end if
          end do
      end do

!     Reconstruct filtered field
!     --------------------------
      if ( grid_type .eq. LATLON ) then
         call SHSEC ( nlat, nlon1, isym, nt, g, nlat, nlon, 
     .                even, odd, mtrunc, ntrunc, 
     .                wshs, lwshs, work, lwork, ier )
      else 
         call SHSGC ( nlat, nlon1, isym, nt, g, nlat, nlon, 
     .                even, odd, mtrunc, ntrunc, 
     .                wshs, lwshs, work, lwork, ier )
      end if

*     Reorder array back to GrADS ordering
*     ------------------------------------
      do i=1,nlat
         do j=1,nlon1/2
            gfil(j+nlon1/2,nlat-i+1) = g(i,j)
         end do
         do j=nlon1/2+1,nlon1
            gfil(j-nlon1/2,nlat-i+1) = g(i,j)
         end do
         if ( xwrap ) gfil(nlon,i) = gfil(1,i)
      end do

*     All done.
*     --------
      return
      end

*.................................................................

      subroutine SPHPWR ( gpwr, grid, lons, lats, nlon, nlat,
     .                    grid_type, ier)
!
!     Computes the power spectra
!
      implicit NONE

      integer   nlon, nlat

      real      grid(nlon,nlat)       ! input  grid
      real      gpwr(nlon,nlat)       ! output grid
      real      lons(nlon)
      real      lats(nlat)
      integer   grid_type             ! LATLON or GAUSSIAN

      integer   ier

!                                ---

      integer mtrunc, ntrunc, lwsha, lwshs, lwork, l1, l2, maxnt, nlon1
      logical xwrap, xglob

      call shwrap_ ( lons, nlon, xwrap, xglob )
      if ( .not. xglob ) then
           ier = 98
           print *, ">>> spherepak requires global grids, but got"
           print *, ">>> lon(1) = ",lons(1),  "lon(nlon) = ",lons(nlon)
           return
      end if

      if ( xwrap ) then
           nlon1 = nlon - 1
      else         
           nlon1 = nlon
      end if

*     Local parameters defined from basic dimensions
*     -----------------------------------------------
      mtrunc = (nlon1+2)/2
      ntrunc = nlat
      l1     = mtrunc
      l2     = (nlat+1)/2
      maxnt  = 1
      if ( grid_type .eq. LATLON ) then
         lwsha  = 2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15
         lwork  = nlat*(maxnt*nlon+max0(3*l2,nlon))
      else if ( grid_type .eq. GAUSSIAN ) then
         lwsha  = nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
         lwork  = nlat*(nlon*maxnt+max0(3*l2,nlon))
      else
         ier = 97
         return
      end if

      call SPHPWR_ ( gpwr, grid, nlon, nlat, ier,
     .               mtrunc, ntrunc, lwsha, lwork, xwrap, grid_type )

      end

      subroutine SPHPWR_ ( gpwr, grid, nlon, nlat, ier,
     .                     mtrunc, ntrunc, lwsha, lwork, xwrap,
     .                     grid_type )

      implicit NONE

      integer   nlon, nlat
      integer   mtrunc, ntrunc, lwsha, lwshs, lwork 
      logical   xwrap

      real      grid(nlon,nlat)       ! input  grid
      real      gpwr(nlon,nlat)       ! output grid
      integer   grid_type             ! LATLON or GAUSSIAN
      
      integer   ier

c
c     This routine applies a spectral filter retaining total wavenumbers
c     in the range [n1,n2].
c
C  ON INPUT:
c
C     grid(nlon,nlat):       array(s) to be filtered
C
c
C  ON OUTPUT:
c
C     gpwr(nlon,nlat):       power spectra (may share storage with grid)
C   
c  WORK SPACE:
c
c                           Allocated internally
c
c
c  NOTE 1:
c  ------
c
C       For the triangular truncation that SPHEREPACK uses,
C  the even and odd spectral coefficients are arranged as follows:
C
C         **********
C         *********0
C         ********00
C         *******000
C         ******0000    * => NONZERO VALUE
C      N  *****00000
C         ****000000
C         ***0000000
C         **00000000
C         *000000000
C
C             M
C
c   m=0 to mtrunc (the 2-delta-x wave)
c   n=m to ntrunc-1
c   the latitudinal wavenumber (number of zeroes) is n-m=0 to ntrunc-m-1
c
c   The even and odd spectral coefficients associated with a particular
c m and n are stored in even(m+1,n+1) and odd(m+1,n+1). IMASK should
c be set appropriately.
c
c   Select resolutions by editing file 'spec.h'. It is assumed a GEOS
c kind of grid: longitudes from -180 to 180 and latitudes from -90
c to 90 (poles are grid points). Proper reordering to comply with
c SPHEREPACK conventions is performed in this routine.
c
c NOTE 2:
c ------
c
c   The original version of this routine allowed to perform analysis
c of several levels at a time. This can be easily restored.
c
c.......................................................................


*     Local work space
*     ----------------
      real      even(mtrunc,ntrunc) ! even spectral coefficients
      real      odd(mtrunc,ntrunc)    ! odd  spectral coefficients
      real      power(mtrunc,ntrunc)      ! power spectra
      real      spectra(ntrunc)           ! power per total wavenumber
      real      g(nlat,nlon)

      real      wsha(lwsha)
      real*8    work(lwork)

      integer isym, nt, i, j, m, n, nlon1

*     Costant to determine whether is necessary to initialize
*     -------------------------------------------------------
      integer lastlon, lastlat
      save lastlon, lastlat
      data lastlon, lastlat / -1, -1 /

!ams      print *, 'inside SHPFIL'
!ams      print *, 'nlon, nlat = ', nlon, nlat
!ams      print *, 'mtrunc, ntrunc = ', mtrunc, ntrunc
!ams      print *, 'lwsha, lwork = ', lwsha, lwork 

      if ( xwrap ) then
           nlon1 = nlon - 1
      else         
           nlon1 = nlon
      end if

*     Initialize spherical harmonics package
*     --------------------------------------
!!!      if ( nlon .ne. lastlon .or. nlat .ne. lastlat ) then
      if ( grid_type .eq. LATLON ) then
         call SHAECI ( nlat, nlon1, wsha, lwsha, work, lwork, ier )
         if ( ier .ne. 0 ) return
      else if ( grid_type .eq. GAUSSIAN ) then
         call SHAGCI ( nlat, nlon1, wsha, lwsha, work, lwork, ier )
         if ( ier .ne. 0 ) return
      else
         ier = 97
         return
      end if
!!!      lastlon = nlon
!!!      lastlat = nlat
!!!      endif

*     Reorder array from GrADS to spherepak ordering
*     ----------------------------------------------
      do i=1,nlat
         do j=1,nlon1/2
            g(i,j)=grid(j+nlon1/2,nlat-i+1)
         end do
         do j=nlon1/2+1,nlon1
            g(i,j)=grid(j-nlon1/2,nlat-i+1)
         end do
      end do

*     Calculate spectral transform
*     ----------------------------
      isym = 0
      nt = 1
      if ( grid_type .eq. LATLON ) then
         call SHAEC ( nlat, nlon1, isym, nt, 
     .                g, nlat, nlon, 
     .                even, odd, mtrunc, ntrunc, 
     .                wsha, lwsha, work, lwork,ier )
         if ( ier .ne. 0 ) return
      else
         call SHAGC ( nlat, nlon1, isym, nt, 
     .                g, nlat, nlon, 
     .                even, odd, mtrunc, ntrunc, 
     .                wsha, lwsha, work, lwork,ier )
         if ( ier .ne. 0 ) return
      end if

*     Initialize output to zero
*     -------------------------
      do 100 n = 1, ntrunc
         spectra(n) = 0.
         do 100 m = 1, mtrunc
            power(m,n) = 0.
 100  continue

*     Compute power in terms of wavenumbers
*     -------------------------------------
      do 20 m = 1, mtrunc
         do 20 n = m, ntrunc
            power(m,n)=((even(m,n)**2+odd(m,n)**2)/4.)
 20   continue

*     Power as function of total wavenumber
*     -------------------------------------
      do 800 m = 1, mtrunc
         do 800 n = m, ntrunc
            spectra(n) = spectra(n) + power(m,n)
 800  continue

*     Return power at x=1
*     -------------------
      do j = 1, nlat
         gpwr(1,j) = spectra(j)
      end do

*     All done.
*     --------
      return
      end

*.....................................................................

      subroutine shwrap_ ( lons, nlon, xwrap, xglob )
      implicit NONE
      integer  nlon
      real     lons(nlon)
      logical  xwrap, xglob
      real x(2), y(2), pi, d2r, dlon, dist
!                         ----
      xwrap = .false. ! are longitudes wrapped?
      xglob = .false. ! is this a full circle?

!     cartesian coordinates on the unit spehere, at equator
      pi  =  4.0 * atan(1.0)
      d2r =  pi / 180.
      dlon = 2. * pi / (nlon-1) ! assumes it is wrapped first
      x(1) = cos(d2r*lons(1))
      y(1) = sin(d2r*lons(1))
      x(2) = cos(d2r*lons(nlon))
      y(2) = sin(d2r*lons(nlon))
      dist = sqrt( (x(1)-x(2))**2 + (y(1)-y(2))**2 )
!     Is first longitude the same as the last?
      if ( dist .lt. 0.01 * dlon ) then
           xwrap = .true.
           xglob = .true.
      else ! it may still be global
         dlon = 2. * pi / nlon
         x(2) = cos(d2r*lons(nlon)+dlon)
         y(2) = sin(d2r*lons(nlon)+dlon)
         dist = sqrt( (x(1)-x(2))**2 + (y(1)-y(2))**2 )
         if ( dist .lt. 0.01 * dlon ) xglob = .true.
      end if
!ams      print *, "---- Wrapped? ", xwrap
!ams      print *, "---- Global? ", xglob
      end


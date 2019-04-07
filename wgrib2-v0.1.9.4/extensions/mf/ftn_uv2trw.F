      subroutine ftnuv2trw(u, v, lon, lat, dum, undef, 
     $     ni, nj, centlon,centlat, 
     $     iout, rc)

      USE gex
      USE trig_vals

      real(kind=iGaKind) ::  u(ni,nj),v(ni,nj)
      real(kind=iGaKind) ::  dum(ni,nj)
      real(kind=iGaKind) ::  lon(ni),lat(nj)
      
      integer rc


C         
C initialize constants         
C
      pi=4.0*atan(1.0)
      dtr = pi/180.0
C        
C         initialize dum to input field
C
C
      do i=1,ni
        xlon=lon(i)

        do j=1,nj

          xlat=lat(j)
          udat=u(i,j)
          vdat=v(i,j)

          call getvrvt (centlon,centlat,xlon,xlat
     &         ,udat,vdat,vr,vt,igvtret)

          if(iout == 1) dum(i,j)=vt
          if(iout == -1) dum(i,j)=vr
CCCC          print*,'qqqq ',i,j,xlon,xlat,udat,vdat,vr,vt

        enddo
      enddo
      
      do i=1,ni
        do j=1,nj
          u(i,j)=dum(i,j)
        end do
      end do

      rc=0
      return
      end



c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      subroutine getvrvt (centlon,centlat,xlon,xlat
     &                   ,udat,vdat,vr,vt,igvtret)
c
c     ABSTRACT: This subroutine takes as input a u-wind and v-wind value
c     at an input (xlon,xlat) location and returns the tangential and
c     radial wind components relative to the input center lat/lon 
c     position (centlon,centlat).  The only trick to this whole 
c     subroutine is figuring out the angle from the center point to the
c     data point, and we do this by creating a triangle with the leg 
c     from the center point to the data point being the hypotenuse.
c
c     NOTE: All longitudes must be in positive degrees east (0-360) !!!
c
c     INPUT:
c     centlon  Longitude of center point
c     centlat  Latitude  of center point
c     xlon     Longitude of pt at which vr & vt will be computed
c     xlat     Latitude  of pt at which vr & vt will be computed
c     udat     u-value of wind at the point (xlon,xlat) 
c     vdat     v-value of wind at the point (xlon,xlat) 
c
c     OUTPUT:
c     vr      Radial wind component at (xlon,xlat) wrt (centlon,centlat)
c     vt      Tang   wind component at (xlon,xlat) wrt (centlon,centlat)
c     igvtret Return code from this subroutine

      USE trig_vals

      logical verb

      real(kind=iGaKind) :: centlon,centlat,xlon,xlat,udat,vdat,vr,vt

      verb=.false.

      call calcdist(centlon,centlat,xlon,xlat,hyp_dist)

      xlatdiff = abs(centlat - xlat)
      xlondiff = abs(centlon - xlon)

      if (xlondiff == 0 .and. xlatdiff > 0) then

        if (centlat > xlat) angle = 180   ! pt directly south of ctr
        if (centlat < xlat) angle = 0     ! pt directly north of ctr

      else if (xlondiff > 0 .and. xlatdiff == 0) then

        if (centlon > xlon) angle = 270   ! pt directly west of ctr
        if (centlon < xlon) angle = 90    ! pt directly east of ctr
C         
C         exactly on the center
C
      else if (hyp_dist ==  0) then

        vr=0.0
        vt=0.0
        return

      else

        ! This next part figures out the angle from the center point
        ! (centlon,centlat) to the data point (xlon,xlat).  It does 
        ! this by setting up a triangle and then using inverse trig
        ! functions to get the angle.  Since this is a kludgy way to
        ! do it that doesn't account for the curvature of the earth,
        ! we'll do it 2 ways, using asin and then acos, then take the
        ! average of those 2 for the angle.  hyp_dist, calculated just
        ! above, is the distance from the center pt to the data pt.

        opp_dist  = xlatdiff/360. * ecircum
        sin_value = opp_dist / hyp_dist
        if (sin_value > 1.0.and.verb) then
          print *,' '
          print *,'!!! In getvrvt, sin_value > 1, setting to 1.'
          print *,'!!! opp_dist= ',opp_dist,' hyp_dist= ',hyp_dist
          print *,'!!! sin_value = ',sin_value
          print *,'!!! centlon= ',centlon,' centlat= ',centlat
          print *,'!!! xlon=    ',xlon,' xlat=    ',xlat
          print *,' '
          sin_value = 0.99999
        endif
        sin_angle = asin(sin_value) / dtr

        call calcdist(centlon,centlat,xlon,centlat,adj_dist)
        cos_value = adj_dist / hyp_dist
        if (cos_value > 1.0.and.verb) then
          print *,' '
          print *,'!!! In getvrvt, cos_value > 1, setting to 1.'
          print *,'!!! adj_dist= ',adj_dist,' hyp_dist= ',hyp_dist
          print *,'!!! cos_value = ',cos_value
          print *,'!!! centlon= ',centlon,' centlat= ',centlat
          print *,'!!! xlon=    ',xlon,' xlat=    ',xlat
          print *,' '
          cos_value = 0.99999
        endif
        cos_angle = acos(cos_value) / dtr

        tmpangle = 0.5 * (sin_angle + cos_angle)

        ! The previous lines of code just calculated an angle between
        ! 0 and 90.  This next if structure adjusts that angle to 
        ! instead be between 0 and 360.

        if      (centlat <= xlat .and. centlon <= xlon) then
          angle = 90 - tmpangle
        else if (centlat >  xlat .and. centlon <= xlon) then
          angle = 90 + tmpangle
        else if (centlat >= xlat .and. centlon >= xlon) then
          angle = 270 - tmpangle
        else if (centlat <  xlat .and. centlon >= xlon) then
          angle = 270 + tmpangle
        endif

      endif

      uvrcomp = udat * sin(angle * dtr)
      vvrcomp = vdat * cos(angle * dtr)
      vr      = uvrcomp + vvrcomp

      uvtcomp = (-udat) * cos(angle * dtr)
      vvtcomp = vdat    * sin(angle * dtr)
      vt      = uvtcomp + vvtcomp

      return 
      end

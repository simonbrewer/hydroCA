!-------------------------------------------------------------------------------
! WCA2D: Weighted 2D Cellular Automata for modeling floods
!
! Guidolin et al. (2016). A weighted cellular automata 2D inundation model 
! for rapid flood analysis Env. Mod. Soft., 84, 378-394
!
! Ver. 0.5 Removed version number from file name
!          Made some changes to the velocity routine 
!          calculation of alpha now includes hydraulic radius
!
! Ver. 0.4 Add subroutines for variable 
!
! Ver. 0.3 Add subroutines for variable time step updating
!          TBD
!      
! Ver. 0.2 Split the code into subroutines
!          calc_flows: calculates all inflows and outflows
!          update_depth: updates water depth in all cells
!-------------------------------------------------------------------------------

      subroutine wca2d( m, n, dem, ppt, wse, dt, nullcell,
     >                  mannn, cellx, cellem, cella, nbit)
      !-------------------------------------------------------------------------
      ! Input variables
      integer m,n ! Grid sizes
      double precision dem( m, n ) ! Elevation values (m)
      double precision ppt( m, n ) ! PPT grid values (mm)
      double precision dt ! Current time step (s)
      double precision nullcell ! Value used to indicate null/border cells
      double precision mannn ! Manning's roughness coefficient
      double precision cellx ! Distance between cell centers (HC)
      double precision cellem ! Cell edge length (HC)
      double precision cella ! Cell area (HC)
      integer nbit ! Number of iterations

      ! Output variables
      double precision wse( m, n ) ! PPT grid values (mm)
      double precision fluxes( m, n, 4 ) ! Flux of water between cells

      ! Internal variables
      integer i,j,k
      double precision itot( m, n ) ! Cuml. inflow into each cell (m3)
      double precision otot( m, n ) ! Total outflow from each cell (m3)
      
      !-------------------------------------------------------------------------
      ! Iteration loop (currently no stopping condition)

      do 10 i=1,nbit
      !-------------------------------------------------------------------------
      ! Main loop through grid cells
      
      call calc_flows( m, n, dem, ppt, wse, otot, itot, 
     >                 fluxes, dt, nullcell, mannn, 
     >                 cellx, cellem, cella)

      call update_depth( m, n, ppt, wse, itot, otot, dt, cella )

10    continue

      end

      subroutine wca2d_1t( m, n, dem, ppt, wse, otot, itot,
     >                     dt, dtu,
     >                     currt, nullcell,
     >                     mannn, cellx, cellem, cella, 
     >                     vcheck, newdt, vamax,
     >                     tolwd, tolslope)
      !-------------------------------------------------------------------------
      ! Input variables
      integer m,n ! Grid sizes
      double precision dem( m, n ) ! Elevation values (m)
      double precision ppt( m, n ) ! PPT grid values (mm)
      double precision dt ! Current time step (s)
      double precision dtu ! Current update time step (s)
      double precision currt ! Current time (s)
      double precision nullcell ! Value used to indicate null/border cells
      double precision mannn ! Manning's roughness coefficient
      double precision cellx ! Distance between cell centers (HC)
      double precision cellem ! Cell edge length (HC)
      double precision cella ! Cell area (HC)
      double precision tolwd ! Minimum water depth for flow to occur
      double precision tolslope ! Minimum slope for flow to occur
      integer vcheck ! Number of iterations

      ! Output variables
      double precision wse( m, n ) ! water surface elevation (m)
      double precision fluxes( m, n, 4 ) ! Fluxes to neighboring cells
      double precision newdt ! Estimated new dt
      double precision vamax ! Estimated max velocity (m s-1)

      ! Internal variables
      double precision itot( m, n ) ! Cuml. inflow into each cell (m3)
      double precision otot( m, n ) ! Total outflow from each cell (m3)
      double precision sarray( m, n ) ! Grid of vector speeds (m s-1)
      double precision aarray( m, n ) ! Grid of vector angles (radians)
      double precision dtarray( m, n ) ! Grid of minimum dt values (s)
      
      !-------------------------------------------------------------------------
      ! Main loop through grid cells
      
      call calc_flows( m, n, dem, ppt, wse, otot, itot, 
     >                 fluxes, dt, nullcell, mannn, 
     >                 cellx, cellem, cella)

      call update_depth( m, n, ppt, wse, itot, otot, dt, cella )

      !if (mod(currt,dtu).eq.0.0) then
      if (vcheck.eq.1) then
      !write(*,*) currt, "Velocity check"
      call velocity ( m, n, dem, wse, itot, fluxes, dt,
     >                sarray, aarray, dtarray, 
     >                cellx, cellem, cella, nullcell, 
     >                mannn, tolwd, tolslope )
      newdt = minval(dtarray)
      vamax = maxval(sarray)

      endif

10    continue

      end

      !-------------------------------------------------------------------------
      ! Subroutine to calculate flows between all cells
      subroutine calc_flows( m, n, dem, ppt, wse, otot, itot, 
     >                       fluxes, dt, nullcell, mannn, 
     >                       cellx, cellem, cella)

      !-------------------------------------------------------------------------
      ! Input variables
      integer m,n ! Grid sizes
      double precision dem( m, n ) ! Elevation values (m)
      double precision ppt( m, n ) ! PPT grid values (mm)
      double precision dt ! Current time step (s)
      double precision nullcell ! Value used to indicate null/border cells
      double precision mannn ! Manning's roughness coefficient
      double precision cellx ! Distance between cell centers (HC)
      double precision cellem ! Cell edge length (HC)
      double precision cella ! Cell area (HC)

      ! Output variables
      double precision wse( m, n ) ! PPT grid values (mm)
      double precision fluxes( m, n, 4 ) ! Flux of water between cells

      ! Internal variables
      integer i,j,k
      double precision taut ! Threshold for flow to occur
      integer offx(4),offy(4)
      double precision dl0i(4),dl0i1(4)
      double precision dV0i(4)
      double precision dVmin,dVmax,dVtot
      double precision wi(4) ! Weights
      double precision wi0 ! Weights
      integer maxid(1) ! Weights
      double precision d0,sd0g,nd0g ! Water depth (m)
      double precision vmax ! Maximum velocity (m/s)
      double precision im ! Max. volume to neighbor (m3)
      double precision g ! Gravitational acceleration (m/s-2)
      double precision itot( m, n ) ! Cuml. inflow into each cell (m3)
      double precision otot( m, n ) ! Total outflow from each cell (m3)
      double precision itotdt

      ! Set up parameters
      parameter(taut = 1e-16)
      parameter(g = 9.81) ! (m/s-2)
      offx = (/-1,+1,0,0/) ! Offsets for von Neumann neighborhood
      offy = (/0,0,-1,+1/)

      !write(*,*) offx

      !-------------------------------------------------------------------------
      ! Main loop through grid cells
      
      ! Reset all inflows to zero
      itot(:,:) = 0
      ! Reset all outflows to zero
      otot(:,:) = 0
        
      do 10 i=2,(m-1)
      
      do 20 j=2,(n-1)

      if (dem(i,j).ne.nullcell) then ! Check for barrier cells
              ! write(*,*) i,j,dem(i,j),ppt(i,j),itot(i,j),wse(i,j)
      if (wse(i,j).gt.dem(i,j)) then ! Check for standing water in cell
      ! Get neighbors
      do 30 k=1,4
      ! Estimate difference in level
      dl0i(k) = ( wse(i,j) + tau ) - wse((i+offx(k)),(j+offy(k)))
      !write(*,*) 1,k,dl0i(k)

      ! Set negative level differences to 0
      !dl0i(k) = dmax1( myz, dl0i1(k) )
      if (dl0i(k).le.0.0) then
        dl0i(k)=0.0
      endif
      !write(*,*) 2,k,dl0i(k)
      dV0i(k) = max( 0.0, dl0i(k)) * cella ! Equation 2 
      !write(*,*) 3,k,dl0i(k),dV0i(k)

30    continue      
      
      ! Volume changes
      dVmax = maxval(dV0i)
      dVmin = minval(dV0i, mask=dV0i.gt.0.0)
      dVtot = sum(dV0i)
      if (maxval(dV0i).gt.0.0) then
        !write (*,*) dVmin,dVmax,dVtot,"DOSTUFF"

        ! Weights (eqn. 6; NEED TO CHECK PRECISION ERROR HERE, wi > 1)
        wi(:) = dV0i(:) / (dVtot + dVmin)
        w0 = dVmin / (dVtot + dVmin)
        maxid = maxloc(wi)
        
        !write(*,*) "wi",w0,maxid,wi

        ! Maximum permissable velocity (m/s eqn. 9)
        d0 = wse(i,j) - dem(i,j) ! Water depth in central cell (m)
        sd0g = sqrt(d0 * g) ! 
        nd0g = (1/mannn) * d0**(2./3.) * sqrt(dl0i(maxid(1)) / cellx)
        vmax = min(sd0g, nd0g)

        !write(*,*) "d0",d0,sd0g,nd0g,vmax

        ! Maximum volume to neighbor w/ highest weight (m3, eqn. 10)
        im = vmax * d0 * dt * cellem

        !write(*,*) "im",im

        ! Total volume to leave cell (m3, eqn. 11)
        itotdt = min( d0*cella, 
     >              im/maxval(wi), 
     >              dVmin + otot(i,j) ) 

        !write(*,*) "d0*cella",d0*cella
        !write(*,*) "im/maxwi",im/maxval(wi)
        !write(*,*) "dvmin+otot",dVmin+otot(j,k) 
        !write(*,*) i,j,"itotdt",itotdt
        ! Copy outflow value to otot (could really just store this as
        ! otot)

        otot(i,j) = itotdt

        ! Calculate flow into all neighboring cells (eqn. 12)
        do 41 k=1,4
        fluxes(i,j,k) = itotdt * wi(k)
        itot((i+offx(k)),(j+offy(k))) = itot((i+offx(k)),(j+offy(k))) +
     >    itotdt * wi(k)
41      continue      
        ! Add "inflow" for volume that is retained in center cell
        itot(i,j) = itot(i,j) + w0 * itotdt

        
      endif ! Check for positive outflow volume
              
      endif ! Standing water check

      endif ! Barrier cell check

20    continue      
10    continue    
      
      end

      !-------------------------------------------------------------------------
      ! Subroutine to update water levels
      ! This is equation 13 from Guidolin et al
      subroutine update_depth( m, n, ppt, wse, itot, otot, dt, cella )

      !-------------------------------------------------------------------------
      ! Input variables
      integer m,n ! Grid sizes
      double precision ppt( m, n ) ! PPT grid values (mm)
      double precision itot( m, n ) ! Cuml. inflow into each cell (m3)
      double precision otot( m, n ) ! Total outflow from each cell (m3)
      double precision dt ! Current time step (s)
      double precision cella ! Cell area (HC)

      ! Output variables
      double precision wse( m, n ) ! PPT grid values (mm)

      ! Internal variables
      integer i,j

      do 10 i=2,(m-1)
      
      do 20 j=2,(n-1)
      wse(i,j) = wse(i,j) + 
     >           itot(i,j) / cella +
     >           (( ppt(i,j) * 1e-3 ) / ( 60 * 60 )) * dt -
     >           otot(i,j) / cella

20    continue      
10    continue      

      end

      !-------------------------------------------------------------------------
      ! Subroutine to calculate velocity
      ! Most of this was taken from the CADDIES code
      ! in velocityDiffusive.ca
      subroutine velocity( m, n, dem, wse, itot, fluxes, dt,
     >                     sarray, aarray, dtarray, 
     >                     cellx, cellem, cella, nullcell, 
     >                     mannn, tolwd, tolslope )

      !-------------------------------------------------------------------------
      ! Input variables
      integer m,n ! Grid sizes
      double precision dem( m, n ) ! PPT grid values (mm)
      double precision wse( m, n ) ! PPT grid values (mm)
      double precision itot( m, n ) ! Cuml. inflow into each cell (m3)
      double precision fluxes( m, n, 4 ) ! Flux of water between cells
      double precision dt ! Current time step (s)
      double precision sarray( m, n ) ! Grid of vector speeds (m s-1)
      double precision aarray( m, n ) ! Grid of vector angles (radians)
      double precision dtarray( m, n ) ! Grid of minimum dt values (s)
      double precision cellx ! Distance between cell centers (HC)
      double precision cellem ! Cell edge length (HC)
      double precision cella ! Cell area (HC)
      double precision nullcell ! Value used to indicate null/border cells
      double precision mannn ! Manning's roughness coefficient

      ! Internal variables
      integer i,j,k
      double precision d0 ! Average water depth for neighboring cells
      double precision vh ! default velocity (m/s)
      double precision x,y ! x and y of velocity vector
      double precision alpha ! Alpha of time step (from code)
      double precision prevdt ! Previous dt
      double precision dwse ! Difference in water elevation
      double precision hr ! Hydraulic radius (m)
      double precision sg ! Hydraulic gradient (-)
      double precision ndem,nwse
      double precision tolwd ! Minimum water depth for flow to occur
      double precision tolslope ! Minimum slope for flow to occur
      double precision wdmain ! Water depth of central cell
      double precision nwd ! Water depth of neighbor cell
      double precision pi ! 3.1415...
      double precision angle(4) ! Angle betw. cell and neighbors (radians)
      integer offx(4),offy(4)

      !parameter(tolwd = 1e-4,tolslope = 1e-4)
      parameter(pi = 4.0*atan(1.0))
      !angle = (/0.0, pi/2.0, pi, 3.0*pi/2.0/)
      angle(1) = 0.0
      angle(2) = pi/2.0
      angle(3) = pi
      angle(4) = 3.0*pi/2.0
      offx = (/0,+1,0,-1/) ! Offsets for von Neumann neighborhood
      offy = (/+1,0,-1,0/)

      prevdt = dt ! Keep previous dt value (delete?)
      sarray(:,:) = 0.0 ! Initialize speed array
      aarray(:,:) = 0.0 ! Initialize angle array
      dtarray(:,:) = 1000.0 ! Initialize possible dt array

      ! Loop through cells
      do 10 i=2,(m-1)
      
      do 20 j=2,(n-1)

      x = 0.0 ! Set vector coordinates to zero
      y = 0.0 ! Set vector coordinates to zero 
      alpha = 1000.0 ! Alpha set to large value
      wdmain = wse(i,j) - dem(i,j)

      ! Loop through neighbors
      do 30 k=1,4
      d0 = 0.0 ! Set difference in depths to zero
      vh = 0.0 ! Set default velocity to zero

      ! Get values for neighboring cell
      ndem = dem((i+offx(k)),(j+offy(k))) 
      nwse = wse((i+offx(k)),(j+offy(k))) 
      nwd = nwse - ndem

      ! Check for border cells
      if (ndem.eq.nullcell) then
        dwse = -1.0e6
        hr = 0.0
      else
        dwse = abs( wse(i,j) - nwse )
        hr = max( wse(i,j), nwse ) - max( dem(i,j), ndem )
      endif

      if ((fluxes(i,j,k).gt.0) .and. (hr.gt.tolwd)) then
        ! Manning and critical velocity
        sg = dwse / cellx ! Hydraulic gradient (S)
        vh = fluxes(i,j,k) / (wdmain * cellem * prevdt) 
        if (sg.gt.tolslope) then
          !write(*,*) "Updating DT",dwse,hr
          !write(*,*) "Updating DT",sg,vh
          !write(*,*) "Updating DT",mannn,hr**(5.0/3.0),sqrt(sg),
        !>                 ((2*mannn) / (dwse**(5.0/3.0))) *
        !>                 sqrt(sg)
          ! Equation 17 (modified Manning's av. velocity) 
          alpha = min( alpha, ((2*mannn) / (hr**(5.0/3.0))) *
     >                 sqrt(sg) )
        endif
        x = x + vh * cos(angle(k)) ! eqn. 15
        y = y + vh * sin(angle(k)) ! eqn. 15
      endif
      
30    continue
      ! If there is some velocity
      if (abs(x).gt.0.0001 .or. abs(y).gt.0.0001) then
        sarray(i,j) = sqrt(x**2 + y**2) ! eqn. 16
        aarray(i,j) = atan2(y,x) ! eqn. 16
      endif
      ! Calculate possible new dt (eqn. 17)
      dtarray(i,j)= min( dtarray(i,j), (cella / 4) * alpha)
      !possdt = min( possdt, (cella / 4) * alpha)
      !write(*,*) prevdt,dt,alpha
      !write(*,*) "possdt",possdt,(cella/4)*alpha

20    continue      
10    continue
      !write(*,*) dtarray      

      end

!
! User module
!
! This module contains functions that may be altered by the user to define a specific problem. 
! 
! These functions are called by tristanmainloop.F90 and should provide
! initialization of parameters specific to the problem being setup, 
! loading of initial plasma, injection of new plasma, and boundary conditions
! on particles and fields that are specific to user's problem.
! 
! Functions by name:
! read_input_user() -- parse the problem-specific section of the input file for user parameters
! get_external_fields() -- called by mover to provide externally imposed EM fields if external_fields is on in input file
! init_EMfields_user() -- initializes EM fields 
! init_particle_distribution_user() -- loads initial particle distribution
! inject_particles_user() -- injects new particles on every time step as needed
! field_bc_user() -- applies user-specific boundary conditions to fields (beyond periodic or radiation BCs)
! particle_bc_user() -- enforces user-specific BCs on particles (e.g., reflecting walls)
! shift_domain_user() -- needed for some shock problems to remove empty space; ignore. 
!
#ifdef twoD 

module m_user

	use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_particles
	use m_inputparser
	use m_fparser
	use m_domain
	
#else

module m_user_3d

	use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_particles_3d
	use m_inputparser_3d 
	use m_fparser_3d
	use m_domain_3d

#endif

	implicit none
		
	private

!-------------------------------------------------------------------------------
!	PARAMETERS
!-------------------------------------------------------------------------------

	
!-------------------------------------------------------------------------------
!	TYPE DEFINITIONS
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!	VARIABLES
!-------------------------------------------------------------------------------

	real(sprec) :: temperature_ratio, sigma_ext, bext0

        real :: betainj,left_wall_speedup 
	integer :: betainj_interval, rightclean, rightwall         
        integer :: rightclean_interval, rightclean_start


!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	public :: init_EMfields_user, init_particle_distribution_user, &
	inject_particles_user, read_input_user, field_bc_user, get_external_fields, &
	particle_bc_user, shift_domain_user

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

!
! External variables read from the input file earlier; they are available for use here
!
! sigma (electron sigma),
! Binit (B field based on sigma); can be reset here
! me, mi, qe, qi (mass and charge of the species, abs(qe)/me = 1
! ppc0 (fiducial particles per cell for ions and electrons together)
! btheta, bhi -- inclination angles of initial B field, in degrees
! delgam -- k T/ m_e c^2 -- fiducial normalized temperature of electrons 
!
! gamma0, and its beta -- 4-velocity of the upstream flow. 
! mx0,my0,mz0 -- real number of cells in each direction, including the ghost zones. Active domain is from 3 to mx0-3, inclusive, and same for y and z. This is global size of the grid.
! each core knows also about mx, my, mz, which is the size of the local grid on this core. This size includes 5 ghost zones. 
!
! iglob, jglob, kglob -- functions that return global i,j,k based on local i,j,k
! xglob, yglob, zglob -- functions that return global x,y,z based on local x,y,z

	contains

!-------------------------------------------------------------------------------
! 						subroutine read_input_user		
!									
! Reads any variables related to (or needed by) this module from section "problem" in input file
! 							
!-------------------------------------------------------------------------------

subroutine read_input_user()

	implicit none
	integer :: lextflds, luserpartbcs, lwall, ldim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CHANGE THE NAME ON THIS LINE IF YOU ARE CREATING A NEW USER FILE
!This helps to identify which user file is being used. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(rank.eq.0) print *, "Using user file user_shock_pu.F90"

!inputpar_getd_def -- read real parameter; input_geti_def -- read integer parameter

	call inputpar_getd_def("problem", "temperature_ratio", 1._sprec, Temperature_ratio)

	call inputpar_getd_def("problem","sigma_ext",0._sprec,sigma_ext)

	call inputpar_geti_def("problem","external_fields",0,lextflds)

	if(lextflds .eq. 1) then 
	   external_fields =.true.
	else
	   external_fields =.false.
	endif

	if(external_fields) bext0 = sqrt((gamma0)*.5*ppc0*c**2*me*(1+me/mi)*sigma_ext)

	call inputpar_geti_def("problem","user_part_bcs",0,luserpartbcs)

	if(luserpartbcs .eq. 1) then 
	   user_part_bcs = .true.
	else
	   user_part_bcs = .false.
	endif

	call inputpar_geti_def("problem", "wall", 0, lwall)
	
	if (lwall==1) then
		wall=.true.
	else
		wall=.false.
	endif
	
	if(wall) user_part_bcs=.true.

!read_input_user is called last, so any parameter can be overwritten here

       	call inputpar_getd_def("problem","wallgam",1._sprec,wallgam)

	if(wallgam .eq. 0.) wallgam=1.	
	if(wallgam<1) wallgam=sign(1./sqrt(1-wallgam**2),wallgam)

        call inputpar_getd_def("problem", "left_wall_speedup",1._sprec,left_wall_speedup) !extra speed multiplier with which left wall will effectively move; wallgam will set the speed for reflecting particles off the wall; this is to eliminate downstream in upstream frame

	call inputpar_getd_def("problem", "betainj", .99999, betainj)
	call inputpar_geti_def("problem", "betainj_interval",1000, betainj_interval) !how often to check text file in output directory for updates on right wall velocity

 	call inputpar_geti_def("problem", "rightclean",0, rightclean) 
	call inputpar_geti_def("problem", "rightclean_interval",300, rightclean_interval) 
	call inputpar_geti_def("problem", "rightclean_start",500, rightclean_start) 
	call inputpar_geti_def("problem", "rightwall",0, rightwall) 

	call inputpar_geti_def("problem", "distr_dim", 2, ldim)
	if (ldim == 2) then 
	   pcosthmult = 0.
	else if (ldim == 3) then 
	   pcosthmult = 1.
	else
	   pcosthmult = 0.
	endif
	
end subroutine read_input_user

!-------------------------------------------------------------
!     Compute external fields to be added to the mover. 
!     These fields do not evolve via Maxwell Eqs, but can depend on time
!-------------------------------------------------------------
	subroutine get_external_fields(x,y,z,ex_ext, ey_ext, ez_ext, bx_ext,by_ext,bz_ext, qm, n)
	
	real,intent(inout):: bx_ext,by_ext,bz_ext, ex_ext, ey_ext, ez_ext
	real, intent(in):: x,y,z
	real, optional:: qm
        integer, optional:: n
	ex_ext=0.
	ey_ext=0.
	ez_ext=0.
	bx_ext=0.
	by_ext=0.
	bz_ext=0.
        ! can use bext0 as fiducial field.
        !external field can be directed in any way; only invoked for external_fields = 1
        !x,y,z come in as local coordinates, convert to global by calling xglob(x), yglob(y), zglob(z)	
	end subroutine get_external_fields

!-------------------------------------------------------------------------------
! 						subroutine init_EMfields_user		 
!												
! Sets the electromagnetic fields of any specific user purpose
!							
!-------------------------------------------------------------------------------

subroutine init_EMfields_user()
	
	! local variables
	
	integer :: i, j, k

! sigma is defined as
! sqrt(sigma)=omega_c/omega_p = (qe B / gamma0 me c)/sqrt( qe 0.5 ppc0 (1+me/mi)/gamma0)  
! where 4 pi=1, qe/me = 1 
! B = sqrt(gamma0*.5*ppc0*(1+me/mi)*c**2*(me)*sigma)
! (1+me/mi) term comes from omega_p^2 = omega_pe^2+omega_pi^2

       	Binit=sqrt(gamma0*ppc0*.5*c**2*(me*(1+me/mi))*sigma)
	!initialize B field to be set by Binit and the inclination angle 	
        !angles btheta and bphi are read earlier

	if(Binit .ne. 0 .and. pcosthmult .eq. 0) then 
	   pcosthmult = 1.
	   if(rank .eq. 0) print *, "Warning: Resetting dimensionality of particle", &
	   " distribution to 3D for magnetized plasma"
	endif


	btheta=btheta/180.*pi
	bphi=bphi/180.*pi

	do  k=1,mz
		do  j=1,my
			do  i=1,mx
! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
				
				bx(i,j,k)=Binit*cos(btheta) 
				by(i,j,k)=Binit*sin(btheta)*sin(bphi)
				bz(i,j,k)=Binit*sin(btheta)*cos(bphi)

				ex(i,j,k)=0.
				ey(i,j,k)=(-beta)*bz(i,j,k) 
				ez(i,j,k)=-(-beta)*by(i,j,k)

			enddo
		enddo
	enddo

end subroutine init_EMfields_user


!-------------------------------------------------------------------------------
! 						subroutine init_particle_distribution_user()	
!											
! Sets the particle distrubtion for a user defined case
!
!-------------------------------------------------------------------------------

subroutine init_particle_distribution_user()

	implicit none

	! local variables
	
	integer :: i, n, direction, ierr
	real    gamma_drift, delgam_i, delgam_e, ppc, weight
	character (len=20) :: dummychar
	real(sprec) :: x1,x2,y1,y2,z1,z2


	!set initial injection points 
 
        leftwall=20.
	xinject=leftwall
	xinject2=(mx0-2.)/2 
        walloc = leftwall
	! ------------------------------------------

!initialize left going upstream

	gamma_drift=-gamma0 ! negative gamma_drift will send the plasma in the negative direction
	delgam_i=delgam  !delgam is read from input file in particles
	delgam_e=delgam*mi/me*Temperature_ratio

	x1=xinject  
	x2=xinject2 

	y1=3. !in global coordinates
	y2=my0-2.  
	z1=3.
	z2=mz0-2. !if 2D, it will reset automatically
	ppc=ppc0
	weight=1.

	direction=1  !drift along x, + or - determined by the sign of gamma_drift
 
	call inject_plasma_region(x1,x2,y1,y2,z1,z2,ppc,&
             gamma_drift,delgam_i,delgam_e,weight,direction)

	x1in=3 !set the location of planes where particles are removed from 
               !simulation, perpendicular to x. 
	x2in=mx0-2

	if(irestart .eq. 1) then 
	   open(unit=20,file="./restart/adjust_param.txt",action="read")
	   read (20,*) dummychar,betainj
	endif

	if(rank .eq. 0) then
	    open(unit=10,file=trim(output_dir_name)//"/adjust_param.txt",action="write",status="replace")
	    write (10,*) "betainj=",betainj
            write (10,*) "left_wall_speedup=",left_wall_speedup

	    close (10)
	endif
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	call check_overflow()
	call reorder_particles()
	
end subroutine init_particle_distribution_user


!-------------------------------------------------------------------------------
! 				subroutine inject_particles_user()					 
!										
! Injects new particles in the simulation on every step. To skip injection, set ppc=0 below
!
!-------------------------------------------------------------------------------

subroutine inject_particles_user()

	implicit none
	real(sprec) :: x1,x2,y1,y2,z1,z2
	real delgam_i, delgam_e, injector_speed, ppc, gamma_drift, weight
	character (len=20) :: dummychar
        integer ierr, n

	injectedions=0 !set counter to 0 here, so that all possible other injections on this step add up
	injectedlecs=0

 	if(rightclean .eq. 1) then
	     xinject2=xinject2+c*.99999 
	else
	   xinject2=xinject2+c*betainj !for when rightclean=0 and the injector moves slowly
	endif

	if(xinject2 .gt. mx0-15) then !stop expansion of the injector has reached the end of the domain
	   xinject2=mx0-15.
	   betainj=0. !injector hit the wall, stop moving the injector
	endif

 	if(rightclean .eq. 1) then 
	   injector_speed = .99999
	else
	   injector_speed = betainj
	endif

        x1=xinject2
        x2=x1
	y1=3. !in global coordinates
	y2=my0-2.  
	z1=3.
	z2=mz0-2. !if 2D, it will reset automatically
        weight=1
        ppc=ppc0

	gamma_drift=-gamma0
	delgam_i=delgam
	delgam_e=delgam*mi/me*Temperature_ratio

	call inject_from_wall(x1,x2,y1,y2,z1,z2,ppc,gamma_drift,&
	  delgam_i,delgam_e,injector_speed,weight)
!inject_from_wall sends a stream of e and ions along normal to the wall 
!which has x1=x2, or y1=y2, etc. This is how it finds which direction to send the stream
 
    x1in=3.
    if(rightclean .eq. 1) then
       x2in=mx0-2.	              
    else
       x2in=xinject2  !for cases when there is no right wall and we want CRs to leave through the right injector
       if(rightwall >0) x2in=xinject2+10 !not sure about +10, but setting to xinject2 was killing subsonic electrons, 
    					!not allowing them to scatter off the right wall. 
       !Is this still necessary?
    endif

    if(rightclean .eq. 1 .and. lap .gt. rightclean_start &
         .and. modulo(lap, rightclean_interval) .eq. 0) then

       xinject2=(xinject2-rightclean_interval*c*(.99999 - betainj))
       x2in=xinject2+c*.99999 !betainj*c

        !remove particles beyond x2in
        if(1>0) then 
        n=1
        do while (n <= ions)
           if(p(n)%x + mxcum > x2in) then 
              call copyprt(p(ions),p(n))
              ions=ions-1
              n=n-1
           endif
           n=n+1
        enddo

        n=1
        do while (n <= lecs)
           if(p(n+maxhlf)%x + mxcum > x2in) then 
              call copyprt(p(maxhlf+lecs),p(maxhlf+n))
              lecs=lecs-1
              n=n-1
           endif
           n=n+1
        enddo
endif !if 1<0

    endif

	if(modulo(lap,betainj_interval) .eq. 0) then
	    open(unit=20,file=trim(output_dir_name)//"/adjust_param.txt",action="read")
	    read (20,*) dummychar,betainj
            read (20,*) dummychar,left_wall_speedup

	    if(rank .eq. 0) print *, "read parameters from file", betainj 

	    close (20)	
	    call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
	    
	    if(rank .eq. 0) then !save into restart directory
	       open(unit=10,file="./restart/adjust_param.txt",action="write",status="replace")
	       write (10,*) "betainj=",betainj
               write (10,*) "left_wall_speedup=",left_wall_speedup
	       close (10)
	    endif
	 endif
	

end subroutine inject_particles_user



!-------------------------------------------------------------------------------
! 				subroutine field_bc_user()	
!										
! Applies boundary conditions specific to user problem. 
! 
!-------------------------------------------------------------------------------

	subroutine field_bc_user()
	implicit none
        
        real(dprec) xmin, xmax
        integer i1, i2
 ! impose conductor behind the left reflecting wall
        
        xmin=1.
        xmax=xinject-5. !global coordinates
        i1=iloc(int(xmin)) !convert to local index
        i2=iloc(int(xmax))

        if(i1 .ne. i2) then  !this means we are on the right CPU
           ey(i1:i2,:,:)=0.
           ez(i1:i2,:,:)=0.
        endif

        xmin=mx0-2.
        xmax=mx0
        i1=iloc(int(xmin)) !convert to local index
        i2=iloc(int(xmax))

        if(i1 .ne. i2) then 
           bx(i1:i2,:,:)=binit*cos(btheta) 
           by(i1:i2,:,:)=binit*sin(btheta)*sin(bphi)
           bz(i1:i2,:,:)=binit*sin(btheta)*cos(bphi)
           ex(i1:i2,:,:)=0.				
           ey(i1:i2,:,:)=(-beta)*bz(i1:i2,:,:) 
           ez(i1:i2,:,:)=-(-beta)*by(i1:i2,:,:)    
        endif


        if(rightclean .eq. 1 .and. (lap .eq. rightclean_start+1 .or. &
             (lap .gt. rightclean_start+1 .and. modulo(lap, rightclean_interval) .eq. 1))) then

           xmin=xinject2+1 !xinject2+1 !+100 !+100 !hack
           xmax=mx0-3
           i1=iloc(int(xmin))
           i2=iloc(int(xmax))
           if(i1 .ne. i2) then
              bx(i1:i2,:,:)=binit*cos(btheta) 
              by(i1:i2,:,:)=binit*sin(btheta)*sin(bphi)
              bz(i1:i2,:,:)=binit*sin(btheta)*cos(bphi)
              ex(i1:i2,:,:)=0.                         
              ey(i1:i2,:,:)=(-beta)*bz(i1:i2,:,:) 
              ez(i1:i2,:,:)=-(-beta)*by(i1:i2,:,:)    
           endif
        endif !if rightclean

	end subroutine field_bc_user

!------------------------------------------------------------------------------

	subroutine particle_bc_user()
	implicit none
	real invgam, gammawall, betawall, gamma,q0
        real x0,y0,z0, tfrac, xcolis, ycolis, zcolis
	integer n1,i0,i1, iter, i00
        logical in, rw
        real (dprec) :: wallocrgt, walloc0

	!loop over particles to check if they crossed special boundaries, like reflecting walls
	!outflow and periodic conditions are handled automatically in deposit_particles
	!
	!This routine is called after the mover and before deposit, thus allowing to avoid
        ! charge deposition behind walls. 
	
        if(wall) then 
           gammawall=wallgam
           betawall=sqrt(1.-1/gammawall**2)
           walloc=walloc+betawall*c

           if(lap .gt. 20000 .and. modulo(lap, 200) .eq. 0 .and. betawall .ne. 0.) then 
              walloc=walloc +200*((left_wall_speedup-1)*betawall*c) ! speedup of left wall
           endif

	   do iter=1,2
              if(iter.eq.1) then 
                  i0=1
                  i1=ions
                  q0=qi
                  i00=0
	      else
		 i0=maxhlf+1
		 i1=maxhlf+lecs
                 q0=qe
                 i00=maxhlf
	      endif
	      
	      n1=i0
	     do while(n1 <= i1) 
	    
	   if(p(n1)%x+mxcum .lt. walloc) then 
	      gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	      
	      !this algorithm ignores change in y and z coordinates
	      !during the scattering. Including it can result in rare
	      !conditions where a particle gets stuck in the ghost zones. 
	      !This can be improved. 

	      !unwind x location of particle
	      x0=p(n1)%x-p(n1)%u/gamma*c
	      y0=p(n1)%y !-p(n1)%v/gamma*c
	      z0=p(n1)%z !-p(n1)%w/gamma*c

	      !unwind wall location
	      walloc0=walloc-betawall*c-mxcum
	      
	      !where did they meet?
	      tfrac=abs((x0-walloc0)/(betawall*c-p(n1)%u/gamma*c))

	      if(tfrac .lt. 1) then 
		 xcolis=x0+p(n1)%u/gamma*c*tfrac
		 ycolis=y0	!+p(n1)%v/gamma*c*tfrac
		 zcolis=z0	!+p(n1)%w/gamma*c*tfrac

	      !deposit current upto intersection
		 q=p(n1)%ch*q0 !real(splitratio)**(1.-real(p(n1)%splitlev))*q0
		 call zigzag(xcolis,ycolis,zcolis,x0,y0,z0,in)

	      !reset particle momentum, getting a kick from the wall
		 p(n1)%u=gammawall**2*gamma*(2*betawall - p(n1)%u/gamma*(1 + betawall**2))
		 gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	      
		 tfrac=min(abs((p(n1)%x-xcolis)/max(abs(p(n1)%x-x0),1e-9)),1.)
              !move particle from the wall position with the new velocity
 
		 p(n1)%x = xcolis + p(n1)%u/gamma*c * tfrac
		 p(n1)%y = ycolis !+ p(n1)%v/gamma*c * tfrac
		 p(n1)%z = zcolis !+ p(n1)%w/gamma*c * tfrac
	     
!now clean up the piece of trajectory behind the wall, that deposit_particles will be adding when it 
!unwinds the position of the particle by the full timestep. 
	      
		 q=-q
		 call zigzag(xcolis,ycolis,zcolis,p(n1)%x-p(n1)%u/gamma*c, & 
		 p(n1)%y-p(n1)%v/gamma*c,p(n1)%z-p(n1)%w/gamma*c,in)
		      
	      else		!if tfrac < 1? 
!       remove the particle
		 call copyprt(p(i1),p(n1))
		 n1=n1-1
		 i1=i1-1
	      endif		!if tfrac < 1
	   endif
	     n1=n1+1
	enddo			! while n1 (particles)

	if(iter .eq. 1)  ions=i1-i00 
	if(iter .eq. 2)  lecs=i1-i00
	enddo			! Do iter  (species)

	endif	! if(wall)

 !=============================================================
!add a second wall at the injector
!	
! algorithm: right wall is moving with background flow starting from xinject2 on every step
! in the process it reflects whatever crossed it, imparting the momentum from reflection 
! off moving wall. 
 rw=.true.


 if(rightclean .eq. 1 .and. lap .gt. rightclean_start &
      .and. modulo(lap, rightclean_interval) .eq. 1 ) rw=.false.

 if(wall .and. rightwall .gt. 0 .and. rw) then 

!	   gammawall=-gamma0 !1.
	   betawall= -beta !-sqrt(1.-1/gammawall**2) !0.  !minus sign to denote left-moving nature of the wall
	    gammawall = 1./sqrt(1.-betawall**2)

	   wallocrgt=xinject2 - beta*c!-betawall*c  !+ betainj*c !+betainj*c

	   do iter=1,2
	      if(iter.eq.1) then 
		 i0=1
		 i1=ions
		 q0=qi
	      else
		 i0=maxhlf+1
		 i1=maxhlf+lecs
		 q0=qe
	      endif
		 
	   do n1=i0,i1	

	   if(p(n1)%x+mxcum .gt. wallocrgt .and. p(n1)%x+mxcum .lt. x2in) then 

	      gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	     
	      !this algorithm ignores change in y and z coordinates
	      !during the scattering. Including it can result in rare
	      !conditions where a particle gets stuck in the ghost zones. 
	      !This can be improved. 

	      !unwind y location of particle
	      x0=p(n1)%x-p(n1)%u/gamma*c
	      y0=p(n1)%y !-p(n1)%v/gamma*c
	      z0=p(n1)%z !-p(n1)%w/gamma*c

	      !unwind wall location
	      walloc0=wallocrgt-betawall*c-mxcum
	      
	      !where did they meet?
	      tfrac=abs((x0-walloc0)/max( abs(betawall*c-p(n1)%u/gamma*c),1e-9))

!	      tfrac=min(abs((x0-walloc0)/max(abs(betawall*c-p(n1)%u/gamma*c),1e-9)),1.)
	if(tfrac .le. 2) then 
	      xcolis=x0+p(n1)%u/gamma*c*tfrac !+mycum
	      ycolis=y0 !+p(n1)%v/gamma*c*tfrac
	      zcolis=z0 !+p(n1)%w/gamma*c*tfrac

	      !deposit current upto intersection
	      q=p(n1)%ch*q0 !real(splitratio)**(1.-real(p(n1)%splitlev))*q0

	      call zigzag(xcolis,ycolis,zcolis,x0,y0,z0,in)

	      !reset particle momentum, getting a kick from the wall
	      p(n1)%u=gammawall**2*gamma*(2*betawall - p(n1)%u/gamma*(1 + betawall**2))
	      gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	      
	      tfrac=min(abs((p(n1)%x-xcolis)/max(abs(p(n1)%x-x0),1e-9)),1.)
              !move particle from the wall position with the new velocity
 
	      p(n1)%x = xcolis + p(n1)%u/gamma*c * tfrac
	      p(n1)%y = ycolis !+ p(n1)%v/gamma*c * tfrac
	      p(n1)%z = zcolis !+ p(n1)%w/gamma*c * tfrac
	     
!now clean up the piece of trajectory behind the wall, that deposit_particles will be adding when it 
!unwinds the position of the particle by the full timestep. 
	      
	      q=-q
	      call zigzag(xcolis,ycolis,zcolis,p(n1)%x-p(n1)%u/gamma*c, & 
	              p(n1)%y-p(n1)%v/gamma*c,p(n1)%z-p(n1)%w/gamma*c,in)
	      
	   else !if tfrac < 2
	      print *, "Warning: particles far beyond right wall", rank, n1, p(n1)%x+mxcum, wallocrgt 
	   endif

	     endif
	   
	enddo	! Do n1 (particles)
	enddo	! Do iter  (species)
	endif !if wall
	       		
        		

	end subroutine particle_bc_user

!-------------------------------------------------------------------------------
! 				subroutine shift_domain_user
!										
! shift fields and particles backward
! only used for some shock problems. Ignore but don't delete. 
!-------------------------------------------------------------------------------
 subroutine shift_domain_user
   implicit none

	integer :: i, n1, n, i1, i2 ! dummy 
	
	integer :: iL, iR	! the first and last ranks in a row
	
	logical :: in
	
	integer :: error , ierr
	
	real, allocatable, dimension(:,:,:) :: tempbx, tempby, tempbz, tempex, tempey,tempez

	integer mxold, mxave, shiftlength
	
	integer buffsize1
 
    shiftlength=walloc-leftwall	! length shifted backward
    if(mxl(1)-shiftlength .lt. leftwall+3 .or. shiftlength .eq. 0) return

	! This condition can be changed.
 if(walloc .gt. 4*leftwall .or. walloc .gt. mxl(1)-10. ) then
    
    if(rank.eq.0)print*, 'domain shift starts', walloc, 4*leftwall, mxl(1)-10.
   
    
    if(modulo(rank,sizex) .eq. 0) then
       
       if(rank.eq.0)print*, 'old walloc=',walloc, 'new walloc=',walloc-shiftlength,"shiftlength", shiftlength
					
       mxold=mx
       mx=mxold-shiftlength
       ix=1
       iy=mx
       iz=iy*my
       lot=mx*my
       
#ifdef twoD
       iz=0
       lot=mx*my
#endif	
       
       i1=1
       i2=i1+shiftlength
       
       if(rank.eq.0)print*, 'mxold, mx=',mxold, mx
     		
       allocate(tempbx(mx,my,mz),tempby(mx,my,mz),tempbz(mx,my,mz))
       allocate(tempex(mx,my,mz),tempey(mx,my,mz),tempez(mx,my,mz))
       
       tempbx=0. !binit*cos(btheta)
       tempby=0. !binit*sin(btheta)*sin(bphi)
       tempbz=0. !binit*sin(btheta)*cos(bphi)				
       tempex=0. !-(-beta)*tempbz				
       tempey=0.
       tempez=0. !-(beta)*tempbx
       
       tempbx(i1:mx,:,:)=bx(i2:mxold,:,:)
       tempby(i1:mx,:,:)=by(i2:mxold,:,:)
       tempbz(i1:mx,:,:)=bz(i2:mxold,:,:)
       tempex(i1:mx,:,:)=ex(i2:mxold,:,:)
       tempey(i1:mx,:,:)=ey(i2:mxold,:,:)
       tempez(i1:mx,:,:)=ez(i2:mxold,:,:)
   
       i1=int(walloc-shiftlength)
       tempex(1:i1-5,:,:)=0.
       tempey(1:i1-5,:,:)=0.
       tempez(1:i1-5,:,:)=0.
       
			
       p(1:ions)%x=p(1:ions)%x-1.*shiftlength
       p(maxhlf+1:maxhlf+lecs)%x=p(maxhlf+1:maxhlf+lecs)%x-1.*shiftlength

		
       deallocate(bx,by,bz,ex,ey,ez,curx,cury,curz)
       deallocate(bufferin1,bufferin2)
       deallocate(bufferin1y,bufferin2y)
       deallocate(sendbufy,sendbufz)
       deallocate(temp)
#ifdef filter2
       deallocate(yghost, zghost) !for filter2
#endif
       deallocate(poutminus, poutplus,pinminus,pinplus)
       deallocate(poutlft, poutrgt,pinlft,pinrgt)
       deallocate(pall)
			
		
       buffsize1 = max(3*int(ppc0*c*max((mx-5),(my-5))),60000)
       
       allocate(bx(mx,my,mz),by(mx,my,mz),bz(mx,my,mz))
       allocate(ex(mx,my,mz),ey(mx,my,mz),ez(mx,my,mz))
       allocate(curx(mx,my,mz),cury(mx,my,mz),curz(mx,my,mz))
       allocate(bufferin1(mx,my,2),bufferin2(mx,my,2))		
       allocate(bufferin1y(mx,2,mz),bufferin2y(mx,2,mz))
       allocate(sendbufy(mx,2,mz),sendbufz(mx,my,2))
       allocate(temp(mx,my,mz))
#ifdef filter2
       allocate(yghost(mx,2*ntimes,mz))
       allocate(zghost(mx,my,2*ntimes)) !for filter2	
#endif
       allocate(poutminus(buffsize1),poutplus(buffsize1))		
       allocate(pinminus(buffsize1),pinplus(buffsize1))
       allocate(poutrgt(buffsize1),poutlft(buffsize1))		
       allocate(pinrgt(buffsize1),pinlft(buffsize1))
       allocate(pall(lot))
       
       pall=0
       curx=0.
       cury=0.
       curz=0.

       bufferin1=0.
       bufferin2=0.
       bufferin1y=0.
       bufferin2y=0.
       sendbufy=0.
       sendbufz=0.		
       temp=0.
       
       bx=tempbx
       by=tempby
       bz=tempbz
       ex=tempex
       ey=tempey
       ez=tempez
       
       deallocate(tempbx,tempby,tempbz,tempex,tempey,tempez)
	
    endif !if(modulo(rank,sizex) .eq. 0)
     	
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
     	
    ! shift global variables
    walloc=walloc-1.*shiftlength
    xinject2=xinject2-1.*shiftlength
    xinject3=xinject3-1.*shiftlength	
    x2in=x2in-1.*shiftlength

    call mpi_allgather(mx,1,mpi_integer,mxl,1,mpi_integer &
         ,mpi_comm_world, error)
			
    iL=1
    iR=modulo(rank,sizex)+1
     	
    mxcum=sum(mxl(iL:iR)-5)-(mxl(iR)-5)
    mx0=sum(mxl(1:sizex)-5)+5		
		
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(rank .eq. 0) print *,"domain shift done  "
 endif

 end subroutine shift_domain_user
 
#ifdef twoD
end module m_user
#else
end module m_user_3d
#endif


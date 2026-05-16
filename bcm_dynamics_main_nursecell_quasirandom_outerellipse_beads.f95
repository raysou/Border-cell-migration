	!!BCM Dynamics with noise and self-propulsion (semi-random arrangement of nurse cells)

	module radii
	implicit none
	double precision,parameter:: a = 24.0d0  !! Major axis length
	double precision,parameter:: b = 12.0d0  !! Minor axis length
	double precision,parameter:: rn = 0.15d0*a/2.0d0 !! Nurse cell radius
	double precision,parameter:: rb = 0.03d0*a/2.0d0 !! Border cell radius
	double precision,parameter:: rb_cluster = 0.1d0*a/2.0d0 !! Border cell cluster radius
	double precision,parameter:: rp = rb/2.0d0    !! Polar cell radius 
	double precision,parameter:: pi = 4.0d0*datan(1.0d0)
	double precision,parameter:: lo_n = 0.1d0 !! Natural spring length of nurse cells
	double precision,parameter:: lo_e = 0.1d0 !! Natural spring length of nurse cells
	double precision,parameter:: lo_b = 0.048d0 !! Natural spring length of border cells 0.045
	double precision,parameter:: lo_p = 0.046d0 !! Natural spring length of polar cells 0.023
	double precision,parameter:: ovrlp_thrsld = 0.08715d0 !! interpenetration overlap threshold. This is not applied
	integer,parameter:: nrexcl = 2 !! for intraring bead-bead repulsion. This is not applied
	integer:: ellipse_bead_num
	end module radii

	module arrays
	use radii
	implicit none
	!integer,parameter:: mb=6,nb=nint(2*pi*rb/lo_b),sb=nb+1,mp=2,np=nint(2*pi*rp/lo_p),sp=np+1,&
			    !&mn=6,nn=nint(2*pi*rn/lo_n),sn=nn+1
	integer,parameter:: mb=6,nb=40,sb=nb+1,mp=2,np=15,sp=np+1,&
			    &mn=6,nn=100,sn=nn+1, full_ellipse_num=1500
	double precision:: xb(mb,0:sb),yb(mb,0:sb),xp(mp,0:sp),yp(mp,0:sp),xn(mn,0:sn),yn(mn,0:sn),&
	       &xe(full_ellipse_num),ye(full_ellipse_num),theta_n(mn),theta_b(mb),theta_p(mp),&
	       &cmx_n(mn),cmx_b(mb),cmx_p(mp),cmy_n(mn),cmy_b(mb),cmy_p(mp),&
	       &xb_prev(mb,nb), yb_prev(mb,nb),xp_prev(mp,np), yp_prev(mp,np),xn_prev(mn,nn), yn_prev(mn,nn),&
	       &xe_prev(full_ellipse_num),ye_prev(full_ellipse_num)
	double precision:: migforce_theta(mb) !! related to mig. force direction
	double precision:: velx_b(mb,nb),vely_b(mb,nb),velx_p(mp,np),vely_p(mp,np),velx_n(mn,nn),vely_n(mn,nn),&
			   &velx_e(full_ellipse_num),vely_e(full_ellipse_num)
	end module arrays

	module forces         
	use arrays
	implicit none
	double precision:: fxn(mn,nn),fyn(mn,nn),f_intx_n(mn,nn),f_inty_n(mn,nn),&
                          &fxb(mb,nb),fyb(mb,nb),f_intx_b(mb,nb),f_inty_b(mb,nb),&
			  &fxp(mp,np),fyp(mp,np),f_intx_p(mp,np),f_inty_p(mp,np),&
			  &f_intx_e(full_ellipse_num),f_inty_e(full_ellipse_num),&
			  &normal_force_unitvec_bx(mb,nb),normal_force_unitvec_by(mb,nb)
	double precision:: fxe(full_ellipse_num), fye(full_ellipse_num)
	double precision:: f_fricx_b(mb,nb),f_fricy_b(mb,nb),f_fricx_p(mp,np),f_fricy_p(mp,np),& !! friction forces on border & polar
			   &f_fricx_n(mn,nn),f_fricy_n(mn,nn),& !! friction forces on nurse
			   &f_fricx_e(full_ellipse_num),f_fricy_e(full_ellipse_num) !! friction force on outer ellipse beads
	double precision:: migration_forcex(mb,nb), migration_forcey(mb,nb) !! test purpose
	logical:: this_border_gets_mig_force(mb,nb)
	end module forces

	module parameters
	use radii
	implicit none
	double precision,parameter:: kn=170.0d0,pn=22.0d0 ! Spring Const and Pressure for Nurse cells
	double precision,parameter:: kb=200.0d0,pb=27.0d0 ! Spring Const and Pressure for Border cells
	double precision,parameter:: kp=150.0d0,pp=22.0d0 ! Spring Const and Pressure for Polar cells
	double precision,parameter:: ke=700.0d0		  ! Spring Const for boundary ellipse beads
	double precision,parameter:: rc_rep_nn=0.18d0,rc_adh_nn=0.28d0,k_adh_nn=100.0d0,k_rep_nn=500.0d0 !Nurse-Nurse IntParameters
	double precision,parameter:: rc_rep_nb=0.18d0,rc_adh_nb=0.28d0,k_adh_nb=0.1d0,k_rep_nb=500.0d0!Nurse-Border Int Parameters
	double precision,parameter:: rc_rep_bb=2*lo_b,rc_adh_bb=0.15d0,k_adh_bb=70.0d0,k_rep_bb=200.0d0 !Border-Border Int Parameters
	double precision,parameter:: rc_rep_bp=2*lo_b,rc_adh_bp=0.15d0,k_adh_bp=70.0d0,k_rep_bp=200.0d0 !Border-Polar Int Parameters
	double precision,parameter:: rc_rep_pp=2*lo_p,rc_adh_pp=0.15d0,k_adh_pp=65.0d0,k_rep_pp=200.0d0 !polar-Polar Int Parameters
	double precision,parameter:: rc_rep_ne=lo_n,rc_adh_ne=0.15d0,k_adh_ne=40d0,k_rep_ne=200.0d0 !ellipse-nurse Int Parameters
	double precision,parameter:: rc_rep_be=lo_e,rc_adh_be=0.15d0,k_adh_be=8.0d0,k_rep_be=200.0d0 !ellipse-border Int Parameters
	double precision,parameter:: rc_rep_pe=lo_e,rc_adh_pe=0.15d0,k_adh_pe=0.0d0,k_rep_pe=200.0d0 !ellipse-polar Int Parameters
	double precision,parameter:: Vo_n = 0.0d0 ! self propulsion speed for nurse cells
	double precision,parameter:: Vo_b = 0.0d0 ! self propulsion speed for border cells
	double precision,parameter:: Vo_p = 0.0d0 ! self propulsion speed for polar cells
	double precision,parameter:: mean_n = 0.0d0 ! Noise-mean for nurse cells
	double precision,parameter:: mean_b = 0.0d0 ! Noise-mean for border cells
	double precision,parameter:: mean_p = 0.0d0 ! Noise-mean for polar cells
	double precision,parameter:: var_n = 0.5d0 ! Noise-strength for nurse cells
	double precision,parameter:: var_b = 0.5d0 ! Noise-strength for border cells
	double precision,parameter:: var_p = 0.5d0 ! Noise-strength for polar cells
	double precision,parameter:: wn = 0.00d0   ! Constant factor for direct noise in nurse cells
	double precision,parameter:: wb = 0.00d0   ! Constant factor for direct noise in border cells
	double precision,parameter:: wp = 0.00d0   ! Constant factor for direct noise in polar cells
	double precision,parameter:: mignoise_mean_b=0.0d0, mignoise_std_b=0.004d0 !! related to migration dir. noise
	end module parameters

   !!! Main program starts here !!!
	program main

	use radii
	use arrays
	use forces
	
	implicit none
	integer:: l,i,j1,jf
	double precision:: dt,ti,tf
	Character(len=1000):: path

	path = 'all_cells/std0.02/Kn170/rx0.05/test_runs_Kn170_rx0.05/With_outer_ellipse_beads/&
		&High_kadhnn/kadhnn100_kadhne40/' !! This is the directory path where you wish to dump the data
	!filepath = 'cell_to_cell/initangl0/Kn170/test_random/'
	!open(16,file=filepath//'test2_2polarbeads15_Kb200_Kadhbp70_Kn170_rx0.1_cmx0.5_yforce_ry0.08_cmx0.5coy0.1.dat',&
         !          & status='unknown')
	open(16,file=trim(adjustl(path))//'Kadhnn100_kadhne40_rx0.1_cmx0.5_cox0.25_coy0.11_ry0.1_noisestd_0.004_20.dat',&
                   & status='unknown')
	open(18,file=trim(adjustl(path))//'ellipse_Kadhnn100_kadhne40_rx0.1_cmx0.5_cox0.25_coy0.11_ry0.1_noisestd_0.004_20.dat',&
                   & status='unknown')	
	open(20,file=trim(adjustl(path))//'com_vel_Kadhnn100_kadhne40_rx0.1_cmx0.5_cox0.25_coy0.11_ry0.1_noisestd_0.004_20.dat',&
		& status='unknown')
         !open(17,file=filepath//'BorderCells_different_int_forces_P2.5_quasirandom_arrngmnt.dat',&
         !& status='unknown',position='append')
         !open(18,file=filepath//'PloarCells_different_int_forces_P2.5_quasirandom_arrngmnt.dat',&
         !& status='unknown',position='append')

        call cpu_time(ti)
	
	dt=0.001d0 !! Integration time-step
	jf=1000000 !! Total number of iterations 

	call initial   
	call initial_angle

	do j1=1,jf  !! Iteration loop begins here
        
	call force(j1)
	call interaction(j1)
	if(j1.gt.10000) then  
	call move_noise(j1,dt)
	else
	call move_deterministic(dt)
	end if

	if(j1.gt.10000) then  
        if(mod(j1,2000).eq.0) then  !! Data will be stored every 2000 steps
          do l=1,mn 
            do i=1,nn        
            	    write(16,*)j1,l,i,xn(l,i),yn(l,i),0d0,0d0	!! Data writing for nurse cells              
            end do
          end do
          do l=1,mb 
            do i=1,nb        
            	    write(16,*)j1,mn+l,i,xb(l,i),yb(l,i),migration_forcex(l,i),migration_forcey(l,i)  !! Data writing for border cells            
            end do
          end do
          do l=1,mp 
            do i=1,np        
            	    write(16,*)j1,mn+mb+l,i,xp(l,i),yp(l,i),0d0,0d0	    !! Data writing for polar cells     
            end do
          end do
	  		
	  do i=1,ellipse_bead_num
		write(18,*)j1, i, xe(i), ye(i)  !! Fixed outer elliptic boundary
	  end do

        end if
	end if

	end do

	call cpu_time(tf)

	write(*,*)'time=',tf-ti

end program main
!!! Main program ends here !!!


subroutine initial  !! Subroutine for creating initial tissue configuration !!

use parameters
use radii
use arrays

implicit none
double precision:: theta,dx,dy,x_tmp,y_tmp,xc_n(10),yc_n(10),xc_b,yc_b,xc_p,yc_p,theta_bcell_centre,g
integer:: l,k,i,t,d,p,s, idum
double precision:: ran2

CALL SYSTEM_CLOCK(COUNT=idum)
   !! Nurse cells arrangement !!
	s=0
	t=0
	do l=1,mn/2
5		xc_n(l) = (2*rb_cluster + rn) + (0.75d0*a - rn - (2*rb_cluster + rn))*ran2(idum)
		yc_n(l) = rn + rc_rep_nn + (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn-0.0d0 - (rn+rc_rep_nn))*ran2(idum)

		if(l.eq.1) then
		xc_n(l) = (2*rb_cluster + rn) + rc_rep_nb !+ 1.5d0!+ 0.25d0*rand(0) !+ rc_rep_nb !+ 1.5d0
		!yc_n(l) = rn + rc_rep_nn +0.2d0!+ (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
		yc_n(l) = rn + rc_rep_nn + (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn-0.25d0 - (rn+rc_rep_nn))*ran2(idum)
		!!!!### here, above, rn-0.25....the extra -0.25 is used to aviod 'exceeded ellipse boundary case'####!!!!
		else if (l.eq.2) then
		!xc_n(l) = xc_n(1) + 2*rn + rc_rep_nn
		!yc_n(l) = rn + rc_rep_nn !+ (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
!		yc_n(l) = rn + rc_rep_nn + (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
		else 
		!xc_n(l) = 0.75d0*a - rn - rc_rep_nn
		!yc_n(l) = rn + rc_rep_nn !+ (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
!		yc_n(l) = rn + rc_rep_nn + (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
		end if
		
		if(l.eq.1) then
		if((yc_n(l).lt.(rn+0.1d0)).or.&
		&((yc_n(l)+rn+0.25d0).gt.(b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2)))) goto 5		
		end if

		if(l.gt.1) then
		do k=1,l-1
		g = sqrt((yc_n(l)-yc_n(k))**2 + (xc_n(l)-xc_n(k))**2)
		if((g.lt.(2*rn + 0.1d0)).or.(yc_n(l).lt.(rn+0.1d0)).or.&
		&((yc_n(l)+rn+0.0d0).gt.(b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2)))) goto 5
		end do
		end if
		!write(*,*)l,xc_n(l),yc_n(l)
	end do

 	t=0
	do l=mn/2 + 1,mn
6		xc_n(l) = (2*rb_cluster + rn) + (0.75d0*a - rn - (2*rb_cluster + rn))*ran2(idum)
		yc_n(l) = -(rn+rc_rep_nn) - (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn-0.0d0 - (rn+rc_rep_nn))*ran2(idum)

		if(l.eq.(mn/2 + 1)) then
		xc_n(l) = (2*rb_cluster + rn) + rc_rep_nb !+ 1.5d0!+ 0.25d0*rand(0) !+ rc_rep_nb !+ 1.5d0
		!yc_n(l) = -(rn + rc_rep_nn) -0.2d0!- (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
		yc_n(l) = -(rn + rc_rep_nn) - (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn-0.25d0 - (rn+rc_rep_nn))*ran2(idum)
		else if (l.eq.(mn/2 + 2)) then
		!xc_n(l) = xc_n(mn/2 + 1) + 2*rn + rc_rep_nn
		!yc_n(l) = -(rn + rc_rep_nn) !- (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
!		yc_n(l) = -(rn + rc_rep_nn) - (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
		else 
		!xc_n(l) = 0.75d0*a - rn - rc_rep_nn
		!yc_n(l) = -(rn + rc_rep_nn) !- (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
!		yc_n(l) = -(rn + rc_rep_nn) - (b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2) - rn - (rn+rc_rep_nn))*rand(0)
		end if

		if(l.eq.(mn/2 + 1)) then
		if((yc_n(l).gt.(-rn-0.1d0)).or.&
		&((yc_n(l)-rn-0.25d0).lt.(-b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2)))) goto 6		
		end if

		if(l.gt.(mn/2 + 1)) then
		do k=mn/2 + 1,l-1
		g = sqrt((yc_n(l)-yc_n(k))**2 + (xc_n(l)-xc_n(k))**2) 
		if((g.lt.(2*rn + 0.1d0)).or.(yc_n(l).gt.(-rn-0.1d0)).or.&
		&((yc_n(l)-rn-0.0d0).lt.(-b*sqrt(0.25d0 - ((xc_n(l)-a/2.0d0)/a)**2)))) goto 6
		end do
		end if
		!write(*,*)l,xc_n(l),yc_n(l)
	end do

	do l=1,mn
		theta = 0.0d0
	   do i=1,nn
		   xn(l,i) = rn*cos(theta) + 0.001d0*(2.0d0*ran2(idum)-1.0d0) + xc_n(l)
		   yn(l,i) = rn*sin(theta) + 0.001d0*(2.0d0*ran2(idum)-1.0d0) + yc_n(l)
		   theta = theta + 2.0d0*pi/nn
		                                                
		  ! write(12,*)l,i,xn(l,i),yn(l,i)
	   end do
	end do



  !! Border cells arrangement !!
	theta_bcell_centre = 2*pi/6
	do l=1,mb
	           xc_b = (2*rp + rb + 0.1d0)*cos(theta_bcell_centre) + rb_cluster
	  	   yc_b = (2*rp + rb + 0.1d0)*sin(theta_bcell_centre)
		   theta = 0.0d0
		   do i=1,nb
			   xb(l,i) = rb*cos(theta) + 0.001d0*(2.0d0*ran2(idum)-1.0d0) + xc_b
			   yb(l,i) = rb*sin(theta) + 0.001d0*(2.0d0*ran2(idum)-1.0d0) + yc_b
		           theta = theta + 2.0d0*pi/nb 

			  ! write(13,*)l,i,xb(l,i),yb(l,i)                  
		   end do
		   theta_bcell_centre = theta_bcell_centre + 2*pi/6.0d0
	end do

  !! Polar cells arrangement !!

	do l=1,mp
		yc_p = 0.0d0
		xc_p = rb_cluster + rp + 0.05d0
		if(l.eq.1) xc_p = rb_cluster - rp
		theta = 0.0d0
	        do i=1,np
			 xp(l,i) = rp*cos(theta) + 0.001d0*(2.0d0*ran2(idum)-1.0d0) + xc_p
			 yp(l,i) = rp*sin(theta) + 0.001d0*(2.0d0*ran2(idum)-1.0d0) + yc_p
                         theta = theta + 2.0d0*pi/np  

			! write(14,*)l,i,xp(l,i),yp(l,i)                 
		   end do
	end do

  !! The beads for the boundary ellipse up to the oocyte boundary !!
        theta = 0.0d0
	    ellipse_bead_num = 0
	    do i=1,full_ellipse_num
	       x_tmp = (a/2.0d0)*(1+cos(theta)) 
	       y_tmp = (b/2.0d0)*sin(theta)	
	       theta = theta + 2.0d0*pi/full_ellipse_num    
	       if (x_tmp.le.(0.75d0*a)) then
	          ellipse_bead_num = ellipse_bead_num + 1
		  	  xe(ellipse_bead_num) = x_tmp
		      ye(ellipse_bead_num) = y_tmp 
           end if  
	       !write(15,*)i,xe(i),ye(i)
	    end do
 
          
return
end 



 	!!*** Subroutine for the forces of interaction ***!!
	subroutine interaction(j1)

	use radii        
	use arrays
	use forces
    use parameters
    implicit none
    integer:: i,j,j1,l,q
    double precision:: r,factor,frepx,frepy,dx,dy,fadhx,fadhy,rlist,ti,tf
	double precision:: cm_d,cm_dx,cm_dy,overlap
	integer:: icell,jcell,jcell0,nabor 
	logical:: same_ring_beyond_nrexcl
	double precision:: relative_velx_bp,relative_vely_bp,angl_diff_wrt_clk,angl_diff_wrt_antclk, &
			&rel_vel_along_clk_vec,rel_vel_along_antclk_vec,f_fricx,f_fricy,fric_coeff_bp, &
			&relative_velx_bn,relative_vely_bn,fric_coeff_bn, &
			&fric_coeff_eb,fric_coeff_en,relative_velx_en,relative_vely_en, &
			&relative_velx_eb,relative_vely_eb, &
			&fric_coeff_ep,relative_velx_ep,relative_vely_ep

	!!#### friction force related parameter ####!!
	fric_coeff_bp = 0.00d0
	fric_coeff_bn = 0.06d0
	fric_coeff_en = 1.0d0
	fric_coeff_eb = 0.00d0
	fric_coeff_ep = 0.00d0
	!!#### friction force related parameter ####!!

	
    f_intx_b=0.0d0
    f_inty_b=0.0d0
    f_intx_n=0.0d0
    f_inty_n=0.0d0
    f_intx_p=0.0d0
    f_inty_p=0.0d0
	f_intx_e=0.0d0
	f_inty_e=0.0d0

	f_fricx_b = 0d0
	f_fricy_b = 0d0
	f_fricx_p = 0d0
	f_fricy_p = 0d0
	f_fricx_n = 0d0
	f_fricy_n = 0d0
	f_fricx_e = 0d0
	f_fricy_e = 0d0        

	this_border_gets_mig_force(:,:) = .true.     
        
       		 !! Nurse-Nurse interactions!!

					do l=1,mn
					    do i=1,nn
					  	do q=l,mn
						   do j=1,nn
				
						    same_ring_beyond_nrexcl = (l == q) &
						    &.and. (min(abs(i - j), nn - abs(i - j)) > nrexcl)

						    if((l.ne.q).or.(same_ring_beyond_nrexcl)) then
							dx = xn(q,j)-xn(l,i)
							dy = yn(q,j)-yn(l,i)
							r = dsqrt(dx*dx + dy*dy)
						
                      		if(r.lt.rc_rep_nn) then

				          		frepx = -k_rep_nn*(rc_rep_nn-r)*(dx)/r
					  			frepy = -k_rep_nn*(rc_rep_nn-r)*(dy)/r

								cm_dx = cmx_n(q) - cmx_n(l)
								cm_dy = cmy_n(q) - cmy_n(l)
								cm_d = hypot(cm_dx, cm_dy)
				                    ! Overlap metric is the projection of bead-bead vector along centre-centre
				                    !! joining unit vector, scaled by rc_rep to make the metric dimensionless
				                   		overlap = (cm_dx*dx + cm_dy*dy)/(rc_rep_nn*cm_d)
				                    ! i-j joining position vector [dx,dy] makes obtuse angle with the position vector
				                    !! joining the com's [cm_dx,cm_dy] only when i and j penetrate each other's cell.
				                    !! Negativity of the above dot product thus indicates cell-cell overlap.

							      if((l.ne.q).and.(overlap.lt.ovrlp_thrsld)) then 
						                factor = -k_rep_nn*rc_rep_nn/cm_d
						                frepx = frepx + factor*cm_dx
						                frepy = frepy + factor*cm_dy
 	                  				end if

				          			f_intx_n(l,i) = f_intx_n(l,i) + frepx 
				          			f_intx_n(q,j) = f_intx_n(q,j) - frepx

				          			f_inty_n(l,i) = f_inty_n(l,i) + frepy 
				          			f_inty_n(q,j) = f_inty_n(q,j) - frepy

                        		else if((r.le.rc_adh_nn).and.(r.ge.rc_rep_nn)) then
							     if(l.ne.q) then

								fadhx = k_adh_nn*(rc_adh_nn-r)*(dx)/r
								fadhy = k_adh_nn*(rc_adh_nn-r)*(dy)/r

								f_intx_n(l,i) = f_intx_n(l,i) + fadhx
								f_intx_n(q,j) = f_intx_n(q,j) - fadhx

						        f_inty_n(l,i) = f_inty_n(l,i) + fadhy 
								f_inty_n(q,j) = f_inty_n(q,j) - fadhy

							      end if
       							end if
						    end if	
						    end do
						end do
					    end do
					end do

       		 !! Polar-Polar interactions!!

					do l=1,mp
					    do i=1,np
						do q=l,mp
						   do j=1,np
				
						    same_ring_beyond_nrexcl = (l == q) &
						    &.and. (min(abs(i - j), np - abs(i - j)) > nrexcl)

						    if((l.ne.q).or.(same_ring_beyond_nrexcl)) then
							dx = xp(q,j)-xp(l,i)
							dy = yp(q,j)-yp(l,i)
							r = dsqrt(dx*dx + dy*dy)

                      		if(r.lt.rc_rep_pp) then

				          		frepx = -k_rep_pp*(rc_rep_pp-r)*(dx)/r
					  			frepy = -k_rep_pp*(rc_rep_pp-r)*(dy)/r

								cm_dx = cmx_p(q) - cmx_p(l)
								cm_dy = cmy_p(q) - cmy_p(l)
								cm_d = hypot(cm_dx, cm_dy)
				                    ! Overlap metric is the projection of bead-bead vector along centre-centre
				                    !! joining unit vector, scaled by rc_rep to make the metric dimensionless
				                   		overlap = (cm_dx*dx + cm_dy*dy)/(rc_rep_pp*cm_d)
				                    ! i-j joining position vector [dx,dy] makes obtuse angle with the position vector
				                    !! joining the com's [cm_dx,cm_dy] only when i and j penetrate each other's cell.
				                    !! Negativity of the above dot product thus indicates cell-cell overlap.

							      if((l.ne.q).and.(overlap.lt.ovrlp_thrsld)) then 
						                factor = -k_rep_pp*rc_rep_pp/cm_d
						                frepx = frepx + factor*cm_dx
						                frepy = frepy + factor*cm_dy
 	                  			  end if

				          			f_intx_p(l,i) = f_intx_p(l,i) + frepx
				          			f_intx_p(q,j) = f_intx_p(q,j) - frepx

				          			f_inty_p(l,i) = f_inty_p(l,i) + frepy 
				          			f_inty_p(q,j) = f_inty_p(q,j) - frepy


                        		else if((r.le.rc_adh_pp).and.(r.ge.rc_rep_pp)) then
							     if(l.ne.q) then

								fadhx = k_adh_pp*(rc_adh_pp-r)*(dx)/r
								fadhy = k_adh_pp*(rc_adh_pp-r)*(dy)/r

								f_intx_p(l,i) = f_intx_p(l,i) + fadhx
								f_intx_p(q,j) = f_intx_p(q,j) - fadhx

						        f_inty_p(l,i) = f_inty_p(l,i) + fadhy 
								f_inty_p(q,j) = f_inty_p(q,j) - fadhy

							      end if
       							end if
						    end if	
						    end do
						end do
					    end do
					end do
		

       		 !! Border-Border interactions!!

					do l=1,mb
					    do i=1,nb
						do q=l,mb
						   do j=1,nb

						    same_ring_beyond_nrexcl = (l == q) &
						    &.and. (min(abs(i - j), nb - abs(i - j)) > nrexcl)

						    if((l.ne.q).or.(same_ring_beyond_nrexcl)) then
											
							dx = xb(q,j)-xb(l,i)
							dy = yb(q,j)-yb(l,i)
							r = dsqrt(dx*dx + dy*dy)

                      		if(r.lt.rc_rep_bb) then

								this_border_gets_mig_force(q,j) = .false.
								this_border_gets_mig_force(l,i) = .false.

				          		frepx = -k_rep_bb*(rc_rep_bb-r)*(dx)/r
					  			frepy = -k_rep_bb*(rc_rep_bb-r)*(dy)/r

								cm_dx = cmx_b(q) - cmx_b(l)
								cm_dy = cmy_b(q) - cmy_b(l)
								cm_d = hypot(cm_dx, cm_dy)
				                    ! Overlap metric is the projection of bead-bead vector along centre-centre
				                    !! joining unit vector, scaled by rc_rep to make the metric dimensionless
				                   		overlap = (cm_dx*dx + cm_dy*dy)/(rc_rep_bb*cm_d)
				                    ! i-j joining position vector [dx,dy] makes obtuse angle with the position vector
				                    !! joining the com's [cm_dx,cm_dy] only when i and j penetrate each other's cell.
				                    !! Negativity of the above dot product thus indicates cell-cell overlap.

							      if((l.ne.q).and.(overlap.lt.ovrlp_thrsld)) then 
						                factor = -k_rep_bb*rc_rep_bb/cm_d
						                frepx = frepx + factor*cm_dx
						                frepy = frepy + factor*cm_dy
 	                  			  end if

				          			f_intx_b(l,i) = f_intx_b(l,i) + frepx 
				          			f_intx_b(q,j) = f_intx_b(q,j) - frepx

				          			f_inty_b(l,i) = f_inty_b(l,i) + frepy 
				          			f_inty_b(q,j) = f_inty_b(q,j) - frepy


                        		else if((r.le.rc_adh_bb).and.(r.ge.rc_rep_bb)) then

								this_border_gets_mig_force(q,j) = .false.
								this_border_gets_mig_force(l,i) = .false.

							    if(l.ne.q) then
								fadhx = k_adh_bb*(rc_adh_bb-r)*(dx)/r
								fadhy = k_adh_bb*(rc_adh_bb-r)*(dy)/r

								f_intx_b(l,i) = f_intx_b(l,i) + fadhx
								f_intx_b(q,j) = f_intx_b(q,j) - fadhx

						        f_inty_b(l,i) = f_inty_b(l,i) + fadhy 
								f_inty_b(q,j) = f_inty_b(q,j) - fadhy

							      end if
       							end if
						    end if	
						    end do
						end do
					    end do
					end do

       		 !! Border-Polar interactions!!

					do l=1,mb
					    do i=1,nb					    
						do q=1,mp
						   do j=1,np
			
							dx = xp(q,j)-xb(l,i)
							dy = yp(q,j)-yb(l,i)
							r = dsqrt(dx*dx + dy*dy)

                      		   if(r.lt.rc_rep_bp) then

								this_border_gets_mig_force(l,i) = .false.

				          		frepx = -k_rep_bp*(rc_rep_bp-r)*(dx)/r
					  			frepy = -k_rep_bp*(rc_rep_bp-r)*(dy)/r

								cm_dx = cmx_p(q) - cmx_b(l)
								cm_dy = cmy_p(q) - cmy_b(l)
								cm_d = hypot(cm_dx, cm_dy)
				                    ! Overlap metric is the projection of bead-bead vector along centre-centre
				                    !! joining unit vector, scaled by rc_rep to make the metric dimensionless
				                   		overlap = (cm_dx*dx + cm_dy*dy)/(rc_rep_bp*cm_d)
				                    ! i-j joining position vector [dx,dy] makes obtuse angle with the position vector
				                    !! joining the com's [cm_dx,cm_dy] only when i and j penetrate each other's cell.
				                    !! Negativity of the above dot product thus indicates cell-cell overlap.

							      if(overlap.lt.ovrlp_thrsld) then 
						                factor = -k_rep_bp*rc_rep_bp/cm_d
						                frepx = frepx + factor*cm_dx
						                frepy = frepy + factor*cm_dy
 	                  			  end if

				          			f_intx_b(l,i) = f_intx_b(l,i) + frepx 
				          			f_intx_p(q,j) = f_intx_p(q,j) - frepx

				          			f_inty_b(l,i) = f_inty_b(l,i) + frepy 
				          			f_inty_p(q,j) = f_inty_p(q,j) - frepy

						!!!!!********************************************************!!!!!
						!!! Traction-like friction force between border-polar depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_bp = velx_b(l,i)-velx_p(q,j)
								relative_vely_bp = vely_b(l,i)-vely_p(q,j)

								angl_diff_wrt_clk = datan2(-dx,dy) - & !!(dy,-dx) is the clkwise rotated&
								&datan2(relative_vely_bp,relative_velx_bp) !!& int. vector.
								angl_diff_wrt_antclk = datan2(dx,-dy) - & !!(-dy,dx) is the antclkwise&
								&datan2(relative_vely_bp,relative_velx_bp) !!& rotated int. vec.

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_bp*(dy*1.0d0/r)&
								& + relative_vely_bp*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_bp*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_bp*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_b(l,i) = f_fricx_b(l,i) + f_fricx
								f_fricy_b(l,i) = f_fricy_b(l,i) + f_fricy
								f_fricx_p(q,j) = f_fricx_p(q,j) - f_fricx
								f_fricy_p(q,j) = f_fricy_p(q,j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_bp*(-dy*1.0d0/r)&
								& + relative_vely_bp*(dx*1.0d0/r)

								f_fricx = - fric_coeff_bp*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_bp*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_b(l,i) = f_fricx_b(l,i) + f_fricx
								f_fricy_b(l,i) = f_fricy_b(l,i) + f_fricy
								f_fricx_p(q,j) = f_fricx_p(q,j) - f_fricx
								f_fricy_p(q,j) = f_fricy_p(q,j) - f_fricy

								end if
						!!!!!********************************************************!!!!!
                        		else if((r.le.rc_adh_bp).and.(r.ge.rc_rep_bp)) then

								this_border_gets_mig_force(l,i) = .false.

								fadhx = k_adh_bp*(rc_adh_bp-r)*(dx)/r
								fadhy = k_adh_bp*(rc_adh_bp-r)*(dy)/r

								f_intx_b(l,i) = f_intx_b(l,i) + fadhx
								f_intx_p(q,j) = f_intx_p(q,j) - fadhx

						        f_inty_b(l,i) = f_inty_b(l,i) + fadhy 
								f_inty_p(q,j) = f_inty_p(q,j) - fadhy
						!!!!!********************************************************!!!!!
						!!! Traction-like friction force betwn border-polar depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_bp = velx_b(l,i)-velx_p(q,j)
								relative_vely_bp = vely_b(l,i)-vely_p(q,j)

								angl_diff_wrt_clk = datan2(-dx,dy) - &
								&datan2(relative_vely_bp,relative_velx_bp)
								angl_diff_wrt_antclk = datan2(dx,-dy) - &
								&datan2(relative_vely_bp,relative_velx_bp)

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_bp*(dy*1.0d0/r)&
								& + relative_vely_bp*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_bp*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_bp*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_b(l,i) = f_fricx_b(l,i) + f_fricx
								f_fricy_b(l,i) = f_fricy_b(l,i) + f_fricy
								f_fricx_p(q,j) = f_fricx_p(q,j) - f_fricx
								f_fricy_p(q,j) = f_fricy_p(q,j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_bp*(-dy*1.0d0/r)&
								& + relative_vely_bp*(dx*1.0d0/r)

								f_fricx = - fric_coeff_bp*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_bp*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_b(l,i) = f_fricx_b(l,i) + f_fricx
								f_fricy_b(l,i) = f_fricy_b(l,i) + f_fricy
								f_fricx_p(q,j) = f_fricx_p(q,j) - f_fricx
								f_fricy_p(q,j) = f_fricy_p(q,j) - f_fricy

								end if
						!!!!!********************************************************!!!!!
								
       							end if
						      end do
						    end do
						 end do
					       end do

       		 !! Nurse-Border interactions!!
					
					do l=1,mn
					    do i=1,nn					    
						do q=1,mb
						   do j=1,nb
											
							dx = xb(q,j)-xn(l,i)
							dy = yb(q,j)-yn(l,i)
							r = dsqrt(dx*dx + dy*dy)

                      			if(r.lt.rc_rep_nb) then

								this_border_gets_mig_force(q,j) = .true.
								
				          		frepx = -k_rep_nb*(rc_rep_nb-r)*(dx)/r
					  			frepy = -k_rep_nb*(rc_rep_nb-r)*(dy)/r

								cm_dx = cmx_b(q) - cmx_n(l)
								cm_dy = cmy_b(q) - cmy_n(l)
								cm_d = hypot(cm_dx, cm_dy)
				                    ! Overlap metric is the projection of bead-bead vector along centre-centre
				                    !! joining unit vector, scaled by rc_rep to make the metric dimensionless
				                   		overlap = (cm_dx*dx + cm_dy*dy)/(rc_rep_nb*cm_d)
				                    ! i-j joining position vector [dx,dy] makes obtuse angle with the position vector
				                    !! joining the com's [cm_dx,cm_dy] only when i and j penetrate each other's cell.
				                    !! Negativity of the above dot product thus indicates cell-cell overlap.

							      if(overlap.lt.ovrlp_thrsld) then 
						                factor = -k_rep_nb*rc_rep_nb/cm_d
						                frepx = frepx + factor*cm_dx
						                frepy = frepy + factor*cm_dy
 	                  			  end if

				          			f_intx_n(l,i) = f_intx_n(l,i) + frepx 
				          			f_intx_b(q,j) = f_intx_b(q,j) - frepx

				          			f_inty_n(l,i) = f_inty_n(l,i) + frepy 
				          			f_inty_b(q,j) = f_inty_b(q,j) - frepy

						!!!!!********************************************************!!!!!
						!!! Traction-like friction force betwn border-polar depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_bn = velx_n(l,i)-velx_b(q,j)
								relative_vely_bn = vely_n(l,i)-vely_b(q,j)

								angl_diff_wrt_clk = datan2(-dx,dy) - &
								&datan2(relative_vely_bn,relative_velx_bn)
								angl_diff_wrt_antclk = datan2(dx,-dy) - &
								&datan2(relative_vely_bn,relative_velx_bn)

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_bn*(dy*1.0d0/r)&
								& + relative_vely_bn*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_bn*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_bn*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_n(l,i) = f_fricx_n(l,i) + f_fricx
								f_fricy_n(l,i) = f_fricy_n(l,i) + f_fricy
								f_fricx_b(q,j) = f_fricx_b(q,j) - f_fricx
								f_fricy_b(q,j) = f_fricy_b(q,j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_bn*(-dy*1.0d0/r)&
								& + relative_vely_bn*(dx*1.0d0/r)

								f_fricx = - fric_coeff_bn*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_bn*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_n(l,i) = f_fricx_n(l,i) + f_fricx
								f_fricy_n(l,i) = f_fricy_n(l,i) + f_fricy
								f_fricx_b(q,j) = f_fricx_b(q,j) - f_fricx
								f_fricy_b(q,j) = f_fricy_b(q,j) - f_fricy

								end if
						!!!!!********************************************************!!!!!

                        		else if((r.le.rc_adh_nb).and.(r.ge.rc_rep_nb)) then
							
								this_border_gets_mig_force(q,j) = .true.

								fadhx = k_adh_nb*(rc_adh_nb-r)*(dx)/r
								fadhy = k_adh_nb*(rc_adh_nb-r)*(dy)/r

								f_intx_n(l,i) = f_intx_n(l,i) + fadhx
								f_intx_b(q,j) = f_intx_b(q,j) - fadhx

						        f_inty_n(l,i) = f_inty_n(l,i) + fadhy 
								f_inty_b(q,j) = f_inty_b(q,j) - fadhy
						!!!!!********************************************************!!!!!
						!!! Traction-like friction force betwn border-polar depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_bn = velx_n(l,i)-velx_b(q,j)
								relative_vely_bn = vely_n(l,i)-vely_b(q,j)

								angl_diff_wrt_clk = datan2(-dx,dy) - &
								&datan2(relative_vely_bn,relative_velx_bn)
								angl_diff_wrt_antclk = datan2(dx,-dy) - &
								&datan2(relative_vely_bn,relative_velx_bn)

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_bn*(dy*1.0d0/r)&
								& + relative_vely_bn*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_bn*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_bn*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_n(l,i) = f_fricx_n(l,i) + f_fricx
								f_fricy_n(l,i) = f_fricy_n(l,i) + f_fricy
								f_fricx_b(q,j) = f_fricx_b(q,j) - f_fricx
								f_fricy_b(q,j) = f_fricy_b(q,j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_bn*(-dy*1.0d0/r)&
								& + relative_vely_bn*(dx*1.0d0/r)

								f_fricx = - fric_coeff_bn*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_bn*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_n(l,i) = f_fricx_n(l,i) + f_fricx
								f_fricy_n(l,i) = f_fricy_n(l,i) + f_fricy
								f_fricx_b(q,j) = f_fricx_b(q,j) - f_fricx
								f_fricy_b(q,j) = f_fricy_b(q,j) - f_fricy

								end if
						!!!!!********************************************************!!!!!
								
       							end if
						      end do
						    end do
						 end do
					       end do		


       		 !! Nurse - Outer-ellipse interactions!!
					
					do l=1,mn
					    do i=1,nn					    
						 do j=1,ellipse_bead_num
											
							dx = xe(j)-xn(l,i)
							dy = ye(j)-yn(l,i)
							r = dsqrt(dx*dx + dy*dy)

                      			if(r.lt.rc_rep_ne) then
								
				          		frepx = -k_rep_ne*(rc_rep_ne-r)*(dx)/r
					  			frepy = -k_rep_ne*(rc_rep_ne-r)*(dy)/r

				          			f_intx_n(l,i) = f_intx_n(l,i) + frepx 
				          			f_intx_e(j)   = f_intx_e(j) - frepx

				          			f_inty_n(l,i) = f_inty_n(l,i) + frepy 
				          			f_inty_e(j)   = f_inty_e(j) - frepy

						!!!!!********************************************************!!!!!
						!!! Traction-like friction force betwn border-polar depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_en = velx_n(l,i)-velx_e(j)
								relative_vely_en = vely_n(l,i)-vely_e(j)

								angl_diff_wrt_clk = datan2(-dx,dy) - &
								&datan2(relative_vely_en,relative_velx_en)
								angl_diff_wrt_antclk = datan2(dx,-dy) - &
								&datan2(relative_vely_en,relative_velx_en)

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_en*(dy*1.0d0/r)&
								& + relative_vely_en*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_en*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_en*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_n(l,i) = f_fricx_n(l,i) + f_fricx
								f_fricy_n(l,i) = f_fricy_n(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_en*(-dy*1.0d0/r)&
								& + relative_vely_en*(dx*1.0d0/r)

								f_fricx = - fric_coeff_en*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_en*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_n(l,i) = f_fricx_n(l,i) + f_fricx
								f_fricy_n(l,i) = f_fricy_n(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								end if
						!!!!!********************************************************!!!!!

                        		else if((r.le.rc_adh_ne).and.(r.ge.rc_rep_ne)) then

								fadhx = k_adh_ne*(rc_adh_ne-r)*(dx)/r
								fadhy = k_adh_ne*(rc_adh_ne-r)*(dy)/r

								f_intx_n(l,i) = f_intx_n(l,i) + fadhx
								f_intx_e(j)   = f_intx_e(j) - fadhx

						        f_inty_n(l,i) = f_inty_n(l,i) + fadhy 
								f_inty_e(j)   = f_inty_e(j) - fadhy
						!!!!!********************************************************!!!!!
						!!! Traction-like friction force betwn border-polar depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_en = velx_n(l,i)-velx_e(j)
								relative_vely_en = vely_n(l,i)-vely_e(j)

								angl_diff_wrt_clk = datan2(-dx,dy) - &
								&datan2(relative_vely_en,relative_velx_en)
								angl_diff_wrt_antclk = datan2(dx,-dy) - &
								&datan2(relative_vely_en,relative_velx_en)

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_en*(dy*1.0d0/r)&
								& + relative_vely_en*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_en*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_en*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_n(l,i) = f_fricx_n(l,i) + f_fricx
								f_fricy_n(l,i) = f_fricy_n(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_en*(-dy*1.0d0/r)&
								& + relative_vely_en*(dx*1.0d0/r)

								f_fricx = - fric_coeff_en*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_en*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_n(l,i) = f_fricx_n(l,i) + f_fricx
								f_fricy_n(l,i) = f_fricy_n(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								end if
						!!!!!********************************************************!!!!!
								
       							end if
						    end do
						 end do
					    end do		


       		 !! Border - Outer-ellipse interactions!!
					
					do l=1,mb
					    do i=1,nb					    
						 do j=1,ellipse_bead_num
											
							dx = xe(j)-xb(l,i)
							dy = ye(j)-yb(l,i)
							r = dsqrt(dx*dx + dy*dy)

                      		if(r.lt.rc_rep_be) then
								
				          		frepx = -k_rep_be*(rc_rep_be-r)*(dx)/r
					  			frepy = -k_rep_be*(rc_rep_be-r)*(dy)/r

				          			f_intx_b(l,i) = f_intx_b(l,i) + frepx 
				          			f_intx_e(j)   = f_intx_e(j) - frepx

				          			f_inty_b(l,i) = f_inty_b(l,i) + frepy 
				          			f_inty_e(j)   = f_inty_e(j) - frepy

						!!!!!********************************************************!!!!!
						!!! Traction-like friction force betwn border-ellipse depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_eb = velx_b(l,i)-velx_e(j)
								relative_vely_eb = vely_b(l,i)-vely_e(j)

								angl_diff_wrt_clk = datan2(-dx,dy) - &
								&datan2(relative_vely_eb,relative_velx_eb)
								angl_diff_wrt_antclk = datan2(dx,-dy) - &
								&datan2(relative_vely_eb,relative_velx_eb)

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_eb*(dy*1.0d0/r)&
								& + relative_vely_eb*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_eb*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_eb*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_b(l,i) = f_fricx_b(l,i) + f_fricx
								f_fricy_b(l,i) = f_fricy_b(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_eb*(-dy*1.0d0/r)&
								& + relative_vely_eb*(dx*1.0d0/r)

								f_fricx = - fric_coeff_eb*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_eb*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_b(l,i) = f_fricx_b(l,i) + f_fricx
								f_fricy_b(l,i) = f_fricy_b(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								end if
						!!!!!********************************************************!!!!!

                        		else if((r.le.rc_adh_be).and.(r.ge.rc_rep_be)) then
							
								!this_border_gets_mig_force(q,j) = .true.

								fadhx = k_adh_be*(rc_adh_be-r)*(dx)/r
								fadhy = k_adh_be*(rc_adh_be-r)*(dy)/r

								f_intx_b(l,i) = f_intx_b(l,i) + fadhx
								f_intx_e(j)   = f_intx_e(j) - fadhx

						        f_inty_b(l,i) = f_inty_b(l,i) + fadhy 
								f_inty_e(j)   = f_inty_e(j) - fadhy
						!!!!!********************************************************!!!!!
						!!! Traction-like friction force betwn border-polar depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_eb = velx_b(l,i)-velx_e(j)
								relative_vely_eb = vely_b(l,i)-vely_e(j)

								angl_diff_wrt_clk = datan2(-dx,dy) - &
								&datan2(relative_vely_eb,relative_velx_eb)
								angl_diff_wrt_antclk = datan2(dx,-dy) - &
								&datan2(relative_vely_eb,relative_velx_eb)

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_eb*(dy*1.0d0/r)&
								& + relative_vely_eb*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_eb*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_eb*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_b(l,i) = f_fricx_b(l,i) + f_fricx
								f_fricy_b(l,i) = f_fricy_b(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_eb*(-dy*1.0d0/r)&
								& + relative_vely_eb*(dx*1.0d0/r)

								f_fricx = - fric_coeff_eb*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_eb*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_b(l,i) = f_fricx_b(l,i) + f_fricx
								f_fricy_b(l,i) = f_fricy_b(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								end if
						!!!!!********************************************************!!!!!
								
       							end if
						    end do
						 end do
					    end do		


       		 !! Polar - Outer-ellipse interactions!!
					
					do l=1,mp
					    do i=1,np					    
						 do j=1,ellipse_bead_num
											
							dx = xe(j)-xp(l,i)
							dy = ye(j)-yp(l,i)
							r = dsqrt(dx*dx + dy*dy)

                      			if(r.lt.rc_rep_pe) then
								
				          		frepx = -k_rep_pe*(rc_rep_pe-r)*(dx)/r
					  			frepy = -k_rep_pe*(rc_rep_pe-r)*(dy)/r

				          			f_intx_p(l,i) = f_intx_p(l,i) + frepx 
				          			f_intx_e(j)   = f_intx_e(j) - frepx

				          			f_inty_p(l,i) = f_inty_p(l,i) + frepy 
				          			f_inty_e(j)   = f_inty_e(j) - frepy

						!!!!!********************************************************!!!!!
						!!! Traction-like friction force betwn border-ellipse depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_ep = velx_p(l,i)-velx_e(j)
								relative_vely_ep = vely_p(l,i)-vely_e(j)

								angl_diff_wrt_clk = datan2(-dx,dy) - &
								&datan2(relative_vely_ep,relative_velx_ep)
								angl_diff_wrt_antclk = datan2(dx,-dy) - &
								&datan2(relative_vely_ep,relative_velx_ep)

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_ep*(dy*1.0d0/r)&
								& + relative_vely_ep*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_ep*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_ep*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_p(l,i) = f_fricx_p(l,i) + f_fricx
								f_fricy_p(l,i) = f_fricy_p(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_ep*(-dy*1.0d0/r)&
								& + relative_vely_ep*(dx*1.0d0/r)

								f_fricx = - fric_coeff_ep*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_ep*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_p(l,i) = f_fricx_p(l,i) + f_fricx
								f_fricy_p(l,i) = f_fricy_p(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								end if
						!!!!!********************************************************!!!!!

                        		else if((r.le.rc_adh_pe).and.(r.ge.rc_rep_pe)) then
							
								!this_border_gets_mig_force(q,j) = .true.

								fadhx = k_adh_pe*(rc_adh_pe-r)*(dx)/r
								fadhy = k_adh_pe*(rc_adh_pe-r)*(dy)/r

								f_intx_p(l,i) = f_intx_p(l,i) + fadhx
								f_intx_e(j)   = f_intx_e(j) - fadhx

						        f_inty_p(l,i) = f_inty_p(l,i) + fadhy 
								f_inty_e(j)   = f_inty_e(j) - fadhy
						!!!!!********************************************************!!!!!
						!!! Traction-like friction force betwn border-polar depending on &
						!!! & their relative velocity in the closest _|_ dir. of &
						!!! & interaction vector. (SoftMat 2014, Karttunen, eq.6)

								relative_velx_ep = velx_p(l,i)-velx_e(j)
								relative_vely_ep = vely_p(l,i)-vely_e(j)

								angl_diff_wrt_clk = datan2(-dx,dy) - &
								&datan2(relative_vely_ep,relative_velx_ep)
								angl_diff_wrt_antclk = datan2(dx,-dy) - &
								&datan2(relative_vely_ep,relative_velx_ep)

								if (abs(angl_diff_wrt_clk).lt.&
								&abs(angl_diff_wrt_antclk)) then
								
								rel_vel_along_clk_vec = relative_velx_ep*(dy*1.0d0/r)&
								& + relative_vely_ep*(-dx*1.0d0/r)
 
								f_fricx = - fric_coeff_ep*rel_vel_along_clk_vec*(dy*1.0d0/r)
								f_fricy = - fric_coeff_ep*rel_vel_along_clk_vec*(-dx*1.0d0/r)
								f_fricx_p(l,i) = f_fricx_p(l,i) + f_fricx
								f_fricy_p(l,i) = f_fricy_p(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								else
								
								rel_vel_along_antclk_vec = relative_velx_ep*(-dy*1.0d0/r)&
								& + relative_vely_ep*(dx*1.0d0/r)

								f_fricx = - fric_coeff_ep*rel_vel_along_antclk_vec*(-dy*1.0d0/r)
								f_fricy = - fric_coeff_ep*rel_vel_along_antclk_vec*(dx*1.0d0/r)
								f_fricx_p(l,i) = f_fricx_p(l,i) + f_fricx
								f_fricy_p(l,i) = f_fricy_p(l,i) + f_fricy
								f_fricx_e(j)   = f_fricx_e(j) - f_fricx
								f_fricy_e(j)   = f_fricy_e(j) - f_fricy

								end if
						!!!!!********************************************************!!!!!
								
       							end if
						    end do
						 end do
					    end do		

				
	return
        end                
              



       subroutine force(j1)  !! Subroutine for calculating the body-force (intracellular force) for each of the border, polar and nurse cells

		use radii        
		use arrays
		use forces
		use parameters
        implicit none
        integer:: i,j,l,j1
        double precision:: l1,l2,dx1,dx2,dy1,dy2
    
       ! boundary conditions
	
        do l=1,mn
          xn(l,0) = xn(l,nn)
          yn(l,0) = yn(l,nn)
          xn(l,nn+1) = xn(l,1)
          yn(l,nn+1) = yn(l,1)
        end do

        do l=1,mb
          xb(l,0) = xb(l,nb)
          yb(l,0) = yb(l,nb)
          xb(l,nb+1) = xb(l,1)
          yb(l,nb+1) = yb(l,1)
        end do

        do l=1,mp
          xp(l,0) = xp(l,np)
          yp(l,0) = yp(l,np)
          xp(l,np+1) = xp(l,1)
          yp(l,np+1) = yp(l,1)
        end do
         
      !! For Nurse cells !!       
  	  do l=1,mn             ! loop for cell no.
              do i=1,nn         ! loop for beads in each cell

                   dx1 = xn(l,i-1)-xn(l,i)
                   dy1 = yn(l,i-1)-yn(l,i)
                   dx2 = xn(l,i)-xn(l,i+1)
                   dy2 = yn(l,i)-yn(l,i+1)

                   l1 = dsqrt(dx1*dx1 + dy1*dy1)

                   l2 = dsqrt(dx2*dx2 + dy2*dy2)

                   fxn(l,i)=kn*(l1-lo_n)*dx1/l1 - kn*(l2-lo_n)*dx2/l2 - & 
                           0.5d0*pn*lo_n*(dy1/l1 + dy2/l2) 
                     

                   fyn(l,i)=kn*(l1-lo_n)*dy1/l1 - kn*(l2-lo_n)*dy2/l2 + &  
                            0.5d0*pn*lo_n*(dx1/l1 + dx2/l2)

               end do
	       cmx_n(l) = sum(xn(l,1:nn))/nn
	       cmy_n(l) = sum(yn(l,1:nn))/nn
	  end do

  	  !! For Border cells !!  
	  do l=1,mb             ! loop for cell no.
              do i=1,nb         ! loop for beads in each cell

                   dx1 = xb(l,i-1)-xb(l,i)
                   dy1 = yb(l,i-1)-yb(l,i)
                   dx2 = xb(l,i)-xb(l,i+1)
                   dy2 = yb(l,i)-yb(l,i+1)

                   l1 = dsqrt(dx1*dx1 + dy1*dy1)

                   l2 = dsqrt(dx2*dx2 + dy2*dy2)

                   fxb(l,i)=kb*(l1-lo_b)*dx1/l1 - kb*(l2-lo_b)*dx2/l2 - & 
                           0.5d0*pb*lo_b*(dy1/l1 + dy2/l2) 
                     

                   fyb(l,i)=kb*(l1-lo_b)*dy1/l1 - kb*(l2-lo_b)*dy2/l2 + &  
                            0.5d0*pb*lo_b*(dx1/l1 + dx2/l2)


		   normal_force_unitvec_bx(l,i) = -(dy1/l1 + dy2/l2)/hypot(dy1/l1 + dy2/l2, dx1/l1 + dx2/l2)
		   normal_force_unitvec_by(l,i) =  (dx1/l1 + dx2/l2)/hypot(dy1/l1 + dy2/l2, dx1/l1 + dx2/l2)
               end do
	       cmx_b(l) = sum(xb(l,1:nb))/nb
	       cmy_b(l) = sum(yb(l,1:nb))/nb
	  end do

  	 !! For Polar cells !!  
	  do l=1,mp             ! loop for cell no.
              do i=1,np         ! loop for beads in each cell

                   dx1 = xp(l,i-1)-xp(l,i)
                   dy1 = yp(l,i-1)-yp(l,i)
                   dx2 = xp(l,i)-xp(l,i+1)
                   dy2 = yp(l,i)-yp(l,i+1)

                   l1 = dsqrt(dx1*dx1 + dy1*dy1)

                   l2 = dsqrt(dx2*dx2 + dy2*dy2)

                   fxp(l,i)=kp*(l1-lo_p)*dx1/l1 - kp*(l2-lo_p)*dx2/l2 - & 
                           0.5d0*pp*lo_p*(dy1/l1 + dy2/l2) 
                     

                   fyp(l,i)=kp*(l1-lo_p)*dy1/l1 - kp*(l2-lo_p)*dy2/l2 + &  
                            0.5d0*pp*lo_p*(dx1/l1 + dx2/l2)

               end do
	       cmx_p(l) = sum(xp(l,1:np))/np
	       cmy_p(l) = sum(yp(l,1:np))/np
	  end do
		 
!!##################### for the spring force calculation in the boundary ellipse beads ####################
  	  do i=2,ellipse_bead_num-1  ! loop for bead number; excluded for the 1st & last beads as they are supposed to be fixed.

                   dx1 = xe(i-1)-xe(i)
                   dy1 = ye(i-1)-ye(i)
                   dx2 = xe(i)-xe(i+1)
                   dy2 = ye(i)-ye(i+1)

                   l1 = dsqrt(dx1*dx1 + dy1*dy1)

                   l2 = dsqrt(dx2*dx2 + dy2*dy2)

                   fxe(i)=ke*(l1-lo_e)*dx1/l1 - kp*(l2-lo_e)*dx2/l2                     
                   fye(i)=ke*(l1-lo_e)*dy1/l1 - ke*(l2-lo_e)*dy2/l2 

	  end do
!!##########################################################################################################

       return
       end



       subroutine move_deterministic(dt)  !! Subroutine for coordinate updatation of the beads, when no active migration force is there
		use radii        
		use arrays
		use forces
       implicit none
       integer:: i,j,l,idum,j1
       double precision:: c,dt  
       
       c = 1.0d0      ! c is coeff. of viscous damping      

	xb_prev(:,:) = xb(:,1:nb)
	yb_prev(:,:) = yb(:,1:nb)
	xp_prev(:,:) = xp(:,1:np)
	yp_prev(:,:) = yp(:,1:np)
	xn_prev(:,:) = xn(:,1:nn)
	yn_prev(:,:) = yn(:,1:nn)
	xe_prev(:) = xe(1:ellipse_bead_num)
	ye_prev(:) = ye(1:ellipse_bead_num)

      	
       do l=1,mn
         do i=1,nn
            xn(l,i) = xn(l,i) + (fxn(l,i) + f_intx_n(l,i))*dt/c
            yn(l,i) = yn(l,i) + (fyn(l,i) + f_inty_n(l,i))*dt/c

	 if(xn(l,i).ge.(0.75d0*a)) then

            xn(l,i) = xn_prev(l,i)
            yn(l,i) = yn_prev(l,i)
	 end if	  

	 !!!!!************************** Velocity determination ***********************!!!!!
	 velx_n(l,i) = (xn(l,i) - xn_prev(l,i))/dt
	 vely_n(l,i) = (yn(l,i) - yn_prev(l,i))/dt
	 !!!!!*************************************************************************!!!!!
  	          
         end do
      end do

       do l=1,mb
         do i=1,nb
            xb(l,i) = xb(l,i) + (fxb(l,i) + f_intx_b(l,i))*dt/c
            yb(l,i) = yb(l,i) + (fyb(l,i) + f_inty_b(l,i))*dt/c

	 if(xb(l,i).ge.(0.75d0*a)) then

            xb(l,i) = xb_prev(l,i)
            yb(l,i) = yb_prev(l,i)
	 end if     

	 !!!!!************************** Velocity determination ***********************!!!!!
	 velx_b(l,i) = (xb(l,i) - xb_prev(l,i))/dt
	 vely_b(l,i) = (yb(l,i) - yb_prev(l,i))/dt
	 !!!!!*************************************************************************!!!!!
 
         end do
      end do

       do l=1,mp
         do i=1,np
            xp(l,i) = xp(l,i) + (fxp(l,i) + f_intx_p(l,i))*dt/c
            yp(l,i) = yp(l,i) + (fyp(l,i) + f_inty_p(l,i))*dt/c

	 if(xp(l,i).ge.(0.75d0*a)) then
            xp(l,i) = xp_prev(l,i) 
            yp(l,i) = yp_prev(l,i)
	 end if    

	 !!!!!************************** Velocity determination ***********************!!!!!
	 velx_p(l,i) = (xp(l,i) - xp_prev(l,i))/dt
	 vely_p(l,i) = (yp(l,i) - yp_prev(l,i))/dt
	 !!!!!*************************************************************************!!!!!
  
         end do
      end do

	!!!##################### movement of the outer ellipse beads ###########################!!!!
         do i=1,ellipse_bead_num
          !  xe(i) = xe(i) + (fxe(i) + f_intx_e(i))*dt/c
          !  ye(i) = ye(i) + (fye(i) + f_inty_e(i))*dt/c

	   ! if((i.eq.1).or.(i.eq.ellipse_bead_num)) then
	    !xe(i) = xe_prev(i)
	    !ye(i) = ye_prev(i)
	    !end if    

	 !!!!!************************** Velocity determination ***********************!!!!!
	 velx_e(i) = (xe(i) - xe_prev(i))/dt
	 vely_e(i) = (ye(i) - ye_prev(i))/dt
	 !!!!!*************************************************************************!!!!!
  
         end do
	!!!################################################################################!!!! 

      return
      end


       subroutine move_noise(j1,dt)  !! Subroutine for coordinate updatation of the beads, when active migration force is present
		use radii        
		use arrays
		use forces
		use parameters
       implicit none
       integer:: i,j,l,idum
       double precision:: ran2
       integer,intent(in):: j1
       double precision:: c,dt,noise_n,noise_b,noise_p,w,g,cm_x,co_x,r_x,cm_y,co_y,r_y,x_sat,&
			& conc_grad_x,conc_grad_y,epsilon_x,epsilon_y,s_x,s_y,force_coeff_x,force_coeff_y
       double precision:: tot_mig_force_x , mignoise_b
       double precision:: cluster_com_x_present, cluster_com_x_prev, cluster_com_y_present, cluster_com_y_prev,&
			& cluster_com_vel_x, cluster_com_vel_y, cluster_com_vel
       
	c = 1.0d0      ! c is coeff. of viscous damping  
	cm_x = 0.5d0  !0.28d0 !0.5d0
	co_x = 0.25d0 !0.25d0
	r_x = 0.1d0
	cm_y = 0.5d0
	co_y = 0.11d0
	r_y = 0.1d0 !0.1d0
	g = 4.00d0  
	epsilon_x = 1.00d0
	epsilon_y = 1.000d0
    force_coeff_x = 1.0d0
	force_coeff_y = 1.0d0 !0.01d0
	x_sat = 0.375d0*a 		
      
	!this_border_gets_mig_force(:,:) = .true.
	tot_mig_force_x = 0d0  !! for all the border cells

	xb_prev(:,:) = xb(:,1:nb)
	yb_prev(:,:) = yb(:,1:nb)
	xp_prev(:,:) = xp(:,1:np)
	yp_prev(:,:) = yp(:,1:np)
	xn_prev(:,:) = xn(:,1:nn)
	yn_prev(:,:) = yn(:,1:nn)
	xe_prev(:) = xe(1:ellipse_bead_num)
	ye_prev(:) = ye(1:ellipse_bead_num)


	cluster_com_x_prev = (sum(xb(:,1:nb)) + sum(xp(:,1:np)))/(mb*nb + mp*np)
	cluster_com_y_prev = (sum(yb(:,1:nb)) + sum(yp(:,1:np)))/(mb*nb + mp*np)

       do l=1,mn
         do i=1,nn

            xn(l,i) = xn(l,i) + (fxn(l,i) + f_intx_n(l,i))*dt/c + wn*noise_n*dsqrt(dt)/c
            yn(l,i) = yn(l,i) + (fyn(l,i) + f_inty_n(l,i))*dt/c + wn*noise_n*dsqrt(dt)/c

	 !if((xn(l,i).lt.0.0d0).or.(abs(yn(l,i)).ge.(b*dsqrt(0.25d0-((xn(l,i)-a/2.0d0)/a)**2)))) then

            !xn(l,i) = xn(l,i) - (fxn(l,i) + f_intx_n(l,i))*dt/c - wn*noise_n*dsqrt(dt)/c
            !yn(l,i) = yn(l,i) - (fyn(l,i) + f_inty_n(l,i))*dt/c - wn*noise_n*dsqrt(dt)/c
	 !else if ((xn(l,i).ge.(0.75d0*a)).and.(abs(yn(l,i)).lt.(b*dsqrt(0.25d0-((xn(l,i)-a/2.0d0)/a)**2)))) then
	 if(xn(l,i).ge.(0.75d0*a)) then
	    xn(l,i) = xn_prev(l,i) !! only x upgrade stops
	 end if     

	 !!!!!************************** Velocity determination ***********************!!!!!
	 velx_n(l,i) = (xn(l,i) - xn_prev(l,i))/dt
	 vely_n(l,i) = (yn(l,i) - yn_prev(l,i))/dt
	 !!!!!*************************************************************************!!!!!

         end do
      end do

      !CALL SYSTEM_CLOCK(COUNT=idum)
      !migforce_theta = (pi*1d0/1.5)*(2.0d0*ran2(idum) - 1.0d0) ! Angle for the direction of mig. force for all cells
      do l=1,mb
	!CALL SYSTEM_CLOCK(COUNT=idum)
        !CALL gasdev(noise_b,idum,mean_b,var_b)
        !CALL SYSTEM_CLOCK(COUNT=idum)
        !migforce_theta = (pi*1d0/3)*(2.0d0*ran2(idum) - 1.0d0) ! Angle for the direction of mig. force for l-th cell

        do i=1,nb
            !xb(l,i) = xb(l,i) + (fxb(l,i) + f_intx_b(l,i))*dt/c + Vo_b*cos(theta_b(l))*dt/c
            !yb(l,i) = yb(l,i) + (fyb(l,i) + f_inty_b(l,i))*dt/c + Vo_b*sin(theta_b(l))*dt/c
	    !if(xb(l,i).ge.x_sat) goto 5

	 if(this_border_gets_mig_force(l,i)) then
	 	!conc_grad_x = co_x + (cm_x-co_x)*(exp((xb(l,i)-0.75d0*a)/g))
	 	!s_x = epsilon_x * conc_grad_x * 1.0d0 / (1 + epsilon_x*conc_grad_x)  !! response function 
         	!xb(l,i) = xb(l,i) + (fxb(l,i) + f_intx_b(l,i) +(cm-co)*(exp((xb(l,i)-a)/g)) )*dt/c + wb*noise_b*dsqrt(dt)/c
            	!CALL SYSTEM_CLOCK(COUNT=idum)
            	!migforce_theta = (pi*1d0/1.5)*(2.0d0*ran2(idum) - 1.0d0) !Angle for the direction of mig. force for i-th bead of l-th cell

		s_x = co_x * cm_x * exp(r_x * xb(l,i))/(cm_x - co_x + co_x * exp(r_x * xb(l,i)))   !! Logistic function

	    	xb(l,i) = xb(l,i) + (fxb(l,i) + f_intx_b(l,i) + force_coeff_x*s_x*cos(migforce_theta(l)) &
			&+ f_fricx_b(l,i))*dt/c + wb*noise_b*dsqrt(dt)/c
	    	!goto 10 
!5	    	xb(l,i) = xb(l,i) + (fxb(l,i) + f_intx_b(l,i) +(cm-co)*(exp((x_sat-a)/g)) )*dt/c + wb*noise_b*dsqrt(dt)/c
!10          	yb(l,i) = yb(l,i) + (fyb(l,i) + f_inty_b(l,i))*dt/c + wb*noise_b*dsqrt(dt)/c

	  	if(xb(l,i).ge.(2.0d0*(0.75d0*a)/3)) then  !! when xb > 2L/3; L=length of egg chamber upto oocyte boundary=0.75*a
		!if(xb(l,i).ge.(1.0d0*(0.75d0*a)/2)) then  !! when xb > L/2; L=length of egg chamber upto oocyte boundary=0.75*a
	  	!if(xb(l,i).ge.10.0d0) then  !! when xb >= 10
	    	!conc_grad_y = co_y + yb(l,i)*(cm_y-co_y)/(b*dsqrt(0.25d0-((0.75d0*a-a/2.0d0)/a)**2))
	    	!s_y = epsilon_y * conc_grad_y * 1.0d0 / (1 + epsilon_y*conc_grad_y)  !! response function	
		s_y = co_y * cm_y * exp(r_y * yb(l,i))/(cm_y - co_y + co_y * exp(r_y * yb(l,i)))   !! Logistic function    
	    	yb(l,i) = yb(l,i) + (fyb(l,i) + f_inty_b(l,i) + force_coeff_y*s_y + force_coeff_x*s_x*sin(migforce_theta(l)) & 
			&+ f_fricy_b(l,i))*dt/c + wb*noise_b*dsqrt(dt)/c
		else
		yb(l,i) = yb(l,i) + (fyb(l,i) + f_inty_b(l,i) + force_coeff_x*s_x*sin(migforce_theta(l)) &
			&+ f_fricy_b(l,i))*dt/c + wb*noise_b*dsqrt(dt)/c
	  	end if

	  	!if(xb(l,i).lt.0.0d0) then
            		!xb(l,i) = xb(l,i) - (fxb(l,i) + f_intx_b(l,i))*dt/c - Vo_b*cos(theta_b(l))*dt/c
            		!yb(l,i) = yb(l,i) - (fyb(l,i) + f_inty_b(l,i))*dt/c - Vo_b*sin(theta_b(l))*dt/c
            		!xb(l,i) = xb(l,i) - (fxb(l,i) + f_intx_b(l,i) +(cm-co)*(exp((xb(l,i)-a)/g)) )*dt/c - wb*noise_b*dsqrt(dt)/c

	    	!xb(l,i) = xb_prev(l,i)
            	!yb(l,i) = yb_prev(l,i)
	  	!else if ((xb(l,i).lt.(2.0d0*(0.75d0*a)/3)).and.(xb(l,i).ge.0.0d0).and.(abs(yb(l,i)).ge.&
		 ! &(b*dsqrt(0.25d0-((xb(l,i)-a/2.0d0)/a)**2)))) then
	    	!xb(l,i) = xb_prev(l,i)
            	!yb(l,i) = yb_prev(l,i)
	  	!else if ((xb(l,i).ge.(2.0d0*(0.75d0*a)/3)).and.(abs(yb(l,i)).ge.&
		 ! &(b*dsqrt(0.25d0-((xb(l,i)-a/2.0d0)/a)**2)))) then
	    	!xb(l,i) = xb_prev(l,i)
            	!yb(l,i) = yb_prev(l,i)
	  	!else if ((xb(l,i).ge.(0.75d0*a)).and.(abs(yb(l,i)).lt.(b*dsqrt(0.25d0-((xb(l,i)-a/2.0d0)/a)**2)))) then
            	!xb(l,i) = xb(l,i) - (fxb(l,i) + f_intx_b(l,i) +(cm-co)*(exp((x_sat-a)/g)) )*dt/c - wb*noise_b*dsqrt(dt)/c
		if(xb(l,i).ge.(0.75d0*a)) then
	    	xb(l,i) = xb_prev(l,i)
            	!yb(l,i) = yb(l,i) - (fyb(l,i) + f_inty_b(l,i))*dt/c - wb*noise_b*dsqrt(dt)/c ! y may be upgraded though
		!else if ((xb(l,i).ge.(0.75d0*a)).and.(abs(yb(l,i)).ge.(b*dsqrt(0.25d0-((xb(l,i)-a/2.0d0)/a)**2)))) then
	    	!xb(l,i) = xb_prev(l,i)
            	!yb(l,i) = yb_prev(l,i)
	  	end if     
 
	  	!!!!!********* total x migration force and velocity calculation ************!!!!!!!!  
		if(xb_prev(l,i).ne.xb(l,i)) then
			tot_mig_force_x = tot_mig_force_x + force_coeff_x*s_x*cos(migforce_theta(l))
			migration_forcex(l,i) = force_coeff_x*s_x*cos(migforce_theta(l))
			migration_forcey(l,i) = force_coeff_x*s_x*sin(migforce_theta(l))
		else
			tot_mig_force_x = tot_mig_force_x
		end if
		!!!!!********************************************************!!!!!!!!

	 else
		xb(l,i) = xb(l,i) + (fxb(l,i) + f_intx_b(l,i) + f_fricx_b(l,i))*dt/c + wb*noise_b*dsqrt(dt)/c
         	yb(l,i) = yb(l,i) + (fyb(l,i) + f_inty_b(l,i) + f_fricy_b(l,i))*dt/c + wb*noise_b*dsqrt(dt)/c

	  	!if((xb(l,i).lt.0.0d0).or.(abs(yb(l,i)).ge.(b*dsqrt(0.25d0-((xb(l,i)-a/2.0d0)/a)**2)))) then
            		!xb(l,i) = xb(l,i) - (fxb(l,i) + f_intx_b(l,i))*dt/c - Vo_b*cos(theta_b(l))*dt/c
            		!yb(l,i) = yb(l,i) - (fyb(l,i) + f_inty_b(l,i))*dt/c - Vo_b*sin(theta_b(l))*dt/c
            		!xb(l,i) = xb(l,i) - (fxb(l,i) + f_intx_b(l,i) +(cm-co)*(exp((xb(l,i)-a)/g)) )*dt/c - wb*noise_b*dsqrt(dt)/c
	    	!xb(l,i) = xb_prev(l,i)
            	!yb(l,i) = yb_prev(l,i)
		!else if ((xb(l,i).ge.(0.75d0*a)).and.(abs(yb(l,i)).lt.(b*dsqrt(0.25d0-((xb(l,i)-a/2.0d0)/a)**2)))) then
		if(xb(l,i).ge.(0.75d0*a)) then
	    	xb(l,i) = xb_prev(l,i)
            		!yb(l,i) = yb(l,i) - (fyb(l,i) + f_inty_b(l,i))*dt/c - wb*noise_b*dsqrt(dt)/c ! y may be upgraded though	    
	  	end if      
	 end if

	 !!!!!************************** Velocity determination ***********************!!!!!
	 velx_b(l,i) = (xb(l,i) - xb_prev(l,i))/dt
	 vely_b(l,i) = (yb(l,i) - yb_prev(l,i))/dt
	 !!!!!*************************************************************************!!!!!
	 	
        end do
            !CALL SYSTEM_CLOCK(COUNT=idum)
            !CALL gasdev(mignoise_b,idum,mignoise_mean_b,mignoise_std_b)
            !migforce_theta(l) = migforce_theta(l) + mignoise_b*dsqrt(dt)
      end do
            CALL SYSTEM_CLOCK(COUNT=idum)
            CALL gasdev(mignoise_b,idum,mignoise_mean_b,mignoise_std_b)
            migforce_theta = migforce_theta + mignoise_b*dsqrt(dt)


       do l=1,mp
         do i=1,np

            xp(l,i) = xp(l,i) + (fxp(l,i) + f_intx_p(l,i) + f_fricx_p(l,i))*dt/c + wp*noise_p*dsqrt(dt)/c
            yp(l,i) = yp(l,i) + (fyp(l,i) + f_inty_p(l,i) + f_fricy_p(l,i))*dt/c + wp*noise_p*dsqrt(dt)/c

	 !if((xp(l,i).lt.0.0d0).or.(abs(yp(l,i)).ge.(b*dsqrt(0.25d0-((xp(l,i)-a/2.0d0)/a)**2)))) then

          !  xp(l,i) = xp_prev(l,i)
          !  yp(l,i) = yp_prev(l,i)

	 !else if ((xp(l,i).ge.(0.75d0*a)).and.(abs(yp(l,i)).lt.(b*dsqrt(0.25d0-((xp(l,i)-a/2.0d0)/a)**2)))) then
	  if(xp(l,i).ge.(0.75d0*a)) then
	    xp(l,i) = xp_prev(l,i)  !! only x upgrade stops
	 end if    

	 !!!!!************************** Velocity determination ***********************!!!!!
	 velx_p(l,i) = (xp(l,i) - xp_prev(l,i))/dt
	 vely_p(l,i) = (yp(l,i) - yp_prev(l,i))/dt
	 !!!!!*************************************************************************!!!!!
  
         end do
      end do


	!!!##################### movement of the outer ellipse beads ###########################!!!!
         do i=1,ellipse_bead_num
            !xe(i) = xe(i) + (fxe(i) + f_intx_e(i))*dt/c
            !ye(i) = ye(i) + (fye(i) + f_inty_e(i))*dt/c

	    !if((i.eq.1).or.(i.eq.ellipse_bead_num)) then
	    !xe(i) = xe_prev(i)
	    !ye(i) = ye_prev(i)
	    !end if    

	 !!!!!************************** Velocity determination ***********************!!!!!
	 velx_e(i) = (xe(i) - xe_prev(i))/dt
	 vely_e(i) = (ye(i) - ye_prev(i))/dt
	 !!!!!*************************************************************************!!!!!
  
         end do
	!!!################################################################################!!!! 

      !!!!!************ velocity of the cluster-com calculation *************!!!!!!
	cluster_com_x_present = (sum(xb(:,1:nb)) + sum(xp(:,1:np)))/(mb*nb + mp*np)
	cluster_com_y_present = (sum(yb(:,1:nb)) + sum(yp(:,1:np)))/(mb*nb + mp*np)	
	cluster_com_vel_x = (cluster_com_x_present - cluster_com_x_prev)/dt
	cluster_com_vel_y = (cluster_com_y_present - cluster_com_y_prev)/dt	
	cluster_com_vel = hypot(cluster_com_vel_x, cluster_com_vel_y)
      !!!!!******************************************************************!!!!!! 

     ! if(mod(j1,2000).eq.0) then
     ! write(20,*) j1, cluster_com_x_present/(0.75d0*a), cluster_com_vel_x, cluster_com_vel_y,&
	!	& sum(velx_b)*1d0/(mb*nb), sum(vely_b)*1d0/(mb*nb), sum(velx_p)*1d0/(mp*np), &
	!	& sum(vely_p)*1d0/(mp*np), tot_mig_force_x
     ! end if

      if(mod(j1,2000).eq.0) then
      write(20,*) j1, cluster_com_x_present/(0.75d0*a), cluster_com_vel_x, cluster_com_vel_y, cluster_com_vel,&
		& sum(velx_b)*1d0/(mb*nb), sum(vely_b)*1d0/(mb*nb), sum(velx_p)*1d0/(mp*np), &
		& sum(vely_p)*1d0/(mp*np), tot_mig_force_x, cluster_com_y_present
      end if
	

      return
      end


       subroutine initial_angle           !! Subroutine for initialization of noise-angle	
		use radii        
		use arrays
        implicit none
        double precision :: ran2
        integer :: idum,i,l

       !  idum = 564326
        CALL SYSTEM_CLOCK(COUNT=idum)         
        do l=1,mn
            theta_n(l) = pi*(2.0d0*ran2(idum) - 1.0d0) 
        end do
        do l=1,mb
            !theta_b(l) = pi*(2.0d0*ran2(idum) - 1.0d0) 
	     migforce_theta(l) = (pi*0d0/3)*(2.0d0*ran2(idum) - 1.0d0)
        end do
        do l=1,mp
            theta_p(l) = pi*(2.0d0*ran2(idum) - 1.0d0) 
        end do
        return
        end

 
!!!!!!!!/////// Gaussian random no. generator \\\\\\\\\\\\\\\\\\\\\\\\\\\

 Subroutine gasdev(g2,idum,mean,std) 
  !use numz
  Implicit none
      INTEGER,INTENT(INOUT)::idum
      DOUBLE PRECISION,INTENT(IN)::std,mean
      DOUBLE PRECISION,INTENT(OUT)::g2
      DOUBLE PRECISION::ran2,g1
      INTEGER:: iset
      DOUBLE PRECISION:: fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
    !  if (iset.eq.0) then
        DO
        v1=2.0d0*ran2(idum)-1.0d0
        v2=2.0d0*ran2(idum)-1.0d0
        rsq=v1**2+v2**2
        if((rsq<1.0d0).AND.(rsq/=0.0d0))EXIT
        ENDDO
        fac=std*SQRT(-2.0d0*log(rsq)/(rsq))
        g1=v1*fac+mean
        g2=v2*fac+mean      
      END subroutine gasdev


!!!!/////// Uniform Random number generators////////////////////////////////////

FUNCTION ran2(idum)
 ! USE numz
  IMPLICIT NONE
  DOUBLE PRECISION:: ran2
  !INTEGER,INTENT(inout),OPTIONAL::idum
  INTEGER,INTENT(inout)::idum
  !INTEGER :: idum
  INTEGER,PARAMETER::IM1=2147483563,IM2=2147483399,IMM1=IM1-1
  INTEGER,PARAMETER::IA1=40014,IA2=40692,IQ1=53668
  INTEGER,PARAMETER::IQ2=52774,IR1=12211,IR2=3791   
  INTEGER,PARAMETER::NTAB=32,NDIV=1+IMM1/NTAB
  DOUBLE PRECISION,PARAMETER::AM=1.0d0/IM1,EPS=1.2e-7,RNMX=1.0d0-EPS
  INTEGER::idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  IF (idum<0) THEN
     idum=MAX(-idum,1)
     idum2=idum
      DO j=NTAB+8,1,-1
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         IF (idum<0) idum=idum+IM1
         IF (j.LE.NTAB) iv(j)=idum
      ENDDO
      iy=iv(1)
   ENDIF
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   IF (idum<0) idum=idum+IM1
   k=idum2/IQ2
   idum2=IA2*(idum2-k*IQ2)-k*IR2
   IF (idum2<0) idum2=idum2+IM2
   j=1+iy/NDIV
   iy=iv(j)-idum2
   iv(j)=idum
   IF(iy.LT.1)iy=iy+IMM1
   ran2=MIN(AM*iy,RNMX)
   RETURN
 END FUNCTION ran2



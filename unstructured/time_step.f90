module time_step
  use time_step_split
  use time_step_unsplit
  implicit none

  integer :: meshAdapted
  data meshAdapted /0/

contains

  subroutine initialize_timestep
    use basic
    implicit none

    select case(isplitstep)
    case(0)
       call initialize_timestep_unsplit
       call assign_variables_unsplit
    case(1:2)
       call initialize_timestep_split
       call assign_variables_split
    end select

    if(myrank.eq.0 .and. iprint.ge.1) then
       print *, 'Index of U: ', u_i
       print *, 'Index of V: ', vz_i
       print *, 'Index of Chi: ', chi_i
       print *, 'Index of Psi: ', psi_i
       print *, 'Index of Bz: ', bz_i
       print *, 'Index of P: ', p_i
       print *, 'Index of n: ', den_i
       print *, 'Index of Pe: ', pe_i
       print *, 'Index of f: ', bf_i
       print *, 'Index of E: ', e_i
    end if
  end subroutine initialize_timestep

  subroutine finalize_timestep
    use basic
    implicit none

    select case(isplitstep)
    case(0)
       call finalize_timestep_unsplit
    case(1:2)
       call finalize_timestep_split
    end select
  end subroutine finalize_timestep

  subroutine clear_matrices
    use basic
    implicit none

    select case(isplitstep)
    case(0)
       call clear_matrices_unsplit
    case(1:2)
       call clear_matrices_split
    end select
  end subroutine clear_matrices
  
  subroutine finalize_matrices
    use basic
    implicit none

    select case(isplitstep)
    case(0)
       call finalize_matrices_unsplit
    case(1:2)
       call finalize_matrices_split
    end select
  end subroutine finalize_matrices


!============================================================
! ONESTEP
! ~~~~~~~
! advance solution to time ntime+1
!============================================================
subroutine onestep

  use basic
  use diagnostics
  use arrays
  use pellet
  use runaway_mod
  use kprad_m3dc1
  use transport_coefficients
  use particles
  use boundary_conditions
  use newvar_mod
  use auxiliary_fields
  use runaway_advection
  use resistive_wall

  implicit none

  integer :: calc_matrices, ivel_def, irepeat
  logical, save :: first_time = .true.

  real :: tstart, tend

#ifdef USE3D
  integer :: icount, maxiter
#endif

  vectype, dimension(dofs_per_element) :: dofs
   integer :: k, itri, izone, izone_index, j
    type(field_type), save ::   p_v
     integer :: ierr
   vectype, dimension(dofs_per_element, dofs_per_element) :: tempxx
       integer, dimension(dofs_per_element) :: imask
       type(vector_type), save :: b1_phi, phi_vec
       type(matrix_type), target, save  :: jphi_mat
       type(field_type), save :: psi_v
     real :: xmag2, xmag_pred, psim, dpsi
  real, save :: totcur0,xmag0



  ! Determine whether matrices should be re-calculated
  if(first_time &
       .or. (linear.eq.0 .and. mod(ntime,nskip).eq.0) &
       .or. (integrator.eq.1 .and. ntime.eq.1) .or. meshAdapted .eq. 1) then
     calc_matrices = 1
  else
     calc_matrices = 0
  endif

  ! Advance impurity charge states
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call kprad_ionize(dt)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_kprad = t_kprad + tend - tstart
  endif

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  transport coefficients"
  call define_transport_coefficients

  totcur0=totcur
  xmag0=xmag

! start of loop to repeat timestep if max iterations exceeded in 3D
  dt=dt/max_repeat
  do irepeat = 1, max_repeat

    ! calculate matrices for time advance
    if(calc_matrices.eq.1) then
       if(myrank.eq.0 .and. iprint.ge.1) print *, "Defining matrices"
       if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

       ! in linear case, eliminate second-order terms from matrix
       ivel_def = 1
       if(istatic.eq.1 .and. isplitstep.eq.1) ivel_def = 0
       call ludefall(ivel_def, idens, ipres, ipressplit, 1-iestatic)
       if(myrank.eq.0 .and. itimer.eq.1) then
          call second(tend)
          t_ludefall = t_ludefall + tend - tstart
       endif
       if(myrank.eq.0 .and. iprint.ge.1) print *, "Done defining matrices."
    endif

#ifdef USEPARTICLES
    if ((kinetic .eq. 1) .and. (linear .eq. 1)) then
       if (myrank .eq. 0 .and. itimer .eq. 1) call second(tstart)
       call ludefvel_nolin
       if (myrank .eq. 0 .and. itimer .eq. 1) then
          call second(tend)
          t_ludefall = t_ludefall + tend - tstart
       end if
    end if
#endif

    ! copy field data to time-advance vectors
    if(myrank.eq.0 .and. iprint.ge.1) print *,"Importing time advance vectors.."
    call import_time_advance_vectors

    ! advance time
    if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing time..."
    if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
    if(isplitstep.ge.1) then
       call step_split(calc_matrices)
    else
       call step_unsplit(calc_matrices)
    end if
    if(myrank.eq.0 .and. itimer.eq.1) then 
       call second(tend)
       if(iprint.ge.1) print *, "Time spent in *_step: ", tend-tstart
    end if

  ! copy time advance vectors to field data
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Exporting time advance vectors.."

! if(eqsubtract.eq.0) call subtract_axi    !DEBUG
  call export_time_advance_vectors
    !!if(nplanes.le.1) exit
    !! check if time step should be repeated

!#ifdef USE3D
    !maxiter = 0
    !do icount=1,maxnumofsolves
    !   maxiter = max(maxiter,int(kspits(icount)))
    !enddo
    !if(myrank.eq.0) write(*,'(A,2i5)') 'ksp_max, maxiter =', ksp_max, maxiter
    !if((maxiter .lt. ksp_max) .or. dtkecrit .le. 0) exit

    !dt = dt/2.
    !if(myrank.eq.0) write(*,'(A,e12.4)') 'Time step reduced by 2:  dt =', dt
!#else
    !exit
    if (linear.eq.1) calc_matrices=0
!#endif
  enddo
  dt=dt*max_repeat


  time = time + dt
  dtold = dt
  if(ntime.gt.1 .and. linear.eq.0) call variable_timestep

  call pellet_advance

!  call runaway_advance

  !! copy time advance vectors to field data
  !if(myrank.eq.0 .and. iprint.ge.2) print *, "Exporting time advance vectors.."

!! if(eqsubtract.eq.0) call subtract_axi    !DEBUG
  !call export_time_advance_vectors
  !if (irunaway .eq. 2) then
  !   call smooth_runaway
  !endif

#ifdef USEPARTICLES
  if(kinetic.eq.0) then
     if (first_time) call set_s1_0_mat
     call set_parallel_velocity
     call filter_velocity
  endif
#endif

#ifdef USEPARTICLES
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(kinetic.eq.1) call particle_step(dt*t0_norm)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_particle = t_particle + tend - tstart
  endif
#endif

   if (first_time) call create_field(p_v)
      p_v=0.
    if (first_time) then
       call create_vector(phi_vec,      1)
    call associate_field(psi_v,    phi_vec,      1)
       call create_vector(b1_phi,      1)

      call set_matrix_index(jphi_mat, 192)
       call create_mat(jphi_mat, 1, 1, icomplex, 1)
       do itri=1,local_elements()
          call define_element_quadrature(itri, int_pts_main, int_pts_tor)
          call define_fields(itri,FIElD_PSI,1,0)
          call get_flux_mask(itri, imask)
          tempxx = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_GS))
          call apply_boundary_mask(itri, 0, tempxx, imask)
          call insert_block(jphi_mat,itri,1,1,tempxx,MAT_ADD)
       enddo
       call flush(jphi_mat)

       call boundary_mag2(b1_phi, psi_v,jphi_mat)
       call finalize(jphi_mat)
    endif

  do itri=1,local_elements()
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,FIELD_PSI+FIELD_ETA,1,0)
     call get_zone_index(itri, izone_index)
     izone = zone_type(izone_index)
     temp79a=pst79(:,OP_GS)
      if((icd_source.eq.1).and.(izone.eq.ZONE_CONDUCTOR)) then
         do j=1,npoints
        !if (eta79(j,OP_1)<1e-4)  temp79a(j) = temp79a(j)+J_0cd * exp( -(x_79(j)-rzero_rekc-0.28)**2/w_cd**2 &
        !     - (z_79(j)-zzero_rekc)**2/w_cd**2 )
        !if (eta79(j,OP_1)<0.5*eta_wall)  then
        !   temp79a(j) = temp79a(j)*0.99+J_0cd*r_79(j)*(ntime*0.01+1)*exp( -(x_79(j)-rzero_rekc-0.275*cos(theta_rekc-phi_79(j)))**2/w_cd**2 &
        !     - (z_79(j)-zzero_rekc-0.275*sin(theta_rekc-phi_79(j)))**2/w_cd**2 )
       !else
           if ((sqrt((x_79(j)-rzero_rekc)**2+(z_79(j)-zzero_rekc)**2)>0.27).and.(sqrt((x_79(j)-rzero_rekc)**2+(z_79(j)-zzero_rekc)**2)<0.295)) then
           temp79a(j) = temp79a(j)*0.01+J_0cd*r_79(j)*(ntime*0.99)*(exp( -(x_79(j)-rzero_rekc-0.28*cos(theta_rekc-phi_79(j)))**2/w_cd**2 &
             - (z_79(j)-zzero_rekc-0.28*sin(theta_rekc-phi_79(j)))**2/w_cd**2 ))*tanh((0.295-sqrt((x_79(j)-rzero_rekc)**2+(z_79(j)-zzero_rekc)**2))/0.005)&
             *tanh((sqrt((x_79(j)-rzero_rekc)**2+(z_79(j)-zzero_rekc)**2)-0.27)/0.005)
       else
           ! if (sqrt((x_79(j)-rzero_rekc)**2+(z_79(j)-zzero_rekc)**2)>0.27) temp79a(j) = temp79a(j)*0.01
           temp79a(j) = temp79a(j)*0.01
          !temp79a(j) = temp79a(j)*0.98
              !temp79a(j)=temp79a(j)*0.9
              !temp79a(j)=0.1
           !endif
           !if (sqrt((x_79(j)-rzero_rekc)**2+(z_79(j)-zzero_rekc)**2)>0.30) temp79a(j)=temp79a(j)*0.5
       endif
        !if (eta79(j,OP_1)<1e-4) temp79a(j) = temp79a(j)+J_0cd*r_79(j)
        !if (eta79(j,OP_1)<1e-4) temp79a(j) = J_0cd*r_79(j)*ntime
       enddo
    endif
      if((icd_source.eq.1).and.(izone.eq.ZONE_VACUUM)) then
         do j=1,npoints
        !if (eta79(j,OP_1)<1e-4)  temp79a(j) = temp79a(j)+J_0cd * exp( -(x_79(j)-rzero_rekc-0.28)**2/w_cd**2 &
        !     - (z_79(j)-zzero_rekc)**2/w_cd**2 )
           temp79a(j) = temp79a(j)*0.01
        !if (eta79(j,OP_1)<1e-4) temp79a(j) = temp79a(j)+J_0cd*r_79(j)
        !if (eta79(j,OP_1)<1e-4) temp79a(j) = J_0cd*r_79(j)*ntime
       enddo
    endif


     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(p_v%vec,itri,1,dofs,VEC_ADD)
  end do
  call sum_shared(p_v%vec)
  call boundary_mag2(p_v%vec, psi_v)
  call newsolve(jphi_mat, p_v%vec, ierr)
  psi_field(1)=p_v

  call find_lcfs()
  !call magaxis(xmag,zmag,psi_field(1),psim,0,ierr)
  xmag2=xmag
  p_v=0.
  do itri=1,local_elements()
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,FIELD_PSI,1,0)
     !temp79a=(pipar79(:,OP_1)*0+piper79(:,OP_1)*3)/3.
     !temp79a=nfi79(:,OP_1)*te079(:,OP_1)*2
     !temp79a=nfi79(:,OP_1)*te079(:,OP_1)+p179(:,OP_1)-n179(:,OP_1)*te0
     temp79a=pst79(:,OP_1)
     !if (xmag>1.0875) temp79a=temp79a+0.0006*(xmag-1.0875)/0.2*x_79**2
     !if (xmag<1.0875) temp79a=temp79a+0.0006*(xmag-1.0875)/0.2*x_79**2
     !if (xmag>0.94) temp79a=temp79a+0.0001*(xmag-0.94)/0.2*x_79**2
     !if (xmag<0.94) temp79a=temp79a+0.0012*(xmag-0.94)/0.2*x_79**2
     !if (xmag>1.0) temp79a=temp79a+0.0001*(xmag-1.0)/0.2*x_79**2
     !if (xmag<1.0) temp79a=temp79a+0.0012*(xmag-1.0)/0.2*x_79**2
     temp79a=temp79a+1e-8*x_79**2
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(p_v%vec,itri,1,dofs,VEC_ADD)
  end do
  call newvar_solve(p_v%vec,mass_mat_lhs)
  !if(calc_matrices.eq.1) then
  !write(0,*) "111111111111111111"
  psi_field(1)=p_v
  ! Calculate all quantities derived from basic fields
  call find_lcfs()
  !write(0,*) xmag
  call derived_quantities(1)

  !call magaxis(xmag,zmag,p_v,psim,0,ierr)
  !dpsi=(1.0875-xmag2)/(xmag-xmag2)*1e-8
  !if (ntime<160) then
  !   xmag_pred=1.05-0.1/160*ntime
  !else
  !   xmag_pred=1.05-0.1
  !endif
     xmag_pred=1.088
   dpsi=(xmag_pred-xmag2)/(xmag-xmag2)*1e-8
  !dpsi=(-totcur/xmag0+totcur0/xmag0)*0.09
  ! dpsi=(totcur/xmag0-totcur0/xmag0)*0.22
  ! dpsi=(totcur-totcur0)*totcur0*(0.7)*(1./xmag)**2
  ! dpsi=(totcur-totcur0)*totcur0*(0.4)*(1./(xmag-0.87))
  dpsi=(totcur-totcur0)*totcur0*(1.0)
  ! dpsi=(totcur-totcur0)*(-0.05)
  dpsi=(totcur-totcur0)*(-0.04)
  if (nplanes.eq.1) dpsi=-5e-5
  ! if (myrank.eq.0) write(0,*) "333333333" ,dpsi/(totcur-totcur0)/totcur0
  if (myrank.eq.0) write(0,*) "333333333" ,dpsi
  ! write(0,*) xmag,xmag2,dpsi
  p_v=0.
  do itri=1,local_elements()
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,FIELD_PSI,1,0)
     !temp79a=(pipar79(:,OP_1)*0+piper79(:,OP_1)*3)/3.
     !temp79a=nfi79(:,OP_1)*te079(:,OP_1)*2
     !temp79a=nfi79(:,OP_1)*te079(:,OP_1)+p179(:,OP_1)-n179(:,OP_1)*te0
     temp79a=pst79(:,OP_1)
     !if (xmag>1.0875) temp79a=temp79a+0.0006*(xmag-1.0875)/0.2*x_79**2
     !if (xmag<1.0875) temp79a=temp79a+0.0006*(xmag-1.0875)/0.2*x_79**2
     !temp79a=temp79a+0.00006*(-0.0875)/0.2*x_79**2
     !if (xmag>0.94) temp79a=temp79a+0.0001*(xmag-0.94)/0.2*x_79**2
     !if (xmag<0.94) temp79a=temp79a+0.0012*(xmag-0.94)/0.2*x_79**2
     !if (xmag>1.0) temp79a=temp79a+0.0001*(xmag-1.0)/0.2*x_79**2
     !if (xmag<1.0) temp79a=temp79a+0.0012*(xmag-1.0)/0.2*x_79**2
     temp79a=temp79a+dpsi*x_79**2
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(p_v%vec,itri,1,dofs,VEC_ADD)
  end do
  call newvar_solve(p_v%vec,mass_mat_lhs)
  !if(calc_matrices.eq.1) then
  !write(0,*) "111111111111111111"
  psi_field(1)=p_v

  ! Calculate all quantities derived from basic fields
  call find_lcfs()
  call derived_quantities(1)
  if(ipellet_abl.gt.0) call pellet_shrink

  ! Advect impurity charge states
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call kprad_advect(dt)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_kprad = t_kprad + tend - tstart
  endif

  ! Conserve toroidal flux
  if(iconstflux.eq.1 .and. numvar.ge.2) then
     call conserve_flux
     tflux = tflux + gbound*area
  endif


  first_time = .false.
  meshAdapted = 0
end subroutine onestep

!======================================================================
! import_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! construct vectors for time-advance
!======================================================================
subroutine import_time_advance_vectors
  use basic
  implicit none

  select case(isplitstep)
  case(0)
     call import_time_advance_vectors_unsplit
  case(1:2)
     call import_time_advance_vectors_split
  end select

end subroutine import_time_advance_vectors


!======================================================================
! export_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! extract values from time-advance vectors
!======================================================================
subroutine export_time_advance_vectors
  use basic
  implicit none

  select case(isplitstep)
  case(0)
     call export_time_advance_vectors_unsplit
  case(1:2)
     call export_time_advance_vectors_split
  end select
end subroutine export_time_advance_vectors



!=============================
! scaleback
! ~~~~~~~~~
! rescale eigenfunction
!=============================
subroutine scaleback

  use basic
  use arrays
  use diagnostics
  use particles

  implicit none

  vectype, parameter :: scalefac = 1.e-10

  if(ekin.lt.max_ke .or. max_ke.eq.0) return
  if(myrank.eq.0) write(*,*) " =>solution scaled back at time", time

  call mult(field_vec, scalefac)
  if(i3d.eq.1) call mult(bf_field(1), scalefac)
  if(i3d.eq.1) call mult(bfp_field(1), scalefac)

#ifdef USEPARTICLES
  if(kinetic.eq.1) then
     call particle_scaleback(scalefac)
     ! epar=epar*scalefac**2
     ! if (hostrank==0) then
     !    do ipart=ipart_begin,ipart_end
     !       pdata(ipart)%wt=pdata(ipart)%wt*scalefac
     !       pdata(ipart)%wt2=pdata(ipart)%wt2*scalefac**2
     !    end do
     !    do ielm=ielm_min,ielm_max
     !       elfieldcoefs(ielm)%psiv1=elfieldcoefs(ielm)%psiv1*scalefac
     !       elfieldcoefs(ielm)%Bzv1=elfieldcoefs(ielm)%Bzv1*scalefac
     !       elfieldcoefs(ielm)%Bfv=elfieldcoefs(ielm)%Bfv*scalefac
     !    end do
     ! endif
     ! call update_particle_pressure
  end if
#endif
end subroutine scaleback


subroutine variable_timestep

  use basic
  use arrays
  use diagnostics

  implicit none
  include 'mpif.h'
  integer :: ierr

#ifdef USE3D
  integer :: maxiter, icount
#endif

if(dtkecrit.gt.0) then
  if(myrank.eq.0) then
!
! decrease timestep based on kinetic energy or ksp iterations
! but limit change to fraction dtfrac and bound by dtmin and dtmax
!
!   
       if(ekin.gt.dtkecrit) then
         dt = dtold/(1. + dtfrac)
       else

#ifdef USE3D
         maxiter = 0
         do icount=1,maxnumofsolves
           maxiter = max(maxiter,int(kspits(icount)))
         enddo

         if(maxiter .gt. ksp_warn) dt = dtold/(1. + dtfrac)
         if(maxiter .lt. ksp_min)  dt = dtold*(1. + dtfrac)
         if(iprint.ge.1) write(*,'("maxiter, ksp_warn, ksp_min",3i5)') maxiter, ksp_warn, ksp_min
#endif

         dt = max(dt,dtmin)
         dt = min(dt,dtmax)

         if(iprint.ge.1) write(*,'("dtold,dt,dtkecrit,ekin",1p4e12.4)') dtold,dt,dtkecrit,ekin
       endif



  endif ! on myrank.eq.0
  call MPI_bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

endif   ! on dtkecrit

end subroutine variable_timestep


end module time_step

! :::::::::: bc2amr ::::::::::::::::::::::::::::::::::::::::::::::;
!> \callgraph
!! \callergraph
!!  Take a grid patch with mesh widths **hx**,**hy**, of dimensions **nrow** by
!!  **ncol**,  and set the values of any piece of
!!  of the patch which extends outside the physical domain 
!!  using the boundary conditions. 
!!
!!  
!!   Specific to geoclaw:  extrapolates aux(i,j,1) at boundaries
!!   to constant.
!!
!!  ### Standard boundary condition choices for amr2ez in clawpack
!!
!!  At each boundary  k = 1 (left),  2 (right),  3 (bottom), 4 (top):
!!
!!  mthbc(k) =  
!!  * 0  for user-supplied BC's (must be inserted!)
!!  * 1  for zero-order extrapolation
!!  * 2  for periodic boundary conditions
!!  * 3  for solid walls, assuming this can be implemented
!!                   by reflecting the data about the boundary and then
!!                   negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
!!                   component of q.
!!  * 4  for sphere bcs (left half maps to right half of same side, and vice versa), as if domain folded in half
!!
!!  The corners of the grid patch are at 
!!     (xlo_patch,ylo_patch)  --  lower left corner
!!     (xhi_patch,yhi_patch) --  upper right corner
!!
!!  The physical domain itself is a rectangle bounded by
!!     (xlower,ylower)  -- lower left corner
!!     (xupper,yupper)  -- upper right corner
!!  
!   This figure below does not work with doxygen
!   the picture is the following: 
!  ____________________________________________________
! 
!                _____________________ (xupper,yupper)
!               |                     |  
!           ____|____ (xhi_patch,yhi_patch)   
!           |   |    |                |
!           |   |    |                |
!           |   |    |                |
!           |___|____|                |
!  (xlo_patch,ylo_patch) |            |
!               |                     |
!               |_____________________|   
!    (xlower,ylower)
!  ____________________________________________________
!!
!!
!>  Any cells that lie outside the physical domain are ghost cells whose
!!  values should be set in this routine.  This is tested for by comparing
!!  xlo_patch with xlower to see if values need to be set at the left
!   as in the figure above, 
!
!>  and similarly at the other boundaries.
!!  Patches are guaranteed to have at least 1 row of cells filled
!!  with interior values so it is possible to extrapolate. 
!!  Fix [trimbd()](@ref trimbd) if you want more than 1 row pre-set.
!!
!!  Make sure the order the boundaries are specified is correct
!!  so that diagonal corner cells are also properly taken care of.
!!
!!  Periodic boundaries are set before calling this routine, so if the
!!  domain is periodic in one direction only you
!!  can safely extrapolate in the other direction. 
!!
!!  Don't overwrite ghost cells in periodic directions!
!!
!! \param val data array for solution \f$q \f$ (cover the whole grid **msrc**)
!! \param aux data array for auxiliary variables 
!! \param nrow number of cells in *i* direction on this grid
!! \param ncol number of cells in *j* direction on this grid
!! \param meqn number of equations for the system
!! \param naux number of auxiliary variables
!! \param hx spacing (mesh size) in *i* direction
!! \param hy spacing (mesh size) in *j* direction
!! \param level AMR level of this grid
!! \param time setting ghost cell values at time **time**
!! \param xlo_patch left bound of the input grid
!! \param xhi_patch right bound of the input grid 
!! \param ylo_patch lower bound of the input grid 
!! \param yhi_patch upper bound of the input grid 
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

subroutine bc2amr(val,aux,nrow,ncol,meqn,naux, hx, hy, level, time,   &
                  xlo_patch, xhi_patch, ylo_patch, yhi_patch) 

    use amr_module, only: mthbc, xlower, ylower, xupper, yupper
    use amr_module, only: xperdom,yperdom,spheredom

    implicit none

    ! Input/Output
    integer, intent(in) :: nrow, ncol, meqn, naux, level
    real(kind=8), intent(in) :: hx, hy, time
    real(kind=8), intent(in) :: xlo_patch, xhi_patch
    real(kind=8), intent(in) :: ylo_patch, yhi_patch
    real(kind=8), intent(in out) :: val(meqn, nrow, ncol)
    real(kind=8), intent(in out) :: aux(naux, nrow, ncol)
    
    ! Local storage
    integer :: i, ii, j, ibeg, jbeg, nxl, nxr, nyb, nyt
    real(kind=8) :: hxmarg, hymarg

    real(kind=8) :: h_0,hu_0,h1,u_1,y,u_0,g

    g = 9.81d0

    hxmarg = hx * .01d0
    hymarg = hy * .01d0

    ! call interpolation routine
    ! call read_file_interpolate('bc.txt',time,hu_0,hx)
    call interpolate(time,hu_0)

    ! Use periodic boundary condition specialized code only, if only one 
    ! boundary is periodic we still proceed below
    if (xperdom .and. (yperdom .or. spheredom)) then
        return
    end if

    ! Each check has an initial check to ensure that the boundary is a real
    ! boundary condition and otherwise skips the code.  Otherwise 
    !-------------------------------------------------------
    ! Left boundary:
    !-------------------------------------------------------
    if (xlo_patch < xlower-hxmarg) then
        ! number of grid cells from this patch lying outside physical domain:
        nxl = int((xlower + hxmarg - xlo_patch) / hx)

        select case(mthbc(1))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                ! stop "A user defined boundary condition was not provided."

                do j = 1,ncol
                    y = ylo_patch + (j - 0.5d0) * hy
                    ! y  = ylower + (j - 0.5d0) * hy
                    if (abs(y-1950.0d0) <= 50.0d0) then
                        ! val(1,1,j) = max(val(1,1,j), 0.001d0)

                        ! if (hu_0 .ne. 0.0d0) then
                        
                        do i=1,nxl
                            if (val(1,1,j) < 1.d-4) then
                                h_0 = max((hu_0/sqrt(9.81d0))**(2.0d0/3.0d0), 0.001d0) 
                                val(1, i, j) = h_0
                                val(2,i,j) = hu_0
                                val(3,i,j) = 0.0d0
                            else

                                u_1 = val(2,1,j)/val(1,1,j)
                                ! h1 = val(1,1,j)
                                ! u_0 = hu_0/h1

                                if (hu_0 .ne. 0.0d0) then
                                    call newton_raphson(h_0,hu_0,val(1,1,j),u_1)
                                    if (h_0 > val(1,1,j)) then
                                        call two_shock(h_0,hu_0,val(1,1,j),u_1)
                                    end if
                                        
                                    ! do ii = 1,100
                                    !     h_0 = ((u_0 - u_1 + 2*sqrt(g*h1))**2)/(4.0d0*g)
                                    !     if (h_0 .le. 0) then
                                    !         h_0 = 0
                                    !         u_0 = 0
                                    !     else
                                    !         if (abs((hu_0/h_0) - u_0) < 1.0d-6) exit
                                    !             u_0 = hu_0/h_0
                                    !     end if
                                    ! enddo

                                    val(1,i,j) = h_0
                                    val(2,i,j) = hu_0   
                                    val(3,i,j) = 0.d0
                                else
                                    aux(:, i, j) = aux(:, 2 * nxl + 1 - i, j)
                                    val(:, i, j) = val(:, 2 * nxl + 1 - i, j)
                                    val(2, i, j) = -val(2, i, j)
                                end if
                            end if
                        end do
                    else
                        ! do j = 1, ncol
                            do i=1, nxl
                                aux(:, i, j) = aux(:, 2 * nxl + 1 - i, j)
                                val(:, i, j) = val(:, 2 * nxl + 1 - i, j)
                            end do
                            ! negate the normal velocity:
                            do i=1, nxl
                                val(2, i, j) = -val(2, i, j)
                            end do
                        ! end do
                    end if
                
                end do
            case(1) ! Zero-order extrapolation
                do j = 1, ncol
                    do i=1, nxl
                        aux(:, i, j) = aux(:, nxl + 1, j)
                        val(:, i, j) = val(:, nxl + 1, j)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do j = 1, ncol
                    do i=1, nxl
                        aux(:, i, j) = aux(:, 2 * nxl + 1 - i, j)
                        val(:, i, j) = val(:, 2 * nxl + 1 - i, j)
                    end do
                end do
                ! negate the normal velocity:
                do j = 1, ncol
                    do i=1, nxl
                        val(2, i, j) = -val(2, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop
        end select
    end if

    !-------------------------------------------------------
    ! Right boundary:
    !-------------------------------------------------------
    if (xhi_patch > xupper+hxmarg) then

        ! number of grid cells lying outside physical domain:
        nxr = int((xhi_patch - xupper + hxmarg) / hx)
        ibeg = max(nrow - nxr + 1, 1)

        select case(mthbc(2))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            case(1) ! Zero-order extrapolation
                do i = ibeg, nrow
                    do j = 1, ncol
                        aux(:, i, j) = aux(:, ibeg - 1, j)
                        val(:, i, j) = val(:, ibeg - 1, j)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do i=ibeg, nrow
                    do j = 1, ncol
                        aux(:, i, j) = aux(:, 2 * ibeg - 1 - i, j)
                        val(:, i, j) = val(:, 2 * ibeg - 1 - i, j)
                    end do
                end do
                ! negate the normal velocity:
                do i = ibeg, nrow
                    do j = 1, ncol
                        val(2, i, j) = -val(2, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! Bottom boundary:
    !-------------------------------------------------------
    if (ylo_patch < ylower - hymarg) then

        ! number of grid cells lying outside physical domain:
        nyb = int((ylower + hymarg - ylo_patch) / hy)

        select case(mthbc(3))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            
            case(1) ! Zero-order extrapolation
                do j = 1, nyb
                    do i = 1, nrow
                        aux(:,i,j) = aux(:, i, nyb + 1)
                        val(:,i,j) = val(:, i, nyb + 1)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do j = 1, nyb
                    do i = 1, nrow
                        aux(:,i,j) = aux(:, i, 2 * nyb + 1 - j)
                        val(:,i,j) = val(:, i, 2 * nyb + 1 - j)
                    end do
                end do
                ! negate the normal velocity:
                do j = 1, nyb
                    do i = 1, nrow
                        val(3,i,j) = -val(3, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! Top boundary:
    !-------------------------------------------------------
    if (yhi_patch > yupper + hymarg) then

        ! number of grid cells lying outside physical domain:
        nyt = int((yhi_patch - yupper + hymarg) / hy)
        jbeg = max(ncol - nyt + 1, 1)

        select case(mthbc(4))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."

            case(1) ! Zero-order extrapolation
                do j = jbeg, ncol
                    do i = 1, nrow
                        aux(:, i, j) = aux(:, i, jbeg - 1)
                        val(:, i, j) = val(:, i, jbeg - 1)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do j = jbeg, ncol 
                    do i = 1, nrow
                        aux(:, i, j) = aux(:, i, 2 * jbeg - 1 - j)
                        val(:, i, j) = val(:, i, 2 * jbeg - 1 - j)
                    end do
                end do
                ! negate the normal velocity:
                do j = jbeg, ncol
                    do i = 1, nrow
                        val(3, i, j) = -val(3, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

end subroutine bc2amr

subroutine interpolate(t,hu0)
    implicit none

    ! declare variables
    real(kind=8) :: t, hu0,zinterp,b
    real(kind=8), dimension(6) :: tt = [0,300,600,5160,5460,172800]
    real(kind=8), dimension(6) :: z = [0,0,20,20,0,0]

    integer :: i

    ! initialize zinterp to zero
    zinterp = 0.0d0

    ! check if t is within the tt range and set the value of zinterp
    ! write(*,*) 't_size' , size(tt)
    if (t < tt(1)) then
        zinterp = z(1)
    else if (t > tt(size(tt))) then
        zinterp = z(size(z))
    else
        do i = 1,size(tt)-1
            if (t >= tt(i) .and. t <= tt(i+1)) then
                zinterp = z(i) + (((z(i+1) - z(i)) / (tt(i+1) - tt(i))) * (t - tt(i)))
                exit
            end if
        end do
    end if

    ! write(*,*) 'The value of zinterp' , zinterp
    ! ----- end of linear interpolation ------------------------
    b = 100.0d0
    ! hu0 = zinterp/(b+2*dx)
    hu0 = zinterp/b 

end subroutine interpolate

subroutine read_file_interpolate(file_name, t, hu0,dx)

    implicit none

    ! declare variables
    character(len=*), intent(in) :: file_name
    real(kind=8), dimension(:), allocatable :: tt,z
    real(kind=8) :: t, zinterp, h0, h1 , u1,hu0, h,dx,dt
    character(len=100) :: line
    real(kind=8) :: hu_1,b,n,slope
    integer :: i,j,num_rows

    ! ----- read tt and z from a file -----------------------
    !  open the file for reading
    open(10,file=file_name,status='old')

    ! count the number of rows in the file
    num_rows = 0
    do 
        read(10,*,iostat=i) line
        if (i /= 0) exit
        num_rows = num_rows + 1
    end do

    ! allocate memory for tt and z
    allocate(tt(num_rows),z(num_rows))

    ! rewind the file
    rewind(10)

    ! read data
    do i = 1,num_rows
        read(10,*) tt(i), z(i)
        ! write(*,*) tt(i), z(i)
    end do

    ! close the file
    close(10)

    ! ------ Linear interpolation -----------------------------

    ! initialize zinterp to zero
    zinterp = 0.0d0

    ! check if t is within the tt range and set the value of zinterp
    write(*,*) 't_size' , size(tt), 'num_rows' , num_rows
    if (t < tt(1)) then
        zinterp = z(1)
    else if (t > tt(size(tt))) then
        zinterp = z(size(z))
    else
        do i = 1,size(tt)-1
            if (t >= tt(i) .and. t <= tt(i+1)) then
                zinterp = z(i) + (((z(i+1) - z(i)) / (tt(i+1) - tt(i))) * (t - tt(i)))
                exit
            end if
        end do
    end if

    ! write(*,*) 'The value of zinterp' , zinterp
    ! ----- end of linear interpolation ------------------------
    b = 100.0d0
    ! hu0 = zinterp/(b+2*dx)
    hu0 = zinterp/b 

    ! free up memory
    deallocate(tt,z)
end subroutine read_file_interpolate


subroutine newton_raphson(h0,hu0,h1,u1)

    implicit none

    ! declare variables
    real(kind=8) :: h0,h1,u1,x0,xn,tol,hu0
    real(kind=8) :: func,fxn,dfxn,dfunc_h0,F,g

    integer :: i, max_iter

    ! initialize variables
    tol = 1.0e-6 ! tolerance for convergence
    max_iter = 100 ! maximum number of iterations
    x0 = 0.01d0 ! initial guess for the inflow discharge
    F = 0.5d0 ! Froude number
    g = 9.81d0 ! gravitational acceleration

    ! solve Riemann invariants
    xn = (hu0/sqrt(g)*F)**(2.0d0/3.0d0)
    ! xn = h1
    do i = 1, max_iter
        fxn = func(hu0,xn,h1,u1)
        if (abs(fxn) < tol) then
            h0 = xn
            return 
        end if
        dfxn = dfunc_h0(hu0,xn,h1,u1)

        xn = xn - fxn/dfxn
        
    end do
    write(*,*) 'Newton-Raphson did not converge'
    xn = 0.0
    
end subroutine newton_raphson

real(kind=8) function func(hu0,h0,h1,u1)
    implicit none
    real(kind=8) :: hu0,h0,h1,u1
    real(kind=8) :: g
    
    g = 9.81d0 ! gravitational acceleration

    func = hu0/h0 - 2*sqrt(g*h0) - u1 + 2*sqrt(g*h1)

end function func

!  given hu0
real(kind=8) function dfunc_h0(hu0,h0,h1,u1)
    implicit none
    real(kind=8) :: hu0,h0,h1,u1
    real(kind=8) :: g
    
    g = 9.81d0 ! gravitational acceleration

    dfunc_h0 = -hu0/(h0**2) - sqrt(g/h0)

end function dfunc_h0

! given h0
real(kind=8) function dfunc_hu0(hu0,h0,h1,u1)
    implicit none
    real(kind=8) :: hu0,h0,h1,u1
    real(kind=8) :: g
    
    g = 9.81d0 ! gravitational acceleration

    dfunc_hu0 = 1/h0

end function dfunc_hu0

subroutine two_shock(h0,hu0,hr,ur)
    implicit none
    real(kind=8) :: hu0,h0,hr,ur
    real(kind=8) :: two_func,dtwo_func,tol
    real(kind=8) :: fxn,dfxn,xn,x0,epi,F,g

    integer :: i, max_iter

    ! initialize variables
    tol = 1.0e-8 ! tolerance for convergence
    max_iter = 100 ! maximum number of iterations
    ! x0 = 0.1d0 ! initial guess for the inflow depth
    epi = 1.0e-11 ! tolerance for the derivativeF = 0.50d ! Froude number
    g = 9.81d0 ! gravitational acceleration
    F = 1.0d0 ! Froude number

    ! solve Riemann invariants
    x0 = (hu0/sqrt(g)*F)**(2.0d0/3.0d0)
    ! x0 = hr

    ! NRM
    ! xn = x0
    do i = 1,max_iter
        fxn = two_func(hu0,x0,hr,ur)
        dfxn = dtwo_func(hu0,x0,hr,ur)

        if (abs(dfxn) < epi) stop

        xn = x0 - fxn/dfxn

        if (abs(xn-x0) <= tol) then
            h0 = xn
            return
        end if
        
        x0 = xn
    end do
    write (*,*) 'Newton-Raphson did not converge for two-shock solution'
    ! xn = 0.0

end subroutine two_shock

! 2-shock solution qr connects to q*
real(kind=8) function two_func(hu0,h0,hr,ur)
    implicit none
    real(kind=8) :: hu0,h0,hr,ur
    real(kind=8) :: g
    
    g = 9.81d0 ! gravitational acceleration

    two_func = hu0/h0 - ur - (h0 - hr)*sqrt((g/2.0d0)*(1.0d0/h0 + 1.0d0/hr)) 

end function two_func

! 2-shock derivative wrt h0
real(kind=8) function dtwo_func(hu0,h0,hr,ur)
    implicit none
    real(kind=8) :: hu0,h0,hr,ur
    real(kind=8) :: g,deno, num

    g = 9.81d0 ! gravitational acceleration

    num =  sqrt(g*(1.0d0/h0 + 1.0d0/hr))
    deno = 2*sqrt(2.0d0)*(h0**2)*num
    dtwo_func = -hu0/(h0**2) - num/sqrt(2.0d0) + g*(h0 - hr)/deno

end function dtwo_func


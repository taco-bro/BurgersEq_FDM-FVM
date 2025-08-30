MODULE common_data
    save
    integer, parameter          :: inputRe = 10000
    integer, parameter          :: inputCFL = 25
    integer, parameter          :: N = 25600
    real*8, parameter           :: Len=1.d0, u0=1.d0, nu=1.d0/real(inputRe), CFL=real(inputCFL)/1000.d0, t_tilde=3.2d0
    real*8, dimension(N+6)      :: u, uOld, x
    real*8, dimension(N+6)      :: LF_flux, D_flux, u_rec_mh, u_rec_ph, Force
    real*8, dimension(N+1)      :: u_interface, du_interface, x_interface
    real*8                      :: dt, dx, pi, Re
    real*8, dimension(N/2+1)    :: EnergySpectrum, DissipationSpectrum
END MODULE

MODULE common_fft
    save
    integer, parameter          :: nn = 25600
END MODULE

PROGRAM MAIN
    USE common_data
    implicit none
    integer                     :: i, j, it, nsteps
    integer, parameter          :: Ntime = 100
    real*8                      :: t, Tend
    real*8                      :: ke, epsilon, PP, etSAV, dkedt
    integer                     :: timeStart, timeEnd
    real*8                      :: ke_fft, epsilon_fft
    real*8                      :: thickness
    logical                     :: convergent_judge
    character(len=100)          :: fileName
    
    write(*,*) &
    '==================================================================='
    write(*,*)'Description:'
    write(*,*)'The code is used to calculate the non-linear Burgers Equation!'
    write(*,*)'Methedology:'
    write(*,*)'Finite-volume Method'
    write(*,*)'Reconstruction ---> fifth-ordter WENO in Lax-Friedrich flux'
    write(*,*)'Time Scheme:'
    write(*,*)'Total          ---> TVD third-order Runge-kutta scheme'
    write(*,*)'Numerical Integral Method:'
    write(*,*)'Trapezoidal rule'
    write(*,*)' '
    write(*,*)'CopyRight@Weiyuan Liang'
    write(*,*) &
    '==================================================================='

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                            Mesh setting
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dx = Len/real(N)
    dt = CFL*dx/u0
    Re = u0*Len/nu
    Tend = t_tilde*Len/u0

    write(*,*) &
    '==================================================================='
    write(*,*) 'Mesh Initializing !!!'
    write(*,*) 'N                ---> ', N
    write(*,*) 'dx               ---> ', dx
    write(*,*) 'dt               ---> ', dt
    write(*,*) 'Tend             ---> ', Tend
    write(*,*) 'Reynolds numbers ---> ', Re
    write(*,*) 'CFL              ---> ', CFL
    write(*,*) &
    '==================================================================='

    ! Cell center point
    do j = 1, N+6
        x(j) = dx/2.d0 + real(j-4)*dx
    end do

    ! Cell interface
    do j = 1, N+1
        x_interface(j) = x(j+2) + dx/2.d0
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                            file operation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(55, file='uSolution.dat', status='replace', action='write')
    write(fileName, '("spatialAveragedRe", I5.5, "N", I5.5, ".dat")') int(Re), N
    open(56, file=fileName, status='replace', action='write')
    write(fileName, '("uSteadySolutionRe", I5.5, "N", I5.5, ".dat")') int(Re), N
    open(57, file=fileName, status='replace', action='write')
    write(fileName, '("spatialAveragedInkSpaceRe", I5.5, "N", I5.5, ".dat")') int(Re), N
    open(58, file=fileName, status='replace', action='write')
    write(fileName, '("spectrumRe", I5.5, "N", I5.5, ".dat")') int(Re), N
    open(59, file=fileName, status='replace', action='write')
    ! open(60, file='output.dat', status='unknown', position='append', action='write')
    open(61, file='thicknessN1600.dat', status='unknown', position='append', action='write')

    pi = 4.0d0*atan(1.0d0)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                           Initialization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Initializing Velocity Field!'
    ! Initial condition
    u = 0.d0
    uOld = u

    etSAV = 0.d0
    ke=0.d0

    !
    convergent_judge = .False.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                           Time Loop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    t = 0.0d0
    nsteps = Tend/dt
    write(*,*) 'Time loop start!'
    call system_clock(timeStart)
    do it = 1, nsteps

        ! RungeKutta step
        CALL ThirdOrderRungeKuttaMethod()

        uOld = u

        t = t + dt

        ! Result output
        Call averageCal(ke, epsilon, PP, dkedt, etSAV)
        CALL calInFourierSpace(ke_fft, epsilon_fft) ! This subrutine must follow the last subrutine
        if ((mod(it, 20*Ntime) .eq. 0) .or. (it .eq. 1)) then
            ! do j = 1, N+1
            !     write(55, 100)t, x(j), u(j)
            ! end do
            write(56, 101)t, ke, epsilon, PP, (PP - epsilon), dkedt
            write(58, 102)t, ke_fft, epsilon_fft, (PP-epsilon_fft)
            do j = 1, N/2+1
                write(59, 103)t, real(j-1), EnergySpectrum(j), DissipationSpectrum(j) 
            end do

            if ((mod(it, 500*Ntime) .eq. 0) .or. (it .eq. 1)) then
                write(*,*)'Current time=', t, 'steps=', it
                write(*,*)'ke=', ke, 'epsilon=', epsilon, 'P=', PP, 'P-epsilon', (PP - epsilon)
                write(*,*)'***<u>=', etSAV, '***d<ke>/dt', dkedt
                write(*,*)'***Calculated by FFT:'
                write(*,*)'ke=', ke_fft, 'epsilon=', epsilon_fft, 'P-epsilon', (PP - epsilon_fft)
                write(*,*) &
                          '==================================================================='
            end if
        end if

        if (it .eq. nsteps) then
            write(*,*)'Current time=', t, 'steps=', it
                write(*,*)'ke=', ke, 'epsilon=', epsilon, 'P=', PP, 'P-epsilon', (PP - epsilon)
                write(*,*)'***<u>=', etSAV, '***d<ke>/dt', dkedt
                write(*,*)'***Calculated by FFT:'
                write(*,*)'ke=', ke_fft, 'epsilon=', epsilon_fft, 'P-epsilon', (PP - epsilon_fft)
                write(*,*) &
                          '==================================================================='
        end if


        ! Output the steady solution
        ! if (it .eq. nsteps) then
        !     do j = 4, N+3
        !         write(57, 100)t, x(j), u(j)
        !     end do
        ! end if

        if (it .eq. nsteps) then
            do j = 1, N+1
                write(57, 100)t, x_interface(j), u_interface(j)
            end do
            CALL calThickness(thickness)
            write(61, 105)Re, real(N), thickness
        end if

        ! judge whether converge
        if ((abs(PP-epsilon_fft) .LT. 2e-3) .and. (.not.convergent_judge)) then
            if (t .GT. 1) then
                write(*,*)'This case is converged !!!'
                convergent_judge = .True.
                ! write(60, 104)int(Re), N
            end if
        elseif ((abs(PP-epsilon_fft) .GT. 2e-3) .and. it .eq. nsteps) then
            print *, 'Finally, this case is not converged !!!'
        end if

    end do
    write(*,*) &
    '==================================================================='
    write(*,*) 'Time loop end!'
    call system_clock(timeEnd)
    write(*,*)'Run time = ', (timeEnd-timeStart)/1000, 's'

    100     format(2x, 3f16.12)
    101     format(2x, 6f16.12)
    102     format(2x, 4f16.12)
    103     format(2x, 4f20.12)
    104     format(2x, 2I8)
    105     format(2x, 3f20.12)

    close(55)
    close(56)
    close(57)
    close(58)
    close(59)
    ! close(60)
    close(61)


END PROGRAM MAIN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !          Define a some subrutine for calculating
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! u_{j-0.5}^+
subroutine REC_u_MinusHalf()
    USE common_data
    integer                 :: j
    real*8, parameter       :: eps = 1e-6, gamma1 = 3.d0/10.d0, & 
                               gamma2 = 3.d0/5.d0, gamma3 = 1.d0/10.d0
    real*8                  :: beta1, beta2, beta3
    real*8                  :: omega_tilde1, omega_tilde2, omega_tilde3
    real*8                  :: omega1, omega2, omega3
    real*8                  :: uu1, uu2, uu3


    ! initialize
    ! u_rec_mh = 0.d0

    ! reconstruction
    do j = 3, N+4
        beta1 = 13.d0/12.d0*((u(j) - 2.d0*u(j-1) + u(j-2))**2.d0) + 1.d0/4.d0*((u(j) - 4.d0*u(j-1) + 3.d0*u(j-2))**2.d0)
        beta2 = 13.d0/12.d0*((u(j+1) - 2.d0*u(j) + u(j-1))**2.d0) + 1.d0/4.d0*((u(j-1) - u(j+1))**2.d0)
        beta3 = 13.d0/12.d0*((u(j+2) - 2.d0*u(j+1) + u(j))**2.d0) + 1.d0/4.d0*((3.d0*u(j+2) - 4.d0*u(j+1) + u(j))**2.d0)
        omega_tilde1 = gamma1/((eps + beta1)**2.d0)
        omega_tilde2 = gamma2/((eps + beta2)**2.d0)
        omega_tilde3 = gamma3/((eps + beta3)**2.d0)
        omega1 = omega_tilde1 / (omega_tilde1 + omega_tilde2 + omega_tilde3)
        omega2 = omega_tilde2 / (omega_tilde1 + omega_tilde2 + omega_tilde3)
        omega3 = omega_tilde3 / (omega_tilde1 + omega_tilde2 + omega_tilde3)
        uu1 = -1.d0/6.d0*u(j-2) + 5.d0/6.d0*u(j-1) + 1.d0/3.d0*u(j)
        uu2 = 1.d0/3.d0*u(j-1) + 5.d0/6.d0*u(j) - 1.d0/6.d0*u(j+1)
        uu3 = 11.d0/6.d0*u(j) - 7.d0/6.d0*u(j+1) + 1.d0/3.d0*u(j+2)
        u_rec_mh(j) = omega1*uu1 + omega2*uu2 + omega3*uu3
    end do
end subroutine REC_u_MinusHalf

! u_{j+0.5}^-
subroutine REC_u_PlusHalf()
    USE common_data
    integer                 :: j
    real*8, parameter       :: eps = 1e-6, gamma1 = 1.d0/10.d0, & 
                               gamma2 = 3.d0/5.d0, gamma3 = 3.d0/10.d0
    real*8                  :: beta1, beta2, beta3
    real*8                  :: omega_tilde1, omega_tilde2, omega_tilde3
    real*8                  :: omega1, omega2, omega3
    real*8                  :: uu1, uu2, uu3

    ! initialize
    ! u_rec_ph = 0.d0

    ! reconstruction
    do j = 3, N+4
        beta1 = 13.d0/12.d0*((u(j-2) - 2.d0*u(j-1) + u(j))**2.d0) + 1.d0/4.d0*((u(j-2) - 4.d0*u(j-1) + 3.d0*u(j))**2.d0)
        beta2 = 13.d0/12.d0*((u(j-1) - 2.d0*u(j) + u(j+1))**2.d0) + 1.d0/4.d0*((u(j-1) - u(j+1))**2.d0)
        beta3 = 13.d0/12.d0*((u(j) - 2.d0*u(j+1) + u(j+2))**2.d0) + 1.d0/4.d0*((3.d0*u(j) - 4.d0*u(j+1) + u(j+2))**2.d0)
        omega_tilde1 = gamma1/((eps + beta1)**2.d0)
        omega_tilde2 = gamma2/((eps + beta2)**2.d0)
        omega_tilde3 = gamma3/((eps + beta3)**2.d0)
        omega1 = omega_tilde1 / (omega_tilde1 + omega_tilde2 + omega_tilde3)
        omega2 = omega_tilde2 / (omega_tilde1 + omega_tilde2 + omega_tilde3)
        omega3 = omega_tilde3 / (omega_tilde1 + omega_tilde2 + omega_tilde3)
        uu1 = 1.d0/3.d0*u(j-2) - 7.d0/6.d0*u(j-1) + 11.d0/6.d0*u(j)
        uu2 = -1.d0/6.d0*u(j-1) + 5.d0/6.d0*u(j) + 1.d0/3.d0*u(j+1)
        uu3 = 1.d0/3.d0*u(j) + 5.d0/6.d0*u(j+1) - 1.d0/6.d0*u(j+2)
        u_rec_ph(j) = omega1*uu1 + omega2*uu2 + omega3*uu3
    end do
end subroutine REC_u_PlusHalf

! f_{j+0.5}
subroutine REC_advection_flux()
    USE common_data
    integer                 :: j
    real*8                  :: alpha

    ! initialize
    ! LF_flux = 0.d0

    CALL REC_u_MinusHalf()
    CALL REC_u_PlusHalf()

    do j = 3, N+4
        ! In this case alpha = max|f'(u)| = max|u|
        ! alpha = maxval([abs(u_rec_mh(j+1)), abs(u_rec_ph(j))])
        if (abs(u_rec_mh(j+1)) >= abs(u_rec_ph(j))) then
            alpha = abs(u_rec_mh(j+1))
        else
            alpha = abs(u_rec_ph(j))
        endif
        ! Lax-Friedrich flux
        LF_flux(j) = 1.d0/2.d0*((u_rec_ph(j)**2.d0)/2.d0 + (u_rec_mh(j+1)**2.d0)/2.d0) + alpha1/2.d0*(u_rec_ph(j) - u_rec_mh(j+1))
    end do
end subroutine REC_advection_flux

subroutine REC_Diffusion_flux()
    USE common_data
    integer                 :: j

    ! initialize
    ! D_flux = 0.d0

    ! reconstruction
    do j = 3, N+4
        D_flux(j) = (15.d0*(u(j+1) - u(j)) - (u(j+2) - u(j-1)))/(12*Re*dx)
    end do
end subroutine REC_Diffusion_flux

subroutine REC_Forcing_term()
    USE common_data
    integer                 :: j

    ! initialize
    ! Force = 0.d0

    ! reconstruction
    do j = 3, N+4
        Force(j) = u0*u0/Len*(dsin(2.d0*pi*x(j)/Len)*dsin(pi*dx/Len)/(pi*dx/Len) & 
                   + dcos(4.d0*pi*x(j)/Len)*dsin(2.d0*pi*dx/Len)/(2.d0*pi*dx/Len))
    end do
end subroutine REC_Forcing_term

subroutine ThirdOrderRungeKuttaMethod()
    USE common_data
    integer                 :: j
    real*8, dimension(N+6)  :: FF_RK0, FF_RK1, FF_RK2
    real*8, parameter       :: alpha = 1.d0/2.d0, beta = 1.d0/6.d0, gamma = 2.d0/3.d0

    ! Initialize
    ! FF_RK0 = 0.d0
    ! FF_RK1 = 0.d0
    ! FF_RK2 = 0.d0

    ! RK first step
    CALL REC_advection_flux()
    CALL REC_Diffusion_flux()
    CALL REC_Forcing_term()
    do j = 4, N+3
        FF_RK0(j) = -(LF_flux(j) - LF_flux(j-1))/dx + (D_flux(j) - D_flux(j-1))/dx + Force(j)
    end do

    do j = 4, N+3
        u(j) = uOld(j) + dt*FF_RK0(j)
    end do

    u(1) = u(N+1)
    u(2) = u(N+2)
    u(3) = u(N+3)
    u(N+4) = u(4)
    u(N+5) = u(5)
    u(N+6) = u(6)

    ! RK second step
    CALL REC_advection_flux()
    CALL REC_Diffusion_flux()
    CALL REC_Forcing_term()

    do j = 4, N+3
        FF_RK1(j) = -(LF_flux(j) - LF_flux(j-1))/dx + (D_flux(j) - D_flux(j-1))/dx + Force(j)
    end do

    do j = 4, N+3
        u(j) = uOld(j) + dt/2.d0*(alpha*FF_RK0(j) + (1.d0-alpha)*FF_RK1(j))
    end do

    u(1) = u(N+1)
    u(2) = u(N+2)
    u(3) = u(N+3)
    u(N+4) = u(4)
    u(N+5) = u(5)
    u(N+6) = u(6)

    ! RK third step
    CALL REC_advection_flux()
    CALL REC_Diffusion_flux()
    CALL REC_Forcing_term()

    do j = 4, N+3
        FF_RK2(j) = -(LF_flux(j) - LF_flux(j-1))/dx + (D_flux(j) - D_flux(j-1))/dx + Force(j)
    end do

    do j = 4, N+3
        u(j) = uOld(j) + dt*(beta*FF_RK0(j) + gamma*FF_RK2(j) + (1.d0-beta-gamma)*FF_RK1(j))
    end do

    u(1) = u(N+1)
    u(2) = u(N+2)
    u(3) = u(N+3)
    u(N+4) = u(4)
    u(N+5) = u(5)
    u(N+6) = u(6)

end subroutine ThirdOrderRungeKuttaMethod


subroutine averageCal(ke, epsilon, PP, dkedt, etSAV)
    USE common_data
    implicit none
    real*8, intent(inout)         :: ke, epsilon, PP, dkedt, etSAV
    integer                       :: j
    real*8, dimension(N+1)        :: tempke, tempepsilon, tempPP
    real*8                        :: keOld

    ! x_interface = 0.d0
    u_interface = 0.d0
    du_interface = 0.d0

    ! reconstruct
    do j = 3, N+3
        u_interface(j-2) = (-u(j-1) + 7.d0*u(j) + 7.d0*u(j+1) - u(j+2)) /12.d0
        du_interface(j-2) = (u(j-1) - 15.d0*u(j) + 15.d0*u(j+1) - u(j+2)) / 12.d0/dx
    end do


    ! The numerical integral are using the Trapezoidal rule

    keOld = ke

    ! Averaged kinetic energy
    tempke = 0.d0
    do j = 1, N+1
        tempke(j) = u_interface(j)**2.d0 / 2.d0 
    end do
    ke = (tempke(1) + tempke(N+1))/2.d0*dx/Len
    do j = 2, N
        ke = ke + tempke(j)*dx/Len
    end do

    ! Averaged viscous dissipation rate
    tempepsilon = 0.d0
    do j = 1, N+1
        tempepsilon(j) = (du_interface(j))**2.d0
    end do
    epsilon = (tempepsilon(1) + tempepsilon(N+1))/2.d0*dx/Len*nu
    do j = 2, N
        epsilon = epsilon + tempepsilon(j)*dx/Len*nu
    end do

    ! Averaged energy input by the forcing term
    tempPP = 0.d0
    do j = 1, N+1
        tempPP(j) = u0*u0/Len*(dsin(2.d0*pi*x_interface(j)/Len) & 
                    + dcos(4.d0*pi*x_interface(j)/Len)) * u_interface(j)
    end do
    PP = (tempPP(1) + tempPP(N+1))/2.d0*dx/Len
    do j = 2, N
        PP = PP + tempPP(j)*dx/Len
    end do

    ! Every time step Spatial Averaged velocity
    etSAV = (u_interface(1) + u_interface(N+1))/2.d0*dx/Len
    do j = 2, N
        etSAV = etSAV + u_interface(j)*dx/Len
    end do

    dkedt = (ke - keOld)/dt

end subroutine averageCal

! 1-D Fast Fourier Transform
!--------------------------------------------------------------------- 
subroutine fft(cur)
    USE common_fft
    real, intent(inout)             :: cur(nn+2)

    integer                         :: INCX, JUMPX
    INTEGER, dimension(13)          :: IFAXX
	REAL, dimension(3*nn/2+1)       :: TRIGSX
    REAL, dimension(nn+1)           :: WORKX

    ! Define common blocks
	COMMON /TRIG/  TRIGSX
	COMMON /FAX/   IFAXX

    ! Arrays used for FFTs
    INCX = 1
    JUMPX = nn + 2


    ! Initialize arrays for FFT
    CALL FFTFAX(nn,IFAXX,TRIGSX)

    ! Fast Fourier Transform
    CALL FFT991(cur(1),WORKX,TRIGSX,IFAXX,INCX,JUMPX,nn,1,-1)

    do i = 1, nn+2
        if (abs(cur(i)) .LT. .00000003 ) then
            cur(i) = 0
        end if
    end do


end subroutine fft
!--------------------------------------------------------------------- 

! 1-D Inverse Fast Fourier Transform
!--------------------------------------------------------------------- 
subroutine ifft(cur)
    USE common_fft
    real, intent(inout)             :: cur(nn+2)

    integer                         :: INCX, JUMPX
    INTEGER, dimension(13)          :: IFAXX
	REAL, dimension(3*nn/2+1)       :: TRIGSX
    REAL, dimension(nn+1)           :: WORKX

    integer                         :: j

    ! Define common blocks
	COMMON /TRIG/  TRIGSX
	COMMON /FAX/   IFAXX

    ! Arrays used for FFTs
    INCX = 1
    JUMPX = nn + 2

    ! Initialize arrays for FFT
    CALL FFTFAX(nn,IFAXX,TRIGSX)

    ! Fast Fourier Transform
    CALL FFT991(cur(1),WORKX,TRIGSX,IFAXX,INCX,JUMPX,nn,1,+1)

    do i = 1, nn+2
        if (abs(cur(i)) .LT. .00003 ) then
            cur(i) = 0
        end if
    end do

    
end subroutine ifft
!--------------------------------------------------------------------- 

subroutine calSpectrum()
    USE common_data
    implicit none
    real, dimension(N+2)              :: ukArray
    integer                           :: j, i
    real*8                            :: coeff
    real*8, dimension(N+2)            :: doubleukArray(N+2) 
    real*8, dimension(N/2+1)          :: realPart(N/2+1), imagePart(N/2+1)

    ! Convert to single precision for fft
    ukArray(1:N+1) = u_interface(1:N+1)
    
    ! velocity component in Fourier space
    CALL fft(ukArray)
    doubleukArray(1:N+2) = ukArray(1:N+2)

    i = 1
    do j = 1, N+2, 2
        realPart(i) = doubleukArray(j)
        imagePart(i) = doubleukArray(j+1)
        i = i + 1
    end do

    ! Energy spectrum
    do j = 1, N/2+1
        EnergySpectrum(j) = realPart(j)**2.d0 + imagePart(j)**2.d0
    end do

    ! Dissipation spectrum
    coeff = 8.d0*nu*pi*pi/Len/Len
    do j = 1, N/2+1
        DissipationSpectrum(j) = coeff*real(j-1)*real(j-1) & 
                                 *(realPart(j)**2.d0 + imagePart(j)**2.d0)
    end do
    
end subroutine calSpectrum
!--------------------------------------------------------------------- 

! Calculate the kenimatic energy and dissipation rate based on spectrum
!--------------------------------------------------------------------- 
subroutine calInFourierSpace(ke_fft, epsilon_fft)
    USE common_data
    implicit none
    real*8, intent(inout)               :: ke_fft, epsilon_fft

    CALL calSpectrum()

    ke_fft = sum(EnergySpectrum)
    epsilon_fft = sum(DissipationSpectrum)
end subroutine calInFourierSpace
!--------------------------------------------------------------------- 


! Calculate the kenimatic energy and dissipation rate based on spectrum
!--------------------------------------------------------------------- 
subroutine calThickness(thickness)
    USE common_data
    implicit none
    real*8, intent(inout)               :: thickness
    real*8, parameter                   :: alpha = 1.d0
    integer                             :: j

    thickness = maxval(u_interface(2:N+1)-u_interface(1:N))

end subroutine calThickness
!--------------------------------------------------------------------- 

program main
	implicit none

	integer, parameter :: N = 128 !横方向格子数
	integer, parameter :: Nc = N / 2 !粗い格子数
	real(8), parameter :: X = 1.d0, Y = 1.d0 !計算領域サイズ
	real(8), parameter :: hf = X / N !細かい格子間隔
	real(8), parameter :: hc = X / Nc !粗い格子間隔
	real(8), parameter :: e0 = 8.85d-12 !真空の誘電率
	integer, parameter :: nu1 = 10, nu2 = 5 !smoothing のステップ数上限
	
	real(8) :: phi(N,N) !計算する電位
	real(8) :: rho(N,N) !電荷密度
	real(8) :: f(N,N) !右辺の値
	real(8) :: Lf(-1:1,-1:1) !L on fine grid
	real(8) :: Lc(-1:1,-1:1) !L on course grid
	real(8) :: Ifc(-1:1,-1:1) !Restriction operator
	real(8) :: Icf(-1:1,-1:1) !Interpolation operation

	call initialize(phi, rho, f, Lf, Lc, Ifc, Icf, N, hf, hc, e0)

	call smooth(phi, f, Lf, N, hf, nu1)

	call CGC(phi, Lf, Lc, Ifc, Icf, N, Nc, f, hc)

	call smooth(phi, f, Lf, N, hf, nu2)

	write(*,*) "end program"

	stop
contains

	subroutine initialize(phi, rho, f, Lf, Lc, Ifc, Icf, N, hf, hc, e0)
		implicit none

		real(8), intent(inout) :: phi(:,:), rho(:,:), f(:,:), Lf(:,:), Lc(:,:), Ifc(:,:), Icf(:,:)
		real(8), intent(in) :: e0
		integer, intent(in) :: N
		integer :: i, j

		phi(:) = 0.d0
		rho(:) = 0.d0
		
		do j = 1, N
			do i = 1, N
				if (((N/2-i)**2+(N/2-j)**2)*h**2 < 0.05**2) then
					rho(i,j) = 1.d-8
				end if
			end do
		end do

		f(:,:) = - rho(:,:) / e0

		Lf( :, :) = 0.d0
		Lf( 0, 0) = 4.d0
		Lf(-1, 0) = -1.d0
		Lf( 0,-1) = -1.d0
		Lf( 1, 0) = -1.d0
		Lf( 0, 1) = -1.d0
		Lf(:,:) = Lf / hf**2

		Lc( :, :) = 0.d0
		Lc( 0, 0) = 4.d0
		Lc(-1, 0) = -1.d0
		Lc( 0,-1) = -1.d0
		Lc( 1, 0) = -1.d0
		Lc( 0, 1) = -1.d0
		Lc(:,:) = Lc / hc**2

		Ifc( :, :) = 1.d0
		Ifc( 0, 0) = 4.d0
		Ifc(-1, 0) = 2.d0
		Ifc( 0,-1) = 2.d0
		Ifc( 1, 0) = 2.d0
		Ifc( 0, 1) = 2.d0
		Ifc(:,:) = Ifc / 16.d0

		Icf( :, :) = 1.d0
		Icf( 0, 0) = 4.d0
		Icf(-1, 0) = 2.d0
		Icf( 0,-1) = 2.d0
		Icf( 1, 0) = 2.d0
		Icf( 0, 1) = 2.d0
		Icf(:,:) = Icf / 4.d0

	end subroutine initialize

	subroutine smooth(phi, f, Lf, N, hf, nu)
		implicit none

		real(8), intent(in) :: f(:,:), Lf(:,:), hf
		integer, intent(in) :: N, nu
		real(8), intent(inout) :: phi(:,:)

		real(8) :: z(N,N)
		real(8), parameter :: omega = 1.d0 
		integer :: i, j, nt

		do nt = 1, nu
			do j = 2, N-1
				do i = 2, N-1
					if(mod(i+j,2) == 0) then
						z(i,j) = (hf**2 * f(i,j) + phi(i-1,j) + phi(i+1,j) + phi(i,j-1) + phi(i,j+1)) / 4.d0
					end if
				end do
			end do
			phi(:,:) = phi(:,:) + omega * (z(:,:) - phi(:,:))
			do j = 2, N-1
				do i = 2, N-1
					if(mod(i+j,2) /= 0) then
						z(i,j) = (hf**2 * f(i,j) + phi(i-1,j) + phi(i+1,j) + phi(i,j-1) + phi(i,j+1)) / 4.d0
					end if
				end do
			end do
			phi(:,:) = phi(:,:) + omega * (z(:,:) - phi(:,:))
		end do

	end subroutine smooth

	subroutine stencil(M, u, N)
		implicit none

		real(8), intent(in) :: M(:,:)
		real(8), intent(inout) :: u(:,:)

		integer :: i, j
		real(8) :: tmp(N,N)

		tmp(:,:) = u(:,:)
		do j = 2, N-1
			do i = 2, N-1
				tmp(i,j) = sum(M(:,:) * u(i-1:i+1,j-1:j+1)) !バグるかも
			end do
		end do
		u(:,:) = tmp(:,:)

	end subroutine stencil

	subroutine Restriction(M, uf, uc, N, Nc)
		implicit none
		
		real(8), intent(in) :: M(:,:), uf(:,:)
		real(8), intent(inout) :: uc(:,:)
		integer, intent(in) :: N, Nc

		integer :: i, j
		real(8) :: tmp(Nc,Nc)

		tmp(:,:) = uc(:,:)
		do j = 2, Nc-1
			do i = 2, Nc-1
				tmp(i,j) = sum(M(:,:) * uf((i-1)*2:(i-1)*2+2,(j-1)*2:(j-1)*2+2)) !バグるかも　
			end do
		end do
		uc(:,:) = tmp(:,:)
		
	end subroutine Restriction

	subroutine Interpolation(M, uf, uc, N, Nc)
		implicit none

		real(8), intent(in) :: M(:,:), uc(:,:)
		real(8), intent(inout) :: uf(:,:)
		integer, intent(in) :: N, Nc

		integer :: i, j
		real(8) :: tmp(N,N)

		tmp(:,:) = uf(:,:)
		do j = 2, N-1
			do i = 2, N-1
				if( (mod(i,2) /= 0) .and. (mod(j,2) /= 0) ) then
					tmp(i,j) = M(0,0) * uc((i-1)/2+1,(j-1)/2+1)
				elseif ( (mod(i,2) /= 0 ) .and. (mod(j,2) == 0) ) then
					tmp(i,j) = M(0,1) * uc((i-1)/2+1,j/2+1) + M(0,-1) * uc((i-1/2+1,j/2))
				elseif ( (mod(i,2) == 0 ) .and. (mod(j,2) /= 0) ) then
					tmp(i,j) = M(1,0) * uc(i/2+1,(j-1)/2+1) + M(-1,0) * uc(i/2,(j-1/2+1))
				else
					tmp(i,j) = M(1,1) * uc(i/2+1,j/2+1) + M(1,-1) * uc(i/2+1,j/2) &
										+ M(-1,1) * uc(i/2,j/2+1) + M(-1,-1) * uc(i/2,j/2)
				end if
			end do
		end do
		uf(:,:) = tmp(:,:)

	end subroutine Interpolation

	subroutine CGC(phi, Lf, Lc, Ifc, Icf, N, Nc, f, hc)
		implicit none
		
		real(8), intent(in) :: Lf(:,:), Lc(:,:), Ifc(:,:), Icf(:,:), f(:), hc
		real(8), intent(inout) :: phi(:,:)
		integer, intent(in) :: N, Nc

		real(8) :: df(N,N), dc(Nc,Nc), vf(N,N), vc(Nc,Nc), z(Nc,Nc)
		integer :: i, j, nt

		call stencil(Lf, phi, N)
		df(:,:) = f(:,:) - phi(:,:)

		call Restriction(Ifc, df, dc, N, Nc)

		z(:,:) = 0.d0

		do nt = 1, 10
			do j = 2, Nc-1
				do i = 2, Nc-1
					if(mod(i+j,2) == 0) then
						z(i,j) = (hc**2 * dc(i,j) + vc(i-1,j) + vc(i+1,j) + vc(i,j-1) + vc(i,j+1)) / 4.d0
					end if
				end do
			end do
			vc(:,:) = z(:,:)
			do j = 2, Nc-1
				do i = 2, Nc-1
					if(mod(i+j,2) /= 0) then
						z(i,j) = (hc**2 * dc(i,j) + vc(i-1,j) + vc(i+1,j) + vc(i,j-1) + vc(i,j+1)) / 4.d0
					end if
				end do
			end do
			vc(:,:) = z(:,:)
		end do

		call Interpolation(Icf, vf, vc, N, Nc)

		phi(:,:) = uf(:,:) + vf(:,:)

	end subroutine CGC


end program main
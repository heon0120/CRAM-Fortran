program cram_prototype
    implicit none

    ! --- 변수 선언 ---
    integer, parameter :: dp = kind(1.0d0) ! Double precision 정의
    
    ! CRAM-16 상수 (Maria Pusa 기준)
    real(dp), parameter :: alpha_0 = 2.124853754224559e-16_dp
    
    ! 8개의 복소수 극점 (theta)과 계수 (alpha)
    complex(dp), parameter :: theta(8) = [ &
        (-10.843917078656922_dp,  19.341377407004452_dp), &
        ( -5.264971343442646_dp,  16.339284281229984_dp), &
        ( -2.122113943206161_dp,  13.564178651833156_dp), &
        ( -0.081823122146031_dp,  10.963425915577646_dp), &
        (  1.383785461947230_dp,   8.498846321287950_dp), &
        (  2.383963462002334_dp,   6.136709875953683_dp), &
        (  2.986641571477435_dp,   3.847167733575691_dp), &
        (  3.238478440051268_dp,   1.597371510488696_dp) ]
        
    complex(dp), parameter :: alpha(8) = [ &
        (-1.144410113106596e-14_dp,  7.747306263795290e-15_dp), &
        ( 2.370597346889418e-12_dp, -2.898864771507421e-12_dp), &
        (-1.579624647318880e-09_dp,  2.457813735165842e-09_dp), &
        ( 4.619041285265691e-07_dp, -9.171960249712797e-07_dp), &
        (-6.527632612749448e-05_dp,  1.831707323871110e-04_dp), &
        ( 4.372573212169602e-03_dp, -1.884784732174033e-02_dp), &
        (-1.129207038753229e-01_dp,  9.588107759521360e-01_dp), &
        ( 4.305106518116345e+00_dp, -1.545937402031174e+01_dp) ]

    ! 시스템 행렬 및 벡터
    real(dp) :: A_matrix(3, 3)
    real(dp) :: N_init(3), N_next(3)
    real(dp) :: N_exact(3)
    
    ! 물리 변수
    real(dp) :: lambda_A, lambda_B, dt, log2

    ! --- 1. 시스템 설정 ---
    log2 = log(2.0_dp)
    lambda_A = log2 / 1.0_dp
    lambda_B = log2 / 100000.0_dp
    
    A_matrix = reshape([ &
        -lambda_A,  lambda_A,  0.0_dp, &
         0.0_dp,   -lambda_B,  lambda_B, &
         0.0_dp,    0.0_dp,    0.0_dp  &
    ], [3, 3])
    
    N_init = [1000.0_dp, 0.0_dp, 0.0_dp]
    dt = 10.0_dp

    ! --- 2. CRAM 알고리즘 구동 ---
    call cram16(A_matrix, N_init, dt, N_next)

    ! --- 3. 해석적 해 계산 ---
    N_exact(1) = N_init(1) * exp(-lambda_A * dt)
    N_exact(2) = N_init(1) * (lambda_A / (lambda_B - lambda_A)) * (exp(-lambda_A * dt) - exp(-lambda_B * dt))
    N_exact(3) = N_init(1) - N_exact(1) - N_exact(2)

    ! --- 4. 결과 출력 ---
    write(*, '(A, F6.1, A)') '=== ', dt, '초 후 결과 ==='
    write(*, '(A)') '핵종 | CRAM 결과    | 해석적 해    | 상대 오차'
    write(*, '(A, F12.4, A, F12.4, A, E12.4)') 'A (빠름) | ', N_next(1), ' | ', N_exact(1), ' | ', abs((N_next(1)-N_exact(1))/N_exact(1))
    write(*, '(A, F12.4, A, F12.4, A, E12.4)') 'B (느림) | ', N_next(2), ' | ', N_exact(2), ' | ', abs((N_next(2)-N_exact(2))/N_exact(2))
    write(*, '(A, F12.4, A, F12.4, A, E12.4)') 'C (안정) | ', N_next(3), ' | ', N_exact(3), ' | ', abs((N_next(3)-N_exact(3))/N_exact(3))

contains

    subroutine cram16(A, N0, dt, Nt)
        real(dp), intent(in)  :: A(:, :)
        real(dp), intent(in)  :: N0(:)
        real(dp), intent(in)  :: dt
        real(dp), intent(out) :: Nt(:)
        
        integer :: n, j, info, k  ! [수정 1] 루프 인덱스 k 명시적 선언 추가
        complex(dp), allocatable :: M_c(:, :), Z(:, :), X(:)
        integer, allocatable :: ipiv(:)
        
        ! [수정 3] 외부 LAPACK 루틴 명시적 선언
        external :: ZGESV
        
        n = size(N0)
        allocate(M_c(n, n), Z(n, n), X(n), ipiv(n))
        
        ! M = A * dt 를 복소수 행렬로 딱 한 번만 생성 (원본 보존용)
        M_c = cmplx(A * dt, 0.0_dp, kind=dp)
        
        Nt = alpha_0 * N0
        
        do j = 1, 8
            ! [수정 2] 매 루프마다 파괴되지 않은 원본 M_c를 Z로 복사
            Z = M_c
            
            ! 복사한 Z의 대각원소에서 theta 빼기
            do k = 1, n
                Z(k, k) = Z(k, k) - theta(j)
            end do
            
            X = cmplx(N0, 0.0_dp, kind=dp)
            
            ! 복소수 연립방정식 풀이: Z * X = X
            ! 주의: ZGESV 호출 후 Z는 LU 분해 결과로 파괴됨!
            call ZGESV(n, 1, Z, n, ipiv, X, n, info)
            
            if (info /= 0) then
                write(*, '(A, I5)') 'LAPACK ZGESV Error! Info: ', info
                stop
            end if
            
            Nt = Nt + 2.0_dp * real(alpha(j) * X)
            
            ! 루프 끝에 Z를 복구할 필요 없음 -> 다음 루프에서 M_c로 다시 덮어씀
        end do
        
        deallocate(M_c, Z, X, ipiv)
    end subroutine cram16

end program cram_prototype

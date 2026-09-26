! Diagnostic extraction using QE 7.4.1 XC/FFT routines.
! Derived from PW/src/v_of_rho.f90 and gradcorr.f90, Quantum ESPRESSO,
! distributed under the GNU General Public License (QE License, version 2).
! This helper reads the accepted state; it does not replace its potential.
SUBROUTINE stmfit_xc_components()
  USE kinds, ONLY: DP
  USE constants, ONLY: e2
  USE fft_base, ONLY: dfftp
  USE gvect, ONLY: ngm, g
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY: rho, rho_core, rhog_core
  USE xc_lib, ONLY: xc, xc_gcx
  USE dft_setting_params, ONLY: iexch, icorr, igcx, igcc, is_libxc, &
      rho_threshold_lda, rho_threshold_gga, grho_threshold_gga
  USE io_global, ONLY: ionode, stdout
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: rt(:,:), grad(:,:,:), ex(:), ec(:)
  REAL(DP), ALLOCATABLE :: vx(:,:), vc(:,:), v2x(:,:), v2c(:,:)
  REAL(DP), ALLOCATABLE :: h(:,:), dh(:), grad2(:), active(:)
  COMPLEX(DP), ALLOCATABLE :: rg(:)
  INTEGER :: k, ipol, unit
  IF (nspin /= 1) CALL errore('stmfit_xc_components','only nonmagnetic states',1)
  IF (iexch/=1 .OR. icorr/=4 .OR. igcx/=3 .OR. igcc/=4 .OR. ANY(is_libxc(1:4))) &
      CALL errore('stmfit_xc_components','only internal PBE drivers',1)
  ALLOCATE(rt(dfftp%nnr,1), grad(3,dfftp%nnr,1), rg(ngm))
  ALLOCATE(ex(dfftp%nnr), ec(dfftp%nnr), vx(dfftp%nnr,1), vc(dfftp%nnr,1))
  ALLOCATE(v2x(dfftp%nnr,1), v2c(dfftp%nnr,1))
  ALLOCATE(h(3,dfftp%nnr), dh(dfftp%nnr), grad2(dfftp%nnr), active(dfftp%nnr))
  rt(:,1) = rho%of_r(:,1) + rho_core(:)
  rg(:) = rho%of_g(:,1) + rhog_core(:)
  CALL stmfit_xc_write('rho_valence.dat', rho%of_r(:,1), 0)
  CALL stmfit_xc_write('rho_xc_input.dat', rt(:,1), 0)
  CALL xc(dfftp%nnr, 1, 1, rt, ex, ec, vx, vc)
  dh(:) = e2*(vx(:,1) + vc(:,1))
  CALL stmfit_xc_write('lda.dat', dh, 1)
  CALL fft_gradient_g2r(dfftp, rg, g, grad(:,:,1))
  DO k = 1, dfftp%nnr
    grad2(k) = grad(1,k,1)**2 + grad(2,k,1)**2 + grad(3,k,1)**2
    active(k) = 0._DP
    ! The internal gcxc driver compares the SQUARED gradient to its threshold.
    IF (ABS(rt(k,1))>rho_threshold_gga .AND. grad2(k)>grho_threshold_gga) active(k)=1._DP
  ENDDO
  CALL stmfit_xc_write('grad2_xc_input.dat', grad2, 0)
  CALL stmfit_xc_write('gga_active.dat', active, 0)
  CALL xc_gcx(dfftp%nnr, 1, rt, grad, ex, ec, vx, v2x, vc, v2c)
  dh(:) = e2*(vx(:,1) + vc(:,1))
  CALL stmfit_xc_write('gga_local.dat', dh, 1)
  DO k = 1, dfftp%nnr
    DO ipol = 1, 3
      h(ipol,k) = e2*(v2x(k,1)+v2c(k,1))*grad(ipol,k,1)
    ENDDO
  ENDDO
  CALL fft_graddot(dfftp, h, g, dh)
  dh(:) = -dh(:)
  CALL stmfit_xc_write('gga_divergence.dat', dh, 1)
  IF (ionode) THEN
    OPEN(NEWUNIT=unit, FILE='xc_thresholds.toml', STATUS='new', ACTION='write')
    WRITE(unit,'(a,es25.17)') 'rho_threshold_lda = ', rho_threshold_lda
    WRITE(unit,'(a,es25.17)') 'rho_threshold_gga = ', rho_threshold_gga
    WRITE(unit,'(a,es25.17)') 'gradient_squared_threshold_gga = ', grho_threshold_gga
    CLOSE(unit)
    WRITE(stdout,'(a)') 'XC_COMPONENTS_EXPORTED; no state, threshold or potential modified'
  ENDIF
  DEALLOCATE(rt, grad, rg, ex, ec, vx, vc, v2x, v2c, h, dh, grad2, active)
END SUBROUTINE stmfit_xc_components

SUBROUTINE stmfit_xc_write(filename, values, observable)
  USE kinds, ONLY: DP
  USE fft_base, ONLY: dfftp
  USE cell_base, ONLY: at, celldm, ibrav
  USE ions_base, ONLY: nat, ntyp => nsp, ityp, tau, zv, atm
  USE run_info, ONLY: title
  USE gvect, ONLY: gcutm
  USE gvecs, ONLY: dual
  USE gvecw, ONLY: ecutwfc
  USE scatter_mod, ONLY: gather_grid
  USE io_global, ONLY: ionode
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: observable
  REAL(DP), INTENT(IN) :: values(dfftp%nnr)
  REAL(DP), ALLOCATABLE :: gathered(:)
  ALLOCATE(gathered(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  CALL gather_grid(dfftp, values, gathered)
  IF (ionode) CALL plot_io(filename, title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
      dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
      gcutm, dual, ecutwfc, observable, atm, ityp, zv, tau, gathered, +1)
  DEALLOCATE(gathered)
END SUBROUTINE stmfit_xc_write

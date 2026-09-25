! Compiler/linker smoke test only: the Gamma job requires sequential MKL BLAS
! and LAPACK. Linking this on a login node does not run QE or an MPI calculation.
program qe_mkl_link_probe
  implicit none
  double precision :: a(1,1), b(1,1), c(1,1), ap(3), w(2), z(1,2), work(6)
  integer :: info
  a = 2.0d0
  b = 3.0d0
  c = 0.0d0
  call dgemm('N', 'N', 1, 1, 1, 1.0d0, a, 1, b, 1, 0.0d0, c, 1)
  if (abs(c(1,1) - 6.0d0) > 1.0d-12) error stop 'MKL BLAS smoke failed'
  ap = [2.0d0, 0.0d0, 3.0d0]
  call dspev('N', 'U', 2, ap, w, z, 1, work, info)
  if (info /= 0) error stop 'MKL LAPACK smoke returned an error'
  if (any(abs(w - [2.0d0, 3.0d0]) > 1.0d-12)) error stop 'MKL LAPACK smoke failed'
  print *, 'SEQUENTIAL_MKL_BLAS_AND_LAPACK_OK'
end program qe_mkl_link_probe

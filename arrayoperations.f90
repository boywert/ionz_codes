subroutine  fortran_subgrid_reionization_with_xfrac(nh,ngamma,xfrac,nxion,nion,len) bind(C)
  use iso_c_binding
  implicit none
  integer(C_INT) :: len
  real(C_FLOAT) :: nion
  real(C_FLOAT) :: nh(len),ngamma(len),xfrac(len),nxion(len)
  nxion(:) = min(nion*ngamma(:)/(nh(:)*xfrac(:)),1.0)
  return
end subroutine fortran_subgrid_reionization_with_xfrac

subroutine  fortran_subgrid_reionization(nh,ngamma,nxion,nion,len) bind(C)
  use iso_c_binding
  implicit none
  integer(C_INT) :: len
  real(C_FLOAT) :: nion
  real(C_FLOAT) :: nh(len),ngamma(len),nxion(len)
  nxion(:) = min(nion*ngamma(:)/nh(:),1.0)
  return
end subroutine fortran_subgrid_reionization

subroutine fortran_prepar_fftw_real_3d_with_xfrac(nh,nhs,xfrac,len) bind(C)
  use iso_c_binding
  implicit none
  integer(C_INT) :: len
  real(C_FLOAT) :: nh(len),nhs(len),xfrac(len)
  nhs(:) = nh(:)*xfrac(:)
  return
end subroutine fortran_prepar_fftw_real_3d_with_xfrac

subroutine fortran_multiply_constant_fftw_real_3d(input,constant,output,len) bind(C)
  use iso_c_binding
  implicit none
  integer(C_INT) :: len
  real(C_FLOAT) :: constant
  real(C_FLOAT) :: input(len),output(len)
  output(:) = constant*input(:)
  return
end subroutine fortran_multiply_constant_fftw_real_3d

subroutine fortran_condition_ionize(nh,ngamma,nxion,len) bind(C)
  use iso_c_binding
  implicit none
  integer(C_INT) :: len
  real(C_FLOAT) :: nh(len),ngamma(len),nxion(len)
  nxion(:) = max(nxion(:),real(min(int(ngamma(:)/nh(:)),1)))
  return
end subroutine fortran_condition_ionize

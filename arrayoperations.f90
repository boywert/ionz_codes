module arraysoperations
contains
  subroutine  fortran_subgrid_reionization_with_xfrac(nh,ngamma,xfrac,nxion,nion,len) 
    use iso_c_binding
    implicit none
    integer(C_INT) :: len
    real(C_FLOAT) :: nion
    real(C_FLOAT) :: nh(len),ngamma(len),xfrac(len),nxion(len)
    nxion(:) = min(nion*ngamma(:)/(nh(:)*(1.0-xfrac(:))),1.0)
    return
  end subroutine fortran_subgrid_reionization_with_xfrac

  subroutine  fortran_subgrid_reionization(nh,ngamma,nxion,nion,len)
    integer(kind=4) :: len
    real(kind=4) :: nion
    real(kind=4) :: nh(len),ngamma(len),nxion(len)
    nxion(:) = min(nion*ngamma(:)/nh(:),1.0)
    return
  end subroutine fortran_subgrid_reionization
end module arraysoperations

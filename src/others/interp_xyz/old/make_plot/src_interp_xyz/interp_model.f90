! loop each point 
do ipt = 1,npts

  ! parent element in source mesh
  ispec = eid(ipt)
  isdone = fin(ipt)

  ! skip point that is not located or is finised
  if (ispec==0 .or. isdone) then
    cycle
  end if

  ! local coordinate
  uvwi= uvw(:,ipt)

  ! interpolation coefficients
  call lagrange_poly(uvwi(1),hxi)
  call lagrange_poly(uvwi(2),heta)
  call lagrange_poly(uvwi(3),hgamma)

  ! loop each parameter
  do ipar = 1,NPARS
    val = 0.0
    ! weighted-sum of model values at all gll point in host element
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          val = val + hxi(i)*heta(j)*hgamma(k)* &
                      model(i,j,k,ispec,ipar)
        end do ! i
      end do ! j
    end do ! k
    modeli(ipt,ipar) = val
  end do ! ipar

  ! record finish status
  fin(ipt) = .true.

end do ! ipt

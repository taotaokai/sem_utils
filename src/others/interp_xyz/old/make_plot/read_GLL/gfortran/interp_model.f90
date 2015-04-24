write(*,*) '------ do interpolation ...'

! loop each element in mesh 2
do ispec2 = 1,nspec2

  ! loop each gll point in this element
  do igllz = 1,NGLLZ
    do iglly = 1,NGLLY
      do igllx = 1,NGLLX

        ! host element in mesh 1
        ispec1 = eid2(igllx,iglly,igllz,ispec2)
        ! local coordinate
        uvw = uvw2(:,igllx,iglly,igllz,ispec2)
        ! interpolation coefficients
        call lagrange_poly(uvw(1),hxi)
        call lagrange_poly(uvw(2),heta)
        call lagrange_poly(uvw(3),hgamma)

        !print *,'igll',igllx,iglly,igllz
        !print *,'ispec1',ispec1
        !print *,'uvw',uvw
        !print *,'hxi',hxi
        !print *,'heta',heta
        !print *,'hgamma',hgamma

        ! loop each parameter
        do ipar = 1,NPARS
          val = 0.0
          ! weighted-sum of model values at all gll point in host element
          do k1 = 1,NGLLZ
            do j1 = 1,NGLLY
              do i1 = 1,NGLLX
                val = val + hxi(i1)*heta(j1)*hgamma(k1)* &
                            model1(i1,j1,k1,ispec1,ipar)
              end do ! i1
            end do ! j1
          end do ! k1
        end do ! ipar

      end do ! igllx
    end do ! iglly
  end do ! igllz

end do ! ispec2

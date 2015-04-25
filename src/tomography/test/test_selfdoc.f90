subroutine selfdoc()

print '(a)', 'NAME'
print '(a)', ''
print '(a)', '   xsem_interp_xyz - interpolate SEM model on given points '
print '(a)', ''
print '(a)', 'SYNOPSIS'
print '(a)', ''
print '(a)', '   xsem_interp_xyz <mesh_dir> <model_dir> <xyz_list> <out_file>'
print '(a)', ''
print '(a)', 'DESCRIPTION'
print '(a)', ''
print '(a)', '   get model update from the kernel of this iteration and the (dmodel,dkernel)'
print '(a)', 'from the previous iterations using l-BFGS method '
print '(a)', ''
print '(a)', '   mesh_dir:  proc000***_reg1_solver_data.bin'
print '(a)', '   model_dir:  a text file containing lines of'
print '(a)', ''

end subroutine


program test

  call selfdoc()

end program
#### import the simple module from the paraview
from paraview.simple import *

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2024, 1276]

#====== read in model
betavtk = LegacyVTKReader(FileNames=['/Users/taok/NEChina/stage10.structure/iter07/model_smooth/vtk/reg_1_beta.vtk'])

#====== clip model
# clip bottom
clip1 = Clip(Input=betavtk)
clip1.Scalars = ['POINTS', 'beta']
clip1.ClipType = 'Sphere'
clip1.ClipType.Center = [0.0, 0.0, 0.0]
clip1.ClipType.Radius = 0.85

# clip top
clip2 = Clip(Input=clip1)
clip2.Scalars = ['POINTS', 'beta']
clip2.ClipType = 'Sphere'
clip2.ClipType.Center = [0.0, 0.0, 0.0]
clip2.ClipType.Radius = 0.98
clip2.InsideOut = 1

# clip west
clip3 = Clip(Input=clip2)
clip3.Scalars = ['POINTS', 'beta']
clip3.ClipType = 'Plane'
clip3.ClipType.Origin = [0.0, 0.0, 0.0]
clip3.ClipType.Normal = [-0.9358033, -0.21307332, -0.28084151]

# clip south
clip4 = Clip(Input=clip3)
clip4.Scalars = ['POINTS', 'beta']
clip4.ClipType = 'Plane'
clip4.ClipType.Origin = [0.0, 0.0, 0.0]
clip4.ClipType.Normal = [-0.17188596, -0.52488806, 0.83363526]

# clip north
clip5 = Clip(Input=clip4)
clip5.Scalars = ['POINTS', 'beta']
clip5.ClipType = 'Plane'
clip5.ClipType.Origin = [0.0, 0.0, 0.0]
clip5.ClipType.Normal = [-0.03341919, 0.86446756, -0.50157651]

# clip east
clip6 = Clip(Input=clip5)
clip6.Scalars = ['POINTS', 'beta']
clip6.ClipType = 'Plane'
clip6.ClipType.Origin = [0.0, 0.0, 0.0]
clip6.ClipType.Normal = [0.73525213, 0.446262, 0.51015148]

#====== contour dvs=2.5%
contour1 = Contour(Input=clip6)
contour1.ContourBy = ['POINTS', 'beta']
contour1.Isosurfaces = [0.025]

# show data in view
contour1Display = Show(contour1, renderView1)
contour1Display.Representation = 'Surface'

#====== depth slice 410/660-km

# get color transfer function/color map for 'beta'
betaLUT = GetColorTransferFunction('beta')
betaLUT.ApplyPreset('Warm to Cool', True)
betaLUT.RescaleTransferFunction(-0.03, 0.03)

## get opacity transfer function/opacity map for 'beta'
#betaPWF = GetOpacityTransferFunction('beta')
#betaPWF.RescaleTransferFunction(-0.03, 0.03)

# 660-km
slice1 = Slice(Input=clip6)
slice1.SliceType = 'Sphere'
slice1.SliceType.Center = [0.0, 0.0, 0.0]
slice1.SliceType.Radius = 1 - 660./6371

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'beta']
slice1Display.LookupTable = betaLUT
slice1Display.RescaleTransferFunctionToDataRange(False)
slice1Display.Opacity = 0.6
# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# 410-km
slice2 = Slice(Input=clip6)
slice2.SliceType = 'Sphere'
slice2.SliceType.Center = [0.0, 0.0, 0.0]
slice2.SliceType.Radius = 1 - 410./6371

# show data in view
slice2Display = Show(slice2, renderView1)
# trace defaults for the display properties.
slice2Display.Representation = 'Surface'
slice2Display.ColorArrayName = ['POINTS', 'beta']
#slice2Display.LookupTable = betaLUT
slice2Display.RescaleTransferFunctionToDataRange(False)
slice2Display.Opacity = 0.6
# show color bar/color legend
#slice2Display.SetScalarBarVisibility(renderView1, True)

#====== vertical slice

#--- south china 1
slicev1 = Slice(Input=clip1)
slicev1.SliceType = 'Plane'
slicev1.SliceType.Origin = [0.0, 0.0, 0.0]
slicev1.SliceType.Normal = [-0.17188596, -0.52488806, 0.83363526]

# show data in view
slicev1Display = Show(slicev1, renderView1)
# trace defaults for the display properties.
slicev1Display.Representation = 'Surface'
slicev1Display.ColorArrayName = ['POINTS', 'beta']
slicev1Display.LookupTable = betaLUT
slicev1Display.RescaleTransferFunctionToDataRange(False)
slicev1Display.Opacity = 0.6
# show color bar/color legend
#slicev1Display.SetScalarBarVisibility(renderView1, True)

#--- south china 2
slicev2 = Slice(Input=clip1)
slicev2.SliceType = 'Plane'
slicev2.SliceType.Origin = [0.0, 0.0, 0.0]
slicev2.SliceType.Normal = [-0.15295609, -0.56800113, 0.80868977]

# show data in view
slicev2Display = Show(slicev2, renderView1)
# trace defaults for the display properties.
slicev2Display.Representation = 'Surface'
slicev2Display.ColorArrayName = ['POINTS', 'beta']
slicev2Display.LookupTable = betaLUT
slicev2Display.RescaleTransferFunctionToDataRange(False)
slicev2Display.Opacity = 0.6
# show color bar/color legend
#slicev2Display.SetScalarBarVisibility(renderView1, True)

#====== plot plate boundary
pbvtk = LegacyVTKReader(FileNames=['/Users/taok/work_local/NEChina/3d_visualization/paraview/pb.vtk'])
pbvtkDisplay = Show(pbvtk, renderView1)
pbvtkDisplay.Representation = 'Surface'
pbvtkDisplay.LineWidth = 4.0
pbvtkDisplay.DiffuseColor = [1.0, 0.0, 0.0]

#====== plot geologic block 
blockvtk = LegacyVTKReader(FileNames=['/Users/taok/work_local/NEChina/3d_visualization/paraview/zhangpz_block.vtk'])
blockvtkDisplay = Show(blockvtk, renderView1)
blockvtkDisplay.Representation = 'Surface'
blockvtkDisplay.LineWidth = 2.0
blockvtkDisplay.DiffuseColor = [1.0, 1.0, 1.0]

#====== plot national boundaries
boundaryvtk = LegacyVTKReader(FileNames=['/Users/taok/work_local/NEChina/3d_visualization/paraview/boundary.vtk'])
boundaryvtkDisplay = Show(boundaryvtk, renderView1)
# trace defaults for the display properties.
boundaryvtkDisplay.Representation = 'Surface'
boundaryvtkDisplay.LineWidth = 3.0
boundaryvtkDisplay.DiffuseColor = [0.0, 0.0, 0.0]

#====== plot shorelines
shorelinevtk = LegacyVTKReader(FileNames=['/Users/taok/work_local/NEChina/3d_visualization/paraview/shoreline.vtk'])
shorelinevtkDisplay = Show(shorelinevtk, renderView1)
# trace defaults for the display properties.
shorelinevtkDisplay.Representation = 'Surface'
shorelinevtkDisplay.LineWidth = 3.0
shorelinevtkDisplay.DiffuseColor = [0.0, 0.0, 0.0]


SetActiveSource(contour1)
renderView1.ResetCamera()
RenderAllViews()
import sys, glob
# state file generated using paraview version 4.3.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [890, 544]
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.0, 0.0, 0.01200000034077675]
renderView1.StereoType = 0
renderView1.CameraPosition = [-0.031077302490787898, 0.010258316107124188, -0.008277791537963331]
renderView1.CameraFocalPoint = [0.014238981652039268, -0.007015440723142775, 0.02089584182818903]
renderView1.CameraViewUp = [0.27287744498438377, 0.9518485032663522, 0.1397223205089287]
renderView1.CameraParallelScale = 0.014647989944891485
renderView1.Background = [1.0, 1.0, 1.0]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Tecplot Reader'
out_fine000000000001dat = TecplotReader(FileNames=['out_fine-{}-00001.dat'.format(sys.argv[1])])
out_fine000000000001dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Contour'
contour1 = Contour(Input=out_fine000000000001dat)
contour1.ContourBy = ['POINTS', 'node']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Tecplot Reader'
out000000000001dat = TecplotReader(FileNames=['out-{}-00001.dat'.format(sys.argv[1])])
out000000000001dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Contour'
contour2 = Contour(Input=out000000000001dat)
contour2.ContourBy = ['POINTS', 'node']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from out000000000001dat
out000000000001datDisplay = Show(out000000000001dat, renderView1)
# trace defaults for the display properties.
out000000000001datDisplay.Representation = 'Outline'
out000000000001datDisplay.AmbientColor = [0.0, 0.0, 0.0]
out000000000001datDisplay.DiffuseColor = [0.0, 0.0, 0.0]
out000000000001datDisplay.BackfaceDiffuseColor = [0.0, 0.0, 0.0]
out000000000001datDisplay.CubeAxesColor = [0.0, 0.0, 0.0]
out000000000001datDisplay.ScalarOpacityUnitDistance = 0.00023291118055700018

# show data from out_fine000000000001dat
out_fine000000000001datDisplay = Show(out_fine000000000001dat, renderView1)
# trace defaults for the display properties.
out_fine000000000001datDisplay.Representation = 'Outline'
out_fine000000000001datDisplay.AmbientColor = [0.0, 0.0, 0.0]
out_fine000000000001datDisplay.DiffuseColor = [0.0, 0.0, 0.0]
out_fine000000000001datDisplay.Opacity = 0.48
out_fine000000000001datDisplay.BackfaceDiffuseColor = [0.0, 0.0, 0.0]
out_fine000000000001datDisplay.CubeAxesColor = [0.0, 0.0, 0.0]
out_fine000000000001datDisplay.ScalarOpacityUnitDistance = 0.00019628481446662218

# show data from contour1
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.AmbientColor = [0.0, 0.0, 0.0]
contour1Display.DiffuseColor = [0.0, 0.0, 0.0]
contour1Display.BackfaceDiffuseColor = [0.0, 0.0, 0.0]
contour1Display.CubeAxesColor = [0.0, 0.0, 0.0]

# show data from contour2
contour2Display = Show(contour2, renderView1)
# trace defaults for the display properties.
contour2Display.AmbientColor = [0.0, 0.0, 0.0]
contour2Display.DiffuseColor = [0.0, 0.0, 0.0]
contour2Display.Opacity = 0.48
contour2Display.BackfaceDiffuseColor = [0.0, 0.0, 0.0]
contour2Display.CubeAxesColor = [0.0, 0.0, 0.0]

WriteImage('t{}_vizIntestine.jpg'.format(sys.argv[1]))

Render()

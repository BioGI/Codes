# state file generated using paraview version 4.4.0
import sys
# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [890, 831]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.0005229531492432216, 0.0, 0.005985000170767307]
renderView1.StereoType = 0
renderView1.CameraPosition = [0.025476728323262428, -5.995498334801768e-06, 0.005132512405494053]
renderView1.CameraFocalPoint = [-0.007329497752889146, -5.995498334801768e-06, 0.005132512405494053]
renderView1.CameraViewUp = [0.0, 0.9999752697698804, -0.007032769629047244]
renderView1.CameraParallelScale = 0.008490876106446998
renderView1.Background = [0.32, 0.34, 0.43]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Tecplot Reader'
out_fine000070000001dat = TecplotReader(FileNames=['out_fine-{}-00001.dat'.format(sys.argv[1])])
out_fine000070000001dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Tecplot Reader'
out000070000002dat = TecplotReader(FileNames=['out-{}-00002.dat'.format(sys.argv[1])])
out000070000002dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Tecplot Reader'
out_fine000070000003dat = TecplotReader(FileNames=['out_fine-{}-00003.dat'.format(sys.argv[1])])
out_fine000070000003dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Tecplot Reader'
out000070000004dat = TecplotReader(FileNames=['out-{}-00004.dat'.format(sys.argv[1])])
out000070000004dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Tecplot Reader'
out_fine000070000004dat = TecplotReader(FileNames=['out_fine-{}-00004.dat'.format(sys.argv[1])])
out_fine000070000004dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Tecplot Reader'
out_fine000070000002dat = TecplotReader(FileNames=['out_fine-{}-00002.dat'.format(sys.argv[1])])
out_fine000070000002dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Group Datasets'
groupDatasets2 = GroupDatasets(Input=[out_fine000070000001dat, out_fine000070000002dat, out_fine000070000003dat, out_fine000070000004dat])

# create a new 'CSV Reader'
pardat000070000004csv = CSVReader(FileName=['pardat-{}-00004.csv'.format(sys.argv[1])])

# create a new 'CSV Reader'
pardat000070000002csv = CSVReader(FileName=['pardat-{}-00002.csv'.format(sys.argv[1])])

# create a new 'CSV Reader'
pardat000070000003csv = CSVReader(FileName=['pardat-{}-00003.csv'.format(sys.argv[1])])

# create a new 'CSV Reader'
pardat000070000001csv = CSVReader(FileName=['pardat-{}-00001.csv'.format(sys.argv[1])])

# create a new 'Tecplot Reader'
out000070000001dat = TecplotReader(FileNames=['out-{}-00001.dat'.format(sys.argv[1])])
out000070000001dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Tecplot Reader'
out000070000003dat = TecplotReader(FileNames=['out-{}-00003.dat'.format(sys.argv[1])])
out000070000003dat.DataArrayStatus = ['u', 'v', 'w', 'P', 'phi', 'node']

# create a new 'Group Datasets'
groupDatasets1 = GroupDatasets(Input=[out000070000001dat, out000070000002dat, out000070000004dat, out000070000003dat])

# create a new 'Clip'
clip3 = Clip(Input=groupDatasets1)
clip3.ClipType = 'Plane'
clip3.Scalars = ['POINTS', 'P']
clip3.Value = 0.0382100045681
clip3.InsideOut = 1

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [0.0, 0.00072, 0.00594000006094575]
clip3.ClipType.Normal = [0.0, 1.0, 0.0]

# create a new 'Clip'
clip4 = Clip(Input=clip3)
clip4.ClipType = 'Plane'
clip4.Scalars = ['POINTS', 'P']
clip4.Value = 0.0382100045681

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [0.0, -0.00072, 0.00594000006094575]
clip4.ClipType.Normal = [0.0, 1.0, 0.0]

# create a new 'Clip'
clip6 = Clip(Input=clip4)
clip6.ClipType = 'Plane'
clip6.Scalars = ['POINTS', 'P']
clip6.Value = 0.0382100045681
clip6.InsideOut = 1

# init the 'Plane' selected for 'ClipType'
clip6.ClipType.Origin = [-0.00072, 0.0, 0.00594000006094575]

# create a new 'Clip'
clip5 = Clip(Input=clip4)
clip5.ClipType = 'Plane'
clip5.Scalars = ['POINTS', 'P']
clip5.Value = 0.0382100045681

# init the 'Plane' selected for 'ClipType'
clip5.ClipType.Origin = [0.00072, 0.0, 0.00594000006094575]

# create a new 'Clip'
clip1 = Clip(Input=groupDatasets1)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'P']
clip1.Value = 0.0382100045681

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.0, 0.00072, 0.00594000006094575]
clip1.ClipType.Normal = [0.0, 1.0, 0.0]

# create a new 'Clip'
clip2 = Clip(Input=groupDatasets1)
clip2.ClipType = 'Plane'
clip2.Scalars = ['POINTS', 'P']
clip2.Value = 0.0382100045681
clip2.InsideOut = 1

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [0.0, -0.00072, 0.00594000006094575]
clip2.ClipType.Normal = [0.0, 1.0, 0.0]

# create a new 'Group Datasets'
coarseMesh = GroupDatasets(Input=[clip1, clip6, clip5, clip2])

# create a new 'Slice'
slice1 = Slice(Input=coarseMesh)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.000315, 0.0, 0.00594000006094575]

# create a new 'Slice'
slice2 = Slice(Input=groupDatasets2)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [0.000315, 0.0, 0.00598500017076731]

# create a new 'Group Datasets'
groupDatasets3 = GroupDatasets(Input=[slice2, slice1])

# create a new 'Calculator'
calculator1 = Calculator(Input=groupDatasets3)
calculator1.ResultArrayName = 'cOverCs'
calculator1.Function = 'phi/3.1e-6'

# create a new 'Contour'
contour1 = Contour(Input=groupDatasets3)
contour1.ContourBy = ['POINTS', 'node']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Group Datasets'
groupDatasets4 = GroupDatasets(Input=[pardat000070000003csv, pardat000070000002csv, pardat000070000004csv, pardat000070000001csv])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=groupDatasets4)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Glyph'
glyph1 = Glyph(Input=tableToPoints1,
    GlyphType='Sphere')
glyph1.Scalars = ['POINTS', 'rp']
glyph1.Vectors = [None, '']
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 0.000294190048165576
glyph1.GlyphTransform = 'Transform2'

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'cOverCs'
cOverCsLUT = GetColorTransferFunction('cOverCs')
cOverCsLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.75, 0.865003, 0.865003, 0.865003, 1.5, 0.705882, 0.0156863, 0.14902]
cOverCsLUT.ScalarRangeInitialized = 1.0
cOverCsLUT.LockDataRange = 1

# get opacity transfer function/opacity map for 'cOverCs'
cOverCsPWF = GetOpacityTransferFunction('cOverCs')
cOverCsPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.5, 1.0, 0.5, 0.0]
cOverCsPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from out000070000001dat
out000070000001datDisplay = Show(out000070000001dat, renderView1)
# trace defaults for the display properties.
out000070000001datDisplay.Representation = 'Outline'
out000070000001datDisplay.ScalarOpacityUnitDistance = 0.000276986101840083
out000070000001datDisplay.Visibility = 0

# show data from out000070000002dat
out000070000002datDisplay = Show(out000070000002dat, renderView1)
# trace defaults for the display properties.
out000070000002datDisplay.Representation = 'Outline'
out000070000002datDisplay.ScalarOpacityUnitDistance = 0.000276986101840083
out000070000002datDisplay.Visibility = 0

# show data from out000070000003dat
out000070000003datDisplay = Show(out000070000003dat, renderView1)
# trace defaults for the display properties.
out000070000003datDisplay.Representation = 'Outline'
out000070000003datDisplay.ScalarOpacityUnitDistance = 0.000276986100586373
out000070000003datDisplay.Visibility = 0

# show data from out000070000004dat
out000070000004datDisplay = Show(out000070000004dat, renderView1)
# trace defaults for the display properties.
out000070000004datDisplay.Representation = 'Outline'
out000070000004datDisplay.ScalarOpacityUnitDistance = 0.000276986103093793
out000070000004datDisplay.Visibility = 0

# show data from groupDatasets1
groupDatasets1Display = Show(groupDatasets1, renderView1)
# trace defaults for the display properties.
groupDatasets1Display.Representation = 'Outline'
groupDatasets1Display.ScalarOpacityUnitDistance = 0.000209993696163408
groupDatasets1Display.Visibility = 0

# show data from clip1
clip1Display = Show(clip1, renderView1)
# trace defaults for the display properties.
clip1Display.ScalarOpacityUnitDistance = 0.000227368143612786
clip1Display.Visibility = 0

# show data from clip2
clip2Display = Show(clip2, renderView1)
# trace defaults for the display properties.
clip2Display.ScalarOpacityUnitDistance = 0.000223509748959036
clip2Display.Visibility = 0

# show data from clip3
clip3Display = Show(clip3, renderView1)
# trace defaults for the display properties.
clip3Display.ScalarOpacityUnitDistance = 0.000234038049778651
clip3Display.Visibility = 0

# show data from clip4
clip4Display = Show(clip4, renderView1)
# trace defaults for the display properties.
clip4Display.ScalarOpacityUnitDistance = 0.000348295163774774
clip4Display.Visibility = 0

# show data from clip5
clip5Display = Show(clip5, renderView1)
# trace defaults for the display properties.
clip5Display.ScalarOpacityUnitDistance = 0.000350794479886739
clip5Display.Visibility = 0

# show data from clip6
clip6Display = Show(clip6, renderView1)
# trace defaults for the display properties.
clip6Display.ScalarOpacityUnitDistance = 0.000350794479886739
clip6Display.Visibility = 0

# show data from coarseMesh
coarseMeshDisplay = Show(coarseMesh, renderView1)
# trace defaults for the display properties.
coarseMeshDisplay.ScalarOpacityUnitDistance = 0.00020943668109054
coarseMeshDisplay.Visibility = 0

# show data from slice1
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.Visibility = 0

# show data from out_fine000070000002dat
out_fine000070000002datDisplay = Show(out_fine000070000002dat, renderView1)
# trace defaults for the display properties.
out_fine000070000002datDisplay.Representation = 'Outline'
out_fine000070000002datDisplay.ScalarOpacityUnitDistance = 5.89384803915182e-05
out_fine000070000002datDisplay.Visibility = 0


# show data from out_fine000070000003dat
out_fine000070000003datDisplay = Show(out_fine000070000003dat, renderView1)
# trace defaults for the display properties.
out_fine000070000003datDisplay.Representation = 'Outline'
out_fine000070000003datDisplay.ScalarOpacityUnitDistance = 5.89384803915182e-05
out_fine000070000003datDisplay.Visibility = 0

# show data from out_fine000070000004dat
out_fine000070000004datDisplay = Show(out_fine000070000004dat, renderView1)
# trace defaults for the display properties.
out_fine000070000004datDisplay.Representation = 'Outline'
out_fine000070000004datDisplay.ScalarOpacityUnitDistance = 5.89384929627893e-05
out_fine000070000004datDisplay.Visibility = 0

# show data from out_fine000070000001dat
out_fine000070000001datDisplay = Show(out_fine000070000001dat, renderView1)
# trace defaults for the display properties.
out_fine000070000001datDisplay.Representation = 'Outline'
out_fine000070000001datDisplay.ScalarOpacityUnitDistance = 5.89384835343359e-05
out_fine000070000001datDisplay.Visibility = 0

# show data from groupDatasets2
groupDatasets2Display = Show(groupDatasets2, renderView1)
# trace defaults for the display properties.
groupDatasets2Display.Representation = 'Volume'
groupDatasets2Display.ScalarOpacityUnitDistance = 0.000125188449218086
groupDatasets2Display.Visibility = 0

# show data from slice2
slice2Display = Show(slice2, renderView1)
# trace defaults for the display properties.
slice2Display.Visibility = 0

# show data from groupDatasets3
groupDatasets3Display = Show(groupDatasets3, renderView1)
# trace defaults for the display properties.
groupDatasets3Display.Visibility = 0

# show data from contour1
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.ColorArrayName = ['POINTS', '']

# show data from tableToPoints1
tableToPoints1Display = Show(tableToPoints1, renderView1)
# trace defaults for the display properties.
tableToPoints1Display.ColorArrayName = [None, '']
tableToPoints1Display.Visibility = 0

# show data from glyph1
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.ColorArrayName = ['POINTS', '']

# show data from calculator1
calculator1Display = Show(calculator1, renderView1)
# trace defaults for the display properties.
calculator1Display.ColorArrayName = ['POINTS', 'cOverCs']
calculator1Display.LookupTable = cOverCsLUT

# show color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# setup the color legend parameters for each legend in this view

# get color legend/bar for cOverCsLUT in view renderView1
cOverCsLUTColorBar = GetScalarBar(cOverCsLUT, renderView1)
cOverCsLUTColorBar.Position = [0.8876712170237978, 0.30248880377580395]
cOverCsLUTColorBar.Position2 = [0.11999999999999955, 0.4300000000000002]
cOverCsLUTColorBar.Title = 'cOverCs'
cOverCsLUTColorBar.ComponentTitle = ''
cOverCsLUTColorBar.RangeLabelFormat = '%3.1f'

WriteImage('t{}_vizPhiDrugRelease.jpg'.format(sys.argv[1]))

Render()

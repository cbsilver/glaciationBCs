# Adding the glacier (spline curve) to a loaded model.
# The PythonAnimation script is supposed to be run as macro when model is loaded.
# Hint: use paraview python shell (pvpython) for class inspections
# Possibly vtk-graphic objects may be added

from paraview.simple import *
import numpy
import glacierclass   # from Christian Silbermann

def glacier_points(t):    
    L_dom = 120000 #m
    L_max = 0.7*L_dom
    H_max = 200 #m
    x_0 = -0.5*L_dom
    t_0 = 0.0 #a
    t_1 = 1.0 #a
    my_glacier = glacierclass.glacier(L_dom, L_max, H_max, x_0, t_0, t_1)
    glacierLength=my_glacier.length(t)
    x_rel=numpy.array([0 , 0.5, 1])
    if glacierLength>0:
        x=x_0+x_rel*glacierLength
        y=my_glacier.local_height(x, t)*25.0
    else: # minimum length if glacier is not in the picture yet
        x=x_0+x_rel*(L_dom/1000.0)
        y=x_rel*0.0
    z=[0.0, 0.0, 0.0]    
    xyz=numpy.ravel([x,y,z], order='F')   # prepare points list for spline [x0,y0,z0, x1,yz,z1, ...]
    return xyz


PythonAnimationCue1 = PythonAnimationCue()
PythonAnimationCue1.Script= """
def start_cue(self):
    #View1 = GetActiveView()
    xyz=glacier_points(0.0)
    SplineSource1 = SplineSource()
    SplineSource1.ParametricFunction = 'Spline'
    SplineSource1.ParametricFunction.Points = xyz
    Show()
    splineSource1Display = GetDisplayProperties(SplineSource1)
    splineSource1Display.AmbientColor = [0.8, 0.8, 0.8]
    splineSource1Display.DiffuseColor = [0.8, 0.8, 0.8]
    splineSource1Display.LineWidth = 4.0

    xyz[1:8:3]=xyz[1::3]*0.0
    SplineSource2 = SplineSource()
    SplineSource2.ParametricFunction = 'Spline'
    SplineSource2.ParametricFunction.Points = xyz
    Show()
    splineSource2Display = GetDisplayProperties(SplineSource2)
    splineSource2Display.AmbientColor = [0.8, 0.8, 0.8]
    splineSource2Display.DiffuseColor = [0.8, 0.8, 0.8]
    splineSource2Display.LineWidth = 4.0

def tick(self): 
    animationScene1 = GetAnimationScene()
    time = animationScene1.TimeKeeper.Time
    print(time)
    xyz=glacier_points(time)
    Spline1 = FindSource("SplineSource1")
    Spline1.ParametricFunction.Points = xyz
    xyz[1:8:3]=xyz[1::3]*0.0 
    Spline2 = FindSource("SplineSource2")
    Spline2.ParametricFunction.Points = xyz
	
def end_cue(self): pass
"""


geolayers_2d_originalvtu = XMLUnstructuredGridReader(FileName=['/home/dominik/projects/ogs/video_tutorial/geolayers_2d_domain.vtu'])
geolayers_2d_originalvtu.CellArrayStatus = ['gmsh:physical', 'gmsh:geometrical']
geolayers_2d_originalvtu.PointArrayStatus = ['gmsh:dim_tags']
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.ViewSize = [1920, 600]
layout1 = GetLayout()
geolayers_2d_originalvtuDisplay = Show(geolayers_2d_originalvtu, renderView1, 'UnstructuredGridRepresentation')
renderView1.ResetCamera()
materialLibrary1 = GetMaterialLibrary()
ColorBy(geolayers_2d_originalvtuDisplay, ('CELLS', 'gmsh:physical'))
geolayers_2d_originalvtuDisplay.RescaleTransferFunctionToDataRange(True, False)

renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, -3750.0, 10000.0]
renderView1.CameraFocalPoint = [0.0, -3750.0, 0.0]
renderView1.CameraParallelScale = 21000
renderView1.Update()

scene = GetAnimationScene()
# dummy run to have window open
scene.NumberOfFrames = 3
scene.PlayMode = 'Sequence' 
scene.Duration = 0 

scene.Play()

any_key = input("command: ")

# actual run with glacier
scene.Cues.append(PythonAnimationCue1)
scene.NumberOfFrames = 100
scene.PlayMode = 'Real Time'
scene.Duration = 5

scene.Play()

any_key = input("-- end --")


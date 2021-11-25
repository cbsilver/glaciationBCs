# glaciationBCs

Python package for capturing glacial cycles via boundary conditions for hydro-geological simulations, especially with the multi-field Simulator OpenGeoSys

## Classes:
* `glacierclass.py` with glacier properties and glacier induced effects (glacier height evolution, deflection under glacier, ...)
* `crustclass.py` with crustal properties given by external displacement field of the lithosphere
* `airclass.py` with atmospheric properties (evolving air temperature and pressure)

* `pythonBCsOGS.py` contains all the BC objects for OGS using the classes for glacier, crust and air

## Tools:
* `glaciation_test.py` simple tests for checking functionality

## Tools:
* `glacier_animate_paraview.py` animation tool for the glacier's shape in ParaView
* `deform_mesh.py` postprocessing tool deforming the mesh according to some given displacement field

## Branches:
* `main` for BCs only on parts of the boundary w.r.t. glaciation

## TODO:
* extend to multiple glacial cycles with hysteretic effects

### how to Simulate
* (1) blockMesh
* (2) cp -r 0.orig 0
* (3) setUFields
* (4) setRhoFields
* (5) foamToVTK -name VTK-0s
* (6) decomposePar
* (7) mpirun -n 36 ISWFoam -parallel &>log
* (8) reconstructPar -time 10,20,30,40,50
* (9) foamToVTK -name VTK-0s-50s
* (10) ./postSensDensity.py

# GraspPointFinder
Based on the "Force Closure" idea, the algorithm intends to find the optimum grasp point of a virtual 3D object.


## How it works
1. Import 3D object (stl or obj file)
2. Find possible grasp points (center of each mesh) and their wrench values in the matrix form
3. Make GWS (grasp wrench space)
4. Find the best grasp (grasp of the largest grasp quality)


## Brief explanation of each code
1. "gws_to_cone.py" : From 3D object, makes GWS data, and 3D object attached with friction cones of all possible grasp points
2. "searchall.py"   : Finds the best grasp among all possible grasp candidates - Most accurate, yet takes too much time
3. "gpsearch.py"    : Finds the locally best grasp through regressively updating the grasp  - Quicker, yet less accurate (80% accurray on average)

## How to execute
### If you want to find the accurate data taking a lot of time
1. gws_to_cone.py
2. searchall.py

### If you want to find the less accurate data within short period of time
1. gws_to_cone.py
2. gpsearch.py

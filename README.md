# GraspPointFinder
Based on the "Force Closure" idea, the algorithm intends to find the optimum grasp point based on virtual 3D object.


### How it works
1. Import 3D object (stl or obj file)
2. Find possible grasp points (center of each mesh) and its wrench
3. Make GWS (grasp wrench space)
4. Find the best grasp (grasp with largest grasp quality)


### How to use
1. "gws_to_cone.py" : From 3D object, makes GWS data, and 3D object attached with friction cones of all possible grasp points
2. "searchall.py"   : Finds the best grasp among all possible grasp candidates - Most accurate, yet takes too much time
3. "gpsearch.py"    : Finds the locally best grasp through regressively updating the grasp  - Quicker, yet less accurate (80% accurrate on average)

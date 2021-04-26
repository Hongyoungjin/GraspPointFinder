import numpy as np  
import math
import pandas as pd
import trimesh 
import sympy
import scipy.linalg as sp
import timeit
import random
from scipy.spatial import Delaunay
import openpyxl

def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0
start_time = timeit.default_timer()
stl_name = 'complex_2_high'

trial = 10
print(trial)

'''
Searching for the best grasp by finding the initial antipodal grasp pair

'''
# Define the tolerance of degree of two normal vectors
degree = 10
degree = degree*math.pi/180 # change in degrees

mesh = trimesh.load('./basic/' + stl_name+'.stl')
mesh.vertices -= mesh.center_mass

ver = mesh.vertices
tri = mesh.triangles

# Import the wrench 
W_1 = np.load('./gws/config/W_1.npy')
W_2 = np.load('./gws/config/W_2.npy')
W_3 = np.load('./gws/config/W_3.npy')
W_4 = np.load('./gws/config/W_4.npy')

# Import the Torsion
T   = np.load('./gws/config/T.npy')

# Store grasp quality and corresponding grasp pair indices into one array
GraspQuality = np.zeros([1,2]) # [ index, grasp_quality ]
Coefficients = np.zeros([1,1,6,6])
Normals      = np.zeros([1,6,6])
Offsets      = np.zeros([1,6])

# Find the grasp point with minimal torque to regard it as the initial grasp point




InPoint = np.array(random.sample(list(range(tri.shape[0])),5))


GPmax = 0
iteration = 0

Performance_iteration = []
Performance_quality = []
Performance_max = []

while 1:
    # Find the initial grasp point pair
    precise = 0
    convexsuccess  = 0
    for init in InPoint:
        index = []
    # 조건 1: 두 벡터의 내적이 -cos 10보다 작은가
        N = mesh.face_normals
        for i in range(N.shape[0]):
            if init != i:
                if np.dot(N[init],N[i]) < - math.cos(degree):    # first condition
                    index.append(i)

        # 조건 2: 두 파지점 사이의 거리가 그리퍼의 너비 안에 있는가 (80mm)
        # 조건 3: 두 파지점이 10도 이내로 틀어져 있는가

        index_2 = []

        # Find all possible secondary grasp point correponding to the initial one
        for index_1 in index:

            point1 = mesh.triangles_center[init]
            point2 = mesh.triangles_center[index_1]

            vector1 = N[init]  # 메쉬 바깥쪽을 향하는 벡터임을 명심하자
            vector2 = N[index_1]

            displace = point1-point2                             # 이어지는 벡터는 point1 벡터를 향해 있다.
            distance = math.sqrt(sum(displace*displace))

            unit_displace = displace/distance

            if distance < 80:                               # 조건 2 필터링
                
                if np.dot(vector1,unit_displace) >  math.cos(degree) and np.dot(vector2,unit_displace) < - math.cos(degree): # 조건 3 필터링
                    index_2.append(index_1)

        # Find the local optimum grasp of given grasp

        index_2 = np.array(index_2)
        
        
        for point2 in index_2 :

        # Import the wrench and torsion data from the first point of the pair

            w11 = W_1[init].reshape(1,6)[0]
            w12 = W_2[init].reshape(1,6)[0]
            w13 = W_3[init].reshape(1,6)[0]
            w14 = W_4[init].reshape(1,6)[0]

            w15 =  np.append(np.zeros([1,3]),T [init][np.newaxis],axis=1).reshape(1,6)[0]

            # Import the wrench and torsion data from the first point of the pair

            w21 = W_1[point2].reshape(1,6)[0]
            w22 = W_2[point2].reshape(1,6)[0]
            w23 = W_3[point2].reshape(1,6)[0]
            w24 = W_4[point2].reshape(1,6)[0]

            w25 =  np.append(np.zeros([1,3]),T [point2][np.newaxis],axis=1).reshape(1,6)[0]

            # Make 6-dimension convex hull using 10 6D vectors just imported

            delunay = np.zeros([1,6])
            b = 1
            
            while b < 6: 
                exec(f'delunay = np.append(delunay,w1{b}[np.newaxis],axis = 0)')
                exec(f'delunay = np.append(delunay,w2{b}[np.newaxis],axis = 0)')
                b+=1
            delunay = delunay[1:]        

            if in_hull(np.zeros([1,6]),delunay) == False:
                print(False)
                continue

            hull = Delaunay(delunay)
                    

            ver_index = hull.convex_hull # vertex indices of each hyperplane of the convex hull
            point     = hull.points     # input vertices used to make the convex hull
        
            
            Coeff = np.zeros([1,6,6])           # vertices of each hyperplane of the convex hull merged into an array
            PlaneInfo   = np.zeros([1,7])

            # Finding the 6D vertices consisting the hyperplane

            for i in range(ver_index.shape[0]):
                coeff = np.zeros([1,6])
                for j in range(ver_index.shape[1]):
                    index = ver_index[i,j]
                    coeff = np.append(coeff,point[index][np.newaxis],axis=0)

                coeff = coeff[1:] #- point[ver_index[i,0]]         
                coeff = coeff[np.newaxis]
                Coeff = np.append(Coeff,coeff,axis=0)
            Coeff = Coeff[1:] 

            # Calculating the normal vector of each hyperplane

            
            for PlaneNum in range(Coeff.shape[0]):
                        
                mat = np.concatenate((Coeff[PlaneNum],np.ones([6,1])),axis = 1)    # make X into the homogeneous form
                mat = sp.null_space(mat).T[0]
                PlaneInfo = np.append(PlaneInfo,mat[np.newaxis],axis = 0)

            PlaneInfo = PlaneInfo[1:]
            Distance = np.array([1,1])
            
            for k in range(ver_index.shape[0]):

                w = PlaneInfo[k,:-1]
                b = PlaneInfo[k, -1]

                w_mag = math.sqrt(np.sum(w*w))

                distance = abs(b)/w_mag
                Distance = np.append(Distance,np.array([distance]),axis=0)

            Distance = Distance[1:]
            grasp_quality = np.min(Distance)
            point2 = int(point2)
            index_quality = np.array([point2,grasp_quality])[np.newaxis]
            iteration+=1
            Performance_iteration.append(iteration)
            Performance_quality.append(grasp_quality)
            Performance_max.append(GPmax)
            GraspQuality = np.append(GraspQuality,index_quality,axis=0)

        GraspQuality = GraspQuality[1:]
            
        
        gp  = pd.DataFrame({
            'index'   : GraspQuality[:,0].tolist(),
            'quality' : GraspQuality[:,1].tolist()
        })
        gp = gp.sort_values(by = 'quality', ascending = False)
        if gp.shape[0] == 0:
            continue

        if gp.iloc[0,1] > GPmax:
            GPmax = gp.iloc[0,1]
            PointCandidate = np.array(gp.iloc[:,0],dtype = np.int)
            InitCandidate = init
            precise+=1
        
        
        
    if precise == 0:
        break
    grasp_pair = np.array([InitCandidate,PointCandidate[0]])
    InPoint = PointCandidate[0][np.newaxis]
    
    # ONE epoch is over
terminate_time = timeit.default_timer()


print('Iteration:',Performance_iteration[-1])

Graph = pd.DataFrame({
    'iteration' : Performance_iteration,
    'quality': Performance_quality,
    'Maximun Quality': Performance_max
})
Graph.to_excel('./SearchResults/Graph/' + stl_name + '_Graph_' + str(trial)+".xlsx")

Performance = pd.DataFrame({
    'Running Time' : terminate_time-start_time,
    'Grasp Quality': GPmax,
    'Grasp Pair_1'   : grasp_pair[0][np.newaxis],
    'Grasp Pair_2'   : grasp_pair[1][np.newaxis]
    
})

Performance.to_excel('./SearchResults/Performance/' + stl_name + '_Performance_' + str(trial)+".xlsx")
print('Time spent: %f'%(terminate_time-start_time))
    




    
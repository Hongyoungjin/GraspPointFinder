from itertools import combinations
import numpy as np 
from stl import mesh
import h5py
import math
import trimesh 
import sympy
import pandas as pd
import scipy.linalg as sp
import timeit
from scipy.spatial import Delaunay

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



stl_name = 'Link_medium'




start_time = timeit.default_timer()



# Import the wrench 
W_1 = np.load('./gws/config/W_1.npy')
W_2 = np.load('./gws/config/W_2.npy')
W_3 = np.load('./gws/config/W_3.npy')
W_4 = np.load('./gws/config/W_4.npy')

# Import the Torsion
T   = np.load('./gws/config/T.npy')


mesh = trimesh.load('./basic/' + stl_name+'.stl')
mesh.vertices -= mesh.center_mass

ver = mesh.vertices
tri = mesh.triangles


mylist = list(range(mesh.face_normals.shape[0]))
combi = combinations(mylist,2)



# Define the tolerance of degree of two normal vectors
degree = 10
degree = degree*math.pi/180 # change in degrees

# Define the index of grasp points that satisfies each condition
index = []

# 조건 1: 두 벡터의 내적이 -cos 10보다 작은가
N = mesh.face_normals
'''
for i in range(N.shape[0]):
    for j in range(N.shape[0]):
        if i != j:                                        # ignore same index
            if np.dot(N[i],N[j]) < - math.cos(degree):    # first condition
                if [j,i] not in index:                  # remove redundant index couple
                    index.append([i,j])
'''
index = np.array(list(combi)        )


# 조건 2: 두 파지점 사이의 거리가 그리퍼의 너비 안에 있는가 (80mm)
# 조건 3: 두 파지점이 10도 이내로 틀어져 있는가

index_2 = []
for i in range(index.shape[0]):

    point1 = mesh.triangles_center[index[i,0]]
    point2 = mesh.triangles_center[index[i,1]]

    vector1 = N[index[i,0]]  # 메쉬 바깥쪽을 향하는 벡터임을 명심하자
    vector2 = N[index[i,1]]

    displace = point1-point2                             # 이어지는 벡터는 point1 벡터를 향해 있다.
    distance = math.sqrt(sum(displace*displace))

    unit_displace = displace/distance


    if distance < 80:                               # 조건 2 필터링
        
        if np.dot(vector1,unit_displace) >  math.cos(degree) and np.dot(vector2,unit_displace) < - math.cos(degree): # 조건 3 필터링
            index_2.append(list(index[i]))

Index = np.array(index_2) # 조건 2,3 으로 filtering한 후의 index
print("Index Evaluation Over")
# Store grasp quality and corresponding grasp pair indices into one array
GraspQuality = np.zeros([1,2]) # [ index, grasp_quality ]
Coefficients = np.zeros([1,1,6,6])
Normals      = np.zeros([1,6,6])
Offsets      = np.zeros([1,6])

for index_pair in range(Index.shape[0]):

    index = Index[index_pair] # ex. [3, 106]

    # Import the wrench and torsion data from the first point of the pair

    w11 = W_1[index[0]].reshape(1,6)[0]
    w12 = W_2[index[0]].reshape(1,6)[0]
    w13 = W_3[index[0]].reshape(1,6)[0]
    w14 = W_4[index[0]].reshape(1,6)[0]

    w15 =  np.append(np.zeros([1,3]),T [index[0]][np.newaxis],axis=1).reshape(1,6)[0]

    # Import the wrench and torsion data from the first point of the pair

    w21 = W_1[index[1]].reshape(1,6)[0]
    w22 = W_2[index[1]].reshape(1,6)[0]
    w23 = W_3[index[1]].reshape(1,6)[0]
    w24 = W_4[index[1]].reshape(1,6)[0]

    w25 =  np.append(np.zeros([1,3]),T [index[1]][np.newaxis],axis=1).reshape(1,6)[0]

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
    point     = hull.points   # vertices used to make the convex hull


    # Define the center of Convex hull. 
    # If the center and origin are both located in the same side of all hyperplanes, the grasp is deemed robust.
    center = np.zeros([1,6])
    num = 0
    for i in range(10):
        if i in ver_index:

            center += point[i]
            num += 1
    center = center/num                  # center of each hyperplane


    Coeff = np.zeros([1,6,6])           # vertices of each hyperplane of the convex hull merged into an array
    PlaneInfo   = np.zeros([1,7])


    for i in range(ver_index.shape[0]):
        coeff = np.zeros([1,6])
        for j in range(ver_index.shape[1]):
            index = ver_index[i,j]
            coeff = np.append(coeff,point[index][np.newaxis],axis=0)
        coeff = coeff[1:]   
        coeff = coeff[np.newaxis]
        Coeff = np.append(Coeff,coeff,axis=0)
    Coeff = Coeff[1:] 

    num = 0
    for PlaneNum in range(Coeff.shape[0]):
                
        mat = np.concatenate((Coeff[PlaneNum],np.ones([6,1])),axis = 1)    # make X into the homogeneous form
        mat = sp.null_space(mat).T[0]
        
        

        # Origin and the center of convex hull should be on the same side divided by the hyperplane 
        origin_ans = mat[-1]                   # - offset
        center_ans = np.sum(mat[:-1]*center) + mat[-1]        

        if origin_ans * center_ans <=0 :        # 원점이 초평변 밖에 있다. 
            break
        num +=1
        PlaneInfo = np.append(PlaneInfo,mat[np.newaxis],axis = 0)

    PlaneInfo = PlaneInfo[1:]
    
    
    Distance = np.array([1,1])

    # 초평면과의 거리를 구해서 grasp quality 구하기
    if PlaneInfo.shape[0] == ver_index.shape[0]:
        for k in range(ver_index.shape[0]):

            w = PlaneInfo[k,:-1]
            b = PlaneInfo[k, -1]

            w_mag = math.sqrt(np.sum(w*w))

            distance = abs(b)/w_mag
            Distance = np.append(Distance,np.array([distance]),axis=0)

        Distance = Distance[1:]
        grasp_quality = np.min(Distance)
        
        index_pair = int(index_pair)
        index_quality = np.array([index_pair,grasp_quality])[np.newaxis]
        GraspQuality = np.append(GraspQuality,index_quality,axis=0)

GraspQuality = GraspQuality[1:]

# 상위 Grasp Point 구하기

gp  = pd.DataFrame({
    'index'   : GraspQuality[:,0].tolist(),
    'quality' : GraspQuality[:,1].tolist()
})

gp = gp.sort_values(by = 'quality', ascending = False)
gp = np.array(gp)
gquality = gp[0,1]
terminate_time = timeit.default_timer()

final_pair = Index[int(gp[0,0])]


terminate_time = timeit.default_timer()
print(gquality)
print('Time spent: %f'%(terminate_time-start_time))


Performance = pd.DataFrame({
    'Running Time' : terminate_time-start_time,
    'Grasp Quality': gquality,
    'Grasp Pair_1'   : final_pair[0][np.newaxis],
    'Grasp Pair_2'   : final_pair[1][np.newaxis]
    
})

Performance.to_excel('./AllResults/' + stl_name + '_Performance_.xlsx')

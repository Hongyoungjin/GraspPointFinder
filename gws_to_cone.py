import numpy as np  
import h5py
import math
import pandas as pd
import trimesh 
from stl import mesh


'''

Input : stl file of object of interest

Output:
R          : Transformation matrix sets of all triangles to the world frame
W, W_1, W_2, W_3 , W_4 : Wrench dataset of all triangles ini world frame

'''
def unit(x):
    div = math.sqrt(sum(x*x))
    y = x/div
    return y

stl_name = 'complex_2_high'

# Load the mesh to remove inner mesh grasp candidates
mesh_tri = trimesh.load('./basic/' + stl_name+'.stl')
mesh_tri.vertices -= mesh_tri.center_mass

ver = mesh_tri.vertices
tri = mesh_tri.triangles

print(tri.shape[0])


# Each value of triangle array in 'tri' refers to the index of vertex in 'ver'.

#config

K   = -1       # Normal vector of each triangle is basically headed outside (tri.face_normals()) 


m   = 0.5      # static friction coefficient
m_n = 0.5      # torsional friction coefficient


W_1 = np.zeros([1,2,3])
W_2 = np.zeros([1,2,3])
W_3 = np.zeros([1,2,3])
W_4 = np.zeros([1,2,3])
T   = np.zeros([1,3])

rho = 0
for i in range(ver.shape[0]):
    dist = ver[i]
    rho += sum(dist*dist) 
rho = math.sqrt(rho / ver.shape[0])

for i in range(tri.shape[0]): 

    e1 = tri[i][0]
    e2 = tri[i][1]
    e3 = tri[i][2] # defining the coordinates of each vertex in i th triangle
    e4 = mesh_tri.triangles_center[i]  # center of the mesh, which is the origin of object frame

    # define the object frame of ith triangle
    x = unit(e3-e2)
    z = mesh_tri.face_normals[i]
    y = unit(np.cross(z,x))

    # define the homogeneous transformation matrix from object frame to world frame
    x = x.reshape(3,1)
    y = y.reshape(3,1)
    z = z.reshape(3,1)
    e4 = e4.reshape(3,1)
    r = np.hstack((x,y,z))  # concatenate rotational frame and translational frame
    
    rr = np. linalg. inv(r)  # inverse transformation function
    # define the normal force onto the center of i th triangle
    
    # find the forces of frictional cone (approximated to rectantular cone)
    f_n = K / math.sqrt(1+m**2)            # normal     force in vector
    f_t = K * m                            # tangential force in scalar

    f   = np.array([  0 ,  0 ,  K ]).reshape(1,3)  # object frame. visualization에도 사용
    f_1 = np.array([f_t ,0   , f_n]).reshape(1,3)
    f_2 = np.array([0   ,f_t , f_n]).reshape(1,3)
    f_3 = np.array([-f_t,0   , f_n]).reshape(1,3)
    f_4 = np.array([0   ,-f_t, f_n]).reshape(1,3)
    
    # change the forces to the world frame  
    f_1 = f_1.T
    f_2 = f_2.T
    f_3 = f_3.T
    f_4 = f_4.T
    
    f_1 = np.matmul(r,f_1).T
    f_2 = np.matmul(r,f_2).T
    f_3 = np.matmul(r,f_3).T
    f_4 = np.matmul(r,f_4).T

    # find the torque
    torque_1 = np.cross(e4.reshape(1,3),f_1)
    torque_2 = np.cross(e4.reshape(1,3),f_2)
    torque_3 = np.cross(e4.reshape(1,3),f_3)
    torque_4 = np.cross(e4.reshape(1,3),f_4)

    # define the torsion of soft contact model - visualization에 사용되지 않는 변수들
    torsion   = f_n * m_n * mesh_tri.face_normals[i][np.newaxis] # world frame

                   
    # conclude the wrench
    w_1 = np.append(f_1  ,torque_1 /rho,  axis = 0)
    w_2 = np.append(f_2  ,torque_2 /rho,  axis = 0)
    w_3 = np.append(f_3  ,torque_3 /rho,  axis = 0)
    w_4 = np.append(f_4  ,torque_4 /rho,  axis = 0)
    
    
    # Store the tf matrix and wrench data into the i th index
    
    
    w_1 = w_1.reshape(1,2,3)
    w_2 = w_2.reshape(1,2,3)
    w_3 = w_3.reshape(1,2,3)
    w_4 = w_4.reshape(1,2,3)

    
    W_1 = np.append(W_1,w_1, axis = 0)
    W_2 = np.append(W_2,w_2, axis = 0)
    W_3 = np.append(W_3,w_3, axis = 0)
    W_4 = np.append(W_4,w_4, axis = 0)

    
    T   = np.append(T,torsion, axis = 0)
    
W_1 = W_1[1:]
W_2 = W_2[1:]
W_3 = W_3[1:]
W_4 = W_4[1:]
T   =  T [1:]


for x in ['W_1','W_2','W_3','W_4','T']:
    x_mat = eval(x)
    path = './gws/config/' + x
    np.save (path, x_mat)




# Returns the stl file of frictional cone of given object & grasp point

# Import wrench data
W_1 = np.load(str('./gws/config/W_1.npy'))
W_2 = np.load(str('./gws/config/W_2.npy'))
W_3 = np.load(str('./gws/config/W_3.npy'))
W_4 = np.load(str('./gws/config/W_4.npy'))

balance = -5


for i in range(W_1.shape[0]):
    e4 = mesh_tri.triangles_center[i][np.newaxis]

    f_1 = W_1[i,0][np.newaxis]* balance
    f_2 = W_2[i,0][np.newaxis]* balance
    f_3 = W_3[i,0][np.newaxis]* balance
    f_4 = W_4[i,0][np.newaxis]* balance


    f_1 = e4 + f_1
    f_2 = e4 + f_2
    f_3 = e4 + f_3
    f_4 = e4 + f_4

    
       
    # Define the vertices composing the frictional cone
    vertices = np.concatenate((e4, f_1, f_2, f_3, f_4),axis=0)
    
    # Define the triangles composing the frictional cone
    faces = np.array([\
        [0,1,2],
        [0,2,3],
        [0,3,4],
        [0,4,1],
        [1,4,3],
        [1,3,2],])

    # Create the mesh
    cone = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))

    for k, f in enumerate(faces):
        for j in range(3):
            cone.vectors[k][j] = vertices[f[j]]

    # Write the mesh to file 
    path = './gws/FrictionCone/' + stl_name+'_' + str(i)+ '.stl'
    cone.save(path)  

    
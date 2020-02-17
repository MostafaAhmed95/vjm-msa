import numpy as np
#code for the msa

#assumptions
#we rotate around z_axis

#variables to be set at the first of the code
E=1
L=1
G=1
S=1
J=1
I_x=1
I_y=1
I_z=1

q=1
w_list=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] #19


#function for the general K
K=np.zeros((6,6))

K[0,0]=(E*S)/L
K[1,1]=(12*E*I_z)/(L**3)
K[2,2]=(12*E*I_y)/(L**3)
K[3,3]=(G*J)/(L)
K[4,4]=(4*E*I_y)/(L)
K[5,5]=(4*E*I_z)/(L)
K[1,5]=(-6*E*I_z)/(L**2)
K[2,4]=(6*E*I_y)/(L**2)
K[4,2]=(6*E*I_y)/(L**2)
K[5,5]=(-6*E*I_z)/(L**2)

#function for inverse RR robot

#take the angle give the rotation matrix
#or give the angle and axis of rotation yield the rot matrix

#write the general stiffnes model
#w_general=k_general*t_general
w_general=np.array(w_list)

#general stiffeness matrix
#np.zeros((19,19))

#A_rigid matrix
A_r = np.diag((1,1,1,1,0,1))
#the non singular version
A_r=np.delete(A_r,-2,0)

#non singular A-p
A_p = np.array([0,0,0,0,1,0])
print(A_p)
A_p = A_p.reshape(1,-1)
print(A_p.shape)

#apply boundry conditions


#our general stiffness after applying boundry will be 114*114
k_after_reduction=np.zeros((114,114))
#create it element by element or assign a row by row
row=0

def create_row(l_t):
    b=np.zeros((l_t[0][1].shape[0],114))
    for i in l_t:
        b[:, i[0]*6:(i[0]*6)+6]=i[1]
    return b

#function for stacking the row
def mat_stack(x):
    np.vstack((k_after_reduction,x))

#the boundry conditions

# at base connection
r_1  = [(0,A_r)]
r_2=[(6,A_r)]
r_3=[(12,A_r)]

# at the end effector
r_4=[(5,A_r),(18,-A_r)]
r_5=[(11,A_r),(18,-A_r)]
r_6=[(17,A_r),(18,-A_r)]

# 2 middle passive joints for every arm
r_7=[(1,A_r),(2,-A_r)]
r_8=[(3,A_r),(4,-A_r)]
r_9=[(7,A_r),(8,-A_r)]
r_10=[(9,A_r),(10,-A_r)]
r_11=[(13,A_r),(14,-A_r)]
r_12=[(15,A_r),(16,-A_r)]

#Wrench related equations
r_13=[(0,A_p),(1,A_p)] #w0
r_14=[(0,A_p),(1,A_p)] #w1

#append all my rows into a list
r=[r_1,r_2,r_3,r_4,r_5,r_6,r_7,r_8,r_9,r_10,r_11,r_12,r_13,r_14]
for i in r:
    mat_stack(create_row(i))

#split the stiffenes matrix to get my A,B,C,D
A=k_after_reduction[0:113,0:113]
B=k_after_reduction[0:113,113]
C=k_after_reduction[113,0:113]
D=k_after_reduction[113,113]

# calculate our Kc
K_c= D - C @ np.linalg.inv(A) @ B
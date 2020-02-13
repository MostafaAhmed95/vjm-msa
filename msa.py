import numpy as np
#code for the msa

#assumptions
#we rotate around z_axis

#variables to be set at the first of the code
q=1
w_list=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] #19


#function for the general K

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

#the boundry conditions

# at base connection
r_1=[(0,A_r)]
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
import math
import numpy as np

Tbase=np.identity((4,4))
Ttool=np.identity((4,4))

#function for rotation matrices
def Rx(q):
    np.array([
        [1,0,0,0],
        [0,np.cos(q),np.sin(q),0],
        [0,-np.sin(q),np.cos(q),0],
        [0,0,0,1]])
def Ry(q):
    np.array([
        [np.cos(q), 0, np.sin(q), 0],
        [0, 1, 0, 0],
        [0, -np.sin(q), np.cos(q), 0],
        [0,0,0,1]])
def Rz(q):
    np.array([
        [np.cos(q), np.sin(q), 0, 0],
        [-np.sin(q), np.cos(q), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]])
def Tx(x):
    np.array([
        [1, 0, 0, x],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]])
def Ty(y):
    np.array([
        [1, 0, 0, 0],
        [0, 1, 0, y],
        [0, 0, 1, 0],
        [0, 0, 0, 1]])
def Tz(z):
    np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, z],
        [0, 0, 0, 1]])

def local_inv(x,y):
    """calculate the angles given the end effector position"""
    l1,l2=1,1  #link lengths
    #elbow up configuration
    p = math.sqrt(x**2+y**2)
    d = (-l2 ** 2 + l1 ** 2 + p**2)/2*p
    n = math.sqrt(l1**2 - x**2)
    phi = math.atan(y/x)
    beta = math.atan(n/d)
    alpha = math.atan(n/(p-d))
    theta_1_1= phi+beta
    theta_1_2= -(alpha+beta)
    return theta_1_1,theta_1_2
def tot_inv(x,y,z):
    "calculate angles for the three arms"
    q1,q2=local_inv(x,-y) #first arm
    q3,q4=local_inv(3-z,-y) #second arm
    q5,q6=local_inv(3-x,z) #third arm
    q = [q1,q2,q3,q4,q5,q6]
def vjm_total(x,y,z):
    """calculate robot stiffness matrix in one configuration"""
    q=tot_inv(x,y,z)
    k = np.linalg.inv(vjm_1(Tbase, Ttool, q[0], q, t, L)+vjm_1(Tbase, Ttool, q[0], q, t, L))
    return k
def forward_kin(Tbase,Ttool,q0,q,t,L):
    T = Tbase.dot(Tx(-d / 2)).dot(Rz(q0(1))).dot(Rz(t(1))).dot(Tx(L)).dot(Tx(t(2))).dot(Ty(t(3))).dot(Tz(t(4))).dot(
    Rx(t(5))).dot(Ry(t(6))).dot(Rz(t(7))).dot(Rz(q(1))).dot(Tx(l)).dot(
    Tx(t(8))).dot(Ty(t(9))).dot(Tz(t(10))).dot(Rx(t(11))).dot(Ry(t(12))).dot(Rz(t(13))).dot(Ttool)
    T_d=T.copy()
    T[0:3,3]=[0,0,0]
    T=np.transpose(T)
    T_d=T_d*T

    J1 = np.array([
        T_d[0,3],T_d[1,3],T_d[2,3],T_d[2,1],T_d[0,2],T_d[1,0]
    ]).transpose()
    return J1

def vjm_1(Tbase, Ttool, q[0], q, t, L):

    k0 = 1e6 # Actuator stif
    # material and shape parameters
    E = 70e9 # Young's modulus
    G = 25.5e9 # shearmodulus
    d = 50e-3
    L=1

    # for cylinder
    S = math.pi * d ^ 2 / 4
    Iy = math.pi * d ^ 4 / 64
    Iz = math.pi * d ^ 4 / 64
    J=Iy+Iz

    #general stiffness matrix
    K_1=np.zeros((6,6))

    K_1[0,0]=(E*S)/L
    K_1[1,1]=(12*E*Iz)/(L**3)
    K_1[2,2]=(12*E*Iy)/(L**3)
    K_1[3,3]=(G*J)/(L)
    K_1[4,4]=(4*E*Iy)/(L)
    K_1[5,5]=(4*E*Iz)/(L)
    K_1[1,5]=(-6*E*Iz)/(L**2)
    K_1[2,4]=(6*E*Iy)/(L**2)
    K_1[4,2]=(6*E*Iy)/(L**2)
    K_1[5,5]=(-6*E*Iz)/(L**2)

    K_2=K_1.copy()

    x = np.hstack((k0,np.zeros((1,12))))
    y = np.hstack((np.zeros((6,1)),K_1))
    y = np.hstack((y,np.zeros(6,6)))
    z = np.hstack((np.zeros(6,1),np.zeros((6,6))))
    z = np.hstack((z,K_2))

    #K_theta
    K_t = np.vstack((x,y))
    K_t = np.vstack((K_t,z))

    #jacobian
    Jq = forward_kin(Tbase,Ttool,q0,q,t,L,l,d)
    Jt = forward_kin(Tbase,Ttool,q0,q,t,L,l,d)

    #analytical solution
    Kc0=np.linalg.inv(Jt*np.linalg.inv(Kt)*Jt')
    Kc = Kc0 - Kc0*Jq*np.linalg.inv(Jq'*Kc0*Jq)*Jq'*Kc0
#assume that we have the whole transformation


def get_def(f):
    """given the force calculate deflection at different configurations"""
    # limits of our robot workspace in the three spaces are from -200 to 200
    dt=[]
    for x in range(-200,200):
        for y in range(-200,200):
            for z in range(-200, 200):
                k = vjm_total(x,y,z)
                #list of deflections
                deflection = k * f
                #calculate magntiude of the deflection for plotting
                dt.append(np.sqrt(deflection[0]**2+deflection[1]**2+deflection[2]**2+))
    return dt



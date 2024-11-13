import numpy as np
import pandas as pd
from utils import Rod

MONOMER_LIST = ["C5"]
############################汎用関数###########################
def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A1,A2,A3):
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv(r'/home/ohno/Working/alkyl_chain/step1_z=1/monomer/{}.csv'.format(monomer_name))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)
    xyz_array = np.matmul(xyz_array,Rod(ey,A1).T)
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    
    if monomer_name in MONOMER_LIST:
        return np.concatenate([xyz_array,R_array],axis=1)
    
    else:
        raise RuntimeError('invalid monomer_name={}'.format(monomer_name))

def vdw_R_x(A1,A2,theta,monomer_name,vx,vy,vz):##2分子が接しないようにするにはtheta方向にどれだけずらせばいいのか
    monomer_1=get_monomer_xyzR(monomer_name,0.,0.,0.,A1,A2,theta)
    R_clps=0
    for x1,y1,z1,rad1 in monomer_1:
        for x2,y2,z2,rad2 in monomer_1:
            x2+=vx;y2+=vy;z2+=vz
            eR=np.array([1.0,0.0,0.0])
            R_12=np.array([x2-x1,y2-y1,z2-z1])
            R_12a=np.linalg.norm([0,y2-y1,z2-z1])
            if (rad1+rad2)**2-R_12a**2<0:
                continue
            else:
                R_clps=max(R_clps,x1-x2+np.sqrt((rad1+rad2)**2-R_12a**2))
    return R_clps

def vdw_R_y(A1,A2,theta,monomer_name,vx,vy,vz):##2分子が接しないようにするにはtheta方向にどれだけずらせばいいのか
    monomer_1=get_monomer_xyzR(monomer_name,0.,0.,0.,A1,A2,theta)
    R_clps=0
    for x1,y1,z1,rad1 in monomer_1:
        for x2,y2,z2,rad2 in monomer_1:
            x2+=vx;y2+=vy;z2+=vz
            eR=np.array([1.0,0.0,0.0])
            R_12=np.array([x2-x1,y2-y1,z2-z1])
            R_12a=np.linalg.norm([x2-x1,0,z2-z1])
            if (rad1+rad2)**2-R_12a**2<0:
                continue
            else:
                R_clps=max(R_clps,y1-y2+np.sqrt((rad1+rad2)**2-R_12a**2))
    return R_clps

monomer_name='C5'

para_list1=[]
for i in range(19):
    for j in range(26):
        theta=int(i*5);z1=np.round(j*0.1,1)
        a1=vdw_R_x(0,0,theta,monomer_name,0,0,z1)
        para_list1.append([theta,z1,np.round(a1,1),a1])

para_list=[]
for theta,z1,a1,a1_raw in para_list1:
    a2_list=[];z2_list=[];b2_list=[]
    for i in range(int(a1*10)+1):
        for j in range(51):
            a2=np.round(-np.round(a1/2,1)+0.1*i,1);z2=np.round(-2.5+0.1*j,1)
            b2_1=vdw_R_y(0,0,theta,monomer_name,a2,0,z2)
            b2_2=vdw_R_y(0,0,theta,monomer_name,a2+a1,0,z2+z1)
            b2_3=vdw_R_y(0,0,theta,monomer_name,a2-a1,0,z2-z1)
            b2=max(b2_1,b2_2,b2_3)
            para_list.append([theta,a1,z1,a2,np.round(b2,1),z2,a1_raw,b2])

df=pd.DataFrame(para_list,columns=['theta','a1','z1','a2','b2','z2','a1_raw','b2_raw'])
df.to_csv(f'/home/ohno/Working/alkyl_chain/step1_z=1/vdw_calc/{monomer_name}_vdw_result.csv',index=False)
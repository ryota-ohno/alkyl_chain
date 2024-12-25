import os
import numpy as np
import pandas as pd
import subprocess
from utils import Rod, R2atom

MONOMER_LIST = ["C5",'C6','C28']
############################汎用関数###########################
def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3):
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv('~/Working/alkyl_chain/step1_z=1_z/monomer/{}_re.csv'.format(monomer_name))
    atoms_array_xyzR=df_mono[['X','Y','Z','R']].values
    
    ex = np.array([1.,0.,0.]); ey = np.array([0.,1.,0.]); ez = np.array([0.,0.,1.])

    xyz_array = atoms_array_xyzR[:,:3]
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)
    xyz_array = xyz_array + T_vec
    R_array = atoms_array_xyzR[:,3].reshape((-1,1))
    
    if monomer_name in MONOMER_LIST:
        return np.concatenate([xyz_array,R_array],axis=1)
    else:
        raise RuntimeError('invalid monomer_name={}'.format(monomer_name))
        
def get_xyzR_lines(xyzR_array,file_description,machine_type):
    if machine_type==1:
        mp_num = 40
    elif machine_type==2:
        mp_num = 52
    lines = [     
        '%mem=20GB\n',
        f'%nproc={mp_num}\n',
        '#B3LYP/6-311G** EmpiricalDispersion=GD3BJ counterpoise=2\n',
        '\n',
        file_description+'\n',
        '\n',
        '0 1 0 1 0 1\n'
    ]
    mol_len = len(xyzR_array)//2
    atom_index = 0
    mol_index = 0
    for x,y,z,R in xyzR_array:
        atom = R2atom(R)
        mol_index = atom_index//mol_len + 1
        line = '{}(Fragment={}) {} {} {}\n'.format(atom,mol_index,x,y,z)     
        lines.append(line)
        atom_index += 1
    return lines

# 実行ファイル作成
def get_one_exe(file_name,machine_type):
    file_basename = os.path.splitext(file_name)[0]
    #mkdir
    if machine_type==1:
        gr_num = 1; mp_num = 40
    elif machine_type==2:
        gr_num = 2; mp_num = 52
    cc_list=[
        '#!/bin/sh \n',
        '#$ -S /bin/sh \n',
        '#$ -cwd \n',
        '#$ -V \n',
        '#$ -q gr{}.q \n'.format(gr_num),
        '#$ -pe OpenMP {} \n'.format(mp_num),
        '\n',
        'hostname \n',
        '\n',
        'export g16root=/home/g03 \n',
        'source $g16root/g16/bsd/g16.profile \n',
        '\n',
        'export GAUSS_SCRDIR=/scr/$JOB_ID \n',
        'mkdir /scr/$JOB_ID \n',
        '\n',
        'g16 < {}.inp > {}.log \n'.format(file_basename,file_basename),
        '\n',
        'rm -rf /scr/$JOB_ID \n',
        '\n',
        '\n',
        '#sleep 5 \n'
#          '#sleep 500 \n'
            ]

    return cc_list

######################################## 特化関数 ########################################

##################gaussview##################
def make_xyzfile(monomer_name,params_dict):
    a = params_dict.get('a',0.0); b = params_dict.get('b',0.0); z = params_dict.get('z',0.0)
    A2 = params_dict.get('A2',0.0); A3 = params_dict.get('theta',0.0)

    monomer_array_1 = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    monomer_array_2 = get_monomer_xyzR(monomer_name,a,b,z,A2,A3)##1,2がb方向
    
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    monomers_array_4 = np.concatenate([monomer_array_1,monomer_array_2],axis=0)
    
    for x,y,z,R in monomers_array_4:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        if key in ['a','b','z']:
            val = np.round(val,2)
        elif key in ['A1','A2','theta']:
            val = int(val)
        xyzfile_name += '_{}={}'.format(key,val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict,machine_type):
    a = params_dict.get('a',0.0); b = params_dict.get('b',0.0); z = params_dict.get('z',0.0)
    A2 = params_dict.get('A2',0.0); A3 = params_dict.get('theta',0.0)

    monomer_array_1 = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    monomer_array_2 = get_monomer_xyzR(monomer_name,a,b,z,A2,A3)##1,2がb方向
    
    dimer_array = np.concatenate([monomer_array_1,monomer_array_2])
    
    file_description = f'{monomer_name}_a={a}_b={b}_z={z}'
    line_list_dimer = get_xyzR_lines(dimer_array,file_description,machine_type)
    
    gij_xyz_lines = ['$ RunGauss\n'] + line_list_dimer + ['\n\n\n']
    
    
    file_name = get_file_name_from_dict(monomer_name,params_dict)
    os.makedirs(os.path.join(auto_dir,'gaussian'),exist_ok=True)
    gij_xyz_path = os.path.join(auto_dir,'gaussian',file_name)
    with open(gij_xyz_path,'w') as f:
        f.writelines(gij_xyz_lines)
    
    return file_name

def get_file_name_from_dict(monomer_name,params_dict):
    file_name = ''
    file_name += monomer_name
    for key,val in params_dict.items():
        if key in ['a','b','z']:
            val = val
        elif key in ['A1','A2','theta']:
            val = int(val)
        file_name += '_{}={}'.format(key,val)
    return file_name + '.inp'
    
def exec_gjf(auto_dir, monomer_name, params_dict, machine_type,isTest=True):
    inp_dir = os.path.join(auto_dir,'gaussian')
    xyz_dir = os.path.join(auto_dir,'gaussview')
    print(params_dict)
    
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict,machine_type)
    cc_list = get_one_exe(file_name,machine_type)
    sh_filename = os.path.splitext(file_name)[0]+'.r1'
    sh_path = os.path.join(inp_dir,sh_filename)
    with open(sh_path,'w') as f:
        f.writelines(cc_list)
    if not(isTest):
        subprocess.run(['qsub',sh_path])
    log_file_name = os.path.splitext(file_name)[0]+'.log'
    return log_file_name
    
############################################################################################
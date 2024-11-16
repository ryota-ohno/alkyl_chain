##tetracene層内計算 step2の精密化に使用　パラメータにR3,R4も加える←これはむしろmake 6分子
import os
os.environ['HOME'] ='/home/ohno'
import pandas as pd
import time
import sys
sys.path.append(os.path.join(os.environ['HOME'],'Working/interaction/'))
from utils import get_E
import argparse
import numpy as np
import random
from make_step2_auto_1__ import exec_gjf

def submit_process(args):
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    isTest= args.isTest
    isMain= args.isMain
    if isMain:
        return
    isEnd= args.isEnd
    if isEnd:
        return
    os.makedirs(os.path.join(auto_dir,'gaussian'),exist_ok=True)
    init_params_csv=os.path.join(auto_dir, 'step2_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params = df_init_params.iloc[0]
    params_dict = df_init_params.to_dict()
    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    log_file= exec_gjf(auto_dir, monomer_name, params_dict,isTest)
    time.sleep(2)
    print(log_file)


def main_process(args):
    isEnd= args.isEnd
    if isEnd:
        return
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)

    init_params_csv=os.path.join(auto_dir, 'step2_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params = df_init_params.iloc[0]
    params_dict = df_init_params.to_dict()
    file_base_name = get_file_base_name(monomer_name,params_dict)
    file_name_1 = file_base_name
    file_name_2 = file_base_name
    file_name_3 = file_base_name
    file_name_1 += '1.log'
    file_name_2 += '2.log'
    file_name_3 += '3.log'
    log_filepath_1 = os.path.join(*[auto_dir,'gaussian',file_name_1])
    log_filepath_2 = os.path.join(*[auto_dir,'gaussian',file_name_2])
    log_filepath_3 = os.path.join(*[auto_dir,'gaussian',file_name_3])
    

    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(log_filepath_1,log_filepath_2,log_filepath_3)##argsの中身を取る
        time.sleep(1)

def listen(log_filepath_1,log_filepath_2,log_filepath_3):##args自体を引数に取るか中身をばらして取るかの違い    
    E_list1=get_E(log_filepath_1)
    E_list2=get_E(log_filepath_2)
    E_list3=get_E(log_filepath_3)
    if len(E_list1)!=40 or len(E_list2)!=40 or len(E_list3)!=40:##計算する層状の分子数
        isOver =False
    else:
        isOver=True
    return isOver
        

def get_file_base_name(monomer_name,params_dict):
    a_ = params_dict['a']; b_ = params_dict['b']; theta = params_dict['theta']
    file_base_name = ''
    file_base_name += monomer_name
    file_base_name += '_step2_'
    file_base_name += 'a={}_b={}_theta={}_'.format(a_,b_,theta)
    return file_base_name

def end_process(args):
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    init_params_csv=os.path.join(auto_dir, 'step2_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params = df_init_params.iloc[0]
    params_dict = df_init_params.to_dict()
    file_base_name = get_file_base_name(monomer_name,params_dict)
    file_name_1 = file_base_name
    file_name_2 = file_base_name
    file_name_3 = file_base_name
    file_name_1 += '1.log'
    file_name_2 += '2.log'
    file_name_3 += '3.log'
    log_filepath_1 = os.path.join(*[auto_dir,'gaussian',file_name_1])
    log_filepath_2 = os.path.join(*[auto_dir,'gaussian',file_name_2])
    log_filepath_3 = os.path.join(*[auto_dir,'gaussian',file_name_3])
    z_list=[np.round(z,1) for z in np.linspace(np.round(-4.0,1),np.round(-0.1,1),int(np.round(np.round(3.9,1)/0.1))+1)]
    E_list1=get_E(log_filepath_1)##t
    E_list2=get_E(log_filepath_2)##p
    E_list3=get_E(log_filepath_3)##p
    para_list=[]
    for i in range(len(E_list1)):
        z=z_list[i];Et1=E_list1[i];Ep1=E_list2[i];Et4=E_list3[i]
        para_list.append([z,Et1,Ep1,Et4])
    df=pd.DataFrame(para_list,columns=['z','Et1','Ep1','Et4'])
    csv_path=os.path.join(auto_dir, 'step2.csv')
    df.to_csv(csv_path,index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    #parser.add_argument('--init',action='store_true')
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--isEnd',action='store_true')
    parser.add_argument('--isMain',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    
    print("----main process----")
    submit_process(args)##step2はここで実行
    main_process(args)##こっちは確認
    end_process(args)
    ##最後に更新とか極小点周りの話
    print("----finish process----")
    
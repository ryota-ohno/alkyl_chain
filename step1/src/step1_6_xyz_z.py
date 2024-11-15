##tetracene層内計算
import os
os.environ['HOME'] ='/home/ohno'
import pandas as pd
import time
import sys
sys.path.append(os.path.join(os.environ['HOME'],'Working/interaction/'))
from make_6_xyz_z import exec_gjf##計算した点のxyzfileを出す
from utils import get_E
import argparse
import numpy as np
import random

def main_process(args):
    auto_dir = args.auto_dir
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    auto_csv_path = os.path.join(auto_dir,'step1.csv')
    if not os.path.exists(auto_csv_path):        
        df_E = pd.DataFrame(columns = ['a','b','theta','z','E','E_p1','E_p2','E_t1','E_t4','machine_type','status','file_name'])##いじる
        df_E.to_csv(auto_csv_path,index=False)##step3を二段階でやる場合二段階目ではinitをやらないので念のためmainにも組み込んでおく

    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(args.auto_dir,args.monomer_name,args.num_nodes,args.num_m2,args.isTest)##argsの中身を取る
        time.sleep(1)

def listen(auto_dir,monomer_name,num_nodes,num_m2,isTest):##args自体を引数に取るか中身をばらして取るかの違い
    auto_csv = os.path.join(auto_dir,'step1.csv')
    df_E = pd.read_csv(auto_csv)
    df_queue = df_E.loc[df_E['status']=='InProgress',['machine_type','file_name']]
    machine_type_list = df_queue['machine_type'].values.tolist()
    len_queue = len(df_queue)
    maxnum_machine2 = num_m2#int(num_nodes/2) ##多分俺のために空けていてくださったので2 3にする
    
    for idx,row in zip(df_queue.index,df_queue.values):
        machine_type,file_name = row
        log_filepath = os.path.join(*[auto_dir,'gaussian',file_name])
        if not(os.path.exists(log_filepath)):#logファイルが生成される直前だとまずいので
            continue
        E_list=get_E(log_filepath)
        if len(E_list)!=4:##get Eの長さは計算した分子の数
            continue
        else:
            len_queue-=1;machine_type_list.remove(machine_type)
            Et1=float(E_list[0]);Et4=float(E_list[1]);Ep1=float(E_list[2]);Ep2=float(E_list[3])##8分子に向けてep1,ep2作成　ep1:b ep2:a
            E = 2*(Et1+Et4+Ep1+Ep2)##エネルギーの値も変える
            df_E.loc[idx, ['E_t1','E_t4','E_p1','E_p2','E','status']] = [Et1,Et4,Ep1,Ep2,E,'Done']
            df_E.to_csv(auto_csv,index=False)
            break#2つ同時に計算終わったりしたらまずいので一個で切る
    isAvailable = len_queue < num_nodes 
    machine2IsFull = machine_type_list.count(2) >= maxnum_machine2
    machine_type = 1 if machine2IsFull else 2
    if isAvailable:
        params_dict = get_params_dict(auto_dir,num_nodes)
        if len(params_dict)!=0:#終わりがまだ見えないなら
            alreadyCalculated = check_calc_status(auto_dir,params_dict)
            if not(alreadyCalculated):
                file_name = exec_gjf(auto_dir, monomer_name, {**params_dict,'cx':0,'cy':0,'cz':0,'A1':0.,'A2':0.}, machine_type,isInterlayer=False,isTest=isTest)##計算を実行並びにxyzファイルの出力
                df_newline = pd.Series({**params_dict,'E':0.,'E_p1':0.,'E_p2':0.,'E_t1':0.,'E_t4':0.,'machine_type':machine_type,'status':'InProgress','file_name':file_name})
                df_E=df_E.append(df_newline,ignore_index=True)
                df_E.to_csv(auto_csv,index=False)
                time.sleep(1)
    
    init_params_csv=os.path.join(auto_dir, 'step1_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params_done = filter_df(df_init_params,{'status':'Done'})
    isOver = True if len(df_init_params_done)==len(df_init_params) else False
    return isOver

def check_calc_status(auto_dir,params_dict):
    df_E= pd.read_csv(os.path.join(auto_dir,'step1.csv'))
    if len(df_E)==0:
        return False
    df_E_filtered = filter_df(df_E, params_dict)
    df_E_filtered = df_E_filtered.reset_index(drop=True)
    try:
        status = get_values_from_df(df_E_filtered,0,'status')
        return status=='Done'
    except KeyError:
        return False

def get_params_dict(auto_dir, num_nodes):
    """
    前提:
        step1_init_params.csvとstep1.csvがauto_dirの下にある
    """
    init_params_csv=os.path.join(auto_dir, 'step1_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_cur = pd.read_csv(os.path.join(auto_dir, 'step1.csv'))
    df_init_params_inprogress = df_init_params[df_init_params['status']=='InProgress']
    fixed_param_keys = ['theta','z']
    opt_param_keys = ['a','b']

    #最初の立ち上がり時
    if len(df_init_params_inprogress) < num_nodes:
        print(1)
        df_init_params_notyet = df_init_params[df_init_params['status']=='NotYet']
        for index in df_init_params_notyet.index:
            df_init_params = update_value_in_df(df_init_params,index,'status','InProgress')
            df_init_params.to_csv(init_params_csv,index=False)
            params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
            return params_dict
    for index in df_init_params.index:
        df_init_params = pd.read_csv(init_params_csv)
        print('index')
        print(index)
        init_params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
        fixed_params_dict = df_init_params.loc[index,fixed_param_keys].to_dict()
        isDone, opt_params_dict = get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict)
        if isDone:
            # df_init_paramsのstatusをupdate
            df_init_params = update_value_in_df(df_init_params,index,'status','Done')
            if np.max(df_init_params.index) < index+1:
                status = 'Done'
            else:
                status = get_values_from_df(df_init_params,index+1,'status')
            df_init_params.to_csv(init_params_csv,index=False)
            
            if status=='NotYet':                
                opt_params_dict = get_values_from_df(df_init_params,index+1,opt_param_keys)
                df_init_params = update_value_in_df(df_init_params,index+1,'status','InProgress')
                df_init_params.to_csv(init_params_csv,index=False)
                return {**fixed_params_dict,**opt_params_dict}
            else:
                continue

        else:
            df_inprogress = filter_df(df_cur, {**fixed_params_dict,**opt_params_dict,'status':'InProgress'})
            print(df_inprogress)
            if len(df_inprogress)>=1:
                print('continue')
                continue
            return {**fixed_params_dict,**opt_params_dict}
    return {}
        
def get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict):
    df_val = filter_df(df_cur, fixed_params_dict)
    a_init_prev = init_params_dict['a']; b_init_prev = init_params_dict['b']; theta = init_params_dict['theta']; z = init_params_dict['z']
    while True:
        E_list=[];ab_list=[]
        for a in [a_init_prev-0.1,a_init_prev,a_init_prev+0.1]:
            for b in [b_init_prev-0.1,b_init_prev,b_init_prev+0.1]:
                a = np.round(a,1);b = np.round(b,1)
                df_val_ab = df_val[(df_val['a']==a)&(df_val['b']==b)&(df_val['theta']==theta)&(df_val['z']==z)&(df_val['status']=='Done')]
                if len(df_val_ab)==0:
                    return False,{'a':a,'b':b}
                ab_list.append([a,b]);E_list.append(df_val_ab['E'].values[0])
        a_init,b_init = ab_list[np.argmin(np.array(E_list))]
        if a_init==a_init_prev and b_init==b_init_prev:
            return True,{'a':a_init,'b':b_init}
        else:
            a_init_prev=a_init;b_init_prev=b_init

def get_values_from_df(df,index,key):
    return df.loc[index,key]

def update_value_in_df(df,index,key,value):
    df.loc[index,key]=value
    return df

def filter_df(df, dict_filter):
    for k, v in dict_filter.items():
        if type(v)==str:
            df=df[df[k]==v]
        else:
            df=df[df[k]==v]
    df_filtered=df
    return df_filtered

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--init',action='store_true')
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    parser.add_argument('--num-m2',type=int,help='num machine2')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    print("----main process----")
    main_process(args)
    print("----finish process----")
    
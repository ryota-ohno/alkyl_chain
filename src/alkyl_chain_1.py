##tetracene層内計算
import os
os.environ['HOME'] ='/home/ohno'
import pandas as pd
import time
import sys
from tqdm import tqdm
sys.path.append(os.path.join(os.environ['HOME'],'Working/interaction/'))
from make_ac_1 import exec_gjf##計算した点のxyzfileを出す
from utils import get_E
from utils import get_E_mono_1
import argparse
import numpy as np
from scipy import signal
import scipy.spatial.distance as distance
import random

def init_process(args):
    # 数理モデル的に自然な定義の元のparams initリスト: not yet
    # 結晶学的に自然なパラメータへ変換: not yet
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    order = 5
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)

    def get_init_para_csv(auto_dir):
        init_params_csv = os.path.join(auto_dir, 'step1_init_params.csv')
        
        init_para_list = []
        A1 = 0; A2 = 0
        for i in range(11):
            for j in range(11):
                a=6.8+0.1*i
                b=5.4+0.1*j
                a=np.round(a,1)
                b=np.round(b,1)
                init_para_list.append([a,b,27,'NotYet'])
            
        df_init_params = pd.DataFrame(np.array(init_para_list),columns = ['a','b','theta','status'])
        df_init_params.to_csv(init_params_csv,index=False)
    
    get_init_para_csv(auto_dir)
    
    auto_csv_path = os.path.join(auto_dir,'step1.csv')
    if not os.path.exists(auto_csv_path):        
        df_E_init = pd.DataFrame(columns = ['a','b','theta','E','E_p','E_t','machine_type','status','file_name'])##隣接6分子
    else:
        df_E_init = pd.read_csv(auto_csv_path)
        df_E_init = df_E_init[df_E_init['status']!='InProgress']
    df_E_init.to_csv(auto_csv_path,index=False)

    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    df_init['status']='NotYet'
    df_init.to_csv(os.path.join(auto_dir,'step1_init_params.csv'),index=False)

def main_process(args):
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    auto_csv_path = os.path.join(auto_dir,'step1.csv')
    if not os.path.exists(auto_csv_path):        
        df_E = pd.DataFrame(columns = ['phi1','a','b','theta','A1','A2','E','E_p1','E_t1','E_t2','E_m','machine_type','status','file_name'])##いじる
        df_E.to_csv(auto_csv_path,index=False)##step3を二段階でやる場合二段階目ではinitをやらないので念のためmainにも組み込んでおく

    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(args.auto_dir,args.monomer_name,args.num_nodes,args.isTest)##argsの中身を取る
        time.sleep(1)

def listen(auto_dir,monomer_name,num_nodes,isTest):##args自体を引数に取るか中身をばらして取るかの違い
    auto_csv = os.path.join(auto_dir,'step1.csv')
    df_E = pd.read_csv(auto_csv)
    df_queue = df_E.loc[df_E['status']=='InProgress',['machine_type','file_name']]
    machine_type_list = df_queue['machine_type'].values.tolist()
    len_queue = len(df_queue)
    maxnum_machine2 = 1#int(num_nodes/2) ##3にする
    
    for idx,row in zip(df_queue.index,df_queue.values):
        machine_type,file_name = row
        log_filepath = os.path.join(*[auto_dir,'gaussian',file_name])
        if not(os.path.exists(log_filepath)):#logファイルが生成される直前だとまずいので
            continue
        E_list=get_E(log_filepath)
        if len(E_list)!=3:##get Eの長さは計算した分子の数
            continue
        else:
            len_queue-=1;machine_type_list.remove(machine_type)
            Ep1=float(E_list[0]);Et1=float(E_list[1]);Et2=float(E_list[2])##8分子に向けてep1,ep2作成　ep1:b ep2:a
            E_m = get_E_mono_1(log_filepath)
            E = Ep1+Et1+Et2+E_m##エネルギーの値も変える
            df_E.loc[idx, ['E_t1','E_t2','E_p1','E_m','E','status']] = [Et1,Et2,Ep1,E_m,E,'Done']
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
                file_name = exec_gjf(auto_dir, monomer_name, {**params_dict,'cx':0,'cy':0,'cz':0,'R3':0.,'R4':0.}, machine_type,isInterlayer=False,isTest=isTest)##計算を実行並びにxyzファイルの出力
                df_newline = pd.Series({**params_dict,'E':0.,'E_m':0.,'E_p1':0.,'E_t1':0.,'E_t2':0.,'machine_type':machine_type,'status':'InProgress','file_name':file_name})
                df_E=df_E.append(df_newline,ignore_index=True)
                df_E.to_csv(auto_csv,index=False)
    
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
    opt_param_keys = ['phi1']
    fixed_param_keys = ['a','b','theta','A1','A2']

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
    phi1_init_prev = init_params_dict['phi1']
    a = init_params_dict['a']; b = init_params_dict['b']; theta = init_params_dict['theta']
    while True:
        E_list=[];phi1_list=[]
        for phi1 in [phi1_init_prev-1,phi1_init_prev,phi1_init_prev+1]:
            phi1 = int(phi1)
            df_val_phi1 = df_val[(df_val['a']==a)&(df_val['b']==b)&(df_val['theta']==theta)&(df_val['phi1']==phi1)&(df_val['status']=='Done')]
            if len(df_val_phi1)==0:
                return False,{'phi1':phi1}
            phi1_list.append(phi1);E_list.append(df_val_phi1['E'].values[0])
        phi1_init = phi1_list[np.argmin(np.array(E_list))]
        if phi1_init == phi1_init_prev:
            return True,{'phi1':phi1_init}
        else:
            phi1_init_prev = phi1_init

def get_values_from_df(df,index,key):
    return df.loc[index,key]

def update_value_in_df(df,index,key,value):
    df.loc[index,key]=value
    return df

def filter_df(df, dict_filter):
    query = []
    for k, v in dict_filter.items():
        if type(v)==str:
            query.append('{} == "{}"'.format(k,v))
        else:
            query.append('{} == {}'.format(k,v))
    df['A1']=df['A1'].astype(float)
    df['A2']=df['A2'].astype(float)
    df_filtered = df.query(' and '.join(query))
    return df_filtered

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--init',action='store_true')
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    if args.init:
        print("----initial process----")
        init_process(args)
    
    print("----main process----")
    main_process(args)
    print("----finish process----")
    
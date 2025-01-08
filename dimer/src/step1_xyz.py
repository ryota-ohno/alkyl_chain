##tetracene層内計算
import os
os.environ['HOME'] ='/home/ohno'
import pandas as pd
import time
import sys
sys.path.append(os.path.join(os.environ['HOME'],'Working/interaction/'))
from make_xyz import exec_gjf##計算した点のxyzfileを出す
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
        df_E = pd.DataFrame(columns = ['x','y','z','E','status','machine_type','file_name'])##いじる
        df_E.to_csv(auto_csv_path,index=False)##step3を二段階でやる場合二段階目ではinitをやらないので念のためmainにも組み込んでおく

    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(args.auto_dir,args.monomer_name,args.num_nodes,args.max_nodes,args.isTest)##argsの中身を取る
        time.sleep(1)

def listen(auto_dir,monomer_name,num_nodes,max_nodes,isTest):##args自体を引数に取るか中身をばらして取るかの違い
    maxnum_machine2 = 3#int(num_nodes/2) ##多分俺のために空けていてくださったので2 3にする
    fixed_param_keys = [];opt_param_keys = ['x','y','z']
    
    auto_csv_1 = os.path.join(auto_dir,'step1.csv');df_E_1 = pd.read_csv(auto_csv_1)
    df_prg_1 = df_E_1.loc[df_E_1['status']=='InProgress',fixed_param_keys+opt_param_keys+['machine_type','file_name']]
    len_prg_1=len(df_prg_1)
    for idx, row in df_prg_1.iterrows():
        params_dict1_ = row[fixed_param_keys + opt_param_keys + ['file_name']].to_dict()
        file_name1=params_dict1_['file_name']##辞書をつくってそこにopt_1とopt_2でファイル名作成
        log_filepath1 = os.path.join(*[auto_dir,'gaussian',file_name1])
        if not(os.path.exists(log_filepath1)):#logファイルが生成される直前だとまずいので
            continue
        E_list1=get_E(log_filepath1)
        if len(E_list1)!=1 :##get Eの長さは計算した分子の数
            continue
        else:
            len_prg_1-=1
            E1=float(E_list1[0])##8分子に向けてep1,ep2作成　ep1:b ep2:a
            df_E_1.loc[idx, ['E','status']] = [E1,'Done']
            df_E_1.to_csv(auto_csv_1,index=False)
            time.sleep(1)
            break#2つ同時に計算終わったりしたらまずいので一個で切る
    
    
    df_qw_1 = df_E_1[df_E_1['status'] == 'qw']
    len_queue = len_prg_1
    len_qw_1 = len(df_qw_1)
    margin = max_nodes - len_queue;machine_type = 1

    df_inpr_1 = df_E_1.loc[df_E_1['status']=='InProgress']
    machine_counts_1 = df_inpr_1['machine_type'].value_counts().to_dict();machine_counts_1.setdefault(1, 0);machine_counts_1.setdefault(2, 0)
    num_machine2 = machine_counts_1.get(2, 0)

    if len_qw_1 > 0 and margin > 0:# 進行中ジョブのマシンタイプをカウント
        for index, row in df_qw_1.iterrows():
            if margin == 0:
                break
            params_dict = row[fixed_param_keys + opt_param_keys].to_dict()# パラメータの辞書を作成
            if num_machine2 >= maxnum_machine2:
                machine_type = 1             
            else:
                machine_type = 2# マシンタイプの決定
                num_machine2 += 1
            file_name = exec_gjf(auto_dir, monomer_name, {**params_dict}, machine_type,isTest=isTest)# ジョブの実行 structure type
            len_queue += 1;margin -= 1
            df_E_1.at[index, 'machine_type'] = machine_type
            df_E_1.at[index, 'status'] = 'InProgress'
            df_E_1.at[index, 'file_name'] = file_name
        df_E_1.to_csv(auto_csv_1, index=False)
        time.sleep(1)#
    dict_matrix = get_params_dict(auto_dir,num_nodes)##更新分を流す a1~z2まで取得
    if len(dict_matrix)!=0:#終わりがまだ見えないなら
        for i in range(len(dict_matrix)):
            params_dict=dict_matrix[i]#print(params_dict)
            params_dict1 = {k: v for k, v in params_dict.items() if (k in fixed_param_keys) or (k in opt_param_keys)}
            alreadyCalculated = check_calc_status(auto_dir,params_dict)
            if not(alreadyCalculated):##ここで各点について1~3を実行しつつ余ったものをqwにぶち込む
                
                auto_csv_1 = os.path.join(auto_dir,'step1.csv');df_E_1 = pd.read_csv(auto_csv_1)
                df_sub_1 = filter_df(df_E_1, params_dict1)
                if len(df_sub_1) == 0:
                    print('debug1')
                    isAvailable = len_queue < max_nodes
                    if isAvailable:
                        machine2IsFull = num_machine2 >= maxnum_machine2
                        if machine2IsFull:
                            machine_type = 1
                        else:
                            machine_type = 2;num_machine2 += 1
                        file_name = exec_gjf(auto_dir, monomer_name, {**params_dict1}, machine_type,isTest=isTest);len_queue +=1
                        df_newline_1 = pd.Series({**params_dict1,'E1':0.,'machine_type':machine_type,'status':'InProgress','file_name':file_name})
                        df_E_new_1=pd.concat([df_E_1,df_newline_1.to_frame().T],axis=0,ignore_index=True);df_E_new_1.to_csv(auto_csv_1,index=False)
                        time.sleep(1)
                    else:
                        file_name = exec_gjf(auto_dir, monomer_name, {**params_dict1}, machine_type,isTest=True)
                        df_newline_1 = pd.Series({**params_dict1,'E1':0.,'machine_type':1,'status':'qw','file_name':file_name})
                        df_E_new_1=pd.concat([df_E_1,df_newline_1.to_frame().T],axis=0,ignore_index=True);df_E_new_1.to_csv(auto_csv_1,index=False)
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
    fixed_param_keys = [];opt_param_keys = ['x','y','z']
    
    #最初の立ち上がり時
    if len(df_init_params_inprogress) < num_nodes:
        #print(1)
        df_init_params_notyet = df_init_params[df_init_params['status']=='NotYet']
        for index in df_init_params_notyet.index:
            df_init_params = update_value_in_df(df_init_params,index,'status','InProgress')
            df_init_params.to_csv(init_params_csv,index=False)
            params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
            return [params_dict]
    dict_matrix=[]
    for index in df_init_params_inprogress.index:##こちら側はinit_params内のある業に関する探索が終わった際の新しい行での探索を開始するもの ###ここを改良すればよさそう
        df_init_params = pd.read_csv(init_params_csv)
        init_params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
        fixed_params_dict = df_init_params.loc[index,fixed_param_keys].to_dict()
        isDone, opt_params_matrix = get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict)
        if isDone:
            opt_params_dict={'x':np.round(opt_params_matrix[0][0],1),'y':np.round(opt_params_matrix[0][1],1),'z':np.round(opt_params_matrix[0][2],1)}
            # df_init_paramsのstatusをupdate
            df_init_params = update_value_in_df(df_init_params,index,'status','Done')
            if np.max(df_init_params.index) < index+1:##もうこれ以上は新しい計算は進まない
                status = 'Done'
            else:
                status = get_values_from_df(df_init_params,index+1,'status')
            df_init_params.to_csv(init_params_csv,index=False)
            
            if status=='NotYet':##計算が始まっていないものがあったらこの時点で開始する　ここでダメでもまた直にlistenでgrt_params_dictまでいけば新しいのが始まる            
                opt_params_dict = get_values_from_df(df_init_params,index+1,fixed_param_keys+opt_param_keys)
                df_init_params = update_value_in_df(df_init_params,index+1,'status','InProgress')
                df_init_params.to_csv(init_params_csv,index=False)
                dict_matrix.append({**fixed_params_dict,**opt_params_dict})
            else:
                continue

        else:
            for i in range(len(opt_params_matrix)):
                opt_params_dict={'x':np.round(opt_params_matrix[i][0],1),'y':np.round(opt_params_matrix[i][1],1),'z':np.round(opt_params_matrix[i][2],1)}
                d={**fixed_params_dict,**opt_params_dict}
                dict_matrix.append(d)
                    #print(d)
    return dict_matrix
        
def get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict):
    df_val = filter_df(df_cur, fixed_params_dict)
    x_init_prev = init_params_dict['x'];y_init_prev = init_params_dict['y'];z_init_prev = init_params_dict['z']
    while True:
        E_list=[];xyz_list=[]
        para_list=[]
        for x in [x_init_prev]:
            for y in [y_init_prev]:
                for z in [z_init_prev]:
                    x = np.round(x,1);y = np.round(y,1);z = np.round(z,1)
                    df_val_xyz = df_val[(df_val['x']==x)&(df_val['y']==y)&(df_val['z']==z)&(df_val['status']=='Done')]
                    if len(df_val_xyz)==0:
                        para_list.append([x,y,z])
                        continue
                    xyz_list.append([x,y,z]);E_list.append(df_val_xyz['E'].values[0])
        if len(para_list) != 0:
            return False,para_list
        x_init,y_init,z_init = xyz_list[np.argmin(np.array(E_list))]
        if x_init==x_init_prev and y_init==y_init_prev and z_init==z_init_prev:
            return True,[[x_init,y_init,z_init]]
        else:
            x_init_prev=x_init;y_init_prev=y_init;z_init_prev=z_init

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
    
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    parser.add_argument('--max-nodes',type=int,help='max nodes')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    print("----main process----")
    main_process(args)
    print("----finish process----")
    
import os
import subprocess
import time

################## USER INPUT BEGIN #######################
"""
    Referring: https://rtklibexplorer.wordpress.com/2022/01/05/batch-processing-rtklib-solutions-with-rnx2rtkp-and-python/
    Note     :
        该脚本为后处理批处理脚本，能同时跑不同配置的多组数组，具体的：
        会运行CONF_PATH配置文件夹下所有的配置，以及DATA_PATH文件夹下所有的RINEX数据，
        结果会以配置文件的名称（CONF_PATH）作为分类文件夹，存储多组结果数据（DATA_PATH内所有RINEX文件）。
    TBD      :
        1）需要解决DATA_DIR下属直接包含RINEX文件但不能处理的BUG
"""
BIN_DIR  = "../app/rnx2rtkp/msc/Release" # 可执行文件所在目录
DATA_DIR = "../data"                     # 输入文件目录。输入文件要求在data目录内部，且名称必须为rover.obs、base.obs、rover.nav
CONF_DIR = "../app/rnx2rtkp/conf"        # 配置文件夹。配置文件必须以.conf结尾
OUTS_DIR = "./output"                    # 输出文件
TRACE_LEVEL = 0                          # trace等级
##################  USER INPUT END  #######################

# set maximum number of simultaneous occurences of rnx2rtkp to run
max_windows = 14
poll_list = [None] * max_windows
slash_type = "\\" if os.name == "nt" else "/"

def ExecPostProc(bin, src, out, cfg):
    # create command to run solution
    if os.name == "nt":
        rtk_cmd = r'%s\rnx2rtkp.exe -ti 1.0 -x %d -y 0 -k %s -o %s.pos %s\rover.obs %s\base.obs %s\rover.nav' % \
                  (bin, TRACE_LEVEL, cfg, out, src, src, src)
    else:
        rtk_cmd = r'%s/rnx2rtkp -x %d -y 0 -k %s -o %s.pos %s/rover.obs %s/base.obs %s/rover.nav' % \
                  (bin, TRACE_LEVEL, cfg, out, src, src, src)
    # run command
    return subprocess.Popen(rtk_cmd, shell=True)

def GetConfigs(cfg_path):
    cfgs = []
    paths = os.walk(cfg_path)
    for _, _, files in paths:
        for f in files:
            if ".conf" in f:
                cfgs.append(f)
    return cfgs

def GenPostProcParas(data_path):
    input_dirs = []
    output_files = []
    
    paths = os.walk(data_path)
    for path, dir_lst, _ in paths:
        for dir_name in dir_lst:
            if (os.path.exists(os.path.join(path, dir_name)+"/rover.obs")):
                target_dir = os.path.join(path, dir_name)
                target_arr = target_dir.split(slash_type+"data"+slash_type)[-1].split(slash_type)
                pos_name = "_".join(target_arr)
                input_dirs.append(target_dir)
                output_files.append(pos_name)
    return input_dirs, output_files

def GetPollStat(polls):
    stats = []
    for p in polls:
        if p is not None and (p.poll() is None):
            stats.append(1)
        else:
            stats.append(0)
    return stats

def CheckIsFinished(polls):
    is_finished = True
    for p in polls:
        if p is not None and (p.poll() is None):
            is_finished = False
    return is_finished

# 删除output中所有的event文件
def DeleteFilesWithName(folder_path, keyword):
    for root, _, files in os.walk(folder_path):
        for file in files:
            if keyword in file:
                file_path = os.path.join(root, file)
                os.remove(file_path)

if __name__ == "__main__":
    start_time = time.time()
    cur_path = os.path.abspath(os.path.dirname(__file__))+"/"
    bin_path = os.path.abspath(cur_path + BIN_DIR)
    data_path = os.path.abspath(cur_path + DATA_DIR)
    conf_path = os.path.abspath(cur_path + CONF_DIR)
    output_path = os.path.abspath(cur_path + OUTS_DIR)

    input_paths, output_files = GenPostProcParas(data_path)
    cfgs = GetConfigs(conf_path)

    # loop through date folders
    ntasks = 0
    for src, out in zip(input_paths, output_files):
        for cfg in cfgs:
            is_idle = 0
            # if max windows open, wait for one to close
            while is_idle == 0:
                for j, p in enumerate(poll_list):
                    if p is None:
                        is_idle = 1
                    elif p.poll() is not None:
                        poll_list[j] = None
                        is_idle = 1
                time.sleep(1)  # wait here for existing window to close

            for j, p in enumerate(poll_list):
                if p is None:
                    is_idle = 1
                    out_path = os.path.join(output_path, cfg+slash_type+out)
                    poll_list[j] = ExecPostProc(bin_path, src, out_path, conf_path+slash_type+cfg)
                    ntasks += 1
                    print('Process: {:0>3d}/{:0>3d} '.format(ntasks, len(output_files)*len(cfgs)))
                    break            

    # Wait for all solutions to finish
    print('Waiting for solutions to complete ...')
    poll_sum = 0
    while CheckIsFinished(poll_list) == False:
        poll_stat = GetPollStat(poll_list)
        if (poll_sum != sum(poll_stat)):
            poll_sum = sum(poll_stat)
            print('Poll busy status({:0>3d}/{:0>3d}): '.format(poll_sum, max_windows), poll_stat)
        time.sleep(1)
    end_time = time.time()
    DeleteFilesWithName(output_path, "_events.pos")
    print('Batch processing completed in {0:.2f}s'.format(end_time-start_time))

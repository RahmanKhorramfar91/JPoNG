import multiprocessing
import sys
from subprocess import PIPE, Popen
import  os
import subprocess
import numpy as np;

 
def system(command):  
    process = Popen(
        args=command,
        stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
        bufsize=0,
        universal_newlines=True,
        shell=True
    )

    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()

    output = process.communicate()[0]
    exitCode = process.returncode

    if (exitCode == 0):
        return output
    else:
        return -1#raise ProcessException(command, exitCode, output)

    return process.communicate()[0]

if __name__ == "__main__":
    import time;
    from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor;
    from datetime import datetime;
    
    folder = 'joint_CF_with_extreme_days';
    net_size = 17;
    rep_days = [5,10,15];
    #rep_days = np.insert(rep_days,len(rep_days),365);
    case = [4];
    elec_scen = ['HE','HX','ME','RF'];
    emis_reduc_goal = [0.8,0.85,0.9,0.95];
    VRE_share = [0.5];
    solver_gap = 0.001;
    wall_clock_time_lim = 10; #hours;
    UC_active = 1;
    relax_UC_vars = 1;
    relax_int_vars = 0;
    solver_thread_num = 4;
    metal_air_cost = ['no-metal-air'];
    param_list=[];
    SuperClound_Thread = 96;
    for i2 in case:
        for i1 in rep_days:
            for i4 in emis_reduc_goal:
                for i3 in elec_scen:            
                    for i5 in VRE_share:
                        for i6 in metal_air_cost:                            
                            param = 'python Main.py '+folder+' '+str(net_size)+' '+ str(i1)+' '+str(i2)+' '+str(i3)+' '+str(i4)+' '+str(i5)+' '+str(solver_gap)+' '+str(wall_clock_time_lim)+' '+str(UC_active)+' '+str(relax_UC_vars)+' '+str(relax_int_vars)+' '+str(solver_thread_num)+' '+i6;
                            param_list.append(param);
                        #print(param)
    del i1,i2,i3,i4,i5;       
    #print(param_list)
    for i in range(int(np.ceil(solver_thread_num*len(param_list)/SuperClound_Thread))):
        j = int(SuperClound_Thread/solver_thread_num);
        with ThreadPoolExecutor() as executor:
            results = executor.map(system, param_list[i*j:(i+1)*j]);
    # for param in param_list:
    #      tmp = system(param)
    # print(tmp)

#%%
# import time
# from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
# from datetime import datetime

# now = datetime.now().time() # time object

# print("now =", now)
# def sleep_secs(seconds):
#   time.sleep(seconds)
#   print(f'{seconds} has been processed')

# secs_list = [2,4, 6, 8, 10, 12];
# with ThreadPoolExecutor() as executor:
#   results = executor.map(sleep_secs, secs_list)
#   print(results)

# now = datetime.now().time() # time object

# print("now =", now)
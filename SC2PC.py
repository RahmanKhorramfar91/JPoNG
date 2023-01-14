# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 13:55:42 2022

@author: Rahman Khorramfar
"""
import multiprocessing
import sys
from subprocess import PIPE, Popen
import  os
import subprocess;
import numpy as np;
def system(command):
    process = Popen(
        args=command,
        stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
        bufsize=0,
        universal_newlines=True,
        shell=True);

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


folder = 'joint_CF_with_extreme_days';
net_size = 6;
rep_days = [30];
#rep_days = np.insert(rep_days,len(rep_days),365);
case = [1,2];
elec_scen = ['2016'];
emis_reduc_goal = [0.0];
VRE_share = [0.0];
solver_gap = 0.001;
wall_clock_time_lim = 6; #hours;
UC_active = 1;
relax_UC_vars = 1;
relax_int_vars = 0;
solver_thread_num = 48;
metal_air_cost = ['no-metal-air'];
param_list=[];
SuperClound_Thread = 96;

for i2 in case:
    for i1 in rep_days:
        for i3 in elec_scen:
            for i4 in emis_reduc_goal:
                for i5 in VRE_share:
                    param = str(net_size)+'-'+ str(i1)+'-'+str(i3)+'-'+str(i2)+'-'+str(i4)+'-'+str(i5)+'.csv';
                    # if i2==1 and i3=='HM' and i4==0.85:continue;
                    # if i2==1 and i3=='RM' and i4==0.95:continue;
                    # if i2==1 and i3=='RM' and i4==0.85:continue;
                    
                    # if i2==3 and not ((i3=='RM' and i4==0.85) or (i2==3 and i3=='RM' and i4==0.95)):
                    #     continue;
                        
                    param_list.append(param);    
del i2,i3,i1,i4,i5;                    
system('cd '+os.getcwd());
system('scp -r rkhorramfar@txe1-login.mit.edu:/home/gridsan/rkhorramfar/JPoNG/JPoNG_Results.csv ./')

for param in param_list:
    system('scp -r rkhorramfar@txe1-login.mit.edu:/home/gridsan/rkhorramfar/JPonG-Python-Codes/'+param+' ./');
    
    

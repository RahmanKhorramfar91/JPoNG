# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:38:13 2023

@author: Rahman Khorramfar
"""

import time;
import gurobipy as gp;
from gurobipy import GRB, quicksum,LinExpr;
from Setting import Setting;
from SystemClasses import EV,GV;
from Setting import Setting;
from SystemClasses import EV,GV;
import sys;import os;import time;import csv;
from ProblemData import Enodes,Gnodes,Branches,CC_CCS;
from ProblemData import PipeLines,Plants,eStore,Other_input,state2zone_id,plant2sym;
from ProblemData import sym2plant,time_weights,zone_id2state,Exist_SVL,SVLs;
import numpy as np;

e_time_weight, g_time_weight, g_rep_days, e_rep_hrs, days_in_cluster,days2Medoid = time_weights(Setting.num_rep_days,Setting.rep_day_folder);

nE = len(Enodes); 
nG = len(Gnodes);
nBr = len(Branches);
nPipe = len(PipeLines);
nPlt = len(Plants);
neSt = len(eStore);
Te = range(len(e_time_weight));
Tg = range(len(g_time_weight));
nPipe = len(PipeLines);
nSVL = len(Exist_SVL);
nG = len(Gnodes);
FY = np.arange(365);
thermal_units = ["ng","OCGT","CCGT","CCGT-CCS","nuclear","nuclear-new"];
NG_units = ["ng","OCGT","CCGT","CCGT-CCS"];
VRE = ["solar","wind","wind-offshore-new","solar-UPV","wind-new"];# hydro not included



# solve a relaxed problem and get a good approximate solution to feed the SP
Setting.UC_active = False;
Setting.relax_UC_vars = True;
Setting.relax_int_vars = True;
Setting.copper_plate_approx = True; 
Model = gp.Model();
import Module;
Module.Power_System_Model(Model); # Power System
Module.NG_System_Model(Model);# NG System
Module.Coupling_constraints(Model);# Coupling Constraints
Model.modelSense = GRB.MINIMIZE;
if Setting.emis_case==1:Model.setObjective(EV.e_system_cost);    
else:Model.setObjective(GV.g_system_cost+ EV.e_system_cost);
Model.setParam('MIPGap', Setting.solver_gap);
Model.setParam('Timelimit', Setting.wall_clock_time_lim);
Model.setParam('Presolve',2); # -1 to 2
Model.optimize();
print(f"First step objective value: {Model.ObjVal}");
Module.Get_var_vals(Model);













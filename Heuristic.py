# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 17:25:25 2023

@author: Rahman Khorramfar
"""
import time;
import gurobipy as gp;
from gurobipy import GRB, quicksum,LinExpr;
from Setting import Setting;
from SystemClasses import EV,GV;
import sys; 
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


#%% 3-step approach: copper-plate with relaxed integrality and no UC
 
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



# Step 2: set the values of integer variable to the nearest integer to their values in the first step
Setting.UC_active = 1;
Setting.relax_UC_vars = True;
Setting.relax_int_vars = 0;
Setting.copper_plate_approx = 1; 
Model = gp.Model();
import Module;
Module.Power_System_Model(Model); # Power System
Module.NG_System_Model(Model);# NG System
Module.Coupling_constraints(Model);# Coupling Constraints
Model.modelSense = GRB.MINIMIZE;
if Setting.emis_case==1:Model.setObjective(EV.e_system_cost); 
else:Model.setObjective(GV.g_system_cost+ EV.e_system_cost);
for n in range(nE):
    for i in range(nPlt):
        if i in thermal_units:
            Model.addConstr(EV.Xop[n,i] == np.round(EV.Xop_val[n,i]));
            Model.addConstr(EV.Xest[n,i] == np.round(EV.Xest_val[n,i]));

# for b in range(nBr):
#     Model.addConstr(EV.Ze == np.round(EV.Ze_val[b]));
    

Model.setParam('MIPGap', Setting.solver_gap);
Model.setParam('Timelimit', Setting.wall_clock_time_lim);
Model.setParam('Presolve',2); # -1 to 2
Model.optimize();
print(f"Second step objective value: {Model.ObjVal}");
Module.Get_var_vals(Model);



# Step 3: set the values of integer variable to the nearest integer to their values in the first step
Setting.UC_active = 1;
Setting.relax_UC_vars = True;
Setting.relax_int_vars = True;
Setting.copper_plate_approx = 0; 
Model = gp.Model();
import Module;
Module.Power_System_Model(Model); # Power System
Module.NG_System_Model(Model);# NG System
Module.Coupling_constraints(Model);# Coupling Constraints
Model.modelSense = GRB.MINIMIZE;
if Setting.emis_case==1:Model.setObjective(EV.e_system_cost); 
else:Model.setObjective(GV.g_system_cost+ EV.e_system_cost);
for n in range(nE):
    for i in range(nPlt):
        if i in thermal_units:
            Model.addConstr(EV.Xop[n,i] == np.round(EV.Xop_val[n,i]));
            Model.addConstr(EV.Xest[n,i] == np.round(EV.Xest_val[n,i]));
        else: 
            Model.addConstr(EV.Xop[n,i] <= np.round(EV.Xop_val[n,i]));

for b in range(nBr):
    Model.addConstr(EV.Ze == np.round(EV.Ze_val[b]));
    

Model.setParam('MIPGap', Setting.solver_gap);
Model.setParam('Timelimit', Setting.wall_clock_time_lim);
Model.setParam('Presolve',2); # -1 to 2
Model.optimize();
print(f"Third step objective value: {Model.ObjVal}");
Module.Get_var_vals(Model);


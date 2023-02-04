# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:39:41 2023

@author: Rahman Khorramfar
"""

import time;
import gurobipy as gp;
from gurobipy import GRB, quicksum,LinExpr;
from Setting import Setting;
from SystemClasses import EV,GV;
import sys; 
from Setting import Setting;
from SystemClasses import EV,GV;
import os;import time;import csv;
from ProblemData import Enodes,Gnodes,Branches,CC_CCS;
from ProblemData import PipeLines,Plants,eStore,Other_input,state2zone_id,plant2sym;
from ProblemData import sym2plant,time_weights,zone_id2state,Exist_SVL,SVLs;

import gurobipy as gp;
from gurobipy import GRB, quicksum,LinExpr;
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


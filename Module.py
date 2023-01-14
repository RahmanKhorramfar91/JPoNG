# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 14:40:30 2022

Power system on hourly resolution and NG on daily

"""
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
thermal_units = ["ng","CT","CC","CC-CCS","nuclear","nuclear-new"];
NG_units = ["ng","CT","CC","CC-CCS"];
VRE = ["solar","wind","wind-offshore-new","solar-UPV","wind-new"];# hydro not included

def Power_System_Model(Model):    
    #% define decision variables 
    EV.Xop = Model.addVars(nE,nPlt,vtype=GRB.INTEGER);
    EV.Xest = Model.addVars(nE,nPlt,vtype=GRB.INTEGER);
    EV.Xdec = Model.addVars(nE,nPlt,vtype=GRB.INTEGER);
    EV.Ze = Model.addVars(nBr,vtype=GRB.BINARY);
    
    if Setting.UC_active:
        EV.X = Model.addVars(nE,len(Te), nPlt,vtype=GRB.INTEGER);
        EV.Xup= Model.addVars(nE,len(Te), nPlt,vtype=GRB.INTEGER);
        EV.Xdown = Model.addVars(nE,len(Te), nPlt,vtype=GRB.INTEGER);
    
    if Setting.relax_int_vars:
        EV.Xop = Model.addVars(nE,nPlt,vtype=GRB.CONTINUOUS);
        EV.Xest = Model.addVars(nE,nPlt,vtype=GRB.CONTINUOUS);
        EV.Xdec = Model.addVars(nE,nPlt,vtype=GRB.CONTINUOUS);
        EV.Ze = Model.addVars(nBr,vtype=GRB.CONTINUOUS);
        
    if (Setting.relax_UC_vars or Setting.relax_int_vars) and (Setting.UC_active):
        EV.X = Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
        EV.Xup= Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
        EV.Xdown = Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
        
    EV.prod= Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
    EV.theta = Model.addVars(nE,len(Te),lb=np.zeros((nE,len(Te)))-GRB.INFINITY,ub=np.zeros((nE,len(Te)))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    EV.Shed = Model.addVars(nE,len(Te),vtype=GRB.CONTINUOUS);
    
    EV.YeCD = Model.addVars(nE,neSt,vtype=GRB.CONTINUOUS);
    EV.YeLev = Model.addVars(nE,neSt,vtype=GRB.CONTINUOUS);
    #EV.YeStr = Model.addVars(nE,neSt,vtype=GRB.BINARY);
    EV.kappa_capt = Model.addVars(nE,len(Te),vtype=GRB.CONTINUOUS);
    EV.kappa_pipe = Model.addVars(nE,vtype=GRB.CONTINUOUS);
    
    EV.eSch =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSdis =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSlev =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSday = Model.addVars(nE,len(FY),neSt,vtype=GRB.CONTINUOUS);
    EV.eSrem = Model.addVars(nE,len(FY),neSt,lb=np.zeros((nE,len(FY),neSt))-GRB.INFINITY,ub=np.zeros((nE,len(FY),neSt))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    
    EV.flowE =  Model.addVars(nBr,len(Te),lb=np.zeros((nBr,len(Te)))-GRB.INFINITY,ub=np.zeros((nBr,len(Te)))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    
    if Setting.expansion_allowed==False:
        Model.addConstrs(EV.Xest[n,i]==0 for n  in range(nE) for i in range(nPlt) )
        Model.addConstrs(EV.Ze[b]==0 for b in range(nBr));
        
    EV.est_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.est_trans_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.decom_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.FOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.startup_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.VOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.nuc_fuel_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.gas_fuel_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.shedding_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.elec_storage_cost1 = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.elec_storage_cost2 = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.emis_amount = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.CCS_cost =  Model.addVar(vtype=GRB.CONTINUOUS);
    EV.trans_FOM_cost =  Model.addVar(vtype=GRB.CONTINUOUS);
    EV.e_system_cost=Model.addVar(vtype=GRB.CONTINUOUS);
    
    
    #% Set some variables
    #1) existing types can not be established because there are new equivalent types
    #2) new types cannot be decommissioned
    #3) offshore only allowed in certain nodes
    
    Model.addConstrs(EV.Xest[n,i]==0 for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==1)
    Model.addConstrs(EV.Xdec[n,i]==0 for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==0)
    Model.addConstrs(EV.Xest[n,sym2plant['wind-offshore-new']]==0 for n in range(nE) if Enodes[n].offshore_allowd==0);
   
    #% Electricity System Objective Function    
    # cost components
    ccs_cost=0;
    est_cost = LinExpr(quicksum(Plants[i].est_cost_coef[n]*Plants[i].capex*EV.Xest[n,i] for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==0));
    dec_cost = LinExpr(quicksum(Plants[i].decom_cost*EV.Xdec[n,i] for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==1));
    fom_cost = LinExpr(quicksum(Plants[i].FOM*EV.Xop[n,i] for n in range(nE) for i in range(nPlt)));
    vom_cost = LinExpr(quicksum(e_time_weight[t]*Plants[i].VOM*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt)));
    nuc_fuel_cost = LinExpr(quicksum(e_time_weight[t]*Plants[i].heat_rate *Other_input.Nuclear_price*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt) if (Plants[i].Type=="nuclear" or Plants[i].Type=="nuclear-new") ));
    
    ccs_cost = LinExpr(quicksum(CC_CCS.node2str_dis[n]*CC_CCS.pipe_capex*EV.kappa_pipe[n] for n in range(nE)) + quicksum(e_time_weight[t]*CC_CCS.str_capex*EV.kappa_capt[n,t] for n in range(nE) for t in Te)); 
    
    gas_fuel_cost = 0;
    gas_fuel_cost = LinExpr(quicksum(Other_input.NG_price*e_time_weight[t]*Plants[i].heat_rate*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt) if plant2sym[i] in NG_units))
    startup_cost = 0;
    if Setting.UC_active:
        startup_cost = LinExpr(quicksum(e_time_weight[t]*Plants[i].startup_cost*EV.Xup[n,t,i] for n in range(nE) for t in Te for i in range(nPlt)));
    shed_cost = LinExpr(quicksum(e_time_weight[t]*Setting.e_shed_penalty*EV.Shed[n,t] for n in range(nE) for t in Te));
    strg_cost1 = LinExpr(quicksum(eStore[0].est_coef*(eStore[0].power_capex*EV.YeCD[n,0]+eStore[0].energy_capex*EV.YeLev[n,0]) for n in range(nE) ));
    strg_cost1 += LinExpr(quicksum((eStore[0].pFOM*EV.YeCD[n,0]+eStore[0].eFOM*EV.YeLev[n,0]) for n in range(nE)));
    strg_cost2=0;
    if Setting.Metal_air_storage_cost!='no-metal-air':
        strg_cost2 = LinExpr(quicksum(eStore[1].est_coef*(eStore[1].power_capex*EV.YeCD[n,1]+eStore[1].energy_capex*EV.YeLev[n,1]) for n in range(nE)));
        strg_cost2 += LinExpr(quicksum((eStore[1].pFOM*EV.YeCD[n,1]+eStore[1].eFOM*EV.YeLev[n,1]) for n in range(nE)));

    tran_cost = LinExpr(quicksum(Branches[b].est_coef*Other_input.trans_unit_cost*Branches[b].maxFlow*Branches[b].length*EV.Ze[b] for b in range(nBr)));
    trans_Fom = LinExpr(quicksum(Branches[b].trans_FOM*Branches[b].maxFlow*Branches[b].length for b in range(nBr) if Branches[b].is_exist==1));
    trans_Fom += LinExpr(quicksum(Branches[b].trans_FOM*Branches[b].maxFlow*Branches[b].length*EV.Ze[b] for b in range(nBr) if Branches[b].is_exist==0));
    # total power system cost function
    if Setting.emis_case==1: # add gas fuel cost only if Case=0 (power system only)
        e_total_cost = gas_fuel_cost+est_cost+dec_cost+fom_cost+vom_cost+nuc_fuel_cost+ccs_cost+startup_cost+shed_cost+strg_cost1+strg_cost2+tran_cost+trans_Fom;
    else:
        e_total_cost = est_cost+dec_cost+fom_cost+vom_cost+nuc_fuel_cost+ccs_cost+startup_cost+shed_cost+strg_cost1+strg_cost2+tran_cost+trans_Fom;

    #e_total_cost = est_cost+dec_cost+fom_cost+vom_cost+shed_cost;
    
    Model.addConstr(EV.est_cost == est_cost);
    Model.addConstr(EV.decom_cost == dec_cost);
    Model.addConstr(EV.FOM_cost == fom_cost);
    Model.addConstr(EV.VOM_cost == vom_cost);
    Model.addConstr(EV.nuc_fuel_cost == nuc_fuel_cost);
    Model.addConstr(EV.startup_cost == startup_cost);
    Model.addConstr(EV.shedding_cost == shed_cost);
    Model.addConstr(EV.elec_storage_cost1 == strg_cost1);
    Model.addConstr(EV.elec_storage_cost2 == strg_cost2);
    Model.addConstr(EV.est_trans_cost == tran_cost);
    Model.addConstr(EV.trans_FOM_cost == trans_Fom);
    Model.addConstr(EV.gas_fuel_cost == gas_fuel_cost);
    Model.addConstr(EV.CCS_cost == ccs_cost);
    Model.addConstr(EV.e_system_cost == e_total_cost);    
    
    #% Electricity System Constraints
    # C1: number of generation units at each node
    Model.addConstrs(EV.Xop[n,i] == Enodes[n].Init_plt_count[i]+EV.Xest[n,i]-EV.Xdec[n,i] for n in range(nE) for i in range(nPlt));
        
    # C2, C3, C4, C5: UC,  production limit, ramping for thermal units (ng, CT, CC, CC-CCS, nuclear)    
    Model.addConstrs(EV.prod[n,t,i]<= Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
    if Setting.UC_active:
        # UC:
        Model.addConstrs(EV.X[n,t,i]-EV.X[n,t-1,i]==EV.Xup[n,t,i]-EV.Xdown[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units));
        # ramping:
        Model.addConstrs(EV.prod[n,t,i]-EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*(EV.X[n,t,i]-EV.Xup[n,t,i])+max(Plants[i].min_output, Plants[i].ramp_rate)*Plants[i].nameplate_cap*EV.Xup[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
        Model.addConstrs(-EV.prod[n,t,i]+EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*(EV.X[n,t,i]-EV.Xup[n,t,i])+max(Plants[i].min_output, Plants[i].ramp_rate)*Plants[i].nameplate_cap*EV.Xup[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
        # production limit
        Model.addConstrs(EV.prod[n,t,i]>= Plants[i].min_output*Plants[i].nameplate_cap*EV.X[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
        Model.addConstrs(EV.prod[n,t,i]<= Plants[i].nameplate_cap*EV.X[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
        
        Model.addConstrs(EV.X[n,t,i] <= EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
    else:
        Model.addConstrs(EV.prod[n,t,i]-EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
        Model.addConstrs(-EV.prod[n,t,i]+EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
    
    # C5, C6: flow limit for electricity
    Model.addConstrs(EV.flowE[b,t]<=Branches[b].maxFlow for b in range(nBr) for t in Te if(Branches[b].is_exist==1));
    Model.addConstrs(-EV.flowE[b,t]<=Branches[b].maxFlow for b in range(nBr) for t in Te if(Branches[b].is_exist==1))
        
    Model.addConstrs(EV.flowE[b,t]<=Branches[b].maxFlow*EV.Ze[b] for b in range(nBr) for t in Te if(Branches[b].is_exist==0));
    Model.addConstrs(-EV.flowE[b,t]<=Branches[b].maxFlow*EV.Ze[b] for b in range(nBr) for t in Te if(Branches[b].is_exist==0))
    
    # C7: power balance  
    
    if Setting.copper_plate_approx:
        Model.addConstrs(quicksum(EV.prod[n,t,i] for i in range(nPlt) for n in range(nE))+quicksum(EV.eSdis[n,t,r]-EV.eSch[n,t,r] for r in range(neSt) for n in range(nE))+quicksum(EV.Shed[n,t] for n in range(nE)) == quicksum(Enodes[n].demand[e_rep_hrs[t]]+CC_CCS.node2str_dis[n]*CC_CCS.elec_req_pipe*EV.kappa_pipe[n] + (CC_CCS.node2str_dis[n]/CC_CCS.comp_dis)*CC_CCS.elec_req_pump*EV.kappa_capt[n,t] for n in range(nE)) for t in Te);
    else:        
        Model.addConstrs(quicksum(EV.prod[n,t,i] for i in range(nPlt))-quicksum(Enodes[n].arc_sign[b]*EV.flowE[Enodes[n].arcs[b],t] for b in range(len(Enodes[n].arcs)))+quicksum(EV.eSdis[n,t,r]-EV.eSch[n,t,r] for r in range(neSt))+EV.Shed[n,t] == Enodes[n].demand[e_rep_hrs[t]]+CC_CCS.node2str_dis[n]*CC_CCS.elec_req_pipe*EV.kappa_pipe[n] + (CC_CCS.node2str_dis[n]/CC_CCS.comp_dis)*CC_CCS.elec_req_pump*EV.kappa_capt[n,t]  for n in range(nE) for t in Te);
  
    
    # C8: flow equation
    Model.addConstrs(EV.flowE[b,t]==Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t]) for b in range(nBr) for t in Te if Branches[b].is_exist==1);
    Model.addConstrs(EV.flowE[b,t]-Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t])<=10e7*(1-EV.Ze[b])  for b in range(nBr) for t in Te if Branches[b].is_exist==0);
    Model.addConstrs(-EV.flowE[b,t]+Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t])<=10e7*(1-EV.Ze[b])  for b in range(nBr) for t in Te if Branches[b].is_exist==0);
        
    # # C9: phase angle (theta) limits. already applied in the definition of the variable
    # Model.addConstrs(EV.theta[n,t]<=Other_input.pi for n in range(nE) for t in Te);
    # Model.addConstrs(-EV.theta[n,t]<=Other_input.pi for n in range(nE) for t in Te);
    Model.addConstrs(EV.theta[0,t]==0 for t in Te);
    
    # C10: VRE production according to capacity factors
    Model.addConstrs(EV.prod[n,t,i]<=Enodes[n].cap_factors[e_rep_hrs[t],i]* Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt) for t in Te);
    
    # C11: demand curtainlment constraint
    Model.addConstrs(EV.Shed[n,t]<= Enodes[n].demand[e_rep_hrs[t]] for  n in range(nE) for t in Te);
    
    # C12: RPS constraints
    Model.addConstr(quicksum(e_time_weight[t]*EV.prod[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if(plant2sym[i] in VRE)) >= Setting.VRE_share*quicksum(e_time_weight[t]*Enodes[n].demand[e_rep_hrs[t]] for n in range(nE) for t in Te));
    
    # C14,C15,C16 storage constraints,  find the starting hour of each rep day
    start_hours = [24*k for k in range(len(g_rep_days))];
    end_hours = [(k+1)*24-1 for k in range(len(g_rep_days))];
    Model.addConstrs(EV.eSlev[n,t,r]-(1-eStore[r].self_discharge)*EV.eSlev[n,t-1,r]==eStore[r].eff_ch*EV.eSch[n,t,r]-EV.eSdis[n,t,r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt) for t in Te if t not in start_hours);   
    # start and ending of storage should be the same in case of using rep. days
    if len(g_rep_days)<365:   
        Model.addConstrs(EV.eSlev[n,start_hours[k],r]==(1-eStore[r].self_discharge)*(EV.eSlev[n,end_hours[k],r]-EV.eSrem[n,g_rep_days[k],r]) + eStore[r].eff_ch*EV.eSch[n,start_hours[k],r]-EV.eSdis[n,start_hours[k],r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt) for k in range(len(g_rep_days)));
    else:
        Model.addConstrs(EV.eSlev[n,0,r]==eStore[r].eff_ch*EV.eSch[n,0,r]-EV.eSdis[n,0,r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt));

    # storage carry over constraints for each period (day)
    Model.addConstrs(EV.eSday[n,tau+1,r] == (1-24*eStore[r].self_discharge)*EV.eSday[n,tau,r]+EV.eSrem[n,days2Medoid[tau],r] for n in range(nE) for r in range(neSt) for tau in FY if tau!=FY[-1]);
    Model.addConstrs(EV.eSday[n,0,r] == (1-24*eStore[r].self_discharge)*EV.eSday[n,FY[-1],r]+EV.eSrem[n,days2Medoid[FY[-1]],r] for n in range(nE) for r in range(neSt));
    Model.addConstrs(EV.eSday[n,g_rep_days[tau],r] == EV.eSlev[n,end_hours[tau],r]-EV.eSrem[n,g_rep_days[tau],r] for n in range(nE) for r in range(neSt) for tau in range(len(g_rep_days)));
    Model.addConstrs(EV.eSday[n,0,r] == 0  for n in range(nE) for r in range(neSt));
    # no storage carry over for Li-ion battery
    Model.addConstrs(EV.eSrem[n,tau,0]==0 for n in range(nE) for tau in FY);
    # Model.addConstrs(EV.eSrem[n,tau,1]==0 for n in range(nE) for tau in FY);
    
    Model.addConstrs(EV.YeCD[n,r]>= EV.eSdis[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    Model.addConstrs(EV.YeCD[n,r]>= EV.eSch[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    Model.addConstrs(EV.YeLev[n,r]>= EV.eSlev[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    
    # productin capacity for each new type    
    Model.addConstr(quicksum(Plants[sym2plant['solar-UPV']].nameplate_cap*EV.Xop[n,sym2plant['solar-UPV']]+Plants[sym2plant['solar']].nameplate_cap*EV.Xop[n,sym2plant['solar']] for n in range(nE) )<= Other_input.type_prod_lim[0]);
    Model.addConstr(quicksum(Plants[sym2plant['wind-new']].nameplate_cap*EV.Xop[n,sym2plant['wind-new']]+Plants[sym2plant['wind']].nameplate_cap*EV.Xop[n,sym2plant['wind']] for n in range(nE) )<= Other_input.type_prod_lim[1]);
    Model.addConstr(quicksum(Plants[sym2plant['wind-offshore-new']].nameplate_cap*EV.Xop[n,sym2plant['wind-offshore-new']] for n in range(nE) )<= Other_input.type_prod_lim[2]);
    Model.addConstr(quicksum(Plants[sym2plant['nuclear-new']].nameplate_cap*EV.Xop[n,sym2plant['nuclear-new']]+Plants[sym2plant['nuclear']].nameplate_cap*EV.Xop[n,sym2plant['nuclear']] for n in range(nE) )<= Other_input.type_prod_lim[3]);
    # CSS constraints
    Model.addConstrs(EV.kappa_capt[n,t]==Other_input.NG_emission*(Plants[sym2plant['CCGT-CCS']].co2_capture_rate)*Plants[sym2plant['CCGT-CCS']].heat_rate*EV.prod[n,t,sym2plant['CCGT-CCS']] for n in range(nE) for t in Te);    
    Model.addConstrs(EV.kappa_pipe[n] >= EV.kappa_capt[n,t] for n in range(nE) for t in Te);
    Model.addConstr(quicksum(EV.kappa_capt[n,t] for n in range(nE) for t in Te)<= Other_input.carbon_str_cap);

    # capacity reserve margin constraint    
    Model.addConstrs(quicksum(Enodes[n].cap_factors[e_rep_hrs[t],i]*Plants[i].nameplate_cap*EV.Xop[n,i] for n in range(nE) for i in range(nPlt)) >= (1+Setting.CRM_reserve)*quicksum(Enodes[n].demand[e_rep_hrs[t]] for n in range(nE)) for t in Te);
    

def NG_System_Model(Model): # for the full year
    
    GV.Zg = Model.addVars(nPipe,vtype=GRB.BINARY);
    GV.ZgDec = Model.addVars(nPipe,vtype=GRB.BINARY);
    GV.ZgOp = Model.addVars(nPipe,vtype=GRB.BINARY);
    if Setting.relax_int_vars==True:
        GV.Zg = Model.addVars(nPipe,vtype=GRB.CONTINUOUS);
        GV.ZgDec = Model.addVars(nPipe,vtype=GRB.CONTINUOUS);
        GV.ZgOp = Model.addVars(nPipe,vtype=GRB.CONTINUOUS);    
    if Setting.expansion_allowed==False:
        Model.addConstrs(GV.Zg[b]==0 for b in range(nPipe));

    GV.Xvpr = Model.addVars(nSVL,vtype=GRB.CONTINUOUS);
    GV.Xstr = Model.addVars(nSVL,vtype=GRB.CONTINUOUS);
    GV.Sstr = Model.addVars(nSVL,len(FY),vtype = GRB.CONTINUOUS);
    GV.Svpr = Model.addVars(nSVL,len(FY),vtype = GRB.CONTINUOUS);
    GV.Sliq = Model.addVars(nSVL,len(FY),vtype = GRB.CONTINUOUS);
    GV.supply = Model.addVars(nG,len(FY),vtype = GRB.CONTINUOUS);
    GV.Shed =  Model.addVars(nG,len(FY),vtype = GRB.CONTINUOUS);
    GV.RNG_use = Model.addVars(nG,len(FY),vtype = GRB.CONTINUOUS);
    GV.flowGG =  Model.addVars(nPipe,len(FY),vtype = GRB.CONTINUOUS);
    GV.flowGE =  Model.addVars(nG,nE,len(FY),vtype = GRB.CONTINUOUS);
    GV.flowGL =  Model.addVars(nG,nSVL,len(FY),vtype = GRB.CONTINUOUS);
    GV.flowVG =  Model.addVars(nSVL,nG,len(FY),vtype = GRB.CONTINUOUS);
    
    GV.inv_str_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.inv_pipe_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.pipe_FOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.pipe_Decom_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.shed_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.RNG_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.fom_str_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.import_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.emis_amount = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.g_system_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    
    # NG System Objective Function
    inv_pipe = LinExpr(quicksum(PipeLines[b].inv_coef*PipeLines[b].length*Other_input.pipe_per_mile*GV.Zg[b] for b in range(nPipe)));   
    ng_import = LinExpr(quicksum(GV.supply[k,tau]*Other_input.NG_price for k in range(nG) for tau in FY));
    storage_inv = LinExpr(quicksum(SVLs[0].inv_coef*(SVLs[0].capex*GV.Xstr[j]+SVLs[1].capex*GV.Xvpr[j]) for j in range(nSVL)));
    storage_FOM = LinExpr(quicksum(SVLs[1].FOM*(Exist_SVL[j].vap_cap+GV.Xvpr[j])+SVLs[0].FOM*(Exist_SVL[j].str_cap+GV.Xstr[j]) for j in range(nSVL)));
    ng_shed = LinExpr(quicksum(Setting.g_shed_penalty*GV.Shed[k,tau] for k in range(nG) for tau in FY));
    rng_import = LinExpr(quicksum(Other_input.RNG_price*GV.RNG_use[k,tau]for k in range(nG) for tau in FY));
    pipe_FOM = LinExpr(quicksum(PipeLines[b].FOM*PipeLines[b].length*GV.ZgOp[b] for b in range(nPipe)));  
    pipe_Dec = LinExpr(quicksum(PipeLines[b].decom*PipeLines[b].length*GV.ZgDec[b] for b in range(nPipe)));  
    
    
    ng_total_cost = inv_pipe+ng_import+storage_inv+storage_FOM+ng_shed+rng_import+pipe_FOM+pipe_Dec;
    Model.addConstr(GV.inv_pipe_cost==inv_pipe);
    Model.addConstr(GV.import_cost==ng_import);
    Model.addConstr(GV.inv_str_cost==storage_inv);
    Model.addConstr(GV.fom_str_cost==storage_FOM);
    Model.addConstr(GV.shed_cost==ng_shed);
    Model.addConstr(GV.RNG_cost==rng_import);
    Model.addConstr(GV.pipe_FOM_cost == pipe_FOM);
    Model.addConstr(GV.pipe_Decom_cost == pipe_Dec);
    Model.addConstr(GV.g_system_cost==ng_total_cost);
    
    # NG System Constraints
    #C1, C2: flow limit for NG
    # Model.addConstrs(GV.flowGG[i,tau]<=PipeLines[i].Cap for i in range(nPipe) for tau in FY if PipeLines[i].is_exist==1);
    Model.addConstrs(GV.flowGG[i,tau]<=PipeLines[i].Cap*GV.ZgOp[i] for i in range(nPipe) for tau in FY);
    
    # C3: flow balance, NG node (no shedding allowed as RNG is considered)
    Model.addConstrs(GV.supply[k,tau]
                      -quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_exp)
                      +quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_imp)
                      -quicksum(GV.flowGE[k,n,tau] for n in Gnodes[k].adjE)
                      +quicksum((GV.flowVG[j,k,tau]-GV.flowGL[k,j,tau]) for j in Gnodes[k].adjS)
                      +GV.Shed[k,tau]
                      +GV.RNG_use[k,tau]
                     == Gnodes[k].demand[tau] for k in range(nG) for tau in FY);                             
    
    # inforce all flowGE of the same cluster to take the equal values
    Model.addConstrs(GV.flowGE[k,n,tau]==GV.flowGE[k,n,days2Medoid[tau]] for k in range(nG) for n in Gnodes[k].adjE for tau in FY)
    
    # C3,C4: injection (supply) limit and curtailment limit
    Model.addConstrs(GV.supply[k,tau]<= Gnodes[k].injU for k in range(nG) for tau in FY);    
    Model.addConstrs(GV.Shed[k,tau]+GV.RNG_use[k,tau]<= Gnodes[k].demand[tau] for k in range(nG) for tau in FY);
    Model.addConstrs(GV.RNG_use[k,tau]<= Gnodes[k].injU for k in range(nG) for tau in FY);    

    # C5: storage balance (storage strats empty)
    Model.addConstrs(GV.Sstr[j,tau]==Exist_SVL[j].str_cap*0+GV.Sliq[j,tau]-GV.Svpr[j,tau]/SVLs[1].eff_disCh for j in range(nSVL) for tau in FY if tau==0);
    Model.addConstrs(GV.Sstr[j,tau]==(1-SVLs[0].BOG)*GV.Sstr[j,tau-1]+GV.Sliq[j,tau]-GV.Svpr[j,tau]/SVLs[1].eff_disCh for j in range(nSVL) for tau in FY if tau>0);
    
    # C6,8: calculate Sliq, Svpr
    for j in range(nSVL):
        for tau in FY:
            NG_adj = [];
            for k in range(nG): 
                for j2 in Gnodes[k].adjS:
                    if j2==j:
                        NG_adj.append(k);
                        
            Model.addConstr(GV.Sliq[j,tau]==quicksum(GV.flowGL[k,j,tau] for k in NG_adj));
            Model.addConstr(GV.Svpr[j,tau]==quicksum(GV.flowVG[j,k,tau] for k in NG_adj));
                
    # C6: Sliq limit
    Model.addConstrs(GV.Sliq[j,tau]<=Exist_SVL[j].liq_cap for j in range(nSVL) for tau in FY);
    
    # C9: Svpr limit
    Model.addConstrs(GV.Svpr[j,tau]<=Exist_SVL[j].vap_cap+GV.Xvpr[j] for j in range(nSVL) for tau in FY);
    
    # C10: Sstr limit
    Model.addConstrs(GV.Sstr[j,tau]<=Exist_SVL[j].str_cap+GV.Xstr[j] for j in range(nSVL) for tau in FY);
    
    # C11: Operational Pipelines
    Model.addConstrs(GV.ZgOp[b] == PipeLines[b].is_exist-GV.ZgDec[b] for b in range(nPipe) if PipeLines[b].is_exist==1);
    Model.addConstrs(GV.ZgOp[b] == GV.Zg[b] for b in range(nPipe) if PipeLines[b].is_exist==0);

    
def Coupling_constraints(Model):
    if (Setting.emis_case !=1): #power system only
        #Model.addConstrs(GV.flowGE[k,n,g_rep_days[tau]]==quicksum(Plants[i].heat_rate*EV.prod[n,t,i] for t in range(tau*24,(tau+1)*24) for i in range(nPlt) if plant2sym[i] in NG_units) for k in range(nG) for tau in Tg for n in Gnodes[k].adjE);
        Model.addConstrs(quicksum(GV.flowGE[k,n,g_rep_days[tau]] for k in range(nG) if n in Gnodes[k].adjE)==quicksum(Plants[i].heat_rate*EV.prod[n,t,i] for t in range(tau*24,(tau+1)*24) for i in range(nPlt) if plant2sym[i] in NG_units) for n in range(nE) for tau in Tg);
    
    #ex_xi = LinExpr(quicksum(g_time_weight[tau]*GV.flowGE[k,n,g_rep_days[tau]] for k in range(nG) for tau in Tg for n in Gnodes[k].adjE));
    e_emis = LinExpr(quicksum(e_time_weight[t]*Plants[i].emission*Plants[i].heat_rate*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt) if plant2sym[i] in NG_units));
    #g_emis = LinExpr(quicksum(g_time_weight[tau]*Other_input.NG_emission*GV.supply[k,tau] for k in range(nG) for tau in Tg))-Other_input.NG_emission*ex_xi;
    
    g_emis = LinExpr(quicksum(Other_input.NG_emission*(Gnodes[k].demand[tau]-GV.RNG_use[k,tau]-GV.Shed[k,tau]) for k in range(nG) for tau in FY));
    Model.addConstr(EV.emis_amount==e_emis);
    Model.addConstr(GV.emis_amount==g_emis);
    
    if Setting.emis_case==1: #power system only
        Model.addConstr(EV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.e_emis_lim);
    if Setting.emis_case==2:  # JPoNG with no RNG and emis const on power system    
        Model.addConstr(EV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.e_emis_lim);
        Model.addConstrs(GV.RNG_use[k,tau]==0 for k in range(nG) for tau in FY);
    if Setting.emis_case==3:#global emission limit with no RNG
        Model.addConstr(EV.emis_amount+GV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.CO2_emission_1990);
        Model.addConstrs(GV.RNG_use[k,tau]==0 for k in range(nG) for tau in FY);
    if Setting.emis_case==4:  #global emission limit
        Model.addConstr(EV.emis_amount+GV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.CO2_emission_1990);

    if Setting.emis_case==5: # JPoNG with separate emission constraints    
        Model.addConstr(EV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.e_emis_lim);
        Model.addConstr(GV.emis_amount<=(1-Setting.emis_reduc_goal)*Setting.g_emis_lim);
        #Model.addConstrs(GV.RNG_use[k,tau]==0 for k in range(nG) for tau in Tg);
        


def Get_var_vals(Model):
    
    EV.Xop_val = Model.getAttr('x',EV.Xop);
    EV.Xest_val = Model.getAttr('x',EV.Xest);
    EV.Xdec_val = Model.getAttr('x',EV.Xdec);
    # EV.X_val = Model.getAttr('x',EV.X);
    # print(EV.X_val);
    #EV.Xup_val = Model.getAttr('x',EV.Xup);
    
    EV.prod_val = Model.getAttr('x',EV.prod);
    EV.flowE_val = Model.getAttr('x',EV.flowE);
    EV.YeCD_val = Model.getAttr('x',EV.YeCD);
    EV.YeLev_val = Model.getAttr('x',EV.YeLev);
    EV.eSlev_val = Model.getAttr('x',EV.eSlev);
    EV.Shed_val = Model.getAttr('x',EV.Shed);
    EV.Ze_val = Model.getAttr('x',EV.Ze);
    EV.eSdis_val = Model.getAttr('x',EV.eSdis);
    EV.eSch_val = Model.getAttr('x',EV.eSch);
    EV.eSday_val =  Model.getAttr('x',EV.eSday);
    EV.eSrem_val =  Model.getAttr('x',EV.eSrem);
    # for n in range(nE):
    #     for tau in FY:
    #         for r in range(neSt):
    #             # if EV.eSday_val[n,tau,r]>0:
    #             #     print(f"eSday_val[{n},{tau},{r}] = {EV.eSday_val[n,tau,r]}");
    #             if EV.eSrem_val[n,tau,r]>0:
    #                 print(f"eSrem_val[{n},{tau},{r}] = {EV.eSrem_val[n,tau,r]}");

    # for n in range(nE):
    #     for t in Te:
    #         for r in range(neSt):
    #             # if EV.eSday_val[n,tau,r]>0:
    #             #     print(f"eSday_val[{n},{tau},{r}] = {EV.eSday_val[n,tau,r]}");
    #             if EV.eSlev_val[n,t,r]>0:
    #                 print(f"eSlev[{n},{t},{r}] = {EV.eSlev_val[n,t,r]}");

    
    EV.kappa_capt_val = Model.getAttr('x',EV.kappa_capt);
    EV.kappa_pipe_val =  Model.getAttr('x',EV.kappa_pipe);    
    EV.est_cost_val = EV.est_cost.X;
    EV.decom_cost_val = EV.decom_cost.X;
    EV.FOM_cost_val = EV.FOM_cost.X;
    EV.VOM_cost_val = EV.VOM_cost.X;
    EV.nuc_fuel_cost_val = EV.nuc_fuel_cost.X;
    EV.startup_cost_val = EV.startup_cost.X;
    EV.shedding_cost_val = EV.shedding_cost.X;
    EV.elec_storage_cost_val1 = EV.elec_storage_cost1.X;
    EV.elec_storage_cost_val2 = EV.elec_storage_cost2.X;
    EV.est_trans_cost_val = EV.est_trans_cost.X;
    EV.trans_FOM_cost_val = EV.trans_FOM_cost.X;
    EV.emis_amount_val = EV.emis_amount.X;
    EV.gas_fuel_cost_val = EV.gas_fuel_cost.X;
    EV.CCS_cost_val = EV.CCS_cost.X;
    EV.e_system_cost_val = EV.e_system_cost.X;
    #if Setting.emis_case==1:
    EV.e_system_cost_val = EV.e_system_cost.X;#-EV.gas_fuel_cost.X;
    
    # print(Model.getAttr('x',EV.theta));
    # print(Model.getAttr('x',GV.Svpr));
    # print(Model.getAttr('x',GV.flowGG));
    
    GV.Zg_val = Model.getAttr('x',GV.Zg);
    GV.ZgOp_val = Model.getAttr('x',GV.ZgOp);
    GV.ZgDec_val = Model.getAttr('x',GV.ZgDec);
    GV.supply_val = Model.getAttr('x',GV.supply);
    #GV.RNG_supply_val = Model.getAttr('x',GV.supply);
    GV.Shed_val = Model.getAttr('x',GV.Shed);
    GV.RNG_use_val = Model.getAttr('x',GV.RNG_use);
    GV.flowGE_val = Model.getAttr('x',GV.flowGE);
    GV.g_system_cost_val = GV.g_system_cost.X;
    GV.flowGL_val = Model.getAttr('x',GV.flowGL);
    # note this one
    GV.import_cost_val = GV.import_cost.X;
    if Setting.emis_case==1:
        GV.import_cost_val = GV.import_cost.X-EV.gas_fuel_cost.X;

    
    
    GV.RNG_cost_val = GV.RNG_cost.X;
    GV.fom_str_cost_val = GV.fom_str_cost.X;
    GV.shed_cost_val = GV.shed_cost.X;
    GV.inv_pipe_cost_val = GV.inv_pipe_cost.X;
    GV.pipe_FOM_cost_val = GV.pipe_FOM_cost.X;
    GV.pipe_Decom_cost_val = GV.pipe_Decom_cost.X;
    GV.emis_amount_val = GV.emis_amount.X;
    GV.inv_str_cost_val = GV.inv_str_cost.X;
    

def Publish_results(s_time,MIP_gap):
    num_e_str1 = 0; str_lev1 = 0;str_cap1=0;
    num_e_str2 = 0; str_lev2 = 0;str_cap2=0;
    num_ze = 0;total_shed=0;total_flow=0;
    [num_e_str1 := num_e_str1+1 if EV.YeLev_val[n,0]>0 else num_e_str1+0 for n in range(nE)];
    [str_lev1 := str_lev1+EV.YeLev_val[n,0] for n in range(nE)];
    [str_cap1 := str_cap1 +EV.YeCD_val[n,0] for n in range(nE)];
    if Setting.Metal_air_storage_cost!='no-metal-air':
        [num_e_str2 := num_e_str2+1 if EV.YeLev_val[n,1]>0 else num_e_str2+0 for n in range(nE)];
        [str_lev2 := str_lev2+EV.YeLev_val[n,1] for n in range(nE)];
        [str_cap2 := str_cap2 +EV.YeCD_val[n,1] for n in range(nE)];
    
    [num_ze := num_ze+1 if EV.Ze_val[b]>0 else num_ze+0 for b in range(nBr)];
    [total_shed := total_shed+EV.Shed_val[n,t] for n in range(nE) for t in Te];
    [total_flow:= total_flow+abs(EV.flowE_val[b,t]) for b in range(nBr) for t in Te];
    total_flow = total_flow/2;
    
    # for b in range(nPipe):
    #     print(f"line {b}, ets: {GV.Zg_val[b]}, Op: {GV.ZgOp_val[b]}, decom: {GV.ZgDec_val[b]}")
     
    pr = np.zeros(nPlt);
    est = np.zeros(nPlt);
    dec = np.zeros(nPlt);
    
    for i in range(nPlt):
        s1 = 0;
        [s1:= s1+EV.prod_val[n,t,i]*e_time_weight[t] for n in range(nE) for t in Te];
        pr[i] = s1;s1 = 0;
        [s1:= s1+EV.Xest_val[n,i] for n in range(nE)];
        est[i] = s1;s1 = 0;
        [s1:= s1+EV.Xdec_val[n,i] for n in range(nE)];
        dec[i] = s1;
        
      
    
    num_est_pipe = 0;total_ng_shed = 0;total_rng=0;total_fge=0;total_fgl =0;
    num_op_pipe = 0; num_decom_pipe = 0;
    [num_est_pipe := num_est_pipe+1 if GV.Zg_val[p]>0 else num_est_pipe+0 for p in range(nPipe)];
    [num_op_pipe := num_op_pipe+1 if GV.ZgOp_val[p]>0 else num_op_pipe+0 for p in range(nPipe)];
    [num_decom_pipe := num_decom_pipe+1 if GV.ZgDec_val[p]>0 else num_decom_pipe+0 for p in range(nPipe)];

    [total_ng_shed := total_ng_shed+GV.Shed_val[n,tau] for n in range(nG) for tau in FY];
    [total_rng := total_rng+GV.RNG_use_val[n,tau] for n in range(nG) for tau in FY];
    #total_rng = total_rng/Setting.poss_gas_consump;
    [total_fge := total_fge+GV.flowGE_val[j,n,tau] for j in range(nG) for n in range(nE) for tau in FY];
    [total_fgl := total_fgl+GV.flowGL_val[k,j,tau] for j in range(nSVL) for k in range(nG) for tau in FY];
    elapsed = time.time()-s_time;
    header0 = ['Power_network_size','cluster_method','Rep-Days',
               'Emis-case','Elec_scenario', 'reduc-goal','RPS','UC-active?',
              'UC-rlx?','int-vars-rlx?','Metal-air-cost',
              'MI-gap(%)', 'Run time(sec)','Total-cost',
              'Power-cost','est-cost','decom-cost','FOM','VOM',
              'nuc-fuel-cost','gas-fuel-cost','startup-cost','e-shed-cost',
              'storage1-cost','storage2-cost','tran-est-cost','trans_FOM','CCS-cost',
              'emission_e','num-est-tran',
              'total-str1-lev','total-str1-cap','num-est-str1',
              'total-str2-lev','total-str2-cap','num-est-str2',
              'total-flow','NG-cost','NG-import-cost','LCDF-import-cost',
              'inv-storage','FOM-storage','g-shed-cost','pipe-est-cost',
              'pipe_FOM','pipe_Dec_cost',
              'emission_g','num-est-pipe','num_decom_pipe','num_operational_pipe',
              'total-ng-shed',
              'total-LCDF-import','total-flowGE','total-flowGL','Production:',
              'ng','solar','wind','hydro','nuclear','OCGT','CCGT','CCGT-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new','established:', # production
              'ng','solar','wind','hydro','nuclear','OCGT','CCGT','CCGT-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new','decommissioned:', # established
              'ng','solar','wind','hydro','nuclear','OCGT','CCGT','CCGT-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new', # decommissioned
              ];
    row = []; row.append(Setting.Power_network_size); row.append(Setting.rep_day_folder);
    row.append(Setting.num_rep_days);row.append(Setting.emis_case);
    row.append(Setting.electrification_scenario);
    row.append(Setting.emis_reduc_goal);
    row.append(Setting.VRE_share);row.append(Setting.UC_active);
    row.append(Setting.relax_UC_vars);row.append(Setting.relax_int_vars);
    row.append(Setting.Metal_air_storage_cost);
    row.append(MIP_gap);row.append(elapsed);
    row.append(EV.e_system_cost_val+GV.g_system_cost_val);
 
    row.append(EV.e_system_cost_val); row.append(EV.est_cost_val);
    row.append(EV.decom_cost_val); row.append(EV.FOM_cost_val);row.append(EV.VOM_cost_val);
    row.append(EV.nuc_fuel_cost_val);row.append(EV.gas_fuel_cost_val);
    row.append(EV.startup_cost_val);
    row.append(EV.shedding_cost_val);
    row.append(EV.elec_storage_cost_val1);
    row.append(EV.elec_storage_cost_val2);
    row.append(EV.est_trans_cost_val); row.append(EV.trans_FOM_cost_val);
    row.append(EV.CCS_cost_val);
    row.append(EV.emis_amount_val);row.append(num_ze);
    row.append(str_lev1);row.append(str_cap1);
    row.append(num_e_str1); row.append(str_lev2);row.append(str_cap2);
    row.append(num_e_str2);     
    row.append(total_flow);

    row.append(GV.g_system_cost_val);row.append(GV.import_cost_val);
    row.append(GV.RNG_cost_val);
    row.append(GV.inv_str_cost_val);
    row.append(GV.fom_str_cost_val);row.append(GV.shed_cost_val);
    row.append(GV.inv_pipe_cost_val);
    row.append(GV.pipe_FOM_cost_val);row.append(GV.pipe_Decom_cost_val);
    row.append(GV.emis_amount_val);
    row.append(num_est_pipe);row.append(num_decom_pipe); row.append(num_op_pipe);
    row.append(total_ng_shed);
    row.append(total_rng);row.append(total_fge);row.append(total_fgl);
    row.append('Production:');
    #print(f"total fge: {total_fge}");
    for i in range(nPlt):
        row.append(pr[i]);
    row.append('established:');
    for i in range(nPlt):
        row.append(est[i]);
    row.append('decommissioned:');
    for i in range(nPlt):
        row.append(dec[i]);
        
    with open(os.getcwd()+'/JPoNG_Results.csv','a',encoding='UTF8',newline='') as f:
        writer = csv.writer(f);
        if Setting.print_result_header:
            writer.writerow(header0);
        writer.writerow(row);
        f.close();
    
    name =os.getcwd()+'/'+ str(Setting.Power_network_size)+'-'+ str(Setting.num_rep_days)+'-'+Setting.electrification_scenario+'-'+str(Setting.emis_case)+'-'+str(Setting.emis_reduc_goal)+'-'+str(Setting.VRE_share)+'.csv';

    

    if Setting.print_all_vars:
        pls = ['ng','solar','wind','hydro','nuclear','OCGT','CCGT','CCGT-CCS','solar-UPV','wind-new','wind-offshore','nuclear-new'];
        header = []; header.append('Hourly Generation:');
        for p in pls: header.append(p);
        header.append('Charge1'); header.append('Discharge1'); 
        header.append('Charge2'); header.append('Discharge2'); 
        header.append('from_g_node');header.append('to_g_node');
        header.append('ZgEst');header.append('ZgOp');
        header.append('from_e_node');header.append('to_e_node');
        header.append('Ze');
        header.append('NG_supply:');
        for g in range(nG): header.append(str(g));
        header.append('LCDF_supply');
        for g in range(nG): header.append(str(g));
        header.append('nodal_generation:');
        for p in pls: header.append(p);
        header.append('nodal_plt_est:');
        for p in pls: header.append(p);
        header.append('nodal_plt_decom:');
        for p in pls: header.append(p);
        header.append('NG_marginal_price:');
        for g in range(nG): header.append(str(g));
              
        row0 = [['']*(len(header)+10) for i in range(max(368,Setting.num_rep_days*24+3))];
        #row0 = np.zeros((max(96,Setting.num_rep_days*24+3),77))-1;
        #row0[:] = np.nan;
        col = 1; 
        #row0 = [[] for t in range(len(Te))];
        for t in Te:
            for i in range(nPlt):
                r0 = [EV.prod_val[n,t,i]  for n in range(nE)];                
                row0[t][i+col] = np.round(sum(r0));
            r0 = [EV.eSch_val[n,t,0] for n in range(nE)];            
            row0[t][nPlt] = np.round(sum(r0));
            r0 = [EV.eSdis_val[n,t,0] for n in range(nE)];            
            row0[t][nPlt+1] = np.round(sum(r0));
            
            if Setting.Metal_air_storage_cost!='no-metal-air':
                r0 = [EV.eSch_val[n,t,1] for n in range(nE)];            
                row0[t][nPlt+2] = np.round(sum(r0));
                r0 = [EV.eSdis_val[n,t,1] for n in range(nE)];            
                row0[t][nPlt+3] = np.round(sum(r0));
                
        col = nPlt+5;#charge and discharge              
        for p in range(nPipe):            
            row0[p][col] = PipeLines[p].from_node;
            row0[p][col+1] = PipeLines[p].to_node;
            row0[p][col+2] = GV.Zg_val[p];
            row0[p][col+3] = GV.ZgOp_val[p];
            
            
        col += 4;
        for br in range(nBr):
            row0[br][col] = Branches[br].from_node;
            row0[br][col+1] = Branches[br].to_node;
            row0[br][col+2] = EV.Ze_val[br];
            
        col += 4;
        for tau in FY:
            for k in range(nG):
                row0[tau][k+col] = GV.supply_val[k,tau];
            
        col += nG+1;
        for tau in FY:
            for k in range(nG):                
                row0[tau][k+col] = GV.RNG_use_val[k,tau];            
        col += nG+1;               
        
        for n in range(nE):
            for i in range(nPlt):
                r0 = [EV.prod_val[n,t,i]*e_time_weight[t] for t in Te];                
                row0[n][i+col] = np.round(sum(r0));
        col += nPlt+1;    
        for n in range(nE):
            for i in range(nPlt):                
                #row0[n].append(EV.Xest_val[n,i]);
                row0[n][i+col] = EV.Xest_val[n,i];
        col += nPlt+1;
        for n in range(nE):
            for i in range(nPlt):                
                #row0[n].append(EV.Xdec_val[n,i]);
                row0[n][i+col] = EV.Xdec_val[n,i];
        
        col += nPlt+1;       
        # for n in range(nE):
        #     row0[0][col+n]= np.round(EV.kappa_pipe_val[n],1);

        # col += nE;
        # for n in range(nE):
        #     for t in Te:
        #         row0[t][col+n] = np.round(EV.kappa_capt_val[n,t],1);        
        # col += nE;
        for k in range(nG):
            for tau in Tg:
                row0[tau][col+k] = GV.marginal_prices_val[k,g_rep_days[tau]];
                
        with open(name,'a',encoding='UTF8',newline='') as fid:
            writer=csv.writer(fid);
            writer.writerow(header0);
            writer.writerow(row);
            writer.writerow(header);
            for r in range(len(row0)):
                writer.writerow(row0[r]);
            fid.close();

   
def Get_marginal_prices():
    Model = gp.Model();
    EV.Xup = Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
    EV.Xdown = Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
    EV.X = Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
    
    EV.prod= Model.addVars(nE,len(Te), nPlt,vtype=GRB.CONTINUOUS);
    EV.theta = Model.addVars(nE,len(Te),lb=np.zeros((nE,len(Te)))-GRB.INFINITY,ub=np.zeros((nE,len(Te)))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    EV.Shed = Model.addVars(nE,len(Te),vtype=GRB.CONTINUOUS);
    
    EV.YeCD = Model.addVars(nE,neSt,vtype=GRB.CONTINUOUS);
    EV.YeLev = Model.addVars(nE,neSt,vtype=GRB.CONTINUOUS);
    #EV.YeStr = Model.addVars(nE,neSt,vtype=GRB.BINARY);
    EV.kappa_capt = Model.addVars(nE,len(Te),vtype=GRB.CONTINUOUS);
    EV.kappa_pipe = Model.addVars(nE,vtype=GRB.CONTINUOUS);
    
    EV.eSch =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSdis =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSlev =  Model.addVars(nE,len(Te),neSt,vtype=GRB.CONTINUOUS);
    EV.eSday = Model.addVars(nE,len(FY),neSt,vtype=GRB.CONTINUOUS);
    EV.eSrem = Model.addVars(nE,len(FY),neSt,lb=np.zeros((nE,len(FY),neSt))-GRB.INFINITY,ub=np.zeros((nE,len(FY),neSt))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    
    EV.flowE =  Model.addVars(nBr,len(Te),lb=np.zeros((nBr,len(Te)))-GRB.INFINITY,ub=np.zeros((nBr,len(Te)))+GRB.INFINITY,vtype=GRB.CONTINUOUS);
    

    
    EV.est_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.est_trans_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.decom_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.FOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.startup_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.VOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.nuc_fuel_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.gas_fuel_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.shedding_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.elec_storage_cost1 = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.elec_storage_cost2 = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.emis_amount = Model.addVar(vtype=GRB.CONTINUOUS);
    EV.CCS_cost =  Model.addVar(vtype=GRB.CONTINUOUS);
    EV.trans_FOM_cost =  Model.addVar(vtype=GRB.CONTINUOUS);
    EV.e_system_cost=Model.addVar(vtype=GRB.CONTINUOUS);
           
    #% Electricity System Objective Function    
    # cost components
    est_cost = LinExpr(quicksum(Plants[i].est_cost_coef[n]*Plants[i].capex* EV.Xest_val[n,i] for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==0));
    dec_cost = LinExpr(quicksum(Plants[i].decom_cost* EV.Xest_val[n,i] for n in range(nE) for i in range(nPlt) if Plants[i].is_exist==1));
    fom_cost = LinExpr(quicksum(Plants[i].FOM* EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt)));
    vom_cost = LinExpr(quicksum(e_time_weight[t]*Plants[i].VOM*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt)));
    nuc_fuel_cost = LinExpr(quicksum(e_time_weight[t]*Plants[i].heat_rate *Other_input.Nuclear_price*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt) if (Plants[i].Type=="nuclear" or Plants[i].Type=="nuclear-new") ));
    ccs_cost = LinExpr(quicksum(CC_CCS.node2str_dis[n]*CC_CCS.pipe_capex*EV.kappa_pipe[n] for n in range(nE)) + quicksum(e_time_weight[t]*CC_CCS.str_capex*EV.kappa_capt[n,t] for n in range(nE) for t in Te)); 

    gas_fuel_cost = 0;
    gas_fuel_cost = LinExpr(quicksum(Other_input.NG_price*e_time_weight[t]*Plants[i].heat_rate*EV.prod[n,t,i] for n in range(nE) for t in Te for i in range(nPlt) if plant2sym[i] in NG_units))
    startup_cost = 0;
    if Setting.UC_active:
        startup_cost = LinExpr(quicksum(e_time_weight[t]*Plants[i].startup_cost*EV.Xup[n,t,i] for n in range(nE) for t in Te for i in range(nPlt)));
    shed_cost = LinExpr(quicksum(e_time_weight[t]*Setting.e_shed_penalty*EV.Shed[n,t] for n in range(nE) for t in Te));
    strg_cost1 = LinExpr(quicksum(eStore[0].est_coef*(eStore[0].power_capex*EV.YeCD[n,0]+eStore[0].energy_capex*EV.YeLev[n,0]) for n in range(nE) ));
    strg_cost1 += LinExpr(quicksum((eStore[0].pFOM*EV.YeCD[n,0]+eStore[0].eFOM*EV.YeLev[n,0]) for n in range(nE)));
    strg_cost2=0;
    if Setting.Metal_air_storage_cost!='no-metal-air':
        strg_cost2 = LinExpr(quicksum(eStore[1].est_coef*(eStore[1].power_capex*EV.YeCD[n,1]+eStore[1].energy_capex*EV.YeLev[n,1]) for n in range(nE)));
        strg_cost2 += LinExpr(quicksum((eStore[1].pFOM*EV.YeCD[n,1]+eStore[1].eFOM*EV.YeLev[n,1]) for n in range(nE)));

    tran_cost = LinExpr(quicksum(Branches[b].est_coef*Other_input.trans_unit_cost*Branches[b].maxFlow*Branches[b].length*EV.Ze_val[b] for b in range(nBr)));
    trans_Fom = LinExpr(quicksum(Branches[b].trans_FOM*Branches[b].maxFlow*Branches[b].length for b in range(nBr) if Branches[b].is_exist==1));
    trans_Fom += LinExpr(quicksum(Branches[b].trans_FOM*Branches[b].maxFlow*Branches[b].length*EV.Ze_val[b] for b in range(nBr) if Branches[b].is_exist==0));
    # total power system cost function
    if Setting.emis_case==1: # add gas fuel cost only if Case=0 (power system only)
        e_total_cost = gas_fuel_cost+est_cost+dec_cost+fom_cost+vom_cost+nuc_fuel_cost+ccs_cost+startup_cost+shed_cost+strg_cost1+strg_cost2+tran_cost+trans_Fom;
    else:
        e_total_cost = est_cost+dec_cost+fom_cost+vom_cost+nuc_fuel_cost+ccs_cost+startup_cost+shed_cost+strg_cost1+strg_cost2+tran_cost+trans_Fom;
   
    Model.addConstr(EV.est_cost == est_cost);
    Model.addConstr(EV.decom_cost == dec_cost);
    Model.addConstr(EV.FOM_cost == fom_cost);
    Model.addConstr(EV.VOM_cost == vom_cost);
    Model.addConstr(EV.nuc_fuel_cost == nuc_fuel_cost);
    Model.addConstr(EV.startup_cost == startup_cost);
    Model.addConstr(EV.shedding_cost == shed_cost);
    Model.addConstr(EV.elec_storage_cost1 == strg_cost1);
    Model.addConstr(EV.elec_storage_cost2 == strg_cost2);
    Model.addConstr(EV.est_trans_cost == tran_cost);
    Model.addConstr(EV.trans_FOM_cost == trans_Fom);
    Model.addConstr(EV.gas_fuel_cost == gas_fuel_cost);
    Model.addConstr(EV.CCS_cost == ccs_cost);
    Model.addConstr(EV.e_system_cost == e_total_cost);
    
    # #% Electricity System Constraints
    # # C1: number of generation units at each node
    # #Model.addConstrs( EV.Xop_val[n,i] == Enodes[n].Init_plt_count[i]+ EV.Xest_val[n,i]- EV.Xest_val[n,i] for n in range(nE) for i in range(nPlt));
        
    # # C2, C3, C4, C5: UC,  production limit, ramping for thermal units (ng, CT, CC, CC-CCS, nuclear)    
    Model.addConstrs(EV.prod[n,t,i]<= Plants[i].nameplate_cap* EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
    if Setting.UC_active:
        # UC:
        Model.addConstrs(EV.X[n,t,i]- EV.X[n,t-1,i]== EV.Xup[n,t,i]- EV.Xdown[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units));
        # ramping:
        Model.addConstrs(EV.prod[n,t,i]-EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*(EV.X[n,t,i]-EV.Xup[n,t,i])+max(Plants[i].min_output, Plants[i].ramp_rate)*Plants[i].nameplate_cap*EV.Xup[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
        Model.addConstrs(-EV.prod[n,t,i]+EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap*(EV.X[n,t,i]-EV.Xup[n,t,i])+max(Plants[i].min_output, Plants[i].ramp_rate)*Plants[i].nameplate_cap*EV.Xup[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
        # production limit
        Model.addConstrs(EV.prod[n,t,i]>= Plants[i].min_output*Plants[i].nameplate_cap* EV.X[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
        Model.addConstrs(EV.prod[n,t,i]<= Plants[i].nameplate_cap*EV.X[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
        
        Model.addConstrs( EV.X[n,t,i] <=  EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (plant2sym[i] in thermal_units))
    else:
        Model.addConstrs(EV.prod[n,t,i]-EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap* EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
        Model.addConstrs(-EV.prod[n,t,i]+EV.prod[n,t-1,i]<=Plants[i].ramp_rate* Plants[i].nameplate_cap* EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt) for t in Te if (t>0 and plant2sym[i] in thermal_units))
    
    # C5, C6: flow limit for electricity
    Model.addConstrs(EV.flowE[b,t]<=Branches[b].maxFlow for b in range(nBr) for t in Te if(Branches[b].is_exist==1));
    Model.addConstrs(-EV.flowE[b,t]<=Branches[b].maxFlow for b in range(nBr) for t in Te if(Branches[b].is_exist==1))
        
    Model.addConstrs(EV.flowE[b,t]<=Branches[b].maxFlow* EV.Ze_val[b] for b in range(nBr) for t in Te if(Branches[b].is_exist==0));
    Model.addConstrs(-EV.flowE[b,t]<=Branches[b].maxFlow* EV.Ze_val[b] for b in range(nBr) for t in Te if(Branches[b].is_exist==0))
    
    # C7: power balance  
    if Setting.copper_plate_approx:
        Model.addConstrs(quicksum(EV.prod[n,t,i] for i in range(nPlt) for n in range(nE))+quicksum(EV.eSdis[n,t,r]-EV.eSch[n,t,r] for r in range(neSt) for n in range(nE))+quicksum(EV.Shed[n,t] for n in range(nE)) == quicksum(Enodes[n].demand[e_rep_hrs[t]] for n in range(nE)) for t in Te);
    else:        
        Model.addConstrs(quicksum(EV.prod[n,t,i] for i in range(nPlt))+quicksum(Enodes[n].arc_sign[b]*EV.flowE[Enodes[n].arcs[b],t] for b in range(len(Enodes[n].arcs)))+quicksum(EV.eSdis[n,t,r]-EV.eSch[n,t,r] for r in range(neSt))+EV.Shed[n,t] == Enodes[n].demand[e_rep_hrs[t]] for n in range(nE) for t in Te);
    
    # C8: flow equation
    Model.addConstrs(EV.flowE[b,t]==Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t]) for b in range(nBr) for t in Te if Branches[b].is_exist==1);
    Model.addConstrs(EV.flowE[b,t]-Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t])<=10e7*(1- EV.Ze_val[b])  for b in range(nBr) for t in Te if Branches[b].is_exist==0);
    Model.addConstrs(-EV.flowE[b,t]+Branches[b].suscept*(EV.theta[Branches[b].to_node,t]-EV.theta[Branches[b].from_node,t])<=10e7*(1- EV.Ze_val[b])  for b in range(nBr) for t in Te if Branches[b].is_exist==0);
        
    # # C9: phase angle (theta) 
    Model.addConstrs(EV.theta[0,t]==0 for t in Te);
    
    # C10: VRE production according to capacity factors
    Model.addConstrs(EV.prod[n,t,i]<=Enodes[n].cap_factors[e_rep_hrs[t],i]* Plants[i].nameplate_cap* EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt) for t in Te);
    
    # C11: load shedding constraint
    Model.addConstrs(EV.Shed[n,t]<= Enodes[n].demand[e_rep_hrs[t]] for  n in range(nE) for t in Te);
    
    # C12: RPS constraints
    Model.addConstr(quicksum(e_time_weight[t]*EV.prod[n,t,i] for n in range(nE) for i in range(nPlt) for t in Te if(plant2sym[i] in VRE)) >= Setting.VRE_share*quicksum(e_time_weight[t]*Enodes[n].demand[e_rep_hrs[t]] for n in range(nE) for t in Te));
    
    # C14,C15,C16 storage constraints
    start_hours = [24*k for k in range(len(g_rep_days))];
    end_hours = [(k+1)*24-1 for k in range(len(g_rep_days))];
    Model.addConstrs(EV.eSlev[n,t,r]-(1-eStore[r].self_discharge)*EV.eSlev[n,t-1,r]==eStore[r].eff_ch*EV.eSch[n,t,r]-EV.eSdis[n,t,r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt) for t in Te if t not in start_hours);   
    # start and ending of storage should be the same in case of using rep. days
    if len(g_rep_days)<365:   
        Model.addConstrs(EV.eSlev[n,start_hours[k],r]==(1-eStore[r].self_discharge)*(EV.eSlev[n,end_hours[k],r]-EV.eSrem[n,g_rep_days[k],r]) + eStore[r].eff_ch*EV.eSch[n,start_hours[k],r]-EV.eSdis[n,start_hours[k],r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt) for k in range(len(g_rep_days)));
    else:
        Model.addConstrs(EV.eSlev[n,0,r]==eStore[r].eff_ch*EV.eSch[n,0,r]-EV.eSdis[n,0,r]/eStore[r].eff_disCh for n in range(nE) for r in range(neSt));

    # storage carry over constraints for each period (day)
    Model.addConstrs(EV.eSday[n,tau+1,r] == (1-24*eStore[r].self_discharge)*EV.eSday[n,tau,r]+EV.eSrem[n,days2Medoid[tau],r] for n in range(nE) for r in range(neSt) for tau in FY if tau!=FY[-1]);
    Model.addConstrs(EV.eSday[n,0,r] == (1-24*eStore[r].self_discharge)*EV.eSday[n,FY[-1],r]+EV.eSrem[n,days2Medoid[FY[-1]],r] for n in range(nE) for r in range(neSt));
    Model.addConstrs(EV.eSday[n,g_rep_days[tau],r] == EV.eSlev[n,end_hours[tau],r]-EV.eSrem[n,g_rep_days[tau],r] for n in range(nE) for r in range(neSt) for tau in range(len(g_rep_days)));
    Model.addConstrs(EV.eSday[n,0,r] == 0  for n in range(nE) for r in range(neSt));
    # no storage carry over for Li-ion battery
    Model.addConstrs(EV.eSrem[n,tau,0]==0 for n in range(nE) for tau in FY);
    # Model.addConstrs(EV.eSrem[n,tau,1]==0 for n in range(nE) for tau in FY);
    
    Model.addConstrs(EV.YeCD[n,r]>= EV.eSdis[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    Model.addConstrs(EV.YeCD[n,r]>= EV.eSch[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    Model.addConstrs(EV.YeLev[n,r]>= EV.eSlev[n,t,r] for n in range(nE) for r in range(neSt) for t in Te if t>0);
    
    # productin capacity for each new type    
    Model.addConstr(quicksum(Plants[sym2plant['solar-UPV']].nameplate_cap*EV.Xop_val[n,sym2plant['solar-UPV']]+Plants[sym2plant['solar']].nameplate_cap*EV.Xop_val[n,sym2plant['solar']] for n in range(nE) )<= Other_input.type_prod_lim[0]);
    Model.addConstr(quicksum(Plants[sym2plant['wind-new']].nameplate_cap*EV.Xop_val[n,sym2plant['wind-new']]+Plants[sym2plant['wind']].nameplate_cap*EV.Xop_val[n,sym2plant['wind']] for n in range(nE) )<= Other_input.type_prod_lim[1]);
    Model.addConstr(quicksum(Plants[sym2plant['wind-offshore-new']].nameplate_cap*EV.Xop_val[n,sym2plant['wind-offshore-new']] for n in range(nE) )<= Other_input.type_prod_lim[2]);
    Model.addConstr(quicksum(Plants[sym2plant['nuclear-new']].nameplate_cap*EV.Xop_val[n,sym2plant['nuclear-new']]+Plants[sym2plant['nuclear']].nameplate_cap*EV.Xop_val[n,sym2plant['nuclear']] for n in range(nE) )<= Other_input.type_prod_lim[3]);
    # CSS constraints
    Model.addConstrs(EV.kappa_capt[n,t]==Other_input.NG_emission*(Plants[sym2plant['CCGT-CCS']].co2_capture_rate)*Plants[sym2plant['CCGT-CCS']].heat_rate*EV.prod[n,t,sym2plant['CCGT-CCS']] for n in range(nE) for t in Te);    
    Model.addConstrs(EV.kappa_pipe[n] >= EV.kappa_capt[n,t] for n in range(nE) for t in Te);
    Model.addConstr(quicksum(EV.kappa_capt[n,t] for n in range(nE) for t in Te)<= Other_input.carbon_str_cap);

    # capacity reserve margin constraint    
    Model.addConstrs(quicksum(Enodes[n].cap_factors[e_rep_hrs[t],i]*Plants[i].nameplate_cap*EV.Xop_val[n,i] for n in range(nE) for i in range(nPlt)) >= (1+Setting.CRM_reserve)*quicksum(Enodes[n].demand[e_rep_hrs[t]] for n in range(nE)) for t in Te);


    # Ng system 
    GV.Xvpr = Model.addVars(nSVL,vtype=GRB.CONTINUOUS);
    GV.Xstr = Model.addVars(nSVL,vtype=GRB.CONTINUOUS);
    GV.Sstr = Model.addVars(nSVL,len(FY),vtype = GRB.CONTINUOUS);
    GV.Svpr = Model.addVars(nSVL,len(FY),vtype = GRB.CONTINUOUS);
    GV.Sliq = Model.addVars(nSVL,len(FY),vtype = GRB.CONTINUOUS);
    GV.supply = Model.addVars(nG,len(FY),vtype = GRB.CONTINUOUS);
    GV.Shed =  Model.addVars(nG,len(FY),vtype = GRB.CONTINUOUS);
    GV.RNG_use = Model.addVars(nG,len(FY),vtype = GRB.CONTINUOUS);
    GV.flowGG =  Model.addVars(nPipe,len(FY),vtype = GRB.CONTINUOUS);
    GV.flowGE =  Model.addVars(nG,nE,len(FY),vtype = GRB.CONTINUOUS);
    GV.flowGL =  Model.addVars(nG,nSVL,len(FY),vtype = GRB.CONTINUOUS);
    GV.flowVG =  Model.addVars(nSVL,nG,len(FY),vtype = GRB.CONTINUOUS);
    
    GV.inv_str_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.inv_pipe_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.pipe_FOM_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.pipe_Decom_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.shed_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.RNG_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.fom_str_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.import_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.emis_amount = Model.addVar(vtype=GRB.CONTINUOUS);
    GV.g_system_cost = Model.addVar(vtype=GRB.CONTINUOUS);
    
    # NG System Objective Function
    inv_pipe = LinExpr(quicksum(PipeLines[b].inv_coef*PipeLines[b].length*Other_input.pipe_per_mile*GV.Zg_val[b] for b in range(nPipe)));   
    ng_import = LinExpr(quicksum(GV.supply[k,tau]*Other_input.NG_price for k in range(nG) for tau in FY));
    storage_inv = LinExpr(quicksum(SVLs[0].inv_coef*(SVLs[0].capex*GV.Xstr[j]+SVLs[1].capex*GV.Xvpr[j]) for j in range(nSVL)));
    storage_FOM = LinExpr(quicksum(SVLs[1].FOM*(Exist_SVL[j].vap_cap+GV.Xvpr[j])+SVLs[0].FOM*(Exist_SVL[j].str_cap+GV.Xstr[j]) for j in range(nSVL)));
    ng_shed = LinExpr(quicksum(Setting.g_shed_penalty*GV.Shed[k,tau] for k in range(nG) for tau in FY));
    rng_import = LinExpr(quicksum(Other_input.RNG_price*GV.RNG_use[k,tau]for k in range(nG) for tau in FY));
    pipe_FOM = LinExpr(quicksum(PipeLines[b].FOM*PipeLines[b].length*GV.ZgOp_val[b] for b in range(nPipe)));  
    pipe_Dec = LinExpr(quicksum(PipeLines[b].decom*PipeLines[b].length*GV.ZgDec_val[b] for b in range(nPipe)));  
    
    ng_total_cost = inv_pipe+ng_import+storage_inv+storage_FOM+ng_shed+rng_import+pipe_FOM+pipe_Dec;
    Model.addConstr(GV.inv_pipe_cost==inv_pipe);
    Model.addConstr(GV.import_cost==ng_import);
    Model.addConstr(GV.inv_str_cost==storage_inv);
    Model.addConstr(GV.fom_str_cost==storage_FOM);
    Model.addConstr(GV.shed_cost==ng_shed);
    Model.addConstr(GV.RNG_cost==rng_import);
    Model.addConstr(GV.pipe_FOM_cost == pipe_FOM);
    Model.addConstr(GV.pipe_Decom_cost == pipe_Dec);
    Model.addConstr(GV.g_system_cost==ng_total_cost);
    
    # NG System Constraints
    # NG System Constraints
    #C1, C2: flow limit for NG
    Model.addConstrs(GV.flowGG[i,tau]<=PipeLines[i].Cap for i in range(nPipe) for tau in FY if PipeLines[i].is_exist==1);
    Model.addConstrs(GV.flowGG[i,tau]<=PipeLines[i].Cap*GV.Zg_val[i] for i in range(nPipe) for tau in FY if PipeLines[i].is_exist==0);
    
    # C3: flow balance, NG node (no shedding allowed as RNG is considered)
    GV.marginal_prices = Model.addConstrs(GV.supply[k,tau]
                      -quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_exp)
                      +quicksum(GV.flowGG[l,tau] for l in Gnodes[k].L_imp)
                      -quicksum(GV.flowGE[k,n,tau] for n in Gnodes[k].adjE)
                      +quicksum((GV.flowVG[j,k,tau]-GV.flowGL[k,j,tau]) for j in Gnodes[k].adjS)
                      +GV.Shed[k,tau]
                      +GV.RNG_use[k,tau]
                     == Gnodes[k].demand[tau] for k in range(nG) for tau in FY);                             
    
    # inforce all flowGE of the same cluster to take the equal values
    Model.addConstrs(GV.flowGE[k,n,tau]==GV.flowGE[k,n,days2Medoid[tau]] for k in range(nG) for n in Gnodes[k].adjE for tau in FY)
    
    # C3,C4: injection (supply) limit and curtailment limit
    Model.addConstrs(GV.supply[k,tau]<= Gnodes[k].injU for k in range(nG) for tau in FY);    
    Model.addConstrs(GV.Shed[k,tau]+GV.RNG_use[k,tau]<= Gnodes[k].demand[tau] for k in range(nG) for tau in FY);
    Model.addConstrs(GV.RNG_use[k,tau]<= Gnodes[k].injU for k in range(nG) for tau in FY);   

     
    # C5: storage balance (storage strats empty)
    Model.addConstrs(GV.Sstr[j,tau]==Exist_SVL[j].str_cap*0+GV.Sliq[j,tau]-GV.Svpr[j,tau]/SVLs[1].eff_disCh for j in range(nSVL) for tau in FY if tau==0);
    Model.addConstrs(GV.Sstr[j,tau]==(1-SVLs[0].BOG)*GV.Sstr[j,tau-1]+GV.Sliq[j,tau]-GV.Svpr[j,tau]/SVLs[1].eff_disCh for j in range(nSVL) for tau in FY if tau>0);
    
    # C6,8: calculate Sliq, Svpr
    for j in range(nSVL):
        for tau in FY:
            NG_adj = [];
            for k in range(nG): 
                for j2 in Gnodes[k].adjS:
                    if j2==j:
                        NG_adj.append(k);
                        
            Model.addConstr(GV.Sliq[j,tau]==quicksum(GV.flowGL[k,j,tau] for k in NG_adj));
            Model.addConstr(GV.Svpr[j,tau]==quicksum(GV.flowVG[j,k,tau] for k in NG_adj));
                
    # C6: Sliq limit
    Model.addConstrs(GV.Sliq[j,tau]<=Exist_SVL[j].liq_cap for j in range(nSVL) for tau in FY);
    
    # C9: Svpr limit
    Model.addConstrs(GV.Svpr[j,tau]<=Exist_SVL[j].vap_cap+GV.Xvpr[j] for j in range(nSVL) for tau in FY);
    
    # C10: Sstr limit
    Model.addConstrs(GV.Sstr[j,tau]<=Exist_SVL[j].str_cap+GV.Xstr[j] for j in range(nSVL) for tau in FY);
    
    
    # # C11: Operational Pipelines
    # Model.addConstrs(GV.ZgOp[b] == PipeLines[b].is_exist+GV.Zg[b]-GV.ZgDec[b] for b in range(nPipe));
    # Model.addConstrs(GV.ZgDec[b]==0 for b in range(nPipe) if PipeLines[b].is_exist==0);
   
    
# coupling constraints
    Coupling_constraints(Model); 
    
    
    
    # Solve the model
    Model.modelSense = GRB.MINIMIZE;
    #Model.setObjective(EV.e_system_cost);
    Model.setObjective(GV.g_system_cost+ EV.e_system_cost);
    Model.setParam('OutputFlag', 0);
    Model.setParam('MIPGap', Setting.solver_gap);
    Model.setParam('Timelimit', Setting.wall_clock_time_lim);
    Model.optimize();
    
    GV.marginal_prices_val = Model.getAttr('Pi',GV.marginal_prices);
    
    # calculate power systems NG fuel cost
    ng_cost = 0;
    for k in range(nG):
        for tau in FY:
            #if GV.marginal_prices_val[k,tau]>5.45:
                #print(f"k: {k}, tau: {tau}, LMP: {GV.marginal_prices_val[k,tau]}")
            for n in Gnodes[k].adjE:
                #ng_cost += g_time_weight[tau]*GV.flowGE_val[k,n,g_rep_days[tau]]*GV.marginal_prices_val[k,g_rep_days[tau]];
                ng_cost += GV.flowGE_val[k,n,tau]*Other_input.NG_price;

                # if GV.flowGE_val[k,n,tau]>0:
                #     print(f"GV.flowGE[{k},{n},{tau}] = {GV.flowGE_val[k,n,tau]}")
    #print(ng_cost);
    if Setting.emis_case!=1:
        EV.gas_fuel_cost_val = ng_cost;
    
#     ng_costp = 0;NG_units = ["ng","CT","CC","CC-CCS"];
    
#     for t in Te:
#         for i in range(nPlt):
#             for n in range(nE):
#                 if plant2sym[i] in NG_units:
#                     ng_costp += Other_input.NG_price*e_time_weight[t]*Plants[i].heat_rate*EV.prod_val[n,t,i];
#     print(ng_costp);
    
#     # test the coupling constraint 1
# #(quicksum(GV.flowGE[k,n,g_rep_days[tau]] for k in range(nG) if n in Gnodes[k].adjE)==quicksum(Plants[i].heat_rate*EV.prod[n,t,i] for t in range(tau*24,(tau+1)*24) for i in range(nPlt) if plant2sym[i] in NG_units) for n in range(nE) for tau in Tg);
#     for n in range(nE):
#         for tau in Tg:
#             ssr=0;
#             for t in range(tau*24,(tau+1)*24):
#                 for i in range(nPlt):
#                     if plant2sym[i] in NG_units:
#                         ssr += Plants[i].heat_rate*EV.prod_val[n,t,i];
#             ssl=0;
#             for k in range(nG):
#                 if n in Gnodes[k].adjE:
#                     ssl += GV.flowGE_val[k,n,g_rep_days[tau]];
             
#             print(f"lhs: {ssl}, rhs: {ssr}");
    
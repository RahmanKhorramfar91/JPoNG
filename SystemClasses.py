# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 22:51:02 2023

@author: Rahman Khorramfar
"""


class EV(): # treated as a static class
    
    # decision variables
    Xest = [];Xdec=[];Xop=[];    
    X = []; Xup=[]; Xdown = [];Ze=[];
    YeCD=[]; YeLev=[]; YeStr=[];
    theta=[];Shed=[];
    prod=[]; 
    eSch=[]; eSdis=[]; eSlev=[];
    eSrem = []; eSday=[];
    flowE=[]; kappa_capt = []; kappa_pipe = [];
    
    # cost variables
    est_cost=[]; est_trans_cost=[];decom_cost=[];
    FOM_cost=[]; startup_cost=[]; VOM_cost=[];
    nuc_fuel_cost=[]; gas_fuel_cost = [];
    shedding_cost=[];
    elec_storage1_cost=[];
    elec_storage2_cost=[];
    CSS_cost = []; trans_FOM_cost = [];
    e_system_cost=[];
    
    
    # other
    emis_amount=[];
    
    # values
    Xop_val=[]; Xest_val=[];Xdec_val=[];Xup_val = []; X_val=[];
    YeCD_val=[]; YeLev_val=[];    
    prod_val=[]; flowE_val=[];Shed_val=[]; eSlev = [];
    Ze_val = []; eSdis_val = []; eSch_val = [];eSrem_val=[];eSday = [];
    CSS_cost_val = []; trans_FOM_cost_val =[];
    
    
class GV(): # treated as a static class
    
    Xstr = []; Xvpr = []; Sstr = []; Svpr = []; Sliq = [];
    supply = []; Shed = []; RNG_use=[]; flowGG=[];flowGE=[]; flowGL=[];flowVG=[];
    Zg=[]; marginal_prices = [];ZgOp = []; ZgDec=[];
    
    inv_str_cost=[]; inv_pipe_cost=[];
    shed_cost=[]; RNG_cost=[];fom_str_cost=[];
    import_cost=[];g_system_cost=[];pipe_FOM_cost=[];pipe_Decom_cost=[];

    # other
    emis_amount=[];
    
    # values
    supply_val = []; Shed_val = [];RNG_use_val=[]; flowGE_val=[];
    Zg_val=[]; marginal_prices_val = [];ZgOp_val=[];ZgDec_val=[];
    
    inv_str_cost_val=[]; inv_pipe_cost_val=[];
    shed_cost_val=[]; RNG_cost_val=[];fom_str_cost_val=[];
    import_cost_val=[];g_system_cost_val=[];
    pipe_FOM_cost_val=[];pipe_Decom_cost_val=[];
    
    
    
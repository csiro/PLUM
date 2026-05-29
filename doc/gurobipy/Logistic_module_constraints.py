# -*- coding: utf-8 -*-
"""The instance.

This module represents an instance of the scheduling problem.
"""

# --------------------------------------------------------------------#
# Imports

import os
import sys
from gurobipy import *

# PYTHONPATH!
here = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(here, "..", "library")))

import intervals
import Logistic_module_variables as LMV
import Logistic_module_parameters as LMP
import instance as Instc


# ----------------------------------------------------------------------#
###########################Defining CONSTRAINTS for MIP#####################
# CONSTRAINT RELATED TO FIXED TIME IMPACTS.A ACTIVITIES STARTED BRFORE FIXED TIME. FOR EXAMPLE LOADING IS DONE BUT
# UNLOADING IS NOT SO THE LOACTAION IS OCCUPIED. THE SOAKING IS STARTED BUT NOT COMPLETED.
# loading of each c should start exactly once in time horizon : sum_{t}(ls_{ct}) = sum_{t}(uls_{ct}) = 1;
def fixed_time_impact_const(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points):
    door_locid = [idl for idl in inst.iid2loc
                  if inst.iid2loc[idl].get_type( ) == Instc.LocationType.DOOR
                  ]
    max_soak_end = 0

    for idc in range_for_component:
        comp_obj = inst.get_comp_obj(idc)
        # id for valid work location
        id_work_loc = [inst.loc2iid[l] for l in LMP.comp_work_loc[idc - 1]]

        # the work location should be occupied at the start of time horizon if the component is loaded there but not unloaded
        if (comp_obj.get_load( ) is not None
                and comp_obj.get_unload( ) is not None
                and LMP.comp_load_move_fstart[idc - 1][0] > -1
                and LMP.comp_unload_move_fstart[idc - 1][0] <= -1
        ):
            mdl.addConstr(VARS.get_loc_occupy_vars( )[id_work_loc[0], idc, no_of_fixed_points] +
                          VARS.get_unload_end_time_vars( )[idc, no_of_fixed_points] == 1,
                          name = 'work_loc_occup_Start_TH_between_Load_Unload' + str(idc))

        # the door should remain closed  during soaking. Therefore if soaking is started before fixed time and
        # finishing after fixed time then for any time point i with value t_i< soak_end_time, the door should be closed.
        # This is implemented by t_i>= soak_end_time*door_closed_i
        if (LMP.comp_soak_duration[idc - 1][0] > -1
                and LMP.comp_fixed_soak_end[idc - 1][0] >= LMP.get_fixed_time( )
        ):
            if max_soak_end < LMP.comp_fixed_soak_end[idc - 1][0]:
                max_soak_end = LMP.comp_fixed_soak_end[idc - 1][0]

    if max_soak_end > 0:
        sum_load_unload = sum(VARS.get_load_start_time_vars( )[idc, no_of_fixed_points]
                              for idc in range_for_component
                              if (idc, no_of_fixed_points) in VARS.get_load_start_time_vars( )
                              # and LMP.comp_move_est[idc - 1][0]<max_soak_end
                              ) \
                          + sum(VARS.get_unload_end_time_vars( )[idc, no_of_fixed_points]
                                for idc in range_for_component
                                if (idc, no_of_fixed_points) in VARS.get_unload_end_time_vars( )
                                # and LMP.comp_move_lct[idc - 1][0]<max_soak_end
                                )
        mdl.addConstr(max_soak_end * sum_load_unload <= VARS.get_time_point_vars( )[no_of_fixed_points],
                      name = 'No_load_unload_during_soak_' + str(no_of_fixed_points))


# ----------------------------------------------------------------------#
###########################Defining CONSTRAINTS for MIP#####################
# CONSTRAINT RELATED TO START TIME OF COMPONENT LOAD, UNLOAD RESPECTIVELY
# loading of each c should start exactly once in time horizon : sum_{t}(ls_{ct}) = sum_{t}(uls_{ct}) = 1;
def move_start_once_consts(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points):
    print(" no_time_points ", no_time_points)
    for idc in range_for_component:
        comp_obj = inst.get_comp_obj(idc)
        if comp_obj.get_load( ) is not None and LMP.comp_load_move_fstart[idc - 1][0] <= -1:
            sum_load_start = VARS.get_load_start_time_vars( ).sum(idc, '*')
            mdl.addConstr(sum_load_start == 1, name = 'loading_starts_once_' + str(idc))
        if comp_obj.get_unload( ) is not None and LMP.comp_unload_move_fstart[idc - 1][0] <= -1:
            sum_unload_end = VARS.get_unload_end_time_vars( ).sum(idc, '*')
            mdl.addConstr(sum_unload_end == 1, name = 'unloading_ends_once_' + str(idc))


# ----------------------------------------------------------------------#
# soaking start/end exactly once in time horizon : sum_{t}(ss_{ct})= sum_{t}(se_{ct}) = 1;
def soak_start_once_consts(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points):
    for idc in range_for_component:
        if LMP.comp_soak_duration[idc - 1][0] > -1 and LMP.comp_fixed_soak_end[idc - 1][0] <= -1:
            sum_soak_start = VARS.get_soak_start_time_vars( ).sum(idc, '*')
            mdl.addConstr(sum_soak_start == 1, name = 'soaking_starts_once_' + str(idc))


#######################################################ADDING OBJECTIVE#################################################
def add_objective(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points):
    sum_lb = 0
    i = 0
    sum_expr = 0
    print("LMP.TimeHorizon[0]", LMP.TimeHorizon[0])
    for i in range(max(0, no_of_fixed_points), no_time_points - 1):
        mdl.addConstr(VARS.get_slk_vars( )[i] >= (VARS.get_time_point_vars( )[i]
                                                  - sum(VARS.get_soak_start_time_vars( )[idc, i1]
                                                        * (max(LMP.comp_move_est[idc - 1][0] +
                                                           LMP.comp_load_move_duration[idc-1][0], inst.iid2comp[idc].get_forecast_soak_start()))
                                                        for (idc, i1) in VARS.get_soak_start_time_vars( )
                                                        if
                                                        i1 == i and VARS.get_soak_start_time_vars( )[idc, i1].ub > 0.5
                                                        )
                                                  - sum(VARS.get_load_start_time_vars( )[idc, i1]
                                                        * (LMP.comp_move_est[idc - 1][0])
                                                        for (idc, i1) in VARS.get_load_start_time_vars( )
                                                        if
                                                        i1 == i and VARS.get_load_start_time_vars( )[idc, i1].ub > 0.5
                                                        )
                                                  - sum(VARS.get_unload_end_time_vars( )[idc, i1]
                                                        * (LMP.comp_move_lct[idc - 1][1])
                                                        for (idc, i1) in VARS.get_unload_end_time_vars( )
                                                        if
                                                        i1 == i and VARS.get_unload_end_time_vars( )[idc, i1].ub > 0.5
                                                        )
                                                  )
                      )
        mdl.addConstr(VARS.get_slk_vars( )[i] >= -(VARS.get_time_point_vars( )[i]
                                                   #- sum(VARS.get_soak_start_time_vars( )[idc, i1]
                                                   #      * (max(LMP.comp_move_est[idc - 1][0] +
                                                   #        LMP.comp_load_move_duration[idc-1][0], inst.iid2comp[idc].get_forecast_soak_start()))
                                                   #      for (idc, i1) in VARS.get_soak_start_time_vars( )
                                                    #     if
                                                    #     i1 == i and VARS.get_soak_start_time_vars( )[idc, i1].ub > 0.5
                                                    #     )
                                                   - sum(VARS.get_load_start_time_vars( )[idc, i1]
                                                         * (LMP.comp_move_est[idc - 1][0])
                                                         for (idc, i1) in VARS.get_load_start_time_vars( )
                                                         if
                                                         i1 == i and VARS.get_load_start_time_vars( )[idc, i1].ub > 0.5
                                                         )
                                                   - sum(VARS.get_unload_end_time_vars( )[idc, i1]
                                                         * (LMP.comp_move_lct[idc - 1][1])
                                                         for (idc, i1) in VARS.get_unload_end_time_vars( )
                                                         if
                                                         i1 == i and VARS.get_unload_end_time_vars( )[idc, i1].ub > 0.5
                                                         )
                                                   )
                      )

    diff_from_bl_timms = sum(VARS.get_slk_vars( )[i] for i in range(max(0, no_of_fixed_points), no_time_points - 1))


    obj_expr = (diff_from_bl_timms )
    # mdl.addConstr(obj_expr<=3*sum_lb)
    mdl.setObjective(obj_expr, GRB.MINIMIZE)
    mdl.update( )


#######################################################CONSTRUCT CONSTRAINTS AND OBJECTIVE##############################
def construct_constraints_objective(mdl, inst, VARS, LMP):
    print(" VARS.get_no_of_components() ", VARS.get_no_of_components( ), " VARS.get_time_points_req() ",
          VARS.get_time_points_req( ))
    no_of_components = VARS.get_no_of_components( )
    no_time_points = VARS.get_time_points_req( )
    range_for_component = range(1, no_of_components + 1)
    fixed_time_points = VARS.get_fixed_time_points_index( )
    no_of_fixed_points = 0
    if fixed_time_points is not None:
        no_of_fixed_points = len(fixed_time_points)
    print(" no_of_fixed_points ", no_of_fixed_points, " fixed_time_points ", fixed_time_points)
    # create constraints
    fixed_time_impact_const(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ fixed time impact constraints constructed...")

    move_start_once_consts(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ move start once constraints constructed...")

    soak_start_once_consts(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ soak_start_once_consts constructed...")

    door_closed_during_soaking_consts(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ door_closed_during_soaking_consts constructed...")

    precedence_consts(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ precedence_consts constructed...")

    time_duration_consts(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ time_duration_consts constructed...")

    load_unload_duration_impact_constraints(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ load_unload_duration_impact_constraints constructed...")

    driver_avail_during_moves(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ driver_avail_during_moves constructed...")

    all_loc_occpd_of_move_path_const(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ all_loc_occpd_of_move_path_const constructed...")

    atmost_one_comp_at_loc_at_anytime_const(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ atmost_one_comp_at_loc_at_anytime_const constructed...")

    exactly_one_path_for_loading_unloading_const(mdl, inst, VARS, no_time_points, range_for_component,
                                                 no_of_fixed_points)
    print("INFO: ------ door_closed_during_soaking_consts constructed...")

    loc_avail_at_time_point_checking_const(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ loc_avail_at_time_point_checking_const constructed...")

    drvr_avail_at_time_point_checking_const(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ drvr_avail_at_time_point_checking_const constructed...")

    one_activity_at_one_time_point(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ------ one_activity_at_one_time_point constructed...")

    # slack_for_positive_deviation_in_finish_unld_move_from_lct_const(mdl, inst, VARS, no_time_points,
    # range_for_component, no_of_fixed_points)
    # print("INFO: ------ slack_for_positive_deviation_in_finish_unld_move_from_lct_const constructed...")

    # create objective
    add_objective(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points)
    print("INFO: ---add_objective done")

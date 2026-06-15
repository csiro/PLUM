# -*- coding: utf-8 -*-
"""Constraint definitions for logistic scheduling mixed integer program (MIP).

This module defines constraints for a scheduling problem in a metabolic gap-filling
context using flux balance analysis. It handles timing constraints for component
loading, unloading, soaking operations, and resource allocation (locations and drivers)
in a logistic scheduling optimization problem.

The module works with Gurobi optimization model to construct constraints that ensure:
- Proper sequencing of operations (load, soak, unload)
- Resource availability and non-overlapping usage
- Time window compliance for fixed-time activities
- Door closure during critical soaking periods
- Driver and location availability

Examples
--------
>>> import Logistic_module_constraints as LMC
>>> # After setting up model, instance, and variables
>>> LMC.construct_constraints_objective(mdl, inst, VARS, LMP)

Notes
-----
This module requires Gurobi optimizer and depends on:
- Logistic_module_variables (LMV) for variable definitions
- Logistic_module_parameters (LMP) for problem parameters
- instance module for component and location objects
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
    """Add constraints for activities with fixed-time impacts to the MIP model.

    Handles constraints for activities that started before a fixed time point but
    continue beyond it. This includes partially completed loading/unloading operations
    and soaking processes that span across the fixed time boundary.

    Parameters
    ----------
    mdl : gurobipy.Model
        The Gurobi optimization model to which constraints will be added.
    inst : instance.Instance
        Instance object containing component and location information.
    VARS : Logistic_module_variables.Variables
        Variable container with decision variables for the optimization model.
    no_time_points : int
        Total number of time points in the scheduling horizon.
    range_for_component : range
        Range object defining component indices to iterate over.
    no_of_fixed_points : int
        Number of fixed time points that have already occurred.

    Returns
    -------
    None
        Constraints are added directly to the model.

    Notes
    -----
    This function enforces two key constraint types:
    1. Work location occupancy when loading is complete but unloading is not
    2. Door closure during soaking periods that extend beyond the fixed time

    The function prevents loading/unloading operations during active soaking
    by constraining time point values based on maximum soaking end time.

    Examples
    --------
    >>> fixed_time_impact_const(mdl, inst, VARS, 20, range(1, 11), 5)
    """
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
    """Add constraints ensuring each component's load and unload operations start exactly once.

    Enforces that within the scheduling time horizon, each component that requires
    loading must have exactly one loading start event, and each component requiring
    unloading must have exactly one unloading end event.

    Parameters
    ----------
    mdl : gurobipy.Model
        The Gurobi optimization model to which constraints will be added.
    inst : instance.Instance
        Instance object containing component and location information.
    VARS : Logistic_module_variables.Variables
        Variable container with decision variables for the optimization model.
    no_time_points : int
        Total number of time points in the scheduling horizon.
    range_for_component : range
        Range object defining component indices to iterate over.
    no_of_fixed_points : int
        Number of fixed time points that have already occurred.

    Returns
    -------
    None
        Constraints are added directly to the model.

    Notes
    -----
    Only components with load/unload operations not already started (indicated by
    negative fixed start time values) receive these constraints. The constraint
    form is: sum_{t}(ls_{ct}) = 1 for loading and sum_{t}(uls_{ct}) = 1 for unloading.

    Examples
    --------
    >>> move_start_once_consts(mdl, inst, VARS, 20, range(1, 11), 5)
    """
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
    """Add constraints ensuring each component's soaking operation starts exactly once.

    Enforces that within the scheduling time horizon, each component requiring a
    soaking period must have exactly one soaking start event.

    Parameters
    ----------
    mdl : gurobipy.Model
        The Gurobi optimization model to which constraints will be added.
    inst : instance.Instance
        Instance object containing component and location information.
    VARS : Logistic_module_variables.Variables
        Variable container with decision variables for the optimization model.
    no_time_points : int
        Total number of time points in the scheduling horizon.
    range_for_component : range
        Range object defining component indices to iterate over.
    no_of_fixed_points : int
        Number of fixed time points that have already occurred.

    Returns
    -------
    None
        Constraints are added directly to the model.

    Notes
    -----
    Only components with positive soaking duration and no fixed soaking end time
    (indicated by negative fixed soak end values) receive these constraints.
    The constraint form is: sum_{t}(ss_{ct}) = 1.

    Examples
    --------
    >>> soak_start_once_consts(mdl, inst, VARS, 20, range(1, 11), 5)
    """
    for idc in range_for_component:
        if LMP.comp_soak_duration[idc - 1][0] > -1 and LMP.comp_fixed_soak_end[idc - 1][0] <= -1:
            sum_soak_start = VARS.get_soak_start_time_vars( ).sum(idc, '*')
            mdl.addConstr(sum_soak_start == 1, name = 'soaking_starts_once_' + str(idc))


#######################################################ADDING OBJECTIVE#################################################
def add_objective(mdl, inst, VARS, no_time_points, range_for_component, no_of_fixed_points):
    """Define and set the objective function for the MIP model.

    Creates slack variables and constraints to minimize deviations between scheduled
    time points and estimated/forecasted start times. The objective minimizes the
    sum of absolute deviations from baseline timings for load, unload, and soak
    operations.

    Parameters
    ----------
    mdl : gurobipy.Model
        The Gurobi optimization model to which the objective will be added.
    inst : instance.Instance
        Instance object containing component forecast information.
    VARS : Logistic_module_variables.Variables
        Variable container with decision variables including slack variables.
    no_time_points : int
        Total number of time points in the scheduling horizon.
    range_for_component : range
        Range object defining component indices to iterate over.
    no_of_fixed_points : int
        Number of fixed time points that have already occurred.

    Returns
    -------
    None
        The objective function is set directly on the model.

    Notes
    -----
    The function uses slack variables to capture both positive and negative
    deviations from target times. Two constraints per time point ensure slack
    variables represent the absolute value of deviation. The objective minimizes
    total slack across all time points from the fixed point onward.

    The deviation is computed as:
    deviation = time_point - (weighted sum of operation starts at that time point)

    where weights are earliest start times (EST) or latest completion times (LCT).

    Examples
    --------
    >>> add_objective(mdl, inst, VARS, 20, range(1, 11), 5)
    """
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
    """Construct all constraints and objective function for the logistic scheduling MIP.

    Master function that orchestrates the construction of the complete MIP model by
    sequentially adding all constraint types and the objective function. This includes
    timing constraints, resource constraints, precedence relationships, and the
    optimization objective.

    Parameters
    ----------
    mdl : gurobipy.Model
        The Gurobi optimization model to be fully specified.
    inst : instance.Instance
        Instance object containing all component and location data.
    VARS : Logistic_module_variables.Variables
        Variable container with all decision variables for the model.
    LMP : module
        Logistic_module_parameters module containing problem parameters.

    Returns
    -------
    None
        The model is modified in place with all constraints and objective.

    Notes
    -----
    This function coordinates the addition of the following constraint groups:
    - Fixed time impact constraints for partial operations
    - Move start-once constraints for load/unload operations
    - Soak start-once constraints for soaking operations
    - Door closure during soaking constraints
    - Precedence constraints for operation sequencing
    - Time duration constraints
    - Load/unload duration impact constraints
    - Driver availability during moves
    - Location occupancy for move paths
    - At-most-one component per location constraints
    - Path selection constraints for loading/unloading
    - Location and driver availability at time points
    - One activity per time point constraints

    Each constraint group addition is logged with an INFO message for debugging.
    Some constraint functions called may not be defined in this module.

    Examples
    --------
    >>> import Logistic_module_constraints as LMC
    >>> import Logistic_module_parameters as LMP
    >>> LMC.construct_constraints_objective(model, instance, variables, LMP)
    INFO: ------ fixed time impact constraints constructed...
    INFO: ------ move start once constraints constructed...
    ...
    """
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

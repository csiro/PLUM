# -*- coding: utf-8 -*-
"""Logistic scheduling variables module for mixed-integer programming optimization.

This module defines variables for a logistics scheduling problem using Gurobi optimization.
It manages time points, loading/unloading events, soak activities, and resource availability
for component scheduling in a logistic system.

The main class Variables encapsulates all decision variables required for the MIP model,
including time variables, event binary variables, and resource allocation variables.

Classes
-------
Variables
    Container for all decision variables in the logistics scheduling MIP model.

Functions
---------
craete_Variables
    Factory function to create and initialize a Variables instance.

Examples
--------
>>> mdl = gurobi.Model()
>>> inst = instance.Instance()
>>> vars_obj = craete_Variables(mdl, inst)
>>> time_points = vars_obj.get_time_points_req()

Notes
-----
This module is part of a larger logistics scheduling optimization system that uses
Gurobi for solving mixed-integer programming problems related to component loading,
unloading, and soak time management.
"""

# --------------------------------------------------------------------#
# Imports

import functools
import os
import sys
from operator import itemgetter
from enum import Enum
from gurobipy import *

# PYTHONPATH!
here = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(here, "..", "library")))

import intervals
from db_data import DataAccessError
import instance as Instc
import Logistic_module_parameters as LMP

# ----------------------------------------------------------------------#
time_points_req = []

fixed_time = 0


###########################CLASS defining all required variables for teh mdoel#####################
class Variables( ):
    """Container for all decision variables in the logistics scheduling MIP model.

    This class manages the creation, storage, and retrieval of all Gurobi decision variables
    required for the logistics scheduling optimization problem. It handles time point variables,
    binary event variables for loading/unloading/soak activities, and resource availability variables.

    Attributes
    ----------
    _time_points : int
        Total number of time points required for the scheduling problem.
    _time_vars : gurobi.tupledict
        Continuous time variables representing scheduling time points.
    _ls_vars : gurobi.tupledict
        Binary variables for component loading start events.
    _le_vars : gurobi.tupledict
        Binary variables for component loading end events.
    _uls_vars : gurobi.tupledict
        Binary variables for component unloading start events.
    _ule_vars : gurobi.tupledict
        Binary variables for component unloading end events.
    _ss_vars : gurobi.tupledict
        Binary variables for soak activity start events.
    _se_vars : gurobi.tupledict
        Binary variables for soak activity end events.
    _lmp_vars : gurobi.tupledict
        Binary variables for load move path selection.
    _ump_vars : gurobi.tupledict
        Binary variables for unload move path selection.
    _l_avail_tw_vars : gurobi.tupledict
        Binary variables for location availability within time windows.
    _l_avail_vars : gurobi.tupledict
        Binary variables for location availability.
    _d_avail_tw_vars : gurobi.tupledict
        Binary variables for driver availability within time windows.
    _d_avail_vars : gurobi.tupledict
        Binary variables for driver availability.
    _loc_occpy_vars : gurobi.tupledict
        Binary variables for location occupancy.
    _slk_unltime_vars : gurobi.tupledict
        Continuous slack variables for unload time constraints.
    _slk_wl_vars : gurobi.tupledict
        Continuous slack variables for workload constraints.
    _no_of_components : int
        Total number of components to be scheduled.
    _fixed_time_points : list of tuple
        List of fixed time points with format (component_id, time_value, event_type).
    _fixed_time_for_comp_load : list of list
        Fixed time indices for component loading events.
    _fixed_time_for_comp_unload : list of list
        Fixed time indices for component unloading events.
    _fixed_time_for_comp_soak_start : list of list
        Fixed time indices for component soak start events.
    _fixed_time_for_comp_soak_end : list of list
        Fixed time indices for component soak end events.
    _load_indexs : dict
        Dictionary mapping component IDs to valid loading time indices.
    _unload_indexs : dict
        Dictionary mapping component IDs to valid unloading time indices.
    _soak_start_indexs : dict
        Dictionary mapping component IDs to valid soak start time indices.
    _soak_end_indexs : dict
        Dictionary mapping component IDs to valid soak end time indices.

    Methods
    -------
    set_time_points_req(inst, range_for_component, fixed_time_points)
        Calculate the total number of required time points.
    set_time_vars(mdl, inst, range_for_component, fixed_time_points)
        Define time variables for the MIP model.
    set_move_start_end_time_vars(mdl, inst, range_for_component, fixed_time_points)
        Define binary variables for move start and end events.
    set_variables(mdl, inst, range_for_component, fixed_time_points)
        Initialize all required variables for the MIP model.
    get_time_points_req()
        Retrieve the total number of time points.
    fix_time_point_value(indx, val, mdl)
        Fix a specific time point variable to a given value.
    fix_load_start_event(compid, indx, mdl)
        Fix a loading start event for a component at a specific time index.

    Examples
    --------
    >>> mdl = gurobi.Model()
    >>> vars_obj = Variables(mdl)
    >>> vars_obj.set_variables(mdl, inst, range(1, 10), [])
    >>> time_vars = vars_obj.get_time_point_vars()
    """
    def __init__(self, mdl):
        """Initialize the Variables container with default values.

        Parameters
        ----------
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        All variable attributes are initialized to None or empty containers and must be
        populated by calling the appropriate set_* methods before use.
        """
        self._time_points = 0
        self._time_vars = None
        self._ls_vars = None
        self._le_vars = None
        self._uls_vars = None
        self._ule_vars = None
        self._ss_vars = None
        self._se_vars = None
        self._lmp_vars = None
        self._ump_vars = None
        self._l_avail_tw_vars = None
        self._l_avail_vars = None
        self._d_avail_tw_vars = None
        self._d_avail_vars = None
        self._loc_occpy_vars = None
        self._slk_unltime_vars = None
        self._slk_wl_vars = None
        self._no_of_components = 0
        self._fixed_time_points = None
        self._fixed_time_for_comp_load = []
        self._fixed_time_for_comp_unload = []
        self._fixed_time_for_comp_soak_start = []
        self._fixed_time_for_comp_soak_end = []
        self._load_indexs = {}
        self._unload_indexs = {}
        self._soak_start_indexs = {}
        self._soak_end_indexs = {}

    def _get(self, name):
        """Retrieve an attribute value with validation.

        Parameters
        ----------
        name : str
            The name of the attribute to retrieve.

        Returns
        -------
        object
            The value of the requested attribute.

        Raises
        ------
        DataAccessError
            If the attribute has not been set (value is None).

        Notes
        -----
        This is a private helper method used by getter methods to ensure variables
        have been initialized before access.
        """
        value = self.__getattribute__(name)
        if value is None:
            raise DataAccessError("Attempted access to " + str(name) + " before it was set!")
        return value

    def _get_fixed_events(self, mdl, inst, range_for_component):
        """Extract and sort fixed time events from component parameters.

        Parameters
        ----------
        mdl : gurobi.Model
            The Gurobi optimization model instance.
        inst : instance.Instance
            The problem instance containing component data.
        range_for_component : range or iterable
            Range or iterable of component IDs to process.

        Returns
        -------
        list of tuple
            Sorted list of tuples (component_id, time_value, event_type) where:
            - component_id : int - ID of the component
            - time_value : float - Fixed time value for the event
            - event_type : int - Type of event (0=load, 1=unload, 2=soak_start)
            Empty list if no fixed events exist.

        Notes
        -----
        Fixed events are extracted from the Logistic_module_parameters (LMP) and sorted
        by time value in ascending order. Event types are:
        - 0: component load move fixed start
        - 1: component unload move fixed start
        - 2: component fixed soak start
        """
        fixed_time_points = []
        for idc in range_for_component:
            if LMP.comp_load_move_fstart[idc - 1][0] > -1:
                fixed_time_points.append((idc, LMP.comp_load_move_fstart[idc - 1][0], 0))
            if LMP.comp_unload_move_fstart[idc - 1][0] > -1:
                fixed_time_points.append((idc, LMP.comp_unload_move_fstart[idc - 1][0], 1))
            if LMP.comp_fixed_soak_start[idc - 1][0] > -1:
                fixed_time_points.append((idc, LMP.comp_fixed_soak_start[idc - 1][0], 2))
        # if LMP.comp_fixed_soak_end[idc-1][0]>-1:
        # fixed_time_points.append((idc, LMP.comp_fixed_soak_end[idc-1][0],  3))
        if len(fixed_time_points) > 0:
            return sorted(fixed_time_points, key = lambda x: x[1])
        return fixed_time_points

    #####################################CALCULATE REQUIRED TIME POINTS######################################################################
    def set_time_points_req(self, inst, range_for_component, fixed_time_points):
        """Calculate the total number of required time points for the scheduling problem.

        Parameters
        ----------
        inst : instance.Instance
            The problem instance containing component data.
        range_for_component : range or iterable
            Range or iterable of component IDs to process.
        fixed_time_points : list of tuple
            List of fixed time points in format (component_id, time_value, event_type).

        Notes
        -----
        The method counts time points needed for:
        - Each component loading move (1 time point)
        - Each component unloading move (1 time point)
        - Each soak activity (1 time point for start)
        - One additional time point for the end of the time horizon

        The calculated value is stored in self._time_points.
        """
        for idc in range_for_component:
            comp_obj = inst.get_comp_obj(idc)
            # add a time points for each load move to be done
            if comp_obj.get_load( ) is not None:
                self._time_points = self._time_points + 1
            # add a time points for each unload move to be done
            if comp_obj.get_unload( ) is not None:
                self._time_points = self._time_points + 1
            # add 2 time points (one for start and one for end) for each soak activity to be done
            if LMP.comp_soak_duration[idc - 1][0] > -1:
                self._time_points = self._time_points + 1
        self._time_points = self._time_points + 1  # one time point corresponding to end of th

    ###########################Defining time variables for MIP#####################
    def set_time_vars(self, mdl, inst, range_for_component, fixed_time_points):
        """Define continuous time variables for the MIP model.

        Parameters
        ----------
        mdl : gurobi.Model
            The Gurobi optimization model instance.
        inst : instance.Instance
            The problem instance containing component data.
        range_for_component : range or iterable
            Range or iterable of component IDs to process.
        fixed_time_points : list of tuple
            List of fixed time points in format (component_id, time_value, event_type).

        Notes
        -----
        Creates Gurobi integer variables for each time point with:
        - Fixed bounds for predetermined time points (LB = UB = fixed_value)
        - Variable bounds for flexible time points (LB = fixed_time or time_horizon, UB = time_horizon)
        - The last time point has LB = time_horizon to represent the end of scheduling

        The created variables are stored in self._time_vars as a Gurobi tupledict.
        """
        indx = []
        LB = []
        UB = []
        indx_varhint = []
        indx_varhintPri = []
        for i in range(self._time_points):
            indx.append(i)
            if i < len(fixed_time_points):
                LB.append(fixed_time_points[i][1])
                UB.append(fixed_time_points[i][1])
            else:
                if i == self._time_points - 1:
                    LB.append(LMP.TimeHorizon[0])
                else:
                    LB.append(LMP.get_fixed_time( ))
                UB.append(LMP.TimeHorizon[0])
            indx_varhintPri.append(0)
            indx_varhint.append(0)
        self._time_vars = mdl.addVars(indx, lb = LB, ub = UB, VarHintVal = indx_varhint, VarHintPri = indx_varhintPri,
                                      vtype = GRB.INTEGER, name = "time_vars")

    ###########################Defining move time variables for MIP#####################
    def set_move_start_end_time_vars(self, mdl, inst, range_for_component, fixed_time_points):
        """Define binary variables for component loading and unloading events.

        Parameters
        ----------
        mdl : gurobi.Model
            The Gurobi optimization model instance.
        inst : instance.Instance
            The problem instance containing component data.
        range_for_component : range or iterable
            Range or iterable of component IDs to process.
        fixed_time_points : list of tuple
            List of fixed time points in format (component_id, time_value, event_type).

        Notes
        -----
        Creates binary variables for:
        - Loading start events (self._ls_vars): indexed by (component_id, time_index)
        - Unloading end events (self._ule_vars): indexed by (component_id, time_index)

        The method handles:
        - Fixed events: variables with LB=UB=1 for predetermined times
        - Flexible events: binary variables (0 or 1) with bounds adjusted based on:
          * Component requirements (load, unload, soak)
          * Valid time index ranges (preventing scheduling near horizon end)
          * Logical dependencies between events

        Valid time indices for each component are stored in:
        - self._load_indexs[component_id]: list of valid loading time indices
        - self._unload_indexs[component_id]: list of valid unloading time indices
        """
        indx_load_start = []
        LB_indx_load_start = []
        UB_indx_load_start = []
        LB_indx_unload_end = []
        UB_indx_unload_end = []
        indx_unload = []


        # add variables for components that requires loading
        for idc in range_for_component:
            comp_obj = inst.get_comp_obj(idc)
            if comp_obj.get_load( ) is not None:
                if idc not in self._load_indexs:
                    self._load_indexs[idc] = []
                if len(self._fixed_time_for_comp_load[idc - 1]) > 0:
                    indx_load_start.append((idc, self._fixed_time_for_comp_load[idc - 1][0]))
                    LB_indx_load_start.append(1.0)
                    UB_indx_load_start.append(1.0)
                    self._load_indexs[idc].append(self._fixed_time_for_comp_load[idc - 1][0])
                    continue
                for i in range(len(fixed_time_points), self._time_points):
                    self._load_indexs[idc].append(i)
                    indx_load_start.append((idc,i))  # iid2comp is dictionary of (id, component object) of insance that maps ids to component_object of instance

                    LB_indx_load_start.append(0.0)
                    if comp_obj.get_unload( ) is not None and LMP.comp_soak_duration[idc - 1][0] > -1:
                        if i >= self._time_points - 3:
                            UB_indx_load_start.append(0)
                        else:
                            UB_indx_load_start.append(1)
                    else:
                        if LMP.comp_soak_duration[idc - 1][0] > -1 or comp_obj.get_unload( ) is not None:
                            if i >= self._time_points - 2:
                                UB_indx_load_start.append(0)
                            else:
                                UB_indx_load_start.append(1)
                        else:
                            if i >= self._time_points - 1:
                                UB_indx_load_start.append(0)
                            else:
                                UB_indx_load_start.append(1)

        # add variables for components that requires unloading
        for idc in range_for_component:
            if comp_obj.get_unload( ) is not None:
                if idc not in self._unload_indexs:
                    self._unload_indexs[idc] = []
                if len(self._fixed_time_for_comp_unload[idc - 1]) > 0:
                    #print("unload fixed ", (idc, self._fixed_time_for_comp_unload[idc-1]))
                    indx_unload.append((idc, self._fixed_time_for_comp_unload[idc - 1][0]))
                    LB_indx_unload_end.append(1.0)
                    UB_indx_unload_end.append(1.0)
                    self._unload_indexs[idc].append(self._fixed_time_for_comp_unload[idc - 1][0])
                    continue
                for i in range(len(fixed_time_points), self._time_points):
                    #print("unload ", (idc, i))
                    indx_unload.append((idc, i))
                    self._unload_indexs[idc].append(i)
                    LB_indx_unload_end.append(0.0)
                    if comp_obj.get_load( ) is not None and LMP.comp_soak_duration[idc - 1][0] > -1:
                        if i <= 1 or i == self._time_points - 1:
                            UB_indx_unload_end.append(0)
                        else:
                            UB_indx_unload_end.append(1)
                    else:
                        if LMP.comp_soak_duration[idc - 1][0] > -1 or comp_obj.get_load( ) is not None:
                            if i < 1 or i == self._time_points - 1:
                                UB_indx_unload_end.append(0)
                            else:
                                UB_indx_unload_end.append(1)
                        else:
                            if i == self._time_points - 1:
                                UB_indx_unload_end.append(0)
                            else:
                                UB_indx_unload_end.append(1)

        self._ls_vars = mdl.addVars(indx_load_start, lb = LB_indx_load_start, ub = UB_indx_load_start, vtype = GRB.BINARY,
                                    name = "loading_start")
        self._ule_vars = mdl.addVars(indx_unload, lb = LB_indx_unload_end, ub = UB_indx_unload_end,
                                     vtype = GRB.BINARY, name = "unloading_end")


    ###########################SET Variables for MIP#####################
    def set_variables(self, mdl, inst, range_for_component, fixed_time_points):
        """Initialize all required variables for the MIP model.

        Parameters
        ----------
        mdl : gurobi.Model
            The Gurobi optimization model instance.
        inst : instance.Instance
            The problem instance containing component data.
        range_for_component : range or iterable
            Range or iterable of component IDs to process.
        fixed_time_points : list of tuple
            List of fixed time points in format (component_id, time_value, event_type).

        Notes
        -----
        This is the main entry point for variable creation. It calls in sequence:
        1. set_time_points_req() - calculate number of time points needed
        2. set_time_vars() - create time point variables
        3. set_move_start_end_time_vars() - create event binary variables
        """
        self.set_time_points_req(inst, range_for_component, fixed_time_points)
        self.set_time_vars(mdl, inst, range_for_component, fixed_time_points)
        self.set_move_start_end_time_vars(mdl, inst, range_for_component, fixed_time_points)
        
    ###########################Get Variables  defined for MIP#####################
    def get_time_points_req(self):
        """Retrieve the total number of required time points.

        Returns
        -------
        int
            The total number of time points in the scheduling problem.

        Raises
        ------
        DataAccessError
            If time points have not been calculated yet via set_time_points_req().
        """
        return self._get('_time_points')

    def get_time_point_vars(self):
        """Retrieve the time point variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of Gurobi integer variables indexed by time point index.

        Raises
        ------
        DataAccessError
            If time variables have not been created yet via set_time_vars().
        """
        return self._get('_time_vars')

    def get_load_start_time_vars(self):
        """Retrieve the loading start time binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables indexed by (component_id, time_index).

        Raises
        ------
        DataAccessError
            If loading start variables have not been created yet.
        """
        return self._get('_ls_vars')

    def get_load_end_time_vars(self):
        """Retrieve the loading end time binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables indexed by (component_id, time_index).

        Raises
        ------
        DataAccessError
            If loading end variables have not been created yet.
        """
        return self._get('_le_vars')

    def get_unload_start_time_vars(self):
        """Retrieve the unloading start time binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables indexed by (component_id, time_index).

        Raises
        ------
        DataAccessError
            If unloading start variables have not been created yet.
        """
        return self._get('_uls_vars')

    def get_unload_end_time_vars(self):
        """Retrieve the unloading end time binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables indexed by (component_id, time_index).

        Raises
        ------
        DataAccessError
            If unloading end variables have not been created yet.
        """
        return self._get('_ule_vars')

    def get_soak_start_time_vars(self):
        """Retrieve the soak start time binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables indexed by (component_id, time_index).

        Raises
        ------
        DataAccessError
            If soak start variables have not been created yet.
        """
        return self._get('_ss_vars')

    def get_soak_end_time_vars(self):
        """Retrieve the soak end time binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables indexed by (component_id, time_index).

        Raises
        ------
        DataAccessError
            If soak end variables have not been created yet.
        """
        return self._get('_se_vars')

    def get_load_path_vars(self):
        """Retrieve the load move path selection binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables indexed by (component_id, path_id).

        Raises
        ------
        DataAccessError
            If load path variables have not been created yet.
        """
        return self._get('_lmp_vars')

    def get_unload_path_vars(self):
        """Retrieve the unload move path selection binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables indexed by (component_id, path_id).

        Raises
        ------
        DataAccessError
            If unload path variables have not been created yet.
        """
        return self._get('_ump_vars')

    def get_loc_avail_tw_vars(self):
        """Retrieve the location availability within time window binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables for location availability in time windows.

        Raises
        ------
        DataAccessError
            If location availability time window variables have not been created yet.
        """
        return self._get('_l_avail_tw_vars')

    def get_loc_avail_vars(self):
        """Retrieve the location availability binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables for location availability.

        Raises
        ------
        DataAccessError
            If location availability variables have not been created yet.
        """
        return self._get('_l_avail_vars')

    def get_driver_avail_tw_vars(self):
        """Retrieve the driver availability within time window binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables for driver availability in time windows.

        Raises
        ------
        DataAccessError
            If driver availability time window variables have not been created yet.
        """
        return self._get('_d_avail_tw_vars')

    def get_driver_avail_vars(self):
        """Retrieve the driver availability binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables for driver availability.

        Raises
        ------
        DataAccessError
            If driver availability variables have not been created yet.
        """
        return self._get('_d_avail_vars')

    def get_loc_occupy_vars(self):
        """Retrieve the location occupancy binary variables.

        Returns
        -------
        gurobi.tupledict
            Dictionary of binary variables for location occupancy status.

        Raises
        ------
        DataAccessError
            If location occupancy variables have not been created yet.
        """
        return self._get('_loc_occpy_vars')

    def get_slk_vars(self):
        """Retrieve the slack variables for unload time constraints.

        Returns
        -------
        gurobi.tupledict
            Dictionary of continuous slack variables for relaxing unload time constraints.

        Raises
        ------
        DataAccessError
            If slack variables for unload time have not been created yet.
        """
        return self._get('_slk_unltime_vars')

    def get_slk_wl_vars(self):
        """Retrieve the slack variables for workload constraints.

        Returns
        -------
        gurobi.tupledict
            Dictionary of continuous slack variables for relaxing workload constraints.

        Raises
        ------
        DataAccessError
            If slack variables for workload have not been created yet.
        """
        return self._get('_slk_wl_vars')

    def get_no_of_components(self):
        """Retrieve the total number of components in the problem.

        Returns
        -------
        int
            The total number of components to be scheduled.

        Raises
        ------
        DataAccessError
            If the number of components has not been set yet.
        """
        return self._get('_no_of_components')

    def fix_time_point_value(self, indx, val, mdl):
        """Fix a time point variable to a specific value.

        Parameters
        ----------
        indx : int
            The index of the time point variable to fix.
        val : float
            The value to fix the time point to.
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        This method sets both lower and upper bounds to the specified value, effectively
        fixing the variable. The model is updated after the change.
        """
        self._time_vars[indx].lb = val
        self._time_vars[indx].ub = val
        mdl.update( )

    def fix_load_start_event(self, compid, indx, mdl):
        """Fix a loading start event for a component at a specific time index.

        Parameters
        ----------
        compid : int
            The component ID.
        indx : int
            The time index at which to fix the loading start event.
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        Sets the lower bound to 1, forcing the binary variable to 1 (event occurs).
        The model is updated after the change.
        """
        # print(compid, indx)
        self._ls_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_load_start_event ", compid, indx, self._ls_vars[compid, indx].lb)

    def fix_load_end_event(self, compid, indx, mdl):
        """Fix a loading end event for a component at a specific time index.

        Parameters
        ----------
        compid : int
            The component ID.
        indx : int
            The time index at which to fix the loading end event.
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        Sets the lower bound to 1, forcing the binary variable to 1 (event occurs).
        The model is updated after the change.
        """
        self._le_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_load_end_event ", compid, indx, self._le_vars[compid, indx].lb)

    def fix_unload_start_event(self, compid, indx, mdl):
        """Fix an unloading start event for a component at a specific time index.

        Parameters
        ----------
        compid : int
            The component ID.
        indx : int
            The time index at which to fix the unloading start event.
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        Sets the lower bound to 1, forcing the binary variable to 1 (event occurs).
        The model is updated after the change.
        """
        self._uls_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_unload_start_event ", compid, indx, self._uls_vars[compid, indx].lb)

    def fix_unload_end_event(self, compid, indx, mdl):
        """Fix an unloading end event for a component at a specific time index.

        Parameters
        ----------
        compid : int
            The component ID.
        indx : int
            The time index at which to fix the unloading end event.
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        Sets the lower bound to 1, forcing the binary variable to 1 (event occurs).
        The model is updated after the change.
        """
        self._ule_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_unload_end_event ", compid, indx, self._ule_vars[compid, indx].lb)

    def fix_soak_start_event(self, compid, indx, mdl):
        """Fix a soak start event for a component at a specific time index.

        Parameters
        ----------
        compid : int
            The component ID.
        indx : int
            The time index at which to fix the soak start event.
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        Sets the lower bound to 1, forcing the binary variable to 1 (event occurs).
        The model is updated after the change.
        """
        self._ss_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_soak_start_event ", compid, indx, self._ss_vars[compid, indx].lb)

    def fix_soak_end_event(self, compid, indx, mdl):
        """Fix a soak end event for a component at a specific time index.

        Parameters
        ----------
        compid : int
            The component ID.
        indx : int
            The time index at which to fix the soak end event.
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        Sets the lower bound to 1, forcing the binary variable to 1 (event occurs).
        The model is updated after the change.
        """
        self._se_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_soak_end_event ", compid, indx, self._se_vars[compid, indx].lb)

    def fix_load_path(self, compid, pathid, mdl):
        """Fix the load path selection for a component.

        Parameters
        ----------
        compid : int
            The component ID.
        pathid : int
            The path ID to fix for the loading move.
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        Sets the lower bound to 1, forcing the binary variable to 1 (path is selected).
        The model is updated after the change.
        """
        self._lmp_vars[compid, pathid].lb = 1
        mdl.update( )
        #print(" fix_load_path ", self._lmp_vars[compid, pathid].lb)

    def fix_unload_path(self, compid, pathid, mdl):
        """Fix the unload path selection for a component.

        Parameters
        ----------
        compid : int
            The component ID.
        pathid : int
            The path ID to fix for the unloading move.
        mdl : gurobi.Model
            The Gurobi optimization model instance.

        Notes
        -----
        Sets the lower bound to 1, forcing the binary variable to 1 (path is selected).
        The model is updated after the change.
        """
        self._ump_vars[compid, pathid].lb = 1
        mdl.update( )
        #print(" fix_unload_path ", self._ump_vars[compid, pathid].lb)

    def get_loading_index(self, compid, mdl):
        """Retrieve the valid loading time indices for a component.

        Parameters
        ----------
        compid : int
            The component ID.
        mdl : gurobi.Model
            The Gurobi optimization model instance (unused but kept for API consistency).

        Returns
        -------
        list of int
            List of valid time indices at which the component can be loaded.
        """
        return self._load_indexs[compid]

    def get_unloading_index(self, compid, mdl):
        """Retrieve the valid unloading time indices for a component.

        Parameters
        ----------
        compid : int
            The component ID.
        mdl : gurobi.Model
            The Gurobi optimization model instance (unused but kept for API consistency).

        Returns
        -------
        list of int
            List of valid time indices at which the component can be unloaded.
        """
        return self._unload_indexs[compid]

    def get_soak_start_index(self, compid, mdl):
        """Retrieve the valid soak start time indices for a component.

        Parameters
        ----------
        compid : int
            The component ID.
        mdl : gurobi.Model
            The Gurobi optimization model instance (unused but kept for API consistency).

        Returns
        -------
        list of int
            List of valid time indices at which the component soak can start.
        """
        return self._soak_start_indexs[compid]

    def get_soak_end_index(self, compid, mdl):
        """Retrieve the valid soak end time indices for a component.

        Parameters
        ----------
        compid : int
            The component ID.
        mdl : gurobi.Model
            The Gurobi optimization model instance (unused but kept for API consistency).

        Returns
        -------
        list of int
            List of valid time indices at which the component soak can end.
        """
        return self._soak_end_indexs[compid]

    def get_loading_fixed_index(self, compid, mdl):
        """Retrieve the fixed loading time index for a component.

        Parameters
        ----------
        compid : int
            The component ID (1-indexed).
        mdl : gurobi.Model
            The Gurobi optimization model instance (unused but kept for API consistency).

        Returns
        -------
        list of int
            List containing the fixed time index for loading, empty if not fixed.
        """
        return self._fixed_time_for_comp_load[compid - 1]

    def get_unloading_fixed_index(self, compid, mdl):
        """Retrieve the fixed unloading time index for a component.

        Parameters
        ----------
        compid : int
            The component ID (1-indexed).
        mdl : gurobi.Model
            The Gurobi optimization model instance (unused but kept for API consistency).

        Returns
        -------
        list of int
            List containing the fixed time index for unloading, empty if not fixed.

        Notes
        -----
        This method appears to reference _unload_fixed_indexs which may not exist in __init__.
        """
        return self._unload_fixed_indexs[compid - 1]

    def get_soak_start_fixed_index(self, compid, mdl):
        """Retrieve the fixed soak start time index for a component.

        Parameters
        ----------
        compid : int
            The component ID (1-indexed).
        mdl : gurobi.Model
            The Gurobi optimization model instance (unused but kept for API consistency).

        Returns
        -------
        list of int
            List containing the fixed time index for soak start, empty if not fixed.

        Notes
        -----
        This method returns _fixed_time_for_comp_unload which appears to be a mapping error.
        Should likely return _fixed_time_for_comp_soak_start.
        """
        return self._fixed_time_for_comp_unload[compid - 1]

    def get_soak_end_fixed_index(self, compid, mdl):
        """Retrieve the fixed soak end time index for a component.

        Parameters
        ----------
        compid : int
            The component ID (1-indexed).
        mdl : gurobi.Model
            The Gurobi optimization model instance (unused but kept for API consistency).

        Returns
        -------
        list of int
            List containing the fixed time index for soak end, empty if not fixed.

        Notes
        -----
        This method returns _fixed_time_for_comp_soak_start which appears to be a mapping error.
        Should likely return _fixed_time_for_comp_soak_end.
        """
        return self._fixed_time_for_comp_soak_start[compid - 1]

    def get_fixed_time_points_index(self):
        """Retrieve all fixed time points.

        Returns
        -------
        list of tuple
            List of tuples (component_id, time_value, event_type) representing all fixed events.
        """
        return self._fixed_time_points


###########################Create Variables####################
def craete_Variables(mdl, inst):
    """Create and initialize a Variables instance for the scheduling problem.

    This factory function creates a Variables object, extracts fixed events from component
    parameters, organizes them by component and event type, and initializes all decision
    variables required for the MIP model.

    Parameters
    ----------
    mdl : gurobi.Model
        The Gurobi optimization model instance.
    inst : instance.Instance
        The problem instance containing component data and scheduling requirements.

    Returns
    -------
    Variables
        Fully initialized Variables object with all decision variables created.

    Notes
    -----
    The function performs the following steps:
    1. Creates Variables object and sets the number of components
    2. Extracts fixed time points from LMP (Logistic_module_parameters)
    3. Organizes fixed time indices by component ID and event type:
       - Load events (event_type=0)
       - Unload events (event_type=1)
       - Soak start events (event_type=2)
       - Soak end events (event_type=3)
    4. Calls set_variables() to create all Gurobi decision variables
    5. Updates the model to register all new variables

    Examples
    --------
    >>> import gurobipy as gp
    >>> mdl = gp.Model('logistics')
    >>> inst = Instance()
    >>> vars_obj = craete_Variables(mdl, inst)
    >>> print(vars_obj.get_no_of_components())
    10
    """
    VARS = Variables(mdl)
    VARS._no_of_components = len(inst.iid2comp)
    print("VARS.get_no_of_components()+1 ", VARS.get_no_of_components( ) + 1)
    # no_time_points = VARS.get_time_points_req()
    range_of_comps = range(1, VARS.get_no_of_components( ) + 1)
    VARS._fixed_time_points = VARS._get_fixed_events(mdl, inst, range_of_comps)
    print("  VARS._fixed_time_points ", VARS._fixed_time_points, " length ", len(VARS._fixed_time_points))
    for idc in range_of_comps:
        VARS._fixed_time_for_comp_load.append([i for i in range(len(VARS._fixed_time_points))
                                                      if VARS._fixed_time_points[i][0] == idc
                                                         and VARS._fixed_time_points[i][2] == 0])
        VARS._fixed_time_for_comp_unload.append([i for i in range(len(VARS._fixed_time_points))
                                                        if VARS._fixed_time_points[i][0] == idc
                                                           and VARS._fixed_time_points[i][2] == 1])
        VARS._fixed_time_for_comp_soak_start.append([i for i in range(len(VARS._fixed_time_points))
                                                            if VARS._fixed_time_points[i][0] == idc
                                                               and VARS._fixed_time_points[i][2] == 2])
        VARS._fixed_time_for_comp_soak_end.append([i for i in range(len(VARS._fixed_time_points))
                                                          if VARS._fixed_time_points[i][0] == idc
                                                             and VARS._fixed_time_points[i][2] == 3])

    VARS.set_variables(mdl, inst, range_of_comps, VARS._fixed_time_points)
    mdl.update( )
    return VARS

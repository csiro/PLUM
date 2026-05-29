# -*- coding: utf-8 -*-
"""The instance.

This module represents an instance of the scheduling problem.
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
    def __init__(self, mdl):
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
        value = self.__getattribute__(name)
        if value is None:
            raise DataAccessError("Attempted access to " + str(name) + " before it was set!")
        return value

    def _get_fixed_events(self, mdl, inst, range_for_component):
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
        self.set_time_points_req(inst, range_for_component, fixed_time_points)
        self.set_time_vars(mdl, inst, range_for_component, fixed_time_points)
        self.set_move_start_end_time_vars(mdl, inst, range_for_component, fixed_time_points)
        
    ###########################Get Variables  defined for MIP#####################
    def get_time_points_req(self):
        return self._get('_time_points')

    def get_time_point_vars(self):
        return self._get('_time_vars')

    def get_load_start_time_vars(self):
        return self._get('_ls_vars')

    def get_load_end_time_vars(self):
        return self._get('_le_vars')

    def get_unload_start_time_vars(self):
        return self._get('_uls_vars')

    def get_unload_end_time_vars(self):
        return self._get('_ule_vars')

    def get_soak_start_time_vars(self):
        return self._get('_ss_vars')

    def get_soak_end_time_vars(self):
        return self._get('_se_vars')

    def get_load_path_vars(self):
        return self._get('_lmp_vars')

    def get_unload_path_vars(self):
        return self._get('_ump_vars')

    def get_loc_avail_tw_vars(self):
        return self._get('_l_avail_tw_vars')

    def get_loc_avail_vars(self):
        return self._get('_l_avail_vars')

    def get_driver_avail_tw_vars(self):
        return self._get('_d_avail_tw_vars')

    def get_driver_avail_vars(self):
        return self._get('_d_avail_vars')

    def get_loc_occupy_vars(self):
        return self._get('_loc_occpy_vars')

    def get_slk_vars(self):
        return self._get('_slk_unltime_vars')

    def get_slk_wl_vars(self):
        return self._get('_slk_wl_vars')

    def get_no_of_components(self):
        return self._get('_no_of_components')

    def fix_time_point_value(self, indx, val, mdl):
        self._time_vars[indx].lb = val
        self._time_vars[indx].ub = val
        mdl.update( )

    def fix_load_start_event(self, compid, indx, mdl):
        # print(compid, indx)
        self._ls_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_load_start_event ", compid, indx, self._ls_vars[compid, indx].lb)

    def fix_load_end_event(self, compid, indx, mdl):
        self._le_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_load_end_event ", compid, indx, self._le_vars[compid, indx].lb)

    def fix_unload_start_event(self, compid, indx, mdl):
        self._uls_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_unload_start_event ", compid, indx, self._uls_vars[compid, indx].lb)

    def fix_unload_end_event(self, compid, indx, mdl):
        self._ule_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_unload_end_event ", compid, indx, self._ule_vars[compid, indx].lb)

    def fix_soak_start_event(self, compid, indx, mdl):
        self._ss_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_soak_start_event ", compid, indx, self._ss_vars[compid, indx].lb)

    def fix_soak_end_event(self, compid, indx, mdl):
        self._se_vars[compid, indx].lb = 1
        mdl.update( )
        #print(" fix_soak_end_event ", compid, indx, self._se_vars[compid, indx].lb)

    def fix_load_path(self, compid, pathid, mdl):
        self._lmp_vars[compid, pathid].lb = 1
        mdl.update( )
        #print(" fix_load_path ", self._lmp_vars[compid, pathid].lb)

    def fix_unload_path(self, compid, pathid, mdl):
        self._ump_vars[compid, pathid].lb = 1
        mdl.update( )
        #print(" fix_unload_path ", self._ump_vars[compid, pathid].lb)

    def get_loading_index(self, compid, mdl):
        return self._load_indexs[compid]

    def get_unloading_index(self, compid, mdl):
        return self._unload_indexs[compid]

    def get_soak_start_index(self, compid, mdl):
        return self._soak_start_indexs[compid]

    def get_soak_end_index(self, compid, mdl):
        return self._soak_end_indexs[compid]

    def get_loading_fixed_index(self, compid, mdl):
        return self._fixed_time_for_comp_load[compid - 1]

    def get_unloading_fixed_index(self, compid, mdl):
        return self._unload_fixed_indexs[compid - 1]

    def get_soak_start_fixed_index(self, compid, mdl):
        return self._fixed_time_for_comp_unload[compid - 1]

    def get_soak_end_fixed_index(self, compid, mdl):
        return self._fixed_time_for_comp_soak_start[compid - 1]

    def get_fixed_time_points_index(self):
        return self._fixed_time_points


###########################Create Variables####################
def craete_Variables(mdl, inst):
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

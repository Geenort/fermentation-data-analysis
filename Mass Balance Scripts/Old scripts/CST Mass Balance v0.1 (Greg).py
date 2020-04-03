#Potential Errors:
#1) If you input a start or end time that does not match up with the interpolation 'fermentation age hours' times (ie, data exists for 0, 0.5, 0.1, 0.15 hours and you enter 0.12), shit will break.
#1) All round numbers exist in 'fermentation age hours'. If necessary, code could be modified to use a generic integer index, and you could screen for index points to use as start and end
#2) all blank cells are replaced with 0.0 to avoid 'nan' errors. if you try to use a time range that includes blank cells, average data will be all messed up
#3) Assumption: biomass loaded effluent liquid has density of 1000 g/L
#4) acetate not considered

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import math

#primary inputs
tank = 13
start_time = 400.0  #must be floats
end_time = 600.0
run_id = 'XYZ'
assumed_working_mass_g = 425.0

#input and output file strings
instr = 'C:/Users/gturk/desktop/Chemostat/run_set_data.csv'
outstr = 'C:/Users/gturk/desktop/Chemostat/' + run_id + '.csv'
varfile = 'C:/Users/gturk/desktop/Chemostat/variable_tag_list.csv'

#Assumptions/Utilities
conv_od_to_gpL = 0.42
T_ref_degC = 25.0
P_ref_barg = 0.0

water_c_h_o_n_ratio_mass = [0.0, 0.1111, 0.8889, 0.0]
biomass_c_h_o_n_ratio_mass = [0.5424, 0.0714, 0.1778, 0.2084]
bdo_c_h_o_n_ratio_mass = [0.5333, 0.1111, 0.355555, 0.0]

def list_avg(list):
    return sum(list)/len(list)

comps = ['argon', 'carbon dioxide', 'ethane', 'hydrogen', 'methane', 'nitrogen', 'oxygen', 'propane']
comp_mw = [40, 44, 30, 2, 16, 28, 16, 44]
c_dict = [0, 1, 2, 0, 1, 0, 0, 3]
h_dict = [0, 0, 6, 2, 4, 0, 0, 8]
o_dict = [0, 2, 0, 0, 0, 0, 2, 0]
n_dict = [0, 0, 0, 0, 0, 2, 0, 0]

variable_list = ['BDO, mg/L', 'OD, raw', 'Tank Scale, grams', 'F_in, sccm', 'x_argon_in_raw', 'x_carbon dioxide_in_raw', 'x_ethane_in_raw', 'x_hydrogen_in_raw', 'x_methane_in_raw', 'x_nitrogen_in_raw', 'x_oxygen_in_raw', 'x_propane_in_raw', 'x_argon_out_raw', 'x_carbon dioxide_out_raw', 'x_ethane_out_raw', 'x_hydrogen_out_raw', 'x_methane_out_raw', 'x_nitrogen_out_raw', 'x_oxygen_out_raw', 'x_propane_out_raw', ]
variable_tags = []

#set up results dataframe
results = pd.DataFrame(index=[], columns=[])
results.loc[run_id, 'Tank'] = tank
results.loc[run_id, 'Start Time, hours'] = start_time
results.loc[run_id, 'End Time, hours'] = end_time

#match variable list to header tags
var_map = pd.read_csv(varfile, index_col='Tank')
for x in range(0, len(variable_list)):
    variable_tags.append(var_map.at[tank, variable_list[x]])

#load data into a database
raw_data = pd.read_csv(instr)
#trim dataframe to only include data within the time window
raw_data = raw_data[raw_data['fermentation_age_hours'] >= start_time]
raw_data = raw_data[raw_data['fermentation_age_hours'] <= end_time]
#raw_data = raw_data[raw_data['offgas2.rms_flow'] >= 50.0]
raw_data = raw_data.set_index('fermentation_age_hours')

#remove any 'nan' values from dataset
for t in raw_data.index.tolist():
    for x in range(0, len(variable_tags)):
        if math.isnan(raw_data.at[t, variable_tags[x]]):
            raw_data.loc[t, variable_tags[x]] = 0.0

#Warn if steady state is not reached. steady state being +/- 2% of BDO, OD concentrations
bdo_start = raw_data.at[start_time, '2,3-butanediol']
bdo_end = raw_data.at[end_time, '2,3-butanediol']
od_start = raw_data.at[start_time, 'od']
od_end = raw_data.at[end_time, 'od']
if (abs(bdo_end - bdo_start) / bdo_end) > 0.02 or (abs(od_end - od_start) / od_end) > 0.02:
    print('WARNING: This time interval is not steady state!')
    print('Initial BDO in mg/L: ', bdo_start)
    print('Final BDO in mg/L ', bdo_end)
    print('Initial OD: ', od_start)
    print('Final OD: ', od_end )


#load average data into results
for x in range(0, len(variable_list)):
    var = variable_list[x]
    tag = variable_tags[x]
    #check units on tank scale, ensure in grams
    if var == 'Media Scale, grams':
        avg = list_avg(raw_data[tag].tolist())
        if avg < 1000.00:
            results.loc[run_id, var] = avg * 1000.0
        else:
            results.loc[run_id, var] = avg
    #for remainder of variables, just take average
    else:
        results.loc[run_id, var] = list_avg(raw_data[tag].tolist())

##Inlet gas
results.loc[run_id, 'Inlet Gas Rate, mol/hr'] = results.at[run_id, 'F_in, sccm'] * (1/1000) * 60 * (P_ref_barg + 1.0) / (0.08314 * (273 + T_ref_degC))
#normalize inlet gas composition, calculate MW
mass_spec_sum = 0.0
for x in range(0, len(comps)):
    tag = 'x_'+comps[x]+'_in_raw'
    mass_spec_sum += results.at[run_id, tag]
MW_sum = 0.0
for x in range(0, len(comps)):
    tag_raw = 'x_' + comps[x] + '_in_raw'
    tag_norm = 'x_'+comps[x]+'_in_normalized'
    results.loc[run_id, tag_norm] = results.at[run_id, tag_raw] * 100.00 / mass_spec_sum
    MW_sum += results.at[run_id, tag_norm] * comp_mw[x] / 100.0
results.loc[run_id, 'Inlet Gas MW'] = MW_sum
#calculate inlet gas component rates in mol/hr, elemental rates in mol/hr, elemental rates in grams/hr
in_gas_C_rate_mph = 0.0
in_gas_H_rate_mph = 0.0
in_gas_O_rate_mph = 0.0
in_gas_N_rate_mph = 0.0
for x in range(0, len(comps)):
    tag_norm = 'x_'+comps[x]+'_in_normalized'
    tag_rate = 'F_'+comps[x]+'_in_molph'
    results.loc[run_id, tag_rate] = results.at[run_id, 'Inlet Gas Rate, mol/hr'] * results.at[run_id, tag_norm] / 100.0
    in_gas_C_rate_mph += results.at[run_id, tag_rate] * c_dict[x]
    in_gas_H_rate_mph += results.at[run_id, tag_rate] * h_dict[x]
    in_gas_O_rate_mph += results.at[run_id, tag_rate] * o_dict[x]
    in_gas_N_rate_mph += results.at[run_id, tag_rate] * n_dict[x]
results.loc[run_id, 'Inlet C Rate, mol/hr'] = in_gas_C_rate_mph
results.loc[run_id, 'Inlet H Rate, mol/hr'] = in_gas_H_rate_mph
results.loc[run_id, 'Inlet O Rate, mol/hr'] = in_gas_O_rate_mph
results.loc[run_id, 'Inlet N Rate, mol/hr'] = in_gas_N_rate_mph

results.loc[run_id, 'Inlet C Rate, g/hr'] = in_gas_C_rate_mph * 12.0
results.loc[run_id, 'Inlet H Rate, g/hr'] = in_gas_H_rate_mph * 1.0
results.loc[run_id, 'Inlet O Rate, g/hr'] = in_gas_O_rate_mph * 16.0
results.loc[run_id, 'Inlet N Rate, g/hr'] = in_gas_N_rate_mph * 14.0

##Effluent Gas
#normalize inlet gas composition, calculate MW
mass_spec_sum = 0.0
for x in range(0, len(comps)):
    tag = 'x_'+comps[x]+'_out_raw'
    mass_spec_sum += results.at[run_id, tag]
MW_sum = 0.0
for x in range(0, len(comps)):
    tag_raw = 'x_' + comps[x] + '_out_raw'
    tag_norm = 'x_'+comps[x]+'_out_normalized'
    results.loc[run_id, tag_norm] = results.at[run_id, tag_raw] * 100.00 / mass_spec_sum
    MW_sum += results.at[run_id, tag_norm] * comp_mw[x] / 100.0
results.loc[run_id, 'Effluent Gas MW'] = MW_sum
#calculate effluent flow rate in mols/hr
results.loc[run_id, 'Effluent Gas Rate, mol/hr'] = results.at[run_id, 'Inlet Gas Rate, mol/hr'] * results.at[run_id, 'x_argon_in_normalized'] / results.at[run_id, 'x_argon_out_normalized']
#calculate effluent gas component rates in mol/hr, elemental rates in mol/hr, elemental rates in grams/hr
out_gas_C_rate_mph = 0.0
out_gas_H_rate_mph = 0.0
out_gas_O_rate_mph = 0.0
out_gas_N_rate_mph = 0.0
for x in range(0, len(comps)):
    tag_norm = 'x_'+comps[x]+'_out_normalized'
    tag_rate = 'F_'+comps[x]+'_out_molph'
    results.loc[run_id, tag_rate] = results.at[run_id, 'Effluent Gas Rate, mol/hr'] * results.at[run_id, tag_norm] / 100.0
    out_gas_C_rate_mph += results.at[run_id, tag_rate] * c_dict[x]
    out_gas_H_rate_mph += results.at[run_id, tag_rate] * h_dict[x]
    out_gas_O_rate_mph += results.at[run_id, tag_rate] * o_dict[x]
    out_gas_N_rate_mph += results.at[run_id, tag_rate] * n_dict[x]
results.loc[run_id, 'Effluent C Rate, mol/hr'] = out_gas_C_rate_mph
results.loc[run_id, 'Effluent H Rate, mol/hr'] = out_gas_H_rate_mph
results.loc[run_id, 'Effluent O Rate, mol/hr'] = out_gas_O_rate_mph
results.loc[run_id, 'Effluent N Rate, mol/hr'] = out_gas_N_rate_mph

results.loc[run_id, 'Effluent C Rate, g/hr'] = out_gas_C_rate_mph * 12.0
results.loc[run_id, 'Effluent H Rate, g/hr'] = out_gas_H_rate_mph * 1.0
results.loc[run_id, 'Effluent O Rate, g/hr'] = out_gas_O_rate_mph * 16.0
results.loc[run_id, 'Effluent N Rate, g/hr'] = out_gas_N_rate_mph * 14.0


##Liquid Feed
#liquid feed rate is d(feed scale)/dt. This mass balance is designed for a single dilution rate (feed rate / V)
#Step 1: check to see if rate at beginning and end of time period are the same
t_start = start_time
t_start_plus1 = start_time + 1.0
t_end = end_time
t_end_minus1 = end_time - 1.0

feed_rate_start = raw_data.at[t_start, 'scales_pumps.scale_feed'] - raw_data.at[t_start_plus1, 'scales_pumps.scale_feed']
feed_rate_end = raw_data.at[t_end_minus1, 'scales_pumps.scale_feed'] - raw_data.at[t_end, 'scales_pumps.scale_feed']

if ((feed_rate_end - feed_rate_start) / feed_rate_end) >= 0.05:
    print('WARNING: Dilution rates not consistent over feed interval')
# check units. should be either grams or kg. change to grams if given in kg
if feed_rate_end < 1.0:
    results.loc[run_id, 'Liquid Feed Rate, g/hr'] = feed_rate_end * 1000.00
else:
    results.loc[run_id, 'Liquid Feed Rate, g/hr'] = 0.5 * (feed_rate_start + feed_rate_end)
    results.loc[run_id, 'Dilution Rate, hrs^-1'] = results.at[run_id, 'Liquid Feed Rate, g/hr'] / assumed_working_mass_g

results.loc[run_id, 'Liquid C Rate, g/hr'] = results.at[run_id, 'Liquid Feed Rate, g/hr'] * water_c_h_o_n_ratio_mass[0]
results.loc[run_id, 'Liquid H Rate, g/hr'] = results.at[run_id, 'Liquid Feed Rate, g/hr'] * water_c_h_o_n_ratio_mass[1]
results.loc[run_id, 'Liquid O Rate, g/hr'] = results.at[run_id, 'Liquid Feed Rate, g/hr'] * water_c_h_o_n_ratio_mass[2]
results.loc[run_id, 'Liquid N Rate, g/hr'] = results.at[run_id, 'Liquid Feed Rate, g/hr'] * water_c_h_o_n_ratio_mass[3]

## Product Generation
# start with accumulation. Accumulation of C = Mass * (delta_w_bdo * sigma_O_bdo + delta_w_bm * sigma_O_bm + delta_w_h2o * sigma_O_h2o)
delta_w_bdo_gpL = (bdo_end - bdo_start) / 1000.0
delta_w_bm_gpL = (od_end - od_end) * conv_od_to_gpL
delta_w_h2o_gpL = -(delta_w_bdo_gpL + delta_w_bm_gpL)
results.loc[run_id, 'Accumulation C, g/hr'] = (assumed_working_mass_g / 1000.0 ) * ( delta_w_bdo_gpL * bdo_c_h_o_n_ratio_mass[0] + \
                                                                             delta_w_bm_gpL * biomass_c_h_o_n_ratio_mass[0] + \
                                                                             delta_w_h2o_gpL * water_c_h_o_n_ratio_mass[0] ) \
                                      / (end_time - start_time)
results.loc[run_id, 'Accumulation H, g/hr'] = (assumed_working_mass_g / 1000.0 ) * ( delta_w_bdo_gpL * bdo_c_h_o_n_ratio_mass[1] + \
                                                                             delta_w_bm_gpL * biomass_c_h_o_n_ratio_mass[1] + \
                                                                             delta_w_h2o_gpL * water_c_h_o_n_ratio_mass[1] ) \
                                      / (end_time - start_time)
results.loc[run_id, 'Accumulation O, g/hr'] = (assumed_working_mass_g / 1000.0 ) * ( delta_w_bdo_gpL * bdo_c_h_o_n_ratio_mass[2] + \
                                                                             delta_w_bm_gpL * biomass_c_h_o_n_ratio_mass[2] + \
                                                                             delta_w_h2o_gpL * water_c_h_o_n_ratio_mass[2] ) \
                                      / (end_time - start_time)
results.loc[run_id, 'Accumulation N, g/hr'] = (assumed_working_mass_g / 1000.0 ) * ( delta_w_bdo_gpL * bdo_c_h_o_n_ratio_mass[3] + \
                                                                             delta_w_bm_gpL * biomass_c_h_o_n_ratio_mass[3] + \
                                                                             delta_w_h2o_gpL * water_c_h_o_n_ratio_mass[3] ) \
                                      / (end_time - start_time)

# determine liquid mass fraction water
results.loc[run_id, 'w_BDO_liq_g/L'] = results.at[run_id, 'BDO, mg/L'] / 1000.0
results.loc[run_id, 'w_BM_liq_g/L'] = results.at[run_id, 'OD, raw'] * conv_od_to_gpL
results.loc[run_id, 'w_H2O_liq_g/L'] = 1000.0 - results.at[run_id, 'w_BDO_liq_g/L'] - results.at[run_id, 'w_BM_liq_g/L']
# Oxygen balance. In - Out = Acc, so m_in_O_gas + m_in_O_liq - m_acc_O = m_out_O_gas + m_out_O_liq
# m_out_O_liq = v_out_liq * (w_BDO_liq * w_O_BDO + w_BM_liq * w_O_BDO + w_H2O_liq * w_O_H2O)
# another option is to simply set 'Draw Rate, g/hr' = 'Liquid Feed Rate, g/hr'
results.loc[run_id, 'Draw Rate, L/hr'] = (results.at[run_id, 'Liquid O Rate, g/hr'] + results.at[run_id, 'Inlet O Rate, g/hr'] - results.at[run_id, 'Effluent O Rate, g/hr'] - results.at[run_id, 'Accumulation O, g/hr']) \
    / (results.at[run_id, 'w_BDO_liq_g/L'] * bdo_c_h_o_n_ratio_mass[2] + results.at[run_id, 'w_BM_liq_g/L'] * biomass_c_h_o_n_ratio_mass[2] + results.at[run_id, 'w_H2O_liq_g/L'] * water_c_h_o_n_ratio_mass[2])
results.loc[run_id, 'Draw Rate, g/hr'] = results.at[run_id, 'Draw Rate, L/hr'] *1000.0

# determine liquid draw elemental rates
results.loc[run_id, 'Draw C Rate, g/hr'] = results.at[run_id, 'Draw Rate, L/hr'] * (results.at[run_id, 'w_BDO_liq_g/L'] * bdo_c_h_o_n_ratio_mass[0] + \
                                                                                    results.at[run_id, 'w_BM_liq_g/L'] * biomass_c_h_o_n_ratio_mass[0] + \
                                                                                    results.at[run_id, 'w_H2O_liq_g/L'] * water_c_h_o_n_ratio_mass[0])
results.loc[run_id, 'Draw H Rate, g/hr'] = results.at[run_id, 'Draw Rate, L/hr'] * (results.at[run_id, 'w_BDO_liq_g/L'] * bdo_c_h_o_n_ratio_mass[1] + \
                                                                                    results.at[run_id, 'w_BM_liq_g/L'] * biomass_c_h_o_n_ratio_mass[1] + \
                                                                                    results.at[run_id, 'w_H2O_liq_g/L'] * water_c_h_o_n_ratio_mass[1])
results.loc[run_id, 'Draw O Rate, g/hr'] = results.at[run_id, 'Draw Rate, L/hr'] * (results.at[run_id, 'w_BDO_liq_g/L'] * bdo_c_h_o_n_ratio_mass[2] + \
                                                                                    results.at[run_id, 'w_BM_liq_g/L'] * biomass_c_h_o_n_ratio_mass[2] + \
                                                                                    results.at[run_id, 'w_H2O_liq_g/L'] * water_c_h_o_n_ratio_mass[2])
results.loc[run_id, 'Draw N Rate, g/hr'] = results.at[run_id, 'Draw Rate, L/hr'] * (results.at[run_id, 'w_BDO_liq_g/L'] * bdo_c_h_o_n_ratio_mass[3] + \
                                                                                    results.at[run_id, 'w_BM_liq_g/L'] * biomass_c_h_o_n_ratio_mass[3] + \
                                                                                    results.at[run_id, 'w_H2O_liq_g/L'] * water_c_h_o_n_ratio_mass[3])


##Calculate Closures
results.loc[run_id, 'C Closure'] = ( results.at[run_id, 'Effluent C Rate, g/hr'] + results.at[run_id, 'Draw C Rate, g/hr'] + results.at[run_id, 'Accumulation C, g/hr'] ) /\
                                   ( results.at[run_id, 'Inlet C Rate, g/hr'] + results.at[run_id, 'Liquid C Rate, g/hr'] )
results.loc[run_id, 'H Closure'] = ( results.at[run_id, 'Effluent H Rate, g/hr'] + results.at[run_id, 'Draw H Rate, g/hr'] + results.at[run_id, 'Accumulation H, g/hr'] ) /\
                                   ( results.at[run_id, 'Inlet H Rate, g/hr'] + results.at[run_id, 'Liquid H Rate, g/hr'] )
results.loc[run_id, 'O Closure'] = ( results.at[run_id, 'Effluent O Rate, g/hr'] + results.at[run_id, 'Draw O Rate, g/hr'] + results.at[run_id, 'Accumulation O, g/hr'] ) /\
                                   ( results.at[run_id, 'Inlet O Rate, g/hr'] + results.at[run_id, 'Liquid O Rate, g/hr'] )
results.loc[run_id, 'N Closure'] = ( results.at[run_id, 'Effluent N Rate, g/hr'] + results.at[run_id, 'Draw N Rate, g/hr'] + results.at[run_id, 'Accumulation N, g/hr'] ) /\
                                   ( results.at[run_id, 'Inlet N Rate, g/hr'] + results.at[run_id, 'Liquid N Rate, g/hr'] )


results.to_csv(outstr)



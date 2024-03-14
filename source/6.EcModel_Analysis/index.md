# 6.EcModel_Analysis
### Import related functions


```python
import cobra
import re 
import pandas as pd
import numpy as np
import plotly
import plotly.graph_objects as go
#from script.ECMpy_function import *
import sys
sys.path.append(r'./script/')
from ECMpy_function import *
```

##  Cumulative distribution figure

### Input and output files


```python
ecModel_file="./model/eciML1515.json"

method='AutoPACMEN'#DLKcat
reaction_kcat_MW_file = "./analysis/get_kcat_mw_by_%s/reaction_change_by_enzuse.csv"%method
reaction_kcat_MW = pd.read_csv(reaction_kcat_MW_file)
reaction_kcat_MW = round(reaction_kcat_MW,3)
reaction_kcat_dis_file='./analysis/reaction_kcat_distrbution.png'
reaction_mw_dis_file='./analysis/reaction_mw_distrbution.png'
```


```python
reaction_kcat_select = reaction_kcat_MW[reaction_kcat_MW['data_type'] != 'fill']
#Sort values
sorted_data = reaction_kcat_select.sort_values('kcat')
sorted_data = sorted_data.reset_index(drop=True)
y_index = sorted_data.index / (sorted_data.shape[0]  - 1)
data_cdf_data = sorted_data['kcat']
x_name="<b>kcat(1/s)<b>"
y_name="<b>Cummulative distribution<b>"
nticks=1000
fig=draw_cdf_fig(data_cdf_data,reaction_kcat_dis_file,x_name,y_name,y_index,nticks)
fig.show()
```


```python
reaction_kcat_select = reaction_kcat_MW[reaction_kcat_MW['data_type'] != 'fill']
#Sort values
sorted_data = reaction_kcat_select.sort_values('MW')
sorted_data = sorted_data.reset_index(drop=True)
y_index = sorted_data.index / (sorted_data.shape[0]  - 1)
data_cdf_data = sorted_data['MW']
y_index = sorted_data.index / (sorted_data.shape[0]  - 1)
data_cdf_data = data_cdf_data/1000# kDa
x_name="<b>mass(kDa)<b>"
y_name="<b>Cummulative distribution<b>"
nticks=10000
fig=draw_cdf_fig(data_cdf_data,reaction_mw_dis_file,x_name,y_name,y_index,nticks)
fig.show()
```

## Phenotype Phase Plane (PhPP) Analysis

### Input and output files


```python
ecModel_file="./model/eciML1515.json"
obj='BIOMASS_Ec_iML1515_core_75p37M'
z_id='EX_o2_e'
x_id='EX_glc__D_e'
substrate_bound=10
o2_bound=20

GEM_output_file='./analysis/iML1515_glc_o2_df.csv'
ecGEM_output_file='./analysis/eciML1515_glc_o2_df.csv'

PhPP_output_fig_file='./analysis/PhPP_combine.png'
```


```python
GEM_glc_o2_df=get_PhPP_data(ecModel_file, 'GEM', obj, substrate_bound,o2_bound,11, GEM_output_file,x_id,z_id)
GEM_glc_o2_df.drop(0,axis=0,inplace=True)
GEM_glc_o2_df.drop(0,axis=1,inplace=True)

ecGEM_glc_o2_df=get_PhPP_data(ecModel_file, 'ecGEM', obj, substrate_bound,o2_bound,11, GEM_output_file,x_id,z_id)
ecGEM_glc_o2_df.drop(0,axis=0,inplace=True)
ecGEM_glc_o2_df.drop(0,axis=1,inplace=True)
```


```python
fig=draw_3d_rbas(GEM_glc_o2_df,substrate_bound,o2_bound,0.9,11,PhPP_output_fig_file)
fig.show()
```


```python
fig=draw_3d_rbas(ecGEM_glc_o2_df,substrate_bound,o2_bound,0.9,11,PhPP_output_fig_file)
fig.show()
```

## Overflow simulation

### Input and output files


```python
# inputfiles
json_model_file="./model/eciML1515.json"
enz_model=get_enzyme_constraint_model(json_model_file)
norm_model=cobra.io.json.load_json_model(json_model_file)
method='AutoPACMEN'#DLKcat
reaction_kcat_MW_file = "./analysis/get_kcat_mw_by_%s/reaction_change_by_enzuse.csv"%method
reaction_kcat_MW = pd.read_csv(reaction_kcat_MW_file,index_col=0)
reaction_kcat_MW = round(reaction_kcat_MW,3)

# outputfiles
overflow_result_figfile="./analysis/pfba_overflow_result.png"
```


```python
use_substrate='EX_glc__D_e'
substrate_name='glucose'
glc_concentration_list = np.arange(1, 10, 0.5)
columns = ['biomass', 'acetate', 'O2', 'CO2', 'pyruvate', 'ethanol']
correspond_rxn = ['BIOMASS_Ec_iML1515_core_75p37M', 'EX_ac_e', 'EX_o2_e_reverse', 'EX_co2_e', 'EX_pyr_e', 'EX_etoh_e']
GEMyield_list = pd.DataFrame()
ecGEMyield_list = pd.DataFrame()

# Calculate yield for growth_model
with enz_model as growth_model:
    for glc_concentration in glc_concentration_list:
        ecGEMyield_list=calculate_yield(growth_model, use_substrate, substrate_name, glc_concentration,columns, correspond_rxn, ecGEMyield_list)

# Calculate yield for norm_model
with norm_model as growth_model:
    for glc_concentration in glc_concentration_list:
        GEMyield_list=calculate_yield(growth_model, use_substrate, substrate_name, glc_concentration,columns, correspond_rxn, GEMyield_list)

```


```python
substrate_name='glucose'
substrate_bound=10
obj_bound=1
secrate_bound=18
column_list=['biomass','CO2','ethanol','acetate']
y_axis_loc_list=['left','right','right','right']
color_list = generate_random_colors(len(column_list)*2)

fig=draw_overfolw_fig(GEMyield_list,ecGEMyield_list,column_list,y_axis_loc_list,color_list,substrate_name, substrate_bound,obj_bound,secrate_bound,overflow_result_figfile)
fig.show()
```

## Trade-off simulation


```python
# inputfiles
method='AutoPACMEN'#DLKcat
reaction_kcat_MW_file = "./analysis/get_kcat_mw_by_%s/reaction_change_by_enzuse.csv"%method
reaction_kcat_MW = pd.read_csv(reaction_kcat_MW_file,index_col=0)
reaction_kcat_MW = round(reaction_kcat_MW,3)
json_model_file="./model/eciML1515.json"
enz_model=get_enzyme_constraint_model(json_model_file)
glc_concentration_list = np.arange(1, 10, 0.5)
efficiency_file="./analysis/efficiency_pfba.csv"
trade_off_enzyme_efficiency_figfile="./analysis/trade_off_enzyme_efficiency.png"
use_substrate='EX_glc__D_e'
obj='BIOMASS_Ec_iML1515_core_75p37M'

yield_cost_efficiency_df=get_yield_cost_efficiency(enz_model,glc_concentration_list,use_substrate,obj,reaction_kcat_MW,efficiency_file)
yield_cost_efficiency_df.head()

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>biomass</th>
      <th>glucose_simu</th>
      <th>glucose_set</th>
      <th>biomass yield</th>
      <th>min enzyme cost</th>
      <th>enzyme efficiency</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1.0</th>
      <td>0.089537</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.497425</td>
      <td>0.062600</td>
      <td>1.430292</td>
    </tr>
    <tr>
      <th>1.5</th>
      <td>0.134305</td>
      <td>1.5</td>
      <td>1.5</td>
      <td>0.497425</td>
      <td>0.093900</td>
      <td>1.430292</td>
    </tr>
    <tr>
      <th>2.0</th>
      <td>0.179073</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>0.497425</td>
      <td>0.125200</td>
      <td>1.430292</td>
    </tr>
    <tr>
      <th>2.5</th>
      <td>0.223841</td>
      <td>2.5</td>
      <td>2.5</td>
      <td>0.497425</td>
      <td>0.156500</td>
      <td>1.430292</td>
    </tr>
    <tr>
      <th>3.0</th>
      <td>0.268610</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>0.497425</td>
      <td>0.187801</td>
      <td>1.430292</td>
    </tr>
  </tbody>
</table>
</div>




```python
yield_cost_efficiency_df = pd.read_csv(efficiency_file)
fig=draw_trade_off(yield_cost_efficiency_df,trade_off_enzyme_efficiency_figfile)
fig.show()
```

## Phenotype Simulation


```python
growth_exp_file = "./data/growth_exp.csv"
json_model_file="./model/eciML1515.json"
enz_model=get_enzyme_constraint_model(json_model_file)
growth_exp = pd.read_csv(growth_exp_file, index_col=0)

growth_rate_diff_substrate_file = "./analysis/enz_model_growth_pfba.csv"
ECMpy_diff_substate_result_figfile="./analysis/ECMpy_diff_substate_result.png"
diff_model_diff_substate_result_figfile="./analysis/diff_model_diff_substate_result.png"

```


```python
substrates = list(growth_exp.index)

#growth = pd.DataFrame()
for substrate in substrates:
    with enz_model as growth_model: 
        growth_model.reactions.get_by_id('EX_glc__D_e_reverse').bounds =(0.0, 0.0) 
        growth_model.reactions.get_by_id(substrate).bounds = (-10, 0.0)
        pfba_solution = cobra.flux_analysis.pfba(growth_model)
        growth_exp.loc[substrate, 'ECMpy_flux'] = pfba_solution.fluxes['BIOMASS_Ec_iML1515_core_75p37M']
        #growth_exp.loc[substrate, 'sub_flux'] = pfba_solution.fluxes[substrate]
        
for substrate in substrates:
    with norm_model as growth_model: 
        growth_model.reactions.get_by_id('EX_glc__D_e_reverse').bounds =(0.0, 0.0) 
        growth_model.reactions.get_by_id(substrate).bounds = (-10, 0.0)
        pfba_solution = cobra.flux_analysis.pfba(growth_model)
        growth_exp.loc[substrate, 'iML1515_flux'] = pfba_solution.fluxes['BIOMASS_Ec_iML1515_core_75p37M']
        #growth_exp.loc[substrate, 'sub_flux'] = pfba_solution.fluxes[substrate]      

growth_exp.to_csv(growth_rate_diff_substrate_file,index=False)
growth_exp.head(5)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>substrate</th>
      <th>EXP</th>
      <th>ECMpy_flux</th>
      <th>iML1515_flux</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>EX_acgam_e</th>
      <td>&lt;b&gt;N-Acetyl-D-glucosamine&lt;b&gt;</td>
      <td>0.61</td>
      <td>0.618304</td>
      <td>1.121990</td>
    </tr>
    <tr>
      <th>EX_ac_e</th>
      <td>&lt;b&gt;Acetate&lt;b&gt;</td>
      <td>0.29</td>
      <td>0.208818</td>
      <td>0.208818</td>
    </tr>
    <tr>
      <th>EX_akg_e</th>
      <td>&lt;b&gt;2-Oxoglutarate&lt;b&gt;</td>
      <td>0.24</td>
      <td>0.472914</td>
      <td>0.543203</td>
    </tr>
    <tr>
      <th>EX_ala__L_e</th>
      <td>&lt;b&gt;L-Alanine&lt;b&gt;</td>
      <td>0.24</td>
      <td>0.347477</td>
      <td>0.377704</td>
    </tr>
    <tr>
      <th>EX_fru_e</th>
      <td>&lt;b&gt;Fructose&lt;b&gt;</td>
      <td>0.54</td>
      <td>0.676556</td>
      <td>0.869773</td>
    </tr>
  </tbody>
</table>
</div>




```python
phenotypes = pd.read_csv(growth_rate_diff_substrate_file,index_col=0)
exp_col_namw = 'EXP'
sim_col_name = 'ECMpy_flux'
x_max = 1
y_max = 1
out_figfile="./analysis/ECMpy_diff_substate_result.png"
visulization_phenotype_results(phenotypes,exp_col_namw,sim_col_name,x_max,y_max,ECMpy_diff_substate_result_figfile)

```


```python
phenotypes = pd.read_csv(growth_rate_diff_substrate_file,index_col=0)

exp_col_name = 'EXP'
sim1_col_name = 'ECMpy_flux'
sim2_col_name = 'iML1515_flux'
sim1_name = '<b><i>eci<i>ML1515<b>'
sim2_name = '<b><i>i<i>ML1515<b>'

out_figfile="./analysis/diff_model_diff_substate_result.png"
compare_phenotype_results(phenotypes,exp_col_name,sim1_col_name,sim1_name,sim2_col_name,sim2_name,diff_model_diff_substate_result_figfile)

```

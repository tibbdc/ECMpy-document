# 4.Get ecModel using ECMpy
### Import related functions


```python
import cobra
import re 
#from script.ECMpy_function import *
import sys
sys.path.append(r'./script/')
from ECMpy_function import *
```

### Input and output files


```python
# reaction kcat_mw
sbml_path = "./data/iML1515R.xml"
taxonom_id=83333
method='AutoPACMEN'#DLKcat
reaction_kcat_MW_file = "./analysis/get_kcat_mw_by_%s/reaction_kcat_MW.csv"%method
#paxdb丰度数据
gene_abundance_colname='abundance'
gene_abundance_file='./data/gene_abundance.csv' # downolad from https://pax-db.org/download
#The enzyme mass fraction,such as 0.405
#f=calculate_f_v2(sbml_path, gene_abundance_file,gene_abundance_colname,taxonom_id)
f=0.405
#Initial parameters
ptot = 0.56 # The total protein fraction in cell.
sigma = 1 # The approximated saturation of enzyme.e.g.,0.5/1.
lowerbound = 0   # Lowerbound  of enzyme concentration constraint. 
upperbound = round(ptot * f * sigma, 3)#total enzyme
ecModel_output_file="./model/iML1515_irr_enz_constraint.json"
```

## Get ecModel and simulation


```python
#Get ecModel
trans_model2enz_json_model_split_isoenzyme(sbml_path, reaction_kcat_MW_file, f, ptot, sigma, lowerbound, upperbound, ecModel_output_file)

#ecModel Simulation
obj='BIOMASS_Ec_iML1515_core_75p37M'# CG_biomass_cgl_ATCC13032 EX_lys_L_e
fluxes_outfile = './analysis/ECMpy_solution_%s_pfba.csv'%obj
use_substrate='EX_glc__D_e'
concentration=10
enz_model=get_enzyme_constraint_model(ecModel_output_file)
enz_model.objective=obj

#change original substrate in model
[ori_obj_id,ori_substrate_id_list,ori_sub_concentration,ori_ATPM]=get_model_substrate_obj(enz_model)
for eachsubid in ori_substrate_id_list:
    if re.search('_reverse',eachsubid):
        r_id_new=eachsubid.split('_reverse')[0]
        enz_model.reactions.get_by_id(eachsubid).bounds = (0, 0) 
        enz_model.reactions.get_by_id(r_id_new).bounds = (0, 0)  
    else:
        r_id_new=eachsubid+'_reverse'
        enz_model.reactions.get_by_id(eachsubid).bounds = (0, 0) 
        enz_model.reactions.get_by_id(r_id_new).bounds = (0, 0) 
        
enz_model.reactions.get_by_id(use_substrate).bounds = (-concentration, 0)
enz_model.reactions.get_by_id(use_substrate+'_reverse').bounds = (0, 0)

enz_model_pfba_solution = cobra.flux_analysis.pfba(enz_model)
enz_model_pfba_solution = get_fluxes_detail_in_model(enz_model,enz_model_pfba_solution,fluxes_outfile,ecModel_output_file)
print(enz_model_pfba_solution.fluxes[obj])
```

    0.1736430506243408
    

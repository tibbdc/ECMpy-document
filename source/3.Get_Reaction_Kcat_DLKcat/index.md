## Import related functions


```python
import cobra
import datetime 
import pandas as pd
import subprocess
import re
# from script.ECMpy_function import *
import sys
sys.path.append(r'./script/')
from ECMpy_function import *
```

# Input and output files


```python
dlkcat_folder = "./analysis/get_kcat_mw_by_DLkcat/"
create_file(dlkcat_folder)

# input files
sbml_path = "./data/iML1515R.xml"

# output files
gene_subnum_path = "%sgene_subnum.csv"%dlkcat_folder
sub_description_path = '%sget_gene_subunitDescription.csv'%dlkcat_folder
inchikey_list_file='%sinchikey_list.csv'%dlkcat_folder
inchikey_list_smilesfile='%sinchikey_list_smiles.csv'%dlkcat_folder
comdf_file= '%scomdf.csv'%dlkcat_folder
DLouputdf_file = '%sDLoutput.tsv'%dlkcat_folder
metdf_outfile='%smetabolites_reactions_gpr_similes_prosequence_mass_dropna.csv'%dlkcat_folder
metabolites_reactions_gpr_file = '%smetabolites_reactions_gpr.csv'%dlkcat_folder
prodf_file = '%sprodf.csv'%dlkcat_folder
DLinput_file= '%sDLinput.tsv'%dlkcat_folder
DL_reaction_kact_mw_file='%sreaction_kcat_MW.csv'%dlkcat_folder
```

    Path exists


# Get reaction kcat_mw using DLKcat

## Step 0: read GEM


```python
# Step 0: read GEM
if re.search('\.xml',sbml_path):
    model = cobra.io.read_sbml_model(sbml_path)
elif re.search('\.json',sbml_path):
    model = cobra.io.json.load_json_model(sbml_path)
```

## Step 1: subunit number of each reaction


```python
starttime=datetime.datetime.now()
# Step 1: subunit number of each reaction
print("Starting to fetch subunit number of each enzyme")
get_gene_subunitDescription(sub_description_path,model)#Download from the UniProt API, run it once.
subbnumdf = get_subunit_number(sub_description_path,gene_subnum_path)
print("Calculation done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to fetch subunit number of each enzyme
    Start downloading from UniProt...
    [0;38;2;66;227;35m100.00%[0;38;2;186;189;250m|█████████████████████████████|[0;38;2;237;166;178m 0:00:00|1:28:45 [0;38;2;146;52;247m ETC: 07-06 00:02:23[0m[K
    Success downloading! :-)
    Calculation done!
    
    1:28:49.808229


## Step 2: convert metbolites bigg id to smiles 


```python
starttime=datetime.datetime.now()
# Step 2: convert metbolites bigg id to smiles 
print("Starting to convert metbolites bigg id to smiles...")
metdf_name = get_met_bigg_id(model)
inchkeydf = convert_bigg_met_to_inchikey(metdf_name['met'],inchikey_list_file)#from BIGG
# inchkeydf = pd.read_csv('./data/inchikey_list.csv')
smilesdf = convert_inchikey_to_smiles(inchkeydf,inchikey_list_smilesfile)#from pubchem
print("Converting done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to convert metbolites bigg id to smiles...
    Converting...
    [0;38;2;66;227;35m100.00%[0;38;2;190;231;233m|█████████████████████████████|[0;38;2;28;97;15m 0:00:00|1:09:21 [0;38;2;146;52;247m ETC: 07-06 10:52:44[0m[K
    [] try again later
    Fail secure!
    Converting done!
    
    2:01:59.345124


## Step 3: get protein sequence and mass in model 


```python
starttime=datetime.datetime.now()
# Step 3: get protein sequence and mass in model 
print("Starting to get protein sequence and mass in model...")
subbnumdf = pd.read_csv(gene_subnum_path)
prodf = get_model_protein_sequence_and_mass(model,subbnumdf,prodf_file)
print("Getting done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to get protein sequence and mass in model...
    Getting done!
    
    0:22:14.932906


## Step 4: split the substrate of reactions to match the gene


```python
starttime=datetime.datetime.now()
# Step 4: split the substrate of reactions to match the gene
print("Starting to split the substrate of reactions to match the gene...")
spdf = split_substrate_to_match_gene(model,metabolites_reactions_gpr_file)
print("Splitting done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to split the substrate of reactions to match the gene...
    Splitting done!
    
    0:00:10.310155


## Step 5: combine the reaction--substrate--gene--protein_sequnce--mass and formate DLKcat input file


```python
starttime=datetime.datetime.now()
# Step 5: combine the reaction--substrate--gene--protein_sequnce--mass and formate DLKcat input file
print("Starting to combine data...")
metdf_name = get_met_bigg_id(model)
smilesdf = pd.read_csv(inchikey_list_smilesfile)
spdf = pd.read_csv(metabolites_reactions_gpr_file)
prodf = pd.read_csv(prodf_file)
comdf = combine_reactions_simles_sequence(spdf,smilesdf,prodf,comdf_file)
DLinputdf = generate_DLKCAT_input(comdf,metdf_name,metdf_outfile,DLinput_file)
print("Combinning done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to combine data...
    DLKCAT input file generated
    Combinning done!
    
    0:00:03.239874


## Step 6: use DLKcat calculate kcat


```python
starttime=datetime.datetime.now()
# Step 6: use DLKcat calculate kcat
print("Starting to Use DLKcat calculate kcat...")
cmd_str = "python ./script/prediction_for_input.py ./analysis/get_kcat_mw_by_DLkcat/DLinput.tsv ./analysis/get_kcat_mw_by_DLkcat/DLoutput.tsv"
subprocess.run(cmd_str, shell=True)
print("DLKcat done!")
print()

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to Use DLKcat calculate kcat...
    ./analysis/get_kcat_mw_by_DLkcat/DLinput.tsv
    It's time to start the prediction!
    -----------------------------------


    [12:24:08] SMILES Parse Error: syntax error while parsing: None
    [12:24:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:09] SMILES Parse Error: syntax error while parsing: None
    [12:24:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:10] SMILES Parse Error: syntax error while parsing: None
    [12:24:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:12] SMILES Parse Error: syntax error while parsing: None
    [12:24:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:13] SMILES Parse Error: syntax error while parsing: None
    [12:24:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:14] SMILES Parse Error: syntax error while parsing: None
    [12:24:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:15] SMILES Parse Error: syntax error while parsing: None
    [12:24:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:16] SMILES Parse Error: syntax error while parsing: None
    [12:24:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:17] SMILES Parse Error: syntax error while parsing: None
    [12:24:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:18] SMILES Parse Error: syntax error while parsing: None
    [12:24:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:19] SMILES Parse Error: syntax error while parsing: None
    [12:24:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:20] SMILES Parse Error: syntax error while parsing: None
    [12:24:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:22] SMILES Parse Error: syntax error while parsing: None
    [12:24:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:23] SMILES Parse Error: syntax error while parsing: None
    [12:24:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:24] SMILES Parse Error: syntax error while parsing: None
    [12:24:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:26] SMILES Parse Error: syntax error while parsing: None
    [12:24:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:29] SMILES Parse Error: syntax error while parsing: None
    [12:24:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:30] SMILES Parse Error: syntax error while parsing: None
    [12:24:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:31] SMILES Parse Error: syntax error while parsing: None
    [12:24:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:32] SMILES Parse Error: syntax error while parsing: None
    [12:24:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:33] SMILES Parse Error: syntax error while parsing: None
    [12:24:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:35] SMILES Parse Error: syntax error while parsing: None
    [12:24:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:36] SMILES Parse Error: syntax error while parsing: None
    [12:24:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:37] SMILES Parse Error: syntax error while parsing: None
    [12:24:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:38] SMILES Parse Error: syntax error while parsing: None
    [12:24:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:39] SMILES Parse Error: syntax error while parsing: None
    [12:24:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:41] SMILES Parse Error: syntax error while parsing: None
    [12:24:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:42] SMILES Parse Error: syntax error while parsing: None
    [12:24:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:43] SMILES Parse Error: syntax error while parsing: None
    [12:24:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:45] SMILES Parse Error: syntax error while parsing: None
    [12:24:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:46] SMILES Parse Error: syntax error while parsing: None
    [12:24:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:47] SMILES Parse Error: syntax error while parsing: None
    [12:24:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:48] SMILES Parse Error: syntax error while parsing: None
    [12:24:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:49] SMILES Parse Error: syntax error while parsing: None
    [12:24:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:52] SMILES Parse Error: syntax error while parsing: None
    [12:24:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:53] SMILES Parse Error: syntax error while parsing: None
    [12:24:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:58] SMILES Parse Error: syntax error while parsing: None
    [12:24:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:24:59] SMILES Parse Error: syntax error while parsing: None
    [12:24:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:01] SMILES Parse Error: syntax error while parsing: None
    [12:25:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:02] SMILES Parse Error: syntax error while parsing: None
    [12:25:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:03] SMILES Parse Error: syntax error while parsing: None
    [12:25:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:04] SMILES Parse Error: syntax error while parsing: None
    [12:25:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:08] SMILES Parse Error: syntax error while parsing: None
    [12:25:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:09] SMILES Parse Error: syntax error while parsing: None
    [12:25:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:11] SMILES Parse Error: syntax error while parsing: None
    [12:25:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:12] SMILES Parse Error: syntax error while parsing: None
    [12:25:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:13] SMILES Parse Error: syntax error while parsing: None
    [12:25:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:14] SMILES Parse Error: syntax error while parsing: None
    [12:25:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:20] SMILES Parse Error: syntax error while parsing: None
    [12:25:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:21] SMILES Parse Error: syntax error while parsing: None
    [12:25:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:34] SMILES Parse Error: syntax error while parsing: None
    [12:25:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:38] SMILES Parse Error: syntax error while parsing: None
    [12:25:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:40] SMILES Parse Error: syntax error while parsing: None
    [12:25:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:25:41] SMILES Parse Error: syntax error while parsing: None
    [12:25:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:02] SMILES Parse Error: syntax error while parsing: None
    [12:26:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:08] SMILES Parse Error: syntax error while parsing: None
    [12:26:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:09] SMILES Parse Error: syntax error while parsing: None
    [12:26:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:10] SMILES Parse Error: syntax error while parsing: None
    [12:26:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:11] SMILES Parse Error: syntax error while parsing: None
    [12:26:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:12] SMILES Parse Error: syntax error while parsing: None
    [12:26:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:14] SMILES Parse Error: syntax error while parsing: None
    [12:26:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:15] SMILES Parse Error: syntax error while parsing: None
    [12:26:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:16] SMILES Parse Error: syntax error while parsing: None
    [12:26:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:17] SMILES Parse Error: syntax error while parsing: None
    [12:26:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:18] SMILES Parse Error: syntax error while parsing: None
    [12:26:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:22] SMILES Parse Error: syntax error while parsing: None
    [12:26:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:23] SMILES Parse Error: syntax error while parsing: None
    [12:26:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:24] SMILES Parse Error: syntax error while parsing: None
    [12:26:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:25] SMILES Parse Error: syntax error while parsing: None
    [12:26:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:27] SMILES Parse Error: syntax error while parsing: None
    [12:26:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:28] SMILES Parse Error: syntax error while parsing: None
    [12:26:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:29] SMILES Parse Error: syntax error while parsing: None
    [12:26:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:30] SMILES Parse Error: syntax error while parsing: None
    [12:26:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:31] SMILES Parse Error: syntax error while parsing: None
    [12:26:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:32] SMILES Parse Error: syntax error while parsing: None
    [12:26:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:33] SMILES Parse Error: syntax error while parsing: None
    [12:26:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:34] SMILES Parse Error: syntax error while parsing: None
    [12:26:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:37] SMILES Parse Error: syntax error while parsing: None
    [12:26:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:39] SMILES Parse Error: syntax error while parsing: None
    [12:26:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:40] SMILES Parse Error: syntax error while parsing: None
    [12:26:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:41] SMILES Parse Error: syntax error while parsing: None
    [12:26:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:42] SMILES Parse Error: syntax error while parsing: None
    [12:26:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:43] SMILES Parse Error: syntax error while parsing: None
    [12:26:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:44] SMILES Parse Error: syntax error while parsing: None
    [12:26:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:45] SMILES Parse Error: syntax error while parsing: None
    [12:26:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:46] SMILES Parse Error: syntax error while parsing: None
    [12:26:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:48] SMILES Parse Error: syntax error while parsing: None
    [12:26:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:50] SMILES Parse Error: syntax error while parsing: None
    [12:26:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:51] SMILES Parse Error: syntax error while parsing: None
    [12:26:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:52] SMILES Parse Error: syntax error while parsing: None
    [12:26:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:53] SMILES Parse Error: syntax error while parsing: None
    [12:26:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:54] SMILES Parse Error: syntax error while parsing: None
    [12:26:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:55] SMILES Parse Error: syntax error while parsing: None
    [12:26:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:56] SMILES Parse Error: syntax error while parsing: None
    [12:26:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:57] SMILES Parse Error: syntax error while parsing: None
    [12:26:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:58] SMILES Parse Error: syntax error while parsing: None
    [12:26:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:26:59] SMILES Parse Error: syntax error while parsing: None
    [12:26:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:00] SMILES Parse Error: syntax error while parsing: None
    [12:27:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:01] SMILES Parse Error: syntax error while parsing: None
    [12:27:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:02] SMILES Parse Error: syntax error while parsing: None
    [12:27:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:03] SMILES Parse Error: syntax error while parsing: None
    [12:27:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:06] SMILES Parse Error: syntax error while parsing: None
    [12:27:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:11] SMILES Parse Error: syntax error while parsing: None
    [12:27:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:12] SMILES Parse Error: syntax error while parsing: None
    [12:27:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:14] SMILES Parse Error: syntax error while parsing: None
    [12:27:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:15] SMILES Parse Error: syntax error while parsing: None
    [12:27:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:16] SMILES Parse Error: syntax error while parsing: None
    [12:27:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:17] SMILES Parse Error: syntax error while parsing: None
    [12:27:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:18] SMILES Parse Error: syntax error while parsing: None
    [12:27:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:23] SMILES Parse Error: syntax error while parsing: None
    [12:27:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:24] SMILES Parse Error: syntax error while parsing: None
    [12:27:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:26] SMILES Parse Error: syntax error while parsing: None
    [12:27:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:28] SMILES Parse Error: syntax error while parsing: None
    [12:27:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:29] SMILES Parse Error: syntax error while parsing: None
    [12:27:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:30] SMILES Parse Error: syntax error while parsing: None
    [12:27:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:31] SMILES Parse Error: syntax error while parsing: None
    [12:27:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:33] SMILES Parse Error: syntax error while parsing: None
    [12:27:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:34] SMILES Parse Error: syntax error while parsing: None
    [12:27:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:35] SMILES Parse Error: syntax error while parsing: None
    [12:27:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:37] SMILES Parse Error: syntax error while parsing: None
    [12:27:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:38] SMILES Parse Error: syntax error while parsing: None
    [12:27:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:39] SMILES Parse Error: syntax error while parsing: None
    [12:27:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:41] SMILES Parse Error: syntax error while parsing: None
    [12:27:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:48] SMILES Parse Error: syntax error while parsing: None
    [12:27:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:49] SMILES Parse Error: syntax error while parsing: None
    [12:27:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:50] SMILES Parse Error: syntax error while parsing: None
    [12:27:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:51] SMILES Parse Error: syntax error while parsing: None
    [12:27:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:52] SMILES Parse Error: syntax error while parsing: None
    [12:27:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:53] SMILES Parse Error: syntax error while parsing: None
    [12:27:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:54] SMILES Parse Error: syntax error while parsing: None
    [12:27:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:27:55] SMILES Parse Error: syntax error while parsing: None
    [12:27:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:00] SMILES Parse Error: syntax error while parsing: None
    [12:28:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:01] SMILES Parse Error: syntax error while parsing: None
    [12:28:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:03] SMILES Parse Error: syntax error while parsing: None
    [12:28:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:04] SMILES Parse Error: syntax error while parsing: None
    [12:28:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:05] SMILES Parse Error: syntax error while parsing: None
    [12:28:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:06] SMILES Parse Error: syntax error while parsing: None
    [12:28:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:07] SMILES Parse Error: syntax error while parsing: None
    [12:28:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:10] SMILES Parse Error: syntax error while parsing: None
    [12:28:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:11] SMILES Parse Error: syntax error while parsing: None
    [12:28:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:12] SMILES Parse Error: syntax error while parsing: None
    [12:28:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:13] SMILES Parse Error: syntax error while parsing: None
    [12:28:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:14] SMILES Parse Error: syntax error while parsing: None
    [12:28:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:15] SMILES Parse Error: syntax error while parsing: None
    [12:28:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:17] SMILES Parse Error: syntax error while parsing: None
    [12:28:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:18] SMILES Parse Error: syntax error while parsing: None
    [12:28:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:19] SMILES Parse Error: syntax error while parsing: None
    [12:28:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:20] SMILES Parse Error: syntax error while parsing: None
    [12:28:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:21] SMILES Parse Error: syntax error while parsing: None
    [12:28:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:22] SMILES Parse Error: syntax error while parsing: None
    [12:28:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:23] SMILES Parse Error: syntax error while parsing: None
    [12:28:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:25] SMILES Parse Error: syntax error while parsing: None
    [12:28:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:26] SMILES Parse Error: syntax error while parsing: None
    [12:28:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:29] SMILES Parse Error: syntax error while parsing: None
    [12:28:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:31] SMILES Parse Error: syntax error while parsing: None
    [12:28:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:32] SMILES Parse Error: syntax error while parsing: None
    [12:28:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:33] SMILES Parse Error: syntax error while parsing: None
    [12:28:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:35] SMILES Parse Error: syntax error while parsing: None
    [12:28:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:36] SMILES Parse Error: syntax error while parsing: None
    [12:28:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:37] SMILES Parse Error: syntax error while parsing: None
    [12:28:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:38] SMILES Parse Error: syntax error while parsing: None
    [12:28:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:42] SMILES Parse Error: syntax error while parsing: None
    [12:28:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:44] SMILES Parse Error: syntax error while parsing: None
    [12:28:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:45] SMILES Parse Error: syntax error while parsing: None
    [12:28:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:46] SMILES Parse Error: syntax error while parsing: None
    [12:28:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:48] SMILES Parse Error: syntax error while parsing: None
    [12:28:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:49] SMILES Parse Error: syntax error while parsing: None
    [12:28:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:51] SMILES Parse Error: syntax error while parsing: None
    [12:28:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:52] SMILES Parse Error: syntax error while parsing: None
    [12:28:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:53] SMILES Parse Error: syntax error while parsing: None
    [12:28:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:55] SMILES Parse Error: syntax error while parsing: None
    [12:28:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:58] SMILES Parse Error: syntax error while parsing: None
    [12:28:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:28:59] SMILES Parse Error: syntax error while parsing: None
    [12:28:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:00] SMILES Parse Error: syntax error while parsing: None
    [12:29:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:01] SMILES Parse Error: syntax error while parsing: None
    [12:29:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:02] SMILES Parse Error: syntax error while parsing: None
    [12:29:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:04] SMILES Parse Error: syntax error while parsing: None
    [12:29:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:05] SMILES Parse Error: syntax error while parsing: None
    [12:29:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:06] SMILES Parse Error: syntax error while parsing: None
    [12:29:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:07] SMILES Parse Error: syntax error while parsing: None
    [12:29:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:08] SMILES Parse Error: syntax error while parsing: None
    [12:29:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:09] SMILES Parse Error: syntax error while parsing: None
    [12:29:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:10] SMILES Parse Error: syntax error while parsing: None
    [12:29:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:11] SMILES Parse Error: syntax error while parsing: None
    [12:29:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:12] SMILES Parse Error: syntax error while parsing: None
    [12:29:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:14] SMILES Parse Error: syntax error while parsing: None
    [12:29:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:16] SMILES Parse Error: syntax error while parsing: None
    [12:29:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:18] SMILES Parse Error: syntax error while parsing: None
    [12:29:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:19] SMILES Parse Error: syntax error while parsing: None
    [12:29:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:20] SMILES Parse Error: syntax error while parsing: None
    [12:29:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:23] SMILES Parse Error: syntax error while parsing: None
    [12:29:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:24] SMILES Parse Error: syntax error while parsing: None
    [12:29:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:25] SMILES Parse Error: syntax error while parsing: None
    [12:29:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:26] SMILES Parse Error: syntax error while parsing: None
    [12:29:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:27] SMILES Parse Error: syntax error while parsing: None
    [12:29:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:28] SMILES Parse Error: syntax error while parsing: None
    [12:29:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:29] SMILES Parse Error: syntax error while parsing: None
    [12:29:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:30] SMILES Parse Error: syntax error while parsing: None
    [12:29:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:31] SMILES Parse Error: syntax error while parsing: None
    [12:29:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:32] SMILES Parse Error: syntax error while parsing: None
    [12:29:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:33] SMILES Parse Error: syntax error while parsing: None
    [12:29:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:34] SMILES Parse Error: syntax error while parsing: None
    [12:29:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:35] SMILES Parse Error: syntax error while parsing: None
    [12:29:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:38] SMILES Parse Error: syntax error while parsing: None
    [12:29:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:39] SMILES Parse Error: syntax error while parsing: None
    [12:29:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:40] SMILES Parse Error: syntax error while parsing: None
    [12:29:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:41] SMILES Parse Error: syntax error while parsing: None
    [12:29:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:42] SMILES Parse Error: syntax error while parsing: None
    [12:29:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:43] SMILES Parse Error: syntax error while parsing: None
    [12:29:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:44] SMILES Parse Error: syntax error while parsing: None
    [12:29:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:45] SMILES Parse Error: syntax error while parsing: None
    [12:29:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:46] SMILES Parse Error: syntax error while parsing: None
    [12:29:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:47] SMILES Parse Error: syntax error while parsing: None
    [12:29:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:48] SMILES Parse Error: syntax error while parsing: None
    [12:29:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:49] SMILES Parse Error: syntax error while parsing: None
    [12:29:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:50] SMILES Parse Error: syntax error while parsing: None
    [12:29:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:52] SMILES Parse Error: syntax error while parsing: None
    [12:29:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:53] SMILES Parse Error: syntax error while parsing: None
    [12:29:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:54] SMILES Parse Error: syntax error while parsing: None
    [12:29:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:55] SMILES Parse Error: syntax error while parsing: None
    [12:29:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:56] SMILES Parse Error: syntax error while parsing: None
    [12:29:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:29:59] SMILES Parse Error: syntax error while parsing: None
    [12:29:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:01] SMILES Parse Error: syntax error while parsing: None
    [12:30:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:02] SMILES Parse Error: syntax error while parsing: None
    [12:30:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:03] SMILES Parse Error: syntax error while parsing: None
    [12:30:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:04] SMILES Parse Error: syntax error while parsing: None
    [12:30:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:05] SMILES Parse Error: syntax error while parsing: None
    [12:30:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:08] SMILES Parse Error: syntax error while parsing: None
    [12:30:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:09] SMILES Parse Error: syntax error while parsing: None
    [12:30:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:10] SMILES Parse Error: syntax error while parsing: None
    [12:30:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:11] SMILES Parse Error: syntax error while parsing: None
    [12:30:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:12] SMILES Parse Error: syntax error while parsing: None
    [12:30:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:15] SMILES Parse Error: syntax error while parsing: None
    [12:30:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:16] SMILES Parse Error: syntax error while parsing: None
    [12:30:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:17] SMILES Parse Error: syntax error while parsing: None
    [12:30:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:18] SMILES Parse Error: syntax error while parsing: None
    [12:30:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:19] SMILES Parse Error: syntax error while parsing: None
    [12:30:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:22] SMILES Parse Error: syntax error while parsing: None
    [12:30:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:23] SMILES Parse Error: syntax error while parsing: None
    [12:30:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:24] SMILES Parse Error: syntax error while parsing: None
    [12:30:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:25] SMILES Parse Error: syntax error while parsing: None
    [12:30:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:28] SMILES Parse Error: syntax error while parsing: None
    [12:30:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:29] SMILES Parse Error: syntax error while parsing: None
    [12:30:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:30] SMILES Parse Error: syntax error while parsing: None
    [12:30:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:31] SMILES Parse Error: syntax error while parsing: None
    [12:30:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:32] SMILES Parse Error: syntax error while parsing: None
    [12:30:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:33] SMILES Parse Error: syntax error while parsing: None
    [12:30:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:34] SMILES Parse Error: syntax error while parsing: None
    [12:30:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:35] SMILES Parse Error: syntax error while parsing: None
    [12:30:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:36] SMILES Parse Error: syntax error while parsing: None
    [12:30:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:37] SMILES Parse Error: syntax error while parsing: None
    [12:30:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:38] SMILES Parse Error: syntax error while parsing: None
    [12:30:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:39] SMILES Parse Error: syntax error while parsing: None
    [12:30:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:40] SMILES Parse Error: syntax error while parsing: None
    [12:30:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:41] SMILES Parse Error: syntax error while parsing: None
    [12:30:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:44] SMILES Parse Error: syntax error while parsing: None
    [12:30:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:45] SMILES Parse Error: syntax error while parsing: None
    [12:30:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:46] SMILES Parse Error: syntax error while parsing: None
    [12:30:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:47] SMILES Parse Error: syntax error while parsing: None
    [12:30:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:48] SMILES Parse Error: syntax error while parsing: None
    [12:30:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:49] SMILES Parse Error: syntax error while parsing: None
    [12:30:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:51] SMILES Parse Error: syntax error while parsing: None
    [12:30:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:52] SMILES Parse Error: syntax error while parsing: None
    [12:30:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:53] SMILES Parse Error: syntax error while parsing: None
    [12:30:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:55] SMILES Parse Error: syntax error while parsing: None
    [12:30:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:56] SMILES Parse Error: syntax error while parsing: None
    [12:30:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:58] SMILES Parse Error: syntax error while parsing: None
    [12:30:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:30:59] SMILES Parse Error: syntax error while parsing: None
    [12:30:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:00] SMILES Parse Error: syntax error while parsing: None
    [12:31:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:01] SMILES Parse Error: syntax error while parsing: None
    [12:31:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:02] SMILES Parse Error: syntax error while parsing: None
    [12:31:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:03] SMILES Parse Error: syntax error while parsing: None
    [12:31:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:04] SMILES Parse Error: syntax error while parsing: None
    [12:31:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:05] SMILES Parse Error: syntax error while parsing: None
    [12:31:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:06] SMILES Parse Error: syntax error while parsing: None
    [12:31:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:07] SMILES Parse Error: syntax error while parsing: None
    [12:31:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:08] SMILES Parse Error: syntax error while parsing: None
    [12:31:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:09] SMILES Parse Error: syntax error while parsing: None
    [12:31:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:11] SMILES Parse Error: syntax error while parsing: None
    [12:31:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:13] SMILES Parse Error: syntax error while parsing: None
    [12:31:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:15] SMILES Parse Error: syntax error while parsing: None
    [12:31:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:16] SMILES Parse Error: syntax error while parsing: None
    [12:31:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:18] SMILES Parse Error: syntax error while parsing: None
    [12:31:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:19] SMILES Parse Error: syntax error while parsing: None
    [12:31:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:20] SMILES Parse Error: syntax error while parsing: None
    [12:31:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:21] SMILES Parse Error: syntax error while parsing: None
    [12:31:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:27] SMILES Parse Error: syntax error while parsing: None
    [12:31:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:29] SMILES Parse Error: syntax error while parsing: None
    [12:31:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:31] SMILES Parse Error: syntax error while parsing: None
    [12:31:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:32] SMILES Parse Error: syntax error while parsing: None
    [12:31:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:33] SMILES Parse Error: syntax error while parsing: None
    [12:31:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:34] SMILES Parse Error: syntax error while parsing: None
    [12:31:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:35] SMILES Parse Error: syntax error while parsing: None
    [12:31:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:36] SMILES Parse Error: syntax error while parsing: None
    [12:31:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:38] SMILES Parse Error: syntax error while parsing: None
    [12:31:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:39] SMILES Parse Error: syntax error while parsing: None
    [12:31:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:42] SMILES Parse Error: syntax error while parsing: None
    [12:31:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:43] SMILES Parse Error: syntax error while parsing: None
    [12:31:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:44] SMILES Parse Error: syntax error while parsing: None
    [12:31:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:45] SMILES Parse Error: syntax error while parsing: None
    [12:31:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:46] SMILES Parse Error: syntax error while parsing: None
    [12:31:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:47] SMILES Parse Error: syntax error while parsing: None
    [12:31:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:48] SMILES Parse Error: syntax error while parsing: None
    [12:31:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:49] SMILES Parse Error: syntax error while parsing: None
    [12:31:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:50] SMILES Parse Error: syntax error while parsing: None
    [12:31:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:51] SMILES Parse Error: syntax error while parsing: None
    [12:31:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:52] SMILES Parse Error: syntax error while parsing: None
    [12:31:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:53] SMILES Parse Error: syntax error while parsing: None
    [12:31:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:56] SMILES Parse Error: syntax error while parsing: None
    [12:31:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:57] SMILES Parse Error: syntax error while parsing: None
    [12:31:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:58] SMILES Parse Error: syntax error while parsing: None
    [12:31:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:31:59] SMILES Parse Error: syntax error while parsing: None
    [12:31:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:00] SMILES Parse Error: syntax error while parsing: None
    [12:32:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:01] SMILES Parse Error: syntax error while parsing: None
    [12:32:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:02] SMILES Parse Error: syntax error while parsing: None
    [12:32:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:03] SMILES Parse Error: syntax error while parsing: None
    [12:32:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:04] SMILES Parse Error: syntax error while parsing: None
    [12:32:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:06] SMILES Parse Error: syntax error while parsing: None
    [12:32:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:07] SMILES Parse Error: syntax error while parsing: None
    [12:32:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:08] SMILES Parse Error: syntax error while parsing: None
    [12:32:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:09] SMILES Parse Error: syntax error while parsing: None
    [12:32:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:10] SMILES Parse Error: syntax error while parsing: None
    [12:32:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:11] SMILES Parse Error: syntax error while parsing: None
    [12:32:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:12] SMILES Parse Error: syntax error while parsing: None
    [12:32:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:13] SMILES Parse Error: syntax error while parsing: None
    [12:32:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:14] SMILES Parse Error: syntax error while parsing: None
    [12:32:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:15] SMILES Parse Error: syntax error while parsing: None
    [12:32:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:16] SMILES Parse Error: syntax error while parsing: None
    [12:32:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:17] SMILES Parse Error: syntax error while parsing: None
    [12:32:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:18] SMILES Parse Error: syntax error while parsing: None
    [12:32:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:20] SMILES Parse Error: syntax error while parsing: None
    [12:32:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:21] SMILES Parse Error: syntax error while parsing: None
    [12:32:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:22] SMILES Parse Error: syntax error while parsing: None
    [12:32:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:23] SMILES Parse Error: syntax error while parsing: None
    [12:32:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:24] SMILES Parse Error: syntax error while parsing: None
    [12:32:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:25] SMILES Parse Error: syntax error while parsing: None
    [12:32:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:26] SMILES Parse Error: syntax error while parsing: None
    [12:32:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:27] SMILES Parse Error: syntax error while parsing: None
    [12:32:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:28] SMILES Parse Error: syntax error while parsing: None
    [12:32:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:29] SMILES Parse Error: syntax error while parsing: None
    [12:32:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:30] SMILES Parse Error: syntax error while parsing: None
    [12:32:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:31] SMILES Parse Error: syntax error while parsing: None
    [12:32:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:32] SMILES Parse Error: syntax error while parsing: None
    [12:32:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:33] SMILES Parse Error: syntax error while parsing: None
    [12:32:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:34] SMILES Parse Error: syntax error while parsing: None
    [12:32:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:37] SMILES Parse Error: syntax error while parsing: None
    [12:32:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:38] SMILES Parse Error: syntax error while parsing: None
    [12:32:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:39] SMILES Parse Error: syntax error while parsing: None
    [12:32:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:41] SMILES Parse Error: syntax error while parsing: None
    [12:32:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:42] SMILES Parse Error: syntax error while parsing: None
    [12:32:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:45] SMILES Parse Error: syntax error while parsing: None
    [12:32:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:47] SMILES Parse Error: syntax error while parsing: None
    [12:32:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:48] SMILES Parse Error: syntax error while parsing: None
    [12:32:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:50] SMILES Parse Error: syntax error while parsing: None
    [12:32:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:51] SMILES Parse Error: syntax error while parsing: None
    [12:32:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:52] SMILES Parse Error: syntax error while parsing: None
    [12:32:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:54] SMILES Parse Error: syntax error while parsing: None
    [12:32:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:56] SMILES Parse Error: syntax error while parsing: None
    [12:32:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:57] SMILES Parse Error: syntax error while parsing: None
    [12:32:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:32:58] SMILES Parse Error: syntax error while parsing: None
    [12:32:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:00] SMILES Parse Error: syntax error while parsing: None
    [12:33:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:02] SMILES Parse Error: syntax error while parsing: None
    [12:33:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:03] SMILES Parse Error: syntax error while parsing: None
    [12:33:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:04] SMILES Parse Error: syntax error while parsing: None
    [12:33:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:05] SMILES Parse Error: syntax error while parsing: None
    [12:33:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:06] SMILES Parse Error: syntax error while parsing: None
    [12:33:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:10] SMILES Parse Error: syntax error while parsing: None
    [12:33:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:11] SMILES Parse Error: syntax error while parsing: None
    [12:33:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:12] SMILES Parse Error: syntax error while parsing: None
    [12:33:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:13] SMILES Parse Error: syntax error while parsing: None
    [12:33:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:14] SMILES Parse Error: syntax error while parsing: None
    [12:33:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:16] SMILES Parse Error: syntax error while parsing: None
    [12:33:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:18] SMILES Parse Error: syntax error while parsing: None
    [12:33:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:19] SMILES Parse Error: syntax error while parsing: None
    [12:33:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:20] SMILES Parse Error: syntax error while parsing: None
    [12:33:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:22] SMILES Parse Error: syntax error while parsing: None
    [12:33:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:23] SMILES Parse Error: syntax error while parsing: None
    [12:33:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:24] SMILES Parse Error: syntax error while parsing: None
    [12:33:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:25] SMILES Parse Error: syntax error while parsing: None
    [12:33:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:26] SMILES Parse Error: syntax error while parsing: None
    [12:33:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:27] SMILES Parse Error: syntax error while parsing: None
    [12:33:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:28] SMILES Parse Error: syntax error while parsing: None
    [12:33:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:29] SMILES Parse Error: syntax error while parsing: None
    [12:33:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:31] SMILES Parse Error: syntax error while parsing: None
    [12:33:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:32] SMILES Parse Error: syntax error while parsing: None
    [12:33:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:33] SMILES Parse Error: syntax error while parsing: None
    [12:33:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:35] SMILES Parse Error: syntax error while parsing: None
    [12:33:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:36] SMILES Parse Error: syntax error while parsing: None
    [12:33:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:37] SMILES Parse Error: syntax error while parsing: None
    [12:33:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:39] SMILES Parse Error: syntax error while parsing: None
    [12:33:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:40] SMILES Parse Error: syntax error while parsing: None
    [12:33:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:41] SMILES Parse Error: syntax error while parsing: None
    [12:33:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:42] SMILES Parse Error: syntax error while parsing: None
    [12:33:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:43] SMILES Parse Error: syntax error while parsing: None
    [12:33:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:44] SMILES Parse Error: syntax error while parsing: None
    [12:33:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:45] SMILES Parse Error: syntax error while parsing: None
    [12:33:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:46] SMILES Parse Error: syntax error while parsing: None
    [12:33:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:51] SMILES Parse Error: syntax error while parsing: None
    [12:33:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:52] SMILES Parse Error: syntax error while parsing: None
    [12:33:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:53] SMILES Parse Error: syntax error while parsing: None
    [12:33:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:54] SMILES Parse Error: syntax error while parsing: None
    [12:33:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:55] SMILES Parse Error: syntax error while parsing: None
    [12:33:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:56] SMILES Parse Error: syntax error while parsing: None
    [12:33:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:57] SMILES Parse Error: syntax error while parsing: None
    [12:33:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:33:59] SMILES Parse Error: syntax error while parsing: None
    [12:33:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:00] SMILES Parse Error: syntax error while parsing: None
    [12:34:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:01] SMILES Parse Error: syntax error while parsing: None
    [12:34:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:02] SMILES Parse Error: syntax error while parsing: None
    [12:34:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:03] SMILES Parse Error: syntax error while parsing: None
    [12:34:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:04] SMILES Parse Error: syntax error while parsing: None
    [12:34:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:05] SMILES Parse Error: syntax error while parsing: None
    [12:34:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:06] SMILES Parse Error: syntax error while parsing: None
    [12:34:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:08] SMILES Parse Error: syntax error while parsing: None
    [12:34:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:09] SMILES Parse Error: syntax error while parsing: None
    [12:34:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:11] SMILES Parse Error: syntax error while parsing: None
    [12:34:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:12] SMILES Parse Error: syntax error while parsing: None
    [12:34:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:13] SMILES Parse Error: syntax error while parsing: None
    [12:34:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:14] SMILES Parse Error: syntax error while parsing: None
    [12:34:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:15] SMILES Parse Error: syntax error while parsing: None
    [12:34:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:16] SMILES Parse Error: syntax error while parsing: None
    [12:34:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:17] SMILES Parse Error: syntax error while parsing: None
    [12:34:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:18] SMILES Parse Error: syntax error while parsing: None
    [12:34:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:19] SMILES Parse Error: syntax error while parsing: None
    [12:34:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:20] SMILES Parse Error: syntax error while parsing: None
    [12:34:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:21] SMILES Parse Error: syntax error while parsing: None
    [12:34:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:22] SMILES Parse Error: syntax error while parsing: None
    [12:34:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:23] SMILES Parse Error: syntax error while parsing: None
    [12:34:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:24] SMILES Parse Error: syntax error while parsing: None
    [12:34:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:25] SMILES Parse Error: syntax error while parsing: None
    [12:34:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:26] SMILES Parse Error: syntax error while parsing: None
    [12:34:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:27] SMILES Parse Error: syntax error while parsing: None
    [12:34:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:28] SMILES Parse Error: syntax error while parsing: None
    [12:34:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:29] SMILES Parse Error: syntax error while parsing: None
    [12:34:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:31] SMILES Parse Error: syntax error while parsing: None
    [12:34:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:32] SMILES Parse Error: syntax error while parsing: None
    [12:34:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:33] SMILES Parse Error: syntax error while parsing: None
    [12:34:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:34] SMILES Parse Error: syntax error while parsing: None
    [12:34:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:35] SMILES Parse Error: syntax error while parsing: None
    [12:34:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:36] SMILES Parse Error: syntax error while parsing: None
    [12:34:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:37] SMILES Parse Error: syntax error while parsing: None
    [12:34:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:38] SMILES Parse Error: syntax error while parsing: None
    [12:34:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:39] SMILES Parse Error: syntax error while parsing: None
    [12:34:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:40] SMILES Parse Error: syntax error while parsing: None
    [12:34:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:41] SMILES Parse Error: syntax error while parsing: None
    [12:34:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:42] SMILES Parse Error: syntax error while parsing: None
    [12:34:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:43] SMILES Parse Error: syntax error while parsing: None
    [12:34:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:44] SMILES Parse Error: syntax error while parsing: None
    [12:34:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:45] SMILES Parse Error: syntax error while parsing: None
    [12:34:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:46] SMILES Parse Error: syntax error while parsing: None
    [12:34:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:48] SMILES Parse Error: syntax error while parsing: None
    [12:34:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:49] SMILES Parse Error: syntax error while parsing: None
    [12:34:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:50] SMILES Parse Error: syntax error while parsing: None
    [12:34:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:51] SMILES Parse Error: syntax error while parsing: None
    [12:34:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:53] SMILES Parse Error: syntax error while parsing: None
    [12:34:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:54] SMILES Parse Error: syntax error while parsing: None
    [12:34:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:55] SMILES Parse Error: syntax error while parsing: None
    [12:34:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:56] SMILES Parse Error: syntax error while parsing: None
    [12:34:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:57] SMILES Parse Error: syntax error while parsing: None
    [12:34:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:58] SMILES Parse Error: syntax error while parsing: None
    [12:34:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:34:59] SMILES Parse Error: syntax error while parsing: None
    [12:34:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:00] SMILES Parse Error: syntax error while parsing: None
    [12:35:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:01] SMILES Parse Error: syntax error while parsing: None
    [12:35:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:02] SMILES Parse Error: syntax error while parsing: None
    [12:35:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:03] SMILES Parse Error: syntax error while parsing: None
    [12:35:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:04] SMILES Parse Error: syntax error while parsing: None
    [12:35:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:05] SMILES Parse Error: syntax error while parsing: None
    [12:35:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:06] SMILES Parse Error: syntax error while parsing: None
    [12:35:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:07] SMILES Parse Error: syntax error while parsing: None
    [12:35:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:10] SMILES Parse Error: syntax error while parsing: None
    [12:35:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:11] SMILES Parse Error: syntax error while parsing: None
    [12:35:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:12] SMILES Parse Error: syntax error while parsing: None
    [12:35:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:13] SMILES Parse Error: syntax error while parsing: None
    [12:35:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:14] SMILES Parse Error: syntax error while parsing: None
    [12:35:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:15] SMILES Parse Error: syntax error while parsing: None
    [12:35:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:16] SMILES Parse Error: syntax error while parsing: None
    [12:35:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:17] SMILES Parse Error: syntax error while parsing: None
    [12:35:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:18] SMILES Parse Error: syntax error while parsing: None
    [12:35:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:19] SMILES Parse Error: syntax error while parsing: None
    [12:35:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:20] SMILES Parse Error: syntax error while parsing: None
    [12:35:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:22] SMILES Parse Error: syntax error while parsing: None
    [12:35:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:23] SMILES Parse Error: syntax error while parsing: None
    [12:35:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:24] SMILES Parse Error: syntax error while parsing: None
    [12:35:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:25] SMILES Parse Error: syntax error while parsing: None
    [12:35:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:26] SMILES Parse Error: syntax error while parsing: None
    [12:35:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:27] SMILES Parse Error: syntax error while parsing: None
    [12:35:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:28] SMILES Parse Error: syntax error while parsing: None
    [12:35:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:29] SMILES Parse Error: syntax error while parsing: None
    [12:35:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:30] SMILES Parse Error: syntax error while parsing: None
    [12:35:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:31] SMILES Parse Error: syntax error while parsing: None
    [12:35:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:32] SMILES Parse Error: syntax error while parsing: None
    [12:35:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:33] SMILES Parse Error: syntax error while parsing: None
    [12:35:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:34] SMILES Parse Error: syntax error while parsing: None
    [12:35:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:36] SMILES Parse Error: syntax error while parsing: None
    [12:35:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:37] SMILES Parse Error: syntax error while parsing: None
    [12:35:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:38] SMILES Parse Error: syntax error while parsing: None
    [12:35:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:39] SMILES Parse Error: syntax error while parsing: None
    [12:35:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:40] SMILES Parse Error: syntax error while parsing: None
    [12:35:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:41] SMILES Parse Error: syntax error while parsing: None
    [12:35:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:42] SMILES Parse Error: syntax error while parsing: None
    [12:35:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:43] SMILES Parse Error: syntax error while parsing: None
    [12:35:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:44] SMILES Parse Error: syntax error while parsing: None
    [12:35:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:45] SMILES Parse Error: syntax error while parsing: None
    [12:35:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:46] SMILES Parse Error: syntax error while parsing: None
    [12:35:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:47] SMILES Parse Error: syntax error while parsing: None
    [12:35:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:48] SMILES Parse Error: syntax error while parsing: None
    [12:35:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:49] SMILES Parse Error: syntax error while parsing: None
    [12:35:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:50] SMILES Parse Error: syntax error while parsing: None
    [12:35:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:52] SMILES Parse Error: syntax error while parsing: None
    [12:35:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:53] SMILES Parse Error: syntax error while parsing: None
    [12:35:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:54] SMILES Parse Error: syntax error while parsing: None
    [12:35:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:55] SMILES Parse Error: syntax error while parsing: None
    [12:35:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:56] SMILES Parse Error: syntax error while parsing: None
    [12:35:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:57] SMILES Parse Error: syntax error while parsing: None
    [12:35:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:58] SMILES Parse Error: syntax error while parsing: None
    [12:35:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:35:59] SMILES Parse Error: syntax error while parsing: None
    [12:35:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:00] SMILES Parse Error: syntax error while parsing: None
    [12:36:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:01] SMILES Parse Error: syntax error while parsing: None
    [12:36:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:02] SMILES Parse Error: syntax error while parsing: None
    [12:36:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:03] SMILES Parse Error: syntax error while parsing: None
    [12:36:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:04] SMILES Parse Error: syntax error while parsing: None
    [12:36:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:05] SMILES Parse Error: syntax error while parsing: None
    [12:36:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:06] SMILES Parse Error: syntax error while parsing: None
    [12:36:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:07] SMILES Parse Error: syntax error while parsing: None
    [12:36:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:09] SMILES Parse Error: syntax error while parsing: None
    [12:36:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:10] SMILES Parse Error: syntax error while parsing: None
    [12:36:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:11] SMILES Parse Error: syntax error while parsing: None
    [12:36:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:12] SMILES Parse Error: syntax error while parsing: None
    [12:36:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:13] SMILES Parse Error: syntax error while parsing: None
    [12:36:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:14] SMILES Parse Error: syntax error while parsing: None
    [12:36:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:15] SMILES Parse Error: syntax error while parsing: None
    [12:36:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:16] SMILES Parse Error: syntax error while parsing: None
    [12:36:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:17] SMILES Parse Error: syntax error while parsing: None
    [12:36:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:18] SMILES Parse Error: syntax error while parsing: None
    [12:36:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:19] SMILES Parse Error: syntax error while parsing: None
    [12:36:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:20] SMILES Parse Error: syntax error while parsing: None
    [12:36:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:25] SMILES Parse Error: syntax error while parsing: None
    [12:36:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:26] SMILES Parse Error: syntax error while parsing: None
    [12:36:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:27] SMILES Parse Error: syntax error while parsing: None
    [12:36:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:28] SMILES Parse Error: syntax error while parsing: None
    [12:36:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:30] SMILES Parse Error: syntax error while parsing: None
    [12:36:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:32] SMILES Parse Error: syntax error while parsing: None
    [12:36:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:33] SMILES Parse Error: syntax error while parsing: None
    [12:36:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:34] SMILES Parse Error: syntax error while parsing: None
    [12:36:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:35] SMILES Parse Error: syntax error while parsing: None
    [12:36:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:36] SMILES Parse Error: syntax error while parsing: None
    [12:36:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:37] SMILES Parse Error: syntax error while parsing: None
    [12:36:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:38] SMILES Parse Error: syntax error while parsing: None
    [12:36:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:39] SMILES Parse Error: syntax error while parsing: None
    [12:36:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:40] SMILES Parse Error: syntax error while parsing: None
    [12:36:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:42] SMILES Parse Error: syntax error while parsing: None
    [12:36:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:43] SMILES Parse Error: syntax error while parsing: None
    [12:36:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:44] SMILES Parse Error: syntax error while parsing: None
    [12:36:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:45] SMILES Parse Error: syntax error while parsing: None
    [12:36:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:46] SMILES Parse Error: syntax error while parsing: None
    [12:36:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:55] SMILES Parse Error: syntax error while parsing: None
    [12:36:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:58] SMILES Parse Error: syntax error while parsing: None
    [12:36:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:36:59] SMILES Parse Error: syntax error while parsing: None
    [12:36:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:00] SMILES Parse Error: syntax error while parsing: None
    [12:37:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:01] SMILES Parse Error: syntax error while parsing: None
    [12:37:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:02] SMILES Parse Error: syntax error while parsing: None
    [12:37:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:03] SMILES Parse Error: syntax error while parsing: None
    [12:37:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:17] SMILES Parse Error: syntax error while parsing: None
    [12:37:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:18] SMILES Parse Error: syntax error while parsing: None
    [12:37:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:35] SMILES Parse Error: syntax error while parsing: None
    [12:37:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:46] SMILES Parse Error: syntax error while parsing: None
    [12:37:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:47] SMILES Parse Error: syntax error while parsing: None
    [12:37:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:37:50] SMILES Parse Error: syntax error while parsing: None
    [12:37:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:02] SMILES Parse Error: syntax error while parsing: None
    [12:38:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:05] SMILES Parse Error: syntax error while parsing: None
    [12:38:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:06] SMILES Parse Error: syntax error while parsing: None
    [12:38:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:07] SMILES Parse Error: syntax error while parsing: None
    [12:38:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:08] SMILES Parse Error: syntax error while parsing: None
    [12:38:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:09] SMILES Parse Error: syntax error while parsing: None
    [12:38:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:10] SMILES Parse Error: syntax error while parsing: None
    [12:38:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:11] SMILES Parse Error: syntax error while parsing: None
    [12:38:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:12] SMILES Parse Error: syntax error while parsing: None
    [12:38:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:16] SMILES Parse Error: syntax error while parsing: None
    [12:38:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:17] SMILES Parse Error: syntax error while parsing: None
    [12:38:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:19] SMILES Parse Error: syntax error while parsing: None
    [12:38:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:21] SMILES Parse Error: syntax error while parsing: None
    [12:38:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:22] SMILES Parse Error: syntax error while parsing: None
    [12:38:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:24] SMILES Parse Error: syntax error while parsing: None
    [12:38:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:25] SMILES Parse Error: syntax error while parsing: None
    [12:38:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:26] SMILES Parse Error: syntax error while parsing: None
    [12:38:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:27] SMILES Parse Error: syntax error while parsing: None
    [12:38:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:28] SMILES Parse Error: syntax error while parsing: None
    [12:38:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:29] SMILES Parse Error: syntax error while parsing: None
    [12:38:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:30] SMILES Parse Error: syntax error while parsing: None
    [12:38:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:31] SMILES Parse Error: syntax error while parsing: None
    [12:38:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:32] SMILES Parse Error: syntax error while parsing: None
    [12:38:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:37] SMILES Parse Error: syntax error while parsing: None
    [12:38:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:38] SMILES Parse Error: syntax error while parsing: None
    [12:38:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:39] SMILES Parse Error: syntax error while parsing: None
    [12:38:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:41] SMILES Parse Error: syntax error while parsing: None
    [12:38:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:42] SMILES Parse Error: syntax error while parsing: None
    [12:38:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:43] SMILES Parse Error: syntax error while parsing: None
    [12:38:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:44] SMILES Parse Error: syntax error while parsing: None
    [12:38:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:45] SMILES Parse Error: syntax error while parsing: None
    [12:38:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:46] SMILES Parse Error: syntax error while parsing: None
    [12:38:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:47] SMILES Parse Error: syntax error while parsing: None
    [12:38:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:48] SMILES Parse Error: syntax error while parsing: None
    [12:38:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:49] SMILES Parse Error: syntax error while parsing: None
    [12:38:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:50] SMILES Parse Error: syntax error while parsing: None
    [12:38:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:51] SMILES Parse Error: syntax error while parsing: None
    [12:38:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:55] SMILES Parse Error: syntax error while parsing: None
    [12:38:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:56] SMILES Parse Error: syntax error while parsing: None
    [12:38:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:57] SMILES Parse Error: syntax error while parsing: None
    [12:38:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:58] SMILES Parse Error: syntax error while parsing: None
    [12:38:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:38:59] SMILES Parse Error: syntax error while parsing: None
    [12:38:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:00] SMILES Parse Error: syntax error while parsing: None
    [12:39:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:02] SMILES Parse Error: syntax error while parsing: None
    [12:39:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:03] SMILES Parse Error: syntax error while parsing: None
    [12:39:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:04] SMILES Parse Error: syntax error while parsing: None
    [12:39:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:10] SMILES Parse Error: syntax error while parsing: None
    [12:39:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:11] SMILES Parse Error: syntax error while parsing: None
    [12:39:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:12] SMILES Parse Error: syntax error while parsing: None
    [12:39:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:13] SMILES Parse Error: syntax error while parsing: None
    [12:39:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:14] SMILES Parse Error: syntax error while parsing: None
    [12:39:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:15] SMILES Parse Error: syntax error while parsing: None
    [12:39:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:16] SMILES Parse Error: syntax error while parsing: None
    [12:39:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:18] SMILES Parse Error: syntax error while parsing: None
    [12:39:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:19] SMILES Parse Error: syntax error while parsing: None
    [12:39:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:20] SMILES Parse Error: syntax error while parsing: None
    [12:39:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:21] SMILES Parse Error: syntax error while parsing: None
    [12:39:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:22] SMILES Parse Error: syntax error while parsing: None
    [12:39:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:24] SMILES Parse Error: syntax error while parsing: None
    [12:39:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:25] SMILES Parse Error: syntax error while parsing: None
    [12:39:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:26] SMILES Parse Error: syntax error while parsing: None
    [12:39:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:27] SMILES Parse Error: syntax error while parsing: None
    [12:39:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:28] SMILES Parse Error: syntax error while parsing: None
    [12:39:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:29] SMILES Parse Error: syntax error while parsing: None
    [12:39:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:31] SMILES Parse Error: syntax error while parsing: None
    [12:39:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:34] SMILES Parse Error: syntax error while parsing: None
    [12:39:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:35] SMILES Parse Error: syntax error while parsing: None
    [12:39:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:36] SMILES Parse Error: syntax error while parsing: None
    [12:39:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:37] SMILES Parse Error: syntax error while parsing: None
    [12:39:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:38] SMILES Parse Error: syntax error while parsing: None
    [12:39:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:39] SMILES Parse Error: syntax error while parsing: None
    [12:39:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:41] SMILES Parse Error: syntax error while parsing: None
    [12:39:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:42] SMILES Parse Error: syntax error while parsing: None
    [12:39:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:43] SMILES Parse Error: syntax error while parsing: None
    [12:39:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:45] SMILES Parse Error: syntax error while parsing: None
    [12:39:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:46] SMILES Parse Error: syntax error while parsing: None
    [12:39:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:47] SMILES Parse Error: syntax error while parsing: None
    [12:39:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:48] SMILES Parse Error: syntax error while parsing: None
    [12:39:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:49] SMILES Parse Error: syntax error while parsing: None
    [12:39:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:50] SMILES Parse Error: syntax error while parsing: None
    [12:39:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:51] SMILES Parse Error: syntax error while parsing: None
    [12:39:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:52] SMILES Parse Error: syntax error while parsing: None
    [12:39:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:53] SMILES Parse Error: syntax error while parsing: None
    [12:39:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:54] SMILES Parse Error: syntax error while parsing: None
    [12:39:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:55] SMILES Parse Error: syntax error while parsing: None
    [12:39:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:56] SMILES Parse Error: syntax error while parsing: None
    [12:39:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:57] SMILES Parse Error: syntax error while parsing: None
    [12:39:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:39:58] SMILES Parse Error: syntax error while parsing: None
    [12:39:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:00] SMILES Parse Error: syntax error while parsing: None
    [12:40:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:01] SMILES Parse Error: syntax error while parsing: None
    [12:40:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:02] SMILES Parse Error: syntax error while parsing: None
    [12:40:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:03] SMILES Parse Error: syntax error while parsing: None
    [12:40:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:04] SMILES Parse Error: syntax error while parsing: None
    [12:40:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:05] SMILES Parse Error: syntax error while parsing: None
    [12:40:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:06] SMILES Parse Error: syntax error while parsing: None
    [12:40:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:07] SMILES Parse Error: syntax error while parsing: None
    [12:40:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:08] SMILES Parse Error: syntax error while parsing: None
    [12:40:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:09] SMILES Parse Error: syntax error while parsing: None
    [12:40:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:11] SMILES Parse Error: syntax error while parsing: None
    [12:40:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:12] SMILES Parse Error: syntax error while parsing: None
    [12:40:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:13] SMILES Parse Error: syntax error while parsing: None
    [12:40:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:14] SMILES Parse Error: syntax error while parsing: None
    [12:40:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:15] SMILES Parse Error: syntax error while parsing: None
    [12:40:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:17] SMILES Parse Error: syntax error while parsing: None
    [12:40:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:18] SMILES Parse Error: syntax error while parsing: None
    [12:40:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:19] SMILES Parse Error: syntax error while parsing: None
    [12:40:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:20] SMILES Parse Error: syntax error while parsing: None
    [12:40:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:21] SMILES Parse Error: syntax error while parsing: None
    [12:40:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:22] SMILES Parse Error: syntax error while parsing: None
    [12:40:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:24] SMILES Parse Error: syntax error while parsing: None
    [12:40:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:25] SMILES Parse Error: syntax error while parsing: None
    [12:40:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:26] SMILES Parse Error: syntax error while parsing: None
    [12:40:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:28] SMILES Parse Error: syntax error while parsing: None
    [12:40:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:29] SMILES Parse Error: syntax error while parsing: None
    [12:40:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:30] SMILES Parse Error: syntax error while parsing: None
    [12:40:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:31] SMILES Parse Error: syntax error while parsing: None
    [12:40:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:32] SMILES Parse Error: syntax error while parsing: None
    [12:40:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:33] SMILES Parse Error: syntax error while parsing: None
    [12:40:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:34] SMILES Parse Error: syntax error while parsing: None
    [12:40:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:35] SMILES Parse Error: syntax error while parsing: None
    [12:40:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:37] SMILES Parse Error: syntax error while parsing: None
    [12:40:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:38] SMILES Parse Error: syntax error while parsing: None
    [12:40:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:39] SMILES Parse Error: syntax error while parsing: None
    [12:40:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:40] SMILES Parse Error: syntax error while parsing: None
    [12:40:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:41] SMILES Parse Error: syntax error while parsing: None
    [12:40:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:42] SMILES Parse Error: syntax error while parsing: None
    [12:40:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:47] SMILES Parse Error: syntax error while parsing: None
    [12:40:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:49] SMILES Parse Error: syntax error while parsing: None
    [12:40:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:50] SMILES Parse Error: syntax error while parsing: None
    [12:40:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:51] SMILES Parse Error: syntax error while parsing: None
    [12:40:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:53] SMILES Parse Error: syntax error while parsing: None
    [12:40:53] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:54] SMILES Parse Error: syntax error while parsing: None
    [12:40:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:55] SMILES Parse Error: syntax error while parsing: None
    [12:40:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:40:56] SMILES Parse Error: syntax error while parsing: None
    [12:40:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:02] SMILES Parse Error: syntax error while parsing: None
    [12:41:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:03] SMILES Parse Error: syntax error while parsing: None
    [12:41:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:04] SMILES Parse Error: syntax error while parsing: None
    [12:41:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:05] SMILES Parse Error: syntax error while parsing: None
    [12:41:05] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:06] SMILES Parse Error: syntax error while parsing: None
    [12:41:06] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:07] SMILES Parse Error: syntax error while parsing: None
    [12:41:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:09] SMILES Parse Error: syntax error while parsing: None
    [12:41:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:10] SMILES Parse Error: syntax error while parsing: None
    [12:41:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:11] SMILES Parse Error: syntax error while parsing: None
    [12:41:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:12] SMILES Parse Error: syntax error while parsing: None
    [12:41:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:13] SMILES Parse Error: syntax error while parsing: None
    [12:41:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:15] SMILES Parse Error: syntax error while parsing: None
    [12:41:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:16] SMILES Parse Error: syntax error while parsing: None
    [12:41:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:17] SMILES Parse Error: syntax error while parsing: None
    [12:41:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:18] SMILES Parse Error: syntax error while parsing: None
    [12:41:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:25] SMILES Parse Error: syntax error while parsing: None
    [12:41:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:27] SMILES Parse Error: syntax error while parsing: None
    [12:41:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:28] SMILES Parse Error: syntax error while parsing: None
    [12:41:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:30] SMILES Parse Error: syntax error while parsing: None
    [12:41:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:31] SMILES Parse Error: syntax error while parsing: None
    [12:41:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:36] SMILES Parse Error: syntax error while parsing: None
    [12:41:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:37] SMILES Parse Error: syntax error while parsing: None
    [12:41:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:39] SMILES Parse Error: syntax error while parsing: None
    [12:41:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:40] SMILES Parse Error: syntax error while parsing: None
    [12:41:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:41] SMILES Parse Error: syntax error while parsing: None
    [12:41:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:42] SMILES Parse Error: syntax error while parsing: None
    [12:41:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:43] SMILES Parse Error: syntax error while parsing: None
    [12:41:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:44] SMILES Parse Error: syntax error while parsing: None
    [12:41:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:45] SMILES Parse Error: syntax error while parsing: None
    [12:41:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:46] SMILES Parse Error: syntax error while parsing: None
    [12:41:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:47] SMILES Parse Error: syntax error while parsing: None
    [12:41:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:48] SMILES Parse Error: syntax error while parsing: None
    [12:41:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:49] SMILES Parse Error: syntax error while parsing: None
    [12:41:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:50] SMILES Parse Error: syntax error while parsing: None
    [12:41:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:51] SMILES Parse Error: syntax error while parsing: None
    [12:41:51] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:52] SMILES Parse Error: syntax error while parsing: None
    [12:41:52] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:54] SMILES Parse Error: syntax error while parsing: None
    [12:41:54] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:55] SMILES Parse Error: syntax error while parsing: None
    [12:41:55] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:56] SMILES Parse Error: syntax error while parsing: None
    [12:41:56] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:57] SMILES Parse Error: syntax error while parsing: None
    [12:41:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:58] SMILES Parse Error: syntax error while parsing: None
    [12:41:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:41:59] SMILES Parse Error: syntax error while parsing: None
    [12:41:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:00] SMILES Parse Error: syntax error while parsing: None
    [12:42:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:01] SMILES Parse Error: syntax error while parsing: None
    [12:42:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:02] SMILES Parse Error: syntax error while parsing: None
    [12:42:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:10] SMILES Parse Error: syntax error while parsing: None
    [12:42:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:11] SMILES Parse Error: syntax error while parsing: None
    [12:42:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:12] SMILES Parse Error: syntax error while parsing: None
    [12:42:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:13] SMILES Parse Error: syntax error while parsing: None
    [12:42:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:14] SMILES Parse Error: syntax error while parsing: None
    [12:42:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:15] SMILES Parse Error: syntax error while parsing: None
    [12:42:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:16] SMILES Parse Error: syntax error while parsing: None
    [12:42:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:17] SMILES Parse Error: syntax error while parsing: None
    [12:42:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:18] SMILES Parse Error: syntax error while parsing: None
    [12:42:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:19] SMILES Parse Error: syntax error while parsing: None
    [12:42:19] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:20] SMILES Parse Error: syntax error while parsing: None
    [12:42:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:21] SMILES Parse Error: syntax error while parsing: None
    [12:42:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:22] SMILES Parse Error: syntax error while parsing: None
    [12:42:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:23] SMILES Parse Error: syntax error while parsing: None
    [12:42:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:24] SMILES Parse Error: syntax error while parsing: None
    [12:42:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:25] SMILES Parse Error: syntax error while parsing: None
    [12:42:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:27] SMILES Parse Error: syntax error while parsing: None
    [12:42:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:28] SMILES Parse Error: syntax error while parsing: None
    [12:42:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:29] SMILES Parse Error: syntax error while parsing: None
    [12:42:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:30] SMILES Parse Error: syntax error while parsing: None
    [12:42:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:31] SMILES Parse Error: syntax error while parsing: None
    [12:42:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:32] SMILES Parse Error: syntax error while parsing: None
    [12:42:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:33] SMILES Parse Error: syntax error while parsing: None
    [12:42:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:34] SMILES Parse Error: syntax error while parsing: None
    [12:42:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:35] SMILES Parse Error: syntax error while parsing: None
    [12:42:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:36] SMILES Parse Error: syntax error while parsing: None
    [12:42:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:37] SMILES Parse Error: syntax error while parsing: None
    [12:42:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:38] SMILES Parse Error: syntax error while parsing: None
    [12:42:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:39] SMILES Parse Error: syntax error while parsing: None
    [12:42:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:40] SMILES Parse Error: syntax error while parsing: None
    [12:42:40] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:41] SMILES Parse Error: syntax error while parsing: None
    [12:42:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:43] SMILES Parse Error: syntax error while parsing: None
    [12:42:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:44] SMILES Parse Error: syntax error while parsing: None
    [12:42:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:45] SMILES Parse Error: syntax error while parsing: None
    [12:42:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:46] SMILES Parse Error: syntax error while parsing: None
    [12:42:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:47] SMILES Parse Error: syntax error while parsing: None
    [12:42:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:57] SMILES Parse Error: syntax error while parsing: None
    [12:42:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:58] SMILES Parse Error: syntax error while parsing: None
    [12:42:58] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:42:59] SMILES Parse Error: syntax error while parsing: None
    [12:42:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:00] SMILES Parse Error: syntax error while parsing: None
    [12:43:00] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:01] SMILES Parse Error: syntax error while parsing: None
    [12:43:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:02] SMILES Parse Error: syntax error while parsing: None
    [12:43:02] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:07] SMILES Parse Error: syntax error while parsing: None
    [12:43:07] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:08] SMILES Parse Error: syntax error while parsing: None
    [12:43:08] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:09] SMILES Parse Error: syntax error while parsing: None
    [12:43:09] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:10] SMILES Parse Error: syntax error while parsing: None
    [12:43:10] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:11] SMILES Parse Error: syntax error while parsing: None
    [12:43:11] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:12] SMILES Parse Error: syntax error while parsing: None
    [12:43:12] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:13] SMILES Parse Error: syntax error while parsing: None
    [12:43:13] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:14] SMILES Parse Error: syntax error while parsing: None
    [12:43:14] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:15] SMILES Parse Error: syntax error while parsing: None
    [12:43:15] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:16] SMILES Parse Error: syntax error while parsing: None
    [12:43:16] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:17] SMILES Parse Error: syntax error while parsing: None
    [12:43:17] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:18] SMILES Parse Error: syntax error while parsing: None
    [12:43:18] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:20] SMILES Parse Error: syntax error while parsing: None
    [12:43:20] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:21] SMILES Parse Error: syntax error while parsing: None
    [12:43:21] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:22] SMILES Parse Error: syntax error while parsing: None
    [12:43:22] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:23] SMILES Parse Error: syntax error while parsing: None
    [12:43:23] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:24] SMILES Parse Error: syntax error while parsing: None
    [12:43:24] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:25] SMILES Parse Error: syntax error while parsing: None
    [12:43:25] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:26] SMILES Parse Error: syntax error while parsing: None
    [12:43:26] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:27] SMILES Parse Error: syntax error while parsing: None
    [12:43:27] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:28] SMILES Parse Error: syntax error while parsing: None
    [12:43:28] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:29] SMILES Parse Error: syntax error while parsing: None
    [12:43:29] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:30] SMILES Parse Error: syntax error while parsing: None
    [12:43:30] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:31] SMILES Parse Error: syntax error while parsing: None
    [12:43:31] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:32] SMILES Parse Error: syntax error while parsing: None
    [12:43:32] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:33] SMILES Parse Error: syntax error while parsing: None
    [12:43:33] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:34] SMILES Parse Error: syntax error while parsing: None
    [12:43:34] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:35] SMILES Parse Error: syntax error while parsing: None
    [12:43:35] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:36] SMILES Parse Error: syntax error while parsing: None
    [12:43:36] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:37] SMILES Parse Error: syntax error while parsing: None
    [12:43:37] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:38] SMILES Parse Error: syntax error while parsing: None
    [12:43:38] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:39] SMILES Parse Error: syntax error while parsing: None
    [12:43:39] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:41] SMILES Parse Error: syntax error while parsing: None
    [12:43:41] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:42] SMILES Parse Error: syntax error while parsing: None
    [12:43:42] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:43] SMILES Parse Error: syntax error while parsing: None
    [12:43:43] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:44] SMILES Parse Error: syntax error while parsing: None
    [12:43:44] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:45] SMILES Parse Error: syntax error while parsing: None
    [12:43:45] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:46] SMILES Parse Error: syntax error while parsing: None
    [12:43:46] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:47] SMILES Parse Error: syntax error while parsing: None
    [12:43:47] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:48] SMILES Parse Error: syntax error while parsing: None
    [12:43:48] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:49] SMILES Parse Error: syntax error while parsing: None
    [12:43:49] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:50] SMILES Parse Error: syntax error while parsing: None
    [12:43:50] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:57] SMILES Parse Error: syntax error while parsing: None
    [12:43:57] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:43:59] SMILES Parse Error: syntax error while parsing: None
    [12:43:59] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:44:01] SMILES Parse Error: syntax error while parsing: None
    [12:44:01] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:44:03] SMILES Parse Error: syntax error while parsing: None
    [12:44:03] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'
    [12:44:04] SMILES Parse Error: syntax error while parsing: None
    [12:44:04] SMILES Parse Error: Failed parsing SMILES 'None' for input: 'None'


    It takes 1198.0919036865234 seconds to predict Kcat values!
    -----------------------------------
    Prediction success!
    DLKcat done!
    
    0:20:00.334137


## Step 7: get the kcat_mw file


```python
starttime=datetime.datetime.now()
# Step 7: get the kcat_mw file
print("Starting to get reaction kcat_mw for model......")
DLouputdf = pd.read_csv(DLouputdf_file, sep='\t')
comdf = pd.read_csv(comdf_file)
DL_reaction_kact_mw = DL_kcat_mw_calculation(DLouputdf, comdf)
DL_reaction_kact_mw.to_csv(DL_reaction_kact_mw_file, index=False)
print("Reaction kcat_mw done!")

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to get reaction kcat_mw for model......
    DL_reaction_kact_mw generated
    Reaction kcat_mw done!
    0:00:00.084022


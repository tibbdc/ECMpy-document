# 9.BiGG_to_ecGEM
## Batch construction of ecGEM from BiGG models 

### Download BiGG Model


```python
!mkdir ./data/BIGG_model
```


```python
# ignore model
wget http://bigg.ucsd.edu/static/models/iAB_RBC_283.json#Homo sapiens
wget http://bigg.ucsd.edu/static/models/iAT_PLT_636.json#Homo sapiens
wget http://bigg.ucsd.edu/static/models/RECON1.json#Homo sapiens
wget http://bigg.ucsd.edu/static/models/Recon3D.json#Homo sapiens

wget http://bigg.ucsd.edu/static/models/iMM1415.json#Mus musculus

wget http://bigg.ucsd.edu/static/models/e_coli_core.json#Escherichia coli str. K-12 substr. MG1655
wget http://bigg.ucsd.edu/static/models/iAF1260b.json#Escherichia coli str. K-12 substr. MG1655
wget http://bigg.ucsd.edu/static/models/iEC1356_Bl21DE3.json#Escherichia coli BL21(DE3)  
wget http://bigg.ucsd.edu/static/models/iECD_1391.json#	Escherichia coli BL21(DE3)
wget http://bigg.ucsd.edu/static/models/iECDH1ME8569_1439.json#Escherichia coli DH1
wget http://bigg.ucsd.edu/static/models/iAF1260.json#Escherichia coli str. K-12 substr. MG1655
wget http://bigg.ucsd.edu/static/models/iAPECO1_1312.json#Escherichia coli APEC O1  
wget http://bigg.ucsd.edu/static/models/iBWG_1329.json#Escherichia coli BW2952
wget http://bigg.ucsd.edu/static/models/ic_1306.json#Escherichia coli CFT073
wget http://bigg.ucsd.edu/static/models/iE2348C_1286.json#Escherichia coli O127:H6 str. E2348/69
wget http://bigg.ucsd.edu/static/models/iEC042_1314.json#Escherichia coli 042
wget http://bigg.ucsd.edu/static/models/iEC1344_C.json#Escherichia coli C
wget http://bigg.ucsd.edu/static/models/iEC1349_Crooks.json#Escherichia coli ATCC 8739 
wget http://bigg.ucsd.edu/static/models/iEC1364_W.json#Escherichia coli W
wget http://bigg.ucsd.edu/static/models/iEC1368_DH5a.json#Escherichia coli DH5[alpha]
wget http://bigg.ucsd.edu/static/models/iB21_1397.json#Escherichia coli BL21(DE3)
wget http://bigg.ucsd.edu/static/models/iEC55989_1330.json#Escherichia coli 55989
wget http://bigg.ucsd.edu/static/models/iECABU_c1320.json#Escherichia coli ABU 83972
wget http://bigg.ucsd.edu/static/models/iECB_1328.json#Escherichia coli B str. REL606
wget http://bigg.ucsd.edu/static/models/iECBD_1354.json#Escherichia coli 'BL21-Gold(DE3)pLysS AG'
wget http://bigg.ucsd.edu/static/models/iECDH10B_1368.json#Escherichia coli str. K-12 substr. DH10B
wget http://bigg.ucsd.edu/static/models/iEcDH1_1363.json#Escherichia coli DH1
wget http://bigg.ucsd.edu/static/models/iEcE24377_1341.json#Escherichia coli O139:H28 str. E24377A
wget http://bigg.ucsd.edu/static/models/iECED1_1282.json#Escherichia coli ED1a          
wget http://bigg.ucsd.edu/static/models/iECH74115_1262.json#Escherichia coli O157:H7 str. EC4115
wget http://bigg.ucsd.edu/static/models/iEcHS_1320.json#Escherichia coli HS
wget http://bigg.ucsd.edu/static/models/iECIAI1_1343.json#Escherichia coli IAI1
wget http://bigg.ucsd.edu/static/models/iECIAI39_1322.json#Escherichia coli IAI39
wget http://bigg.ucsd.edu/static/models/iECNA114_1301.json#Escherichia coli NA114
wget http://bigg.ucsd.edu/static/models/iECO103_1326.json#Escherichia coli O103:H2 str. 12009
wget http://bigg.ucsd.edu/static/models/iECO111_1330.json#Escherichia coli O111:H- str. 11128
wget http://bigg.ucsd.edu/static/models/iECO26_1355.json#Escherichia coli O26:H11 str. 11368
wget http://bigg.ucsd.edu/static/models/iECOK1_1307.json#Escherichia coli IHE3034
wget http://bigg.ucsd.edu/static/models/iEcolC_1368.json#Escherichia coli ATCC 8739
wget http://bigg.ucsd.edu/static/models/iECP_1309.json#Escherichia coli 536
wget http://bigg.ucsd.edu/static/models/iECs_1301.json#Escherichia coli O157:H7 str. Sakai  
wget http://bigg.ucsd.edu/static/models/iECS88_1305.json#Escherichia coli S88         
wget http://bigg.ucsd.edu/static/models/iECSE_1348.json#Escherichia coli SE11
wget http://bigg.ucsd.edu/static/models/iECSF_1327.json#Escherichia coli SE15
wget http://bigg.ucsd.edu/static/models/iEcSMS35_1347.json#Escherichia coli SMS-3-5
wget http://bigg.ucsd.edu/static/models/iECSP_1301.json#Escherichia coli O157:H7 str. TW14359
wget http://bigg.ucsd.edu/static/models/iECUMN_1333.json#Escherichia coli UMN026
wget http://bigg.ucsd.edu/static/models/iECW_1372.json#Escherichia coli W
wget http://bigg.ucsd.edu/static/models/iEKO11_1354.json#	Escherichia coli KO11FL
wget http://bigg.ucsd.edu/static/models/iETEC_1333.json#Escherichia coli ETEC H10407
wget http://bigg.ucsd.edu/static/models/iG2583_1286.json#	Escherichia coli O55:H7 str. CB9615
wget http://bigg.ucsd.edu/static/models/iJO1366.json#Escherichia coli str. K-12 substr. MG1655
wget http://bigg.ucsd.edu/static/models/iJR904.json#Escherichia coli str. K-12 substr. MG1655
wget http://bigg.ucsd.edu/static/models/iLF82_1304.json#Escherichia coli LF82
wget http://bigg.ucsd.edu/static/models/iML1515.json#Escherichia coli str. K-12 substr. MG1655      
wget http://bigg.ucsd.edu/static/models/iUMN146_1321.json#Escherichia coli UM146
wget http://bigg.ucsd.edu/static/models/iUMNK88_1353.json#Escherichia coli UMNK88
wget http://bigg.ucsd.edu/static/models/iUTI89_1310.json#Escherichia coli UTI89
wget http://bigg.ucsd.edu/static/models/iWFL_1372.json#Escherichia coli W
wget http://bigg.ucsd.edu/static/models/iY75_1357.json#	Escherichia coli str. K-12 substr. W3110
wget http://bigg.ucsd.edu/static/models/iZ_1308.json#	Escherichia coli O157:H7 str. EDL933
wget http://bigg.ucsd.edu/static/models/iNRG857_1313.json#Escherichia coli O83:H1 str. NRG 857C
wget http://bigg.ucsd.edu/static/models/iEC1372_W3110.json#Escherichia coli str. K-12 substr. W3110

wget http://bigg.ucsd.edu/static/models/iIS312_Amastigote.json#Trypanosoma cruzi Dm28c       
wget http://bigg.ucsd.edu/static/models/iIS312_Epimastigote.json#Trypanosoma cruzi Dm28c
wget http://bigg.ucsd.edu/static/models/iIS312_Trypomastigote.json#Trypanosoma cruzi Dm28c

wget http://bigg.ucsd.edu/static/models/iAM_Pb448.json#Plasmodium berghei
wget http://bigg.ucsd.edu/static/models/iAM_Pc455.json#Plasmodium cynomolgi strain B
wget http://bigg.ucsd.edu/static/models/iAM_Pk459.json#Plasmodium knowlesi strain H
wget http://bigg.ucsd.edu/static/models/iAM_Pv461.json#Plasmodium vivax Sal-1


wget http://bigg.ucsd.edu/static/models/iS_1188.json#Shigella flexneri 2a str. 2457T
wget http://bigg.ucsd.edu/static/models/iSbBS512_1146.json#Shigella boydii CDC 3083-94
wget http://bigg.ucsd.edu/static/models/iSBO_1134.json#Shigella boydii Sb227
wget http://bigg.ucsd.edu/static/models/iSDY_1059.json#Shigella dysenteriae Sd197    
wget http://bigg.ucsd.edu/static/models/iSF_1195.json#Shigella flexneri 2a str. 301
wget http://bigg.ucsd.edu/static/models/iSFV_1184.json#Shigella flexneri 5 str. 8401
wget http://bigg.ucsd.edu/static/models/iSFxv_1172.json#Shigella flexneri 2002017

wget http://bigg.ucsd.edu/static/models/iYS1720.json#Salmonella pan-reactome

wget http://bigg.ucsd.edu/static/models/iNJ661.json#Mycobacterium tuberculosis H37Rv

wget http://bigg.ucsd.edu/static/models/iJN678.json#Synechocystis sp. PCC 6803

wget http://bigg.ucsd.edu/static/models/iSB619.json#Staphylococcus aureus subsp. aureus N315

wget http://bigg.ucsd.edu/static/models/iCHOv1_DG44.json#Cricetulus griseus

wget http://bigg.ucsd.edu/static/models/iND750.json#Saccharomyces cerevisiae S288C

wget http://bigg.ucsd.edu/static/models/iJN746.json#Pseudomonas putida KT2440

#Constructing ecGEMs using GEMs
wget http://bigg.ucsd.edu/static/models/iAF692.json#Methanosarcina barkeri str. Fusaro
wget http://bigg.ucsd.edu/static/models/iAF987.json#Geobacter metallireducens GS-15
wget http://bigg.ucsd.edu/static/models/iAM_Pf480.json#Plasmodium falciparum 3D7
wget http://bigg.ucsd.edu/static/models/iCHOv1.json#Cricetulus griseus———can not found UniProt information
wget http://bigg.ucsd.edu/static/models/iCN718.json#Acinetobacter baumannii AYE——can not found UniProt information
wget http://bigg.ucsd.edu/static/models/iCN900.json#Clostridioides difficile 630
wget http://bigg.ucsd.edu/static/models/iEK1008.json#Mycobacterium tuberculosis H37Rv
wget http://bigg.ucsd.edu/static/models/iHN637.json#Clostridium ljungdahlii DSM 13528
wget http://bigg.ucsd.edu/static/models/iIS312.json#Trypanosoma cruzi Dm28c 
wget http://bigg.ucsd.edu/static/models/iIT341.json#Helicobacter pylori 26695———UniProt information not enough
wget http://bigg.ucsd.edu/static/models/iJB785.json#Synechococcus elongatus PCC 7942
wget http://bigg.ucsd.edu/static/models/iJN1463.json#Pseudomonas putida KT2440
wget http://bigg.ucsd.edu/static/models/iLB1027_lipid.json#	Phaeodactylum tricornutum CCAP 1055/1
wget http://bigg.ucsd.edu/static/models/iLJ478.json#Thermotoga maritima MSB8 
wget http://bigg.ucsd.edu/static/models/iMM904.json#Saccharomyces cerevisiae S288C
wget http://bigg.ucsd.edu/static/models/iNF517.json#Lactococcus lactis subsp. cremoris MG1363
wget http://bigg.ucsd.edu/static/models/iPC815.json#Yersinia pestis CO92
wget http://bigg.ucsd.edu/static/models/iRC1080.json#Chlamydomonas reinhardtii——can not found UniProt information
wget http://bigg.ucsd.edu/static/models/iSSON_1240.json#Shigella sonnei Ss046
wget http://bigg.ucsd.edu/static/models/iSynCJ816.json#Synechocystis sp. PCC 6803
wget http://bigg.ucsd.edu/static/models/iYL1228.json#Klebsiella pneumoniae subsp. pneumoniae MGH 78578
wget http://bigg.ucsd.edu/static/models/iYO844.json#Bacillus subtilis subsp. subtilis str. 168         
wget http://bigg.ucsd.edu/static/models/STM_v1_0.json#Salmonella enterica subsp. enterica serovar Typhimurium str. LT2
wget http://bigg.ucsd.edu/static/models/iYS854.json#Staphylococcus aureus subsp. aureus USA300_TCH1516——can not found UniProt information

```

### Model add UniProt ID


```python
import glob
import subprocess
import cobra
import re
import sys
sys.path.append(r'./script/')
from uniprot_id_mapping import *
#iLJ478
files = glob.glob('/hpcfs/fproject/mao_zt/MCModel/ECMpy/data/BIGG_model/*.json') #'iSynCJ816' 'iHN637' 

for eachf in files:
    submit_ids_list = []
    print(eachf)
    model_has_uniprot = False
    if re.search('\.xml',eachf):
        model = cobra.io.read_sbml_model(eachf)
    elif re.search('\.json',eachf):
        model = cobra.io.json.load_json_model(eachf)

    model_name = eachf.split('/')[-1].split('.json')[0]

    for eachg in model.genes:
        try:
            eachg.annotation
        except:
            print('Model does not have annotation!')
            submit_ids_list.append(eachg.id)
        else:
            if 'uniprot' in eachg.annotation.keys():
                model.genes.get_by_id(eachg.id).annotation['uniprot'] = model.genes.get_by_id(eachg.id).annotation['uniprot'][0]
                model_has_uniprot = True
            else:
                #print(eachg.notes)
                try:
                    eachg.notes['original_bigg_ids'][0]
                except:
                    print('Model can not find UniProt ID!')
                else:
                    if model_name=='iYO844':
                        submit_ids_list.append(eachg.id)
                    elif model_name=='iNF517' or model_name=='iJB785' or model_name=='iSynCJ816' or model_name=='iHN637' or model_name=='iLB1027_lipid':
                        submit_ids_list.append(eachg.notes['original_bigg_ids'][0].replace('-','_'))      
                        eachg.id = eachg.notes['original_bigg_ids'][0].replace('-','_')
                    else:
                        submit_ids_list.append(eachg.notes['original_bigg_ids'][0].replace('-','_'))
    submit_ids_list = list(set(submit_ids_list))
    if model_name=='iNF517' or model_name=='iJB785' or model_name=='iSynCJ816' or model_name=='iHN637' or model_name=='iLB1027_lipid':
        json_file_path_tmp = '/hpcfs/fproject/mao_zt/MCModel/ECMpy/data/BIGG_model/%s_change.json'%model_name
        cobra.io.save_json_model(model, json_file_path_tmp)
        model = cobra.io.json.load_json_model(json_file_path_tmp)
               
    #print(len(submit_ids_list))
    if model_has_uniprot:
        json_file_path = '/hpcfs/fproject/mao_zt/MCModel/ECMpy/data/BIGG_model/'+model_name+ "_uniprot.json"
        cobra.io.save_json_model(model, json_file_path)
    else:
        #gene to UniProtKB
        job_id = submit_id_mapping(from_db="Gene_Name", to_db="UniProtKB", ids=submit_ids_list)
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)

            for eachr in results['results']:
                #print(eachr['from'],eachr['to']['primaryAccession'])

                try:
                    eachr['to']['primaryAccession']
                except:
                    #pass
                    print(eachr['from'],eachr['to']['primaryAccession'])
                else:
                    try:
                        model.genes.get_by_id(eachr['from'])
                    except:
                        #iLJ478
                        try:
                            model.genes.get_by_id(eachr['from'].replace('_',''))
                        except:
                            eachr['from']
                        else:
                            model.genes.get_by_id(eachr['from'].replace('_','')).annotation['uniprot'] = eachr['to']['primaryAccession']
                    else:
                        model.genes.get_by_id(eachr['from']).annotation['uniprot'] = eachr['to']['primaryAccession']
        
        json_file_path = '/hpcfs/fproject/mao_zt/MCModel/ECMpy/data/BIGG_model/'+model_name+ "_uniprot.json"
        cobra.io.save_json_model(model, json_file_path)
    #break
    
    #iMM904_uniprot need manual correction GPR 
```

    /hpcfs/fproject/mao_zt/MCModel/ECMpy/data/BIGG_model/iYO844.json
    Retrying in 3s
    Fetched: 500 / 842
    Fetched: 842 / 842



```python
model_dict = {
    'iAF692_uniprot':'Methanosarcina barkeri',
    'iAF987_uniprot':'Geobacter metallireducens',
    'iAM_Pf480_uniprot':'Plasmodium falciparum',
    'iCN900_uniprot':'Clostridioides difficile',    
    'iEK1008_uniprot':'Mycobacterium tuberculosis',
    'iHN637_uniprot':'Clostridium ljungdahlii',
    'iIS312_uniprot':'Trypanosoma cruzi',
    'iJB785_uniprot':'Synechococcus elongatus',
    'iJN1463_uniprot':'Pseudomonas putida',
    'iLB1027_lipid_uniprot':'Phaeodactylum tricornutum',
    'iLJ478_uniprot':'Thermotoga maritima',
    'iMM904_uniprot':'Saccharomyces cerevisiae',
    'iNF517_uniprot':'Lactococcus lactis',
    'iPC815_uniprot':'Yersinia pestis',
    'iSSON_1240_uniprot':'Shigella sonnei',
    'iSynCJ816_uniprot':'Synechocystis sp.',
    'iYL1228_uniprot':'Klebsiella pneumoniae',
    'iYO844_uniprot':'Bacillus subtilis',
    'STM_v1_0_uniprot':'Salmonella enterica' 
}
```


```python
import glob
import subprocess
import re
import datetime 

files = glob.glob('/hpcfs/fproject/mao_zt/MCModel/ECMpy/data/BIGG_model/*_uniprot.json')#_uniprot
bigg_models_metabolites_file = '/hpcfs/fproject/mao_zt/MCModel/ECMpy/data/bigg_models_metabolites.txt'
brenda_file = '/hpcfs/fproject/mao_zt/MCModel/ECMpy/data/brenda_2023_1.txt'
uniprot_file = '/hpcfs/fproject/mao_zt/MCModel/ECMpy/data/uniprot_data_accession_key.json'

for eachf in files:
    print(eachf)
    if re.search('.json',eachf):
        model_name = eachf.split('/')[-1].split('.json')[0]
    elif re.search('.xml',eachf):
        model_name = eachf.split('/')[-1].split('.xml')[0]        
    org_name = "'"+model_dict[model_name]+"'"

    #AutoPACMAN
    ecGEM_file ='./model/BiGG/ec%s_AutoPACMEN.json'%model_name
    work_folder = '/hpcfs/fproject/mao_zt/MCModel/ECMpy/analysis/BiGG/get_kcat_mw_for_%s'%model_name
    cmd_str = "python /hpcfs/fproject/mao_zt/MCModel/ECMpy/script/get_ecGEM_onestop.py -m %s -kcat 'No' -f 0.45 -bigg %s -org %s -sigma 0.5 -ptot 0.56 -kcat_method 'AutoPACMEN' -work_folder %s -brenda %s -uniprot %s -kcat_gap_fill 'mean' -r_gap_fill 'mean' -ecGEM %s" %(eachf,bigg_models_metabolites_file,org_name,work_folder,brenda_file,uniprot_file,ecGEM_file)
    try:
        starttime=datetime.datetime.now()
        subprocess.run(cmd_str, shell=True)
        endtime=datetime.datetime.now()
        print(endtime-starttime)
    except:
        print('Can not construct AutoPACMAN ecGEM for %s!'%model_name)

#'''
    #DLKcat
    ecGEM_file ='./model/BiGG/ec%s_DLKcat.json'%model_name
    work_folder = '/hpcfs/fproject/mao_zt/MCModel/ECMpy/analysis/BiGG/get_kcat_mw_for_%s'%model_name
    cmd_str = "python /hpcfs/fproject/mao_zt/MCModel/ECMpy/script/get_ecGEM_onestop.py -m %s -kcat 'No' -f 0.45 -bigg %s -sigma 0.5 -ptot 0.56 -kcat_method 'DLKcat' -work_folder %s -ecGEM %s"%(eachf,bigg_models_metabolites_file,work_folder,ecGEM_file)

    try:
        starttime=datetime.datetime.now()
        subprocess.run(cmd_str, shell=True)
        endtime=datetime.datetime.now()
        print(endtime-starttime)
    except:
        print('Can not construct DLKcat ecGEM for %s!'%model_name)  

#'''

    #break
```


```python

```

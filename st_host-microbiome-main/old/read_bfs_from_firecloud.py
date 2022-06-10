#!/usr/bin/env python
# coding: utf-8

# In[13]:


import os
import pandas as pd
import numpy as np
from numpy import inf
import sys


# In[14]:


def filter_BFs(st_spec):
    
    """
    A function that filters the BF files from cloud
    # 
    # Inputs:
    #    st_spec             - BF data frame
    #    ST_top_gene_dict   - A pd.DataFrame with fields: age_1, age_2, region_1, region_2, AAR1, AAR2, logsBFs (list), Delta (list), genes (list)  


    """
    print(st_spec.head())
    # do some renaming
    st_spec = st_spec[st_spec['BF'] != inf]

    # Log10 BF
    st_spec['BF'] = np.float64(st_spec['BF'])
    x = np.log(st_spec['BF'])
    x[x == -inf] = sys.float_info.min # makes sure no inf
    x[x == inf] = sys.float_info.max # makes sure no inf
    st_spec['logBF'] = x
    print(st_spec['condition_1'].unique())
    print(st_spec['condition_2'].unique())
    # rename gene names
    st_spec['gene_new'] = [i.split("_")[0] for i in st_spec['gene']]
    st_spec['mouse_1'] = [i.split(" ")[0] for i in st_spec['condition_1']]
    st_spec['type_1'] = [i.split(" ")[1] for i in st_spec['condition_1']]
    st_spec['mouse_2'] = [i.split(" ")[0] for i in st_spec['condition_2']]
    st_spec['type_2'] = [i.split(" ")[1] for i in st_spec['condition_2']]
    print('done clean up')
    ## Top 100 ST genes per condition and per region
    ST_top_gene_dict = pd.DataFrame(columns = ['mouse_1', 'mouse_2', 'type_1', 'type_2', 'AAR1', 'AAR2', 'genes', 'logBFs', 'Delta'])
    counter = 0
    df_group = st_spec.groupby(['mouse_1', 'mouse_2', 'type_1', 'type_2', 'AAR1', 'AAR2'])
    for label, dfs in df_group: # this is for splotch_one_level

        # this gets genes super specific against the whole rest of the datset

        #if (label[5] == 'Rest'):

        #print(counter)

        #dfs = df[(df['logBF'] > 2) & (df['Delta'] > 0)]
        #dfs = df[(df['Delta'] > 0)]
        #dfs = df
        #print(df.sort_values(by='logBF', ascending=False)['gene_new'].head(5).tolist())
        if (len(dfs.sort_values(by='logBF', ascending=False)['gene_new'].head(5).tolist()) == 0):
            continue

        ST_top_gene_dict.at[counter, 'mouse_1'] = label[0]
        ST_top_gene_dict.at[counter, 'mouse_2'] = label[1]
        ST_top_gene_dict.at[counter, 'type_1'] = label[2]
        ST_top_gene_dict.at[counter, 'type_2'] = label[3]
        ST_top_gene_dict.at[counter, 'AAR1'] = label[4]
        ST_top_gene_dict.at[counter, 'AAR2'] = label[5]
        ST_top_gene_dict.at[counter, 'genes'] = dfs.sort_values(by=['logBF', 'Delta'], ascending=[False, False])['gene_new'].head(250).tolist()
        ST_top_gene_dict.at[counter, 'logBFs'] = dfs.sort_values(by=['logBF', 'Delta'], ascending=[False, False])['logBF'].head(250).tolist()
        ST_top_gene_dict.at[counter, 'Delta'] = dfs.sort_values(by=['logBF', 'Delta'], ascending=[False, False])['Delta'].head(250).tolist()
        counter += 1
    
    return ST_top_gene_dict


# In[ ]:


#%%capture cap_out --no-stderr
#!gsutil ls gs://fc-2e7e9da8-0b98-4a74-a1ed-61c5ba5a13a4/3c779645-2301-4389-a0d3-da5f554e148a/splotch_diff_exp_workflow/a3bf8b06-c212-49b7-883d-8a96072c8a02/call-diff_exp/ 


# In[16]:


#get_ipython().run_cell_magic('capture', 'cap_out --no-stderr', '!gsutil ls 
#gs://fc-2e7e9da8-0b98-4a74-a1ed-61c5ba5a13a4/14bf41a0-abca-4bb9-890e-3d482f180aee/splotch_diff_exp_workflow/1149cc6e-b6d9-47aa-9862-83ac96039c61/call-diff_exp/ ')



# In[17]:


shards = [i for i in cap_out.stdout.split("\r\n")]


# In[5]:


shards


# In[ ]:


path = '/home/brittalotstedt/host-microbiome/data/bfs'

df_bf = pd.DataFrame([])
for shard in shards:
    
    if shard.startswith('gs://fc'):
        
    #check file names
        filename = shard+'input_dir/analysis_output/BF-beta_df_output.tsv' #Try cacheCopy/input_dir/... if problematic 
        
        print(filename)
        output_filename = os.path.join(path, 'bf.tsv')

        # copy file due to access permissions
        get_ipython().system('gsutil -m cp $filename $output_filename')
        get_ipython().system('gzip $output_filename')

        # read in as pandas and filter
        bf_tmp = pd.read_csv(output_filename+".gz", compression='gzip', header=0, sep='\t', quotechar='"', error_bad_lines=False, index_col=0)
        bf_tmp['BF'] = bf_tmp['BF'].astype(float)
        bf_tmp['Delta'] = bf_tmp['Delta'].astype(float)
        bf_tmp = bf_tmp[(bf_tmp['BF']>0) & (bf_tmp['Delta']>0)]

        # append to larger df
        df_bf = df_bf.append(bf_tmp,ignore_index=True)

        # delete file from memory 
        output_filename = output_filename+".gz"
        get_ipython().system('rm $output_filename')


# In[7]:


ST_top_gene_dict = filter_BFs(df_bf)

ST_top_gene_dict.to_csv(os.path.join(path, 'ST_top_gene_dict_BF0_longnames.csv'))
# In[9]:


shorten_anns_dict = {"muscularis mucosae and muscularis propria muscularis interna and peyer's patch":"mucosae and peyer's patch",
               "epithelium and peyer's patch":"peyer's patch",
                  'epithelium and lamina propria' : 'epithelium',
    'epithelium and lamina propria and mucosa' : 'epithelium and mucosa',
    'epithelium and lamina propria and mucosa and pellet': 'epithelium and mucosa and pellet',
    'epithelium and lamina propria and muscularis mucosae' : 'epithelium and mucosae',
    'epithelium and lamina propria and muscularis mucosae and muscularis propria muscularis externa and muscularis propria muscularis interna and submucosa all':'epithelium and muscle and submucosa',
    'epithelium and lamina propria and muscularis mucosae and submucosa all':'epithelium and mucosae and submucosa',
    "epithelium and lamina propria and peyer's patch" : "peyer's patch",
    'epithelium apex of crypt and lamina propria' : 'crypt apex',
    'epithelium base of crypt and epithelium mid crypt and lamina propria' : 'crypt base and mid',
    'epithelium base of crypt and lamina propria' : 'crypt base',
    'epithelium base of crypt and lamina propria and muscularis mucosae and muscularis propria muscularis externa and muscularis propria muscularis interna and submucosa all' : 'crypt base',
    'epithelium base of crypt and muscularis mucosae and muscularis propria muscularis interna': 'crypt base',
    'epithelium mid crypt and lamina propria' : 'crypt mid',
    'mucosa and pellet' : 'mucosa and pellet',
    'muscularis all and submucosa all' : 'muscle and submucosa',
    'muscularis mucosae and muscularis propria muscularis externa and muscularis propria muscularis interna' : 'muscle and submucosa',
    'muscularis mucosae and muscularis propria muscularis interna' : 'mucosae and interna',
    "muscularis mucosae and muscularis propria muscularis interna and peyer's patch" : "peyer's patch",
    "muscularis mucosae and peyer's patch" : "peyer's patch",
    'muscularis propria muscularis externa' : 'externa',
    'muscularis propria muscularis externa and muscularis propria muscularis interna':'externa and interna',
    'muscularis propria muscularis interna':'interna',
    'muscularis propria muscularis externa and muscularis propria muscularis interna and muscularis mucosae and submucosa all' : 'muscle and submucosa',
    'pellet':'pellet',
    "peyer's patch":"peyer's patch",
    'mucosa':'mucosa',
    'epithelium apex of crypt and mucosa':'crypt apex and mucosa',
    "Rest": "rest"}

expand_ages_dict = {
    '0d' : '0w',
    '1w' : '1w',
    '2w' : '2w',
    '3w' : '3w',
    '4w' : '4w',
    '6w' : '6w',
    '8w' : '8w',
    '12w' : '12w',
    '6m' : '26w',
    '1yr' : '52w',
    '2yr' : '104w'
}


# In[10]:


print([i for i in ST_top_gene_dict.AAR1.unique() if i not in shorten_anns_dict])
print([i for i in ST_top_gene_dict.AAR2.unique() if i not in shorten_anns_dict])


# In[11]:


ST_top_gene_dict['AAR1'] = ST_top_gene_dict['AAR1'].map(shorten_anns_dict)
ST_top_gene_dict['AAR2'] = ST_top_gene_dict['AAR2'].map(shorten_anns_dict)


# In[ ]:


#ST_top_gene_dict['region_1'] = [i.replace("BL6WT.", "") for i in ST_top_gene_dict['region_1']]
#ST_top_gene_dict['region_2'] = [i.replace("BL6WT.", "") for i in ST_top_gene_dict['region_2']]


# In[ ]:


#ST_top_gene_dict['mouse_1'] = ST_top_gene_dict['age_1'].map(expand_ages_dict)
#ST_top_gene_dict['age_2'] = ST_top_gene_dict['age_2'].map(expand_ages_dict)


# In[12]:


# Saves formated DE genes df

ST_top_gene_dict.to_csv(os.path.join(path, 'ST_top_gene_dict_BF0.csv'))


# In[ ]:





# In[ ]:





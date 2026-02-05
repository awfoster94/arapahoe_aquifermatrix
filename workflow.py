# this workflow sets up pestpp-ies to run foward models of phreeqc to simulate downstream water quality changes from varying geochemical reactions inputs along a specified flowpath

# import necessary libraries
import os
import stat
import re
import pandas as pd
import numpy as np
import pyemu
import sys
# sys.path.insert(0,os.path.abspath(os.path.join('dependencies')))
# sys.path.insert(1,os.path.abspath(os.path.join('dependencies','flopy')))
import platform
import shutil
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import random
from random import randint
import time
import geopandas as gpd
import warnings
import math
import zipfile
import subprocess
import glob
import psutil
warnings.filterwarnings("ignore")
import datetime
import os
import shutil
import pandas as pd
import pyemu
import builtins
#_orig = builtins.open
#builtins.open = lambda fn, mode='r', *a, **kw: _orig(fn, mode, encoding='utf-8', *a, **kw)
#import write_condor


### this function prepares dependency folders & files
def prep_deps(d):
    # copy in deps and exes
    print(d)
    if "window" in platform.platform().lower():
        bd = os.path.join("bin")
    elif "linux" in platform.platform().lower():
        bd = os.path.join("bin")
    else:
        bd = os.path.join("bin")
    for f in os.listdir(bd):
        shutil.copy(os.path.join(bd, f), os.path.join(d, f))
    #try:
    #    shutil.rmtree(os.path.join(d, "flopy"))
    #except:
    #    pass
    #shutil.copytree(os.path.join('dependencies', 'flopy'), os.path.join(d, "flopy"), copy_function=shutil.copy)
    try:
        shutil.rmtree(os.path.join(d, "pyemu"))
    except:
        pass
    shutil.copytree(os.path.join('dependencies', "pyemu"), os.path.join(d, "pyemu"), copy_function=shutil.copy)
	
### this function checks the mult2model_info from pestpp-ies
def check_mult2model_info(template_d):
    bd = os.getcwd()
    os.chdir(template_d)
    # load your array‐parameter spec
    arr = pd.read_csv("mult2model_info.csv")
    mylist = []
    for fn in arr["model_file"].unique():
        # how many parameters PEST expects for this file?
        par = len(arr[arr["model_file"] == fn])
        # how many lines does .ref actually have?
        try:
            with open(fn, "r") as f:
                length = sum(1 for _ in f)
        except FileNotFoundError:
            continue
        mylist.append((fn, par, length))
    for fn, par, length in mylist:
        print(f"{fn}: adj parameters={par}, length of file:{length}")
    os.chdir(bd)

### copy original folder to template folder
def copy_folder(ord_d, temp_d):
    os.makedirs(temp_d, exist_ok=True)

    # copy original directory to template directory
    for filename in os.listdir(org_d):
        org_d_fpth = os.path.join(org_d, filename)
        temp_d_fpth = os.path.join(temp_d, filename)
    
        # copy all the files to the new folder
        if os.path.isfile(org_d_fpth):
            shutil.copy(org_d_fpth, temp_d_fpth)

### this function writes the initial template model file for the forward run
def test_write_input_react_chain(d, param_list):
    
    param_array = []
    for param in param_list:
        print(param)
        param_val = pd.read_csv(os.path.join(d,str(param)+'.txt'), header=None)
        param_array.append([param,param_val])
    #print(param_array)
    param_array_df = pd.DataFrame(param_array, columns=['param', 'value'])
    #test  = param_array_df[param_array_df['param'] == 'dolomite']
    #print(test)

    #test_string = np.array2string(test['value'][1].values, separator=' ', prefix='', suffix='')
    #print(test_string)

    text = f"""
            TITLE Forward Model No. 1 Arapahoe

            KNOBS
                -iterations 1000
                
            EQUILIBRIUM_PHASES 1
                    Calcite	0.0 0.0
                    Dolomite 0.0 0.0
                    Gypsum 0.0 0.0	
                    Strontianite 0.0 0.0 dissolve_only
                    Witherite 0.0 0.0 dissolve_only
                    CO2(g) 0.0 0.0
                    CaX2 0.0 0.0 
                    MgX2 0.0 0.0 
                    KX 0.0 0.0 
                    NaX 0.0 0.0 
                    Quartz 0.0 0.0 dissolve_only
                    Albite 0.0 0.0	dissolve_only
                    Kaolinite 0.0 0.0 dissolve_only
                    Muscovite 0.0 0.0
                    Uraninite 0.0 0.0
                    Arsenopyrite 0.0 0.0 dissolve_only
                    Goethite 0.0 0.0
                    Fe(OH)3 0.0 0.0 
                    Siderite 0.0 0.0
                    Mn(OH)3 0.0 0.0			
                    Fluorapatite 0.0 0.0 
                    Halite 0.0 0.0
                    Illite 0.0 0.0
                    Chromite 0.0 0.0	dissolve_only
            END
                    
            REACTION_TEMPERATURE 1
                25.4
            END
                
            Solution 1 # # USGS-392118104362301 Upstream WQ
                    units		mol/kgw
                    redox       pe
                    temp		30.1
                    pH		    8.3
                    O(0)        0.0
                    Na		    0.003062228
                    Ca		    5.81337E-05
                    Mg		    1.67867E-05
                    K           3.91321E-05
                    Alkalinity	0.002817521	as CaCO3
                    S(6)		0.000108263
                    Cl		    4.62623E-05
                    F		    8.57967E-05
                    Mn		    1.98405E-07
                    Ba          2.66517E-07
                    Sr          5.14723E-07
                    Fe          4.7093E-07
                    Al          7.41235E-08
                    Si          0.000150956 as SiO2
                    As          8.00837E-10
                    U           8.40234E-11
                    Cr          5.76967E-10
                    P           1.16227E-06    as P
                
                water 1

            SAVE Solution 1
            END

            REACTION 1
                    Calcite	{str(param_array_df[param_array_df['param']=='calcite']['value'][0].values)} 
                    Dolomite {str(param_array_df[param_array_df['param']=='dolomite']['value'][1].values)}
                    Gypsum {str(param_array_df[param_array_df['param']=='gypsum']['value'][2].values)}
                    Strontianite {str(param_array_df[param_array_df['param']=='strontianite']['value'][3].values)}
                    Witherite {str(param_array_df[param_array_df['param']=='witherite']['value'][4].values)}
                    CO2(g) {str(param_array_df[param_array_df['param']=='co2g']['value'][5].values)}
                    CaX2 {str(param_array_df[param_array_df['param']=='cax2']['value'][6].values)}
                    MgX2 {str(param_array_df[param_array_df['param']=='mgx2']['value'][7].values)} 
                    KX {str(param_array_df[param_array_df['param']=='kx']['value'][8].values)}
                    NaX {str(param_array_df[param_array_df['param']=='nax']['value'][9].values)}
                    Quartz {str(param_array_df[param_array_df['param']=='quartz']['value'][10].values)}
                    Albite {str(param_array_df[param_array_df['param']=='albite']['value'][11].values)}
                    Kaolinite {str(param_array_df[param_array_df['param']=='kaolinite']['value'][12].values)}
                    Muscovite {str(param_array_df[param_array_df['param']=='muscovite']['value'][13].values)}
                    Uraninite {str(param_array_df[param_array_df['param']=='uraninite']['value'][14].values)}
                    Arsenopyrite {str(param_array_df[param_array_df['param']=='arsenopyrite']['value'][15].values)}
                    Goethite {str(param_array_df[param_array_df['param']=='goethite']['value'][16].values)}
                    Fe(OH)3 {str(param_array_df[param_array_df['param']=='feoh3']['value'][17].values)}
                    Siderite {str(param_array_df[param_array_df['param']=='siderite']['value'][18].values)}
                    Mn(OH)3 {str(param_array_df[param_array_df['param']=='mnoh3']['value'][19].values)}
                    Fluorapatite {str(param_array_df[param_array_df['param']=='fluorapatite']['value'][20].values)}
                    Halite {str(param_array_df[param_array_df['param']=='halite']['value'][21].values)}
                    Illite {str(param_array_df[param_array_df['param']=='illite']['value'][22].values)}
                    Chromite {str(param_array_df[param_array_df['param']=='chromite']['value'][23].values)}
            END

            Solution 2 # USGS-394155104425401 Downstream WQ
                    units		mol/kgw
                    redox       pe
                    temp		25.4
                    pH		    8.7
                    O(0)        9.37559E-06
                    Na		    0.002457612
                    Ca		    0.000183383
                    Mg		    2.03662E-05
                    K           4.45032E-05
                    Alkalinity	0.002437855	as CaCO3
                    S(6)		0.00020106
                    Cl		    5.72638E-05
                    F		    9.2113E-05
                    Mn		    2.56653E-07
                    Ba          5.77454E-07
                    Sr          1.60922E-06
                    Fe          1.52739E-06
                    Al          2.96494E-08
                    Si          0.000131407 as SiO2
                    As          1.33473E-09
                    U           8.40234E-11
                    Cr          7.6929E-09
                    P           9.68561E-08    as P
                
                water 1

            SAVE Solution 2
            END

            TITLE Mix Upstream and Downstream Water Qaulity
            MIX 1 
            1 {str(param_array_df[param_array_df['param']=='mix-ratio-sol-upstream']['value'][24].values)}
            2 {str(param_array_df[param_array_df['param']=='mix-ratio-sol-downstream']['value'][25].values)}
            SAVE Solution 3 

            END

            TITLE REACT, EQUILIBRATE WITH SUBSURFACE MATERIALS
            USE Solution 3
            USE REACTION_TEMPERATURE 1
            USE REACTION 1
            USE EQUILIBRIUM_PHASES 1

            SAVE Solution 4

            SELECTED_OUTPUT
                file Arapahoe_01_template_selected_output.txt
                selected_out      TRUE
                high_precision     TRUE
                pH     TRUE
                charge_balance     TRUE
                percent_error     TRUE
                molalities      HCO3- CO3-2
                totals     O(0) Cl Ca Mg Na K S(6) Alkalinity F Mn Ba Sr Fe As U P Al Si Cr
                -saturation_indices		Aragonite Brucite Calcite Dolomite Dolomite-dis Fe(OH)3 Goethite Gypsum Gibbsite Magnesite Mn(OH)3 Rhodochrosite Siderite Sylvite Strontianite Witherite CO2(g) O2(g)
                -equilibrium_phases      Aragonite Brucite Calcite Dolomite Dolomite-dis Fe(OH)3 Goethite Gypsum Gibbsite Magnesite Mn(OH)3 Rhodochrosite Siderite Sylvite Strontianite Witherite  CO2(g) O2(g)
            END"""

    
    ### write input file to folder
    print(text)

    ### update text file to remove brackets
    updated_text = text.replace('[', '').replace(']', '')
    
    ### output phreeqc input file
    with open(os.path.join(d,"Arapahoe_01_template.txt"), "w") as file:
        file.write(updated_text)

# fxn to write the initial template model file for the forward run
def write_input_react_chain():
    param_list = ['calcite', 'dolomite', 'gypsum', 'strontianite', 'witherite', 'co2g', 
                  'cax2', 'mgx2', 'kx', 'nax', 'quartz', 'albite', 'kaolinite', 'muscovite',
                  'uraninite', 'arsenopyrite', 'goethite', 'feoh3', 'siderite', 'mnoh3',
                  'fluorapatite', 'halite', 'illite', 'chromite', 'mix-ratio-sol-upstream', 'mix-ratio-sol-downstream'
                  ]
    
    param_array = []
    for param in param_list:
        print(param)
        param_val = pd.read_csv(str(param)+'.txt', header=None)
        param_array.append([param,param_val])
    #print(param_array)
    param_array_df = pd.DataFrame(param_array, columns=['param', 'value'])
    #test  = param_array_df[param_array_df['param'] == 'dolomite']
    #print(test)

    #test_string = np.array2string(test['value'][1].values, separator=' ', prefix='', suffix='')
    #print(test_string)

    text = f"""
            TITLE Forward Model No. 1 Arapahoe

            KNOBS
                -iterations 1000
                
            EQUILIBRIUM_PHASES 1
                    Calcite	0.0 0.0
                    Dolomite 0.0 0.0
                    Gypsum 0.0 0.0	
                    Strontianite 0.0 0.0 dissolve_only
                    Witherite 0.0 0.0 dissolve_only
                    CO2(g) 0.0 0.0
                    CaX2 0.0 0.0 
                    MgX2 0.0 0.0 
                    KX 0.0 0.0 
                    NaX 0.0 0.0 
                    Quartz 0.0 0.0 dissolve_only
                    Albite 0.0 0.0	dissolve_only
                    Kaolinite 0.0 0.0 dissolve_only
                    Muscovite 0.0 0.0
                    Uraninite 0.0 0.0
                    Arsenopyrite 0.0 0.0 dissolve_only
                    Goethite 0.0 0.0
                    Fe(OH)3 0.0 0.0 
                    Siderite 0.0 0.0
                    Mn(OH)3 0.0 0.0			
                    Fluorapatite 0.0 0.0 
                    Halite 0.0 0.0
                    Illite 0.0 0.0
                    Chromite 0.0 0.0	dissolve_only
            END
                    
            REACTION_TEMPERATURE 1
                25.4
            END
                
            Solution 1 # # USGS-392118104362301 Upstream WQ
                    units		mol/kgw
                    redox       pe
                    temp		30.1
                    pH		    8.3
                    O(0)        0.0
                    Na		    0.003062228
                    Ca		    5.81337E-05
                    Mg		    1.67867E-05
                    K           3.91321E-05
                    Alkalinity	0.002817521	as CaCO3
                    S(6)		0.000108263
                    Cl		    4.62623E-05
                    F		    8.57967E-05
                    Mn		    1.98405E-07
                    Ba          2.66517E-07
                    Sr          5.14723E-07
                    Fe          4.7093E-07
                    Al          7.41235E-08
                    Si          0.000150956 as SiO2
                    As          8.00837E-10
                    U           8.40234E-11
                    Cr          5.76967E-10
                    P           1.16227E-06    as P
                
                water 1

            SAVE Solution 1
            END

            REACTION 1
                    Calcite	{str(param_array_df[param_array_df['param']=='calcite']['value'][0].values)} 
                    Dolomite {str(param_array_df[param_array_df['param']=='dolomite']['value'][1].values)}
                    Gypsum {str(param_array_df[param_array_df['param']=='gypsum']['value'][2].values)}
                    Strontianite {str(param_array_df[param_array_df['param']=='strontianite']['value'][3].values)}
                    Witherite {str(param_array_df[param_array_df['param']=='witherite']['value'][4].values)}
                    CO2(g) {str(param_array_df[param_array_df['param']=='co2g']['value'][5].values)}
                    CaX2 {str(param_array_df[param_array_df['param']=='cax2']['value'][6].values)}
                    MgX2 {str(param_array_df[param_array_df['param']=='mgx2']['value'][7].values)} 
                    KX {str(param_array_df[param_array_df['param']=='kx']['value'][8].values)}
                    NaX {str(param_array_df[param_array_df['param']=='nax']['value'][9].values)}
                    Quartz {str(param_array_df[param_array_df['param']=='quartz']['value'][10].values)}
                    Albite {str(param_array_df[param_array_df['param']=='albite']['value'][11].values)}
                    Kaolinite {str(param_array_df[param_array_df['param']=='kaolinite']['value'][12].values)}
                    Muscovite {str(param_array_df[param_array_df['param']=='muscovite']['value'][13].values)}
                    Uraninite {str(param_array_df[param_array_df['param']=='uraninite']['value'][14].values)}
                    Arsenopyrite {str(param_array_df[param_array_df['param']=='arsenopyrite']['value'][15].values)}
                    Goethite {str(param_array_df[param_array_df['param']=='goethite']['value'][16].values)}
                    Fe(OH)3 {str(param_array_df[param_array_df['param']=='feoh3']['value'][17].values)}
                    Siderite {str(param_array_df[param_array_df['param']=='siderite']['value'][18].values)}
                    Mn(OH)3 {str(param_array_df[param_array_df['param']=='mnoh3']['value'][19].values)}
                    Fluorapatite {str(param_array_df[param_array_df['param']=='fluorapatite']['value'][20].values)}
                    Halite {str(param_array_df[param_array_df['param']=='halite']['value'][21].values)}
                    Illite {str(param_array_df[param_array_df['param']=='illite']['value'][22].values)}
                    Chromite {str(param_array_df[param_array_df['param']=='chromite']['value'][23].values)}
            END

            Solution 2 # USGS-394155104425401 Downstream WQ
                    units		mol/kgw
                    redox       pe
                    temp		25.4
                    pH		    8.7
                    O(0)        9.37559E-06
                    Na		    0.002457612
                    Ca		    0.000183383
                    Mg		    2.03662E-05
                    K           4.45032E-05
                    Alkalinity	0.002437855	as CaCO3
                    S(6)		0.00020106
                    Cl		    5.72638E-05
                    F		    9.2113E-05
                    Mn		    2.56653E-07
                    Ba          5.77454E-07
                    Sr          1.60922E-06
                    Fe          1.52739E-06
                    Al          2.96494E-08
                    Si          0.000131407 as SiO2
                    As          1.33473E-09
                    U           8.40234E-11
                    Cr          7.6929E-09
                    P           9.68561E-08    as P
                
                water 1

            SAVE Solution 2
            END

            TITLE Mix Upstream and Downstream Water Quality
            MIX 1 
            1 {str(param_array_df[param_array_df['param']=='mix-ratio-sol-upstream']['value'][24].values)}
            2 {str(param_array_df[param_array_df['param']=='mix-ratio-sol-downstream']['value'][25].values)}
            SAVE Solution 3 

            END

            TITLE REACT, EQUILIBRATE WITH SUBSURFACE MATERIALS
            USE Solution 3
            USE REACTION_TEMPERATURE 1
            USE REACTION 1
            USE EQUILIBRIUM_PHASES 1

            SAVE Solution 4

            SELECTED_OUTPUT
                file Arapahoe_01_template_selected_output.txt
                selected_out      TRUE
                high_precision     TRUE
                pH     TRUE
                charge_balance     TRUE
                percent_error     TRUE
                molalities      HCO3- CO3-2
                totals     O(0) Cl Ca Mg Na K S(6) Alkalinity F Mn Ba Sr Fe As U P Al Si Cr
                -saturation_indices		Aragonite Brucite Calcite Dolomite Dolomite-dis Fe(OH)3 Goethite Gypsum Gibbsite Magnesite Mn(OH)3 Rhodochrosite Siderite Sylvite Strontianite Witherite CO2(g) O2(g)
                -equilibrium_phases      Aragonite Brucite Calcite Dolomite Dolomite-dis Fe(OH)3 Goethite Gypsum Gibbsite Magnesite Mn(OH)3 Rhodochrosite Siderite Sylvite Strontianite Witherite  CO2(g) O2(g)
            END"""

    
    ### write input file to folder
    print(text)

    ### update text file to remove brackets
    updated_text = text.replace('[', '').replace(']', '')
    
    ### output phreeqc input file
    with open("Arapahoe_01_template.txt", "w") as file:
        file.write(updated_text)

# fxn to run phreeqc
def run_phreeqc(temp_d, inp_fnm, out_fnm, therm_fnm):

    # Specify the path to phreeqc.exe
    phreeqc_path = os.path.join(temp_d, "phreeqc.exe")  # Replace with the actual relative path to phreeqc.exe

    # Build the command to execute
    command = [
        phreeqc_path,
        os.path.join(temp_d, inp_fnm),
        os.path.join(temp_d, out_fnm),
        os.path.join(temp_d, therm_fnm)
    ]

    try:
        # Use subprocess to run the command
        subprocess.run(command, check=True)
        print("Phreeqc executed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error while running Phreeqc: {e}")
        
    # Pause the code for 0.25 seconds
    time.sleep(0.25)

    print("This message is displayed after a 0.25-second pause.")

# fxn to curate the simulated values for the pest object
def build_simulated_pest_obs_test():
    # load in template for simulated data to update
    sim_data_temp = pd.read_csv(os.path.join('model_ws', 'Arapahoe-01', 'template', 'sim.csv'))
    param_list = sim_data_temp['constituent'].unique()
    updated_sim_data = sim_data_temp.copy()

    # load in the simulated data from selected output file
    sim_data = pd.read_csv(os.path.join('model_ws', 'Arapahoe-01', 'template','Arapahoe_01_template_selected_output.txt'), delimiter='\t')
    # Clean the column names: remove all whitespace (space, tabs, etc.)
    sim_data.columns = sim_data.columns.str.replace(r'\s+', '', regex=True)
    # Clean all string cell values (optional, if needed)
    sim_data = sim_data.applymap(lambda x: re.sub(r'\s+', '', x) if isinstance(x, str) else x)
    
    # Update values in the copied DataFrame
    for constituent in param_list:
        if constituent == 'pH':
            # Convert pH to [H+]
            h_plus = 10 ** (-sim_data['pH'].iloc[0])
            updated_sim_data.loc[updated_sim_data['constituent'] == 'pH', 'conc'] = h_plus

        elif constituent == 'Alkalinity':
            if 'm_HCO3-' in sim_data.columns and 'm_CO3-2' in sim_data.columns:
                alk_caco3 = 1 * sim_data['m_HCO3-'].iloc[0] + sim_data['m_CO3-2'].iloc[0]
                updated_sim_data.loc[updated_sim_data['constituent'] == 'Alkalinity', 'conc'] = alk_caco3

        elif constituent in sim_data.columns:
            # Default case: update directly from sim_data
            updated_sim_data.loc[updated_sim_data['constituent'] == constituent, 'conc'] = sim_data[constituent].iloc[0]
    
    updated_sim_data.to_csv(os.path.join('model_ws', 'Arapahoe-01', 'template', 'sim.csv'), index=False)
    return updated_sim_data

# fxn to load the observed downstream water quality data observations
def build_obs_test():
    obs_vals = pd.read_csv(os.path.join('model_ws', 'Arapahoe-01', 'template','obs.csv'))
    return obs_vals

# fxn to curate the simulated values for the pest object
def build_simulated_pest_obs():
    import re
    # load in template for simulated data to update
    sim_data_temp = pd.read_csv(os.path.join('sim_template.csv'))
    param_list = sim_data_temp['constituent'].unique()
    updated_sim_data = sim_data_temp.copy()

    # load in the simulated data from selected output file
    sim_data = pd.read_csv(os.path.join('Arapahoe_01_template_selected_output.txt'), delimiter='\t')
    # Clean the column names: remove all whitespace (space, tabs, etc.)
    sim_data.columns = sim_data.columns.str.replace(r'\s+', '', regex=True)
    # Clean all string cell values (optional, if needed)
    sim_data = sim_data.applymap(lambda x: re.sub(r'\s+', '', x) if isinstance(x, str) else x)
    
    # Update values in the copied DataFrame
    for constituent in param_list:
        if constituent == 'pH':
            # Convert pH to [H+]
            h_plus = 10 ** (-sim_data['pH'].iloc[0])
            updated_sim_data.loc[updated_sim_data['constituent'] == 'pH', 'conc_raw'] = h_plus

        elif constituent == 'Alkalinity':
            if 'm_HCO3-' in sim_data.columns and 'm_CO3-2' in sim_data.columns:
                alk_caco3 = 1 * sim_data['m_HCO3-'].iloc[0] + sim_data['m_CO3-2'].iloc[0]
                updated_sim_data.loc[updated_sim_data['constituent'] == 'Alkalinity', 'conc_raw'] = alk_caco3

        elif constituent in sim_data.columns:
            # Default case: update directly from sim_data
            updated_sim_data.loc[updated_sim_data['constituent'] == constituent, 'conc_raw'] = sim_data[constituent].iloc[0]
    
    # flip to conc_log10 for easier pest weighting rebalancing, idea from the Jeremy White...
    # Replace zeros with 1e-20
    updated_sim_data.loc[updated_sim_data['conc_raw'] == 0, 'conc_raw'] = 1e-20
    updated_sim_data['conc'] = -1*np.log10(updated_sim_data['conc_raw'])

    updated_sim_data.to_csv(os.path.join('sim.csv'), index=False)
    return updated_sim_data

### this function loads the observed downstream water quality data observations
def build_obs():
    obs_vals = pd.read_csv(os.path.join('obs.csv'))
    return obs_vals

### this function deletes old model output files
def delete_old_model_output():
    print('deleting model output...')
    """Delete old model output file so each realization starts fresh."""
    try:
        os.remove("Arapahoe_01_template.txt.out")
    except FileNotFoundError:
        pass
    except Exception as err:
        print(f"error removing Arapahoe_01_template.txt.out: {err}")

    try:
        os.remove("Arapahoe_01_template_selected_output.txt")
    except FileNotFoundError:
        pass
    except Exception as err:
        print(f"error removing Arapahoe_01_template_selected_output.txt: {err}")

    try:
        os.remove("Arapahoe_01_template.txt")
    except FileNotFoundError:
        pass
    except Exception as err:
        print(f"error removing Arapahoe_01_template.txt: {err}")

## this function sets up obs vals
def setup_obs_val(t_d, obs_df, set_weights_phi_facs_csv=0):
        
    # assign weights to all non-zero weight groups
    weight_mask = obs_df["obgnme"].str.startswith(('oname:sim_otype:lst_usecol:conc_tag:general.ph',
                                                   'oname:sim_otype:lst_usecol:conc_tag:majoranion.alkalinity',
                                                   'oname:sim_otype:lst_usecol:conc_tag:majoranion.cl',
                                                   'oname:sim_otype:lst_usecol:conc_tag:majoranion.s(6)',
                                                   'oname:sim_otype:lst_usecol:conc_tag:majorcation.ca',
                                                   'oname:sim_otype:lst_usecol:conc_tag:majorcation.k',
                                                   'oname:sim_otype:lst_usecol:conc_tag:majorcation.mg',
                                                   'oname:sim_otype:lst_usecol:conc_tag:majorcation.na',
                                                   'oname:sim_otype:lst_usecol:conc_tag:minoranion.f',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.al',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.as',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.ba',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.cr',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.fe',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.mn',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.p',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.si',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.sr',
                                                   'oname:sim_otype:lst_usecol:conc_tag:trace.u'), na=False)
    obs_df.loc[weight_mask, 'weight'] = 1.0
    # accuracy of analytical data assumed
    #decimal_percent = 0.025 # 2.5%
    #obs_df.loc[weight_mask, 'standard deviation'] = obs_df.loc[weight_mask, 'obsval'] * decimal_percent

    if set_weights_phi_facs_csv:
        print(f'copying in phi_fac.csv file to {t_d} with specified desired relative contributions for each observation group...')
        phi_factor_csv_srcfpth = os.path.join(os.getcwd(), 'calibration', 'phi_facs.csv')
        phi_factor_csv_dstfpth = os.path.join(os.getcwd(), f'{t_d}', 'phi_facs.csv')

        # copy phi_fac.csv to template_d
        shutil.copy(phi_factor_csv_srcfpth, phi_factor_csv_dstfpth)

    return obs_df


# fxn to add measurement noise...
def draw_noise_reals(m_d,modnm):
    
    # noise definition
    noise_percentage = 0.025

    # load pst object and obs data
    pst = pyemu.Pst(os.path.join(m_d,f"{modnm}.pst"))
    obs = pst.observation_data
    
    # add to weighted obs
    weighted_obs = obs.loc[obs.weight > 0.0, :]
    
    # pull ensemble size
    pr = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(m_d,pst.pestpp_options["ies_par_en"]))
    num_reals = pr.shape[0]

    # initialize observation ensemble
    #oe = pyemu.ObservationEnsemble.from_observation_data(pst=pst, num_reals=num_reals)
    obsvals = obs.obsval.copy()
    df = pd.DataFrame(np.tile(obsvals.values, (num_reals, 1)), columns=obsvals.index)
    oe = pyemu.ObservationEnsemble.from_dataframe(pst, df)
    
    obsvals = weighted_obs.obsval.values
    # 5% noise addition based on gaussian (normal) distribution
    sigma = np.maximum(noise_percentage * np.abs(obsvals), 1E-8)

    # define noise object
    noise = np.random.normal(loc=0.0, scale=sigma, size=(num_reals, len(weighted_obs)))

    oe.loc[:, weighted_obs.obsnme] = obsvals + noise

    # write the ensemble
    noise_fnm = "noise.jcb"
    oe.to_binary(os.path.join(m_d, noise_fnm))

    # update the pst object and control file
    pst.pestpp_options["ies_obs_en"] = noise_fnm
    #pst.control_data.noptmax = -2
    pst.write(os.path.join(m_d, f"{modnm}.pst"), version=2)


# setup pestpp ies
def setup_pstppies(od, template_d, reals, test=1, check_pstpp_setup=0, auto_weight_rebal=0):
    print("Setting up PESTPP workflow...")

    assert os.path.exists(od)
    temp_d = od + '_temp'
    if os.path.exists(temp_d):
        shutil.rmtree(temp_d)
    shutil.copytree(od, temp_d, copy_function=shutil.copy)
    print(f"Copying: {od}, new dir: {temp_d}")

    prep_deps(temp_d)
    bd = os.path.join("bin")
    for f in os.listdir(bd):
        shutil.copy(os.path.join(bd, f), os.path.join(temp_d, f))

    inp_fnm = "Arapahoe_01_template.txt"
    out_fnm = inp_fnm + ".out"
    therm_fnm = "llnl.dat"
    run_phreeqc(temp_d, inp_fnm, out_fnm, therm_fnm)

    print(f"New PEST directory: {template_d}")
    pf = pyemu.utils.PstFrom(
        original_d=temp_d,
        new_d=template_d,
        remove_existing=True,
        longnames=True,
        zero_based=False
    )

    # add helper functions for pre-processing
    pf.add_py_function('workflow.py', 'build_simulated_pest_obs()', is_pre_cmd=True)
    pf.add_py_function('workflow.py', 'delete_old_model_output()', is_pre_cmd=True)
    pf.add_py_function('workflow.py', 'write_input_react_chain()', is_pre_cmd=True)
    pf.add_py_function('workflow.py', 'build_obs()', is_pre_cmd=True)

    # fxn to write .ins files...
    def write_csv_ins_file(csv_path, ins_path, obs_names):
        with open(ins_path, 'w') as f:
            f.write("pif #\n")
            for obs in obs_names:
                f.write(f"l1 !{obs}!\n")

    # load obs.csv, set tag col as obs name
    obs_csv_path = os.path.join(pf.new_d, 'obs.csv')
    obs_csv = pd.read_csv(obs_csv_path)
    obs_csv["obsnme"] = obs_csv["tag"]
    obs_names = obs_csv["obsnme"].tolist()

    # create sim.csv inside pest directory folder if it doesn't exist ###
    sim_csv_filename = "sim.csv"  # JUST filename here
    sim_csv_path = os.path.join(pf.new_d, sim_csv_filename)
    if not os.path.exists(sim_csv_path):
        dummy_df = obs_csv[["tag", "conc"]].copy()
        dummy_df["conc"] = 1.0  # placeholder dummy values
        dummy_df.to_csv(sim_csv_path, index=False)

    # write the .ins file for sim.csv
    ins_file_path = sim_csv_path + ".ins"
    write_csv_ins_file(sim_csv_path, ins_file_path, obs_names)

    # add simulated observations using fnames
    pf.add_observations(
        sim_csv_filename,      # just filename here (important!)
        insfile=os.path.basename(ins_file_path),  # ins filename only
        prefix='sim',
        index_cols=["tag"],
        use_cols=["conc"],
        ofile_sep=','
    )

    # setup parameters from external csv
    mod_param_bnds = pd.read_csv(os.path.join('calibration', 'model_param_bounds.csv'))
    param_list = mod_param_bnds['param'].tolist()
    print(param_list)
    for param in param_list:
        print(param)
        pf.add_parameters(
            f"{param}.txt",  # filename only
            par_type="constant",
            par_name_base=f"{param}_cn",
            pargp=f"{param}_cn",
            lower_bound=mod_param_bnds[mod_param_bnds['param'] == param]['lower_bound_cn'].iloc[0],
            upper_bound=mod_param_bnds[mod_param_bnds['param'] == param]['upper_bound_cn'].iloc[0],
            ult_lbound=mod_param_bnds[mod_param_bnds['param'] == param]['ult_lbound'].iloc[0],
            ult_ubound=mod_param_bnds[mod_param_bnds['param'] == param]['ult_ubound'].iloc[0],
            transform="none"
        )
        pf.add_observations(f"{param}.txt", insfile=f"{param}.ins", prefix=f'{param}', includes_header=False)

    # define running phreeeqc cmd from comd line
    cmd = f"phreeqc.exe {inp_fnm} {out_fnm} {therm_fnm}"
    pf.mod_sys_cmds.append(cmd)

    # add helper fxn to build simulated pest obs for tracking
    pf.add_py_function('workflow.py', 'build_simulated_pest_obs()', is_pre_cmd=False)

    # build the pest object
    pst = pf.build_pst(filename='AR01.pst')
    print("PST observation names:", pst.observation_data.index[:10])
    pst.observation_data['obgnme'] = pst.observation_data.index

    # helper fxn to extract the tag from obs names
    def extract_tag(obs_name):
        parts = obs_name.split(':')
        if 'conc_tag' in parts:
            idx = parts.index('conc_tag')
            if idx + 1 < len(parts):
                return parts[idx + 1]
        return None

    obs_tags = pst.observation_data.index.to_series().apply(extract_tag)
    print("extracted tags sampling:\n", obs_tags.head(10))

    # # load weights from csv 
    # weights_df = pd.read_csv(os.path.join(pf.new_d, "obs_weights.csv"))
    # weights_df.set_index("obs", inplace=True)
    # print("weights DF index example:\n", weights_df.index.tolist())

    # ### mapping dictionary from extracted tags to weights keys
    # tag_to_weight_key = {
    #     "general.ph": "pH",
    #     "majorcation.na": "Na",
    #     "majorcation.ca": "Ca",
    #     "majorcation.mg": "Mg",
    #     "majorcation.k": "K",
    #     "majoranion.alkalinity": "Alkalinity",
    #     "majoranion.s(6)": "S(6)",
    #     "majoranion.cl": "Cl",
    #     "minoranion.f": "F",
    #     "trace.mn": "Mn",
    #     "trace.ba": "Ba",
    #     "trace.sr": "Sr",
    #     "trace.fe": "Fe",
    #     "trace.al": "Al",
    #     "trace.si": "Si",
    #     "trace.as": "As",
    #     "trace.u": "U",
    #     "trace.cr": "Cr",
    #     "trace.p": "P"
    # }

    # # map extracted tags to weights keys
    # mapped_keys = obs_tags.map(tag_to_weight_key).fillna('')
    # print("mapped keys preview:\n", mapped_keys.head(10))

    # # map weights using mapped keys
    # mapped_weights = mapped_keys.map(weights_df["weight"]).fillna(0.0)
    # print("mapped weights preview:\n", mapped_weights.head(10))

    # # assign weights to pst observation data
    # pst.observation_data["weight"] = mapped_weights.values

    # add post-processing to track stats
    pf.post_py_cmds.append('pyemu.helpers.calc_array_par_summary_stats()')

    # define pst options
    pst.control_data.noptmax = 0

    if check_pstpp_setup:
        obs_df = pst.observation_data.copy()
        obs_df['weight'] = 0.0 #init all weights
        pst.observation_data = obs_df
        pst.pestpp_options["ies_num_reals"] = reals
        pst.pestpp_options["save_binary"] = True
        pst.pestpp_options["panther_agent_freeze_on_fail"] = False
        pst.pestpp_options['ies_multimodal_alpha'] = 0.99
        pst.pestpp_options["ies_init_lam"] = -100 # very sensitive parameter option and important, default is around -10, -100 is rharmon recommendation, test -100 later
        #pst.pestpp_options['ies_drop_conflicts'] = True #keep the option disabled for now
        pst.write(os.path.join(pf.new_d, "AR01.pst"), version=2)
        if test and "phreeqc" in template_d:
            pyemu.os_utils.run(F"\\time -v -o timeit.txt pestpp-ies.exe AR01.pst", cwd=template_d)

        
        # # sometimes Ryan H will make the 1e-4 an assert statement, the observed should be ~= simulated with noptmax=0
        # pst.set_res(os.path.join(template_d,'sb.base.rei'))
        # print('noptmax = 0 pre-updating the obs vals and non-zero weights phi is: ', pst.phi)
        # if pst.phi > 1e-4:
        #     print('phi is greater than 1e-4, returning rei, investigate')

        print('Finished')
        return template_d

    else:
        obs_df = pst.observation_data.copy()
        obs_df['weight'] = 0.0 #init all weights
        pst.observation_data = obs_df
        obs_df = pst.observation_data.copy()
        # auto_weight_rebal code from jeremy white to aid to better group weight balancing with extreme range of estimated values
        if auto_weight_rebal:
            obs_df = setup_obs_val(t_d=template_d, obs_df=obs_df, set_weights_phi_facs_csv=0)
            pst.observation_data = obs_df
            pst.pestpp_options["ies_num_reals"] = reals
            pst.pestpp_options["save_binary"] = True
            pst.pestpp_options["panther_agent_freeze_on_fail"] = False
            pst.pestpp_options['ies_multimodal_alpha'] = 0.99
            pst.pestpp_options["ies_init_lam"] = -100 # very sensitive parameter option and important, default is around -10, -100 is rharmon recommendation, test -100 later
            pst.write(os.path.join(pf.new_d, "AR01.pst"), version=2)
            pst = pyemu.Pst(os.path.join(pf.new_d, "AR01.pst"))
            gnames = pst.nnz_obs_groups
            fac = 1./float(len(gnames))
            df = pd.DataFrame({"factor":fac},index=gnames)
            df.to_csv(os.path.join(pf.new_d, "phi_facs.csv"),index=True,header=False)
            pst.pestpp_options['ies_phi_factor_file'] = 'phi_facs.csv'
        else:
            obs_df = setup_obs_val(t_d=template_d, obs_df=obs_df, set_weights_phi_facs_csv=1)
            pst.observation_data = obs_df
            pst.pestpp_options["ies_num_reals"] = reals
            pst.pestpp_options["save_binary"] = True
            pst.pestpp_options["panther_agent_freeze_on_fail"] = False
            pst.pestpp_options['ies_multimodal_alpha'] = 0.99
            pst.pestpp_options["ies_init_lam"] = -100 # very sensitive parameter option and important, default is around -10, -100 is rharmon recommendation, test -100 later
            pst.pestpp_options['ies_phi_factor_file'] = 'phi_facs.csv'
        
        # re-write the pst control file
        pst.write(os.path.join(pf.new_d, "AR01.pst"), version=2)

        if test and "phreeqc" in template_d:
            pyemu.os_utils.run(F"\\time -v -o timeit.txt pestpp-ies.exe AR01.pst", cwd=template_d)

        print('Finished')

        # define prior jcb
        pe = pf.draw(reals, use_specsim=True)
        pe.to_binary(os.path.join(template_d, 'prior.jcb'))
        pst.pestpp_options['ies_par_en'] = 'prior.jcb'
        pst.write(os.path.join(pf.new_d, "AR01.pst"), version=2)
        
        # now draw obs noise ensemble
        draw_noise_reals(template_d, "AR01")

        return template_d

# fxn to run pstpp-ies locally
def run_pstppies_local(t_d, niters, reals, num_workers, pstFile = "AR01.pst"):

    print(f"running PEST locally with {num_workers} workers, {niters} iterations, and {reals} reals...")

    pst = pyemu.Pst(os.path.join(t_d, pstFile))
    pst.control_data.noptmax = niters
    pst.write(os.path.join(t_d, pstFile),version=2)

    date = datetime.datetime.today().strftime("%Y.%m.%d.%H%M%S")
    m_d = os.path.join(f"arap_ies.{date}.{niters}-{reals}")
    
    pyemu.os_utils.start_workers(t_d, # the folder which contains the "template" PEST dataset
                                'pestpp-ies.exe', #the PEST software version we want to run
                                pstFile, # the control file to use with PEST
                                num_workers=num_workers, #how many agents to deploy
                                worker_root='.', #where to deploy the agent directories; relative to where python is running
                                master_dir=m_d, #the manager directory
                                verbose = True
                                )
    return m_d


### define main function to run all functions above
if __name__ == '__main__':



    ### write phreeqc input file from source files
    #param_list = ['calcite', 'dolomite', 'gypsum', 'strontianite', 'witherite', 'co2g', 
    #              'cax2', 'mgx2', 'kx', 'nax', 'quartz', 'albite', 'kaolinite', 'muscovite',
    #              'uraninite', 'arsenopyrite', 'goethite', 'feoh3', 'siderite', 'mnoh3',
    #              'fluorapatite', 'halite', 'illite', 'chromite', 'mix-ratio-sol-upstream', 'mix-ratio-sol-downstream'
    #              ]
    #test_write_input_react_chain(temp_d, param_list)
    
    ### define input, output, thermodynamic database files
    #inp_fnm = "Arapahoe_01_template.txt"
    #out_fnm = inp_fnm + ".out"
    #therm_fnm = "llnl.dat"

    #### run phreeqc as testing
    #run_phreeqc(temp_d, inp_fnm, out_fnm, therm_fnm)


    print('{0}\n\npreparing flow-IES\n\n{1}'.format('*' * 17, '*' * 17))

    # determine number of workers based on available local resources
    num_workers_total = psutil.cpu_count(logical=False)
    num_workers = int(num_workers_total / 2) 
    setup_final_pstpp = 1
    if setup_final_pstpp:

        # define original directory & tempalte directory
        fpth = os.path.join(os.getcwd(), 'model_ws', 'Arapahoe-01')
        org_d = os.path.join(fpth, 'initial')
        temp_d = os.path.join(fpth,'template')
        temp_d_copy = os.path.join(fpth, 'template_temp')

        # if temp copy directory path exists, delete it

        # helper function to overwrite a windows readonly on the folder, if present
        def remove_readonly(func, path, excinfo):
            os.chmod(path,stat.S_IWRITE)
            func(path)
        # delete the temp folder copy
        if os.path.isdir(temp_d_copy):
            shutil.rmtree(temp_d_copy, onerror=remove_readonly)
            
        ### copy original to template directory
        copy_folder(org_d, temp_d)
    
        # set up pestpp-ies
        t_d = 'forward_model_template'
        reals = 60
        optiter = 5
        modnm = 'AR01'
        setup_pstppies(od=temp_d, template_d=t_d, reals=reals, test=1, check_pstpp_setup=0, auto_weight_rebal=1)
        
        # run pestpp-ies
        m_d = run_pstppies_local(t_d = t_d, niters = optiter, reals = reals, num_workers=num_workers)
        
        # plot results of pestpp-ies
        import plot_results
        plot_results.main(m_d=m_d, modnm=modnm, optiter=optiter, realz=reals)


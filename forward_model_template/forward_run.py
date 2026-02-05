import os
import multiprocessing as mp
import numpy as np
import pandas as pd
import pyemu

# function added thru PstFrom.add_py_function()
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



# function added thru PstFrom.add_py_function()
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



# function added thru PstFrom.add_py_function()
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



# function added thru PstFrom.add_py_function()
def build_obs():
    obs_vals = pd.read_csv(os.path.join('obs.csv'))
    return obs_vals



# function added thru PstFrom.add_py_function()
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


def main():

    try:
       os.remove(r'sim.csv')
    except Exception as e:
       print(r'error removing tmp file:sim.csv')
    try:
       os.remove(r'calcite.txt')
    except Exception as e:
       print(r'error removing tmp file:calcite.txt')
    try:
       os.remove(r'dolomite.txt')
    except Exception as e:
       print(r'error removing tmp file:dolomite.txt')
    try:
       os.remove(r'gypsum.txt')
    except Exception as e:
       print(r'error removing tmp file:gypsum.txt')
    try:
       os.remove(r'strontianite.txt')
    except Exception as e:
       print(r'error removing tmp file:strontianite.txt')
    try:
       os.remove(r'witherite.txt')
    except Exception as e:
       print(r'error removing tmp file:witherite.txt')
    try:
       os.remove(r'co2g.txt')
    except Exception as e:
       print(r'error removing tmp file:co2g.txt')
    try:
       os.remove(r'cax2.txt')
    except Exception as e:
       print(r'error removing tmp file:cax2.txt')
    try:
       os.remove(r'mgx2.txt')
    except Exception as e:
       print(r'error removing tmp file:mgx2.txt')
    try:
       os.remove(r'kx.txt')
    except Exception as e:
       print(r'error removing tmp file:kx.txt')
    try:
       os.remove(r'nax.txt')
    except Exception as e:
       print(r'error removing tmp file:nax.txt')
    try:
       os.remove(r'quartz.txt')
    except Exception as e:
       print(r'error removing tmp file:quartz.txt')
    try:
       os.remove(r'albite.txt')
    except Exception as e:
       print(r'error removing tmp file:albite.txt')
    try:
       os.remove(r'kaolinite.txt')
    except Exception as e:
       print(r'error removing tmp file:kaolinite.txt')
    try:
       os.remove(r'muscovite.txt')
    except Exception as e:
       print(r'error removing tmp file:muscovite.txt')
    try:
       os.remove(r'uraninite.txt')
    except Exception as e:
       print(r'error removing tmp file:uraninite.txt')
    try:
       os.remove(r'arsenopyrite.txt')
    except Exception as e:
       print(r'error removing tmp file:arsenopyrite.txt')
    try:
       os.remove(r'goethite.txt')
    except Exception as e:
       print(r'error removing tmp file:goethite.txt')
    try:
       os.remove(r'feoh3.txt')
    except Exception as e:
       print(r'error removing tmp file:feoh3.txt')
    try:
       os.remove(r'siderite.txt')
    except Exception as e:
       print(r'error removing tmp file:siderite.txt')
    try:
       os.remove(r'mnoh3.txt')
    except Exception as e:
       print(r'error removing tmp file:mnoh3.txt')
    try:
       os.remove(r'fluorapatite.txt')
    except Exception as e:
       print(r'error removing tmp file:fluorapatite.txt')
    try:
       os.remove(r'halite.txt')
    except Exception as e:
       print(r'error removing tmp file:halite.txt')
    try:
       os.remove(r'illite.txt')
    except Exception as e:
       print(r'error removing tmp file:illite.txt')
    try:
       os.remove(r'chromite.txt')
    except Exception as e:
       print(r'error removing tmp file:chromite.txt')
    try:
       os.remove(r'mix-ratio-sol-upstream.txt')
    except Exception as e:
       print(r'error removing tmp file:mix-ratio-sol-upstream.txt')
    try:
       os.remove(r'mix-ratio-sol-downstream.txt')
    except Exception as e:
       print(r'error removing tmp file:mix-ratio-sol-downstream.txt')
    pyemu.helpers.apply_list_and_array_pars(arr_par_file='mult2model_info.csv',chunk_len=50)
    build_simulated_pest_obs()
    delete_old_model_output()
    write_input_react_chain()
    build_obs()
    pyemu.os_utils.run(r'phreeqc.exe Arapahoe_01_template.txt Arapahoe_01_template.txt.out llnl.dat')

    build_simulated_pest_obs()

if __name__ == '__main__':
    mp.freeze_support()
    main()


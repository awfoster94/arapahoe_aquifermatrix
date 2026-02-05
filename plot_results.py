import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.animation import FuncAnimation
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import pyemu
import numpy as np
import pandas as pd
import subprocess
import shutil
import time
#import wqchartpy
import ions
import triangle_piper_mod
import stiff_mod
import schoeller_mod

# fxn to extract ies obs to a dictionary
def get_ies_obs_dict(m_d="master_ies", modnm='AR01'):
    pst = pyemu.Pst(os.path.join(m_d,f"{modnm}.pst"))
    if pst.control_data.noptmax == -1:
        itrs = [0]
    else:
        itrs = range(pst.control_data.noptmax+1)

    obs_df_dict = {}
    for i in itrs:
        print(f'loading itr {i}')
        # if usebin:
        obs_df = pyemu.Matrix.from_binary(os.path.join(m_d, f'{modnm}.{i}.obs.jcb')).to_dataframe()
        obs_df['real_name'] = obs_df.index
        print(f'loaded itr {i} - nreals = {len(obs_df)}')
        obs_df_dict[i] = obs_df

    return obs_df_dict

# fxn plot boxplots of simulated ensemble data, overlays observed and base values, and exports csvs for each optimization iteration
def plot_and_export_obs_boxplots_log(m_d, obsdict, log_input=True, modnm='AR01'):

    # define helper fxn for parameter transformations out of log space
    def convert_value(val, obsnme):

        name_lower = obsnme.lower()

        # pH transformations
        # if pH is logged: [val = -log10(pH)], then convert back
        # else val is already pH, then no change
        if "ph" in name_lower:  
            return 10 ** (-val) if log_input else val

        # transform parameters out of log space
        return 10 ** (-val) if log_input else val

    # setup output directory
    fdir = os.path.join(m_d, 'prelim_figs')
    os.makedirs(fdir, exist_ok=True)

    # define/read in pst objects
    pst = pyemu.Pst(os.path.join(m_d, f"{modnm}.pst"))
    nzobs_all = pst.observation_data.loc[pst.nnz_obs_names, :].copy()
    all_obsnames = nzobs_all.obsnme.unique()
    iters = sorted(obsdict.keys())

    # export to csvs for easy access, review
    for itr in iters:
        csv_rows = []

        for obsnme in all_obsnames:
            if obsnme not in obsdict[itr].columns:
                continue

            # clean up obs name
            clean_obsnme = obsnme
            for tag in ['trace','general','majorion','minorion','majorcation','majoranion','minoranion']:
                clean_obsnme = clean_obsnme.replace(tag, '')
            tag = clean_obsnme.split(":")[-1]
            safe_tag = "".join(c for c in tag if c.isalnum() or c in (' ', '_', '-')).rstrip()

            # observed values
            obsval_raw = nzobs_all.loc[nzobs_all.obsnme == obsnme, "obsval"].iloc[0]
            obsval = convert_value(obsval_raw, obsnme)

            # simulated values
            simvals = obsdict[itr][obsnme].dropna()

            for real, simval in simvals.items():
                csv_rows.append({"realization": real, "obsnme": clean_obsnme,"tag": safe_tag, "observed": obsval, "simulated": convert_value(simval, obsnme), "is_base": real == "base", "iteration": itr})

        if not csv_rows:
            continue

        df_csv = pd.DataFrame(csv_rows)

        # pivot to 1 row per realization
        df_pivot = df_csv.pivot(index="realization", columns="tag", values="simulated")
        df_pivot["is_base"] = df_csv.drop_duplicates("realization").set_index("realization")["is_base"]
        df_pivot["iteration"] = itr

        # add observed row
        obs_row = dict(df_pivot.iloc[0])
        obs_row.update(df_csv.groupby("tag").first()["observed"].to_dict())
        obs_row["realization"] = "observed"
        obs_row["is_base"] = False
        obs_row["iteration"] = itr
        df_pivot.loc["observed"] = obs_row

        # write csv
        csv_out = os.path.join(fdir, f"ensemble_obs_iter{itr}.csv")
        df_pivot.to_csv(csv_out)
        print(f"saved csv for iter {itr}: {csv_out}")

    # plotting
    for obsnme in all_obsnames:

        # clean name again
        clean_obsnme = obsnme
        for tag in ['trace','general','majorion','minorion','majorcation','majoranion','minoranion']:
            clean_obsnme = clean_obsnme.replace(tag, '')
        tag = clean_obsnme.split(":")[-1]
        safe_tag = "".join(c for c in tag if c.isalnum() or c in (' ', '_', '-')).rstrip()

        pdf_path = os.path.join(fdir, f"{safe_tag}.pdf")

        with PdfPages(pdf_path) as pdf:

            plot_data = []

            for itr in iters:
                if obsnme not in obsdict[itr].columns:
                    continue
                obsval_raw = nzobs_all.loc[nzobs_all.obsnme == obsnme, "obsval"].iloc[0]
                obsval = convert_value(obsval_raw, obsnme)
                simvals = obsdict[itr][obsnme].dropna()
                for real, val in simvals.items():
                    plot_data.append({"iteration": f"Iter {itr}", "realization": real, "sim_value": convert_value(val, obsnme),"obs_value": obsval, "is_base": real == "base"})

            if not plot_data:
                continue

            df_plot = pd.DataFrame(plot_data)

            # plot start
            plt.figure(figsize=(10, 6))

            sns.boxplot(x="iteration", y="sim_value", data=df_plot, color="lightblue", fliersize=0)
            sns.stripplot(x="iteration", y="sim_value", data=df_plot[~df_plot["is_base"]], color="blue", alpha=0.5, size=4, jitter=True,label="Realizations")

            # base realization
            df_base = df_plot[df_plot["is_base"]]
            plt.scatter(df_base["iteration"], df_base["sim_value"], color="black", s=50, zorder=10, label="Base")

            # observed data
            for itr in df_plot["iteration"].unique():
                obsval = df_plot[df_plot["iteration"] == itr]["obs_value"].iloc[0]
                plt.scatter(itr, obsval, edgecolors="yellow", facecolors="none", s=70, zorder=11, label="Observed")

            # labels
            ylabel = "Simulated Value (real units)"
            if "ph" in obsnme.lower():
                ylabel = "pH"

            plt.title(f"Simulated Ensemble vs Observed: {safe_tag}")
            plt.ylabel(ylabel)
            plt.xlabel("Iteration")
            plt.grid(True, alpha=0.3)
            plt.yscale("log")

            # duplicate legend
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys())

            pdf.savefig(dpi=150)
            plt.close()

        print(f"saved to pdf: {pdf_path}")

    return fdir

# fxn plots modified trilinear, schoeller, stiff piper diagrams for each csv file in the specified directory
# each realization is color-coded: light blue for realizations, black for base, and gold for observed values
def plot_wq_diagrams(fdir):
    # iterate over each csv file in the directory
    for csv_file in os.listdir(fdir):
        if csv_file.endswith(".csv"):
            csv_path = os.path.join(fdir, csv_file)
            df = pd.read_csv(csv_path)

            # check if the necessary columns exist
            required_columns = [
                'Sample', 'Label', 'Color', 'Marker', 'Size', 'Alpha',
                'Ca2+', 'Mg2+', 'Na+', 'K+', 'HCO3-', 'CO3--', 'Cl-', 'SO4--'
            ]
            if not all(col in df.columns for col in required_columns):
                print(f"Missing required columns, now creating required columns and converting from mol/kgw to mg/L....")
                df['Sample'] = df['realization']
                df['Label'] = df['realization']
                df.loc[~df['realization'].str.contains('base|observed', na=False), 'Label'] = 'realizations'
                
                df['Marker'] = np.where(df['Label'] == 'base', 's',
                    np.where(df['Label'] == 'observed', 'o', 'o'))

                df['Color'] = np.where(df['Label'] == 'base', 'black',
                                np.where(df['Label'] == 'observed', 'yellow', 'lightblue'))

                df['Size'] = np.where(df['Label'] == 'base', 80,
                            np.where(df['Label'] == 'observed', 70, 40))

                df['Alpha'] = np.where(df['Label'] == 'base', 1.0,
                                np.where(df['Label'] == 'observed', 1.0, 0.6))
                
                df['Ca'] = df['ca'] * 1000 * 40.08 # convert from mol perkgw to mg/L
                df['Mg'] = df['mg'] * 1000 * 24.305# convert from mol perkgw to mg/L
                df['Na'] = df['na'] * 1000 * 22.9898 # convert from mol perkgw to mg/L
                df['K'] = df['k'] * 1000 * 39.0983 # convert from mol perkgw to mg/L
                df['HCO3'] = df['alkalinity'] * 1000 * 2 * 50.044 # convert from mol perkgw to mg/L
                df['CO3'] = 0 # convert from mol perkgw to mg/L
                df['Cl'] = df['cl'] * 1000 * 35.45 # convert from mol perkgw to mg/L
                df['SO4'] = df['s6'] * 1000 * 96.986 # convert from mol perkgw to mg/L

            print(df)

            ## generate the Piper, Schoeller, and Stiff Diagrams
            fig_name = os.path.splitext(csv_file)[0]
            fig_fpth_piper = os.path.join(fdir, fig_name+'piper')
            triangle_piper_mod.plot(df, unit='mg/L', figname=fig_fpth_piper, figformat='jpg')
            print(f"saved piper diagram...")

            fig_fpth_schoeller = os.path.join(fdir, fig_name+'schoeller')
            schoeller_mod.plot(df, unit='mg/L', figname=fig_fpth_schoeller, figformat='jpg')
            fig = plt.gcf()
            ax = plt.gca()
            ax.set_ylim(1E-9, 10) # set ylim
            fig.savefig(fig_fpth_schoeller + ".jpg", dpi=300)
            print(f"saved schoeller diagram...")

            plot_stiff = 0
            if plot_stiff:
                if '5' in csv_file and csv_file.endswith('.csv'): 
                    fig_fpth_stiff = os.path.join(fdir, fig_name+'stiff')
                    stiff_mod.plot(df, unit='mg/L', figname=fig_fpth_stiff, figformat='jpg')
                    print(f"saved stiff diagrams...")

# fxn to plot parameter histograms quickly
def plot_simple_par_histo(m_d, modnm='AR01'):
    pst = pyemu.Pst(os.path.join(m_d,f'{modnm}.pst'))
    par = pst.parameter_data
    adjpar = par.loc[par.partrans.apply(lambda x: x in ['none','log']),:]
    groups = adjpar.pargp.unique()
    groups.sort()
    #pr = pst.ies.paren0
    par0_file = os.path.join(m_d, "prior.jcb")
    pr = pyemu.ParameterEnsemble.from_binary(filename=par0_file, pst=pst)

    fout = os.path.join(m_d,'prelim_figs','param_distribs')
    os.makedirs(fout, exist_ok=True)

    # detect available .parjcb iterations (excluding 0, prior)
    itr_files = [
        f for f in os.listdir(m_d)
        if f.startswith(f"{modnm}.") and f.endswith("par.jcb")
    ]
    itrs = sorted([int(f.split('.')[1]) for f in itr_files])  # extract iteration number

    print("detected parameter ensemble iterations:", itrs)

    for itr in itrs:
        #pt = pst.ies.__getattr__('paren{0}'.format(itr))
        itr_file = os.path.join(m_d, f"{modnm}.{itr}.par.jcb")
        if not os.path.exists(itr_file):
            print(f"Skipping missing iteration {itr}")
            continue

        with PdfPages(os.path.join(fout,f'simple_par_histo_iter{itr}.pdf')) as pdf:

            #for pname in pst.adj_par_names:
            for group in groups:
                gpar = adjpar.loc[adjpar.pargp==group,:].copy()

                lb = gpar.parlbnd.min()
                ub = gpar.parubnd.max()
                if gpar.partrans.iloc[0] != 'none':
                    lb = np.log10(lb)
                    ub = np.log10(ub)

                # load ensemble for iteration
                pt = pyemu.ParameterEnsemble.from_binary(filename=itr_file, pst=pst)
                
                fig,ax = plt.subplots(1,1,figsize=(6,6))
                if gpar.partrans.iloc[0] == 'none':
                    ax.hist(pr.loc[:,gpar.parnme].values.flatten(),bins=20,facecolor='0.5',edgecolor='none',alpha=0.5,density=True)
                    ax.hist(pt.loc[:,gpar.parnme].values.flatten(),bins=20,facecolor='b',edgecolor='none',alpha=0.5,density=True)
                    ax.set_xlabel('')
                else:
                    ax.hist(np.log10(pr.loc[:,gpar.parnme].values.flatten()),bins=20,facecolor='0.5',edgecolor='none',alpha=0.5,density=True)
                    ax.hist(np.log10(pt.loc[:,gpar.parnme].values.flatten()),bins=20,facecolor='b',edgecolor='none',alpha=0.5,density=True)
                    ax.set_xlabel('$log_{10}$')
                ylim = ax.get_ylim()
                ax.plot([lb,lb],ylim,'k--',lw=3)
                ax.plot([ub,ub],ylim,'k--',lw=3)
                ax.set_ylim(ylim)
                ax.set_title('iteration:{0} pname:{1}, npar:{2}'.format(itr,group,gpar.shape[0]),loc='left')
                plt.tight_layout()
                pdf.savefig()
                plt.close(fig)
                print(itr,group)

# fxn to plot the relative mole fractions for each realization at an optimization iteration of interest
def plot_aquifermatrix_comps(m_d, modnm='AR01', optiter=5, realz=48):
    
    # first let's define helper fxns
    # fxn to extract mole fractions stored in pst object
    def mineral_mole_fractions(pst, ensemble, cn_params):
        
        # store par data
        par = pst.parameter_data

        # extract param values for moles of aquifer matrix constituents
        df = ensemble.loc[:, cn_params].copy()

        # ensure that the mole parameters
        for p in cn_params:
            if par.loc[p, "partrans"] == "log":
                df[p] = 10.0 ** df[p]

        frac = df.div(df.sum(axis=1), axis=0)
        frac.columns = [c.replace("cn_", "") for c in frac.columns]

        return frac

    # fxn to plot stacked bar plot of each realizations mineral groupings normalized by mole fraction
    def plot_stacked_mineral_realizations_grouped(frac, color_map, max_realizations=realz):

        if frac.shape[0] > max_realizations:
            frac = frac.sample(max_realizations, random_state=0)
        
        x = np.arange(frac.shape[0])
        bottom = np.zeros(frac.shape[0])
        
        fig, ax = plt.subplots(figsize=(16,6))
        
        # stacked bars
        for mineral in frac.columns:
            ax.bar(x, frac[mineral].values, bottom=bottom, width=1.0, color=color_map[mineral])
            bottom += frac[mineral].values
        ax.set_xlim(-0.5, frac.shape[0]-0.5)
        ax.set_ylim(0,1)
        ax.set_ylabel('normalized mole fraction (-)', fontsize=16)
        ax.set_xlabel('model realizations \nsimulated with pestpp-ies', fontsize=16)
        ax.set_title('Aquifer matrix compositions', fontsize=18)
        ax.set_xticks(x)
        ax.set_xticklabels(x, rotation=0)
        
        # define mineral groups
        mineral_groups = [
            ('Gases', gases),
            ('Surface Exchange', surface_exchange),
            ('Oxides/Sulfides', oxides_sulfides),
            ('Silicates', silicates),
            ('Evaporites', evaporite),
            ('Carbonates', carbonate)
            ]
        
        # adjustable legend placement
        start_y = -0.10   # vertical position of first legend row (closer to plot)
        y_step = 0.08     # vertical spacing between legend rows
        
        # create one legend per group, stacked below the plot
        for i, (grp_name, minerals_list) in enumerate(mineral_groups):
            handles = [Rectangle((0,0),1,1,color=color_map[m]) for m in minerals_list]
            labels = minerals_list
            fig.legend(handles, labels, ncol=len(labels),
                    title=grp_name, loc='lower center',
                    bbox_to_anchor=(0.5, start_y - i*y_step), frameon=False,
                    fontsize=14)
        
        # leave more space at the bottom for legends
        fig.subplots_adjust(bottom=0.15)
        plt.tight_layout()
        outpth = os.path.join(m_d, 'prelim_figs', 'aquifermatrix')
        if not os.path.exists(outpth):
            os.makedirs(outpth)
        plt.savefig(os.path.join(outpth, 'aquifermatrix_relativemolefraction_grouped.png'), dpi=400, bbox_inches='tight')
        plt.show()
    
    # define list of minerals to plot
    minerals = [
        "calcite","dolomite","gypsum","strontianite","witherite",
        "co2g","cax2","mgx2","kx","nax","quartz","albite",
        "kaolinite","muscovite","uraninite","arsenopyrite",
        "goethite","feoh3","siderite","mnoh3","fluorapatite",
        "halite","illite","chromite"
    ]

    # define mineral composition groups for plotting
    carbonate = ['calcite', 'dolomite', 'siderite', 'strontianite', 'witherite']
    evaporite = ['gypsum', 'halite', 'fluorapatite']
    silicates = ['quartz', 'albite', 'kaolinite', 'muscovite', 'illite']
    oxides_sulfides = ['goethite', 'arsenopyrite', 'feoh3', 'uraninite', 'chromite', 'mnoh3']  # added mnoh3
    surface_exchange = ['cax2','mgx2','kx','nax']
    gases = ['co2g']

    # create custom color maps for each group
    group_colors = {
        'carbonate': plt.cm.Blues(np.linspace(0.3, 0.8, len(carbonate))),
        'evaporite': plt.cm.Reds(np.linspace(0.3, 0.8, len(evaporite))),
        'silicates': plt.cm.YlGn(np.linspace(0.3, 0.8, len(silicates))),
        'oxides_sulfides': plt.cm.Purples(np.linspace(0.3, 0.8, len(oxides_sulfides))),
        'surface_exchange': plt.cm.copper(np.linspace(0.3, 0.8, len(surface_exchange))),
        'gases': plt.cm.Greys(np.linspace(0.3, 0.8, len(gases))),
    }

    # define group order for plotting the legend
    group_order = [
        ('carbonate', carbonate),
        ('evaporite', evaporite),
        ('silicates', silicates),
        ('oxides_sulfides', oxides_sulfides),
        ('surface_exchange', surface_exchange),
        ('gas', gases)
    ]

    # load pst object
    pst = pyemu.Pst(os.path.join(m_d, f'{modnm}.pst'))

    pe = pyemu.ParameterEnsemble.from_binary(pst, os.path.join(f'{m_d}', f'{modnm}.{optiter}.par.jcb'))

    # select ensemble columns that correspond to your minerals
    cn_params = [c for c in pe.columns if any(m in c for m in minerals)]

    # use helper fxn to extract mineral mole fractions
    frac = mineral_mole_fractions(pst, pe, cn_params)

    # clean-up mineral names from pst observations names
    clean_names = {c: c.split(":")[1].split("_cn")[0] for c in frac.columns}
    clean_names_final = {c: c.split(":")[1].split("_inst")[0] for c in frac.columns}
    frac_clean = frac.rename(columns=clean_names_final)

    # map each mineral to each respective color
    color_map = {}
    for group_name, minerals_list in zip(group_colors.keys(),
        [carbonate, evaporite, silicates, oxides_sulfides, surface_exchange, gases]):
        colors = group_colors[group_name]
        for m, c in zip(minerals_list, colors):
            color_map[m] = c

    # flatten the ordered list of minerals
    ordered_minerals = [m for grp, minerals_list in group_order for m in minerals_list]

    # reorder frac columns
    frac_ordered = frac_clean[ordered_minerals]

    # finally plot the mole fraction
    plot_stacked_mineral_realizations_grouped(frac=frac_ordered, color_map=color_map, max_realizations=realz)

# function to run phreeqc 
def run_phreeqc(input_file, output_file, therm_dat_file):
    # Get the current working directory
    current_directory = os.getcwd()

    # Specify the path to phreeqc.exe
    phreeqc_path = os.path.join(current_directory, "phreeqc.exe")  # Replace with the actual relative path to phreeqc.exe

    # Build the command to execute
    command = [
        phreeqc_path,
        os.path.join(current_directory, input_file),
        os.path.join(current_directory, output_file),
        os.path.join(current_directory, therm_dat_file)
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

# fxn to run phreeqc mixing models based on aquifer matrix compositions
def write_phreeqc_mixing(m_d, modnm='AR01', optiter=5, realz=60):
    
    # fxn to extract moles stored in pst object
    def mineral_moles(pst, ensemble, cn_params):
        
        # store par data
        par = pst.parameter_data

        # extract param values for moles of aquifer matrix constituents
        df = ensemble.loc[:, cn_params].copy()

        # ensure that the mole parameters
        for p in cn_params:
            if par.loc[p, "partrans"] == "log":
                df[p] = 10.0 ** df[p]
        frac = df.div(df.sum(axis=1), axis=0)
        df.columns = [c.replace("cn_", "") for c in df.columns]

        return df
    
    # define list of minerals to plot
    minerals = [
        "calcite","dolomite","gypsum","strontianite","witherite",
        "co2g","cax2","mgx2","kx","nax","quartz","albite",
        "kaolinite","muscovite","uraninite","arsenopyrite",
        "goethite","feoh3","siderite","mnoh3","fluorapatite",
        "halite","illite","chromite"]#, "mix-ratio-sol-upstream", "mix-ratio-sol-downstream"]

    # load pst object
    pst = pyemu.Pst(os.path.join(m_d, f'{modnm}.pst'))

    pe = pyemu.ParameterEnsemble.from_binary(pst, os.path.join(f'{m_d}', f'{modnm}.{optiter}.par.jcb'))

    # select ensemble columns that correspond to your minerals
    cn_params = [c for c in pe.columns if any(m in c for m in minerals)]

    # use helper fxn to extract mineral mole fractions
    moles = mineral_moles(pst, pe, cn_params)

    for real in range(realz):
        
        aquifer_matrix_moles = moles.iloc[real].to_dict()

        input_file_text = f"""

        KNOBS
            -iterations 1000

            EQUILIBRIUM_PHASES 1
                Calcite	0.0 0.0
                Dolomite 0.0 0.0
                Gypsum 0.0 0.0	
                Strontianite 0.0 0.0 #dissolve_only
                Witherite 0.0 0.0 #dissolve_only
                CO2(g) 0.0 0.0
                CaX2 0.0 0.0 
                MgX2 0.0 0.0 
                KX 0.0 0.0 
                NaX 0.0 0.0 
                Quartz 0.0 0.0 #dissolve_only
                Albite 0.0 0.0	#dissolve_only
                Kaolinite 0.0 0.0 #dissolve_only
                Muscovite 0.0 0.0
                Uraninite 0.0 0.0
                Arsenopyrite 0.0 0.0 #dissolve_only
                Goethite 0.0 0.0
                Fe(OH)3 0.0 0.0 
                Siderite 0.0 0.0
                Mn(OH)3 0.0 0.0			
                Fluorapatite 0.0 0.0 
                Halite 0.0 0.0
                Illite 0.0 0.0
                Chromite 0.0 0.0	#dissolve_only
            
            END
            
            EQUILIBRIUM_PHASES 2
            CO2(g)     -3.5     1.0 
            O2(g)      -0.76     1.0 
            
            END
            
            REACTION_TEMPERATURE 1
                14.6
                
            REACTION_PRESSURE 1
                3
            
            REACTION_TEMPERATURE 2
                13.327272727272728
                
            REACTION_PRESSURE 2
                0.82
            
            Solution 1	GS-392400104150601
            units		mg/L
            redox       pe
            temp		14.6
            pH		    7.4
            O(0)        0.1
            Na		    273
            Ca		    66.7
            Mg		    6.3
            K           4.73
            Alkalinity	163	    as CaCO3
            S(6)		 599
            Cl		    4.15
            F		    0.44
            Mn		    0.113
            Ba          0.0212
            Sr          1.1
            Fe          0.418
            As          0.0001
            U           2.1000000000000002e-05
            Li          0.0213
            Cr          0.0004
            N(5)        0.03	   as N
            N(3)        0.004	   as N
            N(-3)       1.15	   as N
            P           0.008    as P
            
            water       1
            SAVE solution 1
            END
            
            Solution	2 GS-385646104504601
            units		mg/L
            redox       pe
            temp		13.327272727272728
            pH		    7.568939393939393
            O(0)        8 
            Na		    34.07741935483871
            Ca		    44.91973684210526
            Mg		    9.723684210526315
            K           3.7077419354838708
            Alkalinity	83.64516129032258	as CaCO3
            S(6)		72.99354838709677
            Cl		    39.20967741935484
            F		    0.6829032258064517
            Mn		    0.042124387755102036
            N(5)        0.8441095890410958	as N
            N(-3)       0.3222222222222222	as N
            
            water       1
            USE EQUILIBRIUM_PHASES 2
            USE REACTION_TEMPERATURE 2
            USE REACTION_PRESSURE 2
            SAVE Solution 2
            END
            
            TITLE Mixing Model, groundwater : recharge source water 1
            MIX 1
            1 0.5
            2 0.5
            SAVE solution 3   
            END

            REACTION 0 Convergence
            H20 0
            END

            REACTION 1
                Calcite        {aquifer_matrix_moles['pname:calcite_cn_inst:0_ptype:cn_pstyle:m']}
                Dolomite       {aquifer_matrix_moles['pname:dolomite_cn_inst:0_ptype:cn_pstyle:m']}
                Gypsum         {aquifer_matrix_moles['pname:gypsum_cn_inst:0_ptype:cn_pstyle:m']}
                Strontianite   {aquifer_matrix_moles['pname:strontianite_cn_inst:0_ptype:cn_pstyle:m']}
                Witherite      {aquifer_matrix_moles['pname:witherite_cn_inst:0_ptype:cn_pstyle:m']}
                CO2(g)         {aquifer_matrix_moles['pname:co2g_cn_inst:0_ptype:cn_pstyle:m']}
                CaX2           {aquifer_matrix_moles['pname:cax2_cn_inst:0_ptype:cn_pstyle:m']}
                MgX2           {aquifer_matrix_moles['pname:mgx2_cn_inst:0_ptype:cn_pstyle:m']}
                KX             {aquifer_matrix_moles['pname:kx_cn_inst:0_ptype:cn_pstyle:m']}
                NaX            {aquifer_matrix_moles['pname:nax_cn_inst:0_ptype:cn_pstyle:m']}
                Quartz         {aquifer_matrix_moles['pname:quartz_cn_inst:0_ptype:cn_pstyle:m']}
                Albite         {aquifer_matrix_moles['pname:albite_cn_inst:0_ptype:cn_pstyle:m']}
                Kaolinite      {aquifer_matrix_moles['pname:kaolinite_cn_inst:0_ptype:cn_pstyle:m']}
                Muscovite      {aquifer_matrix_moles['pname:muscovite_cn_inst:0_ptype:cn_pstyle:m']}
                Uraninite      {aquifer_matrix_moles['pname:uraninite_cn_inst:0_ptype:cn_pstyle:m']}
                Arsenopyrite   {aquifer_matrix_moles['pname:arsenopyrite_cn_inst:0_ptype:cn_pstyle:m']}
                Goethite       {aquifer_matrix_moles['pname:goethite_cn_inst:0_ptype:cn_pstyle:m']}
                Fe(OH)3        {aquifer_matrix_moles['pname:feoh3_cn_inst:0_ptype:cn_pstyle:m']}
                Siderite       {aquifer_matrix_moles['pname:siderite_cn_inst:0_ptype:cn_pstyle:m']}
                Mn(OH)3        {aquifer_matrix_moles['pname:mnoh3_cn_inst:0_ptype:cn_pstyle:m']}
                Fluorapatite   {aquifer_matrix_moles['pname:fluorapatite_cn_inst:0_ptype:cn_pstyle:m']}
                Halite         {aquifer_matrix_moles['pname:halite_cn_inst:0_ptype:cn_pstyle:m']}
                Illite         {aquifer_matrix_moles['pname:illite_cn_inst:0_ptype:cn_pstyle:m']}
                Chromite       {aquifer_matrix_moles['pname:chromite_cn_inst:0_ptype:cn_pstyle:m']}

            END

            
            TITLE React Mixture with Aquifer Matrix & Equilibrate with Subsurface Conditions
            USE solution 3
            USE REACTION 1
            USE EQUILIBRIUM_PHASES 1
            USE REACTION_TEMPERATURE 1
            USE REACTION_PRESSURE 1
            SAVE solution 4

            SELECTED_OUTPUT
                file arapahoe-realization_{real}-mixing.txt
                selected_out      TRUE
                pH     TRUE
                charge_balance     TRUE
                percent_error     TRUE
                totals     O(0)
                -saturation_indices		Aragonite Brucite Calcite Dolomite Dolomite-dis Fe(OH)3 Goethite Gypsum Gibbsite Magnesite Mn(OH)3 Rhodochrosite Siderite Sylvite Strontianite Witherite
                -equilibrium_phases      Aragonite Brucite Calcite Dolomite Dolomite-dis Fe(OH)3 Goethite Gypsum Gibbsite Magnesite Mn(OH)3 Rhodochrosite Siderite Sylvite Strontianite Witherite
                
            END
    """
        
        ### save current directory to switch back to
        current_working_dir = os.getcwd()
        
        ### make charge balance folder to write models to 
        folder_name = "aquifer_matrix_mixing_models"
        os.makedirs(folder_name, exist_ok=True)
        
        ### copy phreeqc executable and llnl.dat files to the new folder
        shutil.copy("bin/phreeqc.exe", folder_name)
        shutil.copy("bin/llnl.dat", folder_name)
        
        ### change directory
        os.chdir(folder_name)
        
        ### output phreeqc input file
        with open("arapahoe-realization_"+str(real)+".txt", "w") as file:
            file.write(input_file_text)
        
        ### define input file, output, file, thermodynamic database file
        input_file =  "arapahoe-realization_"+str(real)+".txt"    
        output_file = "arapahoe-realization_"+str(real)+".txt.out" 
        therm_dat_file = "llnl.dat"      

        ### run phreeqc charge balance models
        run_phreeqc(input_file, output_file, therm_dat_file)
        
        ### change back to original directory
        os.chdir(current_working_dir)

    # now let's load the results from the selected output files

    # helper fxn to load data
    def collect_mineral_precipmass(folder, realz):
        
        all_reals = []
        for real in range(realz):
            fname = os.path.join(folder, f"arapahoe-realization_{real}-mixing.txt")
            if not os.path.exists(fname):
                print(f"Missing: {fname}")
                continue
            df = pd.read_csv(fname, delim_whitespace=True)
            # remove leading/trailing whitespace from column names
            df.columns = df.columns.str.strip()
            # select d_ columns robustly
            d_cols = [c for c in df.columns if c.startswith("d_")]
            if not d_cols:
                raise ValueError(f"No d_ columns found in {fname}")
            # take final row
            d_sum = df.loc[df.index[-1], d_cols]
            row = d_sum.to_dict()
            row["realization"] = real
            all_reals.append(row)
        return pd.DataFrame(all_reals).set_index("realization")
    
    precip_mass = collect_mineral_precipmass(folder=os.path.join("aquifer_matrix_mixing_models"), realz=realz)
    #print(precip_mass)
    
    # sum precipitated mass
    precip_mass_row_sum = precip_mass.sum(axis=1)
    precip_mass_sum = precip_mass_row_sum.to_frame(name="precip_mass_sum")
    
    # normalize precipitated mass
    precip_mass_sum["precip_mass_sum_norm"] = ((precip_mass_sum["precip_mass_sum"] - precip_mass_sum["precip_mass_sum"].min()) / (precip_mass_sum["precip_mass_sum"].max() - precip_mass_sum["precip_mass_sum"].min()))
    precip_mass_sum["precip_mass_sum_norm_score"] = 1 - precip_mass_sum["precip_mass_sum_norm"]
    
    # now create an animated plot! 
    x = precip_mass_sum.index.values
    y = precip_mass_sum["precip_mass_sum_norm_score"].values

    cvals = 2 * (y - 0.5)

    fig, ax = plt.subplots(figsize=(8, 4))

    sc = ax.scatter([], [], c=[], cmap="RdBu", vmin=-1, vmax=1, s=70)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Realization", fontsize=12)
    ax.set_ylabel("Normalized Scoring Value(-)", fontsize=14)
    ax.set_title("Mixing Compatibility Scoring")

    # helper fxn for animation
    def update(frame):
        sc.set_offsets(np.column_stack((x[:frame+1], y[:frame+1])))
        sc.set_array(2 * (y[:frame+1] - 0.5))
        return sc,
    # create animation
    ani = FuncAnimation(fig, update, frames=len(x), interval=80, blit=True)
    # save animation
    ani.save("mixing_compatibility_scoring_withaquifermatrix.gif", writer="pillow", dpi=400)


    plt.show()

# main fxn to call all plotting functions
def main(m_d, modnm, optiter, realz):
    # run all post-processing
    obsdict = get_ies_obs_dict(m_d, modnm)
    fdir = plot_and_export_obs_boxplots_log(m_d, obsdict, log_input=True, modnm=modnm)
    plot_wq_diagrams(fdir)
    plot_simple_par_histo(m_d, modnm)
    plot_aquifermatrix_comps(m_d, modnm, optiter, realz)
    write_phreeqc_mixing(m_d, modnm, optiter, realz)

# fxn to plot log-scale boxplots of simulated ensemble data for each observation (obsnme),
# overlays observed and base values, and exports csvs of all obs per optimization iteration.
# def plot_and_export_obs_boxplots_log(m_d, obsdict):
#     # define output directory
#     fdir = os.path.join(m_d, 'prelim_figs')
#     os.makedirs(fdir, exist_ok=True)

#     # define/read in pest object data
#     pst = pyemu.Pst(os.path.join(m_d, "AR01.pst"))
#     nzobs_all = pst.observation_data.loc[pst.nnz_obs_names, :].copy()
#     all_obsnames = nzobs_all.obsnme.unique()
#     iters = sorted(obsdict.keys())

#     # prepare one combined DataFrame per iteration (for CSV export)
#     for itr in iters:
#         csv_rows = []

#         for obsnme in all_obsnames:
#             if obsnme not in obsdict[itr].columns:
#                 continue

#             # Remove specific tags from obsnme
#             clean_obsnme = obsnme
#             for tag in ['trace', 'general', 'majorion', 'minorion', 'majorcation', 'majoranion', 'minoranion']:
#                 clean_obsnme = clean_obsnme.replace(tag, '')

#             # Extract tag (e.g., 'pH') for filename and legend
#             tag = clean_obsnme.split(":")[-1]
#             safe_tag = "".join(c for c in tag if c.isalnum() or c in (' ', '_', '-')).rstrip()

#             obsval = nzobs_all.loc[nzobs_all.obsnme == obsnme, "obsval"].iloc[0]
#             simvals = obsdict[itr][obsnme].dropna()

#             for real, simval in simvals.items():
#                 row = {"realization": real, "obsnme": clean_obsnme, "tag": safe_tag, "observed": obsval, "simulated": simval, "is_base": real == "base", "iteration": itr}
#                 csv_rows.append(row)

#         if not csv_rows:
#             continue

#         df_csv = pd.DataFrame(csv_rows)
#         df_pivot = df_csv.pivot(index="realization", columns="tag", values="simulated")
#         df_pivot["is_base"] = df_csv.drop_duplicates("realization").set_index("realization")["is_base"]
#         df_pivot["iteration"] = itr

#         # add observed values as separate row
#         obs_row = dict(df_pivot.iloc[0])  # copy structure
#         obs_row.update(df_csv.groupby("tag").first()["observed"].to_dict())
#         obs_row["realization"] = "observed"
#         obs_row["is_base"] = False
#         obs_row["iteration"] = itr
#         df_pivot.loc["observed"] = obs_row

#         # save to csv
#         csv_out = os.path.join(fdir, f"ensemble_obs_iter{itr}.csv")
#         df_pivot.to_csv(csv_out)
#         print(f"saved csv for iter {itr}: {csv_out}")

#     # plotting (one plot per obsnme)
#     for obsnme in all_obsnames:
#         clean_obsnme = obsnme
#         for tag in ['trace', 'general', 'majorion', 'minorion', 'majorcation', 'majoranion', 'minoranion']:
#             clean_obsnme = clean_obsnme.replace(tag, '')

#         tag = clean_obsnme.split(":")[-1]
#         safe_tag = "".join(c for c in tag if c.isalnum() or c in (' ', '_', '-')).rstrip()
#         pdf_path = os.path.join(fdir, f"{safe_tag}.pdf")

#         with PdfPages(pdf_path) as pdf:
#             plot_data = []
#             for itr in iters:
#                 if obsnme not in obsdict[itr].columns:
#                     continue

#                 obsval = nzobs_all.loc[nzobs_all.obsnme == obsnme, "obsval"].iloc[0]
#                 simvals = obsdict[itr][obsnme].dropna()

#                 for real, val in simvals.items():
#                     plot_data.append({"iteration": f"Iter {itr}", "realization": real, "sim_value": val, "obs_value": obsval, "is_base": real == "base"})

#             if not plot_data:
#                 continue

#             df_plot = pd.DataFrame(plot_data)

#             plt.figure(figsize=(10, 6))
#             sns.boxplot(x="iteration", y="sim_value", data=df_plot, color="lightblue", fliersize=0)
#             sns.stripplot(x="iteration", y="sim_value", data=df_plot[df_plot["is_base"] == False], color="blue", alpha=0.5, size=4, jitter=True, label="Realizations")

#             # plot base realization
#             df_base = df_plot[df_plot["is_base"] == True]
#             plt.scatter(df_base["iteration"], df_base["sim_value"], color="black", edgecolor="none", zorder=10, s=50, label="Base")

#             # plot observed value as red dot (once per iteration)
#             for itr in df_plot["iteration"].unique():
#                 obsval = df_plot[df_plot["iteration"] == itr]["obs_value"].iloc[0]
#                 plt.scatter(itr, obsval, edgecolors="yellow", facecolors='none', s=60, zorder=11, label="Observed")

#             plt.yscale("log")
#             plt.title(f"Log-scale Simulated Ensemble vs Observed: {safe_tag}")
#             plt.ylabel("Simulated Value (log scale)")
#             plt.xlabel("Iteration")
#             plt.grid(True, alpha=0.3)

#             # legend duplication
#             handles, labels = plt.gca().get_legend_handles_labels()
#             by_label = dict(zip(labels, handles))
#             plt.legend(by_label.values(), by_label.keys())

#             pdf.savefig(dpi=150)
#             plt.close()

#         print(f"saved to pdf: {pdf_path}")

#     return fdir

if __name__ == '__main__':

    m_d = os.path.join(os.getcwd(), 'arap_ies.2026.02.05.143244.5-60')
    modnm = 'AR01'
    optiter = 5
    realz = 60

    main(m_d, modnm, optiter, realz)
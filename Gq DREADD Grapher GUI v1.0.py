# pysimpleGUI for conversion of multiple Med Associate text files into one Excel file

# Import required libraries
import PySimpleGUI as sg
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import f_oneway
from statsmodels.stats.anova import anova_lm
import math
from pathlib import Path
import os

# validate that the file path is entered correctly
def is_valid_path(filepath):
    if filepath and Path(filepath).exists():
        return True
    sg.popup_error("A selected file path is incorrect or has been left empty.")
    return False

# window appears when the program successfully completes
def nom_window():
    layout = [[sg.Text("\n"
    " All Systems Nominal  \n"
    "\n"
    "")]]
    window = sg.Window((""), layout, modal=True)
    choice = None
    while True:
        event, values = window.read()
        if event == "Exit" or event == sg.WIN_CLOSED:
            break
    window.close()
    
# Define the location of the directory
def graph_data_calc_p_vals(input_file):

    # Change the directory
    os.chdir(os.path.dirname(input_file))

    prog_bar_update_val = 0
    window["-Progress_BAR-"].update(max = 11, current_count=int(prog_bar_update_val + 1))
    
    df = pd.read_excel(input_file)

    # columns designating behaviors to be analyzed 
    tot_tm = 'Time_Sum'
    tot_bt = 'Bout_Sum'
    av_bts = 'Av_Bout_Tm'

    formula = f'{tot_tm} ~ Virus + Shock + Virus:Shock'
    model = ols(formula, df).fit()
    tot_tm_aov_table = anova_lm(model, typ=2)

    formula = f'{tot_bt} ~ Virus + Shock + Virus:Shock'
    model = ols(formula, df).fit()
    tot_bt_aov_table = anova_lm(model, typ=2)


    formula = f'{av_bts} ~ Virus + Shock + Virus:Shock'
    model = ols(formula, df).fit()
    av_bts_aov_table = anova_lm(model, typ=2)

    # dataframes and vars for virus group comparisions

    mcher_grp = df['Virus'] == ('Cn')
    virus_grp = df['Virus'] == ('Gq')

    mC = df[mcher_grp]
    Vi = df[virus_grp]

    vi_df_list = [mC,Vi]

    # dataframes and vars for shock group comparisions

    n_sck_grp = df['Shock'] == ('NS')
    y_sck_grp = df['Shock'] == ('FS')

    nS = df[n_sck_grp]
    fS = df[y_sck_grp]

    sk_df_list = [nS,fS]

    # dataframes and vars for all group comparisions

    virus_Gq_grp = df['Virus'] == 'Gq'
    virus_Cn_grp = df['Virus'] == 'Cn'
    shock_FS_grp = df['Shock'] == 'FS'
    shock_NS_grp = df['Shock'] == 'NS'


    Cn_NS = df[virus_Cn_grp & shock_NS_grp]
    Cn_FS = df[virus_Cn_grp & shock_FS_grp]
    Gq_NS = df[virus_Gq_grp & shock_NS_grp]
    Gq_FS = df[virus_Gq_grp & shock_FS_grp]

    df_list = [Cn_NS, Cn_FS, Gq_NS, Gq_FS]

    ######################### Time by Virus #########################################################

    vi_tm_means_list = []
    vi_tm_count_list = []
    vi_tm_std_d_list = []
    vi_tm_std_e_list = []

    for df in vi_df_list:
        vi_tm_means_list.append(df[tot_tm].mean())
        vi_tm_count_list.append(df[tot_tm].count())
        vi_tm_std_d_list.append(df[tot_tm].std())

    for count in vi_df_list:
        t_statistic, tm_vi_p_value = stats.ttest_ind(Vi[tot_tm], mC[tot_tm])

    for count, std_dev in zip(vi_tm_count_list, vi_tm_std_d_list):
        standard_error = std_dev / math.sqrt(count)
        vi_tm_std_e_list.append(standard_error)

    Groups = ['Gq Virus','mCherry']

    tm_virus_grps = vi_tm_means_list[0:2]
    tm_virus_sem = vi_tm_std_e_list[0:2]

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    ###################### Bouts by Virus ###########################################################

    vi_bt_means_list = []
    vi_bt_count_list = []
    vi_bt_std_d_list = []
    vi_bt_std_e_list = []

    for df in vi_df_list:
        vi_bt_means_list.append(df[tot_bt].mean())
        vi_bt_count_list.append(df[tot_bt].count())
        vi_bt_std_d_list.append(df[tot_bt].std())

    for count in vi_df_list:
        t_statistic, bt_vi_p_value = stats.ttest_ind(Vi[tot_bt], mC[tot_bt])

    for count, std_dev in zip(vi_bt_count_list, vi_bt_std_d_list):
        standard_error = std_dev / math.sqrt(count)
        vi_bt_std_e_list.append(standard_error)

    bt_virus_grps = vi_bt_means_list[0:2]
    bt_virus_sem = vi_bt_std_e_list[0:2]

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    ###################### AVG Bout Time by Virus ###################################################

    vi_av_means_list = []
    vi_av_count_list = []
    vi_av_std_d_list = []
    vi_av_std_e_list = []

    for df in vi_df_list:
        vi_av_means_list.append(df[av_bts].mean())
        vi_av_count_list.append(df[av_bts].count())
        vi_av_std_d_list.append(df[av_bts].std())

    for count in vi_df_list:
        t_statistic, av_vi_p_value = stats.ttest_ind(Vi[av_bts], mC[av_bts])

    # Extracting the data for the ANOVA
    data_for_anova = [df[av_bts] for df in vi_df_list]

    # Performing one-way ANOVA
    f_statistic, av_vi_p_value = f_oneway(*data_for_anova)

    for count, std_dev in zip(vi_av_count_list, vi_av_std_d_list):
        standard_error = std_dev / math.sqrt(count)
        vi_av_std_e_list.append(standard_error)

    av_virus_grps = vi_av_means_list[0:2]
    av_virus_sem = vi_av_std_e_list[0:2]

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    ######################### Time by Shock #########################################################

    sk_tm_means_list = []
    sk_tm_count_list = []
    sk_tm_std_d_list = []
    sk_tm_std_e_list = []

    for df in sk_df_list:
        sk_tm_means_list.append(df[tot_tm].mean())
        sk_tm_count_list.append(df[tot_tm].count())
        sk_tm_std_d_list.append(df[tot_tm].std())

    for count in sk_df_list:
        t_statistic, tm_sk_p_value = stats.ttest_ind(nS[tot_tm], fS[tot_tm])

    for count, std_dev in zip(sk_tm_count_list, sk_tm_std_d_list):
        standard_error = std_dev / math.sqrt(count)
        sk_tm_std_e_list.append(standard_error)

    tm_shock_grps = sk_tm_means_list[0:2]
    tm_shock_sem = sk_tm_std_e_list[0:2]

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    ###################### Bouts by Shock ###########################################################

    sk_bt_means_list = []
    sk_bt_count_list = []
    sk_bt_std_d_list = []
    sk_bt_std_e_list = []

    for df in sk_df_list:
        sk_bt_means_list.append(df[tot_bt].mean())
        sk_bt_count_list.append(df[tot_bt].count())
        sk_bt_std_d_list.append(df[tot_bt].std())

    for count in sk_df_list:
        t_statistic, bt_sk_p_value = stats.ttest_ind(nS[tot_bt], fS[tot_bt])

    for count, std_dev in zip(sk_bt_count_list, sk_bt_std_d_list):
        standard_error = std_dev / math.sqrt(count)
        sk_bt_std_e_list.append(standard_error)

    bt_shock_grps = sk_bt_means_list[0:2]
    bt_shock_sem = sk_bt_std_e_list[0:2]

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    ###################### AVG Bout Time by Shock ###################################################

    sk_av_means_list = []
    sk_av_count_list = []
    sk_av_std_d_list = []
    sk_av_std_e_list = []

    for df in sk_df_list:
        sk_av_means_list.append(df[av_bts].mean())
        sk_av_count_list.append(df[av_bts].count())
        sk_av_std_d_list.append(df[av_bts].std())

    for count in sk_df_list:
        t_statistic, av_sk_p_value = stats.ttest_ind(nS[av_bts], fS[av_bts])

    for count, std_dev in zip(sk_av_count_list, sk_av_std_d_list):
        standard_error = std_dev / math.sqrt(count)
        sk_av_std_e_list.append(standard_error)

    av_shock_grps = sk_av_means_list[0:2]
    av_shock_sem = sk_av_std_e_list[0:2]

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    ###################### Total Time Sum All Groups ############################

    tm_means_list = []
    tm_count_list = []
    tm_std_d_list = []
    tm_std_e_list = []

    for df in df_list:
        tm_means_list.append(df[tot_tm].mean())
        tm_count_list.append(df[tot_tm].count())
        tm_std_d_list.append(df[tot_tm].std())

    for count, std_dev in zip(tm_count_list, tm_std_d_list):
        standard_error = std_dev / math.sqrt(count)
        tm_std_e_list.append(standard_error)

    tm_Cn_grp = tm_means_list[0:2]
    tm_Tx_grp = tm_means_list[2:4]

    tm_Cn_sem = tm_std_e_list[0:2]
    tm_Tx_sem = tm_std_e_list[2:4]

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    ###################### Total Bouts Sum All Groups ###########################

    bt_means_list = []
    bt_count_list = []
    bt_std_d_list = []
    bt_std_e_list = []

    for df in df_list:
        bt_means_list.append(df[tot_bt].mean())
        bt_count_list.append(df[tot_bt].count())
        bt_std_d_list.append(df[tot_bt].std())

    for count, std_dev in zip(bt_count_list, bt_std_d_list):
        standard_error = std_dev / math.sqrt(count)
        bt_std_e_list.append(standard_error)

    bt_Cn_grp = bt_means_list[0:2]
    bt_Tx_grp = bt_means_list[2:4]

    bt_Cn_sem = bt_std_e_list[0:2]
    bt_Tx_sem = bt_std_e_list[2:4]

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    ###################### Total Bout AVG All Groups ##########################

    av_means_list = []
    av_count_list = []
    av_std_d_list = []
    av_std_e_list = []

    for df in df_list:
        av_means_list.append(df[av_bts].mean())
        av_count_list.append(df[av_bts].count())
        av_std_d_list.append(df[av_bts].std())

    for count, std_dev in zip(av_count_list, av_std_d_list):
        standard_error = std_dev / math.sqrt(count)
        av_std_e_list.append(standard_error)

    av_Cn_grp = av_means_list[0:2]
    av_Tx_grp = av_means_list[2:4]

    av_Cn_sem = av_std_e_list[0:2]
    av_Tx_sem = av_std_e_list[2:4]

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    #################### Graph Logic ###############################################################

    Vi_Groups = ['mCherry', 'Gq Virus']
    Sk_Groups = ['No Shock','Ft Shock']

    Cn_Groups = ['Cn_NS', 'Cn_FS']
    Tx_Groups = ['Gq_NS', 'Gq_FS']

    # Create subplots
    fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(12, 9.5))

    # Plot the first row bar graph
    axes[0,0].bar(Vi_Groups, tm_virus_grps, yerr=tm_virus_sem, capsize=5, label='Cn', alpha=0.7)
    axes[0,0].set_title('Total Time')
    # Add the p-value as text to the plot
    axes[0,0].text(0.5, max(tm_virus_grps), f'P = {tm_vi_p_value:.4f}', ha='center', va='bottom')

    # Plot the second bar graph
    axes[0,1].bar(Vi_Groups, bt_virus_grps, yerr=bt_virus_sem, capsize=5, label='Cn', alpha=0.7)
    axes[0,1].set_title('Bouts')
    axes[0,1].text(0.5, max(bt_virus_grps), f'P = {bt_vi_p_value:.4f}', ha='center', va='bottom')

    # Plot the third bar graph
    axes[0,2].bar(Vi_Groups, av_virus_grps, yerr=av_virus_sem, capsize=5, label='Cn', alpha=0.7)
    axes[0,2].set_title('Avg Bout Time')
    axes[0,2].text(0.5, max(av_virus_grps), f'P = {av_vi_p_value:.4f}', ha='center', va='bottom')

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))
    prog_bar_update_val += 1

    #----------------------------------------------------------------------------------

    # Plot the second row bar graph
    axes[1,0].bar(Sk_Groups, tm_shock_grps, yerr=tm_shock_sem, capsize=5, label='Cn', alpha=0.7)
    axes[1,0].set_title('Total Time')
    axes[1,0].text(0.5, max(tm_virus_grps), f'P = {tm_sk_p_value:.4f}', ha='center', va='bottom')

    # Plot the second bar graph
    axes[1,1].bar(Sk_Groups, bt_shock_grps, yerr=bt_shock_sem, capsize=5, label='Cn', alpha=0.7)
    axes[1,1].set_title('Bouts')
    axes[1,1].text(0.5, max(bt_shock_grps), f'P = {bt_sk_p_value:.4f}', ha='center', va='bottom')

    # Plot the third bar graph
    axes[1,2].bar(Sk_Groups, av_shock_grps, yerr=av_shock_sem, capsize=5, label='Cn', alpha=0.7)
    axes[1,2].set_title('Avg Bout Time')
    axes[1,2].text(0.5, max(av_shock_grps), f'P = {av_sk_p_value:.4f}', ha='center', va='bottom')

    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val))

    #----------------------------------------------------------------------------------

    # Plot the third row bar graph

    # Plot the first bar graph
    axes[2,0].bar(Cn_Groups, tm_Cn_grp, yerr=tm_Cn_sem, capsize=5, label='Cn', alpha=0.7)
    axes[2,0].bar(Tx_Groups, tm_Tx_grp, yerr=tm_Tx_sem, capsize=5, label='Gq', alpha=0.7)
    axes[2,0].set_title('Total Time')
    # Display ANOVA table beneath the bar graph
    axes[3,0].axis('off')  # Hide axis for the table
    tm_table_data = tot_tm_aov_table.round(4).to_numpy()
    tm_table_col_labels = tot_tm_aov_table.columns
    tm_table_row_labels = tot_tm_aov_table.index

    tm_table = axes[3,0].table(cellText=tm_table_data,
                        colLabels=tm_table_col_labels,
                        rowLabels=tm_table_row_labels,
                        loc='center')

    tm_table.auto_set_font_size(False)
    tm_table.set_fontsize(8)
    tm_table.scale(0.8, 1)

    # Plot the second bar graph
    axes[2,1].bar(Cn_Groups, bt_Cn_grp, yerr=bt_Cn_sem, capsize=5, label='Cn', alpha=0.7)
    axes[2,1].bar(Tx_Groups, bt_Tx_grp, yerr=bt_Tx_sem, capsize=5, label='Gq', alpha=0.7)
    axes[2,1].set_title('Bouts')
    # Display ANOVA table beneath the bar graph
    axes[3,1].axis('off')  # Hide axis for the table
    bt_table_data = tot_bt_aov_table.round(4).to_numpy()
    bt_table_col_labels = tot_bt_aov_table.columns
    bt_table_row_labels = tot_bt_aov_table.index

    bt_table = axes[3,1].table(cellText=bt_table_data,
                        colLabels=bt_table_col_labels,
                        rowLabels=bt_table_row_labels,
                        loc='center')

    bt_table.auto_set_font_size(False)
    bt_table.set_fontsize(8)
    bt_table.scale(0.8, 1)

    # Plot the third bar graph
    axes[2,2].bar(Cn_Groups, av_Cn_grp, yerr=av_Cn_sem, capsize=5, label='Cn', alpha=0.7)
    axes[2,2].bar(Tx_Groups, av_Tx_grp, yerr=av_Tx_sem, capsize=5, label='Gq', alpha=0.7)
    axes[2,2].set_title('Avg Bout Time')
    # Display ANOVA table beneath the bar graph
    axes[3,2].axis('off')  # Hide axis for the table
    av_table_data = av_bts_aov_table.round(4).to_numpy()
    av_table_col_labels = av_bts_aov_table.columns
    av_table_row_labels = av_bts_aov_table.index

    av_table = axes[3,2].table(cellText=av_table_data,
                        colLabels=av_table_col_labels,
                        rowLabels=av_table_row_labels,
                        loc='center')

    av_table.auto_set_font_size(False)
    av_table.set_fontsize(8)
    av_table.scale(0.8, 1)

    prog_bar_update_val += 1
    window["-Progress_BAR-"].update(current_count=int(prog_bar_update_val + 1))

    # Adjust layout
    plt.tight_layout()

    plt.show()

    # window telling the user the program functioned correctly
    nom_window()   

# main GUI creation and GUI elements
sg.theme('DarkBlue7')

layout = [
    [sg.Text("\nSelect the Excel file containing the\n"
             "data to be graphed and analyzed.\n"),
    sg.Input(key="-IN-"),
    sg.FileBrowse()],

    [sg.Exit(), sg.Button("Press to graph behavioral data with p-values"), 
    sg.Text("eBot's progress..."),
    sg.ProgressBar(20, orientation='horizontal', size=(15,10), 
                border_width=4, bar_color=("Blue", "Grey"),
                key="-Progress_BAR-")]
    
]

# create the window
window = sg.Window("Welcome to eBot's behavioral data grapher and analyzer!", layout)

# create an event loop
while True:
    event, values = window.read()
    # end program if user closes window
    if event == "Exit" or event == sg.WIN_CLOSED:
        break
    if event == "Press to graph behavioral data with p-values":
        # check file selections are valid
        if (is_valid_path(values["-IN-"])):

            graph_data_calc_p_vals(
            input_file  = values["-IN-"])   

window.close
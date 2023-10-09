import sys
sys.path.append('/home/kyle/microbial_ecology/custom_functions/')
import nitrite_toxicity_model as ntm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import bmgdata as bd
import glob
from matplotlib.pyplot import cm
import warnings
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as ticker
from pprint import pprint
#warnings.filterwarnings("ignore")

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['text.usetex'] = False
plt.rcParams['font.size'] = 8

pHs = [6, 7.3]
#plt.style.use('dark_background')
#fig, axs = plt.subplots(2,1,figsize=(10,10), sharex = True, sharey = True)
#fig.tight_layout()
#now look at OD data
met_file = '/home/kyle/microbial_ecology/pathway_splitting/denit_experiments_analysis/CRM_predict/CRM_predict_cycle2_experiments.pkl'
fit_file = '/home/kyle/microbial_ecology/pathway_splitting/denit_experiments_analysis/CRM_predict/fits/NarG_fit_pH=6_no_offset_pH7_yields.npz'

xmin_met = -3
xmax_met = 72
ymin_met = -0.2
ymax_met = 2.1

pfit = np.load(fit_file)['pfit']
print(pfit)
ts = np.linspace(0,75, 256)
fig, axs = plt.subplots(2,2,figsize=(6,3.5), sharex = True, sharey = False)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['text.usetex'] = False
plt.rcParams['font.size'] = 8
plt.rcParams['font.family'] = 'Arial'

#A_marker = '^'
#I_marker = 'v'
A_marker = 'o'
I_marker = 'o'
linewidth = 2
alpha = 0.5
marker_2 = '>'
marker_1 = '<'
generic_marker = 's'
markersize = 5
color_1 = 'deepskyblue'
color_2 = 'tab:blue'
errorbar_color = 'Navy'
inset_width = "20%"
inset_height = "65%"
inset_pad = 0.8
#fig.tight_layout()
#now look at OD data
'''
filenames = glob.glob('/home/kyle/microbial_ecology/pathway_splitting/data/FF/KC_OD600_*_FF_cycle*_OD_endpoint_*.csv')
meta_filename = '/home/kyle/microbial_ecology/pathway_splitting/data/FF/plate1_metadata.csv'
meta = pd.read_csv(meta_filename,index_col=0).dropna(how='all')  #import metadata
pseudo_f0s = np.setdiff1d(np.unique(meta["pseudo_f0"]),['O2_control', 'blank', 'presumed_blank'])
'''
#print(pseudo_f0s)
cmap = cm.get_cmap('viridis')
#colors = ['tab:blue', 'tab:orange']
strains = ['NarG', 'NapA']
I0s = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
print('True I0s are ' + str(np.asarray(I0s)*7.0/8.0))

all_experiments = pd.read_pickle(met_file)
for experiments in all_experiments:
    for experiment in experiments:
        #pprint(vars(experiment))
        #print(len(experiment.A))#, axis = 0))
        if experiment.I0==0 and experiment.ID == 'NarG' and experiment.pH == 6:
            color = 'tab:blue'
            y0 = [experiment.N0, 0, np.mean(experiment.A, axis = 0)[0], np.mean(experiment.I, axis = 0)[0]]
            yh = ntm.denitODE(y0,ts,pfit,1)
            N_1 = yh[:,0]                                                                                                                                     
            N_1_dead = yh[:,1]
            A_vals = yh[:,2]
            I_vals = yh[:,3]
            #axs[0][0].plot(ts+experiment.t[0], A_vals, linewidth = linewidth, color = 'tab:blue', alpha = alpha)#, label = 'mono fit pred.')
            if experiment.A0 == 2:
                axs[0][0].plot(ts+experiment.t[0], A_vals,':', linewidth = linewidth, color = color_2, alpha = alpha+0.2, label='pred.')#, label = 'mono fit pred.')
                axs[0][1].errorbar(experiment.t, np.mean(experiment.A, axis = 0), yerr = np.std(experiment.A, axis=0), linewidth = linewidth, marker = A_marker, markersize = markersize, color = color_2, ecolor = color_2, label = 'PD Nar+', ls = 'none')
                axs[0][1].plot(experiment.t, np.mean(experiment.A, axis = 0), linewidth = linewidth, color = color_2, alpha = alpha)
                axs[0][0].errorbar(experiment.t, np.mean(experiment.A, axis = 0), yerr = np.std(experiment.A, axis=0), linewidth = linewidth, marker = A_marker, markersize = markersize, color = color_2, ecolor = color_2, ls = 'none')#label = 'NO3_0='+str(7.0*experiment.A0/8.0), ls = 'none')
                axs[1][0].plot(ts+experiment.t[0], I_vals,':', linewidth = linewidth, color =  color_2, alpha = alpha+0.2)#, label='pred.')#, label = 'mono fit pred.')
                axs[1][0].errorbar(experiment.t, np.mean(experiment.I, axis = 0), yerr = np.std(experiment.I, axis=0), linewidth = linewidth, marker = I_marker, markersize = markersize, color = color_2, ecolor = color_2, ls = 'none')#, alpha = alpha)#, label = 'NarG')
                axs[1][1].errorbar(experiment.t, np.mean(experiment.I, axis = 0), yerr = np.std(experiment.I, axis=0), linewidth = linewidth, marker = I_marker, markersize = markersize, color = color_2, ecolor = color_2, ls = 'none')#, alpha = alpha)#, label = 'NarG') 
                axs[1][1].plot(experiment.t, np.mean(experiment.I, axis = 0),  linewidth = linewidth, color = color_2, alpha = alpha)#, alpha = alpha)#, label = 'NarG') 
                inset_ax_1 = inset_axes(axs[1][0],
                    width=inset_width, # width = 30% of parent_bbox
                    height=inset_height, # height : 1 inch
                                        loc = 'upper right', borderpad = inset_pad, bbox_to_anchor = (0, 0, 1, 0.94), bbox_transform = axs[1][0].transAxes)
                
                inset_ax_1.bar(1, np.mean(experiment.Nend - experiment.N0), color = color_2, width = 0.5)
                print('dOD 2 mM is ' + str(np.mean(experiment.Nend - experiment.N0)))
                inset_ax_1.errorbar(1, np.mean(experiment.Nend - experiment.N0), yerr = np.std(experiment.Nend - experiment.N0), color =errorbar_color, ecolor =errorbar_color, marker = 'None', capsize = 2)
                inset_ax_2 = inset_axes(axs[1][1],
                    width=inset_width, # width = 30% of parent_bbox
                    height=inset_height, # height : 1 inch
                                        loc = 'upper right', borderpad = inset_pad, bbox_to_anchor = (0, 0, 1, 0.94), bbox_transform = axs[1][1].transAxes)
                
                
                inset_ax_2.bar(1, np.mean(experiment.Nend - experiment.N0), color = color_2, width = 0.5)
                inset_ax_2.errorbar(1, np.mean(experiment.Nend - experiment.N0), yerr = np.std(experiment.Nend - experiment.N0), color =errorbar_color, ecolor =errorbar_color, marker = 'None', capsize = 2)
            elif experiment.A0 == 1:
                axs[0][0].plot(ts+experiment.t[0], A_vals,'-', linewidth = linewidth, color = color_1, alpha = alpha, label=' fit')#, label = 'mono fit pred.')
                axs[1][0].plot(ts+experiment.t[0], I_vals,'-', linewidth = linewidth, color = color_1, alpha = alpha)#, label='NO3_0='+str(7.0*experiment.A0/8.0)+' fit')#, label = 'mono fit pred.')
                axs[1][0].errorbar(experiment.t, np.mean(experiment.I, axis = 0), yerr = np.std(experiment.I, axis=0), linewidth = 0, marker = I_marker, markersize = markersize, color = color_1, ecolor = color_1)
                axs[0][0].errorbar(experiment.t, np.mean(experiment.A, axis = 0), yerr = np.std(experiment.A, axis=0), linewidth = 0, marker = A_marker, markersize = markersize, color = color_1, ecolor = color_1, label='NO3_0='+str(7.0*experiment.A0/8.0))
                '''
                inset_ax_1 = inset_axes(axs[1][0],
                    width="20%", # width = 30% of parent_bbox
                    height=0.6, # height : 1 inch
                                      loc=1)
                '''
                inset_ax_1.bar(0, np.mean(experiment.Nend - experiment.N0), color = color_1, width = 0.5)
                print('dOD 1 mM is ' + str(np.mean(experiment.Nend - experiment.N0)))
                #print(experiment.Nend - experiment.N0)
                inset_ax_1.errorbar(0, np.mean(experiment.Nend - experiment.N0), yerr = np.std(experiment.Nend - experiment.N0), color = errorbar_color, ecolor =errorbar_color, marker = 'None', capsize = 2)
                #inset_ax_1.errorbar(0, np.mean(experiment.Nend - experiment.N0), yerr = np.std(experiment.Nend - experiment.N0), color = color_1, ecolor = color_1, marker = marker_1)

        if experiment.I0==0 and experiment.A0 == 2 and experiment.ID == 'NarG+NapA' and experiment.pH == 6:
                axs[1][1].errorbar(experiment.t, np.mean(experiment.I, axis = 0), yerr = np.std(experiment.I, axis=0), linewidth = linewidth, marker = I_marker, markersize = markersize, color = 'tab:purple', ecolor = 'tab:purple', ls = 'none')#, label = 'NarG+NapA')          
                axs[0][1].errorbar(experiment.t, np.mean(experiment.A, axis = 0), yerr = np.std(experiment.A, axis=0), linewidth = linewidth, marker = A_marker, markersize = markersize, color = 'tab:purple', ecolor = 'tab:purple', label = '1:1 coculture', ls = 'none')          
                axs[1][1].plot(experiment.t, np.mean(experiment.I, axis = 0), linewidth = linewidth, marker = I_marker, markersize = markersize, color = 'tab:purple', alpha = alpha)#, label = 'NarG+NapA')          
                axs[0][1].plot(experiment.t, np.mean(experiment.A, axis = 0), linewidth = linewidth, color = 'tab:purple', alpha = alpha)          
                inset_ax_2.bar(0, np.mean(experiment.Nend - experiment.N0), color = 'tab:purple', width = 0.5)
                inset_ax_2.errorbar(0, np.mean(experiment.Nend - experiment.N0), yerr = np.std(experiment.Nend - experiment.N0), color =errorbar_color, ecolor =errorbar_color, marker = 'None', capsize = 2)

#axs[1][2].legend(title = 'NO2_0 [mM]')
inset_ax_2.set_ylim([0,0.07])
inset_ax_1.set_ylim([0,0.07])
axs[0][0].set_ylim([ymin_met, ymax_met])
axs[1][0].set_ylim([ymin_met, ymax_met])
axs[0][1].set_ylim([ymin_met, ymax_met])
axs[1][1].set_ylim([ymin_met, ymax_met])


axs[0][0].set_xlim([xmin_met, xmax_met])
axs[1][0].set_xlim([xmin_met, xmax_met])
axs[0][1].set_xlim([xmin_met, xmax_met])
axs[1][1].set_xlim([xmin_met, xmax_met])

inset_ax_1.set_ylabel('dOD')#, fontsize = 6)
inset_ax_2.set_ylabel('dOD')
inset_ax_1.set_xlim([-0.5, 1.5])
inset_ax_2.set_xlim([-0.5, 1.5])
inset_ax_1.set_xticks([0,1])
inset_ax_1.tick_params(axis='y', pad=2)
inset_ax_2.tick_params(axis='y', pad=2)
inset_ax_1.set_yticks([0,0.05])
inset_ax_2.set_yticks([0,0.05])
#inset_ax_1.set_xticklabels(['', 'NO30=2'], rotation= 90)#, ha='right')
inset_ax_1.set_xticklabels(['', ''], rotation= 90)#, ha='right')
inset_ax_2.set_xticklabels(['', ''], rotation= 90)#, ha='right')
inset_ax_2.set_xticks([0,1])
#inset_ax_2.set_xticklabels(['1:1 co.', 'PD Nar+'], rotation= 90)#, ha='right')
inset_ax_1.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))#, fontsize = 6)
inset_ax_2.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))#, fontsize = 6)
#inset_ax_1.xtickangle(45)

plt.subplots_adjust(left=0.1, bottom=0.12, right=0.9, top=0.93, wspace=0.4, hspace=0.06)
axs[0][0].legend(prop={'size': 6})
axs[0][0].set_ylabel('NO3 [mM]')
axs[0][1].set_ylabel('NO3 [mM]')
axs[0][0].set_title('PD Nar+')
axs[1][0].set_xlabel('t [h]')
axs[1][1].set_xlabel('t [h]')
axs[0][1].set_title('coculture comparison')
axs[0][1].legend(prop={'size': 6})
axs[1][0].set_ylabel('NO2 [mM]')
axs[1][1].set_ylabel('NO2 [mM]')

fig.add_subplot(111, frameon=False)

# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
#plt.title('pH 7.3')
#plt.xlabel("t [h]", fontsize = 8)
#plt.ylabel(r"met. conc.", fontsize = 26)
#plt.tight_layout()
#plt.legend()
#plt.text(0.43,0.27, 'Killing curve', fontsize = 9)
#plt.text(1.01,0.66, 'NO2_0', fontsize = 9)
plt.savefig('fig_4_met_dynamics.svg')
plt.show()
#plt.cla()

###############
####Panel B#####
################
#print('next experiment')

fig2, axs2 = plt.subplots(2,2,figsize=(6,3.5), sharex = True, sharey = True)

for experiments in all_experiments:
    for experiment in experiments:
        #print(len(experiment.A))#, axis = 0))
        for i in range(len(strains)):
            for j in range(len(I0s)):
                if experiment.I0==I0s[j] and experiment.A0 == 2 and experiment.ID == strains[i] and experiment.pH == 6:
                    color_val = cmap((I0s[j]*7.0/8.0)/(2))
                    axs2[0][i].errorbar(experiment.t, np.mean(experiment.A, axis = 0), yerr = np.std(experiment.A, axis=0), linewidth = linewidth, marker = A_marker, markersize = markersize, color = color_val, ecolor = color_val, label = str(round(I0s[j]*7.0/8.0,3)), ls = 'none')#, markersize = markersize)          
                    axs2[1][i].errorbar(experiment.t, np.mean(experiment.I, axis = 0), yerr = np.std(experiment.I, axis=0), linewidth = linewidth, marker = I_marker, markersize = markersize, color = color_val, ecolor = color_val, label = str(round(I0s[j]*7.0/8.0,3)), ls = 'none')#, markersize = markersize)          

                    axs2[0][i].plot(experiment.t, np.mean(experiment.A, axis = 0), linewidth = linewidth,  color = color_val, alpha = alpha)#, markersize = markersize)          
                    axs2[1][i].plot(experiment.t, np.mean(experiment.I, axis = 0),  linewidth = linewidth,  color = color_val, alpha = alpha)#, markersize = markersize)          

plt.subplots_adjust(left=0.1, bottom=0.12, right=0.9, top=0.93, wspace=0.04, hspace=0.06)
axs2[0][0].set_title('PD Nar+')
axs2[0][1].set_title('RH Nap+')
axs2[0][1].legend(title = 'NO2_0 [mM]', ncol = 2, fontsize = 6)
axs2[0][0].set_ylabel('NO3 [mM]')
axs2[1][0].set_ylabel('NO2 [mM]')
axs[0][0].set_ylim([ymin_met, ymax_met])
axs[1][0].set_ylim([ymin_met, ymax_met])
axs[0][1].set_ylim([ymin_met, ymax_met])
axs[1][1].set_ylim([ymin_met, ymax_met])


axs[0][0].set_xlim([xmin_met, xmax_met])
axs[1][0].set_xlim([xmin_met, xmax_met])
axs[0][1].set_xlim([xmin_met, xmax_met])
axs[1][1].set_xlim([xmin_met, xmax_met])

fig2.add_subplot(111, frameon=False)

# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
#plt.title('pH 7.3')
plt.xlabel("t [h]", fontsize = 8)
#plt.ylabel(r"met. conc.", fontsize = 26)
#fig.tight_layout()
#plt.legend()
#plt.text(0.43,0.27, 'Killing curve', fontsize = 9)
#plt.text(1.01,0.66, 'NO2_0', fontsize = 9)
plt.savefig('fig_4_met_dynamics_tox.svg')
plt.show()
#plt.cla()

fig, axs = plt.subplots(2,1,figsize=(3.5,3.5), sharex = True, sharey = True)

cfu_file = np.load('/home/kyle/microbial_ecology/pathway_splitting/denit_experiments_analysis/MCW_drex/processed_CFU_data.npz')
ts_cfus = cfu_file['ts']
cfus = cfu_file['cfus']
cfu_stds = cfu_file['cfu_stds']
strain_names_cfu = cfu_file['strains']
Is_cfus = cfu_file['Is']
markersize = 5
for i in range(len(strain_names_cfu)):
    for j in range(len(Is_cfus)):
        if Is_cfus[j] == 0.5:
            pass
        else:
            avgs = cfus[i,j,:]
            stds = cfu_stds[i,j,:]
            coef = np.polyfit(ts_cfus[~np.isnan(avgs)],np.log(avgs[~np.isnan(avgs)]),1)
            poly1d_fn = np.poly1d(coef)
            color_val = cmap((Is_cfus[j]*7.0/8.0)/(2))
            if strain_names_cfu[i] == 'NarG':
                axs[0].plot(ts_cfus, np.exp(poly1d_fn(ts_cfus)), '-', color = color_val, linewidth = linewidth, alpha = alpha)#, label = str(round(coef[0],3))+'t+'+str(round(coef[1],2))) #'--k'=black dashed line, 'yo' = yellow circle marker
                axs[0].errorbar(ts_cfus, avgs, yerr = stds, ls = 'none',  marker = 'o', markersize = markersize, color = color_val, ecolor = color_val, linewidth = linewidth)#, label = strain_names[i]+' mon')
                axs[0].set_yscale('log')
                #axs[0].set_ylabel('CFU/mL')
            elif strain_names_cfu[i] == 'NapA':
                axs[1].plot(ts_cfus, np.exp(poly1d_fn(ts_cfus)), '-', color = color_val, linewidth = linewidth, alpha = alpha)#, label = str(round(coef[0],3))+'t+'+str(round(coef[1],2))) #'--k'=black dashed line, 'yo' = yellow circle marker
                axs[1].errorbar(ts_cfus, avgs, yerr = stds, ls = 'none',  marker = 'o', markersize = markersize, color = color_val, ecolor = color_val, linewidth = linewidth, label = str(Is_cfus[j]*7.0/8.0))#, label = strain_names[i]+' mon')
                axs[1].set_yscale('log')

axs[1].legend(title = 'NO2_0 [mM]')

axs[0].set_title('PD Nar+')
axs[1].set_title('RH Nap+')
#plt.suptitle('pH 6')
#plt.subplots_adjust(left=0.15, bottom=0.17, right=0.9, top=0.83, wspace=0.04, hspace=0.)
#plt.colorbar(sm)
#ax = axs[1][3]
#cax = ax.inset_axes([1.04, 0, 0.05, 2.45])
#cbar = fig.colorbar(sm, ax=axs[1:, :], cax = cax)#, shrink=2)
#fig.colorbar(sm, ax=axs[1:, :], , shrink=2)
#cbar.set_clim(0,2)
plt.suptitle('Killing curves')
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
#plt.title('Killing curve')
plt.xlabel("t [h]", fontsize = 8)
plt.ylabel('CFU/mL', fontsize = 8)
#plt.ylabel('CFU/mL', fontsize = 8)
#plt.ylabel(r"met. conc.", fontsize = 26)
fig.tight_layout()
#plt.legend()
#plt.text(0.43,0.27, 'Killing curve', fontsize = 9)
#plt.text(1.01,0.66, 'NO2_0', fontsize = 9)
plt.savefig('fig_4_killing_curve.svg')
plt.show()
'''
sm = plt.cm.ScalarMappable(cmap=cmap)
sm.set_clim(vmin=0, vmax=2)

plt.colorbar(sm, label = 'NO2_0 [mM]', aspect = 20)
plt.savefig('fig_4_colorbar.svg')
plt.show()
'''

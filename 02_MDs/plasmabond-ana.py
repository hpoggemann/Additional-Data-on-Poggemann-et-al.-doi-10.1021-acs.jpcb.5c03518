# Script for anaylsis of plasma species interactions with the enyzme atoms
# The script analysis the lammps trajectory file using the ovito python api
# Plot function might need to be adapted to the personal folder structure
## 3 subfolders numbered 01, 02, and 03 are expected and the results averaged

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import warnings
warnings.filterwarnings("ignore")
import ovito
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.modifiers import LoadTrajectoryModifier
import csv

# PLASAMA BOND ANALYSIS

def plasmaspecies_bonds(path, file, residuelist, Protein_ID, Plasma_ID, frame=[], file_number=[]):
    # Protein ID is the ID of the last protein atoms 
    # Plasma_ID is the ID of the first plasma species atoms
    # Frame tells the programm the frames to analyze
    # File_number is a string number the output files
    pipeline = import_file(path+file, sort_particles=True)
    # Create bonds:
    mod = CreateBondsModifier()
    mod.mode = CreateBondsModifier.Mode.Pairwise
    pipeline.modifiers.append(mod)
    mod.set_pairwise_cutoff('C', 'C', 2.040000057220459)
    mod.set_pairwise_cutoff('Fe', 'C', 2.2200000286102295)
    mod.set_pairwise_cutoff('Fe', 'Fe', 2.4)
    mod.set_pairwise_cutoff('H', 'C', 1.740000057220459)
    mod.set_pairwise_cutoff('H', 'Fe', 1.9200000286102294)
    mod.set_pairwise_cutoff('Mg', 'C', 1.727999997138977)
    mod.set_pairwise_cutoff('Mg', 'Fe', 1.9079999685287474)
    mod.set_pairwise_cutoff('Mg', 'H', 1.427999997138977)
    mod.set_pairwise_cutoff('Mg', 'Mg', 1.4159999370574952)
    mod.set_pairwise_cutoff('N', 'C', 1.95)
    mod.set_pairwise_cutoff('N', 'Fe', 2.1299999713897706)
    mod.set_pairwise_cutoff('N', 'H', 1.65)
    mod.set_pairwise_cutoff('N', 'Mg', 1.637999939918518)
    mod.set_pairwise_cutoff('N', 'N', 1.8599999427795408)
    mod.set_pairwise_cutoff('O', 'C', 1.9320000171661376)
    mod.set_pairwise_cutoff('O', 'Fe', 2.111999988555908)
    mod.set_pairwise_cutoff('O', 'H', 1.6320000171661377)
    mod.set_pairwise_cutoff('O', 'Mg', 1.6199999570846557)
    mod.set_pairwise_cutoff('O', 'N', 1.8419999599456787)
    mod.set_pairwise_cutoff('O', 'O', 1.8239999771118163)
    mod.set_pairwise_cutoff('S', 'C', 2.1)
    mod.set_pairwise_cutoff('S', 'Fe', 2.2799999713897705)
    mod.set_pairwise_cutoff('S', 'H', 1.7999999999999998)
    mod.set_pairwise_cutoff('S', 'Mg', 1.787999939918518)
    mod.set_pairwise_cutoff('S', 'N', 2.0099999427795407)
    mod.set_pairwise_cutoff('S', 'O', 1.9919999599456786)
    mod.set_pairwise_cutoff('S', 'S', 2.159999942779541)
    data = pipeline.compute(frame)
    print("Number of bonds:", data.particles.bonds.count)
    # Pint all bonds between protein and plasma species into list, 0 protein ID, 1 species ID
    bonds ={'Protein_ID':[], 'Type':[], 'Residue':[], 'Plasma_ID':[]}
    for a,b in data.particles.bonds.topology:
        if a<=Protein_ID and b>=Plasma_ID:
            bonds['Protein_ID'].append(a)
            bonds['Plasma_ID'].append(b)
    # Get Residue number and Paricle Type number for protein atoms
    for i in bonds['Protein_ID']:
        bonds['Type'].append(data.particles['Particle Type'].type_by_id(data.particles['Particle Type'][i]).name)
        bonds['Residue'].append(data.particles['Molecule Identifier'][i])
    bonds = pd.DataFrame(bonds)
    # Read residuelist to give residues a human readable name    
    res = pd.read_csv(residuelist, delim_whitespace=True,)
    for item in bonds['Residue']:
        loc= res.isin([item]).any(axis=1).idxmax()
        resn = res.iloc[loc]['ResType']
        bonds['Residue'].replace(item, resn, inplace=True)
    # Save output in csv file
    bonds.to_csv(path+'%s_plasmaspecies_bonds.txt'%file_number, sep='\t', index=False)
    # Count all res and attacked res
    res_count = res['ResType'].value_counts()
    bonds_count = bonds['Residue'].value_counts()
    # Make percentage
    result = {'labels': [], 'total':[], 'percent': [],}
    for item in bonds_count.index:
        result['percent'].append(float(bonds_count[item]/res_count[item]))
        result['total'].append(bonds_count[item])
        result['labels'].append(str(item))
    # Export attack amount vs residue to file
    pd.DataFrame(result).to_csv(path+'/interaction_probability.txt', sep='\t', index=False)
    result = pd.DataFrame(result).sort_values('labels', ascending=True)
    return bonds, result

# PLOT BOND PER ATOM ANALYSIS

def load_bond_file(path):
    bond_file = '1_plasmaspecies_bonds.txt'
    folders = ['01','02','03']
    for i in folders:
        if i=='01':
            df1 = pd.DataFrame(pd.read_csv(os.path.join(path+folder+'/'+i +'/'+bond_file), sep='\t', usecols = ['Type']).value_counts()).reset_index()
        elif i=='02':
            df2 = pd.DataFrame(pd.read_csv(os.path.join(path+folder+'/'+i +'/'+bond_file), sep='\t', usecols = ['Type']).value_counts()).reset_index()
        elif i=='03':
            df3 = pd.DataFrame(pd.read_csv(os.path.join(path+folder+'/'+i +'/'+bond_file), sep='\t', usecols = ['Type']).value_counts()).reset_index()
    all_atm = df1.merge(df2, how='outer', on=['Type'], suffixes=('1','2')).merge(df3, how='outer', on=['Type'],
                        suffixes=('1','3')).fillna(0).sort_values('Type', ascending=True).set_index('Type')
    mean = pd.DataFrame(all_atm.mean(axis=1), columns = ['mean'])
    all_atm = pd.merge(all_atm, mean, on=['Type'])
    all_atm.to_csv(path+folder+'/%s_mean_atom_interaction.txt'%folder, sep='\t')
    return all_atm
            
def atom_int_plot(path, atom_data):
    color={'C': 'grey', 'H': 'whitesmoke', 'N': 'b', 'O': 'r', 'S':'y', 'Fe':'orange', 'Mg': 'm'}
    # Plot attack amount vs residue
    fig, ax = plt.subplots(1)
    #fig.set_size_inches(10, 5)
    ## Plot bars ##
    ax.bar(atom_data.index,atom_data.mean(axis=1),edgecolor='black',color=[color[key] for key in atom_data.index])
    ax.errorbar(atom_data.index,atom_data.mean(axis=1), yerr = atom_data.std(axis=1), fmt='none', ecolor = 'k', elinewidth=1.5, capsize = 5)
    # Layout
    ax.set_ylabel('Total bond number',fontsize=18)
    ax.set_xlabel('Atom type',fontsize=18)
    ## Increase line thickness box and ticks      
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.tick_params(labelsize=16, width=1.5, length=6, bottom=True, left=True, right=True,)
    ## Minor ticks
    #ax.minorticks_on()
    ax.tick_params(which='minor',width=1, length=6, right=True,)
    # Save
    plt.savefig(path+folder+'/attack_total_vs_atom.png', dpi = 300, bbox_inches='tight')
    plt.show()


# PLOT BOND PER RESIDUE ANALYSIS    

def load_file(path):
    ## Percentage data
    file = 'interaction_probability.txt'
    folders = ['01','02','03']
    for i in folders:
        if i=='01':
            df1 = pd.read_csv(os.path.join(path+folder+'/'+i +'/'+file), sep='\t', usecols = [0,2]).sort_values('labels', ascending=True).reset_index(drop=True)
        elif i=='02':
            df2 = pd.read_csv(os.path.join(path+folder+'/'+i +'/'+file), sep='\t', usecols = [0,2]).sort_values('labels', ascending=True).reset_index(drop=True)
        elif i=='03':
            df3 = pd.read_csv(os.path.join(path+folder+'/'+i +'/'+file), sep='\t', usecols = [0,2]).sort_values('labels', ascending=True).reset_index(drop=True)

    full = df1.merge(df2, on=['labels'], how='outer',suffixes=('1','2')).merge(df3, on=['labels'],
                     how='outer',suffixes=('1','3')).fillna(0).sort_values('labels', ascending=True).set_index('labels')
    mean = pd.DataFrame(full.mean(axis=1), columns = ['mean'])
    full = pd.merge(full, mean, on=['labels'])
    full.to_csv(path+folder+'/%s_interaction_percent.txt'%folder, sep='\t')
    ## Total data
    file = 'interaction_probability.txt'
    folders = ['01','02','03']
    for i in folders:
        if i=='01':
            df1 = pd.read_csv(os.path.join(path+folder+'/'+i +'/'+file), sep='\t', usecols = [0,1]).sort_values('labels', ascending=True).reset_index(drop=True)
        elif i=='02':
            df2 = pd.read_csv(os.path.join(path+folder+'/'+i +'/'+file), sep='\t', usecols = [0,1]).sort_values('labels', ascending=True).reset_index(drop=True)
        elif i=='03':
            df3 = pd.read_csv(os.path.join(path+folder+'/'+i +'/'+file), sep='\t', usecols = [0,1]).sort_values('labels', ascending=True).reset_index(drop=True)
    full = df1.merge(df2, on=['labels'], how='outer',suffixes=('1','2')).merge(df3, on=['labels'],
                     how='outer',suffixes=('1','3')).fillna(0).sort_values('labels', ascending=True).set_index('labels')
    full.to_csv(path+folder+'/%s_mean_total_interaction.txt'%folder, sep='\t')
    ## Relativ data
    res_names=res['ResType'].reset_index(drop=True).value_counts()
    relativ = {'labels': [], 'percent':[], 'std': [],}
    for item in full.index:
        relativ['labels'].append(item)
        relativ['percent'].append(float((full.mean(axis = 1)[item]*res_names[item])/res.count()[0]))
        relativ['std'].append(float((full.std(axis = 1)[item]*res_names[item])/res.count()[0]))
    return full, relativ

def residue_analysis_total_plot(path, total_data, relativ_data, temp, folder):
    colors={'ALA': '#d6a090',
            'ARG': '#fe3b1e',
            'ASN': '#a12c32',
            'ASP': '#fa2f7a',
            'CYM': '#fb9fda',
            'CYS': '#e61cf7',
            'GLN': '#992f7c',
            'GLU': '#11963b',
            'GLY': '#051155',
            'HEM': '#4f02ec',
            'HISD': '#2d69cb',
            'HISE': '#00a6ee',
            'HISH': '#6febff',
            'ILE': '#08a29a',
            'LEU': '#2a666a',
            'LYS': '#063619',
            'MET': '#000000',
            'MG': '#4a4957',
            'PHE': '#8e7ba4',
            'PRO': '#b7c0ff',
            'SER': '#ffffff',
            'THR': '#acbe9c',
            'TRP': '#827c70',
            'TYR': '#5a3b1c',
            'VAL': '#ae6507'}
    labels={'ALA': 'Ala',
            'ARG': 'Arg',
            'ASN': 'Asn',
            'ASP': 'Asp',
            'CYM': 'Cym',
            'CYS': 'Cys',
            'GLN': 'Gln',
            'GLU': 'Glu',
            'GLY': 'Gly',
            'HEM': 'Heme',
            'HISD': 'Hisd',
            'HISE': 'Hise',
            'HISH': 'Hish',
            'ILE': 'Ile',
            'LEU': 'Leu',
            'LYS': 'Lys',
            'MET': 'Met',
            'MG': 'Mg',
            'PHE': 'Phe',
            'PRO': 'Pro',
            'SER': 'Ser',
            'THR': 'Thr',
            'TRP': 'Trp',
            'TYR': 'Tyr',
            'VAL': 'Val'}
    # plot attack amount vs residue
    fig, ax = plt.subplots(2, 1, sharex=True)
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)
    fig.set_size_inches(9,6)
    fig.supylabel('Interactions with rmdPGS',fontsize=24, x=-0.02, fontweight='bold')
    ax[1].set_ylim(0,49)
    ## Plot bars ##
    ax[1].bar(total_data.index,total_data.mean(axis = 1), edgecolor='black', color=[colors[key] for key in total_data.index],)
    ax[1].errorbar(total_data.index,total_data.mean(axis=1), yerr = total_data.std(axis=1), fmt='none', ecolor = 'k', elinewidth=1.5, capsize = 5)
    # Layout
    ax[1].set_ylabel('$n$ total bonds',fontsize=18,  labelpad=10)
    ax[1].set_xlabel('residue',fontsize=18)
    ## Increase line thickness box and ticks      
    for axis in ['top','bottom','left','right']:
        ax[1].spines[axis].set_linewidth(1.5)
    ax[1].tick_params(labelsize=16, width=1.5, length=6, bottom=True, left=True, right=True,)
    for tick in ax[1].get_xticklabels():
            tick.set_rotation(45)
    ax[1].set_xticklabels([labels[key] for key in total_data.index])
    ## Minor ticks
    ax[1].tick_params(which='minor',width=1, length=6, right=True,)
    ## Plot 2 ##
    ax[0].set_ylim(0,3)
    ax[0].bar( relativ_data['labels'], relativ_data['percent'], edgecolor='black', color=[colors[key] for key in total_data.index],)
    ax[0].errorbar( relativ_data['labels'], relativ_data['percent'], yerr =  relativ_data['std'], fmt='none', ecolor = 'k', elinewidth=1.5, capsize = 5)
    #layout
    ax[0].set_ylabel('relative',fontsize=18,  labelpad=20)
    ax[0].set_xlabel('residue',fontsize=18)
    ## increase line thickness box and ticks      
    for axis in ['top','bottom','left','right']:
        ax[0].spines[axis].set_linewidth(1.5)
    ax[0].tick_params(labelsize=16, width=1.5, length=6, bottom=True, left=True, right=True,)
    ax[0].tick_params(which='minor',width=1, length=6, right=True,)
    # Save
    plt.savefig(path+folder+'/%s_%s_attack_total_vs_residue_new.png'%(temp,folder), dpi = 300, bbox_inches='tight')
    plt.show()
    
# RUN PLASMA BOND ANALYSIS
path = r'path/to/files'
r = 'path/to/residuelistfile/residuelist.txt'
# Protein ID is the ID of the last protein atoms 
Protein_ID = 3719   #GapA: 5063; CviUPO: 3719; AaeUPO: 5040
# Plasma_ID is the ID of the first plasma species atoms
Plasma_ID = 3720  #GapA_vac: 5064; CviUPO_vac: 3720; AaeUPO_vac: 5041 
                    #GapA_sov: ; CviUPO_solv: ; AaeUPO_solv: 5321
    
# MD trajectory file
file = 'traj.lmp'
bonds, result = plasmaspecies_bonds(path, file, r, Protein_ID, Plasma_ID, frame=11, file_number=1)

# PLOT PLASMA BOND ANALYSIS OUTPUT
temp ='temp-of-the-run'
folder = 'foldername'
path = r'path/to/files'
res = pd.read_csv('path/to/residuelistfile/residuelist.txt', sep='\s+',)

## atom analyis
all_atm = load_bond_file(path)
atom_int_plot(path, all_atm)

## residue analysis
full, relativ = load_file(path)
residue_analysis_total_plot(path, full, relativ, temp, folder)

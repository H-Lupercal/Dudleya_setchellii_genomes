# Visualize PCAngsd covariance matrix as a 2D scatterplot of PCs (Principal Components)

# Evan Hackstadt
# Whittall Lab Dudleya Genomics
# 1 November 2025


# TODO: fix warning:
#   PerformanceWarning: DataFrame is highly fragmented.
#   This is usually the result of calling `frame.insert` many times, which has poor performance.
#   Consider joining all columns at once using pd.concat(axis=1) instead. 
#   To get a de-fragmented frame, use `newframe = frame.copy()

import os
import argparse
import numpy as np
import pandas as pd
from scipy.stats import chi2
import matplotlib.pyplot as plt
import seaborn as sns


# --- Primary Function ---

def plot_pca(cov_path, info_path, out_path, x_pc = 1, y_pc = 2):
    """
    Plots PCA graph from a covariance matrix

    Args:
        cov_path (str): Path to the covariance matrix (.cov)
        info_path (str): Path to the info file containing sample names (.info)
        out_path (str): Path to the output folder where plot will be saved
        x_pc (int): (optional) The principal component to plot on the x-axis (1 recommended)
        y_pc (int): (optional) The principal component to plot on the y-axis (2 recommended)

    Returns:
        None
    """
    
    # --- Process Data ---
    cov = np.loadtxt(cov_path)  # Reads estimated covariance matrix
    info = np.loadtxt(info_path, dtype=str)    # Reads sample names
    species_labels = []     # species e.g. DUSE
    pop_labels = []     # popcodes e.g. DAN
    Du_labels = []      # labels e.g. Du-78
    LP_labels = []      # labels e.g. LP_10

    # parse sample names assuming format: "POPCODE_LP_xxx_Du-xxx"
    for sample in info:
        # split on _LP_ to isolate popcode
        substrings1 = sample.split('_LP_')
        popcode = substrings1[0]
        pop_labels.append(popcode)
        
        # extract Du# and LP# from split on '_'
        Du_label = ""
        LP_num = ""
        
        substrings2 = sample.split('_')
        for i, s in enumerate(substrings2):
            Du_label = s if 'Du-' in s else ""
            if s == 'LP':
                LP_num = substrings2[i+1]
        
        Du_labels.append(Du_label)
        LP_labels.append(LP_num)
        
        # check if we have specified a non-setchelli (i.e. popcode = SPECIES_POPCODE)
        if '_' in popcode:
            popcode_substrings = popcode.split('_')
            species = popcode_substrings[0]
            if species in ['ABAB', 'ABBE', 'ABMU', 'CY']:
                if species == 'CY':
                    species_labels.append('DUCY')
                else:
                    species_labels.append('DUAB')
        else:
            species_labels.append('DUSE')

    # process .cov into eigenvectors & eigenvalues
    evals, evecs = np.linalg.eigh(cov)
    # reverse order of evecs and evals so we have PC1, PC2, ... (descending variance)
    evecs = evecs[:,::-1]
    evals = evals[::-1]
    evals_sum = np.sum(evals)

    # put our data into a DataFrame
    points = pd.DataFrame({   # each row is a point
        'species': species_labels,
        'pop': pop_labels,
        'Du#': Du_labels,
        'LP#': LP_labels
    })

    # also store the variance explained by each PC
    variance_explained = {}
    
    # add info for each PC to the DataFrame and dict
    for i in range(len(evals)):
        points[f'PC{i+1}'] = evecs[:,i]
        variance_explained[f'PC{i+1}'] = round(100*evals[i]/evals_sum, 2)

    # Points DF has cols: species, pop, point, PC1, PC2, PC3, ...

    # --- PCA Plot ---
    markers = {"DUSE": "o", "DUCY": "s", "DUAB": "X"}
    hue_order = ["DUSE", "DUCY", "DUAB"]
    ax = sns.scatterplot(data=points, x=f"PC{x_pc}", y=f"PC{y_pc}",
                    hue="species", style="species", 
                    hue_order=hue_order, markers=markers)                 # optional: legend=False
    ax.legend(fontsize=10)
    # sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set_title(f"{len(evals)}-Sample PCA", fontsize=14)
    ax.set_xlabel(f"PC{x_pc} ({variance_explained[f'PC{x_pc}']}% of variance)")
    ax.set_ylabel(f"PC{y_pc} ({variance_explained[f'PC{y_pc}']}% of variance)")
    
    # --- Save Plots ---
    
    # Unlabeled plot
    png_path = os.path.join(out_path, f"pca.png")   # consider: if allowing modular PCs, include in filename
    plt.savefig(png_path, dpi=600)
    
    # Label and plot again
    for x, y, label in zip(points[f"PC{x_pc}"], points[f"PC{y_pc}"], points["Du#"]):
        plt.annotate(label, xy=(x,y), fontsize=6)

    png_path_labeled = os.path.join(out_path, f"pca_labeled.png")   # consider: if allowing modular PCs, include in filename
    plt.savefig(png_path_labeled, dpi=600)

    # TEMP MANUAL ZOOM-IN ON CLUSTERS -- CONSIDER SOME WAY TO AUTOMATE THIS OR AT LEAST PASS IN AS PARAM
    '''
    # 267-sample zooms
    plt.xlim(-0.035, 0.005)
    plt.ylim(-0.005, 0.01)
    png_path_1 = os.path.join(out_path, f"pca_zoomed_1.png")
    plt.savefig(png_path_1, dpi=600)
    
    plt.xlim(0.10, 0.14)
    plt.ylim(-0.20, -0.07)
    png_path_2 = os.path.join(out_path, f"pca_zoomed_2.png")
    plt.savefig(png_path_2, dpi=600)
    
    plt.xlim(0.10, 0.14)
    plt.ylim(0.10, 0.15)
    png_path_3 = os.path.join(out_path, f"pca_zoomed_3.png")
    plt.savefig(png_path_3, dpi=600)
    
    plt.xlim(0.09, 0.11)
    plt.ylim(0.031, 0.0325)
    png_path_4 = os.path.join(out_path, f"pca_zoomed_4.png")
    plt.savefig(png_path_4, dpi=600)
    '''
    # 262-sample zoom
    plt.xlim(-0.0025, 0.0025)
    plt.ylim(-0.01, 0.04)
    png_path_1 = os.path.join(out_path, f"pca_zoomed_1.png")
    plt.savefig(png_path_1, dpi=600) 
    
    print("Plots saved to ", png_path)
    
    
    # Save point coordinates
    csv_path = os.path.join(out_path, "points.csv")
    # currently - only save PCs 1 and 2 since that's what we're plotting - can change
    points.to_csv(csv_path, columns=["species", "pop", "Du#", "LP#", "PC1", "PC2"])
    print("Points saved to ", csv_path)


    # --- Extra Stuff ---

    '''
    # Obtain p-values from PC-based selection scan, if provided
    if selection_path:
        D = np.loadtxt(selection_path)    # Reads PC based selection statistics
        p = chi2.sf(D, 1)
        print("p-values from PC-based selection scan:\n", p)
    '''


    # Save various info & metrics to a file
    metrics_path = os.path.join(out_path, "metrics.txt")

    with open(metrics_path, 'w') as f:
        # f.write("——————\nScript Info\n——————\n")
        
        f.write("——————\nPopulation Info\n——————\n")
        f.write("_species_\t_population_\t_sample_\n")
        for species, pop, Du, LP in zip(species_labels, pop_labels, Du_labels, LP_labels):
            f.write(f"{species}\t{pop}\t\t{Du}\t{LP}\n")

        f.write("——————\nPCA Summary\n——————\n")
        for i, val in enumerate(evals):
            f.write(f"PC{i+1} explains {100*evals[i]/evals_sum}% of variance\n")

        # this will get really big with many samples. make optional?
        f.write("——————\nPCA Details\n——————\n")
        f.write("Eigenvalues:\t\tEigenvectors:\n")
        for val, vecs in zip(evals, evecs):
            f.write(f"{val}\t{vecs}\n")

        
    print("Metrics file saved to ", metrics_path)
    
    return None


# --- Snakemake Arg Handling ---

# snakemake object is defined and passed in by the Snakemake rule
cov_path = snakemake.input[0]
info_path = snakemake.input[1]
out_path = os.path.dirname(snakemake.output[0])     # output is .../pca.png; need the parent dir

# for now, just plot PCs 1 and 2
x_pc = 1
y_pc = 2

plot_pca(cov_path, info_path, out_path, x_pc=x_pc, y_pc=y_pc)
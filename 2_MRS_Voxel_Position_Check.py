import os
import numpy as np
import pandas as pd
import nibabel as ni
import matplotlib.pyplot as plt
from nilearn import plotting
from matplotlib import ticker
import matplotlib.gridspec as gridspec

# --- Project directory setup ---
# ------------------------------------------ !!! Choose correct folder and amend other relevant lines !!! ------------------------------------------
data_dir = '/data/project/CannTeen/data/Adults/'  # One level up from the subject folders

# Directory to save figures in
fig_dir = os.path.join(data_dir, 'figures_Adults')  # e.g., figures_ACC, figures_Tha, etc.
if not os.path.isdir(fig_dir):
    os.mkdir(fig_dir)

# --- Load participant info from TSV ---
tsv_path = os.path.join(data_dir, 'Adults_Participant_filtered.tsv')
df = pd.read_csv(tsv_path, sep='\t')

# Ensure participant IDs are strings and zero-padded to 6 digits
df['participant'] = df['participant'].astype(str).str.zfill(6)

# Construct full session label (e.g., 016001_V1)
df['participant_session'] = df['participant'] + '_' + df['session'].astype(str)

# Define group names
group_name_1 = 'Placebo'
group_name_2 = 'THC'
group_name_3 = 'THC+CBD'

# Filter participant_sessions by group
group1_sessions = df[df['group'] == group_name_1]['participant_session'].tolist()
group2_sessions = df[df['group'] == group_name_2]['participant_session'].tolist()
group3_sessions = df[df['group'] == group_name_3]['participant_session'].tolist()

# Define voxel mask filename
mask_name = 'bw_metab_mask.nii.gz'  # change for Thalamus, ACC, etc.

# --- Load one example to get affine and dimensions ---
example_sub, example_ses = group1_sessions[0].split('_')
example_path = os.path.join(data_dir, example_sub, example_ses, mask_name)
tmp_img = ni.load(example_path)
mask_aff = tmp_img.affine
mask_dims = tmp_img.shape

# --- Preallocate arrays ---
mask_data_group_1 = np.zeros((len(group1_sessions), *mask_dims))
mask_data_group_2 = np.zeros((len(group2_sessions), *mask_dims))
mask_data_group_3 = np.zeros((len(group3_sessions), *mask_dims))

# --- Load mask data for each group ---
def load_group_data(group_sessions, output_array):
    for i, ses_id in enumerate(group_sessions):
        sub, ses = ses_id.split('_')
        fpath = os.path.join(data_dir, sub, ses, mask_name)
        if os.path.isfile(fpath):
            output_array[i, :, :, :] = ni.load(fpath).get_fdata()
        else:
            print(f"Missing file: {fpath}")
    return output_array

mask_data_group_1 = load_group_data(group1_sessions, mask_data_group_1)
mask_data_group_2 = load_group_data(group2_sessions, mask_data_group_2)
mask_data_group_3 = load_group_data(group3_sessions, mask_data_group_3)

# --- Calculate overlap percentage ---
def compute_density(mask_data):
    return (np.sum(mask_data, axis=0) / mask_data.shape[0]) * 100

density_1 = compute_density(mask_data_group_1)
density_2 = compute_density(mask_data_group_2)
density_3 = compute_density(mask_data_group_3)

# Create NIfTI images for density maps
density_map_1 = ni.Nifti1Image(density_1, mask_aff)
density_map_2 = ni.Nifti1Image(density_2, mask_aff)
density_map_3 = ni.Nifti1Image(density_3, mask_aff)

# --- Plotting ---
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

fig = plt.figure()
fig.set_size_inches(cm2inch(28, 14))

gs = gridspec.GridSpec(1, 4, width_ratios=[8, 1, 1, 1])
ax1 = plt.subplot(gs[0])

display = plotting.plot_glass_brain(None, threshold=0, colorbar=False, axes=ax1,
                                    display_mode='xz', alpha=0.5, annotate=False)
display.add_contours(density_map_1, cmap='Reds', alpha=0.7)
display.add_contours(density_map_2, cmap='Blues', alpha=0.7)
display.add_contours(density_map_3, cmap='Greens', alpha=0.7)

# Colorbars
for idx, cmap in enumerate(["Reds", "Blues", "Greens"], start=1):
    ax = plt.subplot(gs[idx])
    a = np.array([[0, 1]])
    img = plt.imshow(a, cmap=cmap, axes=ax)
    ax.set_visible(False)
    cbar = plt.colorbar(ax=ax)
    cbar.set_label('Percentage')

# Adjust y-axis to percentage
for ax in plt.gcf().axes:
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(1))

# Save figure
fig_file = os.path.join(fig_dir, 'voxel-density-map_multiple-groups.png')
fig.savefig(fig_file, bbox_inches='tight', dpi=1200)

# Save NIFTI maps
density_map_1.to_filename(os.path.join(fig_dir, 'voxel_density_map_group_1.nii.gz'))
density_map_2.to_filename(os.path.join(fig_dir, 'voxel_density_map_group_2.nii.gz'))
density_map_3.to_filename(os.path.join(fig_dir, 'voxel_density_map_group_3.nii.gz'))

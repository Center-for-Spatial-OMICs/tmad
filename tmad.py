import pandas as pd
import numpy as np
import os
import scipy.io
import gzip
import sys
import click

# load the files
# user inputted paths
'''
XR_path = '/mnt/scratch1/Luke/XeniumTMASeparate/JP3084_PANC/'
coord_path = '/mnt/scratch1/Luke/XeniumTMASeparate/TMACoords_JP3084_PANC/'
target_path = 'TOUCHSTONE_DIVORCE_XR_FFPE_PA_WCM_N_R1/'
'''

VERSION = "0.1.0"

@click.command()
@click.version_option(VERSION, message="%(version)s")
@click.argument("paths", nargs=3)

def main(paths):
    XR_path = paths[0]
    coord_path = paths[1]
    target_path = paths[2]

    tx_path = os.path.join(XR_path,'transcripts.csv.gz')
    exp_mat = os.path.join(XR_path,'cell_feature_matrix')
    meta_path = os.path.join(XR_path,'cells.csv.gz')

    tx = pd.read_csv(tx_path, compression='gzip')
    tx['TMA'] = 0
    meta = pd.read_csv(meta_path, compression='gzip')
    meta['TMA'] = 0

    exp = scipy.io.mmread(f'{exp_mat}/matrix.mtx.gz').toarray()
    cols = pd.read_csv(f'{exp_mat}/barcodes.tsv.gz', header=None, sep='\t')
    col_names =  cols.iloc[:, 0].values  # This is the barcodes of cells

    df_coord = load_coords(coord_path=coord_path)
    meta = sep_meta(df_coord=df_coord, coord_files=os.listdir(coord_path), meta=meta)
    tx = sep_tx(df_coord=df_coord, coord_files=os.listdir(coord_path), tx=tx)
    TMA_indices = gen_filter_idx(meta=meta, cols=cols)
    save_files(target_path=target_path, coord_files=os.listdir(coord_path), TMA_indices=TMA_indices, 
               cols=cols, exp=exp, exp_mat=exp_mat, meta=meta, tx=tx)


if __name__ == '__main__':
    main()
    
# load the coordinates
def load_coords(coord_path):
    coord_files = os.listdir(coord_path)
    df_coord = pd.DataFrame()
    for i, file in enumerate(coord_files):
        temp = pd.read_csv(os.path.join(coord_path, file), skiprows=2, skipfooter=1, engine='python')
        df_coord[f'X{i+1}'] = temp['X']
        df_coord[f'Y{i+1}'] = temp['Y']
    # print(df_coord)
    return df_coord

# define the function
def incore(xmax, xmin, ymax, ymin, x, y):
    if (x <= xmax) and (x >= xmin) and (y <= ymax) and (y >= ymin):
        return True
    else:
        return False

# go through the metadata and find the TMA core for each cell
# using this, we can then make a dictionary of each cell's TMA
def sep_meta(df_coord, coord_files, meta):
    for i, row in meta.iterrows():
        for j in range(len(coord_files)):
            xmax = max(df_coord['X'+str(j+1)].values)
            xmin = min(df_coord['X'+str(j+1)].values)
            ymax = max(df_coord['Y'+str(j+1)].values)
            ymin = min(df_coord['Y'+str(j+1)].values)
            if incore(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, x=row['x_centroid'], y=row['y_centroid']):
                meta.at[i, 'TMA'] = j+1
    return meta

# separate the tx file
def sep_tx(df_coord, coord_files, tx):
    for i, row in tx.iterrows():
        for j in range(len(coord_files)):
            xmax = max(df_coord['X'+str(j+1)].values)
            xmin = min(df_coord['X'+str(j+1)].values)
            ymax = max(df_coord['Y'+str(j+1)].values)
            ymin = min(df_coord['Y'+str(j+1)].values)
            if incore(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, x=row['x_location'], y=row['y_location']):
                tx.at[i, 'TMA'] = j+1
    # print(tx.head())
    return tx

# create an index to filter the expmat
def gen_filter_idx(meta, cols):
    cell_id_TMA = meta[['cell_id', 'TMA']]
    cell_id_TMA = cell_id_TMA[cell_id_TMA['TMA'] != 0]
    cell_id_TMA = cell_id_TMA.set_index('cell_id')

    TMA_indices = {}

    for index, row in cols.iterrows():
        cell_id = row[0]
        if cell_id in cell_id_TMA.index:
            tma = cell_id_TMA.loc[cell_id, 'TMA']
            if tma not in TMA_indices:
                TMA_indices[tma] = []
            TMA_indices[tma].append(index)
    
    return TMA_indices

# save the files
# first create the home directory
def save_files(target_path, coord_files, TMA_indices, cols, exp, exp_mat, meta, tx):
    if not os.path.exists(target_path):
        os.makedirs(target_path)
    for i in range(len(coord_files)):
        # create a separate folder for each TMA core
        curr_dir_name = target_path.split('_')
        curr_dir_name.insert(5, str(i+1))
        curr_dir_name = '_'.join(curr_dir_name)

        if not os.path.exists(target_path+curr_dir_name):
            os.makedirs(target_path+curr_dir_name)

        os.makedirs(target_path+curr_dir_name+'cell_feature_matrix/')

        # save the barcodes, mtx, and copy the features
        temp_barcode = cols.loc[TMA_indices[i+1]]
        temp_barcode.to_csv(os.path.join(target_path,curr_dir_name,'cell_feature_matrix/barcodes.tsv.gz'), sep='\t', index=False, header=False, compression='gzip')

        temp_mat = exp[:, TMA_indices[i+1]]
        # np.savetxt(target_path+curr_dir_name+'cell_feature_matrix/matrix.mtx', temp_mat, delimiter='\t', fmt='%.18e', comments='', header='%%MatrixMarket matrix coordinate real general')
        scipy.io.mmwrite(os.path.join(target_path,curr_dir_name,'cell_feature_matrix/matrix.mtx'), temp_mat)

        # Compress the MTX file using gzip
        with gzip.open(os.path.join(target_path,curr_dir_name,'cell_feature_matrix/matrix.mtx.gz'), 'wb') as f:
            f.write(open(os.path.join(target_path,curr_dir_name,'cell_feature_matrix/matrix.mtx'), 'rb').read())

        # Remove the uncompressed MTX file
        os.remove(os.join.path(target_path,curr_dir_name,'cell_feature_matrix/matrix.mtx'))

        # copy the features
        '''if os.name == 'nt':  # Windows
            os.system(f'copy "{f'{exp_mat}/features.tsv.gz'}" "{os.path.join(target_path,curr_dir_name,'cell_feature_matrix/', os.path.basename(f'{exp_mat}/features.tsv.gz'))}"')
        else:  # Unix-based (Linux, macOS)
            os.system(f'cp "{f'{exp_mat}/features.tsv.gz'}" "{os.path.join(target_path,curr_dir_name,'cell_feature_matrix/', os.path.basename(f'{exp_mat}/features.tsv.gz'))}"')'''

        # Construct the source and destination paths
        source_path = os.path.join(exp_mat, 'features.tsv.gz')
        destination_dir = os.path.join(target_path, curr_dir_name, 'cell_feature_matrix')
        destination_path = os.path.join(destination_dir, os.path.basename(source_path))

        # Ensure the destination directory exists
        os.makedirs(destination_dir, exist_ok=True)

        # Determine the appropriate copy command based on the OS
        if os.name == 'nt':  # Windows
            os.system(f'copy "{source_path}" "{destination_path}"')
        else:  # Unix-based (Linux, macOS)
            os.system(f'cp "{source_path}" "{destination_path}"')
            
        # save the metadata
        temp_meta = meta[meta['TMA'] == (i+1)]
        temp_meta = temp_meta.drop(['TMA'], axis=1)
        temp_meta.to_csv(os.join.path(target_path,curr_dir_name,'cells.csv.gz'), index=False, compression='gzip')

        # save the tx
        temp_tx = tx[tx['TMA'] == (i+1)]
        temp_tx = temp_tx.drop(['TMA'], axis=1)
        temp_tx.to_csv(os.join.path(target_path,curr_dir_name,'transcripts.csv.gz'), index=False, compression='gzip')
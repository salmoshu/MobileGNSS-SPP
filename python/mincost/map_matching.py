import glob
import numpy as np
import pandas as pd

import transform

INPUT_PATH = '../input/google-smartphone-decimeter-challenge'

COLLECTION_LIST_ALL = np.array(sorted(path.split('/')[-1] for path in glob.glob(f'{INPUT_PATH}/train/*')))
COLLECTION_LIST_DOWNTOWN = COLLECTION_LIST_ALL[np.array([24,27,29]) - 1]

def get_database_vectors():
    gt_df_list = []
    for collection in COLLECTION_LIST_DOWNTOWN:
        phone_list = [path.split('/')[-1] for path in glob.glob(f'{INPUT_PATH}/train/{collection}/*')]
        for phone in phone_list:
            gt_df = pd.read_csv(f'{INPUT_PATH}/train/{collection}/{phone}/ground_truth.csv')
            gt_df_list.append(gt_df)
    gt_df  = pd.concat(gt_df_list, axis=0)
    gt_df  = gt_df[['latDeg', 'lngDeg']].drop_duplicates(ignore_index=True)
    gt_blh = transform.BLH(
        lat=np.deg2rad(gt_df['latDeg']),
        lng=np.deg2rad(gt_df['lngDeg']),
        hgt=np.zeros(gt_df.shape[0]),
    )
    P  = transform.BLH_to_ECEF(gt_blh).to_numpy() # shape = (3, N_P)
    PP = np.sum(P**2, axis=0) # shape = (N_P, )
    return gt_df, P, PP

gt_df, P, PP = get_database_vectors()

def snap_to_nearest_neighbor(pred_df):
    pred_blh = transform.BLH(
        lat=np.deg2rad(pred_df['latDeg']),
        lng=np.deg2rad(pred_df['lngDeg']),
        hgt=np.zeros(pred_df.shape[0]),
    )
    Q = transform.BLH_to_ECEF(pred_blh).to_numpy() # shape=(3, N_Q)

    QQ  = np.sum(Q**2, axis=0) # shape=(N_Q, )
    PQ  = P.T @ Q              # shape=(N_P, N_Q)
    d2  = (QQ - 2 * PQ).T + PP # shape=(N_Q, N_P)
    idx = np.argmin(d2, axis=1)

    nn_df = gt_df.iloc[idx, :].reset_index(drop=True)
    pp_df = pred_df.copy()
    pp_df['latDeg'] = nn_df['latDeg']
    pp_df['lngDeg'] = nn_df['lngDeg']
    return pp_df

def distance_to_nearest_neighbor(pred_df):
    pred_blh = transform.BLH(
        lat=np.deg2rad(pred_df['latDeg']),
        lng=np.deg2rad(pred_df['lngDeg']),
        hgt=np.zeros(pred_df.shape[0]),
    )
    Q = transform.BLH_to_ECEF(pred_blh).to_numpy() # shape=(3, N_Q)
    N_Q = Q.shape[1]

    QQ  = np.sum(Q**2, axis=0) # shape=(N_Q, )
    PQ  = P.T @ Q              # shape=(N_P, N_Q)
    d2  = (QQ - 2 * PQ).T + PP # shape=(N_Q, N_P)
    idx = np.argmin(d2, axis=1)
    d   = np.zeros(N_Q)
    for i, j in zip(range(N_Q), idx):
        d[i] = np.sqrt(np.maximum(0, d2[i, j]))
    return d
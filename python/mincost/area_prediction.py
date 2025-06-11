import numpy as np
import pandas as pd
from pathlib import Path
from glob import glob
from sklearn.neighbors import KNeighborsClassifier

BASE_DIR = Path('../input/google-smartphone-decimeter-challenge')

train_base = pd.read_csv(BASE_DIR / 'baseline_locations_train.csv')
train_base = train_base.sort_values([
    "collectionName", "phoneName", "millisSinceGpsEpoch"
]).reset_index(drop=True)

train_base['area'] = train_base['collectionName'].map(lambda x: x.split('-')[4])

train_name = np.array(sorted(path.split('/')[-1] for path in glob(f'{BASE_DIR}/train/*')))
train_highway  = train_name[np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]) - 1]
train_tree     = train_name[np.array([22,23,25,26,28]) - 1]
train_downtown = train_name[np.array([24,27,29]) - 1]

train_base['area_target'] = -1
train_base.loc[train_base['collectionName'].isin(train_highway),  'area_target'] = 0
train_base.loc[train_base['collectionName'].isin(train_tree),     'area_target'] = 1
train_base.loc[train_base['collectionName'].isin(train_downtown), 'area_target'] = 2

def processing_downtown(input_df: pd.DataFrame, is_train=False):
    output_df = input_df.groupby('collectionName')[['latDeg', 'lngDeg']].std()
    if is_train:
        output_df = output_df.merge(
            input_df.groupby('collectionName')[['area_target']].first(),
            on='collectionName')
    output_df = output_df.merge(
        input_df.groupby('collectionName')['area'].first(),
        on='collectionName')
    output_df = output_df.merge(
        input_df.groupby('collectionName')['phoneName'].unique().apply(list),
        on='collectionName')
    return output_df

train = processing_downtown(train_base, is_train=True)
train['downtown_target'] = (train['area_target']==2).astype(int)

downtown_model_knn = KNeighborsClassifier(n_neighbors=1)
downtown_model_knn.fit(
    train[['latDeg', 'lngDeg']],
    train['downtown_target'],
)

def processing_highway_tree(input_df: pd.DataFrame, is_train=False):
    output_df = input_df.groupby('collectionName')[['latDeg', 'lngDeg']].min()
    if is_train:
        output_df = output_df.merge(
            input_df.groupby('collectionName')[['area_target']].first(),
            on='collectionName')
    output_df = output_df.merge(
        input_df.groupby('collectionName')['area'].first(),
        on='collectionName')
    output_df = output_df.merge(
        input_df.groupby('collectionName')['phoneName'].unique().apply(list),
        on='collectionName')
    return output_df

train = processing_highway_tree(train_base, is_train=True)

highway_tree_model_knn = KNeighborsClassifier(n_neighbors=1)
highway_tree_model_knn.fit(
    train.loc[train['area_target']!=2, ['latDeg', 'lngDeg']],
    train.loc[train['area_target']!=2, 'area_target'],
)

def predict_area(test_base):
    test_base = test_base.copy()
    test_base = test_base.sort_values([
        "collectionName", "phoneName", "millisSinceGpsEpoch"
    ]).reset_index(drop=True)
    test_base['area'] = test_base['collectionName'].map(lambda x: x.split('-')[4])

    test = processing_downtown(test_base)
    downtown_pred = downtown_model_knn.predict(test[['latDeg', 'lngDeg']])

    test = processing_highway_tree(test_base)
    test.loc[downtown_pred==1, 'area_pred'] = 2
    pred = highway_tree_model_knn.predict(test.loc[test['area_pred'].isnull(), ['latDeg', 'lngDeg']])
    test.loc[test['area_pred'].isnull(), 'area_pred'] = pred
    test['area_pred'] = test['area_pred'].astype(int)
    test['collectionName'] = test.index

    test_highway  = []
    test_tree     = []
    test_downtown = []
    for collection, area_pred in test[['collectionName', 'area_pred']].itertuples(index=False):
        if area_pred == 0:
            test_highway.append(collection)
        elif area_pred == 1:
            test_tree.append(collection)
        else:
            test_downtown.append(collection)
    return (test_highway, test_tree, test_downtown)
import pandas as pd
BASELINE_DF = pd.concat([pd.read_csv('../input/gsdc-improved-raw-gnss-baseline-result/raw_gnss_train.csv'),
                         pd.read_csv('../input/gsdc-improved-raw-gnss-baseline-result/raw_gnss_test.csv'),
                        ], axis=0)
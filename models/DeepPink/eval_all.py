import time
import math
import sys
import os

import numpy as np
import pandas as pd
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.linalg import qr

from sklearn.preprocessing import StandardScaler,LabelEncoder

num_epochs = 200;
batch_size = 10;
filterNum = 1;
bias = True;
activation='relu';

q=0.05;

dataDir = '/media/yanglu/TOSHIBA/data/featuresAndResponseDataframeCSVs/2018_6_15_fiveFoldCSVs'
dataTypeList = ['LRnoMotifs', 'LSnoMotifs', 'LTnoMotifs'];


for dataType in dataTypeList:
    dataURL = os.path.join(dataDir, dataType, dataType+'.csv'); print(dataURL);

    columnHeader = pd.read_csv(dataURL, nrows=1,usecols=lambda x:x != 'Response').columns;
    columnHeader = columnHeader[2:];
    print('columnHeader size: '+str(columnHeader.shape))

    dataSubDir = os.path.join(dataDir, dataType, 'result_2layer_epoch' + str(num_epochs) + '_batch' + str(batch_size) + '_knockoff_all');

    # featImportURL = os.path.join(dataSubDir, 'result_epoch' + str(num_epochs) + '_featImport.csv');
    # featImportMat = pd.read_csv(featImportURL, header=None).values.astype(float);

    featWeightURL = os.path.join(dataSubDir, 'result_epoch' + str(num_epochs) + '_featWeight.csv');
    featWeightMat = pd.read_csv(featWeightURL, header=None).values.astype(float);
    print('featWeightMat: ' + str(featWeightMat.shape));

    pVal = featWeightMat.shape[1] / 2;


    # allSelectedIndices = [];
    # allW = np.empty((0, pVal), float);
    # for recordNum in range(featWeightMat.shape[0]):
    #     resultVec = featWeightMat[recordNum,:];
    #
    #     W = np.fabs(resultVec[0:pVal])-np.fabs(resultVec[pVal:]); print(W);
    #     t = np.concatenate(([0],np.sort(np.fabs(W)))); print(t);
    #     allW = np.vstack((allW, W.reshape((1, pVal))));
    #
    #
    #     ratio = np.zeros(pVal);
    #     for j in range(pVal):
    #         ratio[j] = 1.0*(len(np.where(W <= -t[j])[0])) / np.max((1, len(np.where(W >= t[j])[0]) ));
    #     print(ratio);
    #
    #     T = np.inf;
    #     arr = np.where(ratio <= q)[0];
    #     if len(arr) > 0:
    #         id = np.min(arr);
    #         T = t[id];
    #
    #     S = np.where(np.fabs(W) >= T)[0];
    #     allSelectedIndices = np.concatenate((allSelectedIndices, S));
    # print(np.unique(allSelectedIndices, return_counts=True));


    resultVec = np.sum(featWeightMat, axis=0);
    W = np.fabs(resultVec[0:pVal]) - np.fabs(resultVec[pVal:]);
    t = np.concatenate(([0],np.sort(np.fabs(W))));

    ratio = np.zeros(pVal);
    for j in range(pVal):
        ratio[j] = 1.0*(len(np.where(W <= -t[j])[0])) / np.max((1, len(np.where(W >= t[j])[0]) ));
    print(ratio);

    T = np.inf;
    arr = np.where(ratio <= q)[0];
    if len(arr) > 0:
        id = np.min(arr);
        T = t[id];

    qualifiedIndices = np.where(np.fabs(W) >= T)[0];

    meanW = W;
    meanWeightW = np.mean(featWeightMat, axis=0);

    # result = np.unique(allSelectedIndices, return_counts=True);
    # uniqueIndices = result[0];
    # counts = result[1];
    #
    # qualifiedIndices = np.asarray(uniqueIndices[np.where(counts >= 0.5*featWeightMat.shape[0])[0]], int);
    # print(columnHeader[qualifiedIndices])
    #
    # meanW = np.mean(allW, axis=0);
    # meanWeightW = np.mean(featWeightMat, axis=0);
    #
    with open(os.path.join(dataSubDir, 'featuresSelected_fdr'+str(q)+'_knockoff.txt'), 'w') as featureOutput:
        for feat in columnHeader[qualifiedIndices]:
            featureOutput.write(feat+'\n');

    with open(os.path.join(dataSubDir, 'featuresAll_knockoff_stats.txt'),'w') as statOutput:
        for idx in range(len(meanW)):
            statOutput.write(columnHeader[idx]+'\t'+str(meanW[idx])+'\t'+str(meanWeightW[idx])+'\n');



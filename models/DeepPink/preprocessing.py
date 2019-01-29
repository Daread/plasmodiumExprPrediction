import time
import math
import sys
import os

import numpy as np
import pandas as pd
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.linalg import qr

import keras
from keras.models import Sequential, Model, model_from_json
from keras.layers import Input, Dense, Dropout, BatchNormalization, merge, LocallyConnected1D, Flatten, Conv1D
from keras import backend as K
from keras import regularizers
from keras.objectives import mse
from keras import regularizers, optimizers
from keras.callbacks import EarlyStopping
from keras.initializers import Constant

from sklearn.preprocessing import StandardScaler,LabelEncoder
from sklearn.metrics import auc, roc_curve, roc_auc_score
from sklearn.model_selection import GridSearchCV


dataDir = '/media/yanglu/TOSHIBA/data/featuresAndResponseDataframeCSVs/2018_6_15_fiveFoldCSVs'
dataTypeList = ['LRnoMotifs', 'LRwithMotifs', 'LSnoMotifs', 'LSwithMotifs', 'LTnoMotifs', 'LTwithMotifs'];
testIndicesList = [range(0,685), range(685,1396)];



for dataType in dataTypeList:
    dataURL = os.path.join(dataDir, dataType, dataType+'.csv'); print(dataURL);

    dataMat = pd.read_csv(dataURL, usecols=lambda x: x != 'Response').values.astype(float);
    dataMat = dataMat[:, 2:];
    pVal = dataMat.shape[1];

    labelVec = pd.read_csv(dataURL, usecols=lambda x: x == 'Response');
    labelVec = labelVec['Response'].map({'Low': 0, 'High': 1}).values;

    np.savetxt(os.path.join(dataDir, dataType, 'X.csv'), dataMat, fmt='%.6f', delimiter=",")
    np.savetxt(os.path.join(dataDir, dataType, 'Y.csv'), labelVec, fmt='%d', delimiter=",")


# Python 3.7

"""
File:
Author:         Shengyuan Wang
Date:           Mar 16, 2020
Last update:    Jun 6, 2020

Purpose:        Assign secondary structure by tessellation using Keras.

Input files:    Tessellation results.
Output files:   Keras model to assign protein secondary structure.
"""


import os
import pandas as pd
import numpy as np
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.utils import resample

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
import eli5
from eli5.sklearn import PermutationImportance

from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasClassifier
from keras.utils import np_utils


def data_preparation(input_dir, filename):
    """Data preprocessiong."""
    data = pd.read_csv(os.path.join(input_dir, filename))
    data = data[['Ave', 'Tetrahedrality', 'Volume', 'Rotation', 'Bend', 'Twist', 'Class', 'Residue distance',
                 'Simplex_aa', 'Secondary structure']]

    data_temp = []
    sec_classes = ['BBBB', 'BBBC', 'BBBH', 'BBCC', 'BBCH', 'BBHH', 'BCCC', 'BCCH', 'BCHH', 'BHHH', 'CCCC',
                   'CCCH', 'CCHH', 'CHHH', 'HHHH']

    for sec_class in sec_classes:
        data_temp.append(data[data['Secondary structure'] == sec_class])

    downsize = 1000000
    for sec_data in data_temp:
        sec_data_len = len(sec_data)
        if sec_data_len < downsize:
            downsize = sec_data_len

    data_temp_downsize = []
    for sec_data in data_temp:
        data_temp_downsize.append(resample(sec_data, replace=False, n_samples=downsize, random_state=0))

    data = pd.concat(data_temp_downsize)
    data = data.reset_index()

    data1 = data[['Ave', 'Tetrahedrality', 'Volume', 'Rotation', 'Bend', 'Twist', 'Class', 'Residue distance']]
    data2 = data[['Simplex_aa']]
    data3 = data[['Secondary structure']]

    le = preprocessing.LabelEncoder()
    data2 = data2.apply(le.fit_transform)

    x = np.concatenate((data1, data2), axis=1)
    x = preprocessing.scale(x)

    encoder = preprocessing.LabelEncoder()
    encoder.fit(data3)
    encoded_Y = encoder.transform(data3)

    # convert integers to dummy variables (i.e. one hot encoded)
    dummy_y = np_utils.to_categorical(encoded_Y)

    # One hot encoded multi-classes
    X_train, X_test, y_train, y_test = train_test_split(x, dummy_y, test_size=0.2, random_state=0)

    # Single digit multi-classes
    X_train_sd, X_test_sd, y_train_sd, y_test_sd = train_test_split(x, encoded_Y, test_size=0.2, random_state=0)

    one_hot_data = [X_train, X_test, y_train, y_test]
    single_digit_data = [X_train_sd, X_test_sd, y_train_sd, y_test_sd]

    return x, dummy_y, one_hot_data, single_digit_data


def baseline_model():
    # create model
    model = Sequential()
    model.add(Dense(8, input_dim=9, activation='relu'))             # input_dim is the number of features
    model.add(Dense(15, activation='softmax'))                      # input is the number of classes

    # Compile model
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    return model


if __name__ == '__main__':
    # file directory
    input_dir = 'INDIVIDUAL/sec_str'
    filename = 'edges_t_v angles_all.csv'

    # Cleaning data
    x, dummy_y, one_hot_data, single_digit_data = data_preparation(input_dir, filename)
    X_train, X_test, y_train, y_test = one_hot_data[0], one_hot_data[1], one_hot_data[2], one_hot_data[3]
    X_train_sd, X_test_sd, y_train_sd, y_test_sd = single_digit_data[0], single_digit_data[1], single_digit_data[2], \
                                                   single_digit_data[3]

    # Performing Keras model
    estimator = KerasClassifier(build_fn=baseline_model, epochs=200, batch_size=5, verbose=0)

    # Estimate the performance
    kfold = KFold(n_splits=10, shuffle=True)
    results = cross_val_score(estimator, x, dummy_y, cv=kfold)
    print("Baseline: %.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))

    estimator.fit(X_train, y_train)
    y_pred = estimator.predict(X_test)

    print(confusion_matrix(y_test_sd, y_pred))
    print(classification_report(y_test_sd, y_pred))

    perm = PermutationImportance(estimator, random_state=1).fit(X_train, y_train)
    print(eli5.format_as_text(eli5.explain_weights(perm)))
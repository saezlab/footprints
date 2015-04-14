#!/usr/bin/env python2

# plots the number of SPEED platforms vs pathways as color matrix with added values

import numpy as np
import pandas as pd
import pandas.io.parsers as parsers
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def plotMatrix(indexDF, plotX, plotY, plotZ=None):
    groupsX = set(indexDF[plotX])
    df = pd.DataFrame(index=set(indexDF[plotY]))

    for groupX in groupsX:
        mysub = indexDF[indexDF[plotX]==groupX]
        groupY = mysub.groupby(plotY).groups
        if not plotZ:
            cover = np.array([[k,len(v)] for k,v in groupY.iteritems()])
        else:
            cover = np.array([[k,len(set(mysub.loc[v,plotZ]))] for k,v in groupY.iteritems()])
        df[groupX] = pd.Series(np.array(cover[:,1], dtype=int), index=cover[:,0])

    cols = df.sum(axis=0)
    rows = df.sum(axis=1)

    cols.sort()
    df = df[cols.index]
    collabel = [x+" ("+str(int(y))+")" for x,y in zip(df.columns,cols)]
    rowlabel = [x+" ("+str(int(y))+")" for x,y in zip(df.index,rows)]

    plt.imshow(df, interpolation='nearest', cmap='Blues', norm=colors.LogNorm())
    plt.colorbar()
    plt.yticks(range(len(df.index)), rowlabel, size='small')
    plt.xticks(np.arange(len(df.columns))+.5, collabel, size='small', rotation=45, ha='right')

    for y,idx in enumerate(df.index):
        for x,column in enumerate(df.columns):
            if not np.isnan(df.loc[idx,column]):
                plt.annotate(int(df.loc[idx,column]), xy=(x,y), fontsize=8, ha="center", va="center")


if __name__ == '__main__':
    indexDF = parsers.read_csv('index.oldspeed.txt', sep="\t")
    #validPlatform = ['GPL570', 'GPL571', 'GPL96']

    plotMatrix(indexDF, plotX='GPL', plotY='PATHWAY')#, plotZ='cells')
#    plotMatrix(indexDF, plotX='GPL', plotY='pathway')#, plotZ='cells')
    #plotZ e.g. how many different cell types
    plt.show()


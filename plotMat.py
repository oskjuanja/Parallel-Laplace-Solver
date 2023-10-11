import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys


def loadMat(path = "./filedata.out"):
    mat = np.loadtxt(path, dtype=float)
    return mat


def plotMat(mat):
    # plt.imshow(mat, cmap='hot', interpolation='nearest')
    ax = sns.heatmap(mat, cmap='hot', xticklabels=False, yticklabels=False)
    plt.title("Grid Heat Distribution")
    plt.show()


def main():
    if (len(sys.argv) > 1):
        path = sys.argv[1]
        mat = loadMat(path)
    else:
        mat = loadMat()
    
    plotMat(mat)

if __name__ == '__main__':
    main()
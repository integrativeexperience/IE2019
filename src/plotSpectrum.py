# basic plotting code for spectra
import matplotlib.pyplot as plt
import numpy as np

def openFile(f, headerlines = 17):
    with open(f) as fin:
        datatext = fin.readlines()
    data = []
    for line in datatext[headerlines:]:
        if bool(line) == True:
            d = line.split()
            d = [float(dpt) for dpt in d]
        data.append(d)
    data = np.array(data)
    return np.transpose(data)

def plotData(data):
    plt.plot(data[0], data[1])

if __name__ == "__main__":
    data = openFile("spectrum") # change this to the filename you want to plot
    plotData(data)
    plt.show()
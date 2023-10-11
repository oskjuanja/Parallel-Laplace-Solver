import matplotlib.pyplot as plt
import numpy as np

def plot(xdata, ydata, label):

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle('Execution Time w.r.t. Thread Count', fontsize=18)

    ax.scatter(xdata, ydata, label=label)
    ax.set_xlabel("# of threads", fontsize=16)
    ax.set_ylabel("Speedup", fontsize=16)
    ax.legend(fontsize=16)

    plt.show()


def main():
    times = np.loadtxt("./exec_time.data", dtype=float)

    label = "Size: {}".format(int(times[0]))
    times = [float(x) for x in times[1:]]
    threads = [i+1 for i in range(len(times))]

    plot(threads, times, label)

if __name__ == '__main__':
    main()
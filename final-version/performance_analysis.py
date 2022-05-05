import numpy as np
import matplotlib.pyplot as plt


def generate_strong_scaling_plot(n, time, analysis_num):
    """
    Generates a plot of the strong scaling results.
    """
    t1 = time[0]
    Sn = [t1/t for t in time]

    n = np.array(n)
    Sn = np.array(Sn)

    # Plot the data
    plt.plot(n, Sn, '-o')
    # Set the labels
    plt.xlabel("n (number of MPI ranks)")
    plt.ylabel("S(n)")
    # Set the title
    plt.title("Strong Scaling Analysis")
    # Save the plot
    plt.savefig("./graphs/strong_scaling_{}.png".format(analysis_num))
    plt.clf()

def generate_weak_scaling_plot(n, grind_rates, analysis_num):
    """
    Generates a plot of the weak scaling results.
    """
    g1 = grind_rates[0]
    Sn = [g/g1 for g in grind_rates]

    n = np.array(n)
    Sn = np.array(Sn)

    # Plot the data
    plt.plot(n, Sn, '-o')
    # Set the labels
    plt.xlabel("n (number of MPI ranks)")
    plt.ylabel("S(n)")
    # Set the title
    plt.title("Weak Scaling Analysis")
    # Save the plot
    plt.savefig("./graphs/weak_scaling_{}.png".format(analysis_num))
    plt.clf()

n1 = [1,2,4,8,16,32,64]
time1 = [595781.818, 284939.13, 136533.333, 78959.036, 50412.308, 30340.741, 15603.81]
generate_strong_scaling_plot(n1, time1, "1")

n2 = n1
grind_rates = [80,79,78,67,57,48,39]
generate_weak_scaling_plot(n2, grind_rates, "1")

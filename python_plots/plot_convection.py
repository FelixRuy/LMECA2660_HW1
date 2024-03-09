import numpy as np
import matplotlib.pyplot as plt

time1 = float(input("Time for measure 1 ? (use x.y as format): "))
time2 = float(input("Time for measure 2 ? (use x.y as format): "))


def read_vector_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    y = []
    x = []
    t = []
    for line in lines:
        parts = line.split(":")
        time = float(parts[0].split("=")[1].strip())
        values = [float(val.strip()) for val in parts[1].split(",")[:-1]]
        if time < 0:
            x = np.array(values)
        else:
            y.append(values)
            t.append(time)
    y = np.array(y)
    return x, y, np.array(t)



x, y, time = read_vector_file("../data/comp_sol.txt")
x_t, y_t, time_t = read_vector_file("../data/analytical_solution.txt")
N = len(x)
U05 =  np.array([*y[time==time1][0], *y[time==time1][0]])
U05_t = np.array([*y_t[time_t==time1][0], *y_t[time_t==time1][0]])
U1 =  np.array([*y[time==time2][0], *y[time==time2][0]])
U1_t = np.array([*y_t[time_t==time2][0], *y_t[time_t==time2][0]])
X = np.linspace(-1/2, 3/2, 2*N)
plt.figure(figsize=(6, 5))
plt.plot(X, U1, color="seagreen", label=r"$\frac{ct}{L} = $ "+str(time2))
plt.plot(X, U05, color="blue", label=r"$\frac{ct}{L} = $ "+str(time1))
plt.plot(X, U1_t, color="gray", linestyle="dotted")#, label=r"$\frac{ct}{L} = 1$ (analytical)")
plt.plot(X, U05_t, color="gray", linestyle="dashed")#, label=r"$\frac{ct}{L} = 0.5$ (analytical)")

plt.legend(loc='upper left', ncol=2)
plt.xlabel(r"$\frac{\xi}{L}$")
plt.ylabel(r"$\frac{v}{U}$")
plt.grid()
#plt.savefig("../graphs/non_unif_8_I4_num.pdf", bbox_inches='tight')
plt.show()



try:
    x, y, time = read_vector_file("../data/comp_sol_phys.txt")
    x_t, y_t, time_t = read_vector_file("../data/analytical_solution_phys.txt")
    N = len(x)
    U05 =  np.array([*y[time==time1][0], *y[time==time1][0]])
    U05_t = np.array([*y_t[time_t==time1][0], *y_t[time_t==time1][0]])
    U1 =  np.array([*y[time==time2][0], *y[time==time2][0]])
    U1_t = np.array([*y_t[time_t==time2][0], *y_t[time_t==time2][0]])
    xi = np.linspace(-1/2, 3/2, 2*N)
    L = 1; a = 3/5
    mapping = lambda xi : xi - a * L / (2*np.pi) * np.sin(2*np.pi*xi/L)
    X = mapping(xi)
    plt.figure(figsize=(6, 5))
    plt.plot(X, U1, color="seagreen", label=r"$\frac{ct}{L} = $ " +str(time2))
    plt.plot(X, U05, color="blue", label=r"$\frac{ct}{L} = $ " +str(time1))
    plt.plot(X, U1_t, color="gray", linestyle="dotted")
    plt.plot(X, U05_t, color="gray", linestyle="dashed")
    plt.legend(loc='upper left', ncol=2)
    plt.xlabel(r"$\frac{x}{L}$")
    plt.ylabel(r"$\frac{u}{U}$")
    plt.grid()
    #plt.savefig("../graphs/non_unif_8_I4_phys.pdf", bbox_inches='tight')
    plt.show()
except:
    None 
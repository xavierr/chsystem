
# Debugging module
import ipdb

import matplotlib.pyplot as plt
import matplotlib
from numpy import *
import time
from io import *
import os
from subprocess import call

def addToFile(string_arg):
    file = FileIO("param_input.txt", 'a')
    file.write(string_arg)
    file.close()
    

# Class which contains computation settings
# Explanations of parameters  are given in __init__ function below.
class CondInit:
    def __init__(self, rho0 = 0., times = list(), time_end = 2., dt = 0.1,
                 N = 200, R = 20.0,  is_relabeling = "true"):
        self.rho0 = rho0
        self.times = times # times for which the solution is saved (typically dt is much smaller than intervals of this variable)
        self.time_end = time_end # end time of simulation
        self.dt = dt # time step size
        self.N = N # Number of point for the spatial discretization
        self.R = R # [-R, R] is the interval (in space) in which the solution is computed
        self.dxi = 2*R/N # space step size
        self.is_relabeling = is_relabeling # flag, true -> relabeling is done for each time in times. Otherwise, no relabeling

        self.xi = linspace(-R_arg, R_arg, N_arg)

    def writeToFile(self):
        # Routine which writes the parameters to file which we will be read in  c++ program.
        file = FileIO("param_input.txt", 'w')
        file.write("N = " + str(int(self.N)) + ";\n")
        file.write("R = " + str(self.R) + ";\n")
        file.write("time_end = " + str(self.time_end) + ";\n")
        file.write("dt = " + str(self.dt) + ";\n")
        file.write("time = [")
        for t in self.times[0:-1]:
            file.write(str(t) + ", ")
        file.write(str(self.times[-1]) + "];\n")
        file.write("is_relabeling = " + self.is_relabeling + ";\n")
        file.write("rho0 = " + str(self.rho0) + ";\n")
        file.close()


# Parameters for the current computation

# Parameters for "antipeakon" case
p1 = 1.0 
p2 = -1.0
q1 = -1.0
q2 = 1.0

# Parameter for "peakon" case
c = 1.0

# Parameters for "cubic" and "cubic2" case
gamma = 0.1
x0 = 1.0
umax = 1.0

# Parameter for "cubic2" case
rho0_x = 1.0

# Parameters common for all cases
R_arg = 10.0
time_end_arg = 0.4
dt_arg = 0.01
N_arg = 200
number_relabeling = 20
is_relabeling_arg = "false"
rho0_arg = 0.0
init = "cubic2"
# init = "cubic"
# init = "antipeakon"
# init = "peakon"


cond_init = CondInit(rho0 = rho0_arg, times = linspace(0.0, time_end_arg, number_relabeling), time_end = time_end_arg, dt = dt_arg, N = N_arg, R = R_arg, is_relabeling = is_relabeling_arg)
cond_init.writeToFile()

addToFile("init = \"" + init + "\";\n")
addToFile("p1 = " + str(p1) + ";\n")
addToFile("p2 = " + str(p2) + ";\n")
addToFile("q1 = " + str(q1) + ";\n")
addToFile("q2 = " + str(q2) + ";\n")
addToFile("c = " + str(c) + ";\n")
addToFile("x0 = " + str(x0) + ";\n")
addToFile("umax = " + str(umax) + ";\n")
addToFile("gamma = " + str(gamma) + ";\n") 
addToFile("rho0_x = " + str(rho0_x) + ";\n") 



def plot_result(cond_init_arg, save_fig = False, sol_computed = False, with_dots = True):

    # utility which plots results of computations as a movie (y, U, H, etc...)

    if not(sol_computed):
        print "Start Computation"
        call("./run.out")
        print "Computation done"
        
    y = fromfile("y.txt", sep=",");
    U = fromfile("U.txt", sep=",");
    q = fromfile("q.txt", sep=",");
    w = fromfile("w.txt", sep=",");
    r = fromfile("r.txt", sep=",");
    h = fromfile("h.txt", sep=",");
    t = fromfile("t.txt", sep=",");
    xi = fromfile("xi.txt", sep=",");
    M = fromfile("M.txt", sep=",");
    param = fromfile("param.txt", sep=",");

    U_max = amax(U)
    rho_max = amax(r/q)
    R = param[0]
    N = param[1]
        
    plt.ion() # interactive mode for matplotlib

    j = 0
    pos = 0
    for i, tt in enumerate(t):
        new_pos = pos + M[i]
        yt = y[pos : new_pos]
        Ut = U[pos : new_pos]
        qt = q[pos : new_pos]
        wt = w[pos : new_pos]
        rt = r[pos : new_pos]
        ht = h[pos : new_pos]
        pos = new_pos
        plt.figure(1)
        plt.clf()
        if with_dots:
            plt.plot(yt, Ut,'b*-')
        else:
            plt.plot(yt, Ut,'b-')
        plt.axis([-R, R, -U_max, U_max])
        plt.title(r'$u(t,x)$', fontsize = "xx-large", verticalalignment = "bottom")
        plt.draw()
        if save_fig:
            fname="u_%03d.png"%j
            plt.savefig(fname)
        plt.figure(2)
        plt.clf()
        plt.plot(linspace(-1, 1, M[i]), yt, 'b-')
        plt.axis([-1, 1, amin(y), amax(y)])
        plt.title(r'$y(t,\xi)$', fontsize = "xx-large", verticalalignment = "bottom")
        plt.draw()
        if save_fig:
            fname="y_%03d.png"%j
            plt.savefig(fname)
        plt.figure(3)
        plt.clf()
        plt.plot(linspace(-1, 1, M[i]), qt,'b-')
        plt.axis([-1, 1, amin(q), amax(q)])
        plt.title(r'$q(t,\xi)$', fontsize = "xx-large", verticalalignment = "bottom")
        plt.draw()
        plt.figure(4)
        plt.clf()
        plt.plot(linspace(-1, 1, M[i]), ht,'b-')
        plt.axis([-1, 1, amin(h), amax(h)])
        plt.title(r'$h(t,\xi)$', fontsize = "xx-large", verticalalignment = "bottom")
        plt.figure(5)
        plt.clf()
        plt.plot(linspace(-1, 1, M[i]), Ut,'b-')
        plt.axis([-1, 1, amin(U), amax(U)])
        plt.title(r'$U(t,\xi)$', fontsize = "xx-large", verticalalignment = "bottom")
        plt.figure(6)
        plt.clf()
        plt.plot(linspace(-1, 1, M[i]), wt, 'b-')
        plt.axis([-1, 1, amin(w), amax(w)])
        plt.title(r'$w(t,\xi)$', fontsize = "xx-large", verticalalignment = "bottom")
        plt.figure(7)
        plt.clf()
        plt.plot(linspace(-1, 1, M[i]), qt + ht, 'b-')
        plt.axis([-1, 1, amin(q + h), amax(q + h)])
        plt.title(r'$q(t,\xi)+h(t,\xi)$', fontsize = "xx-large", verticalalignment = "bottom")
        plt.figure(8)
        plt.clf()
        plt.plot(linspace(-1, 1, M[i]), rt, 'b-')
        plt.axis([-1, 1, amin(r), amax(r)])
        plt.title(r'$r(t,\xi)$', fontsize = "xx-large", verticalalignment = "bottom")
        plt.draw()
        plt.figure(9)
        plt.clf()
        if with_dots:
            plt.plot(yt, rt/qt,'b*-')
        else:
            plt.plot(yt, rt/qt,'b-')
        plt.axis([-R, R, -rho_max, rho_max])
        plt.title(r'$\rho(t,x)$', fontsize = "xx-large", verticalalignment = "bottom")
        plt.draw()
        time.sleep(0.01)
            
def plot_charac(cond_init_arg, save_fig = False, file_name = "charac.png",  sol_computed = False, xic = [-10, 0, 10]):

    # utility to plot charachteristics.

    if not(sol_computed):
        print "debut calcul"
        call("./run.out")
        print "fin calcul"
        
    y = fromfile("y.txt", sep=",");
    U = fromfile("U.txt", sep=",");
    q = fromfile("q.txt", sep=",");
    r = fromfile("r.txt", sep=",");
    h = fromfile("h.txt", sep=",");
    t = fromfile("t.txt", sep=",");
    xi = fromfile("xi.txt", sep=",");
    param = fromfile("param.txt", sep=",");

    N = cond_init_arg.N
    R = cond_init_arg.R

    ind_xic = empty_like(xic)
    for i, xici in enumerate(xic):
        ind_xic[i] = round((xici + R)/cond_init_arg.dxi)
    
    plt.ion() # interactive mode for matplotlib
    
    plt.figure(5)
    plt.clf()
    for ind_xii in ind_xic:
        yt = empty_like(t)
        for ind_tt in enumerate(t):
            yt[ind_tt[0]] = y[ind_xii + ind_tt[0]*(N + 1)] 
        plt.plot(yt, t)
    plt.axis([cond_init_arg.xi[0], cond_init_arg.xi[-1], 0.0, t[-1]])
    if save_fig:
        plt.savefig(file_name)

def plot_comp_result(save_fig = False,  
                     with_dots = False,
                     create_movie = False,
                     save_fig_directory = "Computations"):

    # utility which plots different results of computations.

    class Datas:
        def __init__(self, directory = "S1"):
            self.y = fromfile(directory + "/y.txt", sep=",");
            self.U = fromfile(directory + "/U.txt", sep=",");
            self.q = fromfile(directory + "/q.txt", sep=",");
            self.w = fromfile(directory + "/w.txt", sep=",");
            self.r = fromfile(directory + "/r.txt", sep=",");
            self.h = fromfile(directory + "/h.txt", sep=",");
            self.t = fromfile(directory + "/t.txt", sep=",");
            self.xi = fromfile(directory + "/xi.txt", sep=",");
            self.M = fromfile(directory + "/M.txt", sep=",");
            self.param = fromfile(directory + "/param.txt", sep=",");
            self.U_max = amax(abs(self.U))
            self.R = self.param[0]
            self.N = self.param[1]

        def compute_sol(self, index = 0):
            pos = 0
            for k in arange(0, index):
                pos = pos + self.M[k]
            self.yt = self.y[pos : pos + self.M[index]]
            self.Ut = self.U[pos : pos + self.M[index]]
            self.qt = self.q[pos : pos + self.M[index]]
            self.wt = self.w[pos : pos + self.M[index]]
            self.rt = self.r[pos : pos + self.M[index]]
            self.ht = self.h[pos : pos + self.M[index]]

            
    if create_movie:
        save_fig = True

    plt.ion() # interactive mode for matplotlib

    
    sol1 = Datas(directory = "Computations/S1")
    sol2 = Datas(directory = "Computations/S2")
    sol3 = Datas(directory = "Computations/S3")
    U_max = amax([sol1.U_max, sol2.U_max, sol3.U_max])

    for i, tt in enumerate(sol1.t):
        sol1.compute_sol(index = i)
        sol2.compute_sol(index = i)
        sol3.compute_sol(index = i)
        fig = plt.figure(1)
        plt.clf()
        if with_dots:
            plt.plot(sol1.yt, sol1.Ut,'b*-', sol2.yt, sol2.Ut, 'k*-', sol3.yt, sol3.Ut, 'r*-')
        else:
            plt.plot(sol1.yt, sol1.Ut,'b-', sol2.yt, sol2.Ut, 'k-', sol3.yt, sol3.Ut, 'r-')

        plt.axis([-sol1.R, sol1.R, -U_max, U_max])
        plt.title(r'$u(t,x)$', fontsize = "xx-large", verticalalignment = "bottom")
        ax = fig.get_axes()
        plt.text(0.8, 0.8, 'time: %.2f'%(tt),
             horizontalalignment='center',
             verticalalignment='center',
             transform = ax[0].transAxes)
        plt.draw()
        if save_fig:
            fname="%s/_tmp_%03d.png"%(save_fig_directory, i)
            plt.savefig(fname)
        time.sleep(0.5)
    if create_movie:
        os.system('mencoder \"mf://' + save_fig_directory + '/_tmp_*.png\" -mf type=png fps=1 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o animation.avi')

    
# plot_charac(cond_init, sol_computed = False)
# plot_charac(cond_init_1, file_name = "charac1.png")
# plot_result(cond_init, sol_computed = False, save_fig = False)
# plot_result(cond_init_1, True)

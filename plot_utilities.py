
# script which generates pictures for article

# Debugging module
import ipdb

import matplotlib.pyplot as plt
import matplotlib
from numpy import *
import time
from io import *
import os
from subprocess import call

# ipdb.set_trace()

matplotlib.rcParams['lines.linewidth'] = 2

def addToFile(string_arg):
    file = FileIO("param_input.txt", 'a')
    file.write(string_arg)
    file.close()
    

class CondInit:
    def __init__(self, rho0 = 0., rho0_x = 2.0, times = list(), time_end = 2., dt = 0.1,
                 abs_err = 1e-6, rel_err = 1e-6, N = 200, R = 20.0, lower_bound = 1e-3, is_relabeling = 1, gamma = 1, xi = list()):
        self.rho0 = rho0
        self.rho0_x = rho0_x
        self.gamma = gamma
        self.times = times
        self.time_end = time_end
        self.dt = dt
        self.abs_err = abs_err
        self.rel_err = rel_err
        self.N = N
        self.R = R
        self.dxi = 2*R/N
        self.lower_bound = lower_bound
        self.is_relabeling = is_relabeling
        self.xi = xi

    def writeToFile(self, init = "peakon"):
        file = FileIO("param_input.txt", 'w')
        file.write("init = \"" + init + "\";\n")
        file.write("N = " + str(int(self.N)) + ";\n")
        file.write("R = " + str(self.R) + ";\n")
        file.write("lower_bound = " + str(self.lower_bound) + ";\n")
        file.write("time_end = " + str(self.time_end) + ";\n")
        file.write("dt = " + str(self.dt) + ";\n")
        file.write("rho0 = " + str(self.rho0) + ";\n")
        file.write("rho0_x = " + str(self.rho0_x) + ";\n")
        file.write("gamma = " + str(self.gamma) + ";\n")
        file.write("time = [")
        for t in self.times[0:-1]:
            file.write(str(t) + ", ")
        file.write(str(self.times[-1]) + "];\n")
        file.write("abs_err = " + str(self.abs_err) + ";\n")
        file.write("rel_err = " + str(self.rel_err) + ";\n")
        file.write("is_relabeling = " + str(self.is_relabeling) + ";\n")
        file.close()


p1 = 1.0
p2 = -1.0
q1 = -1.0
q2 = 1.0
c = 1.0
x0 = 1.0
umax = 1.0
R_arg = 10.0
col_time = 1.783
time_end_arg = 5.0
dt_arg = 0.01
N_arg = 2000
number_relabeling = 30
lower_bound_arg = 1.0e-4
is_relabeling_arg = 0
rho0_arg = 0.01
rho0_x_arg = 1.0
gamma_arg = 0.1


cond_init = CondInit(rho0 = rho0_arg, rho0_x = rho0_x_arg, times = linspace(0.0, time_end_arg, number_relabeling), time_end = time_end_arg, dt = dt_arg,
                     abs_err = 1e-9, rel_err = 1e-9, N = N_arg, R = R_arg, lower_bound = lower_bound_arg, is_relabeling = is_relabeling_arg,
                     xi = linspace(-R_arg, R_arg, N_arg), gamma = gamma_arg)



init = "cubic2"
# init = "cubic"
# init = "antipeakon"
# init = "peakon"

def plot_result(cond_init_arg, save_fig = False, sol_computed = False, with_dots = True):

    cond_init_arg.writeToFile(init)
    addToFile("p1 = " + str(p1) + ";\n")
    addToFile("p2 = " + str(p2) + ";\n")
    addToFile("q1 = " + str(q1) + ";\n")
    addToFile("q2 = " + str(q2) + ";\n")
    addToFile("c = " + str(c) + ";\n")
    addToFile("x0 = " + str(x0) + ";\n")
    addToFile("umax = " + str(umax) + ";\n")

    if not(sol_computed):
        print "Start Computation"
        call("./a.out")
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

    # rho_max = amax(r/q)
    # ener_dens_max = amax(h/q - (r/q)**2)
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
            
def plot_charac(cond_init_arg, file_name = False,  sol_computed = False):
    cond_init_arg.writeToFile(init)
    addToFile("p1 = " + str(p1) + ";\n")
    addToFile("p2 = " + str(p2) + ";\n")
    addToFile("q1 = " + str(q1) + ";\n")
    addToFile("q2 = " + str(q2) + ";\n")
    addToFile("c = " + str(c) + ";\n")
    addToFile("x0 = " + str(x0) + ";\n")
    addToFile("umax = " + str(umax) + ";\n")

    if not(sol_computed):
        print "debut calcul"
        call("./a.out")
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

    xii = empty_like(cond_init_arg.xi)
    for i, xic in enumerate(cond_init_arg.xi):
        xii[i] = round((xic + R)/cond_init_arg.dxi)
    
    plt.ion() # interactive mode for matplotlib
    
    plt.figure(5)
    plt.clf()
    # for ind_xii in xii:
    for ind_xii in arange(N+1):
        yt = empty_like(t)
        for ind_tt in enumerate(t):
            yt[ind_tt[0]] = y[ind_xii + ind_tt[0]*(N + 1)] 
        plt.plot(yt, t)
    plt.axis([cond_init_arg.xi[0], cond_init_arg.xi[-1], 0.0, t[-1]])
    if file_name:
        plt.savefig(file_name)

def plot_comp_result(save_fig = False,  
                     with_dots = False,
                     create_movie = False,
                     save_fig_directory = "Computations"):

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

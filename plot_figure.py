
# plot the figure in article 

import matplotlib.pyplot as plt
import matplotlib
from numpy import *
import time
from io import *
import os
from subprocess import call


def plot_comp_result(save_fig = False,  
                     with_dots = False):

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
    sol2 = Datas(directory = "Computations/S3")
    U_max = amax([sol1.U_max, sol2.U_max])

    for i, tt in enumerate(sol1.t):
        sol1.compute_sol(index = i)
        sol2.compute_sol(index = i)
        fig = plt.figure(1)
        plt.clf()
        if with_dots:
            plt.plot(sol1.yt, sol1.Ut,'b*-', sol2.yt, sol2.Ut, 'k*-')
        else:
            plt.plot(sol1.yt, sol1.Ut,'b-', sol2.yt, sol2.Ut, 'k-')

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


plot_comp_result()

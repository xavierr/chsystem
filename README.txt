

--------------------------------------------------------------
 Requirements for the compilation of the C++ program chsysdiss
--------------------------------------------------------------

* libconfig library 

  http://www.hyperrealm.com/libconfig/

* boost library

  http://www.boost.org

  (the library consists of only header files)

* odeint library

  http://www.odeint.com

  (the library consists of only header files)

* c++ compilator which can handle c++0x standard.


----------------------------------
Example of configuration (my own!)
----------------------------------

* The library libconfig is installed as a shared library
(in ubuntu run in terminal:  sudo apt-get install libconfig)

* The header files of the boost library are installed under the path  /path/to/boost  (dummy name in the example!)

* The header files of the odeint library are installed under the path  /path/to/odeint  (dummy name in the example!)
 
* The current directory where chsysdiss.cpp and chsysdiss.hpp are is /path/to/chsysdiss (dummy name in the example!)

* Compiler g++ version 4.6.3

Then to compile, run in the current directory:

g++ -std=c++0x -O3 chsysdiss.cpp -lconfig++ -I /path/to/chsysdiss/ -I /path/to/boost/ -I /path/to/odeint/ -o run.out

------------------------------
Parameter file param_input.txt
------------------------------

The executable run.out reads the parameters for the computation in the file
param_input.txt.

The significance of the parameters are documented in the python file
plot_utilities.py, which also have functionnality to write the parameter file.

-----------------------------------------------------------------
Requirement for python files plot_utilities.py and plot_figure.py
-----------------------------------------------------------------

* The following python modules:
  * matplotlib (http://matplotlib.org/)
  * numpy (http://numpy.scipy.org/)

* ipython (http://ipython.org/)

-------------------------
Run your first simulation
-------------------------

* Launch ipython

* In ipython
  * load plot_utilities.py, that is, write:
    %run plot_utilities

  * run plotting function, that is, write:
    plot_result(cond_init)

* Change parameters by editing file plot_utilities.py and rerun the simulation
  as described above.

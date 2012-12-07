// Compilation options
// g++   -std=c++0x chsysdiss.cpp -lconfig++ -I /home/xavier/Notes/NonVanDissCH/ -I /home/xavier/boost-trunk/ -I /home/xavier/odeint-v2/
// g++   -std=c++0x chsysdiss.cpp -lconfig++ -I /home/xavier/Notes/CMA/NonVanDissCH/ -I /home/xavier/lib/boost-trunk/ -I /home/xavier/lib/odeint-v2/

#include <libconfig.h++>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <vector>
#include <chsysdiss.hpp>
#include <fstream>
#include <sstream>
#include <algorithm>

namespace {


  double interpolate(double y1, double y2, double u1, double u2, double xi) {
    double a = (y2 - xi)/(y2 - y1);
    return a*u1 + (1 - a)*u2;
  };
  
  double generateXi(double xi_left, double xi_right, int N,
		    std::vector<double>& xi) {
    double dxi = (xi_right - xi_left)/N;
    for (int i = 0; i < N + 1; ++i) {
      xi[i] = xi_left + i*dxi;
    }
    return dxi;
  }

  void computeDerivative(const double* y,
			 const double* U,
			 const double* H,
			 const double* q,
			 const double* w,
			 const double* h,
			 const double* r,
			 const std::vector<bool>& diss,
			 const std::vector<double>& P,
			 const std::vector<double>& Q,
			 std::vector<double>& dZdt,
			 const int N) {
    for (int i = 0; i < N + 1; ++i) {
      dZdt[i + (N+1)*0] = U[i];
      dZdt[i + (N+1)*1] = -Q[i];
      dZdt[i + (N+1)*2] = std::pow(U[i],3)-2*P[i]*U[i];
      if (diss[i] == false) {
	dZdt[i + (N+1)*3] = w[i];
	dZdt[i + (N+1)*4] = 0.5*h[i]+(0.5*U[i]*U[i]-P[i])*q[i];
	dZdt[i + (N+1)*5] = -2*Q[i]*U[i]*q[i]+(3*U[i]*U[i]-2*P[i])*w[i];
      } else {
	dZdt[i + (N+1)*3] = 0.0;
	dZdt[i + (N+1)*4] = 0.0;
	dZdt[i + (N+1)*5] = 0.0;
      }      
      dZdt[i + (N+1)*6] = 0.0;
    }    
  }
  
  void computePandQ(const double* y,
		    const double* U,
		    const double* H,
		    const double* q,
		    const double* w,
		    const double* h,
		    const double* r,
		    const std::vector<bool>& diss,
		    std::vector<double>& P,
		    std::vector<double>& Q,
		    const double dxi,
		    const int N) {

  for (int j = 0; j < N + 1; ++j) {
    P[j] = 0.0;
    Q[j] = 0.0;
    for (int i = 0; i < j; ++i) {
      if (diss[i] == false) {
	double coef =  exp(-(y[j] - y[i]))/4.0 *(h[i] + U[i]*U[i]*q[i])*dxi;
	P[j] += coef;
	Q[j] -= coef;
      }
    }
    for (int i = j; i < N + 1; ++i) {
      if (diss[i] == false) {
	double coef = exp((y[j] - y[i]))/4.0*(h[i] + U[i]*U[i]*q[i])*dxi;
	P[j] += coef;
	Q[j] += coef;
      }
    }
  }
}


  // double gFunction(double y, double U, double H,
  // 		   double q, double w, double h,
  // 		   double r) {
  //   double g2 = q + h + r*r;
  //   if (r == 0) {
  //     if ( w < 0 ) {
  // 	double g1 = -w + 2*(1 + U*U)*q;
  // 	if (g1 < g2) {
  // 	  return g1;
  // 	} else {
  // 	  return g2;
  // 	} 
  //     } else {
  // 	return g2;
  //     }
  //   } else {
  //     return g2;
  //   }
  // }

  // double yPlusHXi(double y, double U, double H,
  // 		  double q, double w, double h,
  // 		  double r) {
  //   return q + h;
  // }



  void writeVectorToFile(const std::string file_name, const std::vector<double>& vector) {
    std::ofstream file(file_name);
    for (std::vector<double>::const_iterator vec = vector.begin(); vec < vector.end(); ++vec) {
      file << *vec << ",";
    }
    file.close();
  }

  double divide(double x1, double x2) {
    if (fabs(x2) < 1e-16) {
      return 0.;
    } else {
      return x1/x2;
    }
  };

  void interpolateVector(std::vector<double>::const_iterator xi_begin, 
  			 std::vector<double>::const_iterator xi_end,
  			 const double * u,
  			 std::vector<double>::const_iterator phi_begin, 
  			 std::vector<double>::const_iterator phi_end,  
  			 double*  u_phi, 
  			 bool is_u_derivative, 
  			 bool is_u_equal_y = false)  {
    if (phi_begin < phi_end) {
      if (*phi_begin < *xi_begin) {
  	// We are outside the scope of phi. 
  	if (is_u_equal_y) { 
  	  // extrapolate with line of slope one.
  	  *u_phi = *u + (*phi_begin - *xi_begin);
  	} else {
  	  // extrapolate with horizontal line.
  	  *u_phi = *u;
  	}
  	interpolateVector(xi_begin, xi_end, u, phi_begin + 1, phi_end, u_phi + 1, is_u_derivative, is_u_equal_y);
      } else {
  	while ((xi_begin + 1 < xi_end) && (*phi_begin >= *(xi_begin + 1))) {
  	  ++xi_begin;
  	  ++u;
  	}
  	if (xi_begin + 1 < xi_end ) {
  	  if (is_u_derivative) {
  	    *u_phi = *u;
  	  } else {
  	    *u_phi = interpolate(*xi_begin, *(xi_begin + 1),
				 *u, *(u + 1),
				 *phi_begin);
  	  }
  	  interpolateVector(xi_begin, xi_end, u, phi_begin + 1, phi_end, u_phi + 1, is_u_derivative, is_u_equal_y);
  	} else {
  	  // We are outside the scope of phi. We extrapolate as above.
  	  if (is_u_equal_y) {
  	    *u_phi = *u + (*phi_begin - *xi_begin);
  	  } else {
  	    *u_phi = *u;
  	  }
  	  interpolateVector(xi_begin, xi_end, u, phi_begin + 1, phi_end, u_phi + 1, is_u_derivative, is_u_equal_y);
  	}
      }
    } else {
      return;
    }
  };

  void relabelingAndCopy(const std::vector<double>& xi,
  			 double* u,
  			 const std::vector<double>& g_inv,
			 std::vector<double>& new_u, const bool is_u_equal_y = false) 
  {
    interpolateVector(xi.begin(), xi.end(),
  		      u, 
  		      g_inv.begin(), g_inv.end(),  
  		      &new_u[0], 
  		      false, is_u_equal_y);
    std::copy(new_u.begin(), new_u.end(), u);
  };

  void relabelingAndCopyDer(const std::vector<double>& xi,
			    double* q,
			    const std::vector<double>& g_inv,
			    const std::vector<double>& g_xi,
			    std::vector<double>& new_q) 
  {

    const int N = xi.size();
    std::transform(q, q + N, g_xi.begin(), new_q.begin(), divide);  

    interpolateVector(xi.begin(), xi.end(),
  		      &new_q[0], 
  		      g_inv.begin(), g_inv.end(),  
  		      q, 
		      true, false);
  };


  void integrateDer(double*  y,
		    const double* y_xi,
		    const double dxi,
		    const int ind0,
		    const double y0,
		    const int N) 
  {

    double* y_p = y + ind0; 
    double* y_n = y + ind0; 
    const double* y_xi_p = y_xi + ind0; 
    const double* y_xi_n = y_xi + ind0; 

    y_p[0] = y0;
    while(y_p < y + N) {
      ++y_p;
      ++y_xi_p;
      y_p[0] = y_p[-1] + y_xi_p[-1]*dxi;
    }

    while(y_n > y) {
      --y_n;
      --y_xi_n;
      y_n[0] = y_n[1] - y_xi_n[0]*dxi;
    }
  }
} // anonymous namespace.

bool Simulation::checkIfInCollisionRegion(double y1, double U1, double H1, double q1, double w1, double h1, double r1) {
  // if ((r1 == 0) && (w1 <= 0) && (-w1 + 2*(1 + U1*U1)*q1 <= q1 + h1) && (q1 >= 0)) {
  if ((r1 == 0) && (w1 <= 0) && (-w1 + 2*(1 + U1*U1)*q1 <= q1 + h1) && (q1 >= 0)) {
    return true;
  } else {
    return false;
  }
}

bool Simulation::checkIfCollided(double y1, double U1, double H1, double q1, double w1, double h1, double r1) {
  double radius = (q1 + h1)/(2*(U1*U1 + 1));
  if ((r1 == 0) && ((q1 <=0) || ( (q1 >=0) && (w1 >= 0) && (w1 <= collision_limiter*(radius - q1)/radius) ))){
    return true;
  } else {
    return false;
  }
}

Simulation::Simulation()  {

  libconfig::Config cfg;
  cfg.readFile("param_input.txt");
  R = cfg.lookup("R");
  N = cfg.lookup("N");
  time_end = cfg.lookup("time_end");
  dt = cfg.lookup("dt");

  libconfig::Setting& time_setting = cfg.lookup("time");
  int size_time = time_setting.getLength();
  time.resize(size_time);
  for (int i = 0; i < size_time; ++i) {
    time[i] = time_setting[i];
  }
  dxi = 2*R/N;
  xi.resize(N + 1);
  for (int i = 0; i < N+1; ++i) {
    xi[i] = -R + i*dxi;
  }
  new_xi.resize(N + 1);
  track.resize(N + 1);
  old_track.resize(N + 1);
  g.resize(N + 1);
  g_inv.resize(N + 1);
  g_xi.resize(N + 1);
  diss.resize(N + 1);
  P.resize(N + 1);
  Q.resize(N + 1);
  local_var1.resize(N + 1);
  Z.resize(7*(N + 1));
  resetPointers(Z);
}

void Simulation::initPeakon() {
  resetPointers(Z);
  libconfig::Config cfg;
  cfg.readFile("param_input.txt");
  const double rho0 = cfg.lookup("rho0");
  const double c = cfg.lookup("c");
  double q1 = 0.0;
  for (int i = 0; i < N + 1; ++i) {
    y[i] = xi[i];
    U[i] = c*exp(-c*fabs(xi[i]));
    w[i] = -c*sign(xi[i])*exp(-c*fabs(xi[i]));
    q[i] = 1.;
    r[i] = rho0;
    h[i] = U[i]*U[i] + w[i]*w[i] + r[i]*r[i];
  }
  int ind0 = ceil((N - 1)/2);
  integrateDer(H, h, dxi, ind0, 0, N);
}

void Simulation::initCubic() {
  resetPointers(Z);
  // u = alpha*exp(-x^2)*x*(x-x0)*(x+x0) 
  libconfig::Config cfg;
  cfg.readFile("param_input.txt");
  const double rho0 = cfg.lookup("rho0");
  const double rho0_x = cfg.lookup("rho0_x");
  const double x0 = cfg.lookup("x0");
  const double umax = cfg.lookup("umax");
  const double xmax = 1/sqrt(3)*x0; // Point where u' vanishes
  const double alpha = umax/(2*pow(xmax, 3));
  const double beta = 1.0;
  for (int i = 0; i < N + 1; ++i) {
    y[i] = xi[i];
    U[i] = alpha*exp(-beta*xi[i]*xi[i])*xi[i]*(xi[i] - x0)*(xi[i] + x0);
    w[i] = alpha*exp(-beta*xi[i]*xi[i])*(3*xi[i]*xi[i] - x0*x0)-2*beta*xi[i]*U[i];
    q[i] = 1.;
    if (fabs(xi[i]) <= rho0_x) {
      r[i] = rho0;
    }
    h[i] = U[i]*U[i] + w[i]*w[i] + r[i]*r[i];
  }
  int ind0 = ceil((N - 1)/2);
  integrateDer(H, h, dxi, ind0, 0, N);
}

void Simulation::initCubic2() {
  resetPointers(Z);
  // u = alpha*exp(-x^2)*x*(x-x0)*(x+x0) 
  libconfig::Config cfg;
  cfg.readFile("param_input.txt");
  const double rho0 = cfg.lookup("rho0");
  const double gamma = cfg.lookup("gamma");
  const double x0 = cfg.lookup("x0");
  const double umax = cfg.lookup("umax");
  const double xmax = 1/sqrt(3)*x0; // Point where u' vanishes
  const double alpha = umax/(2*pow(xmax, 3));
  const double beta = 1.0;
  for (int i = 0; i < N + 1; ++i) {
    y[i] = xi[i];
    U[i] = alpha*exp(-beta*xi[i]*xi[i])*xi[i]*(xi[i] - x0)*(xi[i] + x0);
    w[i] = alpha*exp(-beta*xi[i]*xi[i])*(3*xi[i]*xi[i] - x0*x0)-2*beta*xi[i]*U[i];
    q[i] = 1.;
    r[i] = rho0*exp(-gamma*xi[i]*xi[i]);
    h[i] = U[i]*U[i] + w[i]*w[i] + r[i]*r[i];
  }
  int ind0 = ceil((N - 1)/2);
  integrateDer(H, h, dxi, ind0, 0, N);
}


void Simulation::initAntipeakon() {
  resetPointers(Z);
  for (int i = 0; i < N + 1; ++i) {
    y[i] = xi[i];
  }
  libconfig::Config cfg;
  cfg.readFile("param_input.txt");
  const double p1 = cfg.lookup("p1");
  const double p2 = cfg.lookup("p2");
  const double q1 = cfg.lookup("q1");
  const double q2 = cfg.lookup("q2");
  const double rho0 = cfg.lookup("rho0");
  const double rho0_x = cfg.lookup("rho0_x");
  
  for (int i = 0; i < N + 1; ++i) {
    U[i] = p1*exp(-fabs(y[i] - q1)) + p2*exp(-fabs(y[i] - q2));
    w[i] = -sign(y[i]-q1)*p1*exp(-fabs(y[i]-q1))
      -sign(y[i]-q2)*p2*exp(-fabs(y[i]-q2));
    r[i] = rho0*exp(-10*std::pow(xi[i]/R, 2.));
    if (fabs(xi[i]) <= rho0_x) {
      r[i] = rho0;
    }
    q[i] = 1.0;
    h[i] = U[i]*U[i] + w[i]*w[i] + r[i]*r[i];
  }
  int ind0 = ceil((N - 1)/2);
  integrateDer(H, h, dxi, ind0, 0, N);
}


void Simulation::relabelingAndCopyAll(){

  resetPointers(Z);
  computeG();
  double new_dxi = generateXi(g[0], g[N], N, new_xi);

  interpolateVector(g.begin(), g.end(),
  		    &xi[0], 
  		    new_xi.begin(), new_xi.end(),  
  		    &g_inv[0], 
  		    false, true);

  relabelingAndCopy(xi, y, g_inv, local_var1, true);
  relabelingAndCopy(xi, U, g_inv, local_var1, false);
  relabelingAndCopy(xi, H, g_inv, local_var1, true);


  relabelingAndCopyDer(xi, q, g_inv, g_xi, local_var1);
  relabelingAndCopyDer(xi, w, g_inv, g_xi, local_var1);
  relabelingAndCopyDer(xi, h, g_inv, g_xi, local_var1);

  std::copy(new_xi.begin(), new_xi.end(), xi.begin());
  
  dxi = new_dxi;

  setTrack();
  setDiss();

}

void Simulation::solve() {

  typedef boost::numeric::odeint::runge_kutta_dopri5< state_type > stepper_type;
  stepper_type stepper = boost::numeric::odeint::runge_kutta_dopri5< state_type >();
  double t = 0;
  solution.push_back(Z);
  solution_diss.push_back(diss);
  solution_time.push_back(0.0);
  for (std::vector<double>::const_iterator time_it = time.begin(); time_it < time.end() - 1; ++time_it) {
    while (t <= time_it[1]) {
      resetPointers(Z);
      integrate_const(stepper, *this, Z, t, t + dt, dt);
      removeEnergy();
      t += dt;
    }
    resetPointers(Z);
    // fixYandH();
    if (is_relabeling) {
      relabelingAndCopyAll();
    } 
    
    solution.push_back(Z);
    solution_diss.push_back(diss);
    solution_time.push_back(*time_it);
  }
}    
  
void Simulation::operator() (const state_type& ZZ, state_type& dZdt, const double /* t */ ) {

  std::vector<double> PP(N);
  std::vector<double> QQ(N);

  const double* yy = &ZZ[0*(N + 1)];
  const double* UU = &ZZ[1*(N + 1)];
  const double* HH = &ZZ[2*(N + 1)];
  const double* qq = &ZZ[3*(N + 1)];
  const double* ww = &ZZ[4*(N + 1)];
  const double* hh = &ZZ[5*(N + 1)];
  const double* rr = &ZZ[6*(N + 1)];
  
  ::computePandQ(yy, UU, HH,
		 qq, ww, hh, rr, 
		 diss, 
		 PP, QQ, dxi, N);
	       
  computeDerivative(yy, UU, HH,
		    qq, ww, hh, rr,
		    diss,
		    PP, QQ, 
		    dZdt,
		    N);
}


void Simulation::computeG() {

  std::transform(q, q + N + 1, h, g_xi.begin(), [](double a, double b) {
      return a + b;
    });

  int ind0 = ceil((N - 1)/2);
  integrateDer(&g[0], &g_xi[0], dxi, ind0, 0.0, N);

}



void Simulation::computePandQ() {

  resetPointers(Z);
  ::computePandQ(y, U, H,
		 q, w, h, r, diss,
		 P, Q, dxi, N);
}

  
void Simulation::resetPointers(state_type& ZZ) {
  y = &ZZ[0*(N + 1)];
  U = &ZZ[1*(N + 1)];
  H = &ZZ[2*(N + 1)];
  q = &ZZ[3*(N + 1)];
  w = &ZZ[4*(N + 1)];
  h = &ZZ[5*(N + 1)];
  r = &ZZ[6*(N + 1)];
}


void Simulation::writeToFile()  {
  std::ofstream param_file("param.txt");
  std::ofstream xi_file("xi.txt");
  std::ofstream t_file("t.txt");
  std::ofstream y_file("y.txt");
  std::ofstream U_file("U.txt");
  std::ofstream H_file("H.txt");
  std::ofstream q_file("q.txt");
  std::ofstream w_file("w.txt");
  std::ofstream M_file("M.txt");
  std::ofstream h_file("h.txt");
  std::ofstream r_file("r.txt");
  double* diss_p;
  param_file << R <<"," << N;
  for (int i = 0; i < N + 1; ++i) {
    xi_file << xi[i] << ",";
  }
  for (int j = 0; j < solution.size(); ++j) {
    t_file << solution_time[j] << ",";
    resetPointers(solution[j]);
    const std::vector<bool>& diss_p = solution_diss[j];
    // setDiss();
    int M = N + 1;
    for (int i = 0; i < N + 1; ++i) {
      if (diss_p[i] == true) {
      	--M;
      } else {
	y_file << y[i] << ",";
	U_file << U[i] << ",";
	H_file << H[i] << ",";
	q_file << q[i] << ",";
	w_file << w[i] << ",";
	h_file << h[i] << ",";
	r_file << r[i] << ",";
      }
    }
    M_file << M <<",";
  }
  xi_file.close();
  t_file.close();
  y_file.close();
  U_file.close();
  H_file.close();
  q_file.close();
  M_file.close();
  w_file.close();
  h_file.close();
  r_file.close();
}

void Simulation::init() {

  libconfig::Config cfg;
  cfg.readFile("param_input.txt");
  const std::string init_type = cfg.lookup("init");
  is_relabeling = cfg.lookup("is_relabeling");
  if (init_type == "peakon") {
    initPeakon();
  } else if (init_type == "antipeakon") {
    initAntipeakon();
  } else if (init_type == "cubic") {
    initCubic();
  } else if (init_type == "cubic2") {
    initCubic2();
  }
  diss.assign(N + 1, false);
  collision_limiter = 0.0;
  if (is_relabeling) {
    relabelingAndCopyAll();
  } else {
    setTrack();
  }
}

void Simulation::setTrack() {
  for (int i = 0; i < N + 1; ++i) {
    track[i] = checkIfInCollisionRegion(y[i], U[i], H[i], q[i], w[i], h[i], r[i]);
  }
}

void Simulation::setDiss() {
  for (int i = 0; i < N + 1; ++i) {
    diss[i] = checkIfCollided(y[i], U[i], H[i], q[i], w[i], h[i], r[i]);
  }
}

void Simulation::removeEnergy() {

  resetPointers(Z);
  std::copy(track.begin(), track.end(), old_track.begin());
  setTrack();

  bool need_set_track = false;
  for (int i = 0; i < N + 1; ++i) {
    if (old_track[i] && not(track[i])) {
      // if ((q[i] < (q[i] + h[i])/(2*(1 + U[i]*U[i]))) && (q[i] >= 0) && (w[i] >= 0)) {
      if ((q[i] < (q[i] + h[i])/(2*(1 + U[i]*U[i])))) {
	// collision_limiter = max(collision_limiter, 10*w[i]);
	diss[i] = true;
	q[i] = 0.0;
	w[i] = 0.0;
	need_set_track = true;
      }
    }
  }
  if (need_set_track) {
    setTrack();
  }
}

int main() {
  Simulation sim;
  sim.init();
  sim.solve();
  sim.writeToFile();
}

    
    


#include <libconfig.h++>
#include <vector>
#include <string>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/dense_output_runge_kutta.hpp>
#include <ostream>

#ifndef CHSYS_HPP_INCLUDED
#define CHSYS_HPP_INCLUDED

double sign(double x) {
  return (0.0 < x) - (0.0 > x);
}

typedef std::vector< double > state_type;
typedef std::vector< bool > diss_type;

class Simulation {

public:

  double R;
  int N;
  double dxi;
  double lower_bound;
  std::vector< double > xi;
  std::vector< double > new_xi;
  double* y;
  double* U;
  double* H;
  double* q;
  double* w;
  double* h;
  double* r;
  double collision_limiter;
  std::vector<double> local_var1; // internal working variable
  std::vector< bool > diss;
  std::vector< bool > track;
  std::vector< bool > old_track;
  std::vector< double > g;
  std::vector< double > g_inv;
  std::vector< double > g_xi;
  std::vector< double > P;
  std::vector< double > Q;
  state_type Z;
  std::vector< state_type > solution;
  std::vector< std::vector<bool> > solution_diss;
  std::vector< double > time;
  std::vector< double > solution_time;
  bool is_relabeling;
  double time_end;
  double dt;

  Simulation();
  void resetPointers(state_type& Z);
  void computeG();
  void setTrack();
  void setDiss();
  void init();
  void initAntipeakon();
  void initPeakon();
  void initCubic();
  void initCubic2();
  void operator()(const state_type& ZZ, state_type& dZdt, const double /* t */ );
  void solve();
  void computePandQ();
  void writeToFile();
  void relabelingAndCopyAll();
  void removeEnergy();
  bool checkIfInCollisionRegion(double y1, double U1, double H1, double q1, double w1, double h1, double r1);
  bool checkIfCollided(double y1, double U1, double H1, double q1, double w1, double h1, double r1);
};

#endif // CHSYS_HPP_INCLUDED

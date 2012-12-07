
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

// class YAndPos {
// public:
//   double y;
//   int pos;
// };

// class PartFunction {
// public:
//   const std::vector<double>& getVal() {
//     return val_;
//   };

//   int getInd_mid() {
//     return ind_mid_;
//   };
  
//   void assign (int ind_left, int ind_mid, int ind_right, int N);

// private:
//   std::vector<double> val_;
//   int N_;
//   int ind_left_;
//   int ind_mid_;
//   int ind_right_;
// };


// typedef std::vector< YAndPos > y_and_pos_type;

// namespace boost {
//   namespace numeric {
//     namespace odeint {

//       template<
// 	class State,
// 	class Value = double ,
// 	class Deriv = State ,
// 	class Time = Value ,
// 	class Algebra = range_algebra ,
// 	class Operations = default_operations ,
// 	class Resizer = initially_resizer
// 	>
//       class runge_kutta_with_diss
// 	: public explicit_error_stepper_fsal_base<
// 	runge_kutta_with_diss< State , Value , Deriv , Time , Algebra , Operations , Resizer > , 5 , 5 , 4 , State , Value , Deriv , Time , Algebra , Operations , Resizer >
//       {

//       public:
// 	typedef explicit_error_stepper_fsal_base<
// 	  runge_kutta_with_diss< State , Value , Deriv , Time , Algebra , Operations , Resizer > , 5 , 5 , 4 , State , Value , Deriv , Time , Algebra , Operations , Resizer > stepper_base_type;

// 	runge_kutta_dopri5< State , Value , Deriv , Time , Algebra , Operations , Resizer > runge_kutta_without_diss;

// 	typedef typename stepper_base_type::state_type state_type;
// 	typedef typename stepper_base_type::wrapped_state_type wrapped_state_type;
// 	typedef typename stepper_base_type::value_type value_type;
// 	typedef typename stepper_base_type::deriv_type deriv_type;
// 	typedef typename stepper_base_type::wrapped_deriv_type wrapped_deriv_type;
// 	typedef typename stepper_base_type::time_type time_type;
// 	typedef typename stepper_base_type::algebra_type algebra_type;
// 	typedef typename stepper_base_type::operations_type operations_type;
// 	typedef typename stepper_base_type::resizer_type resizer_type;
// 	typedef typename stepper_base_type::stepper_type stepper_type;

//    	runge_kutta_with_diss( const algebra_type &algebra = algebra_type() ) : stepper_base_type( algebra ), runge_kutta_without_diss(algebra)
// 	{ 
// 	}


// 	template< class System , class StateIn , class DerivIn , class StateOut , class DerivOut >
// 	void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt_in , time_type t ,
// 			   StateOut &out , DerivOut &dxdt_out , time_type dt )
// 	{
// 	  runge_kutta_without_diss.do_step_impl(system, in, dxdt_in, t, out, dxdt_out, dt);
// 	  system.removeEnergy(out, dxdt_out);
// 	}

// 	template< class System , class StateIn , class DerivIn , class StateOut , class DerivOut , class Err >
// 	void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt_in , time_type t ,
// 			   StateOut &out , DerivOut &dxdt_out , time_type dt , Err &xerr )
// 	{
// 	  runge_kutta_without_diss.do_step_impl(system, in, dxdt_in, t, out, dxdt_out, dt, xerr);
// 	  system.removeEnergy(out, dxdt_out);
// 	}

// 	template< class StateOut , class StateIn1 , class DerivIn1 , class StateIn2 , class DerivIn2 >
// 	void calc_state( time_type t , StateOut &x ,
// 			 const StateIn1 &x_old , const DerivIn1 &deriv_old , time_type t_old ,
// 			 const StateIn2 & x_new , const DerivIn2 &deriv_new , time_type t_new )
// 	{
// 	  runge_kutta_without_diss.calc_state(t, x, x_old, deriv_old, t_old, x_new, deriv_new, t_new);
// 	}

// 	template< class StateIn >
// 	void adjust_size( const StateIn &x )
// 	{
// 	  runge_kutta_without_diss(x);
// 	}
  
//       };
  
//       template< class State , class Value , class Deriv , class Time , class Algebra , class Operations , class Resize >
//       struct get_dense_output< runge_kutta_with_diss< State , Value , Deriv , Time , Algebra , Operations , Resize > >
//       {
// 	typedef runge_kutta_with_diss< State , Value , Deriv , Time , Algebra , Operations , Resize > stepper_type;
// 	typedef controlled_runge_kutta< stepper_type > controller_type;
// 	typedef dense_output_runge_kutta< controller_type > type;
//       };
//     } // odeint
//   } // numeric
// } // boost


// class Observer
// {

// public:


//   Observer( std::vector< state_type > &states , 
// 	    std::vector<double>& times, std::vector< double > &given_times, double dt)
//     : given_times_(given_times), states_(states), times_(times),
//       dt_(dt), diff_times_(given_times.size()) { 
//   }

//   void operator()( const state_type &x , double t );

// private:
//   const std::vector< double >& given_times_;
//   std::vector< state_type >& states_;
//   std::vector< double >& times_;
//   const double dt_;
//   std::vector< double > diff_times_;
// };

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
  std::vector<double> new_y;
  std::vector<double> new_U;
  std::vector<double> new_H;
  // std::vector<double> new_qm1;
  // std::vector<double> new_w;
  // std::vector<double> new_h;
  // std::vector<double> new_r;
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
  int is_relabeling;
  double collision_limiter;

  Simulation();
  void resetPointers(state_type& Z);
  void computeG();
  void fixYandH();
  void setTrack();
  void setDiss();
  void init();
  void initAntipeakon();
  void initPeakon();
  void initCubic();
  void initCubic2();
  void initZeroU();
  void operator()(const state_type& ZZ, state_type& dZdt, const double /* t */ );
  void solve();
  void computePandQ();
  // void concatenate(state_type& ZZ) const;
  // void deconcatenate(const state_type& ZZ);
  void writeToFile();
  void relabelingAndCopyAll();
  void removeEnergy();
  bool checkIfInCollisionRegion(double y1, double U1, double H1, double q1, double w1, double h1, double r1);
  bool checkIfCollided(double y1, double U1, double H1, double q1, double w1, double h1, double r1);
  
  // class convertBoolToDouble {
  // public:
  //   double operator()(bool x) {
  //     if (x == true) {
  // 	return 1.0;
  //     } else {
  // 	return 0.0;
  //     }
  //   }
  // };

  // class convertDoubleToBool {
  // public:
  //   bool operator()(double x) {
  //     if (x == 0.0) {
  // 	return false;
  //     } else if (x == 1.0) {
  // 	return true;
  //     } else {
  // 	std::cout<< "convertDoubleToBool failed!";
  // 	return true;
  //     }	
  //   }
  // };

};

#endif // CHSYS_HPP_INCLUDED

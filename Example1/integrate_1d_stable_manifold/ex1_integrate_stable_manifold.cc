/*
 * ex1_integrate_stable_manifold.cc
 * rigorous integration of (5.16) ((5.7) in backward time)
 *
 * written  November 20, 2019 Akitoshi Takayasu (with kv-0.4.48)
 *
 */
#include <kv/ode-maffine.hpp>
namespace ub = boost::numeric::ublas;
typedef kv::interval<double> intval;
FILE* fp;
//
// Desingularized vector field
template <class TT> struct VF {

  TT rho1, rho2, c, c1, c2;
  VF(TT rho1, TT rho2, TT c, TT c1, TT c2): rho1(rho1), rho2(rho2), c(c), c1(c1), c2(c2) {}

  template <class T> ub::vector<T> operator() (const ub::vector<T>& x, const T& t) {
    ub::vector<T> dx(2);
    T x1, x2;
    x1 = x(0); x2 = x(1);
    dx(0) = -(pow(x1,3) - (rho1+rho2)*pow(x1,2) + rho1*rho2*x1 - c*pow(x1,3)*x2 - c1*pow(x1,2)*x2);
    dx(1) = -(-0.5*pow(x1,2)*x2 + 0.5*rho1*rho2*x2 + c*pow(x1,2)*pow(x2,2) + c2*pow(x1,2)*pow(x2,3));
    return dx;
  }
};
// B1_func
struct B1_func {// T will be kv::interval<double>
  template <class T> T operator() (const T& beta, const T& rho1, const T& rho2) {
    return (beta-rho1)*(beta-rho2)/beta;
  }
};
// B2_func
struct B2_func {// T will be kv::interval<double>
  template <class T> T operator() (const T& beta, const T& rho1, const T& rho2) {
    return (pow(beta,2)-rho1*rho2)/(2*pow(beta,2));
  }
};
// ode_callback function:
// Integrate the passing time 't'
namespace kv{
template <class T> struct ode_callback_sample : ode_callback<T> {
    interval<T> &ptime;

    ode_callback_sample(interval<T>& ptime) : ptime(ptime) {}

    virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
      ptime += eval(integrate(result(1)/pow(result(0),2)),end-start);
      // Stopping criterion if it is needed...
      T data_out;
      data_out = x_e(0).lower();
      fwrite(&data_out, sizeof(T), 1, fp);
      data_out = x_e(0).upper();
      fwrite(&data_out, sizeof(T), 1, fp);
      data_out = mid(x_e(0));
      fwrite(&data_out, sizeof(T), 1, fp);
      data_out = x_e(1).lower();
      fwrite(&data_out, sizeof(T), 1, fp);
      data_out = x_e(1).upper();
      fwrite(&data_out, sizeof(T), 1, fp);
      data_out = mid(x_e(1));
      fwrite(&data_out, sizeof(T), 1, fp);
      data_out = ptime.lower();
      fwrite(&data_out, sizeof(T), 1, fp);
      data_out = ptime.upper();
      fwrite(&data_out, sizeof(T), 1, fp);
      data_out = mid(ptime);
      fwrite(&data_out, sizeof(T), 1, fp);
      return true;
    }
};
};
// Main function:
// 1. Rigorously integrate the trajectory on the local stable manifold in backward time
//
int main()
{
  const int n = 2;    // # of unknowns
  const intval rho1(1.), rho2(2.), beta_R(1.5), v_R(5.), beta_L("1.9"), v_L(4.);// parameters
  intval start, end, c, c1, c2, t_n;
  ub::vector< intval > x, p0;
  std::cout.precision(17);

  B1_func B1;
  B2_func B2;

  c = (v_R*B1(beta_R,rho1,rho2) - v_L*B1(beta_L,rho1,rho2))/(beta_R-beta_L);

  c1 = v_L*B1(beta_L,rho1,rho2) - c*beta_L;
  c2 = pow(v_L,2)*B2(beta_L,rho1,rho2) - c*v_L;

  VF<intval> f(rho1, rho2, c, c1, c2);

  // initial data on the local stable manifold by the parametrization method
  x.resize(n);
  p0.resize(n);// p0: initial data of the desingularized ODE


  // ======================= Data of the parametrization method =======================
  // Point on the local stable manifold
  p0(0) = intval("1.99704842870221","1.99704842870362");
  p0(1) = intval("0.06209042154030","0.06209042154164");
  // Blow up time from the above point (in the original time scale)
  t_n  = intval("0.01945344745624","0.01945344745758");
  // ==================================================================================

  for (int i = 0; i < n; i++) {
    x(i) = p0(i);
  }

  fp = fopen("data.bin", "wb"); // file output
  double data_out;
  data_out = p0(0).lower();
  fwrite(&data_out, sizeof(double), 1, fp);
  data_out = p0(0).upper();
  fwrite(&data_out, sizeof(double), 1, fp);
  data_out = mid(p0(0));
  fwrite(&data_out, sizeof(double), 1, fp);
  data_out = p0(1).lower();
  fwrite(&data_out, sizeof(double), 1, fp);
  data_out = p0(1).upper();
  fwrite(&data_out, sizeof(double), 1, fp);
  data_out = mid(p0(1));
  fwrite(&data_out, sizeof(double), 1, fp);
  data_out = t_n.lower();
  fwrite(&data_out, sizeof(double), 1, fp);
  data_out = t_n.upper();
  fwrite(&data_out, sizeof(double), 1, fp);
  data_out = mid(t_n);
  fwrite(&data_out, sizeof(double), 1, fp);

  // Rigorous integration of the vector field
  start = 0;
  end = 40;
  int r = kv::odelong_maffine(
    f,
    x,
    start,
    end,
    kv::ode_param<double>().set_verbose(0).set_order(12).set_restart_max(10).set_epsilon(pow(2,-68)),
    kv::ode_callback_sample<double>(t_n)
  );
  if (r == 0) {
      std::cout << "Cannot calculate solution.\n";
      return 1;
  } else if (r == 1) {
    std::cout << "Solution incomplitely calculated until t = " << end << ".\n";
    std::cout << x << "\n";
  } else if (r == 3) {
      std::cout << "Solution calculated until t = " << end << ".\n";
      std::cout << x << "\n";
  } else {
      std::cout << "Solution completely calculated until t = " << end << ".\n";
      std::cout << "x_n: " << x << "\n";
      std::cout << "t_xi0: " << t_n << std::endl;
      return 1;
  }
  // std::cout << "T_n: " << t_n << std::endl;
  fclose(fp);
}

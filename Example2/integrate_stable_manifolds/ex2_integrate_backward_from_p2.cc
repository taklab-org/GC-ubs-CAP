/*
 * ex2_integrate_backward_from_p2.cc
 *
 * Validation code of connecting orbit using integration backward time
 * from 2d stable manifold attached p1
 * to the equilibrium p_source of
 *
 * x0' = g0(x)
 * x1' = g1(x)
 * x2' = g2(x)
 *
 * which is the desingularized vector field of
 *
 * y0' = y0(y0^2-1)
 * y1' = y0^2y1 + y0^2y2
 * y2' = y0^2y2 + δ(cy0^2y2 - y1(y1-ay0)(y0-y1) + wy0^2)
 *
 * via the Poincare compactification.
 *
 *
 * written  Feb. 17, 2021 Akitoshi Takayasu (with kv-0.4.48)
 *
 */

#include <kv/ode-maffine.hpp>
#include <kv/kraw-approx.hpp>
#include <kv/eig.hpp>
// #include "eig.hpp"
const int size_of_ode = 3;    // Size of ODEs: fixed
const int k = 2; // quasi-homogeneous vector field of the order k+1: fixed
const int alph[] = {1,1,1}; // {alpha_i}: fixed
const int beta[] = {1,1,1}; //{beta_i}: fixed
const int gamm = 1; // c = alpha_1*beta_1 = ... = alpha_n*beta_n: fixed

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> intval;

ub::matrix< intval > Y;
ub::vector< intval > pb;

// f̃(x)
struct tf_func {
    template <class T> ub::vector<T> operator() (const ub::vector<T>& x) {
        int n = size_of_ode;
        ub::vector<T> tf(n);
        T x0 = x(0); T x1 = x(1); T x2 = x(2);
        // a = 0.3; c = 0.7; delta = 9; w = 0.02;
        static T a = kv::constants<T>::str("0.3");
        static T c = kv::constants<T>::str("0.7");
        static T delta = kv::constants<T>::str("9.0");
        static T w = kv::constants<T>::str("0.02");
        tf(0) = pow(x0,3) - (1. - pow(x0,2) - pow(x1,2) - pow(x2,2)) * x0;
        tf(1) = pow(x0,2) * x1 + pow(x0,2) * x2;
        tf(2) = pow(x0,2) * x2 + (1./delta) * (c * pow(x0,2) * x2 - x1 * (x1 - a * x0) * (x0 - x1) + w * pow(x0,3));
        return tf;
    }
};
// sum = x0.*tf0+x1.*tf1+x2.*tf2;
struct sum_func {
    template <class T> T operator() (const ub::vector<T>& x) {
        tf_func tf;
        return ub::inner_prod(x,tf(x));
    }
};
// Generator of desingularized vector field (backward direction)
struct Func {
    template <class T> ub::vector<T> operator() (const ub::vector<T>& x, const T& t) {
        int n = size_of_ode;
        ub::vector<T> y(n);
        ub::vector<T> tf(n);
        tf_func f;
        sum_func sum;
        tf = f(x); T sumx = sum(x);
        T x0 = x(0); T x1 = x(1); T x2 = x(2);
        y(0) = tf(0) - x0 * sumx;
        y(1) = tf(1) - x1 * sumx;
        y(2) = tf(2) - x2 * sumx;
        return -y;
    }
};
// Func for kraw-approx
struct Func_krawapp {
    template <class T> ub::vector<T> operator() (const ub::vector<T>& x) {
        Func f;
        return f(x,(T) 1.);
    }
};
// ode_callback function:
// Check the validated solution is in the Lyapunov domain.
namespace kv{
template <class T> struct ode_callback_sample : ode_callback<T> {
    T eps;

    ode_callback_sample(T eps) : eps(eps) {}

    virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
        ub::vector< interval<T> > xp;
        kv::interval<T> Lx;
        // s = result.size();
        // psa< interval<T> > px2c; // 1-p(x)^{2c}
        // px2c = 1.;
        // for (i = 0; i<s; i++) {
        //   px2c -= pow(result(i),2*beta[i]);
        // }
        // if (k==1) {
        //   result2 += eval(integrate((1.0-0.5*(2*gamm-1)/gamm*px2c)*px2c),end-start);
        // } else {
        //   result2 += eval(integrate((1.0-0.5*(2*gamm-1)/gamm*px2c)*pow(px2c,k)),end-start);
        // }
        // Stopping criterion
        xp = x_e-pb;
        Lx = inner_prod(prod(Y,xp),xp);
        if (Lx < pow(eps,(kv::interval<T>)2)) return false;
        return true;
    }
};
};
//
int main()
{
    int n = size_of_ode;    // # of unknowns
    ub::vector< intval > x, y0, x0;
    intval start, end, c1, cn_tmp, cn, Tmax, t_n, y, Lx, kappa;
    std::cout.precision(17);

    // Initial value: x(0)
    // x0(0) = -0.1; x0(1) = -0.1;
    // // From p0
    // x0(0) = intval("0.876149668192013","0.876149668193990");
    // x0(1) = intval("0.305009506848445","0.305009506850414");
    // x0(2) = intval("-0.000441867191975","-0.000441867190016");

    // From p1_1
    // x0(0) = intval("0.74051007585989","0.74051009207412");
    // x0(1) = intval("0.56646602203592","0.56646603825015");
    // x0(2) = intval("0.01328520714802","0.01328522336223");

    // // From p1_2
    // x0(0) = intval("0.74051007585989","0.74051009207412");
    // x0(1) = intval("0.56646602203592","0.56646603825015");
    // x0(2) = intval("0.01328520714802","0.01328522336223");


    // Validation of the blow-up solution
    // 1. Get the equilibrium point (originally source equilibrium)
    pb.resize(n);
    pb(0) = 0.7071051816183367;
    pb(1) = 0.001504037399468;
    pb(2) = -0.001504037399468;
    kv::krawczyk_approx(Func_krawapp(), mid(pb), pb, 15, 0);
    std::cout << "p_source: " << pb << "\n";

    // 2. Search the maximum range of the Lyapunov domain:
    // Calculate C1
    Func f;
    ub::vector< intval > y1,y1_tmp,y2,radius;
    ub::matrix< intval > Dfx, Ax;
    ub::vector< kv::autodif< intval > > dy1, dy2;
    // Initial vector of auto diff.
    y1.resize(n);
    for (int i = 0; i < n; i++) {
      // y1(i) = intval::hull(x(i),p(i));
      y1(i) = pb(i);
    }
    dy1 = kv::autodif< intval >::init(y1);
    dy2 = f(dy1,(kv::autodif< intval >) 0); // f(dy1) = [f(y1), df/dx[y1](Jacobian) (=Dfx)]
    kv::autodif< intval >::split(dy2, y2, Dfx);
    // std::cout << "Df[p]: " << mid(Dfx) << std::endl;
    // Compute eigenvalue of Df(p)
    // Y=Re(V^-H*V^-1) s.t. V^-1*Dfx*V = D (diagonal matrix)
    ub::vector< kv::complex< intval > > lambda;
    ub::matrix< kv::complex<double> > V, D, R, RH, Y_comp, I;
    kv::eig(mid(Dfx),V,D);
    std::cout << "approximate eigenvalue:" << '\n';
    for (int i = 0; i < n; i++) {
      std::cout << D(i,i) << '\n';
    }
    // std::cout << "D: " << D << std::endl;
    kv::invert_comp(V,R);
    // std::cout << "R: " << R << std::endl;
    RH.resize(n,n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        RH(j,i) = kv::complex<double>(R(i,j).real(),-R(i,j).imag());
      }
    }
    // std::cout << "R^H: " << RH << std::endl;
    I.resize(n,n);
    for (int i = 0; i < n; i++) {
      I(i,i) = 1;
    }
    // std::cout << "I: " << I << std::endl;
    Y_comp = prod(RH,I);
    Y_comp = prod(Y_comp,R);
    // std::cout << "Y_comp: " << Y_comp << std::endl;
    Y.resize(n,n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        Y(i,j) = Y_comp(i,j).real();
      }
    }
    // std::cout << "Y: " << Y << std::endl;
    kv::veig(Y,lambda);
    // std::cout << "lambda_Y: " << lambda <<  std::endl;
    using std::abs;
    c1 = abs(lambda(0).real());
    for (int i = 1; i < n; i++) {
      if (c1 > abs(lambda(i).real())) {
        c1 = abs(lambda(i).real());
      }
    }
    c1 = 1./c1;
    std::cout << "c_1: " << c1 << std::endl;

    // Calculate Cn by the delta inflation scheme
    double eps_tmp = 1e-10, eps;// Radius of N
    y1_tmp.resize(n);
    radius.resize(n);
    while (true) {
      //------------------------- A(x) by using autodiff -----------------------
      for (int i = 0; i < n; i++) {
        radius(i) = intval(-eps_tmp,eps_tmp);
      }
      y1_tmp = pb + radius*sqrt(c1);
      // std::cout << "y1 = " << y1 << std::endl;
      dy1 = kv::autodif< intval >::init(y1_tmp);
      dy2 = f(dy1,(kv::autodif< intval >) 0); // f(dy1) = [f(y1), df/dx[y1](Jacobian) (=Dfx)]
      kv::autodif< intval >::split(dy2, y2, Dfx);
      Ax = prod(trans(Dfx),Y)+prod(Y,Dfx);
      // std::cout << "Ax: " << Ax << std::endl;
      //------------------------- A(x) by using autodiff -----------------------
      kv::veig(Ax,lambda);
      // std::cout << "lambda_AX: " << lambda << std::endl;
      using std::abs;
      cn_tmp = abs(lambda(0));
      for (int i = 1; i < n; i++) {
        if (cn_tmp > abs(lambda(i))) {
          cn_tmp = abs(lambda(i));
        }
      }
      // std::cout << "c_n: " << cn << std::endl;
      // if (in(0,cn_tmp)) {
      if (in(0,cn_tmp)||rad(cn_tmp)>1e-2) { // Uncomment if you want the tight "cn".
        if (eps_tmp==1e-10) {
          std::cout << "Failed in constructing the Lyapunov function" << '\n';
          return 1;
        }
        break;
      }
      y1 = y1_tmp;
      cn = cn_tmp;
      eps = eps_tmp;
      eps_tmp *= 1.1;
    }
    std::cout << "The Lyapunov domain is validated." << '\n';
    std::cout << "c_n: " << cn << std::endl;
    std::cout << "eps = " << eps << '\n';


    // 3. Validate that x(tau_N) is in the Lyapunov domain
    // t_n  = 0; // Initialize the blow-up time
    // Read initial data from the results of paraterization method
    x0.resize(n);
    int num = 54, num_points = num/n/2;
    double pp[num];
    FILE *fp;
    fp = fopen("initial_points_p2.bin","r");
    fread(pp, sizeof(double), num, fp);
    fclose(fp);
    for (int points = 0; points < num_points; points++) {
        std::cout << "Point: " << points+1 << "\n";
        x0(0) = intval(pp[2*n*points],pp[2*n*points+1]);
        x0(1) = intval(pp[2*n*points+2],pp[2*n*points+3]);
        x0(2) = intval(pp[2*n*points+4],pp[2*n*points+5]);
        // std::cout << "x0" << x0 << '\n';
    x = x0;
    start = 0;
    end = 400;
    // std::cout << "t: " << t_n << std::endl;
    // std::cout << x << std::endl;
    int r = kv::odelong_maffine(
      Func(),
      x,
      start,
      end,
      // (intval)eps,
      // p,
      // Y,
      kv::ode_param<double>().set_verbose(0).set_order(24).set_restart_max(10),//.set_epsilon(pow(2,-68)),
      kv::ode_callback_sample<double>(eps)
    );
    // int r = kv::odelong_maffine0(Func(), x, start, end, kv::ode_param<double>().set_verbose(0).set_order(24).set_restart_max(1),kv::ode_callback_sample<double>(t_n));
    if (r == 0) {
        std::cout << "Cannot calculate solution.\n";
        return 1;
    } else if (r == 1 || r == 3) {
        std::cout << "Lyapunov validation is completed by calculating solution at t = " << end << ".\n";
        std::cout << "The end point of stable manifold is \n" << x << "\n";
    } else {
        // std::cout << "Solution completely calculated until t = " << end << ".\n";
        // std::cout << "x_n: " << x << "\n";
        std::cout << "Cannot valitate that x(tau_N) is in the Lyapunov domain." << '\n';
        return 1;
    }
    }
}

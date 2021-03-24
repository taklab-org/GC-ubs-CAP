/*
 * ex3_integrate_forward_with_blowup_time.cc
 *
 * Validation code of blowups from some intial points crossed over the stable manifolds attached to p_{inf,s}
 *
 * x1' = (x1.^2 - x2).*F(x) - x1.*G(x)
 * x2' = (x1.^3/3 - x1.*(1-px4(x)).^2).*F(x) - 2*x2.*G(x),
 *
 * which is the desingularized vector field of
 *
 * u' = u^2 - v
 * v' = u^3/3 - u
 *
 * via the quasi-parabolic compactification.
 *
 *
 * written  Oct. 14, 2020 Akitoshi Takayasu (with kv-0.4.48)
 *
 */

#include <kv/ode-maffine.hpp>
#include <kv/kraw-approx.hpp>
#include <kv/eig.hpp>
// #include "eig.hpp"
const int size_of_ode = 2;    // Size of ODEs: fixed
const int k = 1; // quasi-homogeneous vector field of the order k+1: fixed
const int alph[] = {1,2}; // {alpha_i}: fixed
const int beta[] = {2,1}; //{beta_i}: fixed
const int gamm = 2; // c = alpha_1*beta_1 = ... = alpha_n*beta_n: fixed

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> intval;
FILE* fp;

ub::matrix< intval > Y;
ub::vector< intval > pb;

// p(x)^4 = x1^4 + x2^2;
struct p_func {
    template <class T> T operator() (const ub::vector<T>& x) {
        return pow(x(0),4) + pow(x(1),2);
    }
};
// F(x) = 0.25*(1. + 3*px4(x));
struct F_func {
    template <class T> T operator() (const ub::vector<T>& x) {
        p_func p;
        return 0.25*(1.0 + 3.0*p(x));
    }
};
// G(x) = x1.^3.*(x1.^2-x2) + 0.5*x2.*(x1.^3/3 - x1.*(1-px4(x)).^2);
struct G_func {
    template <class T> T operator() (const ub::vector<T>& x) {
        T x1 = x(0); T x2 = x(1);
        p_func p;
        return pow(x1,3)*(pow(x1,2) - x2) + 0.5*x2*(pow(x1,3)/3. - x1*pow(1 - p(x),2));
    }
};
// Generator of desingularized vector field (backward direction)
struct Func {
    template <class T> ub::vector<T> operator() (const ub::vector<T>& x, const T& t) {
        int n = size_of_ode;
        ub::vector<T> y(n);
        F_func F;
        G_func G;
        p_func p;
        T x1 = x(0); T x2 = x(1);
        y(0) = (pow(x1,2) - x2)*F(x) - x1*G(x);
        y(1) = (pow(x1,3)/3. - x1*pow(1 - p(x),2))*F(x) - 2.*x2*G(x);
        return y;
    }
};
// Func for kraw-approx (vector field on the "horizon")
struct Func_krawapp {
    template <class T> ub::vector<T> operator() (const ub::vector<T>& x) {
        // Func f;
        // return f(x,(T) 1.);
        int n = size_of_ode;
        ub::vector<T> y(n);
        T x1 = x(0); T x2 = x(1);
        T G = pow(x1,5) - 5.*pow(x1,3)*x2;
        y(0) = (pow(x1,2) - x2) - x1*G;
        y(1) = pow(x1,3)/3. - 2.*x2*G;
        return y;
    }
};
// ode_callback function:
// Check the validated solution is in the Lyapunov domain.
namespace kv{
template <class T> struct ode_callback_sample : ode_callback<T> {
    interval<T> &result2;
    T eps;

    ode_callback_sample(interval<T>& result2, T eps) : result2(result2), eps(eps) {}

    virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
        int i, s;
        ub::vector< interval<T> > xp;
        kv::interval<T> Lx;
        s = result.size();
        psa< interval<T> > px2c; // 1-p(x)^{2c}
        px2c = 1.;
        for (i = 0; i<s; i++) {
          px2c -= pow(result(i),2*beta[i]);
        }
        if (k==1) {
          result2 += eval(integrate((1.0-0.5*(2*gamm-1)/gamm*px2c)*px2c),end-start);
        } else {
          result2 += eval(integrate((1.0-0.5*(2*gamm-1)/gamm*px2c)*pow(px2c,k)),end-start);
        }
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
    ub::vector< intval > x, y0, x0, xp;
    intval start, end, c1, cn_tmp, cn, Tmax, t_n, y, Lx, kappa;
    std::cout.precision(17);
    // Initial value: x(0)
    x0.resize(n);
    ub::vector<double> ps, p_infs, v, v_normal;
    ps.resize(n); p_infs.resize(n); v.resize(n); v_normal.resize(n);
    // x0(0) = 0.1; x0(1) = 0.1;
    // From p0
    // 0.320682556228925   0.319896988077054; mid point
    // 0.320682556228873   0.319896988077002; inf
    // 0.320682556228977   0.319896988077106; sup
    // x0(0) = intval("0.320682556228873","0.320682556228977");
    // x0(1) = intval("0.319896988077002","0.319896988077106");
    // From p_inf_s
    p_infs(0) = 0.8861081289780320;
    p_infs(1) = 0.6192579489210105;
    // 0.836474250605419   0.572882744753758; mid point
    // 0.836474250560393   0.572882744708793; inf
    // 0.836474250650447   0.572882744798721; sup
    // x0(0) = intval("0.836474250560393","0.836474250650447");
    // x0(1) = intval("0.572882744708793","0.572882744798721");
    ps(0) = 0.832267964033316;
    ps(1) = 0.570370294663353;
    for (int i = 0; i < n; i++) {
        v(i) = p_infs(i) - ps(i);
    }
    v_normal(0) = v(1);
    v_normal(1) = -v(0);
    v_normal /= norm_2(v_normal);
    // std::cout << "v_normal:" << v_normal << '\n';
    // // For check
    // kappa = 1.;
    // for (int i = 0; i < n; i++) {
    //   kappa -= pow(x0(i),2*beta[i]);
    // }
    // y0: initial vector of y(t): you must edit for each ODE.
    // y0.resize(n);
    // for (int i = 0; i < n; i++) {
    //   y0(i) = x0(i)/pow(kappa,alph[i]);
    // }
    // std::cout << "y0: " << y0 << '\n';

    // Validation of the blow-up solution
    // 1. Get the equilibrium point (originally source equilibrium)
    pb.resize(n);
    pb(0) = 0.989136995894978;
    pb(1) = 0.206758557005181;
    kv::krawczyk_approx(Func_krawapp(), mid(pb), pb, 15, 0);
    std::cout << "pb: " << pb << "\n";
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
    double p_step = 1e-12, param = 1e-8;
    p_func px4;
    int count = 0;
    fp = fopen("data_more.bin", "wb"); // file output
    double data_out;
    while (true) {
        t_n  = 0; // Initialize the blow-up time
        x0 = ps + param*v_normal;
        if (px4(x0)>1 || param <= 1e-12) break;
        //---- just for figures ---
        data_out = param;
        fwrite(&data_out, sizeof(double), 1, fp);
        //-----------
        param -= p_step;
        count++;
        // std::cout << "count: " << count << '\n';
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
            kv::ode_callback_sample<double>(t_n,eps)
        );
        // int r = kv::odelong_maffine0(Func(), x, start, end, kv::ode_param<double>().set_verbose(0).set_order(24).set_restart_max(1),kv::ode_callback_sample<double>(t_n));
        if (r == 0) {
            std::cout << "Cannot calculate solution.\n";
            return 1;
        } else if (r == 1 || r == 3) {
            std::cout << "Lyapunov validation is completed by calculating solution at t = " << end << ".\n";
            std::cout << "The solution inside the Lyapunov domain is \n" << x << "\n";
        } else {
            // std::cout << "Solution completely calculated until t = " << end << ".\n";
            // std::cout << "x_n: " << x << "\n";
            std::cout << "Cannot valitate that x(tau_N) is in the Lyapunov domain." << '\n';
            return 1;
        }
        std::cout << "T_n: " << t_n << std::endl;
        xp = x-pb;
        Lx = inner_prod(prod(Y,xp),xp);
        std::cout << "L(x): " << Lx << std::endl;

        // 4. Compute the validated blow-up time of the ODE
        //----------------- This part is specialized for 1-p(x)^4. -----------------
        // kappa = max{1, 6|x^*_2|^2}
        if (1.<6*pow(abs(pb(0)),2)) {
            kappa = 6*pow(abs(pb(0)),2);
        } else {
            kappa = 1.;
        }
        std::cout << "T_max-T_n: " << mag((
            4*pow((pow(pb(1),2)+4*pow(pb(0),6))*Lx/c1,0.5)
            +kappa*Lx
            +8./3*abs(pb(0))*pow(c1,0.5)*pow(Lx,1.5)
            +0.5*c1*pow(Lx,2))/cn) << std::endl;
        Tmax = intval::hull(t_n,t_n+mag((
            4*pow((pow(pb(1),2)+4*pow(pb(0),6))*Lx/c1,0.5)
            +kappa*Lx
            +8./3*abs(pb(0))*pow(c1,0.5)*pow(Lx,1.5)
            +0.5*c1*pow(Lx,2))/cn));
        //--------------------------------------------------------------------------
        std::cout << "Tmax: " << Tmax << std::endl;
        //---- just for figures ---
        data_out = mid(Tmax);
        fwrite(&data_out, sizeof(double), 1, fp);
        data_out = rad(Tmax);
        fwrite(&data_out, sizeof(double), 1, fp);
        //-----------
    }
    fclose(fp);
}

/*
 * Implementation of the Kamiya algorithm using the CERES non linear solver.
 *
 *
 */


#include "ceres/ceres.h"
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
#include "CLI11.hpp"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"




// A templated cost functor that implements the residual r = 10 -
// x. The method operator() is templated so that we can then use an
// automatic differentiation wrapper around it to generate its
// derivatives.
struct CostFunctor {
  template <typename T>
  bool operator()(const T* const x, T* residual) const {
    residual[0] = 10.0 - x[0];
    return true;
  }
};


struct F28a {
  double deltap1;
  double deltap2;
  double  k;
  double f0;
  double f1;
  double f2;
  double l0;
  double l1;
  double l2;
  
  template <typename T>
  bool operator()(const T* const x1, const T* const x2, T* residual) const {
    residual[0] =  (deltap1/k)*x1[0]*x1[0]*pow((f0*(((x1[0]*x1[0]*x1[0])/f1)+(x2[0]*x2[0]*x2[0])/f2)), 2.0/3.0)
                 - (f0*l0*x1[0]*x1[0]) - f1 * l1 * pow(f0*((x1[0]*x1[0]*x1[0]/f1)+(x2[0]*x2[0]*x2[0])/f2), 2.0/3.0);
    residual[1] =  (deltap2/k)*x2[0]*x2[0]*pow((f0*(((x1[0]*x1[0]*x1[0])/f1)+(x2[0]*x2[0]*x2[0])/f2)), 2.0/3.0)
                 - (f0*l0*x2[0]*x2[0]) - f2 * l2 * pow(f0*((x1[0]*x1[0]*x1[0]/f1)+(x2[0]*x2[0]*x2[0])/f2), 2.0/3.0);
    
    return true;
  }
};



int main(int argc, char** argv) {
  // parse command line CLI ----------------------------------------------
   CLI::App app;
  double x0 {0};
  double y0 {0};
  
  double x1 {0};
  double y1 {0};

  double x2 {0};
  double y2 {0};

  double xb {0};
  double yb {0};
  double k {3};
  double r_ori {1.0};
  
  
  
  app.description("Test the Kamiya algorithm\n");
  app.add_option("--x0", x0, "x0", true );
  app.add_option("--y0", y0, "y0", true );
  
  app.add_option("--x1", x1, "x1", true );
  app.add_option("--y1", y1, "y1", true );
  
  app.add_option("--x2", x2, "x2", true );
  app.add_option("--y2", y2, "y2", true );


  
  app.add_option("-k", k, "k ", true );
  app.add_option("-r", r_ori, "r_ori", true );

  

   app.get_formatter()->column_width(40);
   CLI11_PARSE(app, argc, argv);
   // END parse command line using CLI ----------------------------------------------

  
  
  // Prepare
  DGtal::Z2i::RealPoint p0 (x0, y0);
  DGtal::Z2i::RealPoint p1 (x1, y1);
  DGtal::Z2i::RealPoint p2 (x2, y2);
 
  double r0 = r_ori;
  double r1 = r_ori;
  double r2 = r_ori;

  double f0 = k * r0*r0*r0;
  double f1 = k * r1*r1*r1;
  double f2 = k * r2*r2*r2;

  DGtal::Z2i::RealPoint pb ((f0*x0+f1*x1+f2*x2)/(2.0*f0), (f0*y0+f1*y1+f2*y2)/(2.0*f0));
  double l0 = (p0 - pb).norm();
  double l1 = (p1 - pb).norm();
  double l2 = (p2 - pb).norm();
  
  double deltaP1 = k*((f0*l0)/(r0*r0*r0*r0)+(f1*l1)/(r1*r1*r1*r1));
  double deltaP2 = k*((f0*l0)/(r0*r0*r0*r0)+(f2*l2)/(r2*r2*r2*r2));

  F28a *f = new F28a();
  f->deltap1 = deltaP1;
  f->deltap2 = deltaP2;
  f->k = k;
  f->f0 = f0;
  f->f1 = f1;
  f->f2 = f2;
  f->l0 = l0;
  f->l1 = l1;
  f->l2 = l2;

  
  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  
  double xx1 = 40;
  double xx2 = 40;

 // const double initial_x = x;
  // Build the problem.
  Problem problem;
  // Set up the only cost function (also known as residual). This uses
  // auto-differentiation to obtain the derivative (jacobian).
  CostFunction* cost_function =
      new AutoDiffCostFunction<F28a, 2, 1, 1>(f);
  problem.AddResidualBlock(cost_function, nullptr, &xx1, &xx2);
  // Run the solver!
  Solver::Options options;
  options.minimizer_progress_to_stdout = true;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  std::cout << summary.BriefReport() << "\n";
  
  f1 = k * xx1*xx1*xx1;
  f2 = k * xx2*xx2*xx2;
  r0 = pow((xx1*xx1*xx1 + xx2*xx2*xx2), 1.0/3.0);
  f0 = k * r0*r0*r0;
  
  DGtal::Z2i::RealPoint pbNew ((f0*x0+f1*x1+f2*x2)/(2.0*f0), (f0*y0+f1*y1+f2*y2)/(2.0*f0));

  std::cout << "xx1 : " << xx1 << " and xx2 " << xx2 << " r0 " << r0 << "\n";
  std::cout << "pbNew[0]  : " << pbNew[0] << " and pbNew[1] " << pbNew[1] << "\n";
  return 0;
}

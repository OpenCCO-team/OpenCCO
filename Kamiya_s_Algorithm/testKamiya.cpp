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

#include "DGtal/io/boards/Board2D.h"

using namespace DGtal;
using namespace DGtal::Z2i;

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



static void kamiyaOpt(double deltaP1, double deltaP2, double f0, double f1, double f2, double k, double l0, double l1, double l2, double &xx1, double &xx2) {
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
  options.minimizer_progress_to_stdout = false;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  //std::cout << summary.BriefReport() << "\n";
}

int main(int argc, char** argv) {
  // parse command line CLI ----------------------------------------------
  CLI::App app;
  double x0 {0};
  double y0 {0};
  
  double x1 {0};
  double y1 {0};
  
  double x2 {0};
  double y2 {0};
  
  double gamma {3.0};
  double k {8*3.6/M_PI};
  double r_ori {1.0};
  std::string file_plot {"res.gnuplot"};
  
  app.description("Test the Kamiya algorithm\n");
  app.add_option("--x0", x0, "x0", true );
  app.add_option("--y0", y0, "y0", true );
  
  app.add_option("--x1", x1, "x1", true );
  app.add_option("--y1", y1, "y1", true );
  
  app.add_option("--x2", x2, "x2", true );
  app.add_option("--y2", y2, "y2", true );
  
  app.add_option("-k", k, "k ", true );
  app.add_option("-r", r_ori, "r_ori", true );
  app.add_option("-g", gamma, "gamma", true );
  
  app.add_option("--genGnuPlot", file_plot, "generate a gnuplot to display the result.", true );
  
  
  
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
  double f1 = 0.5*f0;//k * r1*r1*r1;
  double f2 = 0.5*f0;//k * r2*r2*r2;
  // Starting position from Equation (21) for initialisation as mentionned page 11 [Clara Jaquet et HT]
  //DGtal::Z2i::RealPoint pb ((f0*x0+f1*x1+f2*x2)/(f0+f1+f2), (f0*y0+f1*y1+f2*y2)/(f0+f1+f2));
  DGtal::Z2i::RealPoint pb ((f0*x0+f1*x1+f2*x2)/(2.0*f0), (f0*y0+f1*y1+f2*y2)/(2.0*f0));
  std::cout<<"Init:"<<pb<<std::endl;
  double l0 = (p0 - pb).norm();
  double l1 = (p1 - pb).norm();
  double l2 = (p2 - pb).norm();
  std::cout<<"l0="<<l0<<", l1="<<l1<<", l2="<<l2<<std::endl;
  double deltaP1 = k*((f0*l0)/(r0*r0)+(f1*l1)/(r1*r1));
  double deltaP2 = k*((f0*l0)/(r0*r0)+(f2*l2)/(r2*r2));
  
  double rr1 = r_ori;
  double rr2 = r_ori;
  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  DGtal::Z2i::RealPoint pbInit ((f0*x0+f1*x1+f2*x2)/(2.0*f0), (f0*y0+f1*y1+f2*y2)/(2.0*f0));
  DGtal::trace.info() << pbInit << std::endl;
  
  for (int i = 0; i< 100; i++){
    DGtal::trace.progressBar(i, 100);
    l0 = (p0 - pb).norm();
    l1 = (p1 - pb).norm();
    l2 = (p2 - pb).norm();
    
    kamiyaOpt(deltaP1, deltaP2, f0, f1, f2, k, l0, l1, l2, rr1, rr2);
    ///f1 = k * rr1*rr1*rr1;
    ///f2 = k * rr2*rr2*rr2;
    //r0 = pow((pow(rr1, gamma) + pow(rr2, gamma)), 1.0/gamma);
    //Equation 27
    r0 = pow(f0*(pow(rr1, gamma)/f1 + pow(rr2, gamma)/f2), 1.0/gamma);
    ///f0 = k * r0*r0*r0;
    // Equation (26) page 13
    //DGtal::Z2i::RealPoint pbNew ((f0*x0+f1*x1+f2*x2)/(2.0*f0), (f0*y0+f1*y1+f2*y2)/(2.0*f0));
    ///pb[0] = (x0*r0*r0/l0 + x1*r1*r1/l1 + x2*r2*r2/l2)/(r0*r0/l0+r1*r1/l1+r2*r2/l2);
    ///pb[1] = (y0*r0*r0/l0 + y1*r1*r1/l1 + y2*r2*r2/l2)/(r0*r0/l0+r1*r1/l1+r2*r2/l2);
    pb[0] = (x0*r0/l0 + x1*r1/l1 + x2*r2/l2)/(r0/l0+r1/l1+r2/l2);
    pb[1] = (y0*r0/l0 + y1*r1/l1 + y2*r2/l2)/(r0/l0+r1/l1+r2/l2);
    deltaP1 = k*((f0*l0)/(r0*r0)+(f1*l1)/(r1*r1));
    deltaP2 = k*((f0*l0)/(r0*r0)+(f2*l2)/(r2*r2));
    std::cout << "xx1 : " << rr1 << " and xx2 " << rr2 << " r0 " << r0 << "\n";
    std::cout << "pbNew[0]  : " << pb[0] << " and pbNew[1] " << pb[1] << "\n";
  }
  
  f1 = k * rr1*rr1*rr1;
  f2 = k * rr2*rr2*rr2;
  r0 = pow(f0*(pow(rr1, gamma)/f1 + pow(rr2, gamma)/f2), 1.0/gamma);

//  r0 = pow(pow(rr1, gamma) + pow(rr2, gamma), 1.0/gamma);
  f0 = k * r0*r0*r0;
  
  // DGtal::Z2i::RealPoint pbNew ((f0*x0+f1*x1+f2*x2)/(2.0*f0), (f0*y0+f1*y1+f2*y2)/(2.0*f0));
  
  std::cout << "xx1 : " << rr1 << " and xx2 " << rr2 << " r0 " << r0 << "\n";
  std::cout << "pbNew[0]  : " << pb[0] << " and pbNew[1] " << pb[1] << "\n";
  Board2D board;
  // Draw initial points
  board.setFillColor(DGtal::Color::Red);
  board.drawCircle(p0[0], p0[1], r_ori, 1);
  board.drawCircle(p1[0], p1[1], r_ori, 1);
  board.drawCircle(p2[0], p2[1], r_ori, 1);
  board.setFillColor(DGtal::Color::Yellow);
  board.drawCircle(pbInit[0], pbInit[1], r_ori, 1);
  // Draw initial lines
  board.setPenColor(DGtal::Color(200, 200, 0, 100));
  board.setLineWidth(57.5*r_ori);
  board.drawLine(p0[0], p0[1], pbInit[0], pbInit[1],2);
  board.drawLine(p1[0], p1[1], pbInit[0], pbInit[1],2);
  board.drawLine(p2[0], p2[1], pbInit[0], pbInit[1],2);
  // Draw results
  board.setLineWidth(0);
  board.setFillColor(DGtal::Color::Green);
  board.drawCircle(pb[0], pb[1], r0, 0);
  // Draw initial lines
  board.setPenColor(DGtal::Color(20, 200, 200, 100));
  board.setLineWidth(57.5*r0);
  board.drawLine(p0[0], p0[1], pb[0], pb[1],2);
  board.setLineWidth(57.5*sqrt(rr1));
  board.drawLine(p1[0], p1[1], pb[0], pb[1],2);
  board.setLineWidth(57.5*sqrt(rr2  ));
  board.drawLine(p2[0], p2[1], pb[0], pb[1],2);

  
  
  board.saveEPS("logoDGtal.eps");
  board.saveSVG("logoDGtal.svg");

  
  return 0;
}

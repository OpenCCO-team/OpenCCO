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

#include "geomhelpers.h"

#include "DGtal/io/boards/Board2D.h"

using namespace DGtal;
using namespace DGtal::Z2i;

#include "CLI11.hpp"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"


int main(int argc, char** argv) {
  // parse command line CLI ----------------------------------------------
  CLI::App app;
  double x0 {0};
  double y0 {0};
  
  double x1 {100};
  double y1 {0};
  
  double x2 {0};
  double y2 {100};
  
  double gamma {3.0};
  double r_ori {1.0};
  
  double ratioQ {0.5};
  
  std::string file_eps {"resultFig.eps"};
  std::string file_svg {"resultFig.svg"};

  
  
  app.description("Test the Kamiya algorithm\n");
  app.add_option("--x0", x0, "x0", true );
  app.add_option("--y0", y0, "y0", true );
  
  app.add_option("--x1", x1, "x1", true );
  app.add_option("--y1", y1, "y1", true );
  
  app.add_option("--x2", x2, "x2", true );
  app.add_option("--y2", y2, "y2", true );
  
  app.add_option("-r", r_ori, "r_ori", true );
  app.add_option("-g", gamma, "gamma", true );
  app.add_option("-R", ratioQ, "ratio", true );

  
  app.add_option("--genEpsDisplay,-e", file_eps, "generate an eps display of the result.", true );
  app.add_option("--genSvgDisplay,-s", file_svg, "generate an svg display of the result.", true );

  
  
  
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
  
  double f0 =  r0*r0*r0;
  double f1 = ratioQ*f0;//k * r1*r1*r1;
  double f2 = (1.0-ratioQ)*f0;//k * r2*r2*r2;
  // Starting position from Equation (21) for initialisation as mentionned page 11 [Clara Jaquet et HT]
  //DGtal::Z2i::RealPoint pb ((f0*x0+f1*x1+f2*x2)/(f0+f1+f2), (f0*y0+f1*y1+f2*y2)/(f0+f1+f2));
  DGtal::Z2i::RealPoint pb ((f0*x0+f1*x1+f2*x2)/(2.0*f0), (f0*y0+f1*y1+f2*y2)/(2.0*f0));
  std::cout<<"Init:"<<pb<<std::endl;
  double l0 = (p0 - pb).norm();
  double l1 = (p1 - pb).norm();
  double l2 = (p2 - pb).norm();
  std::cout<<"l0="<<l0<<", l1="<<l1<<", l2="<<l2<<std::endl;
  double deltaP1 = (f0*l0)/(r0*r0)+(f1*l1)/(r1*r1);
  double deltaP2 = (f0*l0)/(r0*r0)+(f2*l2)/(r2*r2);
  
  double rr1 = r_ori*r_ori;
  double rr2 = r_ori*r_ori;
  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  DGtal::Z2i::RealPoint pbInit ((f0*x0+f1*x1+f2*x2)/(2.0*f0), (f0*y0+f1*y1+f2*y2)/(2.0*f0));
  DGtal::trace.info() << pbInit << std::endl;
  double R0, R1, R2;
  for (int i = 0; i< 100; i++){
    DGtal::trace.progressBar(i, 100);
    kamiyaOpt(gamma, deltaP1, deltaP2, f0, f1, f2, l0, l1, l2, rr1, rr2);
    ///f1 = k * rr1*rr1*rr1;
    ///f2 = k * rr2*rr2*rr2;
    //r0 = pow((pow(rr1, gamma) + pow(rr2, gamma)), 1.0/gamma);
    //Equation 27
    R0 = pow(f0*(pow(rr1, gamma)/f1 + pow(rr2, gamma)/f2), 1.0/gamma);
    R1 = rr1;
    R2 = rr2;
    // Equation (26) page 13
    pb[0] = (x0*R0/l0 + x1*R1/l1 + x2*R2/l2)/(R0/l0+R1/l1+R2/l2);
    pb[1] = (x0*R0/l0 + y1*R1/l1 + y2*R2/l2)/(R0/l0+R1/l1+R2/l2);
  
    //Update values
    deltaP1 = (f0*l0)/(R0*R0)+(f1*l1)/(R1*R1);
    deltaP2 = (f0*l0)/(R0*R0)+(f2*l2)/(R2*R2);
    r0 = sqrt(R0);
    r1 = sqrt(R1);
    r2 = sqrt(R2);
    f0 = r0*r0*r0;
    f1 = r1*r1*r1;
    f2 = r2*r2*r2;
    l0 = (p0 - pb).norm();
    l1 = (p1 - pb).norm();
    l2 = (p2 - pb).norm();
   
    std::cout << "R0 : " << R0 << " and R1 " << R1 << " R2 " << R2 << "\n";
    std::cout << "pbNew[0]  : " << pb[0] << " and pbNew[1] " << pb[1] << "\n";
  }
  
  std::cout << "r0 : " << r0 << " and r1 " << r1 << " r2 " << r2 << "\n";
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
  
  board.saveEPS(file_eps.c_str());
  board.saveSVG(file_svg.c_str());
  
  return 0;
}

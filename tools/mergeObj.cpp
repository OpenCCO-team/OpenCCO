/**
 *  mergeObj program (used in OpenCCO implementation) 
 *  Copyright (C) 2023 B. Kerautret;  Phuc Ngo, N. Passat H. Talbot and C. Jaquet
 *
 *  This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **/

#include "CLI11.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


using namespace std;

int main (int argc, char * argv [] ) {
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
   
  std::string inputFileNameBase;
  std::string inputFileNameAdded;
  std::string outputFileName {"resultMerged.obj"};

  std::string nameGrObj1 {""};
  std::string nameGrObj2 {""};
  bool cleanNormals {false};
  std::vector<double> material1;
  std::vector<double> material2;

  app.description("Merge two mesh files given in obj format");
  app.add_option("-i,--input,1", inputFileNameBase, "Input base obj file" )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-a,--add,2", inputFileNameAdded, "Second input obj file" )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output,3", outputFileName, "Output obj file" )
    ->required();
  app.add_option("--materialOne", material1, "define the material (RGBA) of the first mesh" )
    ->expected(4);
  app.add_option("--materialTwo", material2, "define the material (RGBA) of the second mesh" )
    ->expected(4);

   
  app.add_option("--nameGrp1,4", nameGrObj1, "Name to reference the first object in the source file" );
  app.add_option("--nameGrp2,5", nameGrObj2, "Name to reference the seconf object in the source file" );
  app.add_flag("-c,--clean", cleanNormals, "Clean normal information by removing normals references");
   
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
   

    ofstream resultFile;
    resultFile.open (outputFileName, std::ofstream::out);
    if (material1.size() == 4){
        resultFile << "newmtl material1" << std::endl;
        resultFile << "Ka 0.9 0.9 0.9" << std::endl;
        resultFile << "Kd " << material1[0]
        << " " << material1[1]
        << " " << material1[2] << std::endl;
        resultFile << "Ks 0.2 0.2 0.2" << std::endl;
        resultFile << "Ns 20  # shininess" << std::endl;
        resultFile << "d " << material1[3] << " # transparency" << std::endl;
    }
    if (nameGrObj1 != ""){
        resultFile << "g " << nameGrObj1 << std::endl;
    }
    if (material1.size() == 4){
        resultFile << "usemtl material1" << std::endl;
    }
    
  string line;
  std::ifstream myfile;
  myfile.open (inputFileNameBase, std::ifstream::in);
  unsigned int nbVertex = 0;

  
  while ( getline (myfile,line) )
  {
    if (line[0] == 'v' && line[1] == ' '){
      nbVertex++;
       
      if( cleanNormals) {
        // limit to size 3
        std::istringstream ss(line);
        unsigned int n=0;
        std::string token;
        while(std::getline(ss, token, ' ') && n<4) {
          resultFile << token << " "; 
          n++;
        }
        resultFile << "\n";
      }
      else if (line[0] != 'v' || line[1] != 'n') {
        resultFile << line << '\n';
      }
    }
    else if (line [0] != 'o') {
      resultFile << line << '\n';         
    }
      
  }
  myfile.close();
 
    if (material2.size() == 4){
        resultFile << "newmtl material2" << std::endl;
        resultFile << "Ka 0.9 0.9 0.9" << std::endl;
        resultFile << "Kd " << material2[0]
        << " " << material2[1]
        << " " << material2[2] << std::endl;
        resultFile << "Ks 0.2 0.2 0.2" << std::endl;
        resultFile << "Ns 20  # shininess" << std::endl;
        resultFile << "d " << material2[3] << " # transparency" << std::endl;
    }

  // File 2:
  if (nameGrObj2 != ""){
    resultFile << "g " << nameGrObj2 << std::endl;
  }
    if (material2.size() == 4){
        resultFile << "usemtl material2" << std::endl;
    }
  std::ifstream myfile2;
  myfile2.open (inputFileNameAdded, std::ifstream::in);
  
  while ( getline (myfile2,line) )
  {
    if (line[0]=='f'){
      resultFile << "f ";
      std::istringstream ss(line);
      std::string token;
      while(std::getline(ss, token, ' ')) {
        std::istringstream fe (token);
        std::string face = "";
        std::getline(fe, face, '/');
        if (face!="f"){
          resultFile << nbVertex+(stoi(face)) << " ";
        }
      }
      resultFile << "\n";
    }
    else if (line [0] != 'o') {
      resultFile << line << '\n';
    }
    
    
  }
  myfile2.close();
  resultFile.close();
  return 0;
}

// g++ -std=c++11  mergeObj.cpp -o mergeObj
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

   
  app.description("Merge two mesh files given in obj format");
  app.add_option("-i,--input,1", inputFileNameBase, "Input base obj file" )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-a,--add,2", inputFileNameAdded, "Second input obj file" )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output,3", outputFileName, "Output obj file" )
    ->required();
 
   
  app.add_option("--nameGrp1,4", nameGrObj1, "Name to reference the first object in the source file" );
  app.add_option("--nameGrp2,5", nameGrObj2, "Name to reference the seconf object in the source file" );
  app.add_flag("-c,--clean", cleanNormals, "Clean normal information by removing normals references");
   
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
   

  ofstream resultFile;
  resultFile.open (outputFileName, std::ofstream::out);
  if (nameGrObj1 != ""){
    resultFile << "g " << nameGrObj1 << std::endl;
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
    else{
      resultFile << line << '\n';         
    }
  }
  myfile.close();
 
    

  // File 2:
  if (nameGrObj2 != ""){
    resultFile << "g " << nameGrObj2 << std::endl;
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
    else
    {
      resultFile << line << '\n';
    }
    
    
  }
  myfile2.close();
  resultFile.close();
  return 0;
}

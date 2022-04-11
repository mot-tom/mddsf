#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <complex>
#include <math.h>
#include <mkl_dfti.h>
using namespace std;

#define MAXELEMENT 16
typedef enum {Si, Ge, C, Sn, Mo, S, O, X1, X2} AtomIonType;

struct info        // information of atom
{
  string name;     // atom name
  double x, y, z;  // coordinate of atom
  double sumx2, sumy2, sumz2;
};

struct optpara{
	string read;
	string wave;
       int part;
       int cutoff;
	string dsflog;
       int k_part;
       int f_part;
	string polytype;
       double dt;
       double maxhz;
       string range;
};

int xyz_vector(int argc, char* argv[], optpara opts)
{  
  /***********************************/

  int num_atom = 300000;
  const int max_atom = 300000;         // maximum number of atom
  info* zero_atom = new info[max_atom];     // each atom at equilibrium position
  string title, temp_line;
  
  ifstream file_fin(opts.read.c_str());
  
  if(!file_fin)
  {
    cerr << "Can't open " << opts.read << "!" << endl;
    return 1;
  }
  
  cerr << "Read \"" << opts.read << "\"." << endl;

  int i=-1;
  int j=0;
  while(getline(file_fin, temp_line))
  {
    if(i==-1 || j>=num_atom)
    {
      j=0;
      i++;
      if(i==0) num_atom = atoi(temp_line.c_str());
      getline(file_fin, temp_line);
      title = temp_line;
    }
    else
    {
      istringstream line_iss(temp_line.c_str());
      if(i==0)
      {
        line_iss >> zero_atom[j].name >> zero_atom[j].x >> zero_atom[j].y >> zero_atom[j].z;
      }
      else
      {
        info atom;
        line_iss >> atom.name >> atom.x >> atom.y >> atom.z;
        zero_atom[j].sumx2 += pow(atom.x - zero_atom[j].x,2.0);
        zero_atom[j].sumy2 += pow(atom.y - zero_atom[j].y,2.0);
        zero_atom[j].sumz2 += pow(atom.z - zero_atom[j].z,2.0);
      }
      j++;
    }
  }

  ostringstream osfile;
  osfile << opts.read << ".dat";
  ofstream fout(osfile.str().c_str());
  if(!fout)
    {
      cout << "Can't open " << osfile.str() << "!" << endl;
      return 1;
    }

  fout << num_atom << endl;
  fout << title << endl;
  for(int k=0; k<num_atom; ++k)
  {
    zero_atom[j].sumx2 /= num_atom;
    zero_atom[j].sumy2 /= num_atom;
    zero_atom[j].sumz2 /= num_atom;
    fout << zero_atom[k].name << " " << zero_atom[k].x << " " << zero_atom[k].y << " " << zero_atom[k].z << " " << zero_atom[k].sumx2 << " " << zero_atom[k].sumy2 << " " << zero_atom[k].sumz2 << endl;
  }
  
  return 0;
}

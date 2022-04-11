#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <math.h>
using namespace std;

#define MAXELEMENT 16
typedef enum {Si, Ge, C, Sn, Mo, S, O} AtomIonType;

struct info        // information of atom
{
  double mass;     // atom mass
  double x, y, z;  // coordinate of atom
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

int hessian_lmp(int argc, char* argv[], optpara opts)
{  
  /***********************************/

  const int max_atom = 300000;         // maximum number of atom
  const int max_file = 20001;           // maximum number of file 
  
  const double pi = 3.14159265;        // pi

  string output = "hessian";           // base name of the output file
  string ext = ".dat";                 // name of extension for the output file
  
  /***********************************/

  double llo, lhi, ls = 0;
  double leng[3] = {0.0, 0.0, 0.0};    // length of each side of crystal  
  double tilt[3][3] = {{0.0, 0.0, 0.0}, 
                       {0.0, 0.0, 0.0},  // x
                       {0.0, 0.0, 0.0}}; // tilt length of each side of crystal  

  int blank = 0;                       // number of present blank in the file
  int all_atom = 0;                    // number of all atom
  int num_data = 0;                    // number of file 

  int time_step = 0;
  int mass;
  string first;                        // first column of line temporarily
  string temp_line, temp_zero, temp_calc;   // line of string (list/first file/calculation file)
  
  // allocate dynamic memory
  
  info* zero_atom = new info[max_atom];     // each atom at equilibrium position
  info* r_atom = new info[max_atom];
  info* f_atom = new info[max_atom];
  
  string* file = new string[max_file];      // each file in the input list
  
  ifstream file_fin(opts.read.c_str());
  
  if(!file_fin)
  {
    cout << "Can't open " << opts.read << "!" << endl;
    return 1;
  }
  
  cout << "Read \"" << opts.read << "\"." << endl;

  while(getline(file_fin, temp_line))
  {
    istringstream is(temp_line.c_str());
    
    if(temp_line != "" && num_data+1 > max_file)
    {
      cout << "Dynamic memory space is not enough for data!" << endl;
      return 1;
    }
    
    is >> file[num_data];
    //cout << num_data << " " << file[num_data] << endl;
    num_data++;
  }
  file_fin.close();

  ifstream fin(file[0].c_str());

  if(!fin)
  {
    cout << "Can't open " << opts.read << "!" << endl;
    return 1;
  }

  while(getline(fin, temp_line))
  {
    getline(fin, temp_zero);
    istringstream zero_time(temp_zero.c_str());
    zero_time >> time_step;

    if(num_data+1 > max_file)
    {
      cout << "Dynamic memory space is not enough for data!" << endl;
      return 1;
    }
    if(num_data == 0 && time_step != 0)
    {
      cout << "File of first structure does not exist in the list!" << endl;
      return 1;
    }

    getline(fin, temp_zero);
    getline(fin, temp_zero);
    istringstream zero_atom(temp_zero.c_str());
    zero_atom >> all_atom;

    if(all_atom+1 > max_atom)
    {
      cout << "Dynamic memory space is not enough for atom!" << endl;
      return 1;
    }

    getline(fin, temp_zero);
    getline(fin, temp_zero);
      istringstream zero_lat_x(temp_zero.c_str());
      zero_lat_x >> llo >> lhi >> ls;
      leng[0] = lhi - llo;
      tilt[0][1] = ls * 0.5;
      tilt[1][0] = ls * 0.5;
    getline(fin, temp_zero);
      istringstream zero_lat_y(temp_zero.c_str());
      zero_lat_y >> llo >> lhi >> ls;
      leng[1] = lhi - llo;
      tilt[0][2] = ls * 0.5;
      tilt[2][0] = ls * 0.5;
    getline(fin, temp_zero);
      istringstream zero_lat_z(temp_zero.c_str());
      zero_lat_z >> llo >> lhi >> ls;
      leng[2] = lhi - llo;
      tilt[1][2] = ls * 0.5;
      tilt[2][1] = ls * 0.5;
      cout << "|ul00 ul01 ul02| ";
      cout << showpoint << "|" << setw(6) << leng[0]    << " " << tilt[0][1] << " " << tilt[0][2] << "|" << endl;
      cout << "|ul10 ul11 ul12|=";
      cout << showpoint << "|" << setw(6) << tilt[1][0] << " " << leng[1]    << " " << tilt[1][2] << "|" << endl;
      cout << "|ul20 ul21 ul22| ";
      cout << showpoint << "|" << setw(6) << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << "|" << endl;

    for(int i=0; i<all_atom+1; i++)
      getline(fin, temp_zero);
  }

  fin.close();

       string atom_name;
	int x_dir=4, y_dir=4, z_dir=4;
       double lattice[3];
	if(opts.dsflog.length() > 0){
             istringstream size_iss(opts.dsflog.c_str());
             size_iss >> atom_name >> x_dir >> y_dir >> z_dir;
             if(x_dir == 0 || y_dir == 0 || z_dir == 0){
		  cerr << "Structure size is wrong" << endl;
               return 1;
               }
	} else {
		  cerr << "Input option is wrong" << endl;
               return 1;
	}
       if(atom_name == "Si"){
         lattice[0] = 5.431;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else if(atom_name == "Ge"){
         lattice[0] = 5.658;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else if(atom_name == "Sn"){
         lattice[0] = 6.489;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else if(atom_name == "C"){
         lattice[0] = 3.567;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else if(atom_name == "3C-SiC"){
         lattice[0] = 4.311;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else if(atom_name == "3C-SiGe"){
         lattice[0] = 5.431 + 0.2*0.5 + 0.027*0.5*0.5;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else if(atom_name == "SiGe"){
         lattice[0] = (leng[0]/x_dir + leng[1]/y_dir + leng[2]/z_dir)/3.0;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else if(atom_name == "GaAs"){
         lattice[0] = 5.65325;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else if(atom_name == "InAs"){
         lattice[0] = 6.0583;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else if(atom_name == "3C-BN"){
         lattice[0] = 3.615;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
       }else{
         cerr << "Atom name is wrong" << endl;
         return 1;
        }
  
  // read the information of the first structure

  fin.open(file[0].c_str());
  for(int i=0; i<9; i++)
    getline(fin, temp_zero);
  for(int i=0; i<all_atom; i++)
  {
    int label = 0;
    getline(fin, temp_zero);
    istringstream zero_iss(temp_zero.c_str());
    zero_iss >> label;
    zero_iss >> zero_atom[label-1].mass;
    zero_iss >> zero_atom[label-1].x >> zero_atom[label-1].y >> zero_atom[label-1].z;
  }

  fin.close();

  cout << "---------------------------------------------" << endl;
  cout << "all atom:" << all_atom << endl;
  cout << "File:" << num_data << endl;
  cout << "---------------------------------------------" << endl;

  cout << "######################################################" << endl;
  cout << " HESSIAN MATRIX [eV x angstrom^-2 ]" << endl;
  cout << " i j  mI[a.m.u.]  mJ[a.m.u.]  rI[x,y,z]  rJ[x,y,z] Hij[eV x angstrom^-2]" << endl;
  cout << "######################################################" << endl;

  
  // read the information of all atoms in each file
  double da = 0.0;  

  for(int l=1; l<=num_data-1; ++l)
  { 
    info temp_r, temp_f;

    fin.open(file[l].c_str());

    //move to x direction
    for(int i=0; i<9; i++)
      getline(fin, temp_calc);
    for(int i=0; i<all_atom; i++)
    {
      int label = 0;
      getline(fin, temp_calc);
      istringstream iss(temp_calc.c_str());
      iss >> label >> temp_r.mass;
      iss >> temp_r.x >> temp_r.y >> temp_r.z;
      iss >> f_atom[label-1].x >> f_atom[label-1].y >> f_atom[label-1].z;
      if(label == 1 && da == 0) da = temp_r.x;
   }
    for(int i=0; i<9; i++)
      getline(fin, temp_calc);
    for(int i=0; i<all_atom; i++)
    {  
      int label = 0;
      getline(fin, temp_calc);
      istringstream iss(temp_calc.c_str());
      iss >> label >> temp_r.mass;
      iss >> temp_r.x >> temp_r.y >> temp_r.z;
      iss >> temp_f.x >> temp_f.y >> temp_f.z;
      temp_r.x = zero_atom[label-1].x - zero_atom[l-1].x;
      temp_r.y = zero_atom[label-1].y - zero_atom[l-1].y;
      temp_r.z = zero_atom[label-1].z - zero_atom[l-1].z;
      if(temp_r.x > leng[0]/2.0) temp_r.x -= leng[0];
      else if(temp_r.x < -leng[0]/2.0-0.0001) temp_r.x += leng[0];
      if(temp_r.y > leng[1]/2.0) temp_r.y -= leng[1];
      else if(temp_r.y < -leng[1]/2.0-0.0001) temp_r.y += leng[1];
      if(temp_r.z > leng[2]/2.0) temp_r.z -= leng[2];
      else if(temp_r.z < -leng[2]/2.0-0.0001) temp_r.z += leng[2];
      temp_r.x += zero_atom[l-1].x;
      temp_r.y += zero_atom[l-1].y;
      temp_r.z += zero_atom[l-1].z;
      //output x direction
      cout << setprecision(10) << fixed << (l-1)*3 << " " << (label-1)*3 << " " << zero_atom[l-1].mass << " " << zero_atom[label-1].mass << " ";
      cout << setw(12) << zero_atom[l-1].x << " " << zero_atom[l-1].y << " " << zero_atom[l-1].z << " ";
      cout << setw(12) << temp_r.x << " " << temp_r.y << " " << temp_r.z << " ";
      cout << setprecision(27) << setw(30) << fixed << (temp_f.x - f_atom[label-1].x)/(2*da) << endl;
      //output y direction
      cout << setprecision(10) << (l-1)*3 << " " << (label-1)*3+1 << " " << zero_atom[l-1].mass << " " << zero_atom[label-1].mass << " ";
      cout << setw(12) << zero_atom[l-1].x << " " << zero_atom[l-1].y << " " << zero_atom[l-1].z << " ";
      cout << setw(12) << temp_r.x << " " << temp_r.y << " " << temp_r.z << " ";
      cout << setprecision(27) << setw(30) << fixed << (temp_f.y - f_atom[label-1].y)/(2*da) << endl;
      //output z direction
      cout << setprecision(10) << (l-1)*3 << " " << (label-1)*3+2 << " " << zero_atom[l-1].mass << " " << zero_atom[label-1].mass << " ";
      cout << setw(12) << zero_atom[l-1].x << " " << zero_atom[l-1].y << " " << zero_atom[l-1].z << " ";
      cout << setw(12) << temp_r.x << " " << temp_r.y << " " << temp_r.z << " ";
      cout << setprecision(27) << setw(30) << fixed << (temp_f.z - f_atom[label-1].z)/(2*da) << endl;
   }

    //move to y direction
    for(int i=0; i<9; i++)
      getline(fin, temp_calc);
    for(int i=0; i<all_atom; i++)
    {
      int label = 0;
      getline(fin, temp_calc);
      istringstream iss(temp_calc.c_str());
      iss >> label >> temp_r.mass;
      iss >> temp_r.x >> temp_r.y >> temp_r.z;
      iss >> f_atom[label-1].x >> f_atom[label-1].y >> f_atom[label-1].z;
   }
    for(int i=0; i<9; i++)
      getline(fin, temp_calc);
    for(int i=0; i<all_atom; i++)
    {  
      int label = 0;
      getline(fin, temp_calc);
      istringstream iss(temp_calc.c_str());
      iss >> label >> temp_r.mass;
      iss >> temp_r.x >> temp_r.y >> temp_r.z;
      iss >> temp_f.x >> temp_f.y >> temp_f.z;
      temp_r.x = zero_atom[label-1].x - zero_atom[l-1].x;
      temp_r.y = zero_atom[label-1].y - zero_atom[l-1].y;
      temp_r.z = zero_atom[label-1].z - zero_atom[l-1].z;
      if(temp_r.x > leng[0]/2.0) temp_r.x -= leng[0];
      else if(temp_r.x < -leng[0]/2.0-0.0001) temp_r.x += leng[0];
      if(temp_r.y > leng[1]/2.0) temp_r.y -= leng[1];
      else if(temp_r.y < -leng[1]/2.0-0.0001) temp_r.y += leng[1];
      if(temp_r.z > leng[2]/2.0) temp_r.z -= leng[2];
      else if(temp_r.z < -leng[2]/2.0-0.0001) temp_r.z += leng[2];
      temp_r.x += zero_atom[l-1].x;
      temp_r.y += zero_atom[l-1].y;
      temp_r.z += zero_atom[l-1].z;
      //output x direction
      cout << setprecision(10) << (l-1)*3+1 << " " << (label-1)*3 << " " << zero_atom[l-1].mass << " " << zero_atom[label-1].mass << " ";
      cout << setw(12) << zero_atom[l-1].x << " " << zero_atom[l-1].y << " " << zero_atom[l-1].z << " ";
      cout << setw(12) << temp_r.x << " " << temp_r.y << " " << temp_r.z << " ";
      cout << setprecision(27) << setw(30) << fixed << (temp_f.x - f_atom[label-1].x)/(2*da) << endl;
      //output y direction
      cout << setprecision(10) << (l-1)*3+1 << " " << (label-1)*3+1 << " " << zero_atom[l-1].mass << " " << zero_atom[label-1].mass << " ";
      cout << setw(12) << zero_atom[l-1].x << " " << zero_atom[l-1].y << " " << zero_atom[l-1].z << " ";
      cout << setw(12) << temp_r.x << " " << temp_r.y << " " << temp_r.z << " ";
      cout << setprecision(27) << setw(30) << fixed << (temp_f.y - f_atom[label-1].y)/(2*da) << endl;
      //output z direction
      cout << setprecision(10) << (l-1)*3+1 << " " << (label-1)*3+2 << " " << zero_atom[l-1].mass << " " << zero_atom[label-1].mass << " ";
      cout << setw(12) << zero_atom[l-1].x << " " << zero_atom[l-1].y << " " << zero_atom[l-1].z << " ";
      cout << setw(12) << temp_r.x << " " << temp_r.y << " " << temp_r.z << " ";
      cout << setprecision(27) << setw(30) << fixed << (temp_f.z - f_atom[label-1].z)/(2*da) << endl;
   }

    //move to z direction
    for(int i=0; i<9; i++)
      getline(fin, temp_calc);
    for(int i=0; i<all_atom; i++)
    {
      int label = 0;
      getline(fin, temp_calc);
      istringstream iss(temp_calc.c_str());
      iss >> label >> temp_r.mass;
      iss >> temp_r.x >> temp_r.y >> temp_r.z;
      iss >> f_atom[label-1].x >> f_atom[label-1].y >> f_atom[label-1].z;
   }
    for(int i=0; i<9; i++)
      getline(fin, temp_calc);
    for(int i=0; i<all_atom; i++)
    {  
      int label = 0;
      getline(fin, temp_calc);
      istringstream iss(temp_calc.c_str());
      iss >> label >> temp_r.mass;
      iss >> temp_r.x >> temp_r.y >> temp_r.z;
      iss >> temp_f.x >> temp_f.y >> temp_f.z;
      temp_r.x = zero_atom[label-1].x - zero_atom[l-1].x;
      temp_r.y = zero_atom[label-1].y - zero_atom[l-1].y;
      temp_r.z = zero_atom[label-1].z - zero_atom[l-1].z;
      if(temp_r.x > leng[0]/2.0) temp_r.x -= leng[0];
      else if(temp_r.x < -leng[0]/2.0-0.0001) temp_r.x += leng[0];
      if(temp_r.y > leng[1]/2.0) temp_r.y -= leng[1];
      else if(temp_r.y < -leng[1]/2.0-0.0001) temp_r.y += leng[1];
      if(temp_r.z > leng[2]/2.0) temp_r.z -= leng[2];
      else if(temp_r.z < -leng[2]/2.0-0.0001) temp_r.z += leng[2];
      temp_r.x += zero_atom[l-1].x;
      temp_r.y += zero_atom[l-1].y;
      temp_r.z += zero_atom[l-1].z;
      //output x direction
      cout << setprecision(10) << (l-1)*3+2 << " " << (label-1)*3 << " " << zero_atom[l-1].mass << " " << zero_atom[label-1].mass << " ";
      cout << setw(12) << zero_atom[l-1].x << " " << zero_atom[l-1].y << " " << zero_atom[l-1].z << " ";
      cout << setw(12) << temp_r.x << " " << temp_r.y << " " << temp_r.z << " ";
      cout << setprecision(27) << setw(30) << fixed << (temp_f.x - f_atom[label-1].x)/(2*da) << endl;
      //output y direction
      cout << setprecision(10) << (l-1)*3+2 << " " << (label-1)*3+1 << " " << zero_atom[l-1].mass << " " << zero_atom[label-1].mass << " ";
      cout << setw(12) << zero_atom[l-1].x << " " << zero_atom[l-1].y << " " << zero_atom[l-1].z << " ";
      cout << setw(12) << temp_r.x << " " << temp_r.y << " " << temp_r.z << " ";
      cout << setprecision(27) << setw(30) << fixed << (temp_f.y - f_atom[label-1].y)/(2*da) << endl;
      //output z direction
      cout << setprecision(10) << (l-1)*3+2 << " " << (label-1)*3+2 << " " << zero_atom[l-1].mass << " " << zero_atom[label-1].mass << " ";
      cout << setw(12) << zero_atom[l-1].x << " " << zero_atom[l-1].y << " " << zero_atom[l-1].z << " ";
      cout << setw(12) << temp_r.x << " " << temp_r.y << " " << temp_r.z << " ";
      cout << setprecision(27) << setw(30) << fixed << (temp_f.z - f_atom[label-1].z)/(2*da) << endl;
   }

    fin.close();
  }

  cout << "######################################################" << endl;
  
  // free dynamic memory
  
  delete[] r_atom;
  delete[] f_atom;
  delete[] zero_atom;
  delete[] file;
  
  return 0;
}

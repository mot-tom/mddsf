#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <complex>
#include <math.h>
#include <mkl_dfti.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

// Outline of this program
/*
  Calculate dynamic structure factor at each wave vector
  for phonon dispersion relation.
*/

#define MAXELEMENT 16
typedef enum {Si, Ge, C, Sn, Mo, S, O} AtomIonType;

struct info        // information of atom
{
  string name;     // atom name
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

int dsf_omp(int argc, char* argv[], optpara opts, int num_omp)
{  
  /***********************************/

  #ifdef _OPENMP
  if(num_omp > 1)
    omp_set_num_threads(num_omp);
  else
    omp_set_num_threads(1);
  printf("set %d OpenMP threads\n",omp_get_max_threads());
  #endif
  
  int part = 10;                 // value to divide wave vector  
  const int max_atom = 300000;         // maximum number of atom
  const int max_file = 20001;           // maximum number of file 
  
  double dt = 1e-14;             // time step[s]  
  const double pi = 3.14159265;        // pi
  const double tera = 1e-12;           // convert unit into tera order

  int cutoff = 1;          // value to cutoff wave vector
  
  string wave = "100";                 // direction index of wave vector
  string output = "dsf";               // base name of the output file
  string ext = ".dat";                 // name of extension for the output file

  if(opts.wave.c_str()) wave = opts.wave;
  if(opts.part) part = opts.part;
  if(opts.cutoff) cutoff = opts.cutoff;
  if(opts.dt) dt = opts.dt;

  double range[3][2] = {{-1.0, 2.0}, 
                       {-1.0, 2.0},     // x
                       {-1.0, 2.0}};    // read range of crystal
  int range_on = 0;
  if(opts.range.length() > 0) {
    istringstream range_iss(opts.range.c_str());
    range_iss >> range[0][0] >> range[0][1] >> range[1][0] >> range[1][1] >> range[2][0] >> range[2][1];
    if(range[0][0] == range[0][1] || range[1][0] == range[1][1] || range[2][0] == range[2][1]) {
      cerr << "Structure range is wrong" << endl;
      return 1;
     }
    range_on = 1;
   }

  cout << "wave " << wave << endl;
  cout << "ldiv " << part << endl;
  cout << "coff " << cutoff << endl;
  
  /***********************************/
  
  double leng[3] = {0.0, 0.0, 0.0};    // length of each side of crystal  
  double tilt[3][3] = {{0.0, 0.0, 0.0}, 
                       {0.0, 0.0, 0.0},  // x
                       {0.0, 0.0, 0.0}}; // tilt length of each side of crystal  
  double lattice[3];                   // lattice constant[A]
  double dir[3];                       // each direction index of wave vector
  double coe[3];                       // coefficient of each direction in the wave vector
  
  double max_dir = 1;                  // maximum displaying wave vector  
  int max_hertz;                       // maximum displaying hertz  
  int blank = 0;                       // number of present blank in the file
  int all_atom = 0;                    // number of all atom
  int num_data = 0;                    // number of file 

  int atomnum[MAXELEMENT];             // number of atom, respectively
  double atomconc[MAXELEMENT];         // concentration of atom, respectively

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=0;i<MAXELEMENT;i++){
    atomnum[i] = 0;
    atomconc[i] = 0;
  }
  
  string atom_name;                    // atom name(Si or Ge)
  string first;                        // first column of line temporarily
  string temp_file, temp_zero;   // line of string (list/first file/calculation file)
  
  // allocate dynamic memory
  
  info* zero_atom = new info[max_atom];     // each atom at equilibrium position
  
  string* file = new string[max_file];      // each file in the input list
  
  // read direction index of wave vector
  
  for(int l=0; l<=2; ++l)
  {
    string pick;        // store part of string wave
    pick = wave.substr(l,1);
    
    istringstream is_pick(pick.c_str());
    is_pick >> dir[l];
    
    if(dir[l] != 0.0 && dir[l] != 1.0 && dir[l] != 2.0)
    {
      cout << "Setting of wave vector direction is wrong!" << endl;
      return 1;
    }
  }

  double norm = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
  dir[0] /= sqrt(norm); dir[1] /= sqrt(norm); dir[2] /= sqrt(norm);

  double dir_la[3];                    // each direction index of LA phonon
  double dir_lo[3];                    // each direction index of LO phonon
  double dir_t1[3];                    // each direction index of T1 phonon
  double dir_t2[3];                    // each direction index of T2 phonon

  if(wave == "100") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 2.0; dir_lo[1] = 0.0; dir_lo[2] = 0.0;
    dir_t1[0] = 0.0; dir_t1[1] = 2.0; dir_t1[2] = 0.0;
    dir_t2[0] = 0.0; dir_t2[1] = 0.0; dir_t2[2] = 2.0;
  }
  else if(wave == "010") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 0.0; dir_lo[1] = 2.0; dir_lo[2] = 0.0;
    dir_t1[0] = 2.0; dir_t1[1] = 0.0; dir_t1[2] = 0.0;
    dir_t2[0] = 0.0; dir_t2[1] = 0.0; dir_t2[2] = 2.0;
  }
  else if(wave == "001") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 0.0; dir_lo[1] = 0.0; dir_lo[2] = 2.0;
    dir_t1[0] = 0.0; dir_t1[1] = 2.0; dir_t1[2] = 0.0;
    dir_t2[0] = 2.0; dir_t2[1] = 0.0; dir_t2[2] = 0.0;
  }
  else if(wave == "110") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 1.4; dir_lo[1] =-1.4; dir_lo[2] = 0.0;
    dir_t1[0] = 1.4; dir_t1[1] = 1.4; dir_t1[2] = 0.0;
    dir_t2[0] = 0.0; dir_t2[1] = 0.0; dir_t2[2] = 2.0;
  }
  else if(wave == "220") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 1.4; dir_lo[1] = 1.4; dir_lo[2] = 0.0;
    dir_t1[0] = 1.4; dir_t1[1] = 1.4; dir_t1[2] = 0.0;
    dir_t2[0] = 0.0; dir_t2[1] = 0.0; dir_t2[2] = 2.0;
  }
  else if(wave == "111") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 1.4; dir_lo[1] =-0.8; dir_lo[2] = 1.15;
    dir_t1[0] = 1.4; dir_t1[1] = 0.8; dir_t1[2] =-1.15;
    dir_t2[0] = 0.0; dir_t2[1] = 1.6; dir_t2[2] = 1.15;
  }
  else{
      cout << "Setting of wave vector direction is wrong!" << endl;
      return 1;
  }
  
  // read the MD file from the input list 

  ifstream file_fin(opts.read.c_str());
  
  if(!file_fin)
  {
    cout << "Can't open " << opts.read << "!" << endl;
    return 1;
  }
  
  cout << "Read \"" << opts.read << "\"." << endl;

  while(getline(file_fin, temp_file))
  {
    istringstream is(temp_file.c_str());
    
    if(temp_file != "" && num_data+1 > max_file)
    {
      cout << "Dynamic memory space is not enough for data!" << endl;
      return 1;
    }
    
    is >> file[num_data];
    //cout << num_data << " " << file[num_data] << endl;
    num_data++;
  }

  file_fin.close();
  
  // check wheter the file of first structure is in the input list
  
  if(file[0].find("00000", 0) == string::npos)
  {
    cout << "File of first structure does not exist in the list!" << endl;
    return 1;
  }
  
  ifstream zero_fin(file[0].c_str());
  
  if(!zero_fin)
  {
    cout << "Can't open " << file[0] << "!" << endl;
    return 1;
  }
  
  // read the information of the first structure

  int lc = 0;
  while(getline(zero_fin, temp_zero))
  {
    istringstream zero_iss(temp_zero.c_str());

    if(blank == 2)
    {
      zero_iss >> zero_atom[0].x >> zero_atom[0].y >> zero_atom[0].z;

      if(zero_atom[0].x != 0.0 && lc == 0)
        leng[0] = zero_atom[0].x;
      if(zero_atom[0].y != 0.0 && lc == 0)
        tilt[0][1] = zero_atom[0].y;
      if(zero_atom[0].z != 0.0 && lc == 0)
        tilt[0][2] = zero_atom[0].z;

      if(zero_atom[0].x != 0.0 && lc == 1)
        tilt[1][0] = zero_atom[0].x;
      if(zero_atom[0].y != 0.0 && lc == 1)
        leng[1] = zero_atom[0].y;
      if(zero_atom[0].z != 0.0 && lc == 1)
        tilt[1][2] = zero_atom[0].z;

      if(zero_atom[0].x != 0.0 && lc == 2)
        tilt[2][0] = zero_atom[0].x;
      if(zero_atom[0].y != 0.0 && lc == 2)
        tilt[2][1] = zero_atom[0].y;
      if(zero_atom[0].z != 0.0 && lc == 2)      
        leng[2] = zero_atom[0].z;

      if(lc < 3){
        cout << "|ul" << lc << "0 ul" << lc << "1 ul" << lc << "2|";
        if(lc == 1) cout << "=";
        else cout << " ";
        if(lc == 0)
          cout << showpoint << "|" << setw(6) << leng[0]    << " " << tilt[0][1] << " " << tilt[0][2] << "|";
        else if(lc == 1)
          cout << showpoint << "|" << setw(6) << tilt[1][0] << " " << leng[1]    << " " << tilt[1][2] << "|";
        else if(lc == 2)
          cout << showpoint << "|" << setw(6) << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << "|";
        cout << endl;
       }
      lc++;
    }

    if(blank == 3)
    {
      if(all_atom+1 > max_atom)
      {
        cout << "Dynamic memory space is not enough for atom!" << endl;
        return 1;
      }
      
      zero_iss >> first;
      
      // determin lattice constant and maximum hertz for displaying dipersion relation
      
      if(all_atom == 0)
      {
        if(first != "Si" && first != "Ge" && first != "Sn" && first != "C"
        && first != "Mo" && first != "S" && first != "S1" && first != "S2" && first != "O")
        {
          cout << "Setting of atom name is wrong!" << endl;
          return 1;
        }
      }
      
      if(first != "V")
      {
        zero_atom[all_atom].name = first;
        zero_iss >> zero_atom[all_atom].x >> zero_atom[all_atom].y >> zero_atom[all_atom].z;
        all_atom++;

        if(first == "Si")
          atomnum[Si]++;
        if(first == "Ge")
          atomnum[Ge]++;
        if(first == "Sn")
          atomnum[Sn]++;
        if(first == "C")
          atomnum[C]++;
        if(first == "Mo")
          atomnum[Mo]++;
        if(first == "S")
          atomnum[S]++;
        if(first == "S1")
          atomnum[S]++;
        if(first == "S2")
          atomnum[S]++;
        if(first == "O")
          atomnum[O]++;
      }
    }
    
    if(temp_zero == "")
      blank++;
    
  }

  if(atomnum[Mo] && atomnum[S])
  {
    if(atomnum[Si] || atomnum[Ge] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "MoS2";
    max_hertz = 20;
    lattice[0] = 3.127;
    lattice[1] = lattice[0]*cos(-30 * M_PI/180)*2.0;
    lattice[2] = 12.066;
  }
  else if(atomnum[Si] && atomnum[Ge] && atomnum[Sn])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "GeSiSn";
    max_hertz = 20;
    atomconc[Si] = (double)atomnum[Si] / (atomnum[Sn] + atomnum[Ge] + atomnum[Si]);
    atomconc[Sn] = (double)atomnum[Sn] / (atomnum[Sn] + atomnum[Ge] + atomnum[Si]);
    lattice[0] = 5.658 - 0.227*atomconc[Si] + 0.831*atomconc[Sn];
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
  }
  else if(atomnum[Si] && atomnum[Ge])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "SiGe";
    max_hertz = 20;
    atomconc[Ge] = (double)atomnum[Ge] / (atomnum[Ge] + atomnum[Si]);
    lattice[0] = 5.431 + 0.2*atomconc[Ge] + 0.027*atomconc[Ge]*atomconc[Ge];
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
  }
  else if(atomnum[Si] && atomnum[C])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[Ge] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "SiC";
    max_hertz = 70;
    atomconc[C] = (double)atomnum[C] / (atomnum[C] + atomnum[Si]);
    lattice[0] = 5.431 - 2.4239*atomconc[C] + 0.5705*atomconc[C]*atomconc[C];
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
  }
  else if(atomnum[Ge] && atomnum[Sn])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Si]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "GeSn";
    max_hertz = 12;
    atomconc[Sn] = (double)atomnum[Sn] / (atomnum[Sn] + atomnum[Ge]);
    lattice[0] = 5.658 + 0.7322*atomconc[Sn] + 0.0988*atomconc[Sn]*atomconc[Sn];
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
  }
  else if(atomnum[Si] && opts.polytype == "3H")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3H-Si";
    max_hertz = 20;
    lattice[0] = 3.84030;
    lattice[1] = lattice[0] / cos(-30 * M_PI/180) / 2.0;
    lattice[2] = 9.40676 / 3.0;
  }
  else if(atomnum[Si] && opts.polytype == "3T")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3T-Si";
    max_hertz = 20;
    lattice[0] = 3.84030;
    lattice[1] = lattice[0];
    lattice[2] = 5.431;
  }
  else if(atomnum[Si])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Si";
    max_hertz = 20;
    lattice[0] = 5.431;
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
  }
  else if(atomnum[Ge])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Ge";
    max_hertz = 12;
    lattice[0] = 5.658; 
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
  }
  else if(atomnum[Sn])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Si]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Sn";
    max_hertz = 8;
    lattice[0] = 6.489; 
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
  }
  else if(atomnum[C])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[Ge] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "C";
    max_hertz = 70;
    lattice[0] = 3.567;
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
  }
  if(max_hertz > tera*0.5/dt)
  {
    cout << "Setting of range for the frequency is wrong!" << endl;
    return 1;
  }

  zero_fin.close();

  if(opts.maxhz) max_hertz = opts.maxhz;
  
  coe[0] = sqrt(norm) * 2.0 * pi / lattice[0];
  coe[1] = sqrt(norm) * 2.0 * pi / lattice[1];
  coe[2] = sqrt(norm) * 2.0 * pi / lattice[2];
  
  cout << "---------------------------------------------" << endl;
  cout << "Atom:" << atom_name << "  Max Frequency:" << max_hertz << "[THz]" << endl;
  cout << "Lattice:" << lattice[0] << " x " << lattice[1] << " x " << lattice[2] << "[A]" << endl;
  cout << "Total atom:" << all_atom;
  if(atomnum[Si]) cout << "  Si atom:" << atomnum[Si];
  if(atomnum[Ge]) cout << "  Ge atom:" << atomnum[Ge];
  if(atomnum[Sn]) cout << "  Sn atom:" << atomnum[Sn];
  if(atomnum[C])  cout << "  C atom:"  << atomnum[C];
  if(atomnum[Mo]) cout << "  Mo atom:" << atomnum[Mo];
  if(atomnum[S]) cout  << "  S atom:"  << atomnum[S];
  if(atomnum[O]) cout  << "  O atom:"  << atomnum[O];
  cout << endl;
  if(atomconc[Ge]) cout << "Ge conc:" << atomconc[Ge] << endl;
  if(atomconc[Sn]) cout << "Sn conc:" << atomconc[Sn] << endl;
  cout << "File:" << num_data << "  ";
  cout << "MD time:" << (double)(num_data - 1) * dt / tera << "[ps]  ";
  cout << "Wave vector:[" << wave << "] from 0 to " << (double)1.0/cutoff << endl; 
  if(range_on) cout << "Input range: " << range[0][0] << " " << range[0][1] << " xlo xhi "
                                       << range[1][0] << " " << range[1][1] << " ylo yhi "
                                       << range[2][0] << " " << range[2][1] << " zlo zhi" << endl;
  cout << "---------------------------------------------" << endl;
  
  // read the information of all atoms in each file

  double sum_la[num_data];
  double sum_lo[num_data];
  double sum_t1[num_data];
  double sum_t2[num_data];
  double sum_ram[num_data];

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int m=0; m<num_data; ++m)
  {
    sum_la[m] = 0;
    sum_lo[m] = 0;
    sum_t1[m] = 0;
    sum_t2[m] = 0;
    sum_ram[m] = 0;
  }

  int break_omp = 0;

  cout << "Calc time: 0 - " << part << " / " << omp_get_max_threads() << " threads" << endl;  
  if(wave == "220")    
  cout << "Wave vector:("
       << dir[0]*coe[0]/cutoff/2.0 << ", " << dir[1]*coe[1]/cutoff/2.0 << ", " << dir[2]*coe[2]/cutoff/2.0 << ") - (0.00000, 0.00000, 0.00000)" << endl << endl; 
  else if(wave == "111") 
  cout << "Wave vector:(0.00000, 0.00000, 0.00000) - ("
       << dir[0]*coe[0]/cutoff/2.0 << ", " << dir[1]*coe[1]/cutoff/2.0 << ", " << dir[2]*coe[2]/cutoff/2.0 << ")" << endl << endl; 
  else    
  cout << "Wave vector:(0.00000, 0.00000, 0.00000) - ("
       << dir[0]*coe[0]/cutoff << ", " << dir[1]*coe[1]/cutoff << ", " << dir[2]*coe[2]/cutoff << ")" << endl << endl;
  cout << " Calculating ... fourier transform ..." << endl;
  cout << " 0%                  100%" << endl << " ";
  
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int l=0; l<=part; ++l)
  {    
    double n = l;
    if(wave == "220") n = part - n;
    else if(wave == "111") n = n/2.0;    
    double q[3] = {dir[0]*coe[0]*n/part/cutoff,
                   dir[1]*coe[1]*n/part/cutoff,
                   dir[2]*coe[2]*n/part/cutoff};        // wave vector[1/A]
    
    complex <double> i(0.0, 1.0);                // imaginary unit
    
    complex <double> c_la[num_data];             // scattering amplitude of longitudinal acoustic mode
    complex <double> c_lo[num_data];             // scattering amplitude of longitudinal optical mode
    complex <double> c_t1[num_data];             // scattering amplitude of transverse 1 mode
    complex <double> c_t2[num_data];             // scattering amplitude of transverse 2 mode
    complex <double> s_fourier_la[num_data];     // dynamic scattering function of longitudinal acoustic mode
    complex <double> s_fourier_lo[num_data];     // dynamic scattering function of longitudinal optical mode
    complex <double> s_fourier_t1[num_data];     // dynamic scattering function of transverse 1 mode
    complex <double> s_fourier_t2[num_data];     // dynamic scattering function of transverse 2 mode

    if(l%((int)part/10) == 0 ) cout << "XX" << flush;
    
    for(int m=1; m<num_data; ++m)
    {
      int stop_omp;
      #pragma omp atomic read
      stop_omp = break_omp;
      if (stop_omp) continue;

      double atom[3] = {0.0, 0.0, 0.0};     // atom information which is read from each file
      double displace[3] = {0.0, 0.0, 0.0}; // displacement from equilibirium position
      
      int label = 0;                        // number of present label in the input file
      string first_omp;
      
      int blank_omp = 0;
      
      c_la[m] = complex <double> (0.0, 0.0);
      c_lo[m] = complex <double> (0.0, 0.0);
      c_t1[m] = complex <double> (0.0, 0.0);
      c_t2[m] = complex <double> (0.0, 0.0);
      
      ifstream fin(file[m].c_str());
      string temp_omp;
      
      if(!fin)
      {
        cout << "Can't open " << file[m] << "!"<< endl;
        #pragma omp atomic write
        break_omp = 1;
      continue;
      }
      
      while(getline(fin, temp_omp))
      {
        istringstream iss(temp_omp.c_str());
        
        if(blank_omp == 3 && temp_omp != "")
        { 
          iss >> first_omp;
          iss >> atom[0] >> atom[1] >> atom[2];
          
          if(range[0][0]<=atom[0] && atom[0]<=range[0][1]
          && range[1][0]<=atom[1] && atom[1]<=range[1][1]
          && range[2][0]<=atom[2] && atom[2]<=range[2][1]
          &&(first_omp =="Si"||first_omp =="Ge"||first_omp =="Sn"||first_omp =="C"
           ||first_omp =="Mo"||first_omp =="S"))
          {
            // adjust data considering periodic boundary condition
            
            displace[0] = atom[0] - zero_atom[label].x;
            displace[1] = atom[1] - zero_atom[label].y;
            displace[2] = atom[2] - zero_atom[label].z;
            
            if(displace[0] > 0.5)
            {
              displace[0] = displace[0] - 1.0;
              atom[0] = atom[0] - 1.0;
            }
            
            if(displace[0] < -0.5)
            {
              displace[0] = displace[0] + 1.0;
              atom[0] = atom[0] + 1.0;
            }
            
            if(displace[1] > 0.5)
            {
              displace[1] = displace[1] - 1.0;
              atom[1] = atom[1] - 1.0;
            }
            
            if(displace[1] < -0.5)
            {
              displace[1] = displace[1] + 1.0;
              atom[1] = atom[1] + 1.0;
            }
            
            if(displace[2] > 0.5)
            {
              displace[2] = displace[2] - 1.0;
              atom[2] = atom[2] - 1.0;
            }
            
            if(displace[2] < -0.5)
            {
              displace[2] = displace[2] + 1.0;
              atom[2] = atom[2] + 1.0;
            }
            
            // calculate scattering amplitude

            c_la[m] += (complex<double>) (exp(-i * ((q[0]+dir_la[0]*coe[0])*(atom[0]*leng[0]   +atom[1]*tilt[0][1]+atom[2]*tilt[0][2]) + 
                                                    (q[1]+dir_la[1]*coe[1])*(atom[0]*tilt[1][0]+atom[1]*leng[1]   +atom[2]*tilt[1][2]) + 
                                                    (q[2]+dir_la[2]*coe[2])*(atom[0]*tilt[2][0]+atom[1]*tilt[2][1]+atom[2]*leng[2]))));

            c_lo[m] += (complex<double>) (exp(-i * ((q[0]+dir_lo[0]*coe[0])*(atom[0]*leng[0]   +atom[1]*tilt[0][1]+atom[2]*tilt[0][2]) + 
                                                    (q[1]+dir_lo[1]*coe[1])*(atom[0]*tilt[1][0]+atom[1]*leng[1]   +atom[2]*tilt[1][2]) + 
                                                    (q[2]+dir_lo[2]*coe[2])*(atom[0]*tilt[2][0]+atom[1]*tilt[2][1]+atom[2]*leng[2]))));

            c_t1[m] += (complex<double>) (exp(-i * ((q[0]+dir_t1[0]*coe[0])*(atom[0]*leng[0]   +atom[1]*tilt[0][1]+atom[2]*tilt[0][2]) + 
                                                    (q[1]+dir_t1[1]*coe[1])*(atom[0]*tilt[1][0]+atom[1]*leng[1]   +atom[2]*tilt[1][2]) + 
                                                    (q[2]+dir_t1[2]*coe[2])*(atom[0]*tilt[2][0]+atom[1]*tilt[2][1]+atom[2]*leng[2]))));

            c_t2[m] += (complex<double>) (exp(-i * ((q[0]+dir_t2[0]*coe[0])*(atom[0]*leng[0]   +atom[1]*tilt[0][1]+atom[2]*tilt[0][2]) + 
                                                    (q[1]+dir_t2[1]*coe[1])*(atom[0]*tilt[1][0]+atom[1]*leng[1]   +atom[2]*tilt[1][2]) + 
                                                    (q[2]+dir_t2[2]*coe[2])*(atom[0]*tilt[2][0]+atom[1]*tilt[2][1]+atom[2]*leng[2]))));
            
            label++;

/*
            if(l == 0 && m == 1) cerr << label << " " << leng[0]
                                      << " x " << atom[0] << " y " << atom[1] << " z " << atom[2]
                                      << " re " << real(c_lo[m]) << " im " << imag(c_lo[m]) << " abs " << abs(c_lo[m]) << endl;
*/
          } 
        }
        
        if(temp_omp == "")
          blank_omp++;
        
      }
      
      c_la[m] /= (complex <double>) all_atom;
      c_lo[m] /= (complex <double>) all_atom;
      c_t1[m] /= (complex <double>) all_atom;
      c_t2[m] /= (complex <double>) all_atom;
      
      fin.close();
    }
    
    ostringstream osfile_c_la;    // name of the new output file of longitudinal acoustic mode
    ostringstream osfile_c_lo;    // name of the new output file of longitudinal optical mode
    ostringstream osfile_c_t1;    // name of the new output file of transverse 1 mode
    ostringstream osfile_c_t2;    // name of the new output file of transverse 2 mode
 
    osfile_c_la << "c-" << atom_name << "-" << output << "-" << wave << "-LA-" << l << ext;
    osfile_c_lo << "c-" << atom_name << "-" << output << "-" << wave << "-LO-" << l << ext;
    osfile_c_t1 << "c-" << atom_name << "-" << output << "-" << wave << "-T1-" << l << ext;
    osfile_c_t2 << "c-" << atom_name << "-" << output << "-" << wave << "-T2-" << l << ext;

    ofstream fout_c_la(osfile_c_la.str().c_str());
    ofstream fout_c_lo(osfile_c_lo.str().c_str());
    ofstream fout_c_t1(osfile_c_t1.str().c_str());
    ofstream fout_c_t2(osfile_c_t2.str().c_str());
    
    if(!fout_c_la)
    {
      cout << "Can't open " << osfile_c_la.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_c_lo)
    {
      cout << "Can't open " << osfile_c_lo.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    
    if(!fout_c_t1)
    {
      cout << "Can't open " << osfile_c_t1.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_c_t2)
    {
      cout << "Can't open " << osfile_c_t2.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    
    
    for(int m=0; m<num_data; ++m)
    {
      fout_c_la << m << " " << real(c_la[m]) << " " << imag(c_la[m]) << " " << abs(c_la[m]) << endl;
      fout_c_lo << m << " " << real(c_lo[m]) << " " << imag(c_lo[m]) << " " << abs(c_lo[m]) << endl;
      fout_c_t1 << m << " " << real(c_t1[m]) << " " << imag(c_t1[m]) << " " << abs(c_t1[m]) << endl;
      fout_c_t2 << m << " " << real(c_t2[m]) << " " << imag(c_t2[m]) << " " << abs(c_t2[m]) << endl;
    }
    
    fout_c_la.close();
    fout_c_lo.close();
    fout_c_t1.close();
    fout_c_t2.close();
    
    // fourier transform using MKL(Math Kernel Library)
    
    DFTI_DESCRIPTOR *my_desc1_handle_la, *my_desc1_handle_lo, *my_desc1_handle_t1, *my_desc1_handle_t2;
    
    long int status_la, status_lo, status_t1, status_t2;
    
    status_la = DftiCreateDescriptor (&my_desc1_handle_la, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_lo = DftiCreateDescriptor (&my_desc1_handle_lo, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_t1 = DftiCreateDescriptor (&my_desc1_handle_t1, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_t2 = DftiCreateDescriptor (&my_desc1_handle_t2, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_la = DftiSetValue(my_desc1_handle_la, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_lo = DftiSetValue(my_desc1_handle_lo, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_t1 = DftiSetValue(my_desc1_handle_t1, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_t2 = DftiSetValue(my_desc1_handle_t2, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_la = DftiCommitDescriptor (my_desc1_handle_la);
    status_lo = DftiCommitDescriptor (my_desc1_handle_lo);
    status_t1 = DftiCommitDescriptor (my_desc1_handle_t1);
    status_t2 = DftiCommitDescriptor (my_desc1_handle_t2);
    status_la = DftiComputeForward (my_desc1_handle_la, c_la, s_fourier_la);
    status_lo = DftiComputeForward (my_desc1_handle_la, c_lo, s_fourier_lo);
    status_t1 = DftiComputeForward (my_desc1_handle_t1, c_t1, s_fourier_t1);
    status_t2 = DftiComputeForward (my_desc1_handle_t2, c_t2, s_fourier_t2);  
    status_la = DftiFreeDescriptor (&my_desc1_handle_la);
    status_lo = DftiFreeDescriptor (&my_desc1_handle_lo);
    status_t1 = DftiFreeDescriptor (&my_desc1_handle_t1);
    status_t2 = DftiFreeDescriptor (&my_desc1_handle_t2);
    
    // create new file automatically
    
    ostringstream osfile_s_la;    // name of the new output file of longitudinal acoustic mode
    ostringstream osfile_s_lo;    // name of the new output file of longitudinal optical mode
    ostringstream osfile_s_t1;    // name of the new output file of transverse 1 mode
    ostringstream osfile_s_t2;    // name of the new output file of transverse 2 mode 
    ostringstream osfile_s_all;   // name of the new output file of all mode 

    osfile_s_la << atom_name << "-" << output << "-" << wave << "-LA-" << l << ext;
    osfile_s_lo << atom_name << "-" << output << "-" << wave << "-LO-" << l << ext;
    osfile_s_t1 << atom_name << "-" << output << "-" << wave << "-T1-" << l << ext;
    osfile_s_t2 << atom_name << "-" << output << "-" << wave << "-T2-" << l << ext;
    osfile_s_all << "a-" << atom_name << "-" << output << "-" << wave << "-ALL-" << l << ext;

    // write the wave vector[1/A], frequency[THz] and power spectrum
    
    ofstream fout_s_la(osfile_s_la.str().c_str());
    ofstream fout_s_lo(osfile_s_lo.str().c_str());
    ofstream fout_s_t1(osfile_s_t1.str().c_str());
    ofstream fout_s_t2(osfile_s_t2.str().c_str());
    ofstream fout_s_all(osfile_s_all.str().c_str());
    
    if(!fout_s_la)
    {
      cout << "Can't open " << osfile_s_la.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    
    if(!fout_s_lo)
    {
      cout << "Can't open " << osfile_s_lo.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_s_t1)
    {
      cout << "Can't open " << osfile_s_t1.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_s_t2)
    {
      cout << "Can't open " << osfile_s_t2.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_s_all)
    {
      cout << "Can't open " << osfile_s_all.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    
    for(int m=0; m<num_data; ++m)
    {
      if(tera*(double)(m+1)/(dt*(double)num_data) < max_hertz)
      {
        fout_s_la << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                  << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                  << abs(s_fourier_la[m]) * dt << endl;

        fout_s_lo << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                  << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                  << abs(s_fourier_lo[m]) * dt << endl;

        fout_s_t1 << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                  << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                  << abs(s_fourier_t1[m]) * dt << endl;

        fout_s_t2 << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                  << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                  << abs(s_fourier_t2[m]) * dt << endl;

        fout_s_all << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                  << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                  << abs(s_fourier_la[m]) * dt << " "
                  << abs(s_fourier_lo[m]) * dt << " "
                  << abs(s_fourier_t1[m]) * dt << " "
                  << abs(s_fourier_t2[m]) * dt << endl;

        sum_la[m] += abs(s_fourier_la[m]) * dt;
        sum_lo[m] += abs(s_fourier_lo[m]) * dt;
        sum_t1[m] += abs(s_fourier_t1[m]) * dt;
        sum_ram[m] += sum_lo[m] + sum_t1[m] + sum_t2[m];
      }
    }
    
    fout_s_la.close();
    fout_s_lo.close();
    fout_s_t1.close();
    fout_s_t2.close();
    
  }

  if(break_omp) return 1;
  cout << "XX ..... Finished !" << endl << endl;

  // create sum file

    ostringstream osfile_sum_la;    // name of the new output file of longitudinal acoustic mode
    ostringstream osfile_sum_lo;    // name of the new output file of longitudinal optical mode
    ostringstream osfile_sum_t1;    // name of the new output file of transverse 1 mode 
    ostringstream osfile_sum_t2;    // name of the new output file of transverse 2 mode 
    ostringstream osfile_sum_ram;   // name of the new output file of raman intensity 
    ostringstream osfile_log;       // name of the new output file of log 

    osfile_sum_la << "sum-" << atom_name << "-" << output << "-" << wave << "-LA" << ext;
    osfile_sum_lo << "sum-" << atom_name << "-" << output << "-" << wave << "-LO" << ext;
    osfile_sum_t1 << "sum-" << atom_name << "-" << output << "-" << wave << "-T1" << ext;
    osfile_sum_t2 << "sum-" << atom_name << "-" << output << "-" << wave << "-T2" << ext;
    osfile_sum_ram << "sum-" << atom_name << "-" << output << "-" << wave << "-RAMAN"  << ext;
    osfile_log << "log-" << atom_name << "-" << output << "-" << wave << ext;
    
    ofstream fout_sum_la(osfile_sum_la.str().c_str());
    ofstream fout_sum_lo(osfile_sum_lo.str().c_str());
    ofstream fout_sum_t1(osfile_sum_t1.str().c_str());
    ofstream fout_sum_t2(osfile_sum_t2.str().c_str());
    ofstream fout_sum_ram(osfile_sum_ram.str().c_str());
    ofstream fout_log(osfile_log.str().c_str());
    
    if(!fout_sum_la)
    {
      cout << "Can't open " << osfile_sum_la.str() << "!" << endl;
      return 1;
    }
    
    if(!fout_sum_lo)
    {
      cout << "Can't open " << osfile_sum_lo.str() << "!" << endl;
      return 1;
    }

    if(!fout_sum_t1)
    {
      cout << "Can't open " << osfile_sum_t1.str() << "!" << endl;
      return 1;
    }

    if(!fout_sum_t2)
    {
      cout << "Can't open " << osfile_sum_t2.str() << "!" << endl;
      return 1;
    }

    if(!fout_sum_ram)
    {
      cout << "Can't open " << osfile_sum_ram.str() << "!" << endl;
      return 1;

    if(!fout_log)
    {
      cout << "Can't open " << osfile_log.str() << "!" << endl;
      return 1;
    }
    }
    
    for(int m=0; m<num_data; ++m)
    {
      if(tera*(double)(m+1)/(dt*(double)num_data) < max_hertz)
      {
        fout_sum_la << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_la[m] << endl;

        fout_sum_lo << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_lo[m] << endl;

        fout_sum_t1 << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_t1[m] << endl;

        fout_sum_t2 << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_t2[m] << endl;

        fout_sum_ram << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                     << sum_ram[m] << endl;
      }
    }

    fout_log << "ATMNAME" << endl << atom_name << endl;
    fout_log << "LATTICE_X" << endl << lattice[0] << endl;
    fout_log << "LATTICE_Y" << endl << lattice[1] << endl;
    fout_log << "LATTICE_Z" << endl << lattice[2] << endl;
    fout_log << "MAXFREQ" << endl << max_hertz << endl;
    fout_log << "MAXWAVE" << endl << max_dir << endl;
    
    fout_sum_la.close();
    fout_sum_lo.close();
    fout_sum_t1.close();
    fout_sum_t2.close();
    fout_sum_ram.close();
    fout_log.close();
  
  // free dynamic memory
  
  delete[] zero_atom;
  delete[] file;
  
  return 0;
}

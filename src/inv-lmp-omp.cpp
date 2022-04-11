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
typedef enum {Si, Ge, C, Sn, Mo, S, O, X1, X2} AtomIonType;

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

int inv_lmp_omp(int argc, char* argv[], optpara opts, string atom_conc, int num_omp)
{  
  /***********************************/

  #ifdef _OPENMP
  if(num_omp > 1)
    omp_set_num_threads(num_omp);
  else
    omp_set_num_threads(1);
  printf("set %d OpenMP threads\n",omp_get_max_threads());
  #endif
  
  int part = 8;                 // value to divide wave vector  
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

  double bpath_center = 520;
  double bpath_width = 20;
  double dir_lo[3] = {2.0, 0.0, 0.0};
  double out_dir[3] = {1.0, 0.0, 0.0};
  double boost = 1000.0;
  double qz = 0;
  if(opts.range.length() > 0) {
    istringstream range_iss(opts.range.c_str());
    range_iss >> part >> bpath_center >> bpath_width >> boost >> qz >> dir_lo[0] >> dir_lo[1] >> dir_lo[2] >> out_dir[0] >> out_dir[1] >> out_dir[2];
   }

  cout << "wave " << wave << endl;
  cout << "ldiv " << part << endl;
  cout << "coff " << cutoff << endl;
  
  /***********************************/

  double llo, lhi, ls = 0;
  double leng[3] = {0.0, 0.0, 0.0};    // length of each side of crystal  
  double tilt[3][3] = {{0.0, 0.0, 0.0}, 
                       {0.0, 0.0, 0.0},  // x
                       {0.0, 0.0, 0.0}}; // tilt length of each side of crystal  
  double lattice[3] = {0.0, 0.0, 0.0}; // lattice constant[A]
  double dir[3];                       // each direction index of wave vector
  double coe[3];                       // coefficient of each direction in the wave vector
  double mass[4];
  
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

  int time_step = 0;
  int atom_type;
  int atom_type_O = 0;
  string atom_name = atom_conc;        // atom name(Si or Ge)
  string first;                        // first column of line temporarily
  string temp_line, temp_zero, temp_calc;   // line of string (list/first file/calculation file)
  
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
  
  // read the MD file from the input list 

  ifstream fin(opts.read.c_str());
  
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
    if(num_data == 0) {
      istringstream zero_atom(temp_zero.c_str());
      zero_atom >> all_atom;
    }

    if(all_atom+1 > max_atom)
    {
      cout << "Dynamic memory space is not enough for atom!" << endl;
      return 1;
    }

    getline(fin, temp_zero);
    getline(fin, temp_zero);
    if(num_data == 0) {
      istringstream zero_lat(temp_zero.c_str());
      zero_lat >> llo >> lhi >> ls;
      leng[0] = lhi - llo;
      tilt[0][1] = ls * 0.5;
      tilt[1][0] = ls * 0.5;
    }
    getline(fin, temp_zero);
    if(num_data == 0) {
      istringstream zero_lat(temp_zero.c_str());
      zero_lat >> llo >> lhi >> ls;
      leng[1] = lhi - llo;
      tilt[0][2] = ls * 0.5;
      tilt[2][0] = ls * 0.5;
    }
    getline(fin, temp_zero);
    if(num_data == 0) {
      istringstream zero_lat(temp_zero.c_str());
      zero_lat >> llo >> lhi >> ls;
      leng[2] = lhi - llo;
      tilt[1][2] = ls * 0.5;
      tilt[2][1] = ls * 0.5;
      cout << "|ul00 ul01 ul02| ";
      cout << showpoint << "|" << setw(6) << leng[0]    << " " << tilt[0][1] << " " << tilt[0][2] << "|" << endl;
      cout << "|ul10 ul11 ul12|=";
      cout << showpoint << "|" << setw(6) << tilt[1][0] << " " << leng[1]    << " " << tilt[1][2] << "|" << endl;
      cout << "|ul20 ul21 ul22| ";
      cout << showpoint << "|" << setw(6) << tilt[2][0] << " " << tilt[2][1] << " " << leng[2] << "|" << endl;
    }

    for(int i=0; i<all_atom+1; i++)
      getline(fin, temp_zero);

    num_data++;
  }

  fin.close();
  
  // read the information of the first structure

  fin.open(opts.read.c_str());
  for(int i=0; i<9; i++)
    getline(fin, temp_zero);
  for(int i=0; i<all_atom; i++)
  {
    int label = 0;
    getline(fin, temp_zero);
    istringstream zero_iss(temp_zero.c_str());
    zero_iss >> label >> atom_type;
    zero_iss >> zero_atom[label-1].x >> zero_atom[label-1].y >> zero_atom[label-1].z;

    if(atom_name == "XX") {
      if(atom_type == 1) atomnum[X1]++;
      if(atom_type == 2) atomnum[X2]++;
    }
    else if(atom_name == "MoS2") {
      if(atom_type == 1) atomnum[S]++;
      if(atom_type == 2) atomnum[S]++;
      if(atom_type == 3) atomnum[Mo]++;
    }
    else if(atom_name == "SiGe") {
      if(atom_type == 1) atomnum[Si]++;
      if(atom_type == 2) atomnum[Ge]++;
    }
    else if(atom_name == "GeSn") {
      if(atom_type == 1) atomnum[Ge]++;
      if(atom_type == 2) atomnum[Sn]++;
    }
    else if(atom_name == "XSiO2") {
      if(atom_type == 1) atomnum[X1]++;
      if(atom_type == 2) atomnum[Si]++;
      if(atom_type == 3) atomnum[O]++;
    }
    else if(atom_name == "SiO2") {
      if(atom_type == 1) atomnum[Si]++;
      if(atom_type == 2) atomnum[O]++;
    }
    else if(atom_name == "Si") {
      if(atom_type == 1) atomnum[Si]++;
    }
    else if(atom_name == "Ge") {
      if(atom_type == 1) atomnum[Ge]++;
    }
    else if(atom_name == "Sn") {
      if(atom_type == 1) atomnum[Sn]++;
    }
    else if(atom_name == "C") {
      if(atom_type == 1) atomnum[C]++;
    }
  }

  fin.close();

  if(atomnum[X1] && atomnum[X2])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[Si] || atomnum[Ge] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Custom-XX";
    atomconc[X2] = (double)atomnum[X2] / (atomnum[X2] + atomnum[X1]);
    while(!max_hertz){
      cout << "Please enter a max frequency [THz] (ex. Si : 20 THz)" << endl;
      cin >> max_hertz;
      if(!max_hertz) { cout << "input is wrong!!" << endl; cin.clear(); cin.ignore(); }
    }
    while(!lattice[0]){
      cout << "Please enter a lattice constant of a [A] (ex. Si : 5.431 A)" << endl;
      cin >> lattice[0];
      if(!lattice[0]) { cout << "input is wrong!!" << endl; cin.clear(); cin.ignore(); }
    }
    while(!lattice[1]){
      cout << "Please enter a lattice constant of b [A] (ex. Si : 5.431 A)" << endl;
      cin >> lattice[1];
      if(!lattice[1]) { cout << "input is wrong!!" << endl; cin.clear(); cin.ignore(); }
    }
    while(!lattice[2]){
      cout << "Please enter a lattice constant of c [A] (ex. Si : 5.431 A)" << endl;
      cin >> lattice[2];
      if(!lattice[2]) { cout << "input is wrong!!" << endl; cin.clear(); cin.ignore(); }
    }
  }
  else if(atomnum[Mo] && atomnum[S])
  {
    if(atomnum[Si] || atomnum[Ge] || atomnum[C] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "MoS2";
    max_hertz = 20;
    lattice[0] = 3.127;
    lattice[1] = lattice[0]*cos(-30 * M_PI/180)*2.0;
    lattice[2] = 12.066;
    mass[0] = 32.065;
    mass[1] = 32.065;
    mass[2] = 95.94;
  }
  else if(atomnum[Si] && atomnum[Ge] && atomnum[Sn] && opts.polytype == "3H")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3H-GeSiSn";
    max_hertz = 20;
    atomconc[Si] = (double)atomnum[Si] / (atomnum[Sn] + atomnum[Ge] + atomnum[Si]);
    atomconc[Sn] = (double)atomnum[Sn] / (atomnum[Sn] + atomnum[Ge] + atomnum[Si]);
    double lattice_temp;
    lattice_temp = 5.658 - 0.227*atomconc[Si] + 0.831*atomconc[Sn];
    lattice[0] = lattice_temp * sqrt(3);
    lattice[1] = lattice_temp * sqrt(2) / 2.0;
    lattice[2] = lattice[1];
    mass[0] = 72.63;
    mass[1] = 28.0855;
    mass[2] = 118.710;
  }
  else if(atomnum[Si] && atomnum[Ge] && atomnum[Sn] && opts.polytype == "3T")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3T-GeSiSn";
    max_hertz = 20;
    atomconc[Si] = (double)atomnum[Si] / (atomnum[Sn] + atomnum[Ge] + atomnum[Si]);
    atomconc[Sn] = (double)atomnum[Sn] / (atomnum[Sn] + atomnum[Ge] + atomnum[Si]);
    lattice[2] = 5.658 - 0.227*atomconc[Si] + 0.831*atomconc[Sn];
    lattice[1] = lattice[2] * sqrt(2) / 2.0;
    lattice[0] = lattice[1];
    mass[0] = 72.63;
    mass[1] = 28.0855;
    mass[2] = 118.710;
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
    mass[0] = 72.63;
    mass[1] = 28.0855;
    mass[2] = 118.710;
  }
  else if(atomnum[Si] && atomnum[Ge] && opts.polytype == "3H")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3H-SiGe";
    max_hertz = 20;
    atomconc[Ge] = (double)atomnum[Ge] / (atomnum[Ge] + atomnum[Si]);
    double lattice_temp;
    lattice_temp = 5.431 + 0.2*atomconc[Ge] + 0.027*atomconc[Ge]*atomconc[Ge];
    lattice[0] = lattice_temp * sqrt(3);
    lattice[1] = lattice_temp * sqrt(2) / 2.0;
    lattice[2] = lattice[1];
    mass[0] = 28.0855;
    mass[1] = 72.63;
  }
  else if(atomnum[Si] && atomnum[Ge] && opts.polytype == "3T")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3T-SiGe";
    max_hertz = 20;
    atomconc[Ge] = (double)atomnum[Ge] / (atomnum[Ge] + atomnum[Si]);
    lattice[2] = 5.431 + 0.2*atomconc[Ge] + 0.027*atomconc[Ge]*atomconc[Ge];
    lattice[1] = lattice[2] * sqrt(2) / 2.0;
    lattice[0] = lattice[1];
    mass[0] = 28.0855;
    mass[1] = 72.63;
  }
  else if(atomnum[Si] && atomnum[Ge])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "SiGe";
    max_hertz = 20;
    atomconc[Ge] = (double)atomnum[Ge] / (atomnum[Ge] + atomnum[Si]);
    lattice[0] = 5.431 + 0.2*atomconc[Ge] + 0.027*atomconc[Ge]*atomconc[Ge];
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
    mass[0] = 28.0855;
    mass[1] = 72.63;
  }
  else if(atomnum[Ge] && atomnum[Sn])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Si] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "GeSn";
    max_hertz = 12;
    atomconc[Sn] = (double)atomnum[Sn] / (atomnum[Sn] + atomnum[Ge]);
    lattice[0] = 5.658 + 0.7322*atomconc[Sn] + 0.0988*atomconc[Sn]*atomconc[Sn];
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
    mass[0] = 72.63;
    mass[1] = 118.710;
  }
  else if(atomnum[X1] && atomnum[Si] && atomnum[O])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "XSiO2";
    max_hertz = 20;
    lattice[0] = 5.431;
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
    atom_type_O = 3;
    mass[0] = 28.0855;
    mass[1] = 28.0855;
    mass[2] = 15.9994;
  }
  else if(atomnum[Si] && atomnum[O])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "SiO2";
    max_hertz = 20;
    lattice[0] = 5.431;
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
    atom_type_O = 3;
    mass[0] = 28.0855;
    mass[1] = 15.9994;
  }
  else if(atomnum[Si] && opts.polytype == "3H")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3H-Si";
    max_hertz = 20;
    lattice[0] = 9.40676;
    lattice[1] = 3.84030 / cos(-30 * M_PI/180) / 2.0;
    lattice[2] = 3.84030;
    mass[0] = 28.0855;
  }
  else if(atomnum[Si] && opts.polytype == "3T")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3T-Si";
    max_hertz = 20;
    lattice[0] = 3.84030;
    lattice[1] = lattice[0];
    lattice[2] = 5.431;
    mass[0] = 28.0855;
  }
  else if(atomnum[Si])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Si";
    max_hertz = 20;
    lattice[0] = 5.431;
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
    mass[0] = 28.0855;
  }
  else if(atomnum[Ge])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Ge";
    max_hertz = 12;
    lattice[0] = 5.658; 
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
    mass[0] = 72.63;
  }
  else if(atomnum[Sn])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Si] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Sn";
    max_hertz = 8;
    lattice[0] = 6.489; 
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
    mass[0] = 118.710;
  }
  else if(atomnum[C])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[Ge] || atomnum[Sn] || atomnum[X1] || atomnum[X2]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "C";
    max_hertz = 56;
    lattice[0] = 3.567;
    lattice[1] = lattice[0];
    lattice[2] = lattice[0];
    mass[0] = 12.0107;
  }
  if(max_hertz > tera*0.5/dt)
  {
    cout << "Setting of range for the frequency is wrong!" << endl;
    return 1;
  }
  
  coe[0] = sqrt(norm) * 2.0 * pi / lattice[0];
  coe[1] = sqrt(norm) * 2.0 * pi / lattice[1];
  coe[2] = sqrt(norm) * 2.0 * pi / lattice[2];
  
  cout << "---------------------------------------------" << endl;
  cout << "Atom:" << atom_name << "  Max Frequency:" << max_hertz << "[THz]" << endl;
  cout << "Lattice:" << lattice[0] << " x " << lattice[1] << " x " << lattice[2] << "[A]" << endl;
  cout << "all atom:" << all_atom;
  if(atomnum[Si]) cout << "  Si atom:" << atomnum[Si];
  if(atomnum[Ge]) cout << "  Ge atom:" << atomnum[Ge];
  if(atomnum[Sn]) cout << "  Sn atom:" << atomnum[Sn];
  if(atomnum[C])  cout << "  C atom:"  << atomnum[C];
  if(atomnum[Mo]) cout << "  Mo atom:" << atomnum[Mo];
  if(atomnum[S]) cout  << "  S atom:"  << atomnum[S];
  if(atomnum[O]) cout  << "  O atom:"  << atomnum[O];
  if(atomnum[X1]) cout << "  X1 atom:" << atomnum[X1];
  if(atomnum[X2]) cout << "  X2 atom:" << atomnum[X2];
  cout << endl;
  if(atomconc[Ge]) cout << "Ge conc:" << atomconc[Ge] << endl;
  if(atomconc[Sn]) cout << "Sn conc:" << atomconc[Sn] << endl;
  if(atomconc[X2]) cout << "X2 conc:" << atomconc[X2] << endl;
  cout << "File:" << num_data << "  ";
  cout << "MD time:" << (double)(num_data - 1) * dt / tera << "[ps]  ";
  cout << "Wave vector:[" << wave << "] of " << (double)qz/cutoff << endl;
  cout << "Extruct vector:[" << dir_lo[0] << " " << dir_lo[1] << " " << dir_lo[2] << "]" << endl;
  cout << "Vibration vector:[" << out_dir[0] << " " << out_dir[1] << " " << out_dir[2] << "]" << endl;
  cout << "---------------------------------------------" << endl;
  
  // read the information of all atoms in each file

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
  cout << " Extracting ... atomic trajectory ..." << endl;
  cout << " 0%                  100%" << endl << " ";


    ostringstream* osfile_v = new ostringstream[part+1];    // name of the new output file of longitudinal acoustic mode
    ofstream* fout_v = new ofstream[part+1];
 
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int l=1; l<=part; ++l)
  {
    osfile_v[l] << "v-" << atom_name << "-" << output << "-" << wave << "-" << l << ext;

    fout_v[l].open(osfile_v[l].str().c_str());

    if(!fout_v[l])
    {
      cout << "Can't open " << osfile_v[l].str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
  }
    
    ifstream fin_omp(opts.read.c_str());
    string temp_omp;
  
    for(int m=1; m<num_data; ++m)
    {
      double atom[3] = {0.0, 0.0, 0.0};     // atom information which is read from each file
      double vel[3] = {0.0, 0.0, 0.0};
      int label = 0;                        // number of present label in the input file
      int atom_type_omp;

      for(int i=0; i<9; i++)
        getline(fin_omp, temp_omp);
      for(int i=0; i<all_atom; i++)
      {
        label = 0;
        getline(fin_omp, temp_omp);
        istringstream iss(temp_omp.c_str());
        iss >> label >> atom_type_omp;
        iss >> atom[0] >> atom[1] >> atom[2];
        iss >> vel[0] >> vel[1] >> vel[2];
        if(label<=part) fout_v[label] << label << " " << atom_type_omp << " " << (m-1) * dt / tera << " " << atom[0] << " " << atom[1] << " " << atom[2] << " " << vel[0] << " " << vel[1] << " " << vel[2] << endl;
      }
    }
      fin_omp.close();

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int l=1; l<=part; ++l)
    fout_v[l].close();

  cout << "Finished !" << endl << endl;

  cout << " Calculating ... fourier transform ..." << endl;
  cout << " 0%                  100%" << endl << " ";

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int l=1; l<=part; ++l)
  {
    ostringstream osfile_v;    // name of the new output file of longitudinal acoustic mode
 
    osfile_v << "v-" << atom_name << "-" << output << "-" << wave << "-" << l << ext;

    ifstream fin_v(osfile_v.str().c_str());

    string temp_omp;

    complex <double> j(0.0, 1.0);                // imaginary unit
    double zero_atom[3] = {0.0, 0.0, 0.0};
    double atom[3] = {0.0, 0.0, 0.0};     // atom information which is read from each file
    double vel[3] = {0.0, 0.0, 0.0};
    double displace[3] = {0.0, 0.0, 0.0}; // displacement from equilibirium position
    int atom_type_omp;
    double time;
    int label;

    double q[3] = {dir[0]*coe[0]*qz/cutoff,
                   dir[1]*coe[1]*qz/cutoff,
                   dir[2]*coe[2]*qz/cutoff};        // wave vector[1/A]

    complex <double> c_vx[num_data];
    complex <double> c_vy[num_data];
    complex <double> c_vz[num_data];
    complex <double> s_fourier_cx[num_data];
    complex <double> s_fourier_cy[num_data];
    complex <double> s_fourier_cz[num_data];
    complex <double> inv_vx[num_data];
    complex <double> inv_vy[num_data];
    complex <double> inv_vz[num_data];

    if(part<10) cout << "XX" << flush;
    else if(l%((int)part/10) == 0 ) cout << "XX" << flush;
    
    for(int m=0; m<num_data; ++m)
    {
       getline(fin_v, temp_omp);
       istringstream iss(temp_omp.c_str());
       iss >> label >> atom_type_omp >> time;
       iss >> atom[0] >> atom[1] >> atom[2];
       iss >> vel[0] >> vel[1] >> vel[2];

       if(m == 0){zero_atom[0]=atom[0]; zero_atom[1]=atom[1]; zero_atom[2]=atom[2];}

          displace[0] = atom[0] - zero_atom[0];
          displace[1] = atom[1] - zero_atom[1];
          displace[2] = atom[2] - zero_atom[2];
            
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

          c_vx[m+1] = complex <double> (0.0, 0.0);
          c_vx[m+1] = (complex<double>) displace[0] *
                      (complex<double>) (exp(-j * ((q[0]+dir_lo[0]*coe[0])*(zero_atom[0]*leng[0]   +zero_atom[1]*tilt[0][1]+zero_atom[2]*tilt[0][2]) + 
                                                   (q[1]+dir_lo[1]*coe[1])*(zero_atom[0]*tilt[1][0]+zero_atom[1]*leng[1]   +zero_atom[2]*tilt[1][2]) + 
                                                   (q[2]+dir_lo[2]*coe[2])*(zero_atom[0]*tilt[2][0]+zero_atom[1]*tilt[2][1]+zero_atom[2]*leng[2]))));

          c_vy[m+1] = complex <double> (0.0, 0.0);
          c_vy[m+1] = (complex<double>) displace[1] *
                      (complex<double>) (exp(-j * ((q[0]+dir_lo[0]*coe[0])*(zero_atom[0]*leng[0]   +zero_atom[1]*tilt[0][1]+zero_atom[2]*tilt[0][2]) + 
                                                   (q[1]+dir_lo[1]*coe[1])*(zero_atom[0]*tilt[1][0]+zero_atom[1]*leng[1]   +zero_atom[2]*tilt[1][2]) + 
                                                   (q[2]+dir_lo[2]*coe[2])*(zero_atom[0]*tilt[2][0]+zero_atom[1]*tilt[2][1]+zero_atom[2]*leng[2]))));

          c_vz[m+1] = complex <double> (0.0, 0.0);
          c_vz[m+1] = (complex<double>) displace[2] *
                      (complex<double>) (exp(-j * ((q[0]+dir_lo[0]*coe[0])*(zero_atom[0]*leng[0]   +zero_atom[1]*tilt[0][1]+zero_atom[2]*tilt[0][2]) + 
                                                   (q[1]+dir_lo[1]*coe[1])*(zero_atom[0]*tilt[1][0]+zero_atom[1]*leng[1]   +zero_atom[2]*tilt[1][2]) + 
                                                   (q[2]+dir_lo[2]*coe[2])*(zero_atom[0]*tilt[2][0]+zero_atom[1]*tilt[2][1]+zero_atom[2]*leng[2]))));
    }

    fin_v.close();

    ostringstream osfile_cx;
    ostringstream osfile_cy;
    ostringstream osfile_cz;
    ostringstream osfile_fx;
    ostringstream osfile_fy;
    ostringstream osfile_fz;
    ostringstream osfile_cix;
    ostringstream osfile_ciy;
    ostringstream osfile_ciz;
    ostringstream osfile_i;
 
    osfile_cx << "c-" << atom_name << "-" << output << "-" << wave << "-X-" << l << ext;
    osfile_cy << "c-" << atom_name << "-" << output << "-" << wave << "-Y-" << l << ext;
    osfile_cz << "c-" << atom_name << "-" << output << "-" << wave << "-Z-" << l << ext;
    osfile_fx << "f-" << atom_name << "-" << output << "-" << wave << "-X-" << l << ext;
    osfile_fy << "f-" << atom_name << "-" << output << "-" << wave << "-Y-" << l << ext;
    osfile_fz << "f-" << atom_name << "-" << output << "-" << wave << "-Z-" << l << ext;
    osfile_cix << "ci-" << atom_name << "-" << output << "-" << wave << "-X-" << l << ext;
    osfile_ciy << "ci-" << atom_name << "-" << output << "-" << wave << "-Y-" << l << ext;
    osfile_ciz << "ci-" << atom_name << "-" << output << "-" << wave << "-Z-" << l << ext;
    osfile_i << "vi-" << atom_name << "-" << output << "-" << wave << "-" << l << ext;

    ofstream fout_cx(osfile_cx.str().c_str());
    ofstream fout_cy(osfile_cy.str().c_str());
    ofstream fout_cz(osfile_cz.str().c_str());
    ofstream fout_fx(osfile_fx.str().c_str());
    ofstream fout_fy(osfile_fy.str().c_str());
    ofstream fout_fz(osfile_fz.str().c_str());
    ofstream fout_cix(osfile_cix.str().c_str());
    ofstream fout_ciy(osfile_ciy.str().c_str());
    ofstream fout_ciz(osfile_ciz.str().c_str());
    ofstream fout_i(osfile_i.str().c_str());
    
    if(!fout_cx)
    {
      cout << "Can't open " << osfile_cx.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    if(!fout_cy)
    {
      cout << "Can't open " << osfile_cy.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    if(!fout_cz)
    {
      cout << "Can't open " << osfile_cz.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    if(!fout_fx)
    {
      cout << "Can't open " << osfile_fx.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    if(!fout_fy)
    {
      cout << "Can't open " << osfile_fy.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    if(!fout_fz)
    {
      cout << "Can't open " << osfile_fz.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    if(!fout_cix)
    {
      cout << "Can't open " << osfile_cix.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    if(!fout_ciy)
    {
      cout << "Can't open " << osfile_ciy.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    if(!fout_ciz)
    {
      cout << "Can't open " << osfile_ciz.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    if(!fout_i)
    {
      cout << "Can't open " << osfile_i.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    for(int m=1; m<num_data; ++m)
    {
      fout_cx << l << " " << atom_type_omp << " " << m << " " << (m-1) * dt / tera << " " << real(c_vx[m]) << " " << imag(c_vx[m]) << " " << abs(c_vx[m]) << endl;
      fout_cy << l << " " << atom_type_omp << " " << m << " " << (m-1) * dt / tera << " " << real(c_vy[m]) << " " << imag(c_vy[m]) << " " << abs(c_vy[m]) << endl;
      fout_cz << l << " " << atom_type_omp << " " << m << " " << (m-1) * dt / tera << " " << real(c_vz[m]) << " " << imag(c_vz[m]) << " " << abs(c_vz[m]) << endl;
    }

    fout_cx.close();
    fout_cy.close();
    fout_cz.close();

   // fourier transform using MKL(Math Kernel Library)
    DFTI_DESCRIPTOR *my_desc1_handle_x, *my_desc1_handle_y, *my_desc1_handle_z;
    
    long int status_x, status_y, status_z;
    
    status_x = DftiCreateDescriptor (&my_desc1_handle_x, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_y = DftiCreateDescriptor (&my_desc1_handle_y, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_z = DftiCreateDescriptor (&my_desc1_handle_z, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);

    status_x = DftiSetValue(my_desc1_handle_x, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_y = DftiSetValue(my_desc1_handle_y, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_z = DftiSetValue(my_desc1_handle_z, DFTI_PLACEMENT, DFTI_NOT_INPLACE);

    status_x = DftiCommitDescriptor (my_desc1_handle_x);
    status_y = DftiCommitDescriptor (my_desc1_handle_y);
    status_z = DftiCommitDescriptor (my_desc1_handle_z);
 
    status_x = DftiComputeForward (my_desc1_handle_x, c_vx, s_fourier_cx);
    status_y = DftiComputeForward (my_desc1_handle_y, c_vy, s_fourier_cy);
    status_z = DftiComputeForward (my_desc1_handle_z, c_vz, s_fourier_cz);
 
    status_x = DftiFreeDescriptor (&my_desc1_handle_x);
    status_y = DftiFreeDescriptor (&my_desc1_handle_y);
    status_z = DftiFreeDescriptor (&my_desc1_handle_z);

    for(int m=0; m<num_data; ++m)
    {
      double wf = tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458;
      if(wf < (bpath_center - bpath_width) || wf > (bpath_center + bpath_width))
      {
        s_fourier_cx[m] = complex <double> (0.0, 0.0);
        s_fourier_cy[m] = complex <double> (0.0, 0.0);
        s_fourier_cz[m] = complex <double> (0.0, 0.0);
      }
    }
   
    for(int m=0; m<num_data; ++m)
    {
        fout_fx << l << " " << atom_type_omp << " " << m << " "
                << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                << real(s_fourier_cx[m]) << " " << imag(s_fourier_cx[m]) << " "
                << abs(s_fourier_cx[m]) * abs(s_fourier_cx[m]) / (4.0 * pi * num_data) << endl;
        fout_fy << l << " " << atom_type_omp << " " << m << " "
                << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                << real(s_fourier_cy[m]) << " " << imag(s_fourier_cy[m]) << " "
                << abs(s_fourier_cy[m]) * abs(s_fourier_cy[m]) / (4.0 * pi * num_data) << endl;
        fout_fz << l << " " << atom_type_omp << " " << m << " "
                << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                << real(s_fourier_cz[m]) << " " << imag(s_fourier_cz[m]) << " "
                << abs(s_fourier_cz[m]) * abs(s_fourier_cz[m]) / (4.0 * pi * num_data) << endl;
    }

    fout_fx.close();
    fout_fy.close();
    fout_fz.close();

    status_x = DftiCreateDescriptor (&my_desc1_handle_x, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_y = DftiCreateDescriptor (&my_desc1_handle_y, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_z = DftiCreateDescriptor (&my_desc1_handle_z, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);

    status_x = DftiSetValue(my_desc1_handle_x, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_y = DftiSetValue(my_desc1_handle_y, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_z = DftiSetValue(my_desc1_handle_z, DFTI_PLACEMENT, DFTI_NOT_INPLACE);

    status_x = DftiCommitDescriptor (my_desc1_handle_x);
    status_y = DftiCommitDescriptor (my_desc1_handle_y);
    status_z = DftiCommitDescriptor (my_desc1_handle_z);
 
    status_x = DftiComputeBackward (my_desc1_handle_x, s_fourier_cx, inv_vx);
    status_y = DftiComputeBackward (my_desc1_handle_y, s_fourier_cy, inv_vy);
    status_z = DftiComputeBackward (my_desc1_handle_z, s_fourier_cz, inv_vz);
 
    status_x = DftiFreeDescriptor (&my_desc1_handle_x);
    status_y = DftiFreeDescriptor (&my_desc1_handle_y);
    status_z = DftiFreeDescriptor (&my_desc1_handle_z);

    for(int m=1; m<num_data; ++m)
    {
          inv_vx[m] /= complex <double> (boost, 0.0);
          inv_vy[m] /= complex <double> (boost, 0.0);
          inv_vz[m] /= complex <double> (boost, 0.0);
      fout_cix << l << " " << atom_type_omp << " " << m << " " << (m-1) * dt / tera << " " << real(inv_vx[m]) << " " << imag(inv_vx[m]) << " " << abs(inv_vx[m]) << endl;
      fout_ciy << l << " " << atom_type_omp << " " << m << " " << (m-1) * dt / tera << " " << real(inv_vy[m]) << " " << imag(inv_vy[m]) << " " << abs(inv_vy[m]) << endl;
      fout_ciz << l << " " << atom_type_omp << " " << m << " " << (m-1) * dt / tera << " " << real(inv_vz[m]) << " " << imag(inv_vz[m]) << " " << abs(inv_vz[m]) << endl;
    }

    fout_cix.close();
    fout_ciy.close();
    fout_ciz.close();

    for(int m=1; m<num_data; ++m)
    {
          inv_vx[m] /= (complex<double>) (exp(-j * ((q[0]+dir_lo[0]*coe[0])*(zero_atom[0]*leng[0]   +zero_atom[1]*tilt[0][1]+zero_atom[2]*tilt[0][2]) + 
                                                    (q[1]+dir_lo[1]*coe[1])*(zero_atom[0]*tilt[1][0]+zero_atom[1]*leng[1]   +zero_atom[2]*tilt[1][2]) + 
                                                    (q[2]+dir_lo[2]*coe[2])*(zero_atom[0]*tilt[2][0]+zero_atom[1]*tilt[2][1]+zero_atom[2]*leng[2]))));
          inv_vx[m] *= (complex<double>)out_dir[0];
          inv_vx[m] += (complex<double>)zero_atom[0];

          inv_vy[m] /= (complex<double>) (exp(-j * ((q[0]+dir_lo[0]*coe[0])*(zero_atom[0]*leng[0]   +zero_atom[1]*tilt[0][1]+zero_atom[2]*tilt[0][2]) + 
                                                    (q[1]+dir_lo[1]*coe[1])*(zero_atom[0]*tilt[1][0]+zero_atom[1]*leng[1]   +zero_atom[2]*tilt[1][2]) + 
                                                    (q[2]+dir_lo[2]*coe[2])*(zero_atom[0]*tilt[2][0]+zero_atom[1]*tilt[2][1]+zero_atom[2]*leng[2]))));
          inv_vy[m] *= (complex<double>)out_dir[1];
          inv_vy[m] += (complex<double>)zero_atom[1];

          inv_vz[m] /= (complex<double>) (exp(-j * ((q[0]+dir_lo[0]*coe[0])*(zero_atom[0]*leng[0]   +zero_atom[1]*tilt[0][1]+zero_atom[2]*tilt[0][2]) + 
                                                    (q[1]+dir_lo[1]*coe[1])*(zero_atom[0]*tilt[1][0]+zero_atom[1]*leng[1]   +zero_atom[2]*tilt[1][2]) + 
                                                    (q[2]+dir_lo[2]*coe[2])*(zero_atom[0]*tilt[2][0]+zero_atom[1]*tilt[2][1]+zero_atom[2]*leng[2]))));
          inv_vz[m] *= (complex<double>)out_dir[2];
          inv_vz[m] += (complex<double>)zero_atom[2];

      fout_i << l << " " << atom_type_omp << " " << (m-1) * dt / tera << " " << real(inv_vx[m]) << " " << real(inv_vy[m]) << " " << real(inv_vz[m]) << endl;
    }

    fout_i.close();
  }

  if(break_omp) return 1;
  cout << "XX ..... Finished !" << endl << endl;

  cout << " Output ... XYZ file ..." << endl;
  cout << " 0%                  100%" << endl << " ";

    ostringstream osfile_xyz;
 
    osfile_xyz << atom_name << "-" << output << "-" << wave << ".xyz";

    ofstream fout_xyz(osfile_xyz.str().c_str());

    if(!fout_xyz)
    {
      cout << "Can't open " << osfile_xyz.str() << "!" << endl;
    }

  for(int m=1; m<num_data; ++m)
  {
    if(num_data<10) cout << "XX" << flush;
    else if(m%((int)num_data/10) == 0 ) cout << "XX" << flush;

    fout_xyz << part << endl;
    fout_xyz << "Frequency (cm^-1) = " << bpath_center << " Time step = " << m << endl;

    for(int l=1; l<=part; ++l)
    {
      ostringstream osfile_v;    // name of the new output file of longitudinal acoustic mode
 
      osfile_v << "vi-" << atom_name << "-" << output << "-" << wave << "-" << l << ext;

      ifstream fin_v(osfile_v.str().c_str());

      double atom[3] = {0.0, 0.0, 0.0};     // atom information which is read from each file
      int label = 0;                        // number of present label in the input file
      int atom_type_omp;
      double time;
      string temp_omp;

      for(int i=0; i<m-1; i++)
        getline(fin_v, temp_omp);

      label = 0;
      getline(fin_v, temp_omp);
      istringstream iss(temp_omp.c_str());
      iss >> label >> atom_type_omp >> time;
      iss >> atom[0] >> atom[1] >> atom[2];

          if(atom[0] > 0.5 && l < all_atom/2.0)
          {
            atom[0] = atom[0] - 1.0;
          }
          if(atom[1] > 0.5 && l < all_atom/2.0)
          {
            atom[1] = atom[1] - 1.0;
          }
          if(atom[2] > 0.5 && l < all_atom/2.0)
          {
            atom[2] = atom[2] - 1.0;
          }

    string atom_type_name = "Si";

    if(atom_name == "XX") {
      if(atom_type_omp == 1) atom_type_name = "X1";
      if(atom_type_omp == 2) atom_type_name = "X1";
    }
    else if(atom_name == "MoS2") {
      if(atom_type_omp == 1) atom_type_name = "S ";
      if(atom_type_omp == 2) atom_type_name = "S ";
      if(atom_type_omp == 3) atom_type_name = "Mo";
    }
    else if(atom_name == "SiGe") {
      if(atom_type_omp == 1) atom_type_name = "Si";
      if(atom_type_omp == 2) atom_type_name = "Ge";
    }
    else if(atom_name == "GeSn") {
      if(atom_type_omp == 1) atom_type_name = "Ge";
      if(atom_type_omp == 2) atom_type_name = "Sn";
    }
    else if(atom_name == "XSiO2") {
      if(atom_type_omp == 1) atom_type_name = "X1";
      if(atom_type_omp == 2) atom_type_name = "Si";
      if(atom_type_omp == 3) atom_type_name = "O ";
    }
    else if(atom_name == "SiO2") {
      if(atom_type_omp == 1) atom_type_name = "Si";
      if(atom_type_omp == 2) atom_type_name = "O ";
    }
    else if(atom_name == "Si") {
      if(atom_type_omp == 1) atom_type_name = "Si";
    }
    else if(atom_name == "Ge") {
      if(atom_type_omp == 1) atom_type_name = "Ge";
    }
    else if(atom_name == "Sn") {
      if(atom_type_omp == 1) atom_type_name = "Sn";
    }
    else if(atom_name == "C") {
      if(atom_type_omp == 1) atom_type_name = "C";
    }

      fout_xyz << atom_type_name << " " << atom[0]*leng[0] << " " << atom[1]*leng[1] << " " << atom[2]*leng[2] << endl;

      fin_v.close();
    }
  }

    fout_xyz.close();

  if(break_omp) return 1;
  cout << "XX ..... Finished !" << endl << endl;

    ostringstream osfile_log;       // name of the new output file of log 

    osfile_log << "log-" << atom_name << "-" << output << "-" << wave << ext;
    
    ofstream fout_log(osfile_log.str().c_str());

    if(!fout_log)
    {
      cout << "Can't open " << osfile_log.str() << "!" << endl;
      return 1;
    }

    fout_log << "ATMNAME" << endl << atom_name << endl;
    fout_log << "LATTICE_X" << endl << lattice[0] << endl;
    fout_log << "LATTICE_Y" << endl << lattice[1] << endl;
    fout_log << "LATTICE_Z" << endl << lattice[2] << endl;
    fout_log << "MAXFREQ" << endl << max_hertz << endl;
    fout_log << "MAXWAVE" << endl << max_dir << endl;
  
  // free dynamic memory
  
  delete[] zero_atom;
  delete[] file;
  
  return 0;
}

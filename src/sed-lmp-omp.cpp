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

int sed_lmp_omp(int argc, char* argv[], optpara opts, string atom_conc, int num_omp)
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
  const int max_file = 40005;           // maximum number of file 
  
  double dt = 1e-14;             // time step[s]  
  const double pi = 3.14159265;        // pi
  const double tera = 1e-12;           // convert unit into tera order

  int cutoff = 1;          // value to cutoff wave vector
  
  string wave = "100";                 // direction index of wave vector
  string output = "sed";               // base name of the output file
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
  cout << "sed calculation" << endl;
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
    
    if(dir[l] != 0.0 && dir[l] != 1.0 && dir[l] != 2.0 && dir[l] != 4.0)
    {
      cout << "Setting of wave vector direction is wrong!" << endl;
      return 1;
    }
  }

  double norm = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
  dir[0] /= sqrt(norm); dir[1] /= sqrt(norm); dir[2] /= sqrt(norm);

  double dir_la[3];                    // each direction index of LA phonon
  double dir_lo[3];                    // each direction index of LO phonon

  if(wave == "100" && opts.polytype == "3H") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 3.0; dir_lo[1] = 0.0; dir_lo[2] = 0.0;
  }
  else if(wave == "100" && opts.polytype == "3T") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 0.0; dir_lo[1] = 0.0; dir_lo[2] = 2.0;
  }
  else if(wave == "100") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 2.0; dir_lo[1] = 0.0; dir_lo[2] = 0.0;
  }
  else if(wave == "400") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 0.0; dir_lo[1] = 0.0; dir_lo[2] = 0.0;
    max_dir = 4;
  }
  else if(wave == "200" && opts.polytype == "3H") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 3.0; dir_lo[1] = 0.0; dir_lo[2] = 0.0;
  }
  else if(wave == "200") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 2.0; dir_lo[1] = 0.0; dir_lo[2] = 0.0;
  }
  else if(wave == "010") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 0.0; dir_lo[1] = 2.0; dir_lo[2] = 0.0;
  }
  else if(wave == "001") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 0.0; dir_lo[1] = 0.0; dir_lo[2] = 2.0;
  }
  else if(wave == "110") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 2.0; dir_lo[1] = 2.0; dir_lo[2] = 0.0;
  }
  else if(wave == "220") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 1.4; dir_lo[1] =-1.4; dir_lo[2] = 0.0;
  }
  else if(wave == "111") {
    dir_la[0] = 0.0; dir_la[1] = 0.0; dir_la[2] = 0.0;
    dir_lo[0] = 1.4; dir_lo[1] =-0.8; dir_lo[2] = 1.15;
  }
  else{
      cout << "Setting of wave vector direction is wrong!" << endl;
      return 1;
  }
  
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
      cout << num_data+1 << " steps over than "  << max_file << " steps!" << endl;
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
      cout << all_atom+1 << " atoms over than "  << max_atom << " atoms!" << endl;
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
    else if(atom_name == "Si-of-SiGe") {
      if(atom_type == 1) atomnum[Si]++;
      if(atom_type == 2) atomnum[Ge]++;
    }
    else if(atom_name == "Ge-of-SiGe") {
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
    if(atom_name == "Si-of-SiGe") atom_type_O = 2;
    if(atom_name == "Ge-of-SiGe") atom_type_O = 1;
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
  cout << "Wave vector:[" << wave << "] from 0 to " << (double)1.0/cutoff << endl;
  if(range_on) cout << "Input range: " << range[0][0] << " " << range[0][1] << " xlo xhi "
                                       << range[1][0] << " " << range[1][1] << " ylo yhi "
                                       << range[2][0] << " " << range[2][1] << " zlo zhi" << endl;
  cout << "---------------------------------------------" << endl;
  
  // read the information of all atoms in each file

  double sum_la_x[num_data];
  double sum_la_y[num_data];
  double sum_la_z[num_data];
  double sum_lo_x[num_data];
  double sum_lo_y[num_data];
  double sum_lo_z[num_data];
  double sum_la_all[num_data];
  double sum_lo_all[num_data];
  double sum_all[num_data];

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int m=0; m<num_data; ++m)
  {
    sum_la_x[m] = 0;
    sum_la_y[m] = 0;
    sum_la_z[m] = 0;
    sum_lo_x[m] = 0;
    sum_lo_y[m] = 0;
    sum_lo_z[m] = 0;
    sum_la_all[m] = 0;
    sum_lo_all[m] = 0;
    sum_all[m] = 0;
  }

  int break_omp = 0;

  cout << "Calc time: 0 - " << part << " / " << omp_get_max_threads() << " threads" << endl;  
  if(wave == "220")    
  cout << "Wave vector:("
       << dir[0]*coe[0]/cutoff/2.0 << ", " << dir[1]*coe[1]/cutoff/2.0 << ", " << dir[2]*coe[2]/cutoff/2.0 << ") - (0.00000, 0.00000, 0.00000)" << endl << endl; 
  else if(wave == "111") 
  cout << "Wave vector:(0.00000, 0.00000, 0.00000) - ("
       << dir[0]*coe[0]/cutoff/2.0 << ", " << dir[1]*coe[1]/cutoff/2.0 << ", " << dir[2]*coe[2]/cutoff/2.0 << ")" << endl << endl; 
  else if(wave == "400") 
  cout << "Wave vector:(0.00000, 0.00000, 0.00000) - ("
       << dir[0]*coe[0]/cutoff*4.0 << ", " << dir[1]*coe[1]/cutoff*4.0 << ", " << dir[2]*coe[2]/cutoff*4.0 << ")" << endl << endl;  
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
    double z = l;
    if(wave == "220") z = (part - l)/2.0;
    else if(wave == "111") z = l/2.0;  
    double q[3] = {dir[0]*coe[0]*z/part/cutoff,
                   dir[1]*coe[1]*z/part/cutoff,
                   dir[2]*coe[2]*z/part/cutoff};        // wave vector[1/A]
    
    complex <double> j(0.0, 1.0);                // imaginary unit
    
    complex <double> c_la_x[num_data];           // scattering amplitude of acoustic mode
    complex <double> c_la_y[num_data];           // scattering amplitude of acoustic mode
    complex <double> c_la_z[num_data];           // scattering amplitude of acoustic mode
    complex <double> c_lo_x[num_data];           // scattering amplitude of optical mode
    complex <double> c_lo_y[num_data];           // scattering amplitude of optical mode
    complex <double> c_lo_z[num_data];           // scattering amplitude of optical mode
    complex <double> s_fourier_la_x[num_data];   // dynamic scattering function of acoustic mode
    complex <double> s_fourier_la_y[num_data];   // dynamic scattering function of acoustic mode
    complex <double> s_fourier_la_z[num_data];   // dynamic scattering function of acoustic mode
    complex <double> s_fourier_lo_x[num_data];   // dynamic scattering function of optical mode
    complex <double> s_fourier_lo_y[num_data];   // dynamic scattering function of optical mode
    complex <double> s_fourier_lo_z[num_data];   // dynamic scattering function of optical mode
    
    if(l%((int)part/10) == 0 ) cout << "XX" << flush;

    ifstream fin_omp(opts.read.c_str());
    string temp_omp;
    for(int i=0; i<all_atom+9; i++)
      getline(fin_omp, temp_omp);
    
    for(int m=1; m<num_data; ++m)
    {
      int stop_omp;
      #pragma omp atomic read
      stop_omp = break_omp;
      if (stop_omp) continue;

      double atom[3] = {0.0, 0.0, 0.0};     // atom information which is read from each file
      double displace[3] = {0.0, 0.0, 0.0}; // displacement from equilibirium position
      double vel[3] = {0.0, 0.0, 0.0};
      
      int label = 0;                        // number of present label in the input file
      int atom_type_omp;
      c_la_x[m] = complex <double> (0.0, 0.0);
      c_la_y[m] = complex <double> (0.0, 0.0);
      c_la_z[m] = complex <double> (0.0, 0.0);
      c_lo_x[m] = complex <double> (0.0, 0.0);
      c_lo_y[m] = complex <double> (0.0, 0.0);
      c_lo_z[m] = complex <double> (0.0, 0.0);

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

        if(range[0][0]<=atom[0] && atom[0]<=range[0][1]
        && range[1][0]<=atom[1] && atom[1]<=range[1][1]
        && range[2][0]<=atom[2] && atom[2]<=range[2][1]
        && atom_type_omp != atom_type_O)
        {
            
          // calculate scattering amplitude

          c_la_x[m] += (complex<double>) vel[0] * mass[atom_type_omp-1] *
                       (complex<double>) (exp(-j * ((q[0]+dir_la[0]*coe[0])*(zero_atom[label-1].x*leng[0]   +zero_atom[label-1].y*tilt[0][1]+zero_atom[label-1].z*tilt[0][2]) + 
                                                    (q[1]+dir_la[1]*coe[1])*(zero_atom[label-1].x*tilt[1][0]+zero_atom[label-1].y*leng[1]   +zero_atom[label-1].z*tilt[1][2]) + 
                                                    (q[2]+dir_la[2]*coe[2])*(zero_atom[label-1].x*tilt[2][0]+zero_atom[label-1].y*tilt[2][1]+zero_atom[label-1].z*leng[2]))));

          c_la_y[m] += (complex<double>) vel[1] * mass[atom_type_omp-1] *
                       (complex<double>) (exp(-j * ((q[0]+dir_la[0]*coe[0])*(zero_atom[label-1].x*leng[0]   +zero_atom[label-1].y*tilt[0][1]+zero_atom[label-1].z*tilt[0][2]) + 
                                                    (q[1]+dir_la[1]*coe[1])*(zero_atom[label-1].x*tilt[1][0]+zero_atom[label-1].y*leng[1]   +zero_atom[label-1].z*tilt[1][2]) + 
                                                    (q[2]+dir_la[2]*coe[2])*(zero_atom[label-1].x*tilt[2][0]+zero_atom[label-1].y*tilt[2][1]+zero_atom[label-1].z*leng[2]))));

          c_la_z[m] += (complex<double>) vel[2] * mass[atom_type_omp-1] *
                       (complex<double>) (exp(-j * ((q[0]+dir_la[0]*coe[0])*(zero_atom[label-1].x*leng[0]   +zero_atom[label-1].y*tilt[0][1]+zero_atom[label-1].z*tilt[0][2]) + 
                                                    (q[1]+dir_la[1]*coe[1])*(zero_atom[label-1].x*tilt[1][0]+zero_atom[label-1].y*leng[1]   +zero_atom[label-1].z*tilt[1][2]) + 
                                                    (q[2]+dir_la[2]*coe[2])*(zero_atom[label-1].x*tilt[2][0]+zero_atom[label-1].y*tilt[2][1]+zero_atom[label-1].z*leng[2]))));

          c_lo_x[m] += (complex<double>) vel[0] * mass[atom_type_omp-1] *
                       (complex<double>) (exp(-j * ((q[0]+dir_lo[0]*coe[0])*(zero_atom[label-1].x*leng[0]   +zero_atom[label-1].y*tilt[0][1]+zero_atom[label-1].z*tilt[0][2]) + 
                                                    (q[1]+dir_lo[1]*coe[1])*(zero_atom[label-1].x*tilt[1][0]+zero_atom[label-1].y*leng[1]   +zero_atom[label-1].z*tilt[1][2]) + 
                                                    (q[2]+dir_lo[2]*coe[2])*(zero_atom[label-1].x*tilt[2][0]+zero_atom[label-1].y*tilt[2][1]+zero_atom[label-1].z*leng[2]))));

          c_lo_y[m] += (complex<double>) vel[1] * mass[atom_type_omp-1] *
                       (complex<double>) (exp(-j * ((q[0]+dir_lo[0]*coe[0])*(zero_atom[label-1].x*leng[0]   +zero_atom[label-1].y*tilt[0][1]+zero_atom[label-1].z*tilt[0][2]) + 
                                                    (q[1]+dir_lo[1]*coe[1])*(zero_atom[label-1].x*tilt[1][0]+zero_atom[label-1].y*leng[1]   +zero_atom[label-1].z*tilt[1][2]) + 
                                                    (q[2]+dir_lo[2]*coe[2])*(zero_atom[label-1].x*tilt[2][0]+zero_atom[label-1].y*tilt[2][1]+zero_atom[label-1].z*leng[2]))));

          c_lo_z[m] += (complex<double>) vel[2] * mass[atom_type_omp-1] *
                       (complex<double>) (exp(-j * ((q[0]+dir_lo[0]*coe[0])*(zero_atom[label-1].x*leng[0]   +zero_atom[label-1].y*tilt[0][1]+zero_atom[label-1].z*tilt[0][2]) + 
                                                    (q[1]+dir_lo[1]*coe[1])*(zero_atom[label-1].x*tilt[1][0]+zero_atom[label-1].y*leng[1]   +zero_atom[label-1].z*tilt[1][2]) + 
                                                    (q[2]+dir_lo[2]*coe[2])*(zero_atom[label-1].x*tilt[2][0]+zero_atom[label-1].y*tilt[2][1]+zero_atom[label-1].z*leng[2]))));
        }
      }
      c_la_x[m] /= (complex <double>) all_atom;
      c_la_y[m] /= (complex <double>) all_atom;
      c_la_z[m] /= (complex <double>) all_atom;
      c_lo_x[m] /= (complex <double>) all_atom;
      c_lo_y[m] /= (complex <double>) all_atom;
      c_lo_z[m] /= (complex <double>) all_atom;
    }

    fin.close();

    ostringstream osfile_c_la_x;    // name of the new output file of acoustic mode
    ostringstream osfile_c_la_y;    // name of the new output file of acoustic mode
    ostringstream osfile_c_la_z;    // name of the new output file of acoustic mode
    ostringstream osfile_c_lo_x;    // name of the new output file of optical mode
    ostringstream osfile_c_lo_y;    // name of the new output file of optical mode
    ostringstream osfile_c_lo_z;    // name of the new output file of optical mode
 
    osfile_c_la_x << "c-" << atom_name << "-" << output << "-" << wave << "-AX-" << l << ext;
    osfile_c_la_y << "c-" << atom_name << "-" << output << "-" << wave << "-AY-" << l << ext;
    osfile_c_la_z << "c-" << atom_name << "-" << output << "-" << wave << "-AZ-" << l << ext;
    osfile_c_lo_x << "c-" << atom_name << "-" << output << "-" << wave << "-OX-" << l << ext;
    osfile_c_lo_y << "c-" << atom_name << "-" << output << "-" << wave << "-OY-" << l << ext;
    osfile_c_lo_z << "c-" << atom_name << "-" << output << "-" << wave << "-OZ-" << l << ext;

    ofstream fout_c_la_x(osfile_c_la_x.str().c_str());
    ofstream fout_c_la_y(osfile_c_la_y.str().c_str());
    ofstream fout_c_la_z(osfile_c_la_z.str().c_str());
    ofstream fout_c_lo_x(osfile_c_lo_x.str().c_str());
    ofstream fout_c_lo_y(osfile_c_lo_y.str().c_str());
    ofstream fout_c_lo_z(osfile_c_lo_z.str().c_str());
    
    if(!fout_c_la_x)
    {
      cout << "Can't open " << osfile_c_la_x.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_c_la_y)
    {
      cout << "Can't open " << osfile_c_la_y.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_c_la_z)
    {
      cout << "Can't open " << osfile_c_la_z.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_c_lo_x)
    {
      cout << "Can't open " << osfile_c_lo_x.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_c_lo_y)
    {
      cout << "Can't open " << osfile_c_lo_y.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_c_lo_z)
    {
      cout << "Can't open " << osfile_c_lo_z.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    
    for(int m=0; m<num_data; ++m)
    {
      fout_c_la_x << m << " " << real(c_la_x[m]) << " " << imag(c_la_x[m]) << " " << abs(c_la_x[m]) << endl;
      fout_c_la_y << m << " " << real(c_la_y[m]) << " " << imag(c_la_y[m]) << " " << abs(c_la_y[m]) << endl;
      fout_c_la_z << m << " " << real(c_la_z[m]) << " " << imag(c_la_z[m]) << " " << abs(c_la_z[m]) << endl;
      fout_c_lo_x << m << " " << real(c_lo_x[m]) << " " << imag(c_lo_x[m]) << " " << abs(c_lo_x[m]) << endl;
      fout_c_lo_y << m << " " << real(c_lo_y[m]) << " " << imag(c_lo_y[m]) << " " << abs(c_lo_y[m]) << endl;
      fout_c_lo_z << m << " " << real(c_lo_z[m]) << " " << imag(c_lo_z[m]) << " " << abs(c_lo_z[m]) << endl;
    }
    
    fout_c_la_x.close();
    fout_c_la_y.close();
    fout_c_la_z.close();
    fout_c_lo_x.close();
    fout_c_lo_y.close();
    fout_c_lo_z.close();
    
    // fourier transform using MKL(Math Kernel Library)
    
    DFTI_DESCRIPTOR *my_desc1_handle_la_x, *my_desc1_handle_la_y, *my_desc1_handle_la_z, *my_desc1_handle_lo_x, *my_desc1_handle_lo_y, *my_desc1_handle_lo_z;
    
    long int status_la_x, status_la_y, status_la_z, status_lo_x, status_lo_y, status_lo_z;
    
    status_la_x = DftiCreateDescriptor (&my_desc1_handle_la_x, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_la_y = DftiCreateDescriptor (&my_desc1_handle_la_y, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_la_z = DftiCreateDescriptor (&my_desc1_handle_la_z, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_lo_x = DftiCreateDescriptor (&my_desc1_handle_lo_x, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_lo_y = DftiCreateDescriptor (&my_desc1_handle_lo_y, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_lo_z = DftiCreateDescriptor (&my_desc1_handle_lo_z, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
    status_la_x = DftiSetValue(my_desc1_handle_la_x, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_la_y = DftiSetValue(my_desc1_handle_la_y, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_la_z = DftiSetValue(my_desc1_handle_la_z, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_lo_x = DftiSetValue(my_desc1_handle_lo_x, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_lo_y = DftiSetValue(my_desc1_handle_lo_y, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_lo_z = DftiSetValue(my_desc1_handle_lo_z, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status_la_x = DftiCommitDescriptor (my_desc1_handle_la_x);
    status_la_y = DftiCommitDescriptor (my_desc1_handle_la_y);
    status_la_z = DftiCommitDescriptor (my_desc1_handle_la_z);
    status_lo_x = DftiCommitDescriptor (my_desc1_handle_lo_x);
    status_lo_y = DftiCommitDescriptor (my_desc1_handle_lo_y);
    status_lo_z = DftiCommitDescriptor (my_desc1_handle_lo_z);
    status_la_x = DftiComputeForward (my_desc1_handle_la_x, c_la_x, s_fourier_la_x);
    status_la_y = DftiComputeForward (my_desc1_handle_la_y, c_la_y, s_fourier_la_y);
    status_la_z = DftiComputeForward (my_desc1_handle_la_z, c_la_z, s_fourier_la_z);
    status_lo_x = DftiComputeForward (my_desc1_handle_lo_x, c_lo_x, s_fourier_lo_x);
    status_lo_y = DftiComputeForward (my_desc1_handle_lo_y, c_lo_y, s_fourier_lo_y);
    status_lo_z = DftiComputeForward (my_desc1_handle_lo_z, c_lo_z, s_fourier_lo_z);
    status_la_x = DftiFreeDescriptor (&my_desc1_handle_la_x);
    status_la_y = DftiFreeDescriptor (&my_desc1_handle_la_y);
    status_la_z = DftiFreeDescriptor (&my_desc1_handle_la_z);
    status_lo_x = DftiFreeDescriptor (&my_desc1_handle_lo_x);
    status_lo_y = DftiFreeDescriptor (&my_desc1_handle_lo_y);
    status_lo_z = DftiFreeDescriptor (&my_desc1_handle_lo_z);
    
    // create new file automatically
    
    ostringstream osfile_s_la_x;    // name of the new output file of acoustic mode
    ostringstream osfile_s_la_y;    // name of the new output file of acoustic mode
    ostringstream osfile_s_la_z;    // name of the new output file of acoustic mode
    ostringstream osfile_s_lo_x;    // name of the new output file of optical mode
    ostringstream osfile_s_lo_y;    // name of the new output file of optical mode
    ostringstream osfile_s_lo_z;    // name of the new output file of optical mode

    osfile_s_la_x << atom_name << "-" << output << "-" << wave << "-AX-" << l << ext;
    osfile_s_la_y << atom_name << "-" << output << "-" << wave << "-AY-" << l << ext;
    osfile_s_la_z << atom_name << "-" << output << "-" << wave << "-AZ-" << l << ext;
    osfile_s_lo_x << atom_name << "-" << output << "-" << wave << "-OX-" << l << ext;
    osfile_s_lo_y << atom_name << "-" << output << "-" << wave << "-OY-" << l << ext;
    osfile_s_lo_z << atom_name << "-" << output << "-" << wave << "-OZ-" << l << ext;

    // write the wave vector[1/A], frequency[THz] and power spectrum
    
    ofstream fout_s_la_x(osfile_s_la_x.str().c_str());
    ofstream fout_s_la_y(osfile_s_la_y.str().c_str());
    ofstream fout_s_la_z(osfile_s_la_z.str().c_str());
    ofstream fout_s_lo_x(osfile_s_lo_x.str().c_str());
    ofstream fout_s_lo_y(osfile_s_lo_y.str().c_str());
    ofstream fout_s_lo_z(osfile_s_lo_z.str().c_str());
    
    if(!fout_s_la_x)
    {
      cout << "Can't open " << osfile_s_la_x.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_s_la_y)
    {
      cout << "Can't open " << osfile_s_la_y.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_s_la_z)
    {
      cout << "Can't open " << osfile_s_la_z.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    
    if(!fout_s_lo_x)
    {
      cout << "Can't open " << osfile_s_lo_x.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_s_lo_y)
    {
      cout << "Can't open " << osfile_s_lo_y.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }

    if(!fout_s_lo_z)
    {
      cout << "Can't open " << osfile_s_lo_z.str() << "!" << endl;
      #pragma omp atomic write
      break_omp = 1;
      continue;
    }
    
    for(int m=0; m<num_data; ++m)
    {
      if(tera*(double)(m+1)/(dt*(double)num_data) < max_hertz)
      {
        fout_s_la_x << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                    << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << abs(s_fourier_la_x[m]) * abs(s_fourier_la_x[m]) / (4.0 * pi * num_data) << endl;

        fout_s_la_y << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                    << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << abs(s_fourier_la_y[m]) * abs(s_fourier_la_y[m]) / (4.0 * pi * num_data) << endl;

        fout_s_la_z << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                    << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << abs(s_fourier_la_z[m]) * abs(s_fourier_la_z[m]) / (4.0 * pi * num_data) << endl;

        fout_s_lo_x << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                    << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << abs(s_fourier_lo_x[m]) * abs(s_fourier_lo_x[m]) / (4.0 * pi * num_data) << endl;

        fout_s_lo_y << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                    << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << abs(s_fourier_lo_y[m]) * abs(s_fourier_lo_y[m]) / (4.0 * pi * num_data) << endl;

        fout_s_lo_z << sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]) << " "
                    << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << abs(s_fourier_lo_z[m]) * abs(s_fourier_lo_z[m]) / (4.0 * pi * num_data) << endl;

        sum_la_x[m] += abs(s_fourier_la_x[m]) * abs(s_fourier_la_x[m]) / (4.0 * pi * num_data);
        sum_la_y[m] += abs(s_fourier_la_y[m]) * abs(s_fourier_la_y[m]) / (4.0 * pi * num_data);
        sum_la_z[m] += abs(s_fourier_la_z[m]) * abs(s_fourier_la_z[m]) / (4.0 * pi * num_data);
        sum_lo_x[m] += abs(s_fourier_lo_x[m]) * abs(s_fourier_lo_x[m]) / (4.0 * pi * num_data);
        sum_lo_y[m] += abs(s_fourier_lo_y[m]) * abs(s_fourier_lo_y[m]) / (4.0 * pi * num_data);
        sum_lo_z[m] += abs(s_fourier_lo_z[m]) * abs(s_fourier_lo_z[m]) / (4.0 * pi * num_data);
        sum_la_all[m] += sum_la_x[m] + sum_la_y[m] + sum_la_z[m];
        sum_lo_all[m] += sum_lo_x[m] + sum_lo_y[m] + sum_lo_z[m];
        sum_all[m] += sum_la_all[m] + sum_lo_all[m];
      }
    }
    
    fout_s_la_x.close();
    fout_s_la_y.close();
    fout_s_la_z.close();
    fout_s_lo_x.close();
    fout_s_lo_y.close();
    fout_s_lo_z.close();
  }

  if(break_omp) return 1;
  cout << "XX ..... Finished !" << endl << endl;

  // create sum file

    ostringstream osfile_sum_la_x;    // name of the new output file of acoustic mode
    ostringstream osfile_sum_la_y;    // name of the new output file of acoustic mode
    ostringstream osfile_sum_la_z;    // name of the new output file of acoustic mode
    ostringstream osfile_sum_lo_x;    // name of the new output file of optical mode
    ostringstream osfile_sum_lo_y;    // name of the new output file of optical mode
    ostringstream osfile_sum_lo_z;    // name of the new output file of optical mode
    ostringstream osfile_sum_la_all;    // name of the new output file of acoustic mode
    ostringstream osfile_sum_lo_all;    // name of the new output file of optical mode
    ostringstream osfile_sum_all;    // name of the new output file of all mode

    osfile_sum_la_x << "sum-" << atom_name << "-" << output << "-" << wave << "-AX" << ext;
    osfile_sum_la_y << "sum-" << atom_name << "-" << output << "-" << wave << "-AY" << ext;
    osfile_sum_la_z << "sum-" << atom_name << "-" << output << "-" << wave << "-AZ" << ext;
    osfile_sum_lo_x << "sum-" << atom_name << "-" << output << "-" << wave << "-OX" << ext;
    osfile_sum_lo_y << "sum-" << atom_name << "-" << output << "-" << wave << "-OY" << ext;
    osfile_sum_lo_z << "sum-" << atom_name << "-" << output << "-" << wave << "-OZ" << ext;
    osfile_sum_la_all << "sum-" << atom_name << "-" << output << "-" << wave << "-A-ALL" << ext;
    osfile_sum_lo_all << "sum-" << atom_name << "-" << output << "-" << wave << "-O-ALL" << ext;
    osfile_sum_all << "sum-" << atom_name << "-" << output << "-" << wave << "-ALL" << ext;

    ofstream fout_sum_la_x(osfile_sum_la_x.str().c_str());
    ofstream fout_sum_la_y(osfile_sum_la_y.str().c_str());
    ofstream fout_sum_la_z(osfile_sum_la_z.str().c_str());
    ofstream fout_sum_lo_x(osfile_sum_lo_x.str().c_str());
    ofstream fout_sum_lo_y(osfile_sum_lo_y.str().c_str());
    ofstream fout_sum_lo_z(osfile_sum_lo_z.str().c_str());
    ofstream fout_sum_la_all(osfile_sum_la_all.str().c_str());
    ofstream fout_sum_lo_all(osfile_sum_lo_all.str().c_str());
    ofstream fout_sum_all(osfile_sum_all.str().c_str());

    if(!fout_sum_la_x)
    {
      cout << "Can't open " << osfile_sum_la_x.str() << "!" << endl;
      return 1;
    }
    if(!fout_sum_la_y)
    {
      cout << "Can't open " << osfile_sum_la_y.str() << "!" << endl;
      return 1;
    }
    if(!fout_sum_la_z)
    {
      cout << "Can't open " << osfile_sum_la_z.str() << "!" << endl;
      return 1;
    }
    if(!fout_sum_lo_x)
    {
      cout << "Can't open " << osfile_sum_lo_x.str() << "!" << endl;
      return 1;
    }
    if(!fout_sum_lo_y)
    {
      cout << "Can't open " << osfile_sum_lo_y.str() << "!" << endl;
      return 1;
    }
    if(!fout_sum_lo_z)
    {
      cout << "Can't open " << osfile_sum_lo_z.str() << "!" << endl;
      return 1;
    }
    if(!fout_sum_la_all)
    {
      cout << "Can't open " << osfile_sum_la_all.str() << "!" << endl;
      return 1;
    }
    if(!fout_sum_lo_all)
    {
      cout << "Can't open " << osfile_sum_lo_all.str() << "!" << endl;
      return 1;
    }
    if(!fout_sum_all)
    {
      cout << "Can't open " << osfile_sum_all.str() << "!" << endl;
      return 1;
    }

    for(int m=0; m<num_data; ++m)
    {
      if(tera*(double)(m+1)/(dt*(double)num_data) < max_hertz)
      {
        fout_sum_la_x << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_la_x[m] << endl;
        fout_sum_la_y << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_la_y[m] << endl;
        fout_sum_la_z << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_la_z[m] << endl;
        fout_sum_lo_x << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_lo_x[m] << endl;
        fout_sum_lo_y << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_lo_y[m] << endl;
        fout_sum_lo_z << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_lo_z[m] << endl;
        fout_sum_la_all << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_la_all[m] << endl;
        fout_sum_lo_all << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_lo_all[m] << endl;
        fout_sum_all << tera * (double) (m+1) / (dt * (double) num_data) / 0.0299792458 << " "
                    << sum_all[m] << endl;
      }
    }

  // create log file

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
    
    fout_log.close();
  
  // free dynamic memory
  
  delete[] zero_atom;
  delete[] file;
  
  return 0;
}

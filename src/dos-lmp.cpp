#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <math.h>
#include <mkl_dfti.h>

#define FILE_MAX 20000

#define DT 1e-14
#define TERA 1e-12

#define MAXELEMENT 16
using namespace std;

struct AtomInfo{
  string atom_type;
  double vx, vy, vz;
};

typedef enum {Si, Ge, C, Sn, Mo, S, O, X1, X2} AtomIonType;

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

/** Input: MD Scratch File List */

int dos_lmp(int iargc, char* argv[], optpara opts, string atom_conc){
  double dt = DT;
  double tera = TERA;  
  int max_data = FILE_MAX;

  if(opts.dt) dt = opts.dt;

  string output = "dos.dat";
  string first;
 
  int a, l, m, n;
  int num_data;
  int total_atom;
  int atomnum[MAXELEMENT];             // number of atom, respectively
  double atomconc[MAXELEMENT];         // concentration of atom, respectively
  for(int i=0;i<MAXELEMENT;i++){
    atomnum[i] = 0;
    atomconc[i] = 0;
  }

  double vx0, vy0, vz0, max_hertz; 
  //douteki new
  int **array;
  array = new int*[ num_data ];

  int time_step = 0;
  int atom_type;
  string atom_name = atom_conc; 
  string temp_file, temp_zero, temp_calc;
  
  ifstream file_fin(opts.read.c_str());
  if(!file_fin){
    cout << "Can't Open " << opts.read << " File!" << endl;
    return 1;
  }
  cout << " Reading ... \"" << opts.read << "\" File." << endl;
  num_data = 0;  
  while(getline(file_fin, temp_file)){
    getline(file_fin, temp_zero);
    istringstream zero_time(temp_zero.c_str());
    zero_time >> time_step;

    if(num_data+1 > max_data){
      cout << "Dynamic Memory Space is NOT Enough for Data!" << endl;
      return -1;
    }

    getline(file_fin, temp_zero);
    getline(file_fin, temp_zero);
    if(num_data == 0) {
      istringstream zero_atom(temp_zero.c_str());
      zero_atom >> total_atom;
    }

    getline(file_fin, temp_zero);
    getline(file_fin, temp_zero);
    getline(file_fin, temp_zero);
    getline(file_fin, temp_zero);

    for(int i=0; i<total_atom+1; i++)
      getline(file_fin, temp_zero);

    if(time_step != 0)
      num_data++;
  }
  file_fin.close();

  ifstream fin(opts.read.c_str());

  for(int i=0; i<9; i++)
    getline(fin, temp_zero);
  for(int i=0; i<total_atom; i++){
    int label = 0;
    getline(fin, temp_zero);
    istringstream zero_iss(temp_zero.c_str());
    zero_iss >> label >> atom_type;
    if(atom_name == "SiGe") {
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

  if(atomnum[Mo] && atomnum[S])
  {
    if(atomnum[Si] || atomnum[Ge] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "MoS2";
    max_hertz = 20;
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
  }
  else if(atomnum[Si] && atomnum[C])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[Ge] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "SiC";
    max_hertz = 49;
    atomconc[C] = (double)atomnum[C] / (atomnum[C] + atomnum[Si]);
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
  }
  else if(atomnum[Si] && atomnum[O])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "SiO2";
    max_hertz = 20;
  }
  else if(atomnum[Si] && opts.polytype == "3H")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3H-Si";
    max_hertz = 20;
  }
  else if(atomnum[Si] && opts.polytype == "3T")
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "3T-Si";
    max_hertz = 20;
  }
  else if(atomnum[Si])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Si";
    max_hertz = 20;
  }
  else if(atomnum[Ge])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Ge";
    max_hertz = 12;
  }
  else if(atomnum[Sn])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[C] || atomnum[Si]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "Sn";
    max_hertz = 8;
  }
  else if(atomnum[C])
  {
    if(atomnum[S] || atomnum[Mo] || atomnum[Ge] || atomnum[Sn]){
      cout << "atomconposition is wrong!!" << endl;
      return 1;
    }
    atom_name = "C";
    max_hertz = 49;
  }
  
  cout << " ***************************************************** " << endl;
  cout << "  Atom:" << atom_name << "  Max Frequency:" << max_hertz << "[THz]" << endl;
  cout << "  Total Atom:" << total_atom;
  if(atomnum[Si]) cout << " | Si atom:" << atomnum[Si];
  if(atomnum[Ge]) cout << " | Ge atom:" << atomnum[Ge];
  if(atomnum[Sn]) cout << " | Sn atom:" << atomnum[Sn];
  if(atomnum[C])  cout << " | C atom:"  << atomnum[C];
  if(atomnum[Mo]) cout << " | Mo atom:" << atomnum[Mo];
  if(atomnum[S]) cout  << " | S atom:"  << atomnum[S];
  if(atomnum[O]) cout  << " | O atom:"  << atomnum[O];
  if(atomnum[X1]) cout << " | Si atom:" << atomnum[X1];
  if(atomnum[X2]) cout << " | Si atom:" << atomnum[X2];
  cout << endl;
  if(atomconc[Ge]) cout << "  Ge conc:" << atomconc[Ge] << endl;
  if(atomconc[Sn]) cout << "  Sn conc:" << atomconc[Sn] << endl;
  cout << "  File Data:" << num_data << " | "
       << "MD Time:" << (double)num_data*dt/tera << "[ps]" << endl;
  cout << " ***************************************************** " << endl << endl;
  
  complex <double> cor[num_data];
  complex <double> dos[num_data]; 

  AtomInfo **vel=new AtomInfo*[total_atom];
  for(int i = 0; i < total_atom; ++i ) {
        vel[i] = new AtomInfo[ num_data ];
   }

  cout << " Reading ... Atomic Velocity ..."<< endl;
  cout << " 0%                  100%" << endl << " ";
  int bar[10];
  for(int i=0; i<10; i++)
    bar[i] = (int)(num_data*i/10);

  for(a=0; a<num_data; a++){
    int label;
    double x0, y0, z0;

    if(a%((int)num_data/10) == 0) cout << "XX" << flush;

    for(int i=0; i<9; i++)
      getline(fin, temp_calc);
    for(int i=0; i<total_atom; i++){
      label = 0;
      getline(fin, temp_calc);
      istringstream iss(temp_calc.c_str());
      iss >> label >> atom_type;
      iss >> x0 >> y0 >> z0;
      iss >> vx0 >> vy0 >> vz0;

      if(atom_name == "SiGe") {
        if(atom_type == 1) first = "Si";
        if(atom_type == 2) first = "Ge";
      }
      else if(atom_name == "GeSn") {
        if(atom_type == 1) first = "Ge";
        if(atom_type == 2) first = "Sn";
      }
      else if(atom_name == "XSiO2") {
        if(atom_type == 1) first = "X1";
        if(atom_type == 2) first = "Si";
        if(atom_type == 3) first = "O";
      }
      else if(atom_name == "SiO2") {
        if(atom_type == 1) first = "Si";
        if(atom_type == 2) first = "O";
      }
      else if(atom_name == "Si") {
        if(atom_type == 1) first = "Si";
      }
      else if(atom_name == "Ge") {
        if(atom_type == 1) first = "Ge";
      }
      else if(atom_name == "Sn") {
        if(atom_type == 1) first = "Sn";
      }
      else if(atom_name == "C") {
        if(atom_type == 1) first = "C";
      }

      vel[label-1][a].atom_type = first;
      vel[label-1][a].vx = 1.0e3 * vx0;
      vel[label-1][a].vy = 1.0e3 * vy0;
      vel[label-1][a].vz = 1.0e3 * vz0;
    }
  }

  fin.close();

  cout << "XX ..... Finished !" << endl << endl;

  cout << " Creating ... Autocorrelation Function ..." << endl;
  cout << " 0%                  10%" << endl << " ";

  int ii, ij, ik;
  for(ij=0; ij<num_data; ij++){
    if(ij%((int)num_data/10) == 0 && ij)
      cout << "XX" << endl << (int)(ij*100/num_data) << "%                  " << (int)(ij*100/num_data)+10 << "%" << endl << " ";
    if(ij%((int)num_data/100) == 0)
      cout << "XX" << flush;
    cor[ij] = complex <double> (0.0, 0.0);
    for(ik=0; ik<total_atom; ik++){
      for(ii=0; ii<num_data-ij; ii++){
        if(vel[ik][ii].atom_type != "O"){
          cor[ij] += ( vel[ik][ii].vx*vel[ik][ii+ij].vx +
                       vel[ik][ii].vy*vel[ik][ii+ij].vy +
                       vel[ik][ii].vz*vel[ik][ii+ij].vz );
        }
      }
    }
  }

  cout << "XX ..... Finished !" << endl << endl;

  for(ij=0; ij<num_data; ij++){
    if(abs(cor[0]) > 0.0){ cor[ij] = cor[ij]/cor[0];}
    else{ cout << "Error: Correlation Function Division by Zero. " << endl; return -1;}
  }
 
  cout << " Creating ... Density of States for All Direction ....." << endl;
  DFTI_DESCRIPTOR *my_desc1_handle;
  long int status;
  status = DftiCreateDescriptor (&my_desc1_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, num_data);
  status = DftiSetValue(my_desc1_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
  status = DftiCommitDescriptor (my_desc1_handle);
  status = DftiComputeForward (my_desc1_handle, cor, dos);
  status = DftiFreeDescriptor (&my_desc1_handle);
  ofstream fout(output.c_str());
  if(!fout){ cout << "Can't Open " << output << " File!" << endl; return -1;}
  double freq;
  for(a=0; a<num_data; a++){
    freq = (double)(a+1)/((double)num_data*dt)*tera;
    if(freq < max_hertz)
      fout << freq / 0.0299792458 << " " << real(dos[a]) << " " << imag(dos[a]) << " " << abs(dos[a]) << " " << endl;
  }
  fout.close();
  cout << "..... Finished !" << endl;

  for(int i = 0; i < total_atom; ++i ) {     
    delete[] vel[i];
  }

  delete[] vel;

  return 0;

}


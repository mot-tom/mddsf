#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
using namespace std;

static option options[] =
  {
    {"help", no_argument, NULL, 'h'},
    {"wave", required_argument, NULL, 'w'},
    {"div", required_argument, NULL, 'd'},
    {"cutoff", required_argument, NULL, 'c'},
    {"dispersion", no_argument, NULL, 'v'},
    {"vdos", no_argument, NULL, 'u'},
    {"dmat", required_argument, NULL, 'D'},
    {"dmat2", required_argument, NULL, 'T'},
    {"evec", required_argument, NULL, 'E'},
    {"dinam", required_argument, NULL, 'X'},
    {"loginput", required_argument, NULL, 'l'},
    {"kdiv", required_argument, NULL, 'k'},
    {"fdiv", required_argument, NULL, 'f'},
    {"poly", required_argument, NULL, 'p'},
    {"dt", required_argument, NULL, 't'},
    {"maxhz", required_argument, NULL, 'h'},
    {"hessian", required_argument, NULL, 'H'},
    {"lammps", required_argument, NULL, 'L'},
    {"range", required_argument, NULL, 'R'},
    {"hist", no_argument, NULL, 's'},
    {"sed", no_argument, NULL, 'S'},
    {"omp", required_argument, NULL, 'O'},
    {"inv", required_argument, NULL, 'I'},
    {"xyz", no_argument, NULL, 'Z'},
    {0, 0, 0, 0}
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

int dispersion(int argc, char* argv[], optpara opts);
int dos(int argc, char* argv[], optpara opts);
int dos_lmp(int argc, char* argv[], optpara opts, string atom_name);
int dos_lmp_omp(int argc, char* argv[], optpara opts, string atom_name, int num_omp);
int dsf(int argc, char* argv[], optpara opts);
int dsf_omp(int argc, char* argv[], optpara opts, int num_omp);
int dsf_lmp(int argc, char* argv[], optpara opts, string atom_name);
int dsf_lmp_omp(int argc, char* argv[], optpara opts, string atom_name, int num_omp);
int sed_lmp_omp(int argc, char* argv[], optpara opts, string atom_name, int num_omp);
int inv_lmp_omp(int argc, char* argv[], optpara opts, string atom_name, int num_omp);
int dmatrix(int argc, char* argv[], optpara opts);
int dmatrix2(int argc, char* argv[], optpara opts);
int hessian_lmp(int argc, char* argv[], optpara opts);
int hist(int argc, char* argv[], optpara opts);
int xyz_vector(int argc, char* argv[], optpara opts);

int main(int argc, char **argv)
{
       cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
       cout << "            Phonon DFT analizer for Si,Ge,Mo,S version 2.01" << endl << endl;
       cout << "                        Motohiro Tomita 2015" << endl;
       cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  int opt, index;
  int disp = 0;
  int vdos = 0;
  int dmat = 0;
  int lmp = 0;
  int hessian = 0;
  int his = 0;
  int omp = 0;
  int sed = 0;
  int inv = 0;
  int xyz = 0;
  string atom_name = "Si";
  optpara opts;
  opts.read = "";
  opts.wave = "100";
  opts.part = 30;
  opts.cutoff = 1;
  opts.dsflog = "";
  opts.k_part = 29;
  opts.f_part = 100;
  opts.dt = 1e-14;
  opts.maxhz = 0;
  opts.range = "";
  string inv_range = "";

  while((opt = getopt_long(argc, argv, "hw:d:c:vuD:T:E:X:l:k:f:p:t:z:H:L:R:sSO:I:Z", options, &index)) != -1){
    switch(opt){
      case 'h':
        cout << "Help of options." << endl;
        cout << "           -w: Input wave direction. \"100\"" << endl;
        cout << "           -d: Input division number(Int). def. 30" << endl;
        cout << "           -c: Input cutoff number(Int)." << endl;
        cout << "           -p: Input crystal poly-type. ex. 2H" << endl;
        cout << "           -t: Input timestep in MD calculation. def. 1e-14" << endl;
        cout << "           -z: Input max frequency." << endl;
        cout << "      --range: Input structure range. ex. xlo xhi ylo yhi zlo zhi" << endl;
        cout << "       --vdos: Make VDOS program." << endl;
        cout << " --dispersion: Make dispersion curve program." << endl;
        cout << "Help of --dispersion options." << endl;
        cout << "           -l: Input log file of dsf." << endl;
        cout << "           -k: Input division number of k(Int). def. 29" << endl;
        cout << "           -f: Input division number of f(Int). def. 100" << endl;
        return 0;
      case 'w':
        opts.wave = optarg;
        break;
      case 'd':
        opts.part = atoi(optarg);
        break;
      case 'c':
        opts.cutoff = atoi(optarg);
        break;
      case 'v':
        disp = 1;
        break;
      case 'u':
        vdos = 1;
        break;
      case 'T':
        opts.dsflog = optarg;
        dmat = 2;
        break;
      case 'D':
        opts.dsflog = optarg;
        dmat = 1;
        break;
      case 'E':
        opts.dsflog = optarg;
        opts.range = "evec";
        dmat = 1;
        break;
      case 'X':
        opts.dsflog = optarg;
        opts.range = "dinam";
        dmat = 1;
        break;
      case 'l':
        opts.dsflog = optarg;
        break;
      case 'k':
        opts.k_part = atoi(optarg);
        break;
      case 'f':
        opts.f_part = atoi(optarg);
        break;
      case 'p':
        opts.polytype = optarg;
        break;
      case 't':
        opts.dt = atof(optarg);
        break;
      case 'z':
        opts.maxhz = atof(optarg);
        break;
      case 'H':
        opts.dsflog = optarg;
        hessian = 1;
        break;
      case 'L':
        atom_name = optarg;
        lmp = 1;
        break;
      case 'R':
        opts.range = optarg;
        break;
      case 's':
        his = 1;
        break;
      case 'S':
        sed = 1;
        break;
      case 'O':
        omp = atoi(optarg);
        break;
      case 'I':
        inv = 1;
        inv_range = optarg;
        break;
      case 'Z':
        xyz = 1;
        break;
      default:
        cout << "Wrong option" << endl;
        return 1;
        break;
     }
  }

  if(argv[optind]) opts.read = argv[optind];

  if(opts.read == "") {
    cout << "no input files" << endl;
    return 1;
  }

  if(disp) dispersion(argc, argv, opts);
  else if(dmat == 2) dmatrix2(argc, argv, opts);
  else if(dmat) dmatrix(argc, argv, opts);
  else if(vdos && lmp && omp) dos_lmp_omp(argc, argv, opts, atom_name, omp);
  else if(vdos && lmp) dos_lmp(argc, argv, opts, atom_name);
  else if(vdos) dos(argc, argv, opts);
  else if(his) hist(argc, argv, opts);
  else if(hessian) hessian_lmp(argc, argv, opts);
  else if(xyz) xyz_vector(argc, argv, opts);
  else if(inv && lmp && omp) {opts.range = inv_range; inv_lmp_omp(argc, argv, opts, atom_name, omp);}
  else if(sed && lmp && omp) sed_lmp_omp(argc, argv, opts, atom_name, omp);
  else if(sed && lmp) sed_lmp_omp(argc, argv, opts, atom_name, omp);
  else if(lmp && omp) dsf_lmp_omp(argc, argv, opts, atom_name, 1);
  else if(lmp) dsf_lmp(argc, argv, opts, atom_name);
  else if(omp) dsf_omp(argc, argv, opts, omp);
  else dsf(argc, argv, opts);
  return 0;
}

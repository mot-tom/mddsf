#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <complex>
#include <math.h>
#include <mkl.h>
using namespace std;

const int X = 0;
const int Y = 1;
const int Z = 2;

struct coordinate{
	double x, y, z;
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

int dmatrix2(int argc, char **argv, optpara opts){
	
	double mass = 1.66054021*pow(10.0, -27.0);

       string atom_name;
	int x_dir=4, y_dir=4, z_dir=4;
       double conc=0.5;
	if(opts.dsflog.length() > 0){
             istringstream size_iss(opts.dsflog.c_str());
             size_iss >> atom_name >> conc >> x_dir >> y_dir >> z_dir;
             if(x_dir == 0 || y_dir == 0 || z_dir == 0){
		  cerr << "Structure size is wrong" << endl;
               return 1;
               }
	} else {
		  cerr << "Input option is wrong" << endl;
               return 1;
	}
	int cellnum = x_dir*y_dir*z_dir;

	int atomnum = cellnum*8;
	int inum = 3*atomnum;
	int jnum = 3*atomnum;
	int N = inum*jnum;

	double* hessian = new double[N];
	double* mi = new double[N];
	double* mj = new double[N];
	coordinate* iatom = new coordinate[N];
	coordinate* jatom = new coordinate[N];

	int UCatomnum = 2;
	const int UCinum = 3*UCatomnum;
	const int UCjnum = 3*UCatomnum;
	const int UCN = UCinum*UCjnum;

	int UCatomnum2 = UCatomnum*2;
	const int UCinum2 = 3*UCatomnum2;
	const int UCjnum2 = 3*UCatomnum2;
	const int UCN2 = UCinum2*UCjnum2;

	int UCatomnum4 = UCatomnum*4;
	const int UCinum4 = 3*UCatomnum4;
	const int UCjnum4 = 3*UCatomnum4;
	const int UCN4 = UCinum4*UCjnum4;

        int lwork = 504;

	MKL_Complex16* dinam = new MKL_Complex16[UCN];
	MKL_Complex16* dinam2 = new MKL_Complex16[UCN2];
	MKL_Complex16* dinam4 = new MKL_Complex16[UCN4];

	double* eigenval = new double[3*UCatomnum];
	MKL_Complex16* work = new MKL_Complex16[2*UCN];
	double* rwork = new double[3*UCN];

	double* eigenval2 = new double[3*UCatomnum2];
	MKL_Complex16* work2 = new MKL_Complex16[2*UCN2];
	double* rwork2 = new double[3*UCN2];

	double* eigenval4 = new double[3*UCatomnum4];
	MKL_Complex16* work4 = new MKL_Complex16[2*UCN4];
	double* rwork4 = new double[3*UCN4];
 
	double* axis = new double[3];
	axis[X] = 1.0; axis[Y] = 1.0; axis[Z] = 1.0;

	double norm = axis[X]*axis[X] + axis[Y]*axis[Y] + axis[Z]*axis[Z];
	axis[X] /= sqrt(norm); axis[Y] /= sqrt(norm); axis[Z] /= sqrt(norm);

       int division = 10;
       double lattice[3] = {0.0,0.0,0.0};
       double atom_mass = 0, atom_mass2 = 0;
       if(atom_name == "X"){
         while(!lattice[0]){
           cerr << "Please enter a lattice constant of a [A] (ex. Si : 5.431 A)" << endl;
           cin >> lattice[0];
           if(!lattice[0]) { cerr << "input is wrong!!" << endl; cin.clear(); cin.ignore(); }
          }
         while(!lattice[1]){
           cerr << "Please enter a lattice constant of b [A] (ex. Si : 5.431 A)" << endl;
           cin >> lattice[1];
           if(!lattice[1]) { cerr << "input is wrong!!" << endl; cin.clear(); cin.ignore(); }
          }
         while(!lattice[2]){
           cerr << "Please enter a lattice constant of c [A] (ex. Si : 5.431 A)" << endl;
           cin >> lattice[2];
           if(!lattice[2]) { cerr << "input is wrong!!" << endl; cin.clear(); cin.ignore(); }
          }
         while(!atom_mass){
           cerr << "Please enter an atomic mass [amu] (ex. Si : 28.0855 amu)" << endl;
           cin >> atom_mass;
           if(!atom_mass) { cerr << "input is wrong!!" << endl; cin.clear(); cin.ignore(); }
          }
         atom_mass2 = atom_mass;
       }else if(atom_name == "SiGe"){
         lattice[0] = 5.431 + 0.2*conc + 0.027*conc*conc;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
         atom_mass = 28.0855;
         atom_mass2 = 72.63;
       }else if(atom_name == "Si"){
         lattice[0] = 5.431;
         lattice[1] = lattice[0];
         lattice[2] = lattice[0];
         atom_mass = 28.0855;
         atom_mass2 = 28.0855;
       }else{
         cerr << "Atom name is wrong" << endl;
         return 1;
        }
	double space[3];
       for(int p=0;p<3;p++) space[p] = sqrt(norm)*2*M_PI/lattice[p]/division;

	double* ki = new double[3];
	ki[X] = 0.0; ki[Y] = 0.0; ki[Z] = 0.0;

        double* dk = new double[3];
        dk[X] = axis[X]*space[0]; dk[Y] = axis[Y]*space[1]; dk[Z] = axis[Z]*space[2];

	double* k = new double[3];

  cout << "---------------------------------------------" << endl;
  cout << "Atom:" << atom_name << " Mass1:" << atom_mass << " Mass2:" << atom_mass2 << endl;
  cout << "Lattice:" << lattice[0] << " x " << lattice[1] << " x " << lattice[2] << "[A]" << endl;
  cout << "---------------------------------------------" << endl;

	ifstream fin1(opts.read.c_str());
	vector<string> data;
	string hess;

	while(getline(fin1,hess)){
		data.push_back(hess);
	}
	int num1, num2;

	for(vector<string>::size_type n1=0; n1<data.size(); n1++){
		if(data[n1].find("HESSIAN MATRIX") != string::npos){num1 = n1;}
	}

	for(vector<string>::size_type n2=num1+3; n2<data.size(); n2++){
         	if(data[n2].find("########################") != string::npos){num2 = n2;}
        }
	fin1.close();

	num2 = num1+N+3;

	ofstream fout("dyn_mtx_hessian.txt");
	for(vector<string>::size_type n=num1+3; n<num2; n++){
		fout << data[n] << endl;
        }
	fout.close();

	ifstream fin2("dyn_mtx_hessian.txt",ios::in);
	if(!fin2){ return 1;}
	string temp_string;
	int l = 0;

	while(getline(fin2,temp_string)){
		istringstream iss(temp_string.c_str());
		int line, row;
		iss >> line >> row;
		if(line<inum && row<jnum){
			iss >> mi[l] >> mj[l];
			iss >> iatom[l].x >> iatom[l].y >> iatom[l].z;
	                iss >> jatom[l].x >> jatom[l].y >> jatom[l].z;
			iss >> hessian[l];
			l++;
		}
		if(l>N){
			break;
			cout << "atomnumber is wrong!!" << endl;
		}
	}
	fin2.close();

	for(int SI=0; SI<N; SI++){
		hessian[SI]*=16.0;
        }

       double smi = 1.0/sqrt(atom_mass*mass);
       double smj = 1.0/sqrt(atom_mass2*mass);

	const char jobz='V';
	const char uplo='L';
 	int info;

	complex<double> re(1.0, 0);
	complex<double> im(0, 1.0);

	complex<double>* expart = new complex<double>[N];
	complex<double> arg = 0.0;

       //i-i 2 atom mode
       cout << "1-1 2atom mode";
	for(int form=0; form<UCN; form++){
		dinam[form].real = 0;
		dinam[form].imag = 0;
	}
	for(int num=0; num<=division; num++){//100 direction
	       k[X] = ki[X]+num*dk[X];
              k[Y] = ki[Y]+num*dk[Y]*0.0;
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smi*exp(arg);
		}

		for(int i=0; i<3*UCatomnum; i++){
		  for(int j=0; j<3*UCatomnum; j++){
		    for(int k=0; k<8/UCatomnum*cellnum; k++){
			dinam[3*UCatomnum*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
			dinam[3*UCatomnum*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum, dinam, &UCinum, eigenval, work, &lwork, rwork, &info);

		for(int SI=0; SI<3*UCatomnum; SI++){
			eigenval[SI] = eigenval[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval[SI] = sqrt(eigenval[SI]);
			eigenval[SI] = eigenval[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum; num1++){
			cout << setprecision(6) << eigenval[num1] << " ";
                }
	}
	for(int form=0; form<UCN; form++){
		dinam[form].real = 0;
		dinam[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//110 direction
	       k[X] = ki[X]+(division-num)*dk[X];
              k[Y] = ki[Y]+(division-num)*dk[Y];
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smi*exp(arg);
		}

		for(int i=0; i<3*UCatomnum; i++){
		  for(int j=0; j<3*UCatomnum; j++){
		    for(int k=0; k<8/UCatomnum*cellnum; k++){
			dinam[3*UCatomnum*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
			dinam[3*UCatomnum*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum, dinam, &UCinum, eigenval, work, &lwork, rwork, &info);

		for(int SI=0; SI<3*UCatomnum; SI++){
			eigenval[SI] = eigenval[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval[SI] = sqrt(eigenval[SI]);
			eigenval[SI] = eigenval[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum; num1++){
			cout << setprecision(6) << eigenval[num1] << " ";
                }
	}
	for(int form=0; form<UCN; form++){
		dinam[form].real = 0;
		dinam[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//111 direction
	       k[X] = ki[X]+num*dk[X]*0.5;
              k[Y] = ki[Y]+num*dk[Y]*0.5;
              k[Z] = ki[Z]+num*dk[Z]*0.5;

		cout << setprecision(2) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smi*exp(arg);
		}

		for(int i=0; i<3*UCatomnum; i++){
		  for(int j=0; j<3*UCatomnum; j++){
		    for(int k=0; k<8/UCatomnum*cellnum; k++){
			dinam[3*UCatomnum*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
			dinam[3*UCatomnum*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum, dinam, &UCinum, eigenval, work, &lwork, rwork, &info);

		for(int SI=0; SI<3*UCatomnum; SI++){
			eigenval[SI] = eigenval[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval[SI] = sqrt(eigenval[SI]);
			eigenval[SI] = eigenval[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum; num1++){
			cout << setprecision(6) << eigenval[num1] << " ";
                }
	}
      cout << endl << endl;

       //i-j 2 atom mode
       cout << "1-2 2atom mode";
	for(int form=0; form<UCN; form++){
		dinam[form].real = 0;
		dinam[form].imag = 0;
	}
	for(int num=0; num<=division; num++){//100 direction
	       k[X] = ki[X]+num*dk[X];
              k[Y] = ki[Y]+num*dk[Y]*0.0;
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smj*exp(arg);
		}

		for(int i=0; i<3*UCatomnum; i++){
		  for(int j=0; j<3*UCatomnum; j++){
		    for(int k=0; k<8/UCatomnum*cellnum; k++){
			dinam[3*UCatomnum*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
			dinam[3*UCatomnum*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum, dinam, &UCinum, eigenval, work, &lwork, rwork, &info);

		for(int SI=0; SI<3*UCatomnum; SI++){
			eigenval[SI] = eigenval[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval[SI] = sqrt(eigenval[SI]);
			eigenval[SI] = eigenval[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum; num1++){
			cout << setprecision(6) << eigenval[num1] << " ";
                }
	}
	for(int form=0; form<UCN; form++){
		dinam[form].real = 0;
		dinam[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//110 direction
	       k[X] = ki[X]+(division-num)*dk[X];
              k[Y] = ki[Y]+(division-num)*dk[Y];
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smj*exp(arg);
		}

		for(int i=0; i<3*UCatomnum; i++){
		  for(int j=0; j<3*UCatomnum; j++){
		    for(int k=0; k<8/UCatomnum*cellnum; k++){
			dinam[3*UCatomnum*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
			dinam[3*UCatomnum*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum, dinam, &UCinum, eigenval, work, &lwork, rwork, &info);

		for(int SI=0; SI<3*UCatomnum; SI++){
			eigenval[SI] = eigenval[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval[SI] = sqrt(eigenval[SI]);
			eigenval[SI] = eigenval[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum; num1++){
			cout << setprecision(6) << eigenval[num1] << " ";
                }
	}
	for(int form=0; form<UCN; form++){
		dinam[form].real = 0;
		dinam[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//111 direction
	       k[X] = ki[X]+num*dk[X]*0.5;
              k[Y] = ki[Y]+num*dk[Y]*0.5;
              k[Z] = ki[Z]+num*dk[Z]*0.5;

		cout << setprecision(2) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smj*exp(arg);
		}

		for(int i=0; i<3*UCatomnum; i++){
		  for(int j=0; j<3*UCatomnum; j++){
		    for(int k=0; k<8/UCatomnum*cellnum; k++){
			dinam[3*UCatomnum*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
			dinam[3*UCatomnum*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum, dinam, &UCinum, eigenval, work, &lwork, rwork, &info);

		for(int SI=0; SI<3*UCatomnum; SI++){
			eigenval[SI] = eigenval[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval[SI] = sqrt(eigenval[SI]);
			eigenval[SI] = eigenval[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum; num1++){
			cout << setprecision(6) << eigenval[num1] << " ";
                }
	}
      cout << endl << endl;

       //j-j 2 atom mode
       cout << "2-2 2atom mode";
	for(int form=0; form<UCN; form++){
		dinam[form].real = 0;
		dinam[form].imag = 0;
	}
	for(int num=0; num<=division; num++){//100 direction
	       k[X] = ki[X]+num*dk[X];
              k[Y] = ki[Y]+num*dk[Y]*0.0;
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smj*smj*exp(arg);
		}

		for(int i=0; i<3*UCatomnum; i++){
		  for(int j=0; j<3*UCatomnum; j++){
		    for(int k=0; k<8/UCatomnum*cellnum; k++){
			dinam[3*UCatomnum*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
			dinam[3*UCatomnum*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum, dinam, &UCinum, eigenval, work, &lwork, rwork, &info);

		for(int SI=0; SI<3*UCatomnum; SI++){
			eigenval[SI] = eigenval[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval[SI] = sqrt(eigenval[SI]);
			eigenval[SI] = eigenval[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum; num1++){
			cout << setprecision(6) << eigenval[num1] << " ";
                }
	}
	for(int form=0; form<UCN; form++){
		dinam[form].real = 0;
		dinam[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//110 direction
	       k[X] = ki[X]+(division-num)*dk[X];
              k[Y] = ki[Y]+(division-num)*dk[Y];
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smj*smj*exp(arg);
		}

		for(int i=0; i<3*UCatomnum; i++){
		  for(int j=0; j<3*UCatomnum; j++){
		    for(int k=0; k<8/UCatomnum*cellnum; k++){
			dinam[3*UCatomnum*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
			dinam[3*UCatomnum*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum, dinam, &UCinum, eigenval, work, &lwork, rwork, &info);

		for(int SI=0; SI<3*UCatomnum; SI++){
			eigenval[SI] = eigenval[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval[SI] = sqrt(eigenval[SI]);
			eigenval[SI] = eigenval[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum; num1++){
			cout << setprecision(6) << eigenval[num1] << " ";
                }
	}
	for(int form=0; form<UCN; form++){
		dinam[form].real = 0;
		dinam[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//111 direction
	       k[X] = ki[X]+num*dk[X]*0.5;
              k[Y] = ki[Y]+num*dk[Y]*0.5;
              k[Z] = ki[Z]+num*dk[Z]*0.5;

		cout << setprecision(2) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smj*smj*exp(arg);
		}

		for(int i=0; i<3*UCatomnum; i++){
		  for(int j=0; j<3*UCatomnum; j++){
		    for(int k=0; k<8/UCatomnum*cellnum; k++){
			dinam[3*UCatomnum*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
			dinam[3*UCatomnum*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum*k]*hessian[3*atomnum*i+j+3*UCatomnum*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum, dinam, &UCinum, eigenval, work, &lwork, rwork, &info);

		for(int SI=0; SI<3*UCatomnum; SI++){
			eigenval[SI] = eigenval[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval[SI] = sqrt(eigenval[SI]);
			eigenval[SI] = eigenval[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum; num1++){
			cout << setprecision(6) << eigenval[num1] << " ";
                }
	}
      cout << endl << endl;

       //i-i-i-i 4 atom mode
       cout << "1-1-1-1 4atom mode";
	for(int form=0; form<UCN2; form++){
		dinam2[form].real = 0;
		dinam2[form].imag = 0;
	}
	for(int num=0; num<=division; num++){//100 direction
	       k[X] = ki[X]+num*dk[X];
              k[Y] = ki[Y]+num*dk[Y]*0.0;
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smi*exp(arg);
		}

		for(int i=0; i<3*UCatomnum2; i++){
		  for(int j=0; j<3*UCatomnum2; j++){
		    for(int k=0; k<8/UCatomnum2*cellnum; k++){
			dinam2[3*UCatomnum2*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum2*k]*hessian[3*atomnum*i+j+3*UCatomnum2*k]);
			dinam2[3*UCatomnum2*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum2*k]*hessian[3*atomnum*i+j+3*UCatomnum2*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum2, dinam2, &UCinum2, eigenval2, work2, &lwork, rwork2, &info);

		for(int SI=0; SI<3*UCatomnum2; SI++){
			eigenval2[SI] = eigenval2[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval2[SI] = sqrt(eigenval2[SI]);
			eigenval2[SI] = eigenval2[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum2; num1++){
			cout << setprecision(6) << eigenval2[num1] << " ";
                }
	}
	for(int form=0; form<UCN2; form++){
		dinam2[form].real = 0;
		dinam2[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//110 direction
	       k[X] = ki[X]+(division-num)*dk[X];
              k[Y] = ki[Y]+(division-num)*dk[Y];
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smi*exp(arg);
		}

		for(int i=0; i<3*UCatomnum2; i++){
		  for(int j=0; j<3*UCatomnum2; j++){
		    for(int k=0; k<8/UCatomnum2*cellnum; k++){
			dinam2[3*UCatomnum2*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum2*k]*hessian[3*atomnum*i+j+3*UCatomnum2*k]);
			dinam2[3*UCatomnum2*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum2*k]*hessian[3*atomnum*i+j+3*UCatomnum2*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum2, dinam2, &UCinum2, eigenval2, work2, &lwork, rwork2, &info);

		for(int SI=0; SI<3*UCatomnum2; SI++){
			eigenval2[SI] = eigenval2[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval2[SI] = sqrt(eigenval2[SI]);
			eigenval2[SI] = eigenval2[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum2; num1++){
			cout << setprecision(6) << eigenval2[num1] << " ";
                }
	}
	for(int form=0; form<UCN2; form++){
		dinam2[form].real = 0;
		dinam2[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//111 direction
	       k[X] = ki[X]+num*dk[X]*0.5;
              k[Y] = ki[Y]+num*dk[Y]*0.5;
              k[Z] = ki[Z]+num*dk[Z]*0.5;

		cout << setprecision(2) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smi*exp(arg);
		}

		for(int i=0; i<3*UCatomnum2; i++){
		  for(int j=0; j<3*UCatomnum2; j++){
		    for(int k=0; k<8/UCatomnum2*cellnum; k++){
			dinam2[3*UCatomnum2*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum2*k]*hessian[3*atomnum*i+j+3*UCatomnum2*k]);
			dinam2[3*UCatomnum2*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum2*k]*hessian[3*atomnum*i+j+3*UCatomnum2*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum2, dinam2, &UCinum2, eigenval2, work2, &lwork, rwork2, &info);

		for(int SI=0; SI<3*UCatomnum2; SI++){
			eigenval2[SI] = eigenval2[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval2[SI] = sqrt(eigenval2[SI]);
			eigenval2[SI] = eigenval2[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum2; num1++){
			cout << setprecision(6) << eigenval2[num1] << " ";
                }
	}
      cout << endl << endl;

       //i-i-i-i-i-i-i-i 8 atom mode
       cout << "1-1-1-1-1-1-1-1 8atom mode";
	for(int form=0; form<UCN4; form++){
		dinam4[form].real = 0;
		dinam4[form].imag = 0;
	}
	for(int num=0; num<=division; num++){//100 direction
	       k[X] = ki[X]+num*dk[X];
              k[Y] = ki[Y]+num*dk[Y]*0.0;
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smi*exp(arg);
		}

		for(int i=0; i<3*UCatomnum4; i++){
		  for(int j=0; j<3*UCatomnum4; j++){
		    for(int k=0; k<8/UCatomnum4*cellnum; k++){
			dinam4[3*UCatomnum4*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum4*k]*hessian[3*atomnum*i+j+3*UCatomnum4*k]);
			dinam4[3*UCatomnum4*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum4*k]*hessian[3*atomnum*i+j+3*UCatomnum4*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum4, dinam4, &UCinum4, eigenval4, work4, &lwork, rwork4, &info);

		for(int SI=0; SI<3*UCatomnum4; SI++){
			eigenval4[SI] = eigenval4[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval4[SI] = sqrt(eigenval4[SI]);
			eigenval4[SI] = eigenval4[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum4; num1++){
			cout << setprecision(6) << eigenval4[num1] << " ";
                }
	}
	for(int form=0; form<UCN4; form++){
		dinam4[form].real = 0;
		dinam4[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//110 direction
	       k[X] = ki[X]+(division-num)*dk[X];
              k[Y] = ki[Y]+(division-num)*dk[Y];
              k[Z] = ki[Z]+num*dk[Z]*0.0;

		cout << setprecision(1) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smi*exp(arg);
		}

		for(int i=0; i<3*UCatomnum4; i++){
		  for(int j=0; j<3*UCatomnum4; j++){
		    for(int k=0; k<8/UCatomnum4*cellnum; k++){
			dinam4[3*UCatomnum4*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum4*k]*hessian[3*atomnum*i+j+3*UCatomnum4*k]);
			dinam4[3*UCatomnum4*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum4*k]*hessian[3*atomnum*i+j+3*UCatomnum4*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum4, dinam4, &UCinum4, eigenval4, work4, &lwork, rwork4, &info);

		for(int SI=0; SI<3*UCatomnum4; SI++){
			eigenval4[SI] = eigenval4[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval4[SI] = sqrt(eigenval4[SI]);
			eigenval4[SI] = eigenval4[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum4; num1++){
			cout << setprecision(6) << eigenval4[num1] << " ";
                }
	}
	for(int form=0; form<UCN4; form++){
		dinam4[form].real = 0;
		dinam4[form].imag = 0;
	}
	for(int num=1; num<=division; num++){//111 direction
	       k[X] = ki[X]+num*dk[X]*0.5;
              k[Y] = ki[Y]+num*dk[Y]*0.5;
              k[Z] = ki[Z]+num*dk[Z]*0.5;

		cout << setprecision(2) << endl;
		cout << " " << k[X]/(2*M_PI/lattice[0]) << " ";
		cout << " " << k[Y]/(2*M_PI/lattice[1]) << " ";
		cout << " " << k[Z]/(2*M_PI/lattice[2]) << " ";

		for(int dline=0; dline<N; dline++){
			arg = -im*(k[X]*(iatom[dline].x-jatom[dline].x)+k[Y]*(iatom[dline].y-jatom[dline].y)+k[Z]*(iatom[dline].z-jatom[dline].z));
			expart[dline] = smi*smi*exp(arg);
		}

		for(int i=0; i<3*UCatomnum4; i++){
		  for(int j=0; j<3*UCatomnum4; j++){
		    for(int k=0; k<8/UCatomnum4*cellnum; k++){
			dinam4[3*UCatomnum4*i+j].real += real(expart[3*atomnum*i+j+3*UCatomnum4*k]*hessian[3*atomnum*i+j+3*UCatomnum4*k]);
			dinam4[3*UCatomnum4*i+j].imag += -imag(expart[3*atomnum*i+j+3*UCatomnum4*k]*hessian[3*atomnum*i+j+3*UCatomnum4*k]);
		    }
		  }
		}

		zheev(&jobz, &uplo, &UCinum4, dinam4, &UCinum4, eigenval4, work4, &lwork, rwork4, &info);

		for(int SI=0; SI<3*UCatomnum4; SI++){
			eigenval4[SI] = eigenval4[SI]*pow(10.0,-24.0)/(2*2)/(M_PI*M_PI);
			eigenval4[SI] = sqrt(eigenval4[SI]);
			eigenval4[SI] = eigenval4[SI]/0.0299792458;
		}

		for(int num1=0; num1<3*UCatomnum4; num1++){
			cout << setprecision(6) << eigenval4[num1] << " ";
                }
	}
      cout << endl << endl;

	delete[] hessian;
	delete[] mi;
	delete[] mj;
        delete[] iatom;
        delete[] jatom;
        delete[] axis;
        delete[] ki;
        delete[] dk;
        delete[] k;
        delete[] expart;
        delete[] dinam;
        delete[] dinam2;
        delete[] dinam4;
        delete[] work;
        delete[] rwork;
        delete[] eigenval;
        delete[] eigenval2;
        delete[] eigenval4;

	return 0;
}


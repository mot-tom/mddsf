#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
using namespace std;

// declare struct info

struct info
{
  double k;   // wave vector
  double f;   // frequency
  double d;   // powerspectrum of dynamical structure factor
  double g;   // gauss powerspectrum of dynamical structure factor
  int a;      // number of atom
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

int hist(int argc, char* argv[], optpara opts)
{
  /*********************************************/
  
  const double color = pow(10.0, 0.0);     // value for color compensation
  const double pi = 3.14159265;             // pi
  const double sigma = 2.0*0.0299792458;   // sigma of gauss
  double gau_den = 1.0/sqrt(2*pi*sigma*sigma);
  const double x_st = 0.0;                  // limit of start edge in x range
  const double y_st = 0.2;                  // limit of start edge in y range
  
  int k_part = 0;                    // value to divide wave vector
  int f_part = 100;                   // value to divide frequency  
  const int max_file = 100;           // maximum number of file
  int wave = 1;                       // direction of wave vector; "1":[100]; "2":[110]; "3":[111]

  string output = "hist";                    // base name of the output file
  string ext = ".dat";                      // name of extension for the output file

  for(int l=0; l<=2; ++l)
  {
    string swave, pick;        // store part of string wave
    int dir;
    swave = opts.wave;
    pick = swave.substr(l,1);
    
    istringstream is_pick(pick.c_str());
    is_pick >> dir;
    
    if(dir != 0.0 && dir != 1.0)
    {
      cout << "Setting of wave vector direction is wrong!" << endl;
      return 1;
    }
  }
  if(opts.wave == "100" || opts.wave == "010"|| opts.wave == "001") wave = 1;
  if(opts.wave == "110" || opts.wave == "101"|| opts.wave == "011") wave = 2;
  if(opts.wave == "111") wave = 3;
  
  /*********************************************/
  
  double lattice[3];                        // lattice constant[A]
  double dk;                                // mesh width of wave vector
  double df;                                // mesh width of frequency
  
  double k_min = 0.0;                       // minimum of wave vector
  double k_max = 1.0;                       // maximum of wave vector
  
  double f_min = 0.0;                       // minimum of frequency
  double f_max = 20.0;                       // maximum of frequency
  
  int sec;                                  // section number
  
  string atom_name;                         // atom name(Si or Ge)
  string temp_ar, temp_file, temp_data;     // line of string (list/file/data)
   
  if(wave != 1 && wave != 2 && wave != 3)
  {
    cout << "Setting of direction of wave vector is wrong!" << endl;
    return 1;
  }

  if(opts.dsflog.length() > 0)
  {
    ifstream fin_log(opts.dsflog.c_str());
  
    if(!fin_log)
    {
      cout << "Can't open " << opts.dsflog << "!" << endl;
      return 1;
    }
  
    cout << "Read \"" << opts.dsflog << "\"." << endl;

    while(getline(fin_log, temp_ar))
    {
      istringstream is_file(temp_ar.c_str());
      is_file >> temp_file;

      if(temp_file.find("ATMNAME") != string::npos)
      {
        getline(fin_log, temp_ar);
        istringstream is_read(temp_ar.c_str());
        is_read >> atom_name;
      }
      else if(temp_file.find("LATTICE_X") != string::npos)
      {
        getline(fin_log, temp_ar);
        istringstream is_read(temp_ar.c_str());
        is_read >> lattice[0];
      }
      else if(temp_file.find("LATTICE_Y") != string::npos)
      {
        getline(fin_log, temp_ar);
        istringstream is_read(temp_ar.c_str());
        is_read >> lattice[1];
      }
      else if(temp_file.find("LATTICE_Z") != string::npos)
      {
        getline(fin_log, temp_ar);
        istringstream is_read(temp_ar.c_str());
        is_read >> lattice[2];
      }
      else if(temp_file.find("MAXFREQ") != string::npos)
      {
        getline(fin_log, temp_ar);
        istringstream is_read(temp_ar.c_str());
        is_read >> f_max;
      }
      else if(temp_file.find("MAXWAVE") != string::npos)
      {
        getline(fin_log, temp_ar);
        istringstream is_read(temp_ar.c_str());
        is_read >> k_max;
      }
      else getline(fin_log, temp_ar);
    }

    fin_log.close();
  }

  if(opts.maxhz) f_max = opts.maxhz;

  ifstream fin_ar(opts.read.c_str());
  while(getline(fin_ar, temp_ar))
  {
    if(temp_ar.find("#") == string::npos){k_part++;}
  }
  fin_ar.close();

  k_part--;

  if(opts.k_part) k_part = opts.k_part;
  if(opts.f_part) f_part = opts.f_part;
  cout << "fdiv " << f_part << endl;
  cout << "kdiv " << k_part << endl;

  // allocate dynamic memory
  info* pdr = new info[k_part*f_part];

  dk = (k_max - k_min) / (double)k_part;
  df = (f_max - f_min) / (double)f_part;
      
  // initiallization
  for(int l=0; l<k_part*f_part; ++l)
  {
    pdr[l].k = ((double)(l % k_part) + 0.5) * dk;
    pdr[l].f = ((int)(l / k_part) + 0.5) * df;
    pdr[l].d = 0.0;
    pdr[l].g = 0.0;
    pdr[l].a = 0;
  }

  fin_ar.open(opts.read.c_str());
  for(int m=0; m<k_part; ++m)
  {
    getline(fin_ar, temp_ar);
    while(temp_ar.find("#") != string::npos) getline(fin_ar, temp_ar);
    istringstream is_data(temp_ar.c_str());
    stringstream ss_data;
    double data_f;
    is_data >> data_f;
    while(is_data.str().find(ss_data.str()) != string::npos && !is_data.eof())
    {
      is_data >> data_f;
      ss_data.clear();
      ss_data.str("");
      ss_data << scientific << setprecision(6) << data_f;
      int l=1;
      data_f *= 0.0299792458;
      while(data_f > pdr[l*k_part+m].f && l < f_part) l++;
      if(l >= f_part) {continue;}
      double temp_f = pdr[l*k_part+m].f - data_f;
      if(temp_f > pdr[(l-1)*k_part+m].f - data_f)
        {for(int i=1; i<f_part; ++i){pdr[i*k_part+m].g += gau_den*exp(-pow(pdr[(i-1)*k_part+m].f-data_f,2)/(2*sigma*sigma));}
         pdr[(l-1)*k_part+m].d++; pdr[(l-1)*k_part+m].a++;}
      else {for(int i=1; i<f_part; ++i){pdr[i*k_part+m].g += gau_den*exp(-pow(pdr[i*k_part+m].f-data_f,2)/(2*sigma*sigma));} 
            pdr[l*k_part+m].d++; pdr[l*k_part+m].a++;}
    }
  }
  fin_ar.close();

  cout << "------------------------------" << endl;
  cout << "Atom:" << atom_name << "  Max Frequency:" << f_max << "[THz]" << endl;
  cout << "Lattice:" << lattice[0] << " x " << lattice[1] << " x " << lattice[2] << "[A]" << endl;
  cout << "X mesh:" << k_part << "  Y mesh:" << f_part << endl;
  cout << "------------------------------" << endl;
  
  // write the wave vector[1/A], frequency[THz] and corrected power spectrum in each section
  
  ostringstream ost;         // name of the new output file of transverse mode

  if(atom_name.c_str())
    ost << atom_name << "-";
  ost << output;
  if(opts.maxhz) ost << "-fmax" << f_max;
  ost << ext;
  
  ofstream fout(ost.str().c_str());
  
  if(!fout)
  {
    cout << "Can't open " << output << "!" << endl;
    return 1;
  }
  
  for(int l=0; l<f_part; ++l)
  {
    for(int m=0; m<k_part; ++m)
    {
      sec = k_part*l + m;
      
      if(pdr[sec].k > x_st && pdr[sec].f > y_st)
        fout << pdr[sec].k << " " << pdr[sec].f / 0.0299792458 << " " << pdr[sec].g*color << endl;
      
      else
        fout << pdr[sec].k << " " << pdr[sec].f / 0.0299792458 << " " << 0.0 << endl;
      
    }
    
    fout << endl;
  }
  
  fout.close();
  
  cout << "\"" << ost.str() << "\" is made." << " "<< pdr[sec].a << endl;
  
  // free dynamic memory
  
  delete[] pdr;
  
  return 0;
}


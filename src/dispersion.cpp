#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
using namespace std;

// Outline of this program
/*
  Calculate dynamic structure factor in each section
  to display phonon dispersion relation.
*/

double round(double round, double digit);   // function to round off the number

// declare struct info

struct info
{
  double k;   // wave vector
  double f;   // frequency
  double d;   // powerspectrum of dynamical structure factor
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

int dispersion(int iargc, char* argv[], optpara opts)
{
  /*********************************************/
  
  const double color = pow(10.0, 0.0);     // value for color compensation
  const double pi = 3.14159265;             // pi
  const double x_st = 0.0;                  // limit of start edge in x range
  const double y_st = 0.2;                  // limit of start edge in y range
  
  int k_part = 29;                    // value to divide wave vector
  int f_part = 100;                   // value to divide frequency  
  const int max_file = 100;           // maximum number of file
  int wave = 1;                       // direction of wave vector; "1":[100]; "2":[110]; "3":[111]

  string output = "pdr";                    // base name of the output file
  string ext = ".dat";                      // name of extension for the output file

  for(int l=0; l<=2; ++l)
  {
    string swave, pick;        // store part of string wave
    int dir;
    swave = opts.wave;
    pick = swave.substr(l,1);
    
    istringstream is_pick(pick.c_str());
    is_pick >> dir;
    
    if(dir != 0.0 && dir != 1.0 && dir != 2.0 && dir != 4.0)
    {
      cout << "Setting of wave vector direction is wrong!" << endl;
      return 1;
    }
  }
  if(opts.wave == "100" || opts.wave == "010" || opts.wave == "001" || opts.wave == "200" || opts.wave == "400") wave = 1;
  if(opts.wave == "110" || opts.wave == "101" || opts.wave == "011" || opts.wave == "220") wave = 2;
  if(opts.wave == "111") wave = 3;

  if(opts.k_part) k_part = opts.k_part;
  if(opts.f_part) f_part = opts.f_part;

  cout << "wave " << opts.wave << endl;
  cout << "kdiv " << k_part << endl;
  cout << "fdiv " << f_part << endl;
  
  /*********************************************/
  
  double lattice[3];                        // lattice constant[A]
  double dk;                                // mesh width of wave vector
  double df;                                // mesh width of frequency
  
  double k_min = 0.0;                       // minimum of wave vector
  double k_max = 1.0;                       // maximum of wave vector
  
  double f_min = 0.0;                       // minimum of frequency
  double f_max = 0.0;                       // maximum of frequency
  
  int num_data = 0;                         // number of data
  int num_file = 0;                         // number of file
  int sec;                                  // section number
  
  string atom_name;                         // atom name(Si or Ge)
  string temp_ar, temp_file, temp_data;     // line of string (list/file/data)
  
  // allocate dynamic memory
  
  info* pdr = new info[k_part*f_part];  
   
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
  
  if(!fin_ar)
  {
    cout << "Can't open " << opts.read << "!" << endl;
    return 1;
  }
  
  cout << "Read \"" << opts.read << "\"." << endl;
  
  // read the information of data in each file
  
  while(getline(fin_ar, temp_ar))
  {
    istringstream is_file(temp_ar.c_str());
    is_file >> temp_file;
    
    if(num_file == 0 && f_max == 0)
    {
      if(temp_file.find("Si") != string::npos)
      {
        atom_name = "Si";
        lattice[0] = 5.43;
        lattice[1] = lattice[0];
        lattice[2] = lattice[0];
        f_max = 20.0;
      }
      else if(temp_file.find("Ge") != string::npos)
      {
        atom_name = "Ge";
        lattice[0] = 5.65;
        lattice[1] = lattice[0];
        lattice[2] = lattice[0];
        f_max = 12.0;
      }

      else      
      {
        cout << "Set the atom name (1:Si, 2:Ge) : ";
        
        int input_num;       // input from the keyboard
        cin >> input_num;
        
        if(input_num == 1)
        {
          atom_name = "Si";
          lattice[0] = 5.43;
          lattice[1] = lattice[0];
          lattice[2] = lattice[0];
          f_max = 20.0;
        }
        
        else if(input_num == 2)
        {
          atom_name = "Ge";
          lattice[0] = 5.65;
          lattice[1] = lattice[0];
          lattice[2] = lattice[0];
          f_max = 12.0;
        }
        
        else
        {
          cout << "Fail to set the atom name!" << endl;
          return 1;
        }
      }
    }  
    if(num_file == 0)
    {
      dk = (k_max - k_min) / (double)k_part;
      df = (f_max - f_min) / (double)f_part;
      
      // initiallization

      for(int l=0; l<k_part*f_part; ++l)
      {
        pdr[l].k = ((double)(l % k_part) + 0.5) * dk;
        pdr[l].f = ((int)(l / k_part) + 0.5) * df;
        pdr[l].d = 0.0;
        pdr[l].a = 0;
      }
    }  
    
    ifstream fin_file(temp_file.c_str());
    
    if(!fin_file)
    {
      cout << "Can't open " << temp_file << "!" << endl;
      return 1;
    }
    
    // define the area(mesh) of the data 
     
    while(getline(fin_file, temp_data))
    {
      info data;   // data stored temporarily
      
      istringstream is_data(temp_data.c_str());
      
      if(temp_data != "")
      {
        is_data >> data.k >> data.f >> data.d;
        data.f = data.f * 0.0299792458;

        // normalized the wave vector 
        
        //if(opts.wave == "100" && atom_name == "3T-Si")
          //data.k /= 2.0*pi/sqrt(2.0)/lattice[0];

        //if(opts.wave == "100" && atom_name == "3H-Si")
          //data.k /= 2.0*pi/sqrt(2.0)/lattice[0];

        if(opts.wave == "100")
          data.k /= 2.0*pi/lattice[0];
        else if(opts.wave == "010")
          data.k /= 2.0*pi/lattice[1];
        else if(opts.wave == "001")
          data.k /= 2.0*pi/lattice[2];

        else if(opts.wave == "200")
          data.k /= 4.0*pi/lattice[0];

        else if(opts.wave == "400")
          data.k /= 8.0*pi/lattice[0];
        
        else if(opts.wave == "110")
          data.k /= 2.0*pi*sqrt(1.0/(lattice[0]*lattice[0]) + 1.0/(lattice[1]*lattice[1]));
        else if(opts.wave == "011")
          data.k /= 2.0*pi*sqrt(1.0/(lattice[1]*lattice[1]) + 1.0/(lattice[2]*lattice[2]));
        else if(opts.wave == "101")
          data.k /= 2.0*pi*sqrt(1.0/(lattice[0]*lattice[0]) + 1.0/(lattice[2]*lattice[2]));

        else if(opts.wave == "220")
          data.k /= 2.0*pi*sqrt(1.0/(lattice[0]*lattice[0]) + 1.0/(lattice[1]*lattice[1]));
        
        else
          data.k /= 2.0*pi*sqrt(1.0/(lattice[0]*lattice[0]) + 1.0/(lattice[1]*lattice[1]) + 1.0/(lattice[2]*lattice[2]));
        
        data.k = round(data.k, 5.0);
        
        if(data.k > k_max)
        {
          cout << "Setting of wave vector value is wrong!" << endl;
          cout.precision(10);
          cout << data.k << endl;
          return 1;
        }
         
        if(data.f <= f_max)
        {
          int k_num = (int)(data.k/dk);     // number of x direction of the data
          int f_num = (int)(data.f/df);     // number of y direction of the data
          
          // correct the area of the data when the wave vector is [110] and define the boundary condition
          /*
          if(wave == 2)
          {
            k_num += k_part-1;
            
            if(k_num < 0)
              k_num++;
          }
          */
          if(data.k == k_num*dk && data.k != 0.0)
            k_num--;
          
          if(data.f == f_num*df && data.f != 0.0)
            f_num -= k_part;
          
          sec = k_num + k_part*f_num;
          
          if(sec < 0 || sec >= k_part*f_part)
          {
            cout << "Assessment of the mesh is wrong!" << endl;
            cout << k_num << " " << f_num << endl;
            cout << data.k << " " << data.f << " " << sec << endl;
            return 1;
          }
          
          pdr[sec].d += data.d;
          pdr[sec].a++; 
          
          if(num_file == 0)
            num_data++; 
          
        }
      }
    }
    
    num_file++;
    
    fin_file.close();
  }

  fin_ar.close();
  
  cout << "------------------------------" << endl;
  cout << "Atom:" << atom_name << "  Max Frequency:" << f_max << "[THz]" << endl;
  cout << "Lattice:" << lattice[0] << " x " << lattice[1] << " x " << lattice[2] << "[A]" << endl;
  cout << "File:" << num_file << "  Data of one file:" << num_data << endl;
  cout << "X mesh:" << k_part << "  Y mesh:" << f_part << endl;
  cout << "------------------------------" << endl;
  
  if(num_file <= k_part)
    cout << "X mesh is not appropriate" << endl;
  
  if(num_data <= f_part)
    cout << "Y mesh is not appropriate" << endl;
  
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
      
      if(pdr[sec].a != 0 && pdr[sec].k > x_st && pdr[sec].f > y_st)
        fout << pdr[sec].k << " " << pdr[sec].f / 0.0299792458 << " " << pdr[sec].d*color/(double)pdr[sec].a << endl;
      
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

// round off the number to "digit" places of decimals

double round(double round, double digit)
{
  return (double)((int)(round * pow(10.0, digit) + 0.5)) / pow(10.0, digit);
}

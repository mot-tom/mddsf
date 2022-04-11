#include <Python.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
using namespace std;

const char *init_Py =
"#!/usr/local/env python\n"
"# -*- coding: utf-8 -*-\n"

//#Import Library

"import warnings\n"
"import sys\n"
"import matplotlib.pyplot as plt\n"
"import re\n"

"from pylab import *\n"
"from numpy import *\n"

"warnings.simplefilter(\"ignore\", DeprecationWarning)\n";

const char *density_Py =
//# Outline of this program #
//# Plot the dispersion relation.

//#Main Routine START

//#Set parameter

//####################

"MESH_MAX = 1000\n"

"fig = ''\n"           //# figure name

"x_st = 0.0\n"         //# minimum of x range
"x_ed = 1.0\n"         //# maximum of x range

"y_st = 0.0\n"         //# minimum of y range

"c_limit = 'N'\n"      //# limit of contrasting; 'Y'(Yes) or 'N'(No)
"c_st = 0.0\n"         //# minimum of contrasting
"c_ed = 1.0\n"         //# maximum of contrasting
"d_ed = 1.0\n"         //# maximum of density

//####################

"density = [[0 for j in range(MESH_MAX)] for i in range (MESH_MAX)]\n"   //# density

"argv = sys.argv\n"    //# command line arguments
"iargc = len(argv)\n"  //# number of command line arguments

"if(iargc < 2):\n"
"  print 'Input file is nothing!'\n"
"  quit()\n"

"if len(argv[4]) > 0:\n"
"  c_ed = float(argv[4])\n"
"  d_ed = c_ed*15\n"
"  c_limit = 'Y'\n"

"if len(argv[5]) > 0:\n"
"  c_st = float(argv[5])\n"

"if len(argv[6]) > 0:\n"
"  d_ed = float(argv[6])\n"

/*
"if argv[1].find('MoS2') != -1:\n"
"  print 'Atom:MoS2'\n"
"  y_ed = 20.0 / 0.0299792458\n"

"elif argv[1].find('GeSiSn') != -1:\n"
"  print 'Atom:GeSiSn'\n"
"  y_ed = 20.0 / 0.0299792458\n"

"elif argv[1].find('SiGe') != -1:\n"
"  print 'Atom:SiGe'\n"
"  y_ed = 20.0 / 0.0299792458\n"

"elif argv[1].find('SiC') != -1:\n"
"  print 'Atom:SiGe'\n"
"  y_ed = 70.0 / 0.0299792458\n"

"elif argv[1].find('GeSn') != -1:\n"
"  print 'Atom:GeSn'\n"
"  y_ed = 12.0 / 0.0299792458\n"

"elif argv[1].find('Si') != -1:\n"
"  print 'Atom:Si'\n"
"  y_ed = 20.0 / 0.0299792458\n"

"elif argv[1].find('Ge') != -1:\n"
"  print 'Atom:Ge'\n"
"  y_ed = 12.0 / 0.0299792458\n"

"elif argv[1].find('Sn') != -1:\n"
"  print 'Atom:Sn'\n"
"  y_ed = 8.0 / 0.0299792458\n"

"elif argv[1].find('C') != -1:\n"
"  print 'Atom:C'\n"
"  y_ed = 70.0 / 0.0299792458\n"
*/

"if argv[1].find('MoS2') != -1:\n"
"  print 'Atom:MoS2'\n"
"  y_ed = 20.0\n"

"elif argv[1].find('GeSiSn') != -1:\n"
"  print 'Atom:GeSiSn'\n"
"  y_ed = 20.0\n"

"elif argv[1].find('SiGe') != -1:\n"
"  print 'Atom:SiGe'\n"
"  y_ed = 20.0\n"

"elif argv[1].find('SiC') != -1:\n"
"  print 'Atom:SiGe'\n"
"  y_ed = 70.0\n"

"elif argv[1].find('GeSn') != -1:\n"
"  print 'Atom:GeSn'\n"
"  y_ed = 12.0\n"

"elif argv[1].find('Si') != -1:\n"
"  print 'Atom:Si'\n"
"  y_ed = 20.0\n"

"elif argv[1].find('Ge') != -1:\n"
"  print 'Atom:Ge'\n"
"  y_ed = 12.0\n"

"elif argv[1].find('Sn') != -1:\n"
"  print 'Atom:Sn'\n"
"  y_ed = 8.0\n"

"elif argv[1].find('C') != -1:\n"
"  print 'Atom:C'\n"
"  y_ed = 70.0\n"

"else:\n"
"  input_num = raw_input('Please enter a number(1:Si, 2:Ge, 3:Sn, 4:SiGe, 5:GeSn):')\n"
/*
"  if input_num == '1':\n"
"    print 'Atom:Si'\n"
"    y_ed = 20.0 / 0.0299792458\n"           //# maximum of y range

"  elif input_num == '2':\n"
"    print 'Atom:Ge'\n"
"    y_ed = 12.0 / 0.0299792458\n"

"  elif input_num == '3':\n"
"    print 'Atom:Sn'\n"
"    y_ed = 8.0 / 0.0299792458\n"

"  elif input_num == '4':\n"
"    print 'Atom:SiGe'\n"
"    y_ed = 20.0 / 0.0299792458\n"

"  elif input_num == '5':\n"
"    print 'Atom:GeSn'\n"
"    y_ed = 12.0 / 0.0299792458\n"
*/

"  if input_num == '1':\n"
"    print 'Atom:Si'\n"
"    y_ed = 20.0\n"           //# maximum of y range

"  elif input_num == '2':\n"
"    print 'Atom:Ge'\n"
"    y_ed = 12.0\n"

"  elif input_num == '3':\n"
"    print 'Atom:Sn'\n"
"    y_ed = 8.0\n"

"  elif input_num == '4':\n"
"    print 'Atom:SiGe'\n"
"    y_ed = 20.0\n"

"  elif input_num == '5':\n"
"    print 'Atom:GeSn'\n"
"    y_ed = 12.0\n"

"  else:\n"
"    print 'Input number is wrong!'\n"
"    quit()\n"

"if argv[1].find('fmax') != -1:\n"
"  ftemp = re.findall(r'fmax[0-9]+', argv[1])\n"
"  fmax = re.findall(r'[0-9]+', ftemp[0])\n"
"  y_ed = float(fmax[0]) / 0.0299792458\n"
"  print 'F Max:%f' % y_ed\n"

"if len(argv[7]) > 0:\n"
"  y_ed = float(argv[7]) / 0.0299792458\n"

"if argv[3] != 'Y' and  argv[3] != 'N':\n"
"  print 'Setting of color condition is wrong!'\n"
"  quit()\n"

"print 'Read \"%s\".' % argv[1]\n"

"fi = open(argv[1], 'r')\n"

"x_index = 0\n"        //# number of elements in x axis
"y_index = 0\n"        //# number of elements in y axis
"temp = 0\n"           //# store "x_index" temporarily

//#Read the input data

"for line in fi:\n"
"  item = line.split();\n"

"  if len(item) == 3:\n"
"    if argv[3] == 'Y':\n"
"      density[y_index][x_index] = float(item[2])\n"

"    else:\n"
"      density[y_index][x_index] = - float(item[2])\n"

"    x_index += 1\n"

"  if len(item) == 0:\n"
"    temp = x_index\n"
"    y_index += 1\n"
"    x_index = 0\n"

"fi.close()\n"

"x_index = temp\n"

"x_index = float(x_index)\n"
"y_index = float(y_index)\n"

"dx = (x_ed-x_st)/x_index\n"     //# value of x range
"dy = (y_ed-y_st)/y_index\n"     //# value of y range

"if x_index == 0 or y_index == 0:\n"
"  print \"Setting of index is wrong!\"\n"
"  quit()\n"

"elif x_ed-x_st == 0 or y_ed-y_st == 0:\n"
"  print \"Setting of range is wrong!\"\n"
"  quit()\n"

"else:\n"
"  x = arange(x_st, x_ed+dx, dx)\n"
"  y = arange(y_st, y_ed+dy, dy)\n"

"print 'X Mesh:%d  Y Mesh:%d' % (x_index, y_index)\n"

//# Set option of Matplotlib

"x_index = int(x_index)\n"
"y_index = int(y_index)\n"

"figure(figsize=figaspect(0.45), facecolor=\"w\")\n"
"axes(axisbg='k')\n"
"subplot(1,3,1)\n"
"title('Dispersion', fontsize = 16, fontname = 'serif')\n"
"X, Y = meshgrid(x, y)\n"
"pcolor(X, Y, density)\n"
"colorbar()\n"
"grid(True)\n"
"set_cmap(cmap=plt.get_cmap(argv[2]))\n"

"if c_limit == 'Y':\n"
"  if argv[3] == 'Y':\n"
"    clim(c_st, c_ed)\n"
"  else:\n"
"    clim(-c_ed, c_st)\n"
"else:\n"
"  dmin = 1.0\n"
"  dmax = 0.0\n"
"  for i in range(1, y_index):\n"
"    for j in range(0, x_index):\n"
"      if dmin > density[i][j]:\n"
"        dmin = density[i][j]\n"
"      if dmax < density[i][j]:\n"
"        dmax = density[i][j]\n"
"  c_st = dmin\n"
"  c_ed = dmax*1.05\n"
"  clim(c_st, c_ed)\n"

"print 'Z Min:%e  Z Max:%e' % (c_st, c_ed)\n"

"xlim(x_st, x_ed)\n"
"ylim(y_st, y_ed)\n"

"xticks(fontsize = 14, fontname = 'serif')\n"
"yticks(fontsize = 14, fontname = 'serif')\n"
"xlabel('Wave Vector [($\pi$/a)]', fontsize = 16, fontname = 'serif')\n"
//"ylabel('Wave Number [cm-1]', fontsize = 16, fontname = 'serif')\n"
"ylabel('Phonon Frequency [THz]', fontsize = 16, fontname = 'serif')\n"

"ram = [0 for i in range (y_index+1)]\n"
"sum = [0 for i in range (y_index+1)]\n"
"for i in range(0, y_index):\n"
"  ram[i] = density[i][0]\n"
"  for j in range(0, x_index):\n"
"    sum[i] += density[i][j]\n"
"d_st = c_st\n"
"r_st = c_st\n"
"r_ed = c_ed\n"

"if argv[2].find('Reds') != -1:\n"
"  sumcolor = 'r'\n"
"  ramcolor = 'r'\n"
"elif argv[2].find('Blues') != -1:\n"
"  sumcolor = 'b'\n"
"  ramcolor = 'b'\n"
"elif argv[2].find('Greens') != -1:\n"
"  sumcolor = 'g'\n"
"  ramcolor = 'g'\n"
"elif argv[2].find('Purples') != -1:\n"
"  sumcolor = 'purple'\n"
"  ramcolor = 'purple'\n"
"elif argv[2].find('Oranges') != -1:\n"
"  sumcolor = 'orange'\n"
"  ramcolor = 'orange'\n"
"elif argv[2].find('Greys') != -1:\n"
"  sumcolor = 'grey'\n"
"  ramcolor = 'grey'\n"
"else:\n"
"  sumcolor = 'b'\n"
"  ramcolor = 'r'\n"

"subplot(1,3,2)\n"
"title('DOS Gamma to X', fontsize = 16, fontname = 'serif')\n"
"plot(sum, Y, sumcolor)\n"
"if c_limit == 'Y':\n"
"  smin = 1.0\n"
"  for i in range(1, y_index):\n"
"    if smin > sum[i] and sum[i] > 0.0:\n"
"      smin = sum[i]\n"
"  d_st = smin\n"
"  xlim(d_st, d_ed)\n"
"else:\n"
"  smin = 1.0\n"
"  smax = 0.0\n"
"  for i in range(1, y_index):\n"
"    if smin > sum[i] and sum[i] > 0.0:\n"
"      smin = sum[i]\n"
"    if smax < sum[i]:\n"
"      smax = sum[i]\n"
"  d_st = smin\n"
"  d_ed = smax*1.1\n"
"  xlim(d_st, d_ed)\n"
"print 'D Min:%e  D Max:%e' % (d_st, d_ed)\n"
"ylim(y_st, y_ed)\n"
"yticks(visible=0)\n"
"xlabel('Intensity', fontsize = 16, fontname = 'serif')\n"

"subplot(1,3,3)\n"
"title('DOS only Gammma point', fontsize = 16, fontname = 'serif')\n"
"plot(ram, Y, ramcolor)\n"
"if c_limit == 'Y':\n"
"  rmin = 1.0\n"
"  for i in range(1, y_index):\n"
"    if rmin > ram[i] and ram[i] > 0.0:\n"
"      rmin = ram[i]\n"
"  r_st = rmin\n"
"  xlim(r_st, c_ed)\n"
"else:\n"
"  rmin = 1.0\n"
"  rmax = 0.0\n"
"  for i in range(1, y_index):\n"
"    if rmin > ram[i] and ram[i] > 0.0:\n"
"      rmin = ram[i]\n"
"    if rmax < ram[i]:\n"
"      rmax = ram[i]\n"
"  r_st = rmin\n"
"  r_ed = rmax*1.05\n"
"  xlim(r_st, r_ed)\n"
"print 'R Min:%e  R Max:%e' % (r_st, r_ed)\n"
"ylim(y_st, y_ed)\n"
"yticks(visible=0)\n"
"xlabel('Intensity', fontsize = 16, fontname = 'serif')\n"

"if fig != '':\n"
"  savefig(fig)\n"

"show()\n";

const char *cmap_Py =
"cmaps = [[1, 'rainbow'], [2, 'gist_rainbow'],\n"
"[3, 'spectral'], [4, 'RdBu'],  [5, 'hot'], [6, 'cool'],\n"
"[7, 'Reds'], [8, 'Greens'], [9, 'Blues'], [10, 'Purples'], [11, 'Oranges'], [12, 'Greys'], [13, 'bone']]\n"
"figure(facecolor=\"w\")\n"
"axes(axisbg='k')\n"
"gradient = linspace(0, 1, 256)\n"
"gradient = vstack((gradient, gradient))\n"

"subplot(13,1,1)\n"
"text(-2, 1.2, '(default)', va='center', ha='right', fontsize=10)\n"

"for n, name in cmaps:\n"
"  subplot(13,1,n)\n"
"  imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))\n"
"  axis('off')\n"
"  text(-2, 0.5, name, va='center', ha='right', fontsize=10)\n"

"show()\n";

#define MAXOPT 8

static option options[] =
  {
    {"help", no_argument, NULL, 'h'},
    {"sample", no_argument, NULL, 's'},
    {"reverse", no_argument, NULL, 'r'},
    {"color", required_argument, NULL, 'c'},
    {"limit", required_argument, NULL, 'l'},
    {"minimamlimit", required_argument, NULL, 'm'},
    {"densitylimit", required_argument, NULL, 'd'},
    {"maxhz", required_argument, NULL, 'z'},
    {0, 0, 0, 0}
  };

int main(int argc, char **argv) {

       cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
       cout << "            Phonon DFT plot for Si,Ge,Mo,S version 1.00" << endl << endl;
       cout << "                        Motohiro Tomita 2015" << endl;
       cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  int opt, index;
  int sample = 0;
  char *opts[MAXOPT];
  opts[0] = "";         //name of this program
  opts[1] = "";         //name of graph data
  opts[2] = "rainbow";  //color name
  opts[3] = "Y";        //reverse color
  opts[4] = "";         //color scale limit
  opts[5] = "";         //color scale min limit
  opts[6] = "";         //density limit
  opts[7] = "";         //max hz

  while((opt = getopt_long(argc, argv, "hsrc:m:l:d:z:", options, &index)) != -1){
    switch(opt){
      case 'h':
        cout << "Help of options." << endl;
        cout << "           -s: Show sample color." << endl;
        cout << "           -r: Set reverse color." << endl;
        cout << "           -c: Set color name." << endl;
        cout << "           -l: Set color scale and raman limit." << endl;
        cout << "           -m: Set color scale min limit with -l option." << endl;
        cout << "           -d: Set density limit with -l option." << endl;
        cout << "           -z: Set max frequency." << endl;
        return 0;
      case 's':
        sample = 1;
        break;
      case 'c':
        opts[2] = optarg;
        break;
      case 'r':
        opts[3] = "N";
        break;
      case 'l':
        opts[4] = optarg;
        break;
      case 'm':
        opts[5] = optarg;
        break;
      case 'd':
        opts[6] = optarg;
        break;
      case 'z':
        opts[7] = optarg;
        break;
      default:
        cout << "Wrong option" << endl;
        return 1;
        break;
     }
  }

  opts[0] = argv[0];
  if(argv[optind]) opts[1] = argv[optind];
  if(opts[1] == "" && sample == 0) {
    cout << "no input files" << endl;
    return 1;
  }

  Py_SetProgramName(argv[0]);
  Py_Initialize();
  PySys_SetArgv(MAXOPT, opts);
  PyRun_SimpleString(init_Py);
  if(sample) PyRun_SimpleString(cmap_Py);
  else PyRun_SimpleString(density_Py);
  Py_Finalize();
  return 0;
}

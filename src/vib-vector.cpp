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
"import math\n"
"import matplotlib.pyplot as plt\n"
"import re\n"

"from pylab import *\n"
"from numpy import *\n"

"warnings.simplefilter(\"ignore\", DeprecationWarning)\n";

const char *vib2d_Py =
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

"C = [0 for i in range(MESH_MAX)]\n"
"X = [-255 for i in range(MESH_MAX)]\n"
"Y = [-255 for i in range(MESH_MAX)]\n"
"Z = [-255 for i in range(MESH_MAX)]\n"
"YpZ = [-255 for i in range(MESH_MAX)]\n"
"ZpY = [-255 for i in range(MESH_MAX)]\n"
"U = [0 for i in range(MESH_MAX)]\n"
"V = [0 for i in range(MESH_MAX)]\n"
"W = [0 for i in range(MESH_MAX)]\n"
"R = [0 for i in range(MESH_MAX)]\n"

"argv = sys.argv\n"    //# command line arguments
"iargc = len(argv)\n"  //# number of command line arguments

//#Read the input 1st data

"print 'Read \"%s\".' % argv[1]\n"

"i_index = 0\n"        //# number of elements

"fi = open(argv[1], 'r')\n"

"xmax = 0\n"
"ymax = 0\n"
"zmax = 0\n"
"limit = 2000\n"
"Rlim = 5\n"
"Ylimh = 0.5\n"
"Yliml = -0.3\n"
"Zlimh = 0.5\n"
"Zliml = -0.3\n"
"deg = 1\n"

"for line in fi:\n"
"  item = line.split();\n"
"  if len(item) == 7:\n"
"    if item[0].find('Si') != -1:\n"
"      C[i_index] = 10\n"
"    elif item[0].find('Ge') != -1:\n"
"      C[i_index] = 20\n"
"    else:\n"
"      C[i_index] = 0\n"
"    X[i_index] = float(item[1])\n"
"    Y[i_index] = float(item[2])\n"
"    Z[i_index] = float(item[3])\n"
"    if Y[i_index] < Yliml or Y[i_index] > Ylimh:\n"
"      ZpY[i_index] = -255\n"
"    else:\n"
"      ZpY[i_index] = Z[i_index] + Y[i_index]*math.sin(math.radians(deg))\n"
"    if Z[i_index] < Zliml or Z[i_index] > Zlimh:\n"
"      YpZ[i_index] = -255\n"
"    else:\n"
"      YpZ[i_index] = Y[i_index] + Z[i_index]*math.sin(math.radians(deg))\n"
"    U[i_index] = float(item[4])/limit\n"
"    V[i_index] = float(item[5])/limit\n"
"    W[i_index] = float(item[6])/limit\n"
"    if U[i_index] > Rlim:\n"
"      U[i_index] = Rlim\n"
"    if V[i_index] > Rlim:\n"
"      V[i_index] = Rlim\n"
"    if W[i_index] > Rlim:\n"
"      W[i_index] = Rlim\n"
"    R[i_index] = math.sqrt(U[i_index]*U[i_index] + V[i_index]*V[i_index] + W[i_index]*W[i_index])\n"
"    if xmax < X[i_index]:\n"
"      xmax = X[i_index]\n"
"    if ymax < Y[i_index]:\n"
"      ymax = Y[i_index]\n"
"    if zmax < Z[i_index]:\n"
"      zmax = Z[i_index]\n"
"    i_index += 1\n"

"fi.close()\n"

"i_index = float(i_index)\n"

"print 'Atom number:%d' % (i_index)\n"

//# Set option of Matplotlib

"i_index = int(i_index)\n"

"ratio = (ymax*2)/xmax\n"
"figure(figsize=figaspect(ratio), facecolor=\"w\")\n"
"axes(axisbg='k')\n"
"subplot(2,1,1)\n"
"scatter(X, YpZ, c=C, cmap=\"rainbow\", s=300)\n"
"quiver(X, YpZ, U, V, R, cmap=\"rainbow\", width=0.01)\n"
"xlabel('X(A)', fontsize = 16, fontname = 'serif')\n"
"ylabel('Y(A)', fontsize = 16, fontname = 'serif')\n"
"xlim(-xmax*0.1, xmax*1.1)\n"
"ylim(-ymax*1.1, ymax*1.1)\n"

"subplot(2,1,2)\n"
"scatter(X, ZpY, c=C, cmap=\"rainbow\", s=300)\n"
"quiver(X, ZpY, U, W, R, cmap=\"rainbow\", width=0.01)\n"
"xlabel('X(A)', fontsize = 16, fontname = 'serif')\n"
"ylabel('Z(A)', fontsize = 16, fontname = 'serif')\n"
"xlim(-xmax*0.1, xmax*1.1)\n"
"ylim(-zmax*1.1, zmax*1.1)\n"

"show()\n";

const char *test_Py =
"fig = plt.figure()\n"
"ax = fig.gca(projection='3d')\n"

"x, y, z = np.meshgrid(np.arange(-0.8, 1, 0.2),\n"
"                      np.arange(-0.8, 1, 0.2),\n"
"                      np.arange(-0.8, 1, 0.8))\n"

"u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)\n"
"v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)\n"
"w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *\n"
"     np.sin(np.pi * z))\n"

"ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)\n"

"plt.show()\n";

#define MAXOPT 9

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
  opts[8] = "";         //name of graph data2

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
  //PyRun_SimpleString(test_Py);
  PyRun_SimpleString(vib2d_Py);
  Py_Finalize();
  return 0;
}

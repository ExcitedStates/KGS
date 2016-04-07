/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

        Permission is hereby granted, free of charge, to any person obtaining a copy of
        this software and associated documentation files (the "Software"), to deal in
        the Software without restriction, including without limitation the rights to
        use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
        of the Software, and to permit persons to whom the Software is furnished to do
        so, subject to the following conditions:

        This entire text, including the above copyright notice and this permission notice
        shall be included in all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
        OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
        FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
        IN THE SOFTWARE.


*/
#include "RunFirst.h"
//#include "Froda.h"
#include <iostream>
#include <cstdlib>
#include <stdio.h>

using namespace std;


/* Use the unix 'which' to locate the FIRST executable and return the path. Exit with an error if it can't be found. */
string RunFirst::FirstPath(){
	//From http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c 
	FILE* pipe = popen("which FIRST", "r"); 
	char buffer[128]; 
	string result = "";
	while(!feof(pipe)) { 
		if(fgets(buffer, 128, pipe) != NULL) result += buffer; 
	} 
	pclose(pipe); 

	if(result==""){ 
		cerr<<"FIRST executable could not be located. Make sure its specified in PATH"<<endl;
		exit(-1); 
	} 
	return result.substr(0,result.rfind("/")+1); 
}

void RunFirst::RigidityAnalysis (string protein_file_name) {
	cout<<"RunFirst::RigidityAnalysis(...)"<<endl;
//        string command = "/home/vdbedem/DEV/ModelBuildTool/MBT/CTK/KGSX/trunk/FIRST-6.2.1-bin-64-gcc3.4.3-O3/myfirst "+protein_file_name+" -non -L /home/vdbedem/DEV/ModelBuildTool/MBT/CTK/FIRST-6.2.1-bin-64-gcc3.4.3-O3 -E -5.00 -covin -hbin -phin -minRC 1";
//        string command = "/Users/vdbedem/DEV/KGSX/trunk/FIRST-6.2.1-bin-64-gcc3.4.3-O3/tmp/FIRST_peggy/src/FIRST "+protein_file_name+" -non -L ../FIRST-6.2.1-bin-64-gcc3.4.3-O3 -E -5.00 -covin -hbin -phin -minRC 1";
//        string command = "/Users/rfonseca/Documents/RNAKGS/KGSX/trunk/FIRST-macbook/FIRST "+protein_file_name+" -non -L ../FIRST-macbook -E -5.00 -covin -hbin -phin -minRC 1";
//        string command = "/home/dpachov/mymodules/FIRST/6.2.1-GCC3.4.3/myfirst "+protein_file_name+" -non -L /home/dpachov/mymodules/FIRST/6.2.1-GCC3.4.3 -E -5.00 -covin -hbin -phin -minRC 1";

	string firstPath = FirstPath();
	string command = firstPath+"FIRST "+protein_file_name+" -non -L "+firstPath+" -E -1.00 -covin -hbin  -minRC 1";
	int status = system(command.c_str());

	if(status){ cerr<<"FIRST exited with error"<<endl;	exit(-1);	}

/*	// argv = FIRST protein_file_name -non -E 100 -covin -hbin -phin -minRC 1
	char* argv[10];
	argv[0] = "FIRST";
	argv[1] = (char*) protein_file_name.c_str();
	argv[2] = "-non";
	argv[3] = "-E";
	argv[4] = "100";
	argv[5] = "-covin";
	argv[6] = "-hbin";
	argv[7] = "-phin";
	argv[8] = "-minRC";
	argv[9] = "1";
	main_FIRST(10,argv);
*/
}

void RunFirst::Froda (string protein_file_name) {
	cout<<"RunFirst::Froda(...)"<<endl;
//	string command = "/home/vdbedem/DEV/ModelBuildTool/MBT/CTK/KGSX/trunk/FIRST-6.2.1-bin-64-gcc3.4.3-O3/myfirst "+protein_file_name+" -non -L /home/vdbedem/DEV/ModelBuildTool/MBT/CTK/FIRST-6.2.1-bin-64-gcc3.4.3-O3 -E -5.00 -covin -hbin -phin -minRC 1 -FRODA -body -CMforward 0.2 -freq 100 -vdw 0.65 -totconf 1000";
// string command = "/Users/vdbedem/DEV/KGSX/trunk/FIRST-6.2.1-bin-64-gcc3.4.3-O3/tmp/FIRST_peggy/src/FIRST "+protein_file_name+" -non -L ../FIRST-6.2.1-bin-64-gcc3.4.3-O3 -E -5.00 -covin -hbin -phin -minRC 1 -FRODA -body -CMforward 0.2 -freq 100 -vdw 0.65 -totconf 1000";
// string command = "/Users/rfonseca/Documents/RNAKGS/KGSX/trunk/FIRST-macbook/FIRST "+protein_file_name+" -non -L ../FIRST-macbook -E -5.00 -covin -hbin -phin -minRC 1 -FRODA -body -CMforward 0.2 -freq 100 -vdw 0.65 -totconf 1000";
//	string command = "/home/dpachov/mymodules/FIRST/6.2.1-GCC3.4.3/myfirst "+protein_file_name+" -non -L /home/dpachov/mymodules/FIRST/6.2.1-GCC3.4.3 -E -5.00 -covin -hbin -phin -minRC 1 -FRODA -body -CMforward 0.2 -freq 100 -vdw 0.65 -totconf 1000";

	string firstPath = FirstPath();
	string command = firstPath+"FIRST "+protein_file_name+" -non -L "+firstPath+" -E -1.00 -covin -hbin -minRC 1 -FRODA -body -CMforward 0.2 -freq 100 -vdw 0.65 -totconf 1000";
	int status = system(command.c_str());

	if(status){ cerr<<"FIRST exited with error"<<endl;	exit(-1);	}

/*	// argv = FIRST protein_file_name -non -E 100 -covin -hbin -phin -minRC 1 
	//        -FRODA -body -CMforward 0.2 -freq 100 -vdw 0.65 -totconf 1000
	char* argv[20];
	argv[0] = "FIRST";
	argv[1] = (char*) protein_file_name.c_str();
	argv[2] = "-non";
	argv[3] = "-E";
	argv[4] = "100";
	argv[5] = "-covin";
	argv[6] = "-hbin";
	argv[7] = "-phin";
	argv[8] = "-minRC";
	argv[9] = "1";
	argv[10] = "-FRODA";
	argv[11] = "-body";
	argv[12] = "-CMforward";
	argv[13] = "0.2"; // default value in FIRST
	argv[14] = "-freq";
	argv[15] = "100";
	argv[16] = "-vdw";
	argv[17] = "0.65";
	argv[18] = "-totconf";
	argv[19] = "1000";
	main_FIRST(20,argv);
*/
}

void RunFirst::ComputeHbondEnergy (string path, string protein_name, string input_file, string output_file) {
	cout<<"RunFirst::ComputeHbondEnergy(...)"<<endl;
	string command;

	if (input_file!="hbonds.in") {
		command = "cp "+path+input_file+" "+path+"hbonds.in";
		system(command.c_str());
	}

//	command = "/home/vdbedem/DEV/ModelBuildTool/MBT/CTK/KGSX/trunk/FIRST-6.2.1-bin-64-gcc3.4.3-O3/FIRST "+path+protein_name+".pdb -non -L home/vdbedem/DEV/ModelBuildTool/MBT/CTK/FIRST-6.2.1-bin-64-gcc3.4.3-O3 -E -1 -covin -phin -hbin -hbout";
// command = "/Users/vdbedem/DEV/KGSX/trunk/FIRST-6.2.1-bin-64-gcc3.4.3-O3/tmp/FIRST_peggy/src/FIRST "+path+protein_name+".pdb -non -L ../FIRST-6.2.1-bin-64-gcc3.4.3-O3 -E -1 -covin -phin -hbin -hbout";
// command = "/Users/rfonseca/Documents/RNAKGS/KGSX/trunk/FIRST-macbook/FIRST "+path+protein_name+".pdb -non -L ../FIRST-macbook -E -1 -covin -phin -hbin -hbout";
//	command = "/home/dpachov/mymodules/FIRST/6.2.1-GCC3.4.3/FIRST "+path+protein_name+".pdb -non -L /home/dpachov/mymodules/FIRST/6.2.1-GCC3.4.3 -E -1 -covin -phin -hbin -hbout";

	string firstPath = FirstPath();
	command = firstPath+"FIRST "+path+protein_name+".pdb -non -L "+firstPath+" -E -1 -covin -hbin -hbout";
	int status = system(command.c_str());

	if(status){ cerr<<"FIRST exited with error"<<endl;	exit(-1);	}

	if (output_file!="hbonds.out") {
		command = "mv "+path+"hbonds.out "+path+output_file;
		system(command.c_str());
	}
}

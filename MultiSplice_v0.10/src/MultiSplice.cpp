/*    
 *    MultiSplice.cpp		
 *    MultiSplice
 *
 *    Copyright (C) 2012 University of Kentucky and
 *                       Yan Huang
 *
 *    Authors: Yan Huang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define UNIX

#ifdef UNIX


#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <sstream>
#include <cmath>
#include <ctime>
#include <assert.h>
#include <algorithm>

#else
#include <fstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <assert.h>
#include <sstream>
#include <time.h>
#include <assert.h>
#include <algorithm>

#endif


using namespace std;


const int default_dataset_num = 100;

string itostr(long t);
void inputChrName(string namefile);
vector<string> chrName;


void processResult(string infilename, string outfilename, string datapath)
{
	string chr, genename, ensemlIden, translen, info, filename;
	vector<double> abundance;
	ifstream infile, file;
	infile.open(infilename.c_str());
	ofstream outfile;
	outfile.open(outfilename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

	while(infile >> chr)
	{
		infile >> genename;
		infile >> ensemlIden;
		infile >> translen;
		getline(infile, info);

		outfile << chr << "\t" << genename << "\t" << ensemlIden << "\t" << translen << "\t";
		if (abundance.size() > 0)
		{
			outfile << abundance[0] << endl;
			abundance.erase(abundance.begin());
		}
		else
		{
			filename = datapath + "result_" + genename + ".txt";
			file.open(filename.c_str());
			while(file >> info)
			{
				abundance.push_back(atof(info.c_str()));
				getline(file, info);
			}
			file.close();
			outfile << abundance[0] << endl;
			abundance.erase(abundance.begin());
		}
	}

	infile.close();
	outfile.close();

	return;
}


int main(int argc, char* argv[])
{
	if (argc != 6)
	{
		cout << argv[0] << "<MultiSplicePath>\t<GTFFileName>\t<SAMFileName>\t<OutputFilePath>\t<FragmentLength>" << endl;
		return 1;
	}

	string comd, srcpath, filepath, filename;
	int i;
	srcpath = argv[1];


	/************************************************************************/
	// temporary file folder	
	filepath = argv[4];
	comd = "mkdir " + filepath + "tmp_foler";
	i = system(comd.c_str());


	/************************************************************************/
	// separate GTF file by chromosomes
	comd = "mkdir " + filepath + "tmp_foler/GTF";
	i = system(comd.c_str());

	filename = argv[2];
	comd = srcpath + "bin/SeparateChrFromGTF " + filename + " " + filepath + "tmp_foler/GTF/";
	i = system(comd.c_str());


	/************************************************************************/
	// separate GTF by gene
	comd = "mkdir " + filepath + "tmp_foler/GTF/Genes";
	i = system(comd.c_str());

	filename = filepath + "tmp_foler/GTF/ChromosomeName.txt";
	comd = srcpath + "bin/SeparateGeneFromGTF " + filename + " " + filepath + "tmp_foler/GTF/";
	i = system(comd.c_str());

	cout << "Parse GTFile done......" << endl;

	/************************************************************************/
	// parse SAM file
	comd = "mkdir " + filepath + "tmp_foler/SAM";
	i = system(comd.c_str());

	filename = argv[3];
	comd = srcpath + "bin/SeparateChrFromSAM " + filename + " " + filepath + "tmp_foler/SAM/";
	i = system(comd.c_str());

	/************************************************************************/
	// parse fragments
	comd = "mkdir " + filepath + "tmp_foler/Frag";
	i = system(comd.c_str());
	filename = filepath + "tmp_foler/SAM/ChromosomeName.txt";
	inputChrName(filename);

	for (int chrIndex = 0; chrIndex < chrName.size(); chrIndex++)
	{
		comd = "mkdir " + filepath + "tmp_foler/Frag/" + chrName[chrIndex];
		i = system(comd.c_str());
		comd = srcpath + "bin/fragment " + filepath + "tmp_foler/SAM/ " + chrName[chrIndex] + " 1 noDB " + filepath + "tmp_foler/Frag/" + chrName[chrIndex] + "/";
	    i = system(comd.c_str());
	}
	chrName.clear();

	cout << "Parse SAMFile done......" << endl;

	/************************************************************************/
	// run MultiSplice
	comd = "mkdir " + filepath + "tmp_foler/Result";
	i = system(comd.c_str());
	inputChrName(filename);
	string fraglength = argv[5];

	for (int chrIndex = 0; chrIndex < chrName.size(); chrIndex++)
	{
		comd = "mkdir " + filepath + "tmp_foler/Result/" + chrName[chrIndex];
		i = system(comd.c_str());
		comd = "mkdir " + filepath + "tmp_foler/Result/" + chrName[chrIndex] + "/temp";
		i = system(comd.c_str());
		
		comd = srcpath + "bin/Estimate " + filepath + "tmp_foler/GTF/Genes/" + chrName[chrIndex] + "/" + chrName[chrIndex] + "W/ ";
		comd += filepath + "tmp_foler/Frag/" + chrName[chrIndex] + "/ ";
		comd += filepath + "tmp_foler/Result/" + chrName[chrIndex] + "/ " + fraglength + " " + fraglength + " ";
		comd += chrName[chrIndex] + " noSeq " + fraglength + " 1";
		i = system(comd.c_str());

		comd = srcpath + "bin/Estimate " + filepath + "tmp_foler/GTF/Genes/" + chrName[chrIndex] + "/" + chrName[chrIndex] + "C/ ";
		comd += filepath + "tmp_foler/Frag/" + chrName[chrIndex] + "/ ";
		comd += filepath + "tmp_foler/Result/" + chrName[chrIndex] + "/ " + fraglength + " " + fraglength + " ";
		comd += chrName[chrIndex] + " noSeq " + fraglength + " 1";
		i = system(comd.c_str());
	}
	chrName.clear();

	/************************************************************************/
	// run Matlab
	char chrcomd[5000];
	string tmpstr = srcpath + "bin/";
	const char *srcpath_matlab = tmpstr.c_str();
	tmpstr = filepath + "tmp_foler/Result/";
	const char *OutputFilePath = tmpstr.c_str();
	tmpstr = filepath + "tmp_foler/GTF/";
	const char *path = tmpstr.c_str();
	sprintf(chrcomd,  "matlab -r \"%sEstimate('%sChromosomeName.txt', '%sGenes/', '%s');exit;\"", srcpath_matlab, path, path, OutputFilePath);
	i = system(chrcomd);

	/************************************************************************/
	// finalize results
	filename = filepath + "tmp_foler/GTF/ChromosomeName.txt";
	inputChrName(filename);
	string outfilename = filepath + "MultiSplice_estimation.txt", tmpfilepath;
	for (int chrIndex = 0; chrIndex < chrName.size(); chrIndex++)
	{
		filename = filepath + "tmp_foler/Result/" + chrName[chrIndex] + "/" + chrName[chrIndex] + ".pos";
		tmpfilepath = filepath + "tmp_foler/Result/" + chrName[chrIndex] + "/temp/";
		processResult(filename, outfilename, tmpfilepath);
	}
	chrName.clear();

	/************************************************************************/
	// clean temporary files
	comd = "rm -r " + filepath + "tmp_foler";
	i = system(comd.c_str());

	return 0;
}

void inputChrName(string namefile)
{
	ifstream chrFile;
	string chrNm, info;
	chrFile.open(namefile.c_str());

	while(chrFile >> chrNm)
	{
		getline(chrFile, info);
		if (chrNm[0] == 'c' && chrNm[1] == 'h' && chrNm[2] == 'r')
		{
			if (chrName.size() >= chrName.capacity())
				chrName.reserve(default_dataset_num + chrName.capacity());
			chrName.push_back(chrNm);
		}

	}
	chrFile.close();

	return;
}

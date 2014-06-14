/*    
 *    SeparateGeneFromGTF.cpp		
 *    SeparateGeneFromGTF
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

#endif

using namespace std;

// separate GTF file into genes
void processGTF(string chrfilename, string path)
{
	// info on a line of GTF
	string chr, seqname, source, feature, genename, genesymbol, start, end, strand, info, comd;
	int i;

	ifstream chrfile, gtffile;
	ofstream genelistfile, outfile;
	string GeneList, GTFfilename, outfilename;
	long transCnt = 0;
	size_t found;

	chrfile.open(chrfilename.c_str());
	while(chrfile >> chr)
	{
		getline(chrfile, info);
		comd = "mkdir " + path + "Genes/" + chr;
		i = system(comd.c_str());
		GeneList = path + "Genes/" + chr + "/GeneList.txt";
		genelistfile.open(GeneList.c_str());

		transCnt = 0;
		GTFfilename = path + chr + ".gtf";
		gtffile.open(GTFfilename.c_str());
		while(gtffile >> seqname)
		{
			gtffile >> source;
			gtffile >> feature;
			if (feature.compare("gene") == 0)
			{
				// a new gene
				if (transCnt > 0)
				{
					genelistfile << transCnt << endl;
					outfile.close();
				}

				gtffile >> start;
				gtffile >> end;
				gtffile >> info;
				gtffile >> strand;
				if (strand.compare("+") == 0)
				{
					outfilename = path + "Genes/" + chr + "/" + start + "-" + end + "W";
					genelistfile << start << "-" << end << "W\t";
				}
				else
				{
					outfilename = path + "Genes/" + chr + "/" + start + "-" + end + "C";
					genelistfile << start << "-" << end << "C\t";
				} 
				outfile.open(outfilename.c_str());
				outfile << chr << "\t" << source << "\t" << feature << "\t" << start << "\t" << end << "\t" << info << "\t" << strand;
				getline(gtffile, info);
				outfile << info << endl;

				// extract the gene names and gene symbols
				genename = "NA";
				genesymbol = "NA";
				found = info.find("\"");
				if (found != string::npos)
				{
					info = info.substr(found + 1);
				}
				found = info.find("\"");
				if (found != string::npos)
				{
					genename = info.substr(0, found);
				}
				info = info.substr(found + 1);
				found = info.find("\"");
				if (found != string::npos)
				{
					info = info.substr(found + 1);
				}
				found = info.find("\"");
				if (found != string::npos)
				{
					genesymbol = info.substr(0, found);
				}
				genelistfile << genename << "\t" << genesymbol << "\t";
				transCnt = 0;
			}
			else if (feature.compare("transcript") == 0)
			{
				getline(gtffile, info);
				outfile << chr << "\t" << source << "\t" << feature << info << endl;
				transCnt++;
			}
			else
			{
				getline(gtffile, info);
				outfile << chr << "\t" << source << "\t" << feature << info << endl;
			}
		}
		// deal with the last gene
		genelistfile << transCnt << endl;
		outfile.close();
		gtffile.close();
		genelistfile.close();
	}
	chrfile.close();

	return;
}

// put positive and negative strand gtf in separate folders
void separateGTF(string chrfilename, string path)
{
	string chr, info, genename, comd, filename;
	int i;

	ifstream chrfile, genefile;
	fstream inputfile;
	chrfile.open(chrfilename.c_str());
	while(chrfile >> chr)
	{
		comd = "mkdir " + path + "Genes/" + chr + "/" + chr + "W";
		i = system(comd.c_str());
		comd = "mkdir " + path + "Genes/" + chr + "/" + chr + "C";
		i = system(comd.c_str());

		filename = path + "Genes/" + chr + "/GeneList.txt";
		genefile.open(filename.c_str());
		while(genefile >> genename)
		{
			getline(genefile, info);
			if (genename[genename.length()-1] == 'W')
			{
				comd = "mv " + path + "Genes/" + chr + "/" + genename + " " + path + "Genes/" + chr + "/" + chr + "W/";
				filename = path + "Genes/" + chr + "/" + chr + "W/GeneList.txt";
				inputfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
				inputfile << genename << endl;
				inputfile.close();
			}
			else
			{
				comd = "mv " + path + "Genes/" + chr + "/" + genename + " " + path + "Genes/" + chr + "/" + chr + "C/";
				filename = path + "Genes/" + chr + "/" + chr + "C/GeneList.txt";
				inputfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
				inputfile << genename << endl;
				inputfile.close();
			}
			
			i = system(comd.c_str());
		}
		genefile.close();

		getline(chrfile, info);

	}

	return;
}


int main (int argc, char *argv[])
{
	if (argc != 3)
	{
		cout << argv[0] << "\t<chrFile>" << "\t<Input_Path>" << endl;
		return 1;
	}
	processGTF(argv[1], argv[2]);
	separateGTF(argv[1], argv[2]);

	return 0;
}
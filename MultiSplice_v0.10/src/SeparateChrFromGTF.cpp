/*    
 *    SeparateChrFromGTF.cpp		
 *    SeparateChrFromGTF
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
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <assert.h>
#include <algorithm>

#else

#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
#include <assert.h>
#include <algorithm>

#endif


using namespace std;


/************************************************************************/
/* Separate GTF file by chromosome */
/************************************************************************/

// Input a GTF file, separate it by chromsome names
void Parse(char *inputfilename, char *outputfile_path)
{
	ifstream inputfile;
	inputfile.open(inputfilename);
	fstream outputfile_gtf;
	ofstream outputfile;
	
	char chromsome[100], outputfilename_gtf[1000], outputfilename[1000], chrName[1000][100];
	string info;
	int tmp, iLoop, chrNm;
	bool found;
	
	chrNm = 0;
	for (tmp = 0; tmp < 100; tmp++)
	{
		chromsome[tmp] = '\0';
		for (iLoop = 0; iLoop < 100; iLoop++)
		{
			chrName[tmp][iLoop] = '\0';
		}
	}
	inputfile >> chromsome;
	while(chromsome[0] != '\0')
	{
		getline(inputfile, info);
		
		// Separate GTF file by its chromosome
		if (chromsome[0] != '#')
		{
			sprintf(outputfilename_gtf, "%schr%s.gtf", outputfile_path, chromsome);
			outputfile_gtf.open (outputfilename_gtf, fstream::in | fstream::out | fstream::app);
			outputfile_gtf << chromsome << info << endl;
			outputfile_gtf.close();

			found = false;
			for (tmp = 1; tmp <= chrNm; tmp++)
			{
				if (strcmp(chrName[tmp], chromsome) == 0)
				{
					found = true;
					break;
				}
			}
			if (found == false)
			{
				chrNm++;
				strcpy(chrName[chrNm], chromsome);
			}
		}

		chromsome[0] = '\0';
		inputfile >> chromsome;
	}

	sprintf(outputfilename, "%sChromosomeName.txt", outputfile_path);
	outputfile.open(outputfilename);
	for (tmp = 1; tmp <= chrNm; tmp++)
	{
		outputfile << "chr" << chrName[tmp] << endl;
	}

	inputfile.close();
	outputfile.close();
	return;
}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		cout << argv[0] << "\t<gtfFile>" << "\t<Output_Path>" << endl;
		return 1;
	}
	Parse(argv[1], argv[2]);
	return 0;
}


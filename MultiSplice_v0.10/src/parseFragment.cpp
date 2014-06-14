/*    
 *    parseFragment.cpp		
 *    fragment
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

string chrName;
const long MAX = 2000000000;

ofstream exonFile;
ofstream junctionFile;
ofstream readFile;
ofstream blockFile;
fstream statFile;
fstream datastatFile;
fstream crossGeneFile;

fstream allExonFile;
fstream allJunctionFile;


class transcript
{
public:
	char transID[100];
	long start;
	long end;

	transcript *next;

	transcript();
};

class gene
{
public:
	char geneID[100];
	char chrNm[100];
	long start;
	long end;

	long db_trans_cnt;

	transcript *transList;

	gene();
};

gene* geneList[1000000];
long geneListNum = 0;

long sortKey[1000000];
long mergeSort_Larray[1000000];
long mergeSort_Rarray[1000000];
gene* mergeSort_LorderedList[1000000];
gene* mergeSort_RorderedList[1000000];

char dirPrefix[400];


enum read_type {read_ingene, read_crossgene, read_betweengene, read_crossbetweengene, read_crossIngeneBetweengene};


transcript::transcript()
{
	for (int i = 0; i < 100; i++)
	{
		transID[i] = '\0';
	}

	start = MAX;
	end = 0;

	next = NULL;
}

gene::gene()
{
	for (int i = 0; i < 100; i++)
	{
		geneID[i] = '\0';
		chrNm[i] = '\0';
	}

	start = MAX;
	end = 0;

	db_trans_cnt = 0;

	transList = NULL;
}


long compute_endpoint(long startPosition, char *end)
{
	long i, tmp = 0, endPosition;

	endPosition = startPosition - 1;

	for (i = 0; end[i] != '\0'; i++)
	{
		if (end[i] == 'M')
		{
			endPosition = endPosition + tmp;
			tmp = 0;
		} 
		else if (end[i] == 'N')
		{
			endPosition = endPosition + tmp;
			tmp = 0;
		}
		else if (end[i] >= '0' && end[i] <= '9')
		{
			tmp = tmp * 10 + end[i] - 48;
		}
		else
		{
			tmp = 0;
		}
	}

	return endPosition;
}

//sort fragment
void merge(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortKey[p + i - 1];
		mergeSort_LorderedList[i] = geneList[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortKey[q + j];
		mergeSort_RorderedList[j] = geneList[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX;
	mergeSort_Rarray[n2 + 1] = MAX;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortKey[k] = mergeSort_Larray[i];
			geneList[k] = mergeSort_LorderedList[i];

			i++;
		} 
		else
		{
			sortKey[k] = mergeSort_Rarray[j];
			geneList[k] = mergeSort_RorderedList[j];

			j++;
		}
	}

	return;
}


void mergeSort(long sortList_size)
{
	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = sortList_size;

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}

void parseDB(char* inputfilename)
{
	gene *newGene;
	newGene = NULL;

	ifstream inputfile;
	inputfile.open(inputfilename);

	char chrNm[100], lineCategory[100], IDinfo[500], tmpChar[100];
	int tmp;
	long start, end, iLoop;
	string otherInfo;
	transcript *newTrans;

	for (tmp = 0; tmp < 100; tmp++)
	{
		chrNm[tmp] = '\0';
	}

	inputfile >> chrNm;
	while (chrNm[0] != '\0')
	{
		inputfile >> tmpChar;
		inputfile >> lineCategory;
		inputfile >> start;
		inputfile >> end;
		inputfile >> tmpChar;
		inputfile >> tmpChar;
		inputfile >> tmpChar;
		inputfile >> IDinfo;
		getline(inputfile, otherInfo);

		if (strcmp(lineCategory, "gene") == 0)
		{
			newGene = new gene;
			strcpy(newGene->chrNm, chrNm);
			newGene->start = start;
			newGene->end = end;
			newGene->db_trans_cnt = 0;

			if (IDinfo[0] == 'I' && IDinfo[1] == 'D' && IDinfo[2] == '=')
			{
				for (tmp = 0; IDinfo[tmp + 3] != '|' && IDinfo[tmp + 3] != '\0'; tmp++)
				{
					newGene->geneID[tmp] = IDinfo[tmp + 3];
				}
				newGene->geneID[tmp] = '\0';
			} 
			else
			{
				cout << "Error: abnormal gene ID. Please confirm... ";
				cin >> tmpChar;
				exit(1);
			}

			geneListNum++;
			geneList[geneListNum] = newGene;
		}
		else if (strcmp(lineCategory, "transcript") == 0)
		{
			newTrans = new transcript;

			newGene->db_trans_cnt += 1;

			newTrans->start = start;
			newTrans->end = end;
			newTrans->next = newGene->transList;
			newGene->transList = newTrans; 

			if (IDinfo[0] == 'I' && IDinfo[1] == 'D' && IDinfo[2] == '=')
			{
				for (tmp = 0; IDinfo[tmp + 3] != ';' && IDinfo[tmp + 3] != '|' && IDinfo[tmp + 3] != '\0'; tmp++)
				{
					newTrans->transID[tmp] = IDinfo[tmp + 3];
				}
				newTrans->transID[tmp] = '\0';
			} 
			else
			{
				cout << "Error: abnormal gene ID. Please confirm... ";
				cin >> tmpChar;
				exit(1);
			}
		}
		else
		{
			//do nothing
		}

		chrNm[0] = '\0';
		inputfile >> chrNm;		
	}

	inputfile.close();

	for (iLoop = 1; iLoop <= geneListNum; iLoop++)
	{
		sortKey[iLoop] = geneList[iLoop]->end;
	}
	mergeSort(geneListNum);

	for (iLoop = 1; iLoop <= geneListNum; iLoop++)
	{
		sortKey[iLoop] = geneList[iLoop]->start;
	}
	mergeSort(geneListNum);

	return;
}

double getGeneIndex(long position)
{
	double posIndex = 0;
	long iLoop;

	if (position < geneList[1]->start)
	{
		posIndex = 0.5;
	}
	else if (position > geneList[geneListNum]->end)
	{
		posIndex = geneListNum + 0.5;
	}
	else
	{
		for (iLoop = 1; iLoop < geneListNum; iLoop++)
		{
			if (geneList[iLoop]->start <= position && position <= geneList[iLoop]->end)
			{
				posIndex = iLoop;
			}
			else if (geneList[iLoop]->end < position && position < geneList[iLoop + 1]->start)
			{
				posIndex = iLoop + 0.5;
			}
		}
		if (geneList[geneListNum]->start <= position && position <= geneList[geneListNum]->end)
		{
			posIndex = geneListNum;
		}
	}

	return posIndex;
}

read_type crossGene(long startPosition, long endPosition)
{
	//return true if the read is cross gene
	long iLoop;
	double startIndex = 0, endIndex = 0;

	startIndex = getGeneIndex(startPosition);
	endIndex = getGeneIndex(endPosition);

	if (startIndex == endIndex)
	{
		if (startIndex == int(startIndex))
		{
			return read_ingene;
		} 
		else
		{
			return read_betweengene;
		}
	}
	else
	{
		if (startIndex == int(startIndex) && endIndex == int(endIndex))
		{
			return read_crossgene;
		}
		else if (startIndex != int(startIndex) && endIndex != int(endIndex))
		{
			return read_crossbetweengene;
		}
		else
		{
			return read_crossIngeneBetweengene;
		}
	}
}

void parse(char* inputfilename)
{
	ifstream inputfile;
	inputfile.open(inputfilename);

	string info, curLine, XS_field;

	char end[2000];
	long startPoint, endPoint, tmp, i, sign, totalReadNum = 0, crossGeneReadNum = 0, cnt_read_ingene = 0, cnt_read_crossgene = 0, cnt_read_betweengene = 0, cnt_read_crossbetweengene = 0, cnt_read_crossIngeneBetweengene = 0;
	bool spliced, flag, transDirection;
	read_type curReadType;

	long stat_MinStart = MAX, stat_MaxEnd = 0;


	flag = false; // indicates whether we need to start a new line for blocks
	while (inputfile >> info)
	{
		startPoint = atol(info.c_str());
		inputfile >> end;
		inputfile >> XS_field;
		getline(inputfile, info);

		if (end[0] != '*')
		{
//			curReadType = crossGene(startPoint, compute_endpoint(startPoint, end));
			curReadType = read_ingene;

			totalReadNum++;

// 			if (curReadType == read_ingene)
// 			{
// 				cnt_read_ingene++;
// 			}
// 			else if (curReadType == read_crossgene)
// 			{
// 				cnt_read_crossgene++;
// 			}
// 			else if (curReadType == read_betweengene)
// 			{
// 				cnt_read_betweengene++;
// 			}
// 			else if (curReadType == read_crossbetweengene)
// 			{
// 				cnt_read_crossbetweengene++;
// 			}
// 			else if (curReadType == read_crossIngeneBetweengene)
// 			{
// 				cnt_read_crossIngeneBetweengene++;
// 			}

			if (curReadType == read_crossgene || curReadType == read_crossbetweengene || curReadType == read_crossIngeneBetweengene)
			{
				crossGeneReadNum++;
//				crossGeneFile << name << "\t" << field1 << "\t" << field2 << "\t" << startPoint << "\t" << field4 << "\t" << end << endl;
			}
			else if (curReadType == read_ingene || curReadType == read_betweengene)
			{
				if (flag == true)
				{
					blockFile << endl;
					flag = false;
				}

				tmp = 0;
				sign = 1;
				spliced = false;

				if (XS_field.compare("-") == 0)
				{
					transDirection = false;
				} 
				else
				{
					transDirection = true;
				}

				if (startPoint < stat_MinStart)
					stat_MinStart = startPoint;

				for (i = 0; end[i] != '\0'; i++)
				{
					if (end[i] == 'N')
					{
						spliced = true;
						blockFile << chrName << "\t";
						break;
					}
				}

				readFile << chrName << "\t" << startPoint << '\t' << startPoint << '\t' << spliced << endl;

				for (i = 0; end[i] != '\0'; i++)
				{
					if (end[i] == 'M')
					{
						endPoint = startPoint + tmp * sign;
						if (startPoint < 0)
						{
							exit(1);
						}
						exonFile << chrName << "\t" << startPoint << '\t' << endPoint - 1 << '\t' << spliced << endl;
						allExonFile << chrName << "\t" << startPoint << '\t' << endPoint - 1 << '\t' << spliced << endl;
						startPoint = endPoint;
						tmp = 0;
						sign = 1;
					} 
					else if (end[i] == 'N')
					{
						endPoint = startPoint + tmp * sign;
						if (startPoint < 0)
						{
							exit(1);
						}
						junctionFile << chrName << "\t" << startPoint - 1 << '\t' << endPoint << "\t" << transDirection << endl;
						blockFile << startPoint - 1 << '\t' << endPoint << '\t'; 
						flag = true;
						allJunctionFile << chrName << "\t" << startPoint - 1 << '\t' << endPoint << "\t" << transDirection << endl;
						startPoint = endPoint;
						tmp = 0;
						sign = 1;
					}
					else if (end[i] == '-')
					{
						exit(1);
					}
					else if (end[i] >= '0' && end[i] <= '9')
					{
						tmp = tmp * 10 + end[i] - 48;
					}
					else
					{
						tmp = 0;
						sign = 1;
					}
				}

				if (endPoint > stat_MaxEnd)
					stat_MaxEnd = endPoint;
			}
		}
	}
	
	statFile << stat_MinStart << '\t' << stat_MaxEnd << endl;
	datastatFile << totalReadNum << "\t" << crossGeneReadNum << "\t" << cnt_read_ingene << "\t" << cnt_read_crossgene << "\t" << cnt_read_betweengene << "\t" << cnt_read_crossbetweengene << "\t" << cnt_read_crossIngeneBetweengene << endl;

	return;
}



int main(int argc, char* argv[])
{
	if (argc == 5)
	{
		sprintf(dirPrefix, "./tmp/");
	}
	else if (argc == 6)
	{
		strcpy(dirPrefix, argv[5]);
	}
	else
	{
		cout << argv[0] << "\t<readfilepath>\t<chr>\t<outputfile_suffix>\t<geneDB_filename>\t<output_folder>" << endl;
		return 1;
	}
	chrName = argv[2];
	char inputfilename[1000], outputfilename[1000], comd[2000];
	
	sprintf(outputfilename, "%sjunction.txt", dirPrefix);
	junctionFile.open(outputfilename);
	sprintf(outputfilename, "%sexon.txt", dirPrefix);
	exonFile.open(outputfilename);
	sprintf(outputfilename, "%sread.txt", dirPrefix);
	readFile.open(outputfilename);
	sprintf(outputfilename, "%sblock.txt", dirPrefix);
	blockFile.open(outputfilename);
	sprintf(outputfilename, "%sstat.txt", dirPrefix);
	statFile.open(outputfilename, fstream::in | fstream::out | fstream::app);
	sprintf(outputfilename, "%sdatastat.txt", dirPrefix);
	datastatFile.open(outputfilename, fstream::in | fstream::out | fstream::app);
	sprintf(outputfilename, "%scrossgenedata.txt", dirPrefix);
	crossGeneFile.open(outputfilename, fstream::in | fstream::out | fstream::app);
	
	sprintf(outputfilename, "%sallJunction.txt", dirPrefix);
	allJunctionFile.open(outputfilename, fstream::in | fstream::out | fstream::app);
	sprintf(outputfilename, "%sallExon.txt", dirPrefix);
	allExonFile.open(outputfilename, fstream::in | fstream::out | fstream::app);

	strcpy(inputfilename, argv[4]);
	parseDB(inputfilename);

	datastatFile << argv[3] << "\t";
	sprintf(inputfilename, "%s%s.txt", argv[1], argv[2]);
	parse(inputfilename);

	junctionFile.close();
	exonFile.close();
	readFile.close();
	blockFile.close();
	statFile.close();
	datastatFile.close();
	crossGeneFile.close();

	allJunctionFile.close();
	allExonFile.close();

	sprintf(comd, "sort -n +1 -3 %sjunction.txt > %sjunction_1.txt", dirPrefix, dirPrefix);
	system(comd);
	sprintf(comd, "sort -n +1 -3 %sexon.txt > %sexon_1.txt", dirPrefix, dirPrefix);
	system(comd);
	sprintf(comd, "sort -n +1 -3 %sread.txt > %sread_1.txt", dirPrefix, dirPrefix);
	system(comd);
	sprintf(comd, "sort -n +1 -3 %sblock.txt > %sblock_1.txt", dirPrefix, dirPrefix);
	system(comd);
	
//	cout << "sort sort sort!!!!!" << endl;
	sprintf(comd, "rm -r %sjunction.txt", dirPrefix);
	system(comd);
	sprintf(comd, "rm -r %sexon.txt", dirPrefix);
	system(comd);
	sprintf(comd, "rm -r %sread.txt", dirPrefix);
	system(comd);
	sprintf(comd, "rm -r %sblock.txt", dirPrefix);
	system(comd);
	sprintf(comd, "rm -r %sallJunction.txt", dirPrefix);
	system(comd);
	sprintf(comd, "rm -r %sallExon.txt", dirPrefix);
	system(comd);
//	sprintf(comd, "rm -r %sstat.txt", dirPrefix);
//	system(comd);
//	sprintf(comd, "rm -r %sdatastat.txt", dirPrefix);
//	system(comd);
	sprintf(comd, "rm -r %scrossgenedata.txt", dirPrefix);
	system(comd);

	string readCount_total;
	ifstream inputDataStat;
	sprintf(inputfilename, "%sDatasetStat.txt", argv[1]);
	inputDataStat.open(inputfilename);
	getline(inputDataStat, readCount_total);
	inputDataStat.close();

	ofstream outputDataStat;
	sprintf(outputfilename, "%sDatasetStat.txt", dirPrefix);
	outputDataStat.open(outputfilename, fstream::app);
	outputDataStat << readCount_total << endl;
	outputDataStat.close();


	return 0;
}


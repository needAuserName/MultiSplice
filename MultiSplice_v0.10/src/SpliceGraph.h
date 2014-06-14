/*    
 *    SpliceGraph.cpp		
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

const long MAX = 2000000000;
const long MAX_CHR_LENGTH = 1000000000;
const long MAX_Gene = 10000;
const long MAX_Block_NUM = 100000;
const int MAX_Trans_NUM = 100;
const int MAX_TransNameLength = 30;
const int MAX_Read_Type = 10;

class Fragment;
class Transcript;

class parentTranscript
{
public:
	Transcript *parent;
	parentTranscript *next;

	parentTranscript();
};

class Exon
{
public:
	long start_exon;
	long end_exon;
	char strand[10];
	long read_count; // total read counts on this exon
	// for Poisson model
	double read_count_Poisson;
	long *coverage_base;
	parentTranscript *parents;
	Exon *next_exon;

	Exon *clone();
	Exon();
	~Exon();
};

class Junction
{
public:
	long start_juntion;
	long end_juntion;
	long support;

	parentTranscript *parents;
	Junction *next_junction;
//	long connected_exonIdx;

	Junction *clone();
	Junction();
	~Junction();
};

class Transcript
{
public:
	char transID[MAX_TransNameLength];
	long start_trans;
	long end_trans;
	long exonNum;
	long junNum;
	double abundance;
	double samplingProb;
	long transLength;
	Exon *exons_tran;
	Junction *junctions_tran;
	Fragment *fragments_tran;

	Transcript *next_transcript;

	Transcript();
	~Transcript();
};

class Gene
{
public:
	char chromosome[20];
	long startPoint; // gene start point
	long endPoint; // gene end point
	long geneLength;
	char strand[10];
	long tranNum;
	long exonNum;
	long junctionsNum;
	double readcount;
	Transcript *transcripts;
	Exon *exons_gene;
	Junction *junctions_gene;

	Gene *clone();
	Gene();
	~Gene();

};

class geneExonBoundary
{
public:
	long position;
	geneExonBoundary *next;

	geneExonBoundary();
};

class Read
{
public:
	char readNm[100];
	int flag; // 0 for '+' and 16 for '-'
	char chromosome[10];
	long start_Read;
	long start_MatePair; // for matepair end
	long distance; // for matepair end
	char CIGAR[500];
	char TAG[50];

	Read();
};

/************************************************************************/
/* For Block Enumerate */
/************************************************************************/
int BlockNum;
class Block
{
public:
	char* transNameQueue[MAX_Trans_NUM + 1];
	double* ReadOrientedProb; // For each transcript, given all the possible reads, weight to adjust the M matrix

	Fragment *fragsInBlock;
	long consumedEffectLength; // effect exonic length which has already been consumed
	long blocklength; // block length (all the exonic length this block covered)
	int junctionCount;
	long blockCount; // read count on this block

	Block* clone(); // clone a new path same as this one
	Block* next; 
	Block();
	~Block();
};

class Fragment
{
public:
	bool type; // 0 for junction; 1 for exon
	long startPoint;
	long endPoint;
	Fragment* next_frag;
	Fragment* clone();
	Fragment();
};

ifstream inputfile_exon, inputfile_junction, inputfile_read, inputfile_block, inputfile_seq;
bool input_exon = true, input_junction = true, input_read = true, input_block = true, input_gene = true;

/************************************************************************/
/* For bias correction */
/************************************************************************/
const int TRAINSETNUM = 100;
const int MAXITER = 5;
const int epsilon = 10e-4;

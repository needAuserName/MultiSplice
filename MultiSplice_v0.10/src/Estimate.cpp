/*    
 *    Estimate.cpp		
 *    Estimate
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

#include "SpliceGraph.h"

double *EstimatedExonCov;
parentTranscript::parentTranscript()
{
	parent = NULL;
	next = NULL;
}

Exon::Exon()
{
	start_exon = MAX;
	end_exon = 0;
	strand[0] = '\0';
	read_count = 0;
	coverage_base = NULL;
	read_count_Poisson = 0;

	parents = NULL;
	next_exon = NULL;
}

Exon::~Exon()
{
	if (coverage_base != NULL)
	{
		delete[] coverage_base;
	}

	parentTranscript *delparent, *curparent;
	curparent = parents;
	while(curparent != NULL)
	{
		delparent = curparent;
		curparent = delparent->next;
		delete delparent;
	}
	parents = NULL;
}

Exon *Exon::clone()
{
	long i;
	Exon *resultexon;
	resultexon = new Exon;

	resultexon->start_exon = start_exon;
	resultexon->end_exon = end_exon;
	strcpy(resultexon->strand, strand);
	resultexon->read_count = read_count;
	resultexon->coverage_base = new long[end_exon-start_exon+1];
	for (i = 0; i < end_exon-start_exon+1; i++)
	{
		resultexon->coverage_base[i] = 0;
	}

	resultexon->read_count_Poisson = read_count_Poisson;

	resultexon->next_exon = NULL;

	return resultexon;
}

Junction::Junction()
{
	start_juntion = MAX;
	end_juntion = 0;
	support = 0;

	parents = NULL;
	next_junction = NULL;
}

Junction::~Junction()
{
	parentTranscript *delparent, *curparent;
	curparent = parents;
	while(curparent != NULL)
	{
		delparent = curparent;
		curparent = delparent->next;
		delete delparent;
	}
	parents = NULL;
}

Junction * Junction::clone()
{
	Junction *resultjun;
	resultjun = new Junction;

	resultjun->start_juntion = start_juntion;
	resultjun->end_juntion = end_juntion;
	resultjun->support = support;

	resultjun->next_junction = NULL;

	return resultjun;
}


Transcript::Transcript()
{
	transID[0] = '\0';
	start_trans = MAX;
	end_trans = 0;
	exonNum = junNum = 0;
	abundance = 0.0;
	samplingProb = 0.0;
	transLength = 0;
	exons_tran = NULL;
	junctions_tran = NULL;
	fragments_tran = NULL;
	next_transcript = NULL;
}

Transcript::~Transcript()
{
	Fragment *delfrag, *curfrag;
	curfrag = fragments_tran;
	while(curfrag != NULL)
	{
		delfrag = curfrag;
		curfrag = delfrag->next_frag;
		delete delfrag;
	}
	fragments_tran = NULL;

	Exon *delexon, *curexon;
	curexon = exons_tran;
	while(curexon != NULL)
	{
		delexon = curexon;
		curexon = delexon->next_exon;
		delete delexon;
	}
	exons_tran = NULL;

	Junction *deljunction, *curjunction;
	curjunction = junctions_tran;
	while(curjunction != NULL)
	{
		deljunction = curjunction;
		curjunction = deljunction->next_junction;
		delete deljunction;
	}
	junctions_tran = NULL;
}

Gene::Gene()
{
	chromosome[0] = '\0';
	startPoint = MAX;
	endPoint = 0;
	strand[0] = '\0';
	tranNum = 0;
	exonNum = 0;
	junctionsNum = 0;

	transcripts = NULL;
	exons_gene = NULL;
	junctions_gene = NULL;
}

Gene::~Gene()
{
	Exon *delexon, *curexon;
	curexon = exons_gene;
	while(curexon != NULL)
	{
		delexon = curexon;
		curexon = delexon->next_exon;
		delete delexon;
	}
	exons_gene = NULL;

	Junction *deljunction, *curjunction;
	curjunction = junctions_gene;
	while(curjunction != NULL)
	{
		deljunction = curjunction;
		curjunction = deljunction->next_junction;
		delete deljunction;
	}
	junctions_gene = NULL;

	Transcript *deltrans, *curtrans;
	curtrans = transcripts;
	while(curtrans != NULL)
	{
		deltrans = curtrans;
		curtrans = deltrans->next_transcript;
		delete deltrans;
	}
	transcripts = NULL;
}

geneExonBoundary::geneExonBoundary()
{
	position = 0;
	next = NULL;
}

Read::Read()
{
	readNm[0] = '\0';
	chromosome[0] = '\0';
	start_Read = MAX;
	CIGAR[0] = '\0';
	TAG[0] = '\0';
}

Block::Block()
{
	int iLoop, jLoop;
	for (iLoop = 0; iLoop <= MAX_Trans_NUM; iLoop++)
	{
		transNameQueue[iLoop] = new char[MAX_TransNameLength];
		for (jLoop = 0; jLoop < MAX_TransNameLength; jLoop++)
		{
			transNameQueue[iLoop][jLoop] = '\0';
		}
	}
	ReadOrientedProb = new double[MAX_Trans_NUM + 1];
	for (iLoop = 0; iLoop <= MAX_Trans_NUM; iLoop++)
	{
		ReadOrientedProb[iLoop] = 0;
	}
	fragsInBlock = NULL;
	consumedEffectLength = 0;
	blocklength = 0;
	junctionCount = 0;
	blockCount = 0;
	next = NULL;
}

Block::~Block()
{
	if (ReadOrientedProb != NULL)
	{
		delete [] ReadOrientedProb;
	}
	for (int iLoop = 0; iLoop <= MAX_Trans_NUM; iLoop++)
	{
		if (transNameQueue[iLoop] != NULL)
		{
			delete [] transNameQueue[iLoop];
		}
	}
	Fragment *delfrag, *curfrag;
	curfrag = fragsInBlock;
	while(curfrag != NULL)
	{
		delfrag = curfrag;
		curfrag = delfrag->next_frag;
		delete delfrag;
	}
	fragsInBlock = NULL;
}

Block* Block::clone()
{
	Block *resultBlock = new Block;

	for (int iLoop = 0; transNameQueue[iLoop][0] != '\0'; iLoop++)
	{
		strcpy(resultBlock->transNameQueue[iLoop], transNameQueue[iLoop]);
		resultBlock->ReadOrientedProb[iLoop] = ReadOrientedProb[iLoop];
	}

	// copy the information of pathJunction
	Fragment *listTail, *curList, *newlist;
	listTail = NULL;
	curList = fragsInBlock;
	while(curList != NULL)
	{
		//newlist = new Fragment;
		newlist = curList->clone();

		if (resultBlock->fragsInBlock == NULL)
		{
			resultBlock->fragsInBlock = newlist;
			listTail = newlist;
		}
		else
		{
			listTail->next_frag = newlist;
			listTail = newlist;
		}

		curList = curList->next_frag;

	}

	resultBlock->consumedEffectLength = consumedEffectLength;
	resultBlock->blocklength = blocklength;
	resultBlock->junctionCount = junctionCount;
	resultBlock->blockCount = blockCount;

	return resultBlock;
}

Fragment::Fragment()
{
	type = -1;
	startPoint = MAX;
	endPoint = 0;
	next_frag = NULL;
}

Fragment* Fragment::clone()
{
	Fragment *resultFragment;

	resultFragment = new Fragment;
	resultFragment->type = type;
	resultFragment->startPoint = startPoint;
	resultFragment->endPoint = endPoint;
	resultFragment->next_frag = NULL;

	return resultFragment;
}


long ExonNum = 0;
long ReadNum = 0;
Read *readList[10000000];

double mergeSort_Larray[10000000];
double mergeSort_Rarray[10000000];
void* mergeSort_LorderedList[10000000];
void* mergeSort_RorderedList[10000000];
double sortKey[10000000]; 
void* sortList[10000000];

void merge(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortKey[p + i - 1];
		mergeSort_LorderedList[i] = sortList[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortKey[q + j];
		mergeSort_RorderedList[j] = sortList[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_CHR_LENGTH * 2;
	mergeSort_Rarray[n2 + 1] = MAX_CHR_LENGTH * 2;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortKey[k] = mergeSort_Larray[i];
			sortList[k] = mergeSort_LorderedList[i];

			i++;
		} 
		else
		{
			sortKey[k] = mergeSort_Rarray[j];
			sortList[k] = mergeSort_RorderedList[j];

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

// Get gene ends from geneID
void getEnd(char *geneName, long &start, long &end)
{
	int i, tmp1, tmp2;
	char startPoint[1000], endPoint[1000];
	i = tmp1 = tmp2 = 0;
	while(geneName[i] != '-')
	{
		startPoint[tmp1] = geneName[i];
		i++;
		tmp1++;
	}
	startPoint[tmp1] = '\0';
	i++;
	while(geneName[i] != '\0')
	{
		endPoint[tmp2] = geneName[i];
		i++;
		tmp2++;
	}
	endPoint[tmp2] = '\0';

	start = atol(startPoint);
	end = atol(endPoint);

}

long geneListNum;
Gene *geneList[MAX_Gene];


// From like "chr2C" get the chromsome name "chr2" and the strand "C"
void GetChrStrand(char *name, char *chromosome, char *strand)
{
	int iLoop, tmp;

	for (iLoop = 0, tmp = 0; name[iLoop] != 'C' && name[iLoop] != 'W'; iLoop++, tmp++)
	{
		chromosome[tmp] = name[iLoop];
	}
	chromosome[tmp] = '\0';
	strand[0] = name[iLoop];
	strand[1] = '\0';

	return;
}


void getTranscriptsFromGTF(char *gtf_path, char *geneName)
	//this is for GTF format!!
{
	char inputfilename[500];
	sprintf(inputfilename, "%s%s", gtf_path, geneName);
	ifstream inputfile;
	inputfile.open(inputfilename);

	bool newTranscript;
	char chrNm[100], chromosome[50], strand[5], lineCategory[100], geneID[300], transID[300], tmpChar[100];
	int tmp;
	long start, end, iLoop;
	string otherInfo;
	Transcript *curTrans, *newTrans;
	Exon *newExon;
	Gene *curGene;
	curGene = new Gene;
	getEnd(geneName, start, end);
	curGene->startPoint = start;
	curGene->endPoint = end;

	for (tmp = 0; tmp < 100; tmp++)
	{
		chrNm[tmp] = '\0';
	}

	inputfile >> chrNm;
	if (chrNm[0] != '\0')
	{
		//		GetChrStrand(chrNm, chromosome, strand);
		strcpy(curGene->chromosome, chrNm);
		/*		strcpy(curGene->strand, strand);*/
	}
	while (chrNm[0] != '\0')
	{
		inputfile >> tmpChar;
		inputfile >> lineCategory;
		inputfile >> start;
		inputfile >> end;
		inputfile >> tmpChar;
		inputfile >> strand;
		inputfile >> tmpChar;
		inputfile >> tmpChar;
		inputfile >> geneID;
		inputfile >> tmpChar;
		inputfile >> transID;
		getline(inputfile, otherInfo);

		if (strcmp(lineCategory, "exon") == 0)
		{
			for (tmp = 0; geneID[tmp+1] != '\"'; tmp++)
			{
				geneID[tmp] = geneID[tmp+1];
			}
			geneID[tmp] = '\0';

			for (tmp = 0; transID[tmp+1] != '\"'; tmp++)
			{
				transID[tmp] = transID[tmp+1];
			}
			transID[tmp] = '\0';

			// insert this transcript into curGene if it is new
			curTrans = curGene->transcripts;
			newTranscript = true;
			if (curTrans == NULL)
			{
				newTranscript = true;
				newTrans = new Transcript;
				strcpy(newTrans->transID, transID);
				curTrans = curGene->transcripts = newTrans;
				curGene->tranNum++;
			}
			else
			{
				while(curTrans != NULL)
				{
					if (strcmp(curTrans->transID, transID) == 0)
					{
						newTranscript = false;
						break;
					}
					curTrans = curTrans->next_transcript;
				}
				if (newTranscript == true)
				{
					// this is a new transcript, insert it into the geneList
					newTrans = new Transcript;
					strcpy(newTrans->transID, transID);
					curTrans = curGene->transcripts;
					while (curTrans->next_transcript != NULL)
					{
						curTrans = curTrans->next_transcript;
					}
					curTrans->next_transcript = newTrans;
					curTrans = newTrans;
					curGene->tranNum++;
				}
			}

			//make a new exon for the transcript
			newExon = new Exon;

			if (strcmp(strand, "+") == 0)
			{
				strcpy(newExon->strand, "+");
			} 
			else if (strcmp(strand, "-") == 0)
			{
				strcpy(newExon->strand, "-");
			}
			else
			{
				cout << "Error: unrecognized strand. Please confirm... ";
				cin >> tmpChar;
				exit(1);
			}

			newExon->start_exon = start;
			newExon->end_exon = end;
			newExon->next_exon = curTrans->exons_tran;
			curTrans->exons_tran = newExon;

			if (start < curTrans->start_trans)
			{
				curTrans->start_trans = start;
			}
			if (end	> curTrans->end_trans)
			{
				curTrans->end_trans = end;
			}
		}
		else
		{
			//do nothing
		}

		chrNm[0] = '\0';
		inputfile >> chrNm;		
	}

	geneListNum++;
	geneList[geneListNum] = curGene;

	inputfile.close();

	return;
}

// Decide whether or not the exon is the start or end exon of one transcript
// 1 indicates the start exon; 2 indicates the end exon;  3 indicates both; and 0 indicates neither
int ExonType(Gene *curGene, Exon *exon)
{
	int type = 0;
	Transcript *curTrans;
	curTrans = curGene->transcripts;
	Exon *curExon;
	while(curTrans != NULL)
	{
		curExon = curTrans->exons_tran;
		while(curExon->next_exon != NULL)
		{
			curExon = curExon->next_exon;
		}
		if (curTrans->exons_tran->start_exon == exon->start_exon && curTrans->exons_tran->end_exon >= exon->end_exon)
		{
			if (type == 2 || type == 3)
			{
				type = 3;
			}
			else
				type = 1;
		}
		else if (curExon->start_exon <= exon->start_exon && curExon->end_exon == exon->end_exon)
		{
			if (type == 1 || type == 3)
			{
				type = 3;
			}
			else
				type = 2;
		}

		if (type == 3)
		{
			break;
		}
		curTrans = curTrans->next_transcript;
	}
	return type;
}

void makeGeneExons(long geneIndex)
{
	long iLoop, jLoop, kLoop, exonCnt, junCnt, splicesite_cnt, *splicesiteList, curExon_start, curExon_end;
	Gene *curGene;
	Transcript *curTrans;
	Exon *curExon, *curGeneExon, *newExon, *deleteExon;
	Junction *newJunc, *curJunc, *curGeneJunc, *deleteJunc;
	bool exonExist, juncExist, splicesiteExist;
	geneExonBoundary *boundary_head, *boundary_2ndlast, *boundary_new; //2nd last: where new ones should be inserted
	parentTranscript *newParent, *curParent;

	curGene = geneList[geneIndex];
	curGene->exons_gene = NULL;
	curGene->junctions_gene = NULL;

	//build junction list
	curTrans = curGene->transcripts;
	while (curTrans != NULL)
	{
		//sort the exon list 
		exonCnt = 0;
		curExon = curTrans->exons_tran;
		while (curExon != NULL)
		{
			exonCnt++;
			sortKey[exonCnt]  = double(curExon->start_exon);
			sortList[exonCnt] = (void*)(curExon);

			curExon = curExon->next_exon;
		}
		mergeSort(exonCnt);

		curTrans->exons_tran = NULL;
		for (jLoop = exonCnt; jLoop >= 1 ; jLoop--)
		{
			((Exon*)sortList[jLoop])->next_exon = curTrans->exons_tran;
			curTrans->exons_tran = (Exon*)sortList[jLoop];
		}
		curTrans->exonNum = exonCnt;

		//build junction list
		curTrans->junctions_tran = NULL;
		for (jLoop = exonCnt; jLoop > 1 ; jLoop--)
		{
			if ((((Exon*)sortList[jLoop])->start_exon - ((Exon*)sortList[jLoop-1])->end_exon) > 1)
			{
				newJunc = new Junction;
				newJunc->start_juntion = ((Exon*)sortList[jLoop-1])->end_exon;
				newJunc->end_juntion   = ((Exon*)sortList[jLoop])->start_exon;
				newJunc->next_junction  = curTrans->junctions_tran;
				curTrans->junctions_tran = newJunc;
				curTrans->junNum += 1;
			}
		}

		curTrans = curTrans->next_transcript;
	}

	//build gene junction list
	curTrans = curGene->transcripts;
	while (curTrans != NULL)
	{
		curJunc = curTrans->junctions_tran;
		while (curJunc != NULL)
		{
			juncExist = false;
			curGeneJunc = curGene->junctions_gene;
			while (curGeneJunc != NULL)
			{
				if (curJunc->start_juntion == curGeneJunc->start_juntion && curJunc->end_juntion == curGeneJunc->end_juntion)
				{
					juncExist = true;

					newParent = new parentTranscript;
					newParent->parent = curTrans;
					newParent->next = curGeneJunc->parents;
					curGeneJunc->parents = newParent;


					break;
				}
				curGeneJunc = curGeneJunc->next_junction;
			}

			if (juncExist == false)
			{
				newJunc = new Junction;
				newJunc->start_juntion = curJunc->start_juntion;
				newJunc->end_juntion   = curJunc->end_juntion;
				newJunc->next_junction  = curGene->junctions_gene;

				newParent = new parentTranscript;
				newParent->parent = curTrans;
				newJunc->parents = newParent;

				curGene->junctions_gene = newJunc;
				curGene->junctionsNum += 1;
			}

			curJunc = curJunc->next_junction;
		}

		curTrans = curTrans->next_transcript;
	}

	//build splice sites list !! ADD exon start and end... vegfa is an example 
	jLoop = 0;
	curGeneJunc = curGene->junctions_gene;
	while (curGeneJunc != NULL)
	{
		splicesiteExist = false;
		for (kLoop = 1; kLoop <= jLoop; kLoop++)
		{
			if (double(curGeneJunc->start_juntion) == sortKey[kLoop])
			{
				splicesiteExist = true;
				break;
			}
		}
		if (splicesiteExist == false)
		{
			jLoop++;
			sortList[jLoop] = NULL;
			sortKey[jLoop]  = double(curGeneJunc->start_juntion);
		}

		splicesiteExist = false;
		for (kLoop = 1; kLoop <= jLoop; kLoop++)
		{
			if (double(curGeneJunc->end_juntion) == sortKey[kLoop])
			{
				splicesiteExist = true;
				break;
			}
		}
		if (splicesiteExist == false)
		{
			jLoop++;
			sortList[jLoop] = NULL;
			sortKey[jLoop]  = double(curGeneJunc->end_juntion);
		}

		curGeneJunc = curGeneJunc->next_junction;
	}

	curTrans = curGene->transcripts;
	while (curTrans != NULL)
	{
		curExon = curTrans->exons_tran;
		while (curExon != NULL)
		{
			splicesiteExist = false;
			for (kLoop = 1; kLoop <= jLoop; kLoop++)
			{
				if (double(curExon->start_exon) == sortKey[kLoop])
				{
					splicesiteExist = true;
					break;
				}
			}
			if (splicesiteExist == false)
			{
				jLoop++;
				sortList[jLoop] = NULL;
				sortKey[jLoop]  = double(curExon->start_exon);
			}

			splicesiteExist = false;
			for (kLoop = 1; kLoop <= jLoop; kLoop++)
			{
				if (double(curExon->end_exon) == sortKey[kLoop])
				{
					splicesiteExist = true;
					break;
				}
			}
			if (splicesiteExist == false)
			{
				jLoop++;
				sortList[jLoop] = NULL;
				sortKey[jLoop]  = double(curExon->end_exon);
			}

			curExon = curExon->next_exon;
		}

		curTrans = curTrans->next_transcript;
	}

	splicesite_cnt = jLoop;
	mergeSort(splicesite_cnt);

	splicesiteList = new long [splicesite_cnt+1];
	for (jLoop = 1; jLoop <= splicesite_cnt; jLoop++)
	{
		splicesiteList[jLoop] = long(sortKey[jLoop]);
	}


	//build gene exon list: cut exons & throw away redundant ones
	curTrans = curGene->transcripts;
	while (curTrans != NULL)
	{
		jLoop = 1;
		curExon = curTrans->exons_tran;
		while (curExon != NULL)
		{
			boundary_head = new geneExonBoundary;
			boundary_head->position = curExon->start_exon;
			boundary_new = new geneExonBoundary;
			boundary_new->position = curExon->end_exon;
			boundary_head->next = boundary_new;
			boundary_2ndlast = boundary_head;

			for (; jLoop <= splicesite_cnt && splicesiteList[jLoop] <= curExon->end_exon; jLoop++)
			{
				if (splicesiteList[jLoop] > curExon->start_exon && splicesiteList[jLoop] < curExon->end_exon)
				{
					boundary_new = new geneExonBoundary;
					boundary_new->position = splicesiteList[jLoop];
					boundary_new->next = boundary_2ndlast->next;
					boundary_2ndlast->next = boundary_new;
					boundary_2ndlast = boundary_new;
				}
			}

			boundary_new = boundary_head;
			while (boundary_new->next != NULL)
			{
				curExon_start = boundary_new->position;
				curExon_end   = boundary_new->next->position;

				exonExist = false;
				curGeneExon = curGene->exons_gene;
				while (curGeneExon != NULL)
				{
					if (curExon_start == curGeneExon->start_exon && curExon_end == curGeneExon->end_exon)
					{
						exonExist = true;

						newParent = new parentTranscript;
						newParent->parent = curTrans;
						newParent->next = curGeneExon->parents;
						curGeneExon->parents = newParent;

						break;
					}
					curGeneExon = curGeneExon->next_exon;
				}

				if (exonExist == false)
				{
					newExon = new Exon;
					newExon->start_exon = curExon_start;
					newExon->end_exon   = curExon_end;
					newExon->next_exon  = curGene->exons_gene;

					newParent = new parentTranscript;
					newParent->parent = curTrans;
					newExon->parents = newParent;

					curGene->exons_gene = newExon;
					curGene->exonNum += 1;
				}

				boundary_new = boundary_new->next;
			}

			curExon = curExon->next_exon;
		}

		curTrans = curTrans->next_transcript;
	}

	//clean up
	delete [] splicesiteList;

	//sort the exon list for this gene
	exonCnt = 0;
	curExon = curGene->exons_gene;
	while (curExon != NULL)
	{
		exonCnt++;
		sortKey[exonCnt]  = double(curExon->start_exon);
		sortList[exonCnt] = (void*)(curExon);
		curExon = curExon->next_exon;
	}
	mergeSort(exonCnt);

	curGene->exons_gene = NULL;
	for (jLoop = exonCnt; jLoop >= 1 ; jLoop--)
	{
		((Exon*)sortList[jLoop])->next_exon = curGene->exons_gene;
		curGene->exons_gene = (Exon*)sortList[jLoop];
	}

	// 	curExon = curGene->exons_gene;
	// 	while(curExon != NULL)
	// 	{
	// 		if (curExon->start_exon == 54360112 && curExon->end_exon == 54360111)
	// 		{
	// 			cout << 111 << endl;
	// 		}
	// 		curExon = curExon->next_exon;
	// 	}

	// First use junction to adjust the adjacent bases of exons
	Exon *prevExon;
	prevExon = curExon = curGene->exons_gene;
	while(curExon->next_exon != NULL)
	{
		if (curExon->end_exon == curExon->next_exon->start_exon)
		{
			curJunc = curGene->junctions_gene;
			while(curJunc != NULL)
			{
				if (curJunc->start_juntion == curExon->end_exon)
				{
					curExon->next_exon->start_exon = curExon->next_exon->start_exon + 1;
					if (curExon->next_exon->start_exon > curExon->next_exon->end_exon)
					{
						// delete curExon->next_exon
						curExon->next_exon = curExon->next_exon->next_exon;
					}
					break;
				}
				else if (curJunc->end_juntion == curExon->end_exon)
				{
					curExon->end_exon = curExon->end_exon - 1;
					if (curExon->end_exon < curExon->start_exon)
					{
						// delete curExon
						prevExon->next_exon = curExon->next_exon;
						curExon = prevExon;
					}
					break;
				}
				curJunc = curJunc->next_junction;
			}

		}
		prevExon = curExon;
		curExon = curExon->next_exon;
	}

	// 	curExon = curGene->exons_gene;
	// 	while(curExon != NULL)
	// 	{
	// 		if (curExon->start_exon == 54360112 && curExon->end_exon == 54360111)
	// 		{
	// 			cout << 111 << endl;
	// 		}
	// 		curExon = curExon->next_exon;
	// 	}

	// Use exons to adjust the adjacent bases of exons
	int type;
	prevExon = curExon = curGene->exons_gene;
	while(curExon->next_exon != NULL)
	{
		if (curExon->end_exon == curExon->next_exon->start_exon)
		{
			type = ExonType(curGene, curExon->next_exon);
			if (type == 1 || type == 3)
			{
				// the start exon of one transcript
				curExon->end_exon = curExon->end_exon - 1;
				if (curExon->end_exon < curExon->start_exon)
				{
					// delete curExon
					prevExon->next_exon = curExon->next_exon;
					curExon = prevExon;
				}
			}
			else
			{
				curExon->next_exon->start_exon = curExon->next_exon->start_exon + 1;
				if (curExon->next_exon->start_exon > curExon->next_exon->end_exon)
				{
					// delete curExon->next_exon
					curExon->next_exon = curExon->next_exon->next_exon;
				}
			}
		}
		prevExon = curExon;
		curExon = curExon->next_exon;
	}

	// 	curExon = curGene->exons_gene;
	// 	while(curExon != NULL)
	// 	{
	// 		if (curExon->start_exon == 54360112 && curExon->end_exon == 54360111)
	// 		{
	// 			cout << 111 << endl;
	// 		}
	// 		curExon = curExon->next_exon;
	// 	}

	// rebuild the exonlist and fragment list for transcripts
	curTrans = curGene->transcripts;
	while(curTrans != NULL)
	{
		// build fragment
		Fragment *newFrag, *curFrag;
		curExon = curTrans->exons_tran;
		while (curExon != NULL)
		{
			newFrag = new Fragment;
			newFrag->type = 1;
			newFrag->startPoint = curExon->start_exon;
			newFrag->endPoint = curExon->end_exon;
			curTrans->transLength += newFrag->endPoint - newFrag->startPoint + 1;
			if(curTrans->fragments_tran == NULL)
			{
				curFrag = curTrans->fragments_tran = newFrag;
			}
			else
			{
				curFrag->next_frag = newFrag;
				curFrag = newFrag;
			}
			if (curExon->next_exon != NULL)
			{
				if ((curExon->next_exon->start_exon - curExon->end_exon) > 1)
				{
					newFrag = new Fragment;
					newFrag->type = 0;
					newFrag->startPoint = curExon->end_exon;
					newFrag->endPoint = curExon->next_exon->start_exon;
					curFrag->next_frag = newFrag;
					curFrag = newFrag;
				}
			}

			curExon = curExon->next_exon;
		}


		curExon = curTrans->exons_tran;
		while (curExon != NULL)
		{
			deleteExon = curExon;
			curExon = deleteExon->next_exon;
			delete deleteExon;

		}
		curTrans->exons_tran = NULL;

		curJunc = curTrans->junctions_tran;
		while (curJunc != NULL)
		{
			deleteJunc = curJunc;
			curJunc = deleteJunc->next_junction;
			delete deleteJunc;

		}
		curTrans->junctions_tran = NULL;		

		// insert exons into transcripts
		curExon = curGene->exons_gene;
		while (curExon != NULL)
		{
			curFrag = curTrans->fragments_tran;
			while(curFrag != NULL)
			{
				if (curFrag->type == 1)
				{
					if (curExon->start_exon >= curFrag->startPoint && curExon->end_exon <= curFrag->endPoint)
					{
						newExon = new Exon;
						newExon->start_exon = curExon->start_exon;
						newExon->end_exon   = curExon->end_exon;
						newExon->next_exon  = curTrans->exons_tran;

						curTrans->exons_tran = newExon;
						break;
					}
				}
				curFrag = curFrag->next_frag;
			}
			// 			curParent = curExon->parents;
			// 			while(curParent != NULL)
			// 			{
			// 				if (strcmp(curTrans->transID, curParent->parent->transID) == 0)
			// 				{
			// 					newExon = new Exon;
			// 					newExon->start_exon = curExon->start_exon;
			// 					newExon->end_exon   = curExon->end_exon;
			// 					newExon->next_exon  = curTrans->exons_tran;
			// 
			// 					curTrans->exons_tran = newExon;
			// 					break;
			// 				}
			// 				curParent = curParent->next;
			// 			}
			curExon = curExon->next_exon;
		}

		//sort the exon list 
		exonCnt = 0;
		curExon = curTrans->exons_tran;
		while (curExon != NULL)
		{
			exonCnt++;
			sortKey[exonCnt]  = double(curExon->start_exon);
			sortList[exonCnt] = (void*)(curExon);

			curExon = curExon->next_exon;
		}
		mergeSort(exonCnt);

		curTrans->exons_tran = NULL;
		for (jLoop = exonCnt; jLoop >= 1 ; jLoop--)
		{
			((Exon*)sortList[jLoop])->next_exon = curTrans->exons_tran;
			curTrans->exons_tran = (Exon*)sortList[jLoop];
		}
		curTrans->exonNum = exonCnt;

		curExon = curTrans->exons_tran;
		while(curExon->next_exon != NULL)
		{
			if ((curExon->next_exon->start_exon - curExon->end_exon) > 1)
			{
				newJunc = new Junction;
				newJunc->start_juntion = curExon->end_exon;
				newJunc->end_juntion   = curExon->next_exon->start_exon;
				newJunc->next_junction  = curTrans->junctions_tran;

				curTrans->junctions_tran = newJunc;
			}
			curExon = curExon->next_exon;
		}

		// sort the junction list 
		junCnt = 0;
		curJunc = curTrans->junctions_tran;
		while (curJunc != NULL)
		{
			junCnt++;
			sortKey[junCnt]  = double(curJunc->start_juntion);
			sortList[junCnt] = (void*)(curJunc);

			curJunc = curJunc->next_junction;
		}
		mergeSort(junCnt);

		curTrans->junctions_tran = NULL;
		for (jLoop = junCnt; jLoop >= 1 ; jLoop--)
		{
			((Junction*)sortList[jLoop])->next_junction = curTrans->junctions_tran;
			curTrans->junctions_tran = (Junction*)sortList[jLoop];
		}
		curTrans->junNum = junCnt;

		curTrans = curTrans->next_transcript;
	}

}


int ReadTypeNum;
long *Readlength; // an array saves all the possible read length;
double *ReadProb; // an array saves the probability for all possible read length;

// Get all the possible read length and probability


// Compare two fragments
bool CompareFragments(Fragment *frag1, Fragment *frag2)
{
	bool Equal = false;
	if (frag1->type == frag2->type && frag1->startPoint == frag2->startPoint && frag1->endPoint == frag2->endPoint)
	{
		Equal = true;
	}
	return Equal;
}

// Compare whether or not two blocks having the same junction sets
bool CompareTransNames(Block *existBlock, Block *newBlock)
{
	bool EqualTransName = true;
	Fragment *existFrags, *newFrags;
	int iLoop, jLoop, junctionCount;

	if (existBlock->junctionCount == newBlock->junctionCount)
	{
		existFrags = existBlock->fragsInBlock;
		newFrags = newBlock->fragsInBlock;
		junctionCount = 0;
		while(existFrags != NULL && newFrags != NULL)
		{
			if (existFrags->type == 0 && newFrags->type == 0)
			{
				junctionCount++;
				if (CompareFragments(existFrags, newFrags) == false)
				{
					EqualTransName = false;
					break;
				}
				existFrags = existFrags->next_frag;
				newFrags = newFrags->next_frag;
			}	
			else if (existFrags->type == 1 && newFrags->type == 0)
			{
				existFrags = existFrags->next_frag;
			}
			else if (existFrags->type == 0 && newFrags->type == 1)
			{
				newFrags = newFrags->next_frag;
			}
			else
			{
				existFrags = existFrags->next_frag;
				newFrags = newFrags->next_frag;
			}

		}
		if (EqualTransName == true && junctionCount == existBlock->junctionCount)
		{
			// Block has already exists, copy the transcripts names from newBlock to existBlock
			bool flag = false;
			for (jLoop = 0; newBlock->transNameQueue[jLoop][0] != '\0'; jLoop++)
			{
				flag = false;
				for (iLoop = 0; existBlock->transNameQueue[iLoop][0] != '\0'; iLoop++)
				{
					if (strcmp(existBlock->transNameQueue[iLoop], newBlock->transNameQueue[jLoop]) == 0)
					{
						flag = true;
						break;
					}
				}
				if (flag == false)
				{
					strcpy(existBlock->transNameQueue[iLoop], newBlock->transNameQueue[jLoop]);
				}
			}

		}
		else
			EqualTransName = false;
	}
	else
	{
		int i, j;
		EqualTransName = false;

		bool flag = true;
		for (j = 0; newBlock->transNameQueue[j][0] != '\0'; j++)
		{
			for (i = 0; existBlock->transNameQueue[i][0] != '\0'; i++)
			{
				if (strcmp(existBlock->transNameQueue[i], newBlock->transNameQueue[j]) == 0)
				{
					break;
				}
			}
			if (existBlock->transNameQueue[i][0] == '\0')
			{
				flag = false;
				break;
			}
		}
		if (flag == true && newBlock->transNameQueue[j][0] == '\0')
		{
			// do nothing
		}
		else
		{
			if (existBlock->junctionCount < newBlock->junctionCount)
			{
				if (existBlock->junctionCount >= 1)
				{
					Block *subBlock;
					Fragment *curFrag, *tempFrag, *listTail, *curList, *newlist;
					long junCount, index, count;
					curFrag = tempFrag = newBlock->fragsInBlock;
					count = 0;
					Fragment *junctionfrag;
					junctionfrag = existBlock->fragsInBlock;
					while(junctionfrag != NULL)
					{
						if (junctionfrag->type == 0)
						{
							break;
						}
						junctionfrag = junctionfrag->next_frag;
					}
					while(curFrag != NULL)
					{
						if (curFrag->type == 0)
						{
							if (curFrag->startPoint != junctionfrag->startPoint)
							{
								// do nothing
							}
							else
							{
								junCount = 0;
								for (junCount = 1; junCount <= newBlock->junctionCount - count; junCount++)
								{
									index = 0;
									tempFrag = curFrag;
									while(tempFrag != NULL)
									{
										if (tempFrag->type == 0)
										{
											index++;
											if (index == junCount)
											{
												break;
											}
										}
										tempFrag = tempFrag->next_frag;
									}
									if (tempFrag->next_frag != NULL)
									{
										if (tempFrag->next_frag->type == true)
										{
											tempFrag = tempFrag->next_frag;
										}
									}

									// create new prefix subBlock and compare to the resultBlocks
									subBlock = new Block;
									listTail = NULL;
									curList = curFrag;
									while(curList != tempFrag->next_frag)
									{
										//							newlist = new Fragment;
										newlist = curList->clone();
										if (curList->type == 0)
										{
											subBlock->junctionCount++;
										}
										if (subBlock->fragsInBlock == NULL)
										{
											subBlock->fragsInBlock = newlist;
											listTail = newlist;
										}
										else
										{
											listTail->next_frag = newlist;
											listTail = newlist;
										}
										curList = curList->next_frag;
									}
									for (i = 0; newBlock->transNameQueue[i][0] != '\0'; i++)
									{
										strcpy(subBlock->transNameQueue[i], newBlock->transNameQueue[i]);
									}
									if (existBlock->junctionCount == subBlock->junctionCount)
									{
										if (CompareTransNames(existBlock, subBlock) == true)
										{
											// do nothing
										}
									}
									delete subBlock;
								}
							}
							count++;
						}
						curFrag = curFrag->next_frag;
					}
				}

			}
			else if (existBlock->junctionCount > newBlock->junctionCount)
			{
				Block *subBlock;
				Fragment *curFrag, *tempFrag, *listTail, *curList, *newlist;
				long junCount, index, count;
				curFrag = tempFrag = existBlock->fragsInBlock;
				count = 0;
				Fragment *junctionfrag;
				junctionfrag = newBlock->fragsInBlock;
				while(junctionfrag != NULL)
				{
					if (junctionfrag->type == 0)
					{
						break;
					}
					junctionfrag = junctionfrag->next_frag;
				}
				while(curFrag != NULL)
				{
					if (curFrag->type == 0)
					{
						if (curFrag->startPoint != junctionfrag->startPoint)
						{
							// do nothing
						}
						else
						{
							junCount = 0;
							for (junCount = 1; junCount <= existBlock->junctionCount - count; junCount++)
							{
								index = 0;
								tempFrag = curFrag;
								while(tempFrag != NULL)
								{
									if (tempFrag->type == 0)
									{
										index++;
										if (index == junCount)
										{
											break;
										}
									}
									tempFrag = tempFrag->next_frag;
								}
								if (tempFrag->next_frag != NULL)
								{
									if (tempFrag->next_frag->type == true)
									{
										tempFrag = tempFrag->next_frag;
									}
								}

								// create new prefix subBlock and compare to the resultBlocks
								subBlock = new Block;
								listTail = NULL;
								curList = curFrag;
								while(curList != tempFrag->next_frag)
								{
									//							newlist = new Fragment;
									newlist = curList->clone();
									if (curList->type == 0)
									{
										subBlock->junctionCount++;
									}
									if (subBlock->fragsInBlock == NULL)
									{
										subBlock->fragsInBlock = newlist;
										listTail = newlist;
									}
									else
									{
										listTail->next_frag = newlist;
										listTail = newlist;
									}
									curList = curList->next_frag;
								}
								for (i = 0; existBlock->transNameQueue[i][0] != '\0'; i++)
								{
									strcpy(subBlock->transNameQueue[i], existBlock->transNameQueue[i]);
								}
								if (newBlock->junctionCount == subBlock->junctionCount)
								{
									if (CompareTransNames(newBlock, subBlock) == true)
									{
										// do nothing
									}
								}
								delete subBlock;
							}
						}
						count++;
					}
					curFrag = curFrag->next_frag;
				}
			}
		}

	}
	return EqualTransName;
}

double ComputeBlockProb(Block *block, Block *parent, long readlength)
{
	double prob;
	long effectiveLength, lastfragLength;
	Fragment *curFrag, *curFrag_temp;
	curFrag = block->fragsInBlock;
	while(curFrag->next_frag != NULL)
	{
		curFrag = curFrag->next_frag;
	}
	lastfragLength = curFrag->endPoint - curFrag->startPoint + 1;

	curFrag = block->fragsInBlock;
	curFrag_temp = parent->fragsInBlock;

	while (curFrag->next_frag != NULL && curFrag_temp->next_frag != NULL)
	{
		curFrag = curFrag->next_frag;
		curFrag_temp = curFrag_temp->next_frag;
	}

	effectiveLength = 0;
	while(curFrag != NULL)
	{
		if (curFrag->type == 1)
		{
			effectiveLength += curFrag->endPoint - curFrag->startPoint + 1;
		}
		curFrag = curFrag->next_frag;
	}
	effectiveLength = effectiveLength - lastfragLength;
	prob = double (readlength - effectiveLength)/readlength;
	return prob;
}

/* Enumerate all the junction blocks given the maximum read length and the minimum read length*/
Block* Enumerate_Junction_Block(long geneIndex, long maxi_length, long mini_length)
{
	Block *resultBlocks, *newBlock, *curBlock, *tempBlock, *parent;
	Gene *curGene;
	Transcript *curTrans;
	Fragment *curFrag;
	curTrans = NULL;
	Fragment *listTail, *curList, *newlist;
	listTail = curList = newlist = NULL;

	resultBlocks = newBlock = curBlock = NULL;
	bool EnqueueSuccess;
	char tmpchar;


	// Find all the blocks
	BlockNum = 0;

	curGene = geneList[geneIndex];
	curTrans = curGene->transcripts;
	while(curTrans != NULL)
	{
		if (curTrans->transLength <= mini_length && curTrans->junNum != 0)
		{
			// this transcript can only generate one block
			newBlock = new Block;
			listTail = NULL;
			curList = curTrans->fragments_tran;
			while(curList != NULL)
			{
				//				newlist = new Fragment;
				newlist = curList->clone();
				if (curList->type == 0)
				{
					newBlock->junctionCount++;
				}
				if (newBlock->fragsInBlock == NULL)
				{
					newBlock->fragsInBlock = newlist;
					listTail = newlist;
				}
				else
				{
					listTail->next_frag = newlist;
					listTail = newlist;
				}

				curList = curList->next_frag;

			}
			newBlock->blocklength = curTrans->transLength;
			newBlock->consumedEffectLength = newBlock->blocklength - (newBlock->fragsInBlock->endPoint - newBlock->fragsInBlock->startPoint);
			strcpy(newBlock->transNameQueue[0], curTrans->transID);
			newBlock->ReadOrientedProb[0] = 1.0;
			BlockNum++;
			if (resultBlocks == NULL)
			{
				resultBlocks = newBlock;
			}
			else
			{
				newBlock->next = resultBlocks;
				resultBlocks = newBlock;
			}
		}
		else
		{
			// this transcript can generate multiple blocks 
			// enumerate blocks from the first exon
			int junctionCount = 0;
			long tempLength = curTrans->transLength;
			long currentLength, currentConsumedLength;
			Fragment *nextJun, *tempFrag, *tempFrag2;
			nextJun = tempFrag = tempFrag2 = NULL;
			curFrag = curTrans->fragments_tran;
			double prob;
			while(curFrag!= NULL && tempLength >= mini_length)
			{
				if (curFrag->type == 1)
				{
					// this fragment is an exon
					junctionCount = 0;
					currentLength = curFrag->endPoint - curFrag->startPoint + 1;
					currentConsumedLength = 1;
					// try to find first junction after this exon
					nextJun = curFrag->next_frag;

					if (nextJun != NULL && nextJun->type == 0)
					{
						// found the first junction after this exon
						tempFrag = nextJun;
						while(tempFrag != NULL && (currentLength < mini_length || junctionCount == 0))
						{
							if (tempFrag->type == 0)
							{
								junctionCount++;
							}
							if (tempFrag->type == 1)
							{
								currentLength += tempFrag->endPoint - tempFrag->startPoint + 1;
								currentConsumedLength += tempFrag->endPoint - tempFrag->startPoint + 1;
							}
							tempFrag = tempFrag->next_frag;
						}

						if (tempFrag != NULL)
						{
							if (tempFrag->type == 1)
							{
								currentLength += tempFrag->endPoint - tempFrag->startPoint + 1;
								currentConsumedLength += tempFrag->endPoint - tempFrag->startPoint + 1;
								tempFrag = tempFrag->next_frag;
							}
						}
						newBlock = new Block;
						strcpy(newBlock->transNameQueue[0], curTrans->transID);
						listTail = NULL;
						curList = curFrag;
						while(curList != tempFrag)
						{
							//								newlist = new Fragment;
							newlist = curList->clone();
							if (curList->type == 0)
							{
								newBlock->junctionCount++;
							}
							if (newBlock->fragsInBlock == NULL)
							{
								newBlock->fragsInBlock = newlist;
								listTail = newlist;
							}
							else
							{
								listTail->next_frag = newlist;
								listTail = newlist;
							}

							curList = curList->next_frag;
						}

						newBlock->blocklength = currentLength;
						newBlock->consumedEffectLength = newBlock->blocklength - (newBlock->fragsInBlock->endPoint - newBlock->fragsInBlock->startPoint);
						newBlock->ReadOrientedProb[0] = 1.0;
						parent = newBlock->clone();
						//							prob = ComputeBlockProb(newBlock, curTrans->transLength);
						// insert this block into resultBlocks

						if (newBlock->junctionCount >= 1)
						{
							bool flag = false;
							if (resultBlocks == NULL)
							{
								BlockNum++;
								resultBlocks = newBlock;
							}
							else
							{
								curBlock = resultBlocks;
								while(curBlock != NULL)
								{
									if (CompareTransNames(curBlock, newBlock) == true)
									{
										flag = true;
									}
									curBlock = curBlock->next;
								}
								if (flag == false)
								{
									BlockNum++;
									newBlock->next = resultBlocks;
									resultBlocks = newBlock;
								}
								else
								{
									delete newBlock;
								}
							}
						}
						else
						{
							delete newBlock;
						}

						if (tempFrag != NULL)
						{
							tempFrag2 = tempFrag->next_frag;
							while(tempFrag2 != NULL && currentConsumedLength < maxi_length)
							{
								if (tempFrag2->type == 1)
								{
									// first find all the exons directly connected to this exon
									if ((tempFrag2->next_frag != NULL && tempFrag2->next_frag->type == 0) || tempFrag2->next_frag == NULL)
									{
										newBlock = new Block;
										strcpy(newBlock->transNameQueue[0], curTrans->transID);
										listTail = NULL;
										curList = curFrag;
										while(curList != tempFrag2->next_frag)
										{
											//												newlist = new Fragment;
											newlist = curList->clone();
											if (curList->type == 0)
											{
												newBlock->junctionCount++;
											}
											if (newBlock->fragsInBlock == NULL)
											{
												newBlock->fragsInBlock = newlist;
												listTail = newlist;
											}
											else
											{
												listTail->next_frag = newlist;
												listTail = newlist;
											}

											curList = curList->next_frag;
										}

										currentConsumedLength += tempFrag2->endPoint - tempFrag2->startPoint + 1;
										newBlock->consumedEffectLength = currentConsumedLength;
										newBlock->blocklength = newBlock->consumedEffectLength + (newBlock->fragsInBlock->endPoint - newBlock->fragsInBlock->startPoint);
										newBlock->ReadOrientedProb[0] = ComputeBlockProb(newBlock, parent, maxi_length);

										if (newBlock->junctionCount >= 1)
										{
											bool flag = false;
											if (resultBlocks == NULL)
											{
												BlockNum++;
												resultBlocks = newBlock;
											}
											else
											{
												curBlock = resultBlocks;
												while(curBlock != NULL)
												{
													if (CompareTransNames(curBlock, newBlock) == true)
													{
														flag = true;
													}
													curBlock = curBlock->next;
												}
												if (flag == false)
												{
													BlockNum++;
													newBlock->next = resultBlocks;
													resultBlocks = newBlock;
												}
												else
												{
													delete newBlock;
												}
											}
										}
										else
											delete newBlock;

										tempFrag2 = tempFrag2->next_frag;
									}
									else
									{
										Fragment *tempFrag3, *prev;
										tempFrag3 = tempFrag2->next_frag;
										while(tempFrag3 != NULL && tempFrag3->type != 0)
										{
											currentConsumedLength += tempFrag3->endPoint - tempFrag3->startPoint + 1;
											prev = tempFrag3;
											tempFrag3 = tempFrag3->next_frag;
										}
										newBlock = new Block;
										strcpy(newBlock->transNameQueue[0], curTrans->transID);
										listTail = NULL;
										curList = curFrag;
										while(curList != prev->next_frag)
										{
											//												newlist = new Fragment;
											newlist = curList->clone();
											if (curList->type == 0)
											{
												newBlock->junctionCount++;
											}
											if (newBlock->fragsInBlock == NULL)
											{
												newBlock->fragsInBlock = newlist;
												listTail = newlist;
											}
											else
											{
												listTail->next_frag = newlist;
												listTail = newlist;
											}

											curList = curList->next_frag;
										}

										currentConsumedLength += tempFrag2->endPoint - tempFrag2->startPoint + 1;
										newBlock->consumedEffectLength = currentConsumedLength;
										newBlock->blocklength = newBlock->consumedEffectLength + (newBlock->fragsInBlock->endPoint - newBlock->fragsInBlock->startPoint); 
										newBlock->ReadOrientedProb[0] = ComputeBlockProb(newBlock, parent, maxi_length);	

										if (newBlock->junctionCount >= 1)
										{
											bool flag = false;
											if (resultBlocks == NULL)
											{
												BlockNum++;
												resultBlocks = newBlock;
											}
											else
											{
												curBlock = resultBlocks;
												while(curBlock != NULL)
												{
													if (CompareTransNames(curBlock, newBlock) == true)
													{
														flag = true;

													}
													curBlock = curBlock->next;
												}
												if (flag == false)
												{
													BlockNum++;
													newBlock->next = resultBlocks;
													resultBlocks = newBlock;
												}
												else
												{
													delete newBlock;
												}
											}
										}
										else 
											delete newBlock;
										tempFrag2 = prev->next_frag;
									}

								}
								else
								{
									tempFrag2 = tempFrag2->next_frag;
								}

							}
						}
						delete parent;

					}
					tempLength = tempLength - (curFrag->endPoint - curFrag->startPoint + 1);

				}// end if (curFrag is an exon)

				curFrag = curFrag->next_frag;
			}

		}
		curTrans = curTrans->next_transcript;
	}

// 	ofstream outputfile;
// 	outputfile.open("./test_block.txt");
// 	curBlock = resultBlocks;
// 	while(curBlock != NULL)
// 	{
// 		outputfile << curBlock->blocklength << endl;
// 		for (int iLoop = 0; curBlock->transNameQueue[iLoop][0] != '\0'; iLoop++)
// 		{
// 			outputfile << curBlock->transNameQueue[iLoop] << '\t' << curBlock->ReadOrientedProb[iLoop] << endl;
// 		}
// 		curFrag = curBlock->fragsInBlock;
// 		while(curFrag != NULL)
// 		{
// 			if (curFrag->type == 1)
// 			{
// 				outputfile << curFrag->startPoint << "-" << curFrag->endPoint << '\t';
// 			}
// 
// 			curFrag = curFrag->next_frag;
// 		}
// 		outputfile << endl;
// 		curBlock = curBlock->next;
// 	}
// 	outputfile.close();

	return resultBlocks;

}

// Compute block count from simulation
void ComputeBlockCount(Block *resultBlocks, Block *newBlock)
{
	Block *curBlock;
	if (newBlock->junctionCount == 1)
	{
		curBlock = resultBlocks;
		while(curBlock != NULL)
		{
			if (CompareTransNames(curBlock, newBlock) == true)
			{
				curBlock->blockCount++;
				break;
			}
			curBlock = curBlock->next;
		}
	}
	else
	{
		// build subBlocks
		Block *subBlock;
		Fragment *curFrag, *tempFrag, *listTail, *curList, *newlist;
		long junctionCount, index, count;
		curFrag = tempFrag = newBlock->fragsInBlock;
		count = 0;

		while(curFrag != NULL)
		{
			if (curFrag->type == 0)
			{
				junctionCount = 0;
				for (junctionCount = 1; junctionCount <= newBlock->junctionCount - count; junctionCount++)
				{
					index = 0;
					tempFrag = curFrag;
					while(tempFrag != NULL)
					{
						if (tempFrag->type == 0)
						{
							index++;
							if (index == junctionCount)
							{
								break;
							}
						}
						tempFrag = tempFrag->next_frag;
					}
					if (tempFrag->next_frag != NULL)
					{
						if (tempFrag->next_frag->type == true)
						{
							tempFrag = tempFrag->next_frag;
						}
					}

					// create new prefix subBlock and compare to the resultBlocks
					subBlock = new Block;
					listTail = NULL;
					curList = curFrag;
					while(curList != tempFrag->next_frag)
					{
						//						newlist = new Fragment;
						newlist = curList->clone();
						if (curList->type == 0)
						{
							subBlock->junctionCount++;
						}
						if (subBlock->fragsInBlock == NULL)
						{
							subBlock->fragsInBlock = newlist;
							listTail = newlist;
						}
						else
						{
							listTail->next_frag = newlist;
							listTail = newlist;
						}
						curList = curList->next_frag;
					}
					curBlock = resultBlocks;
					while(curBlock != NULL)
					{
						if (CompareTransNames(curBlock, subBlock) == true)
						{
							curBlock->blockCount++;
							break;
						}
						curBlock = curBlock->next;
					}
					delete subBlock;
				}
				count++;
			}
			curFrag = curFrag->next_frag;
		}

	}

	return;
}



double *PoissonEstimation;
void Poisson_Estimate(long geneIndex)
{
	long row, col, iLoop, jLoop, kLoop, tmp1, tmp2;
	char tmpChar;
	Exon * currExon, *exonToComp;
	Transcript *currenttrans;
	long **TransFrag;
	currExon = exonToComp = NULL;
	currenttrans = NULL;
	TransFrag = NULL;

	row = geneList[geneIndex]->exonNum;
	col = geneList[geneIndex]->tranNum;

	PoissonEstimation = new double[col];
	// Calculate transcript-exon matrix
	TransFrag = new long *[col];
	for (iLoop = 0; iLoop < col; iLoop++)
	{
		TransFrag[iLoop] = new long[row];
		for (jLoop = 0; jLoop < row; jLoop++)
		{
			TransFrag[iLoop][jLoop] = 0;
		}
	}

	currExon = geneList[geneIndex]->exons_gene;
	tmp1 = tmp2 = 0;
	// tmp1 keeps the index of exon
	// tmp2 keeps the index of transcript
	while(currExon != NULL)
	{
		currenttrans = geneList[geneIndex]->transcripts;
		tmp2 = 0;
		while(currenttrans != NULL)
		{
			exonToComp = currenttrans->exons_tran;
			while (exonToComp != NULL)
			{
				if ((currExon->start_exon == exonToComp->start_exon) && (currExon->end_exon == exonToComp->end_exon))
				{
					// This transcript covers this exon
					TransFrag[tmp2][tmp1] = 1;
					break;
				}

				exonToComp = exonToComp->next_exon;
			}

			tmp2++;
			currenttrans = currenttrans->next_transcript;
		}
		tmp1++;
		currExon = currExon->next_exon;
	}

	//estimation
	double *P = new double [row]; // probability that reads fall into each exon
	double *Y = new double [row]; // observed exon expression
	double *Q = new double [col]; //proportions of alternative transcripts
	double *Qold = new double [col]; //new Q vector
	double *Qtmp = new double [col];
	double *Qdelta = new double [col]; //change on Q
	double *Z = new double [col]; //expectation of read count
	double lambda; 

	//calculate P
	long lengthSum = 0;
	currExon = geneList[geneIndex]->exons_gene;
	while(currExon != NULL)
	{
		lengthSum += currExon->end_exon - currExon->start_exon + 1;
		currExon = currExon->next_exon;
	}

	currExon = geneList[geneIndex]->exons_gene;
	iLoop = 0;
	while(currExon != NULL)
	{
		P[iLoop] = (double)(currExon->end_exon - currExon->start_exon + 1)/lengthSum;
		iLoop++;
		currExon = currExon->next_exon;
	}


	//fill in Y
	currExon = geneList[geneIndex]->exons_gene;
	iLoop = 0;
	while(currExon != NULL)
	{
		Y[iLoop] = currExon->read_count_Poisson;
		iLoop++;
		currExon = currExon->next_exon;
	}

	//initiate Q, Z, transExonCnt}
	for (iLoop = 0; iLoop < col; iLoop++)
	{
		Q[iLoop] = 1/col;
		Qtmp[iLoop] = 0.0;
		Qold[iLoop] = 0.0;
		Qdelta[iLoop] = 0.0;

		Z[iLoop] = 0.0;
	}

	//EM
	double sum_count, sum_denominator, sum_numerator;
	bool stoppable = false, tmpFlag;
	long MAX_ITER = 100, iterCnt = 0;

	while (iterCnt < MAX_ITER && stoppable == false)
	{
		for (iLoop = 0; iLoop < col; iLoop++)
		{
			Qold[iLoop] = Q[iLoop];
		}

		//E-step
		for (iLoop = 0; iLoop < col; iLoop++)
		{
			sum_count = 0.0;

			for (jLoop = 0; jLoop < row; jLoop++)
			{
				if (TransFrag[iLoop][jLoop] == 1)
				{
					sum_denominator = 0.0;
					for (kLoop = 0; kLoop < col; kLoop++)
					{
						if (TransFrag[kLoop][jLoop] == 1)
						{
							sum_denominator += P[jLoop] * Q[kLoop];
						}
					}

					if (sum_denominator > 0.0)
					{
						sum_count += P[jLoop] * Q[iLoop] * Y[jLoop] / sum_denominator;
					}
				}
			}

			Z[iLoop] = sum_count;
		}

		//M-step
		//calculate lambda
		lambda = 0.0;
		for (iLoop = 0; iLoop < col; iLoop++)
		{
			lambda += Z[iLoop];
		}

		//calculate Q
		//get Q iteratively until stable
		if (col > 1)
		{
			tmpFlag = false;
			while (tmpFlag == false)
			{
				for (iLoop = 0; iLoop < col; iLoop++)
				{
					Qtmp[iLoop] = Q[iLoop];
				}
				for (iLoop = 0; iLoop < col; iLoop++)
				{
					sum_denominator = 0.0;
					for (jLoop = 0; jLoop < row; jLoop++)
					{
						if (TransFrag[iLoop][jLoop] == 1)
						{
							sum_denominator += P[jLoop];
						}
					}
					sum_denominator *= (lambda - Z[iLoop]);

					sum_numerator = 0.0;
					for (jLoop = 0; jLoop < col; jLoop++)
					{
						if (jLoop == iLoop)
						{
							continue;
						}

						sum_count = 0.0;
						for (kLoop = 0; kLoop < row; kLoop++)
						{
							if (TransFrag[jLoop][kLoop] == 1)
							{
								sum_count += P[kLoop];
							}
						}
						sum_count *= Qtmp[jLoop];

						sum_numerator += sum_count;
					}
					sum_numerator *= Z[iLoop];

					if (sum_denominator > 0.0)
					{
						Qtmp[iLoop] = sum_numerator / sum_denominator;
					} 
					else
					{
						if (lambda == Z[iLoop])
						{
							Qtmp[iLoop] = 1.0;
						}
						else
						{
							Qtmp[iLoop] = 0.0;
						}
					}
				}

				sum_count = 0.0;
				for (iLoop = 0; iLoop < col; iLoop++)
				{
					sum_count += Qtmp[iLoop];
				}
				for (iLoop = 0; iLoop < col; iLoop++)
				{
					//calculate change on Q
					if (sum_count > 0.0)
					{
						Qtmp[iLoop] = Qtmp[iLoop] / sum_count;
					} 
					else
					{
						Qtmp[iLoop] = 1.0 / col;
					}
					Qdelta[iLoop] = fabs(Qtmp[iLoop] - Q[iLoop]);
					Q[iLoop] = Qtmp[iLoop];
				}

				tmpFlag = true;
				for (iLoop = 0; iLoop < col; iLoop++)
				{
					if (Q[iLoop] > 0.0 && fabs(Qdelta[iLoop] / Q[iLoop]) > 0.001)
					{
						tmpFlag = false;
						break;
					}
				}
			}

			for (iLoop = 0; iLoop < col; iLoop++)
			{
				//calculate change on Q
				Qdelta[iLoop] = fabs(Qold[iLoop] - Q[iLoop]);
			}

			iterCnt++;
			tmpFlag = true;
			for (iLoop = 0; iLoop < col; iLoop++)
			{
				if (Q[iLoop] > 0.0 && fabs(Qdelta[iLoop] / Q[iLoop]) > 0.001)
				{
					tmpFlag = false;
					break;
				}
			}
			if (iterCnt >= 3 && tmpFlag == true)
			{
				stoppable = true;
			}	
		} 
		else
		{
			stoppable = true;
		}

	}

	//one more E-step
	for (iLoop = 0; iLoop < col; iLoop++)
	{
		sum_count = 0.0;

		for (jLoop = 0; jLoop < row; jLoop++)
		{
			if (TransFrag[iLoop][jLoop] == 1)
			{
				sum_denominator = 0.0;
				for (kLoop = 0; kLoop < col; kLoop++)
				{
					if (TransFrag[kLoop][jLoop] == true)
					{
						sum_denominator += P[jLoop] * Q[kLoop];
					}
				}

				if (sum_denominator > 0.0)
				{
					sum_count += P[jLoop] * Q[iLoop] * Y[jLoop] / sum_denominator;
				}
			}
		}

		Z[iLoop] = sum_count;
	}

	// Copy the proportion result
	for (iLoop = 0; iLoop < col; iLoop++)
	{
		//		PoissonEstimation[iLoop] = Q[iLoop];
		PoissonEstimation[iLoop] = Z[iLoop];
	}

	return;
}

/************************************************************************/
/* Process the exon, read and block files */
/* Note that both GTF and SAM format are 1-based; BED format are 0-based */
/************************************************************************/
// information on a line of the inputfile
char chromosome_exon[10];
long start_exon;
long end_exon;
string s_exon;

// process the exon file to get the read coverage for the exon.
void process_exons()
{
	long j;

	Exon *currentexon;

	if (input_exon == true)
	{
		chromosome_exon[0] = '\0';
		inputfile_exon >> chromosome_exon;
	}
	while(chromosome_exon[0] != '\0')
	{
		if (input_exon == true)
		{
			inputfile_exon >> start_exon;
			inputfile_exon >> end_exon;
			end_exon = end_exon - 1;
			getline(inputfile_exon, s_exon);
		}

		if ((start_exon >= geneList[1]->startPoint) && (end_exon <= geneList[1]->endPoint))
		{
			// this read is contained in this gene
			input_exon = true;

			// count the coverage per base on exons
			currentexon = geneList[1]->exons_gene;
			while(currentexon != NULL)
			{
				if (start_exon >= currentexon->start_exon && end_exon <= currentexon->end_exon)
				{
					// this read is completely falling on this exon 
					for (j = start_exon; j <= end_exon; j++)
					{
						currentexon->coverage_base[j - currentexon->start_exon]++;
						currentexon->read_count++;
					}
					break;
				}
				else if ((start_exon >= currentexon->start_exon && start_exon <= currentexon->end_exon) && end_exon > currentexon->end_exon)
				{
					// the downstream part of the read is in this exon
					for(j = start_exon; j <= currentexon->end_exon; j++)
					{
						currentexon->coverage_base[j - currentexon->start_exon]++;
						currentexon->read_count++;
					}
				}
				else if (start_exon < currentexon->start_exon && (end_exon >= currentexon->start_exon && end_exon <= currentexon->end_exon))
				{
					// the upstream part of the read is in this exon
					for (j = currentexon->start_exon; j <= end_exon; j++)
					{
						currentexon->coverage_base[j - currentexon->start_exon]++;
						currentexon->read_count++;
					}
					break;
				}	
				else if (start_exon < currentexon->start_exon && end_exon > currentexon->end_exon)
				{
					// this exon is completely covered within this read
					for (j = currentexon->start_exon; j <= currentexon->end_exon; j++)
					{
						currentexon->coverage_base[j - currentexon->start_exon]++;
						currentexon->read_count++;
					}
				}
				currentexon = currentexon->next_exon;
			}
		}
		else if ((start_exon < geneList[1]->startPoint && end_exon > geneList[1]->startPoint) || (start_exon >= geneList[1]->startPoint && start_exon < geneList[1]->endPoint))
		{
			/*			cout << "This kind of read sampling can never happen!!!(parts are the gene and parts not)" << endl;*/
			input_exon = true;

			// count the coverage per base of the part reads falling on exons
			currentexon = geneList[1]->exons_gene;
			while(currentexon != NULL)
			{

				if ((start_exon >= currentexon->start_exon && start_exon <= currentexon->end_exon) && end_exon > currentexon->end_exon)
				{
					// the downstream part of the read is in this exon
					for(j = start_exon; j <= currentexon->end_exon; j++)
					{
						currentexon->coverage_base[j - currentexon->start_exon]++;
						currentexon->read_count++;
					}
				}
				else if (start_exon < currentexon->start_exon && (end_exon >= currentexon->start_exon && end_exon <= currentexon->end_exon))
				{
					// the upstream part of the read is in this exon
					for (j = currentexon->start_exon; j <= end_exon; j++)
					{
						currentexon->coverage_base[j - currentexon->start_exon]++;
						currentexon->read_count++;
					}
					break;
				}	
				else if (start_exon < currentexon->start_exon && end_exon > currentexon->end_exon)
				{
					// this exon is completely covered within this read
					for (j = currentexon->start_exon; j <= currentexon->end_exon; j++)
					{
						currentexon->coverage_base[j - currentexon->start_exon]++;
						currentexon->read_count++;
					}
				}
				currentexon = currentexon->next_exon;
			}
		}
		else if (end_exon <= geneList[1]->startPoint)
		{
			// this read is sampled downstream of the gene
			input_exon = true;

		}
		else if (start_exon >= geneList[1]->endPoint)
		{
			// this read is sampled upstream of the gene
			input_exon = false;
			break;
		}		
		if (input_exon == true)
		{
			chromosome_exon[0] = '\0';
			inputfile_exon >> chromosome_exon;
		}

	}
	return;
}

// information on a line of the inputfile
char chromosome_junc[10];
long start_junc;
long end_junc;
string s_junc;

// process the exon file to get the read coverage for the exon.
void process_junctions()
{
	long j;

	Junction *curJunc;

	if (input_junction == true)
	{
		chromosome_junc[0] = '\0';
		inputfile_junction >> chromosome_junc;
	}
	while(chromosome_junc[0] != '\0')
	{
		if (input_junction == true)
		{
			inputfile_junction >> start_junc;
			inputfile_junction >> end_junc;
			getline(inputfile_junction, s_junc);
		}

		if ((start_junc >= geneList[1]->startPoint) && (end_junc <= geneList[1]->endPoint))
		{
			// this read is contained in this gene
			input_junction = true;

			// count the coverage per base on exons
			curJunc = geneList[1]->junctions_gene;
			while(curJunc != NULL)
			{
				if (curJunc->start_juntion == start_junc && curJunc->end_juntion == end_junc)
				{
					curJunc->support++;
					break;
				}
				curJunc = curJunc->next_junction;
			}
		}
		else if ((start_junc < geneList[1]->startPoint && (end_junc > geneList[1]->startPoint && end_junc <= geneList[1]->endPoint)) || ((start_junc >= geneList[1]->startPoint && end_junc < geneList[1]->endPoint) && end_junc > geneList[1]->endPoint))
		{
			input_junction = true;
		}
		else if (end_junc <= geneList[1]->startPoint)
		{
			// this read is sampled downstream of the gene
			input_junction = true;

		}
		else if (start_junc >= geneList[1]->endPoint)
		{
			// this read is sampled upstream of the gene
			input_junction = false;
			break;
		}	
		else
		{
			input_junction = true;
		}
		if (input_junction == true)
		{
			chromosome_junc[0] = '\0';
			inputfile_junction >> chromosome_junc;
		}

	}
	return;
}

// information on a line of the inputfile
char chromosome[10];
long start;
long end;
string s;
// recount the read count on exon based on Poisson model
void process_read()
{
	long j;

	Exon *currentexon;

	if (input_read == true)
	{
		chromosome[0] = '\0';
		inputfile_read >> chromosome;
	}
	while(chromosome[0] != '\0')
	{
		if (input_read == true)
		{
			inputfile_read >> start;
			inputfile_read >> start;
			getline(inputfile_read, s);
		}

		if ((start >= geneList[1]->startPoint) && (start < geneList[1]->endPoint))
		{
			// this read is contained in this gene
			input_read = true;

			// count the read count starting from that base on a specific
			currentexon = geneList[1]->exons_gene;
			while(currentexon != NULL)
			{
				if (start >= currentexon->start_exon && start <= currentexon->end_exon)
				{
					currentexon->read_count_Poisson++;
					break;
				}

				currentexon = currentexon->next_exon;
			}
		}
		else if (start < geneList[1]->startPoint)
		{
			// this read is sampled downstream of the gene
			input_read = true;

		}
		else if (start >= geneList[1]->endPoint)
		{
			// this read is sampled upstream of the gene
			input_read = false;
			break;
		}
		if (input_read == true)
		{
			chromosome[0] = '\0';
			inputfile_read >> chromosome;
		}

	}

	return;
}


Block *block;
// information on a line of the inputfile
char buffer[100];

// process the block file and get the support on each multi-junction block
void process_block(Block *resultBlocks)
{
	long j, start, end, point;
	bool pair = false, flag = false;
	Exon *currentexon;
	Fragment *newFrag, *FragCursor, *curFrag;
	newFrag = FragCursor = curFrag = NULL;

	if (input_block == true)
	{
		buffer[0] = '\0';
		inputfile_block >> buffer;
	}
	while(buffer[0] != '\0')
	{
		if (input_block == true)
		{
			if (buffer[0] == 'c')
			{
				// a new line
				buffer[0] = '\0';
				inputfile_block >> buffer;
				start = atol(buffer);
				if (start <= geneList[1]->startPoint)
				{
					flag = true;

					while(buffer[0] >= '0' && buffer[0] <= '9')
					{
						buffer[0] = '\0';
						inputfile_block >> buffer;
					}
				}
				else
				{
					flag = false;

					block = new Block;
					pair = false;
					while(buffer[0] >= '0' && buffer[0] <= '9')
					{
						point = atol(buffer);
						if (pair == false)
						{
							start = point;
							pair = true;
						}
						else
						{
							// find the start and end point of the junction
							end = point;
							pair = false;
							block->junctionCount++;
							newFrag = new Fragment;
							newFrag->type = 0;
							newFrag->startPoint = start;
							newFrag->endPoint = end;
							if (block->fragsInBlock == NULL)
							{
								FragCursor = block->fragsInBlock = newFrag;
							}
							else
							{
								FragCursor->next_frag = newFrag;
								FragCursor = newFrag;
							}
						}
						buffer[0] = '\0';
						inputfile_block >> buffer;
					}
				}
			}
		}
		if (flag == true)
		{
			input_block = true;
			if (block != NULL)
			{
				delete block;
				block = NULL;
			}
		}
		else
		{
			start = block->fragsInBlock->startPoint;
			curFrag = block->fragsInBlock;
			while(curFrag->next_frag != NULL)
			{
				curFrag = curFrag->next_frag;
			}
			end = curFrag->endPoint;

			if ((start > geneList[1]->startPoint) && (end < geneList[1]->endPoint))
			{
				// this multi-junction block is contained in this gene
				input_block = true;

				// add support to the corresponding block in all the blocks
				ComputeBlockCount(resultBlocks, block);
				delete block;
				block = NULL;
			}
			else if (start >= geneList[1]->endPoint)
			{
				// this read is sampled upstream of the gene
				input_block = false;
				break;
			}
			else
			{
				input_block = true;
				delete block;
				block = NULL;
			}
		}
	}

	return;
}



/************************************************************************/
/* Bias along the transcripts */
/************************************************************************/
void InitiaTransAbundance(double **structure, double *coverage, long transnum, long segnum, char *outputfile_path, double *abundance, Gene *gene, int type, char *genenames)
{
	char inputfilename[500], outputfilename_struc[500], outputfilename_cov[500], comd[1000];
	sprintf(outputfilename_struc, "%sstructure.txt", outputfile_path);
	sprintf(outputfilename_cov, "%scoverage.txt", outputfile_path);
	ofstream outputfile_structure, outputfile_coverage;

	outputfile_structure.open(outputfilename_struc);
	outputfile_coverage.open(outputfilename_cov);
	long transIndex, segmentIndex;
	for (segmentIndex = 0; segmentIndex < segnum; segmentIndex++)
	{
		for (transIndex = 0; transIndex < transnum; transIndex++)
		{
			outputfile_structure << structure[segmentIndex][transIndex] << '\t';
		}
		outputfile_structure << endl;
		outputfile_coverage << coverage[segmentIndex] << endl;
	}

	outputfile_structure.close();
	outputfile_coverage.close();

	sprintf(comd,  "mv %sstructure.txt %stemp/structure_%s.txt", outputfile_path, outputfile_path, genenames);
	system(comd);
	sprintf(comd,  "mv %scoverage.txt %stemp/coverage_%s.txt", outputfile_path, outputfile_path, genenames);
	system(comd);

// 	sprintf(comd,  "matlab -r \"Estimate('%s', %ld, %d);exit;\"", outputfile_path, gene->exonNum, type);
// 	system(comd);
// 
// 	sprintf(inputfilename, "%sresult.txt", outputfile_path);
// 	ifstream inputfile;
// 	inputfile.open(inputfilename);
// 
// 	// read abundances of transcripts
// 	transIndex = 0;
// 	double number = -1;
// 	string info;
// 	inputfile >> number;
// 	while(number != -1)
// 	{
// 		getline(inputfile, info);
// 		abundance[transIndex] = number;
// 		transIndex++;
// 		number = -1;
// 		inputfile >> number;
// 	}
// 	inputfile.close();

	// read estimate coverage of the exons within this gene
// 	if (segnum > gene->exonNum)
// 	{
		// 		sprintf(inputfilename, "%stemp.txt", outputfile_path);
		// 		inputfile.open(inputfilename);
		// 		transIndex = 0;
		// 		number = -1;
		// 		inputfile >> number;
		// 		while(number != -1)
		// 		{
		// 			getline(inputfile, info);
		// 			EstimatedExonCov[transIndex] = number;
		// 			transIndex++;
		// 			if (transIndex >= gene->exonNum)
		// 			{
		// 				break;
		// 			}
		// 			number = -1;
		// 			inputfile >> number;
		// 		}
		// 
		// 		inputfile.close();

// 		sprintf(comd,  "cp %sstructure.txt %stemp/", outputfile_path, outputfile_path);
// 		system(comd);
// 		sprintf(comd,  "mv %stemp/structure.txt %stemp/structure_%ld-%ld%s.txt", outputfile_path, outputfile_path, gene->startPoint, gene->endPoint, gene->strand);
// 		system(comd);
// 		sprintf(comd,  "cp %scoverage.txt %stemp/", outputfile_path, outputfile_path);
// 		system(comd);
// 		sprintf(comd,  "mv %stemp/coverage.txt %stemp/coverage_%ld-%ld%s.txt", outputfile_path, outputfile_path, gene->startPoint, gene->endPoint, gene->strand);
// 		system(comd);
// 		sprintf(comd,  "cp %sresult.txt %stemp/%s/", outputfile_path, outputfile_path, gene->chromosome);
// 		system(comd);
// 		sprintf(comd,  "mv %stemp/%s/result.txt %stemp/%s/result_%ld-%ld%s.txt", outputfile_path, gene->chromosome, outputfile_path, gene->chromosome, gene->startPoint, gene->endPoint, gene->strand);
// 		system(comd);

//	}
	return;
}

/************************************************************************/
/*  Get sequence-specific bias */
/************************************************************************/
// Get gene ends from geneID
void getEnd_Strand(char *geneName, long &start, long &end)
{
	int i, tmp1, tmp2;
	char startPoint[1000], endPoint[1000];
	i = tmp1 = tmp2 = 0;
	while(geneName[i] != '-')
	{
		startPoint[tmp1] = geneName[i];
		i++;
		tmp1++;
	}
	startPoint[tmp1] = '\0';
	i++;
	while(geneName[i] != 'W' && geneName[i] != 'C')
	{
		endPoint[tmp2] = geneName[i];
		i++;
		tmp2++;
	}
	endPoint[tmp2] = '\0';

	start = atol(startPoint);
	end = atol(endPoint);

}

// Get the sequence bias from learned data;

void process_seqbias(char *inputdata_seqbias, double **sequencebias, char **transName)
{
	char category[20], Id[30];
	long number;
	long transNum, transLength, transIndex;
	string info;
	long start, end;
	char strand[10];
	double prob;
	bool flag;

	ifstream inputfile;
	inputfile.open(inputdata_seqbias);

	category[0] = '\0';
	inputfile >> category;
	while(category[0] != '\0')
	{
		inputfile >> Id;			
		inputfile >> number;

		if (strcmp(category, "Gene") == 0)
		{
			// This line provides gene information
			getEnd_Strand(Id, start, end);
			if (start >= geneList[1]->startPoint && end <= geneList[1]->endPoint)
			{
				flag = true;
			}
			else
			{
				flag = false;
			}

			getline(inputfile, info);
		}
		else if (strcmp(category, "Transcript") == 0)
		{
			// This line provides the probability for sampling bias of each transcript
			if (flag == true)
			{
				transLength = number;
				for (transIndex = 0; transIndex < geneList[1]->tranNum; transIndex++)
				{
					if (strcmp(transName[transIndex], Id) == 0)
					{
						for (long iLoop = 1; iLoop <= transLength; iLoop++)
						{
							inputfile >> prob;
							sequencebias[transIndex][iLoop] = prob;
						}
					}
				}
			}

			getline(inputfile, info);
		}

		category[0] = '\0';
		inputfile >> category;

	}

	inputfile.close();
	return;
}


void FindMapLength(Transcript *trans, Block *block, long &start, long &end, int readlength)
{
	Fragment *curFrag, *FragTocomp;

	// Find first junction in the block
	curFrag = block->fragsInBlock;
	while(curFrag->type != false)
	{
		curFrag = curFrag->next_frag;
	}
	// 1. Find the exonic length of this transcript before the first junction in the block
	FragTocomp = trans->fragments_tran;
	long tmplength = 0;
	while(FragTocomp != NULL)
	{
		if (curFrag->startPoint == FragTocomp->startPoint)
		{
			break;
		}
		else
		{
			if (FragTocomp->type == 1)
			{
				tmplength += FragTocomp->endPoint - FragTocomp->startPoint + 1;
			}
		}
		FragTocomp = FragTocomp->next_frag;
	}

	// 2. Find the exonic length of this transcript after the spliced junctions in this block
	// curFrag now points to the first junction, find the middle exons' length
	int intervallength = 0;
	int juncount = 1;
	while(juncount <= block->junctionCount)
	{
		if (curFrag->type == 1)
		{
			intervallength += curFrag->endPoint - curFrag->startPoint + 1;
		}
		else if (curFrag->type == 0)
		{
			juncount++;
			if (juncount > block->junctionCount)
			{
				break;
			}
		}
		curFrag = curFrag->next_frag;
	}
	// curFrag points to the last junction of this block
	long downstreamlength = 0; // 3'end
	bool flag = false;
	FragTocomp = trans->fragments_tran;
	while(FragTocomp != NULL)
	{
		if (FragTocomp->startPoint == curFrag->startPoint && FragTocomp->endPoint == curFrag->endPoint)
		{
			flag = true;
		}
		if (flag == true && FragTocomp->type == 1)
		{
			downstreamlength += FragTocomp->endPoint - FragTocomp->startPoint + 1;
		}
		FragTocomp = FragTocomp->next_frag;
	}
	start = (tmplength - (readlength-intervallength - 1) + 1) > 1? (tmplength - (readlength-intervallength - 1) + 1):1;
	end = (readlength - downstreamlength - intervallength) > 0? (tmplength - (readlength - downstreamlength - intervallength) + 1):tmplength;

	return;
}


void Bias(double *phi, double *xi, int readlength, Gene *gene, Block *blocks, int blockcount, double **sequencebias, double ***Segbias, double **structure, int fragmentlength)
{
	double **biasweight;
	Transcript *trans = gene->transcripts;
	long transIndex = 0;
	while(trans != NULL)
	{
		biasweight = new double *[2];
		biasweight[0] = new double[trans->transLength + 1];
		biasweight[1] = new double[trans->transLength + 1];
		if (trans->transLength > fragmentlength)
		{
			// Get the bias weight per base
			long index = 1;
			if (trans->transLength >= 2 * fragmentlength)
			{
				while (index <= fragmentlength)
				{
					biasweight[0][index] = (double)fragmentlength*index/fragmentlength * sequencebias[transIndex][index];
					biasweight[1][index] = (double)index/fragmentlength * sequencebias[transIndex][index];

					index++;
				}
				while(index <= trans->transLength - fragmentlength)
				{
					biasweight[0][index] = (double)index * sequencebias[transIndex][index];
					biasweight[1][index] = (double)sequencebias[transIndex][index];

					index++;
				}
				while(index <= trans->transLength)
				{
					biasweight[0][index] = (double)(trans->transLength-fragmentlength)*(trans->transLength-index+1)/fragmentlength * sequencebias[transIndex][index];
					biasweight[1][index] = (double)(trans->transLength-index+1)/fragmentlength * sequencebias[transIndex][index];

					index++;
				}
			}
			else
			{
				while (index <= trans->transLength - fragmentlength + 1)
				{
					biasweight[0][index] = (double)index * (trans->transLength - fragmentlength + 1)/fragmentlength * sequencebias[transIndex][index];
					biasweight[1][index] = (double)index/(trans->transLength - fragmentlength + 1) * (trans->transLength - fragmentlength + 1)/fragmentlength * sequencebias[transIndex][index];

					index++;
				}
				while(index <= fragmentlength)
				{
					biasweight[0][index] = (double)index * (trans->transLength - fragmentlength + 1)/fragmentlength * sequencebias[transIndex][index];
					biasweight[1][index] = (double)(trans->transLength - fragmentlength + 1)/fragmentlength * sequencebias[transIndex][index];

					index++;
				}
				while(index <= trans->transLength)
				{
					biasweight[0][index] = (double)(trans->transLength - fragmentlength + 1)/fragmentlength * fragmentlength*(trans->transLength-index+1)/(trans->transLength - fragmentlength + 1) * sequencebias[transIndex][index];
					biasweight[1][index] = (double)(trans->transLength - fragmentlength + 1)/fragmentlength * 1/(trans->transLength - fragmentlength + 1) *(trans->transLength-index+1) * sequencebias[transIndex][index];

					index++;
				}
			}
			double total0, total1;
			Exon *prevExon = NULL, *curExon = gene->exons_gene, *exonTocomp;
			long tmplength = 0, effectivelength = 0; 
			long segmentIndex = 0;
			while(curExon != NULL)
			{
				exonTocomp = trans->exons_tran;
				while(exonTocomp != NULL)
				{
					if (curExon->start_exon == exonTocomp->start_exon && curExon->end_exon == exonTocomp->end_exon)
					{
						if (prevExon == NULL || (prevExon != NULL && prevExon->end_exon != curExon->start_exon))
						{
							effectivelength += curExon->end_exon - curExon->start_exon + 1;
						}
						else
						{
							effectivelength += curExon->end_exon - curExon->start_exon;
						}

						total0 = total1 = 0;
						for (index = tmplength+1; index <= effectivelength; index++)
						{
							total0 += biasweight[0][index];
							total1 += biasweight[1][index];
						}
						Segbias[0][transIndex][segmentIndex] = total0/(effectivelength - tmplength);
						Segbias[1][transIndex][segmentIndex] = total1/(effectivelength - tmplength);

						break;
					}
					exonTocomp = exonTocomp->next_exon;
				}
				segmentIndex++;
				tmplength = effectivelength;
				prevExon = curExon;
				curExon = curExon->next_exon;
			}

			// Get the bias weight per base
			index = 1;
			while (index <= trans->transLength)
			{
				biasweight[0][index] = (double)index * sequencebias[transIndex][index];
				biasweight[1][index] = (double)sequencebias[transIndex][index];

				index++;
			}

			Block *curBlock;
			curBlock = blocks;
			long start, end;
			while(curBlock != NULL)
			{
				for (int iLoop = 0; curBlock->transNameQueue[iLoop][0] != '\0'; iLoop++)
				{
					if (strcmp(trans->transID, curBlock->transNameQueue[iLoop]) == 0)
					{
						// This block is within this transcript, find the mappable length
						FindMapLength(trans, curBlock, start, end, readlength);
						total0 = total1 = 0;
						for (index = start; index <= end; index++)
						{
							total0 += biasweight[0][index];
							total1 += biasweight[1][index];
						}
						Segbias[0][transIndex][segmentIndex] = total0/readlength;
						Segbias[1][transIndex][segmentIndex] = total1/readlength;
						// 						Segbias[0][transIndex][segmentIndex] = (double)(end - start + 1)/readlength;
						// 						Segbias[1][transIndex][segmentIndex] = (double)(end - start + 1)/readlength;

						break;
					}
				}
				segmentIndex++;
				curBlock = curBlock->next;
			}
		}
		else
		{

			long index = 1;
			while(index <= trans->transLength)
			{
				biasweight[0][index] = sequencebias[transIndex][index];
				biasweight[1][index] = sequencebias[transIndex][index];

				index++;
			}

			double total0, total1;
			Exon *prevExon = NULL, *curExon = gene->exons_gene, *exonTocomp;
			long tmplength = 0, effectivelength = 0; 
			long segmentIndex = 0;
			while(curExon != NULL)
			{
				exonTocomp = trans->exons_tran;
				while(exonTocomp != NULL)
				{
					if (curExon->start_exon == exonTocomp->start_exon && curExon->end_exon == exonTocomp->end_exon)
					{
						if (prevExon == NULL || (prevExon != NULL && prevExon->end_exon != curExon->start_exon))
						{
							effectivelength += curExon->end_exon - curExon->start_exon + 1;
						}
						else
						{
							effectivelength += curExon->end_exon - curExon->start_exon;
						}

						// calculate the accumulative bias on each exon

						total0 = total1 = 0;
						for (index = tmplength+1; index <= effectivelength; index++)
						{
							total0 += biasweight[0][index];
							total1 += biasweight[1][index];
						}
						Segbias[0][transIndex][segmentIndex] = total0/(effectivelength - tmplength);
						Segbias[1][transIndex][segmentIndex] = total1/(effectivelength - tmplength);
						break;
					}
					exonTocomp = exonTocomp->next_exon;
				}
				segmentIndex++;
				tmplength = effectivelength;
				prevExon = curExon;
				curExon = curExon->next_exon;
			}
			Block *curBlock;
			curBlock = blocks;
			long start, end;
			while(curBlock != NULL)
			{
				for (int iLoop = 0; curBlock->transNameQueue[iLoop][0] != '\0'; iLoop++)
				{
					if (strcmp(trans->transID, curBlock->transNameQueue[iLoop]) == 0)
					{
						// This block is within this transcript, find the mappable length
						FindMapLength(trans, curBlock, start, end, readlength);
						total0 = total1 = 0;
						for (index = start; index <= end; index++)
						{
							total0 += biasweight[0][index];
							total1 += biasweight[1][index];
						}
						Segbias[0][transIndex][segmentIndex] = total0/readlength;
						Segbias[1][transIndex][segmentIndex] = total1/readlength;
						// 						Segbias[0][transIndex][segmentIndex] = (double)(end - start + 1)/readlength;
						// 						Segbias[1][transIndex][segmentIndex] = (double)(end - start + 1)/readlength;
						break;
					}
				}
				segmentIndex++;
				curBlock = curBlock->next;
			}
		}
		delete [] biasweight[0];
		delete [] biasweight[1];
		delete [] biasweight;
		transIndex++;
		trans = trans->next_transcript;
	}
	return;
}


double objectiveVal(long segNum, long transNum, double ***Segbias, double *phi, double *xi, double **structure, double *coverage, double *abundance, bool &flag)
{
	double G = 0, temp;
	long segmentIndex, transIndex;
	for (segmentIndex = 0; segmentIndex < segNum; segmentIndex++)
	{
		temp = 0;
		for (transIndex = 0; transIndex < transNum; transIndex++)
		{
			if ((Segbias[0][transIndex][segmentIndex]*phi[transIndex]+Segbias[1][transIndex][segmentIndex]*xi[transIndex])*structure[segmentIndex][transIndex] < 0)
			{
				flag = false;
			}
			temp += (Segbias[0][transIndex][segmentIndex]*phi[transIndex]+Segbias[1][transIndex][segmentIndex]*xi[transIndex])*structure[segmentIndex][transIndex]*abundance[transIndex];
		}
		G += pow((coverage[segmentIndex] - temp), 2);
	}

	return G;
}

void Maximization(int geneIndex, Block *blocks, int readlength, char *outputfile_path, double **sequencebias, double *&abundance, double *&exonabundance, int fragmentlength, char *genenames)
{
	Gene *curGene = geneList[geneIndex];
	Transcript *curTrans;

	// Initial abundance estimation
	abundance = new double[curGene->tranNum];
	exonabundance = new double[curGene->tranNum];
	// find the number of blocks within this gene
	Block *curBlock = blocks;
	long blockcount = 0;
	while(curBlock != NULL)
	{
		blockcount++;
		curBlock = curBlock->next;
	}
	// Gene structure: matrix M &  Coverage on each segment

	long segmentIndex, transIndex;

	// parameter initialize
	double *phi, *xi;
	phi = new double[curGene->tranNum]; // slope
	xi = new double[curGene->tranNum]; // intercept
	for (transIndex = 0; transIndex < curGene->tranNum; transIndex++)
	{
		phi[transIndex] = 0;
		xi[transIndex] = 1;
	}

	double ***Segbias = NULL;
	// Get the bias weight per segment
	Segbias = new double**[2];
	Segbias[0] = new double *[curGene->tranNum];
	Segbias[1] = new double *[curGene->tranNum];
	transIndex = 0, segmentIndex = 0;
	for (transIndex = 0; transIndex < curGene->tranNum; transIndex++)
	{
		Segbias[0][transIndex] = new double[curGene->exonNum+blockcount];
		Segbias[1][transIndex] = new double[curGene->exonNum+blockcount];
		for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
		{
			Segbias[0][transIndex][segmentIndex] = 0;
			Segbias[1][transIndex][segmentIndex] = 0;
		}
	}

	double **structure = NULL, **structureTemp = NULL, **exonstructure = NULL;
	structure = new double *[curGene->exonNum+blockcount];
	structureTemp = new double *[curGene->exonNum+blockcount];
	exonstructure = new double *[curGene->exonNum];
	double *coverage = NULL, *exoncoverage = NULL;
	coverage = new double[curGene->exonNum+blockcount];
	exoncoverage = new double[curGene->exonNum];
	for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
	{
		structure[segmentIndex] = new double[curGene->tranNum];
		structureTemp[segmentIndex] = new double[curGene->tranNum];
		if (segmentIndex < curGene->exonNum)
		{
			exonstructure[segmentIndex] = new double[curGene->tranNum];
		}
		for (transIndex = 0; transIndex < curGene->tranNum; transIndex++)
		{
			structure[segmentIndex][transIndex] = 0;
			if (segmentIndex < curGene->exonNum)
			{
				exonstructure[segmentIndex][transIndex] = 0;
			}
		}
		coverage[segmentIndex] = 0;
		if (segmentIndex < curGene->exonNum)
		{
			coverage[segmentIndex] = 0;
		}
	}
	segmentIndex = 0;
	Exon *curExon = curGene->exons_gene, *exonTocomp;
	while(curExon != NULL)
	{
		curTrans = curGene->transcripts;
		transIndex = 0;
		while(curTrans != NULL)
		{
			exonTocomp = curTrans->exons_tran;
			while(exonTocomp != NULL)
			{
				if (curExon->start_exon == exonTocomp->start_exon && curExon->end_exon == exonTocomp->end_exon)
				{
					structure[segmentIndex][transIndex] = 1;
					exonstructure[segmentIndex][transIndex] = 1;
					break;
				}
				exonTocomp = exonTocomp->next_exon;
			}
			curTrans = curTrans->next_transcript;
			transIndex++;
		}
		coverage[segmentIndex] = (double)curExon->read_count/(curExon->end_exon - curExon->start_exon + 1);
		exoncoverage[segmentIndex] = (double)curExon->read_count/(curExon->end_exon - curExon->start_exon + 1);
		segmentIndex++;
		curExon = curExon->next_exon;
	}

	curBlock = blocks;
	while(curBlock != NULL)
	{
		curTrans = curGene->transcripts;
		transIndex = 0;
		while(curTrans != NULL)
		{
			for (int iLoop = 0; curBlock->transNameQueue[iLoop][0] != '\0'; iLoop++)
			{
				if (strcmp(curTrans->transID, curBlock->transNameQueue[iLoop]) == 0)
				{
					structure[segmentIndex][transIndex] = 1;
					break;
				}
			}
			transIndex++;
			curTrans = curTrans->next_transcript;
		}
		coverage[segmentIndex] = curBlock->blockCount;
		segmentIndex++;
		curBlock = curBlock->next;
	}
	int type = 0;
	//	InitiaTransAbundance(exonstructure, exoncoverage, curGene->tranNum, curGene->exonNum, outputfile_path, exonabundance, curGene);

	// Get the initial value of the objective function
	Bias(phi, xi, readlength, curGene, blocks, blockcount, sequencebias, Segbias, structure, fragmentlength);
	for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
	{
		for (transIndex = 0; transIndex < curGene->tranNum; transIndex++)
		{
			structureTemp[segmentIndex][transIndex] = structure[segmentIndex][transIndex] * (Segbias[0][transIndex][segmentIndex]*phi[transIndex] + Segbias[1][transIndex][segmentIndex]*xi[transIndex]); 
		}
	}
	type = 1;
	InitiaTransAbundance(structureTemp, coverage, curGene->tranNum, curGene->exonNum+blockcount, outputfile_path, abundance, curGene, type, genenames);

// 	double G, G_new, G_old, G_iter, temp;
// 	bool flag = true;
// 	G = objectiveVal(curGene->exonNum+blockcount, curGene->tranNum, Segbias, phi, xi, structure, coverage, abundance, flag);
// 
// 	double norminator, denominator;
// 	for (int iteration = 0; iteration < MAXITER; iteration++)
// 	{
// 		G_iter = G;
// 		curTrans = curGene->transcripts;
// 		for (transIndex = 0; transIndex < curGene->tranNum; transIndex++)
// 		{
// 			if (curTrans->transLength > 2*fragmentlength)
// 			{
// 				bool temflag = false; int inneriter = 0;
// 				while(temflag == false)
// 				{
// 					G_old = G;
// 					// Maximize abundance
// 					norminator = denominator = 0;
// 					for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
// 					{
// 						norminator += coverage[segmentIndex]*(Segbias[0][transIndex][segmentIndex]*phi[transIndex]+Segbias[1][transIndex][segmentIndex]*xi[transIndex])*structure[segmentIndex][transIndex];
// 						denominator += (Segbias[0][transIndex][segmentIndex]*phi[transIndex]+Segbias[1][transIndex][segmentIndex]*xi[transIndex])*structure[segmentIndex][transIndex]*(Segbias[0][transIndex][segmentIndex]*phi[transIndex]+Segbias[1][transIndex][segmentIndex]*xi[transIndex])*structure[segmentIndex][transIndex];
// 					}
// 					//					norminator -= lambda/2;
// 					long transIndextemp;
// 					for (transIndextemp = 0; transIndextemp < curGene->tranNum; transIndextemp++)
// 					{
// 						if (transIndextemp != transIndex)
// 						{
// 							for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
// 							{
// 								norminator -= abundance[transIndextemp]*(Segbias[0][transIndextemp][segmentIndex]*phi[transIndextemp]+Segbias[1][transIndextemp][segmentIndex]*xi[transIndextemp])*structure[segmentIndex][transIndextemp]*(Segbias[0][transIndex][segmentIndex]*phi[transIndex]+Segbias[1][transIndex][segmentIndex]*xi[transIndex])*structure[segmentIndex][transIndex];
// 							}
// 						}
// 					}
// 					temp = abundance[transIndex];
// 					abundance[transIndex] = norminator/denominator;
// 					if (abundance[transIndex] < 0)
// 					{
// 						abundance[transIndex] = 0;
// 					}
// 					flag = true;
// 					G_new = objectiveVal(curGene->exonNum+blockcount, curGene->tranNum, Segbias, phi, xi, structure, coverage, abundance, flag);
// 					if (G_new <= G && flag == true)
// 					{
// 						G = G_new;
// 					}
// 					else
// 					{
// 						abundance[transIndex] = temp;
// 					}
// 
// 					// Maximize phi
// 					if (abundance[transIndex] > 0)
// 					{
// 						norminator = denominator = 0;
// 						for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
// 						{
// 							norminator += coverage[segmentIndex]*Segbias[0][transIndex][segmentIndex]*structure[segmentIndex][transIndex];
// 							norminator -= abundance[transIndex]*Segbias[1][transIndex][segmentIndex]*xi[transIndex]*structure[segmentIndex][transIndex]*Segbias[0][transIndex][segmentIndex]*structure[segmentIndex][transIndex];
// 							denominator += abundance[transIndex]*Segbias[0][transIndex][segmentIndex]*structure[segmentIndex][transIndex]*Segbias[0][transIndex][segmentIndex]*structure[segmentIndex][transIndex];
// 						}
// 						for (transIndextemp = 0; transIndextemp < curGene->tranNum; transIndextemp++)
// 						{
// 							if (transIndextemp != transIndex)
// 							{
// 								for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
// 								{
// 									norminator -= abundance[transIndextemp]*(Segbias[0][transIndextemp][segmentIndex]*phi[transIndextemp]+Segbias[1][transIndextemp][segmentIndex]*xi[transIndextemp])*structure[segmentIndex][transIndextemp]*Segbias[0][transIndex][segmentIndex]*structure[segmentIndex][transIndex];
// 								}
// 							}
// 						}
// 						temp = phi[transIndex];
// 						phi[transIndex] = norminator/denominator;
// 						Bias(phi, xi, readlength, curGene, blocks, blockcount, sequencebias, Segbias, structure, fragmentlength);
// 						flag = true;
// 						G_new = objectiveVal(curGene->exonNum+blockcount, curGene->tranNum, Segbias, phi, xi, structure, coverage, abundance, flag);
// 						if (G_new <= G && flag == true)
// 						{
// 							G = G_new;
// 						}
// 						else
// 						{
// 							phi[transIndex] = temp;
// 						}
// 					}
// 
// 					// Maximize xi
// 					if (abundance[transIndex] > 0)
// 					{
// 						norminator = denominator = 0;
// 						for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
// 						{
// 							norminator += coverage[segmentIndex]*Segbias[1][transIndex][segmentIndex]*structure[segmentIndex][transIndex];
// 							norminator -= abundance[transIndex]*Segbias[0][transIndex][segmentIndex]*phi[transIndex]*structure[segmentIndex][transIndex]*Segbias[1][transIndex][segmentIndex]*structure[segmentIndex][transIndex];
// 							denominator += abundance[transIndex]*Segbias[1][transIndex][segmentIndex]*structure[segmentIndex][transIndex]*Segbias[1][transIndex][segmentIndex]*structure[segmentIndex][transIndex];
// 						}
// 						for (transIndextemp = 0; transIndextemp < curGene->tranNum && transIndextemp != transIndex; transIndextemp++)
// 						{
// 							if (transIndextemp != transIndex)
// 							{
// 								for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
// 								{
// 									norminator -= abundance[transIndextemp]*(Segbias[0][transIndextemp][segmentIndex]*phi[transIndextemp]+Segbias[1][transIndextemp][segmentIndex]*xi[transIndextemp])*structure[segmentIndex][transIndextemp]*Segbias[1][transIndex][segmentIndex]*structure[segmentIndex][transIndex];
// 								}
// 							}
// 						}
// 						temp = xi[transIndex];
// 						xi[transIndex] = norminator/denominator;
// 						if (xi[transIndex] < 0)
// 						{
// 							xi[transIndex] = 0;
// 						}
// 						Bias(phi, xi, readlength, curGene, blocks, blockcount, sequencebias, Segbias, structure, fragmentlength);
// 						flag = true;
// 						G_new = objectiveVal(curGene->exonNum+blockcount, curGene->tranNum, Segbias, phi, xi, structure, coverage, abundance, flag);
// 						if (G_new <= G && flag == true)
// 						{
// 							G = G_new;
// 						}
// 						else
// 						{
// 							xi[transIndex] = temp;
// 						}
// 					}	
// 					Bias(phi, xi, readlength, curGene, blocks, blockcount, sequencebias, Segbias, structure, fragmentlength);
// 
// 					temflag = true;
// 					inneriter++;
// 					if ((G_old - G)/G > epsilon && inneriter < MAXITER*2)
// 					{
// 						temflag = false;
// 					}
// 					else
// 					{
// 						temflag = true;
// 					}
// 				}
// 
// 			}
// 			curTrans = curTrans->next_transcript;
// 
// 		}
// 		if ((G_iter - G)/G < epsilon)
// 		{
// 			break;
// 		}
// 	}
// 	double maxi = 0;
// 	for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
// 	{
// 		for (transIndex = 0; transIndex < curGene->tranNum; transIndex++)
// 		{
// 			structureTemp[segmentIndex][transIndex] = structure[segmentIndex][transIndex] * (Segbias[0][transIndex][segmentIndex]*phi[transIndex] + Segbias[1][transIndex][segmentIndex]*xi[transIndex]); 
// 			if (structureTemp[segmentIndex][transIndex] > maxi)
// 			{
// 				maxi = structureTemp[segmentIndex][transIndex];
// 			}
// 		}
// 	}
// 	for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
// 	{
// 		for (transIndex = 0; transIndex < curGene->tranNum; transIndex++)
// 		{
// 			if (maxi != 0)
// 			{
// 				structureTemp[segmentIndex][transIndex] = structureTemp[segmentIndex][transIndex]/maxi;
// 			}
// 		}
// 	}
// 	InitiaTransAbundance(structureTemp, coverage, curGene->tranNum, curGene->exonNum+blockcount, outputfile_path, abundance, curGene, type);


	for (segmentIndex = 0; segmentIndex < curGene->exonNum+blockcount; segmentIndex++)
	{
		delete [] structure[segmentIndex];
		delete [] structureTemp[segmentIndex];
		if (segmentIndex < curGene->exonNum)
		{
			delete [] exonstructure[segmentIndex];
		}
	}
	delete [] structure;
	delete [] structureTemp;
	delete [] exonstructure;

	delete [] coverage;
	delete [] exoncoverage;
	for (transIndex = 0; transIndex < curGene->tranNum; transIndex++)
	{
		delete [] Segbias[0][transIndex];
		delete [] Segbias[1][transIndex];
	}
	delete [] Segbias[0];
	delete [] Segbias[1];
	delete [] Segbias;

	delete [] phi;
	delete [] xi;

	return;
}

int main(int argc, char* argv[])
{
	if (argc != 10)
	{
		cout << argv[0] << "\t<inputfile_gtf_path>" << "\t<inputfile_data_path>" << "\t<outputfile_path>" << "\t<MaxiReadLength>" << "\t<MiniReadLength>" << "\tchromosome" << "\tSequnceBiasPath" << "\tfragmentlength" << "\treadtype" << endl;
		return 1;
	}


	Block *resultBlocks;
	long iLoop, index;
	double *abundance = NULL, total, *exonabundance = NULL;
	char inputfilename[1000], outputfilename[1000], geneNames[500], inputdata_exon[1000], inputdata_junction[1000], inputdata_read[1000], inputdata_block[1000], inputdata_seqbias[1000], outputfilename_coverage[1000];
	string info;
	double **sequencebias = NULL;


	ifstream inputfile;
	ofstream outputfile_estimation, outputfile_cov;
	sprintf(outputfilename, "%s%s.pos", argv[3], argv[6]);
	outputfile_estimation.open (outputfilename, fstream::in | fstream::out | fstream::app);
	sprintf(inputdata_exon, "%sexon_%s.txt", argv[2], argv[9]);
	inputfile_exon.open(inputdata_exon);
	sprintf(inputdata_read, "%sread_%s.txt", argv[2], argv[9]);
	inputfile_read.open(inputdata_read);
	sprintf(inputdata_block, "%sblock_%s.txt", argv[2], argv[9]);
	inputfile_block.open(inputdata_block);
	sprintf(inputdata_junction, "%sjunction_%s.txt", argv[2], argv[9]);
	inputfile_junction.open(inputdata_junction);

	geneListNum = 0;
	input_exon = true, input_junction = true, input_read = true, input_block = true, input_gene = true;
	block = NULL;

	sprintf(inputfilename, "%sGeneList.txt", argv[1]);
	inputfile.open(inputfilename);

	geneNames[0] = '\0';
	inputfile >> geneNames;
	while(geneNames[0] != '\0')
	{
		getline(inputfile, info);

		getTranscriptsFromGTF(argv[1], geneNames);
		cout << geneList[1]->chromosome << '\t' << geneNames << endl;
		makeGeneExons(1);
		resultBlocks = NULL;
		if (geneList[1]->exonNum > 1)
		{
			resultBlocks = Enumerate_Junction_Block(1, atol(argv[4]), atol(argv[5]));
		}
		Exon *curExon = geneList[1]->exons_gene;
		while(curExon != NULL)
		{
			curExon->coverage_base = new long[curExon->end_exon - curExon->start_exon + 1];
			for (long index = 0; index < curExon->end_exon - curExon->start_exon + 1; index++)
			{
				curExon->coverage_base[index] = 0;
			}
			curExon = curExon->next_exon;
		}
		process_exons();
		if (resultBlocks != NULL)
		{
			process_block(resultBlocks);
		}
		sprintf(inputdata_seqbias, "%sseqbias.txt", argv[7]);

		long index = 0;
		char **transName;
		sequencebias = new double *[geneList[1]->tranNum];
		transName = new char *[geneList[1]->tranNum];
		Transcript *trans = geneList[1]->transcripts;
		for (index = 0; index < geneList[1]->tranNum; index++)
		{
			sequencebias[index] = new double[trans->transLength + 1];
			transName[index] = new char[100];
			strcpy(transName[index], trans->transID);
			for (long iLoop = 0; iLoop <= trans->transLength; iLoop++)
			{
				sequencebias[index][iLoop] = 1;
			}
			trans = trans->next_transcript;
		}
		Maximization(1, resultBlocks, atol(argv[4]), argv[3], sequencebias, abundance, exonabundance, atol(argv[8]), geneNames);

		Transcript *curTrans = geneList[1]->transcripts;
		long transIndex = 0;
		while(curTrans != NULL)
		{
			outputfile_estimation << geneList[1]->chromosome << "\t" << geneNames << '\t' << curTrans->transID << '\t' << curTrans->transLength << endl;
			delete [] sequencebias[transIndex];
			delete [] transName[transIndex];
			transIndex++;
			curTrans = curTrans->next_transcript;
		}
		delete [] sequencebias;
		delete [] transName;

		delete [] abundance;
		abundance = NULL;
		delete [] exonabundance;
		exonabundance = NULL;

		Block *delblock;
		while(resultBlocks)
		{
			delblock = resultBlocks;
			resultBlocks = resultBlocks->next;
			delete delblock;
		}

		delete geneList[1];
		geneList[1] = NULL;
		geneListNum = 0;

		geneNames[0] = '\0';
		inputfile >> geneNames;

	}

	return 0;
}


// int main(int argc, char* argv[])
// {
// 
// 	Block *resultBlocks;
// 	long iLoop, index;
// 	double *abundance = NULL, total, *exonabundance = NULL;
// 	char inputfilename[1000], outputfilename[1000], geneNames[500], inputdata_exon[1000], inputdata_junction[1000], inputdata_read[1000], inputdata_block[1000], inputdata_seqbias[1000], outputfilename_coverage[1000];
// 	string info;
// 	double **sequencebias = NULL;
// 
// 
// 	ifstream inputfile;
// 	ofstream outputfile_estimation, outputfile_cov;
// 	sprintf(outputfilename, "chrX.pos");
// 	outputfile_estimation.open (outputfilename, fstream::in | fstream::out | fstream::app);
// 	sprintf(inputdata_exon, "exon_1.txt");
// 	inputfile_exon.open(inputdata_exon);
// 	sprintf(inputdata_read, "read_1.txt");
// 	inputfile_read.open(inputdata_read);
// 	sprintf(inputdata_block, "block_1.txt");
// 	inputfile_block.open(inputdata_block);
// 	sprintf(inputdata_junction, "junction_1.txt");
// 	inputfile_junction.open(inputdata_junction);
// 
// 	geneListNum = 0;
// 	input_exon = true, input_junction = true, input_read = true, input_block = true, input_gene = true;
// 	block = NULL;
// 
// 	sprintf(inputfilename, "GeneList.txt");
// 	inputfile.open(inputfilename);
// 
// 	geneNames[0] = '\0';
// 	inputfile >> geneNames;
// 	while(geneNames[0] != '\0')
// 	{
// 		getline(inputfile, info);
// 
// 		getTranscriptsFromGTF("", geneNames);
// 		cout << geneList[1]->chromosome << '\t' << geneNames << endl;
// 		cout << 111 << endl;
// 		makeGeneExons(1);
// 		resultBlocks = NULL;
// 		cout << 222 << endl;
// 		if (geneList[1]->exonNum > 1)
// 		{
// 			resultBlocks = Enumerate_Junction_Block(1, 190, 190);
// 		}
// 		Exon *curExon = geneList[1]->exons_gene;
// 		while(curExon != NULL)
// 		{
// 			curExon->coverage_base = new long[curExon->end_exon - curExon->start_exon + 1];
// 			for (long index = 0; index < curExon->end_exon - curExon->start_exon + 1; index++)
// 			{
// 				curExon->coverage_base[index] = 0;
// 			}
// 			curExon = curExon->next_exon;
// 		}
// 		cout << 333 << endl;
// 		process_exons();
// 		cout << 333333 << endl;
// 		if (resultBlocks != NULL)
// 		{
// 			process_block(resultBlocks);
// 		}
// 		cout << 444 << endl;
// 		sprintf(inputdata_seqbias, "seqbias.txt");
// 
// 		long index = 0;
// 		char **transName;
// 		sequencebias = new double *[geneList[1]->tranNum];
// 		transName = new char *[geneList[1]->tranNum];
// 		Transcript *trans = geneList[1]->transcripts;
// 		for (index = 0; index < geneList[1]->tranNum; index++)
// 		{
// 			sequencebias[index] = new double[trans->transLength + 1];
// 			transName[index] = new char[100];
// 			strcpy(transName[index], trans->transID);
// 			for (long iLoop = 0; iLoop <= trans->transLength; iLoop++)
// 			{
// 				sequencebias[index][iLoop] = 1;
// 			}
// 			trans = trans->next_transcript;
// 		}
// 		cout << 555 << endl;
// 		Maximization(1, resultBlocks, 190, "", sequencebias, abundance, exonabundance, 190, geneNames);
// 
// 		Transcript *curTrans = geneList[1]->transcripts;
// 		long transIndex = 0;
// 		while(curTrans != NULL)
// 		{
// 			outputfile_estimation << geneList[1]->chromosome << "\t" << geneNames << '\t' << curTrans->transID << '\t' << curTrans->transLength << endl;
// 			delete [] sequencebias[transIndex];
// 			delete [] transName[transIndex];
// 			transIndex++;
// 			curTrans = curTrans->next_transcript;
// 		}
// 		delete [] sequencebias;
// 		delete [] transName;
// 
// 		delete [] abundance;
// 		abundance = NULL;
// 		delete [] exonabundance;
// 		exonabundance = NULL;
// 
// 		Block *delblock;
// 		while(resultBlocks)
// 		{
// 			delblock = resultBlocks;
// 			resultBlocks = resultBlocks->next;
// 			delete delblock;
// 		}
// 
// 		delete geneList[1];
// 		geneList[1] = NULL;
// 		geneListNum = 0;
// 
// 		geneNames[0] = '\0';
// 		inputfile >> geneNames;
// 
// 	}
// 
// 	return 0;
// }



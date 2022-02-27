//============================================================================
// Name        : DNA.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <chrono>
#include <fstream>
#include <cassert>
#include <algorithm>
#include "Clone.hpp"
#include "Cluster2.hpp"
#include "LongestPath.hpp"
#include "EditDistance.hpp"
#include "DividerBMA.hpp"

using namespace std;

void GetAllCopies(ifstream &cluster_file, ifstream &ref_file, string &original, vector<string> &copies)
{
	string line;
	original.clear();
	copies.clear();

	// read ref from file
	getline(ref_file, original);

    // read cluster
    getline(cluster_file, line); // line containing number
    for(int i = 0; i < 10; i++) {
        getline(cluster_file, line);
        copies.push_back(line);
    }
}

struct compare
{
	bool operator()(const pair<int, int> &first, const pair<int, int> &second)
	{
		return first.first < second.first;
	}
};

void TestFromFile(const string &clusterFileName, const string &refFileName, int testNum, int strandLen, int delPatternLen,
				  const int subPriority, const int delPriority, const int insPriority, const int maxReps, const string &outputFileName)
{
	ifstream cluster_file;
	cluster_file.open(clusterFileName.c_str());

    ifstream ref_file;
	ref_file.open(refFileName.c_str());

	std::ofstream output;
    output.open(outputFileName);


	if (!cluster_file.is_open() || !ref_file.is_open())
	{
		cout << "Failed opening input file!" << endl;
		return;
	}

	string original;
	vector<string> copies;
	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);

	int cumTotalFinalGuessEditDist = 0, roundFinalGuessEditDist = 0;
	int cumFinalGuessSubstitutions = 0, cumFinalGuessInsertions = 0, cumFinalGuessDeletions = 0;
	map<int, int> editDistanceHist;

	int countFiltered = 0;
	for (int i = 1; i <= testNum; i++)
	{
		GetAllCopies(cluster_file, ref_file, original, copies);
		
		if (original.empty())
		{
			break;
		}
		Cluster2 cluster(original, copies);
		// computing the avg length of the sequence;
		double avglen = 0;
		for (string &cp : copies)
		{
			avglen += cp.length();
		}
		avglen = avglen / copies.size();

		string finalGuess;
		if (copies.size() == 1)
		{ // if only 1 copy return copy.
			finalGuess = copies[0];
		}
		else if (copies.size() >= 20 || (avglen / original.length() > 0.99 &&
										 avglen / original.length() < 1.01))
		{
			finalGuess = dividerBMA(copies, original.length());
		}
		else
		{
			finalGuess = cluster.TestBest(delPatternLen, roundFinalGuessEditDist, subPriority, delPriority, insPriority,
										  generator, maxReps);
		}

		output << finalGuess << endl;

		roundFinalGuessEditDist = ComputeEditDistanceNum(cluster.Original(), finalGuess);
		editDistanceHist[roundFinalGuessEditDist]++;
		cumTotalFinalGuessEditDist += roundFinalGuessEditDist;

		vector<LetterOps> result = ComputeEditDistancePriority(finalGuess, cluster.Original(), 0, generator);
		map<string, double> countOperations = CountOperations(result);
		assert(countOperations["I"] + countOperations["D"] + countOperations["R"] == roundFinalGuessEditDist);

		cumFinalGuessSubstitutions += countOperations["R"];
		cumFinalGuessInsertions += countOperations["I"];
		cumFinalGuessDeletions += countOperations["D"];

	}


	cluster_file.close();
    ref_file.close();
}

int getTestNum(const string &refFileName) {
    ifstream ref_file;
	ref_file.open(refFileName.c_str());

    	if (!ref_file.is_open())
	{
		cout << "Failed opening input file!" << endl;
		return 0;
	}

    int i = 0;
    while(!ref_file.eof()) {
        string line; 
        getline(ref_file, line);
        i++;
    }


    ref_file.close();

    return i;
}

// TODO:	decide fix order by copy len. too long -> prioritize fix inserts, too short -> prioritize fix deletions

int main(int argc, char *argv[])
{
	clock_t begin = clock();

	int delPatternLen = 3;

	int subPriority = 0;
	int delPriority = 0;
	int insPriority = 0;
	int maxReps = 2;


	string clusterFileName = "evyat.txt";
	string refFileName = "evyat.txt";
	string outputFileName = "evyat.txt";
	int strandLen = 120;

	if (argc > 1)
		clusterFileName = argv[1];
    if (argc > 2)
		refFileName = argv[2];
	if (argc > 3)
        outputFileName = argv[3];
	if (argc > 4)
		strandLen = atoi(argv[4]);

    
    int testNum = getTestNum(refFileName);

	TestFromFile(clusterFileName, refFileName, testNum, strandLen, delPatternLen, subPriority, delPriority, insPriority,
					maxReps, outputFileName);

	cout << "Done" << endl;
	return 0;
}



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

void TestFixAll(int testNum, int strandLen, int cloneNum, int delPatternLen, const int subPriority,
				const int delPriority, const int insPriority, const int maxReps, const double delProb, const double insProb,
				const double subProb)
{
	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);

	int cumTotalFinalGuessEditDist = 0, roundFinalGuessEditDist = 0;
	int cumFinalGuessSubstitutions = 0, cumFinalGuessInsertions = 0, cumFinalGuessDeletions = 0;
	map<int, int> editDistanceHist;

	for (int i = 0; i < testNum; i++)
	{
		Cluster2 cluster(strandLen, cloneNum, delProb, insProb, subProb, generator);

		//		cluster.TestFixAll(subPatternLen, delPatternLen, insThreshold, roundCloneOriginalEDS, roundAvgCloneEditDist,
		//				roundCumCorrectSizeCloneEditDist, roundCorrectSizeCloneNum, roundBingoNum, round1stHalfBingoNum, round2ndHalfBingoNum,roundMinED, roundMaxED,
		//				roundFinalGuessEditDist, finalGuessNum, finalStage,subPriority, delPriority, insPriority, generator, minEqual, maxReps);

		//				cluster.TestFixAllClonesLast(subPatternLen, delPatternLen, insThreshold, roundCloneOriginalEDS,
		//				roundAvgCloneEditDist, roundCumCorrectSizeCloneEditDist, roundCorrectSizeCloneNum, roundBingoNum,round1stHalfBingoNum, round2ndHalfBingoNum,
		//				roundMinED, roundMaxED, roundFinalGuessEditDist, finalGuessNum, finalStage, subPriority, delPriority,
		//				insPriority, generator, minEqual, maxReps);

		//		cluster.TestFixAllClonesLastInitial(subPatternLen, delPatternLen, insThreshold, roundCloneOriginalEDS,
		//				roundAvgCloneEditDist, roundCumCorrectSizeCloneEditDist, roundCorrectSizeCloneNum, roundBingoNum,round1stHalfBingoNum, round2ndHalfBingoNum,
		//				roundMinED, roundMaxED, roundFinalGuessEditDist, finalGuessNum, finalStage, subPriority, delPriority,
		//				insPriority, generator, minEqual, maxReps);

		//		cluster.TestFixAllByErrorType(subPatternLen, delPatternLen, insThreshold, roundCloneOriginalEDS,
		//				roundAvgCloneEditDist, roundCumCorrectSizeCloneEditDist, roundCorrectSizeCloneNum, roundBingoNum,round1stHalfBingoNum, round2ndHalfBingoNum,
		//				roundMinED, roundMaxED, roundFinalGuessEditDist, finalGuessNum, finalStage, subPriority, delPriority,
		//				insPriority, generator, minEqual, maxReps, convReps);

		string finalGuess = cluster.TestBest(delPatternLen, roundFinalGuessEditDist, subPriority, delPriority,
											 insPriority, generator, maxReps);
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
	map<int, int>::reverse_iterator rit = editDistanceHist.rbegin(); // points to last element in map
	int highestED = rit->first;
	int cumDist = 0;
	cout << "Edit distance hist:" << endl;
	for (int i = 0; i <= highestED; i++)
	{
		cumDist += editDistanceHist[i];
		cout << i << "\t" << cumDist << endl;
	}
	cout << "Avg. guess substitutions:\t" << 1000 * (double)cumFinalGuessSubstitutions / (testNum * strandLen) << endl;
	cout << "Avg. guess deletions:\t" << 1000 * (double)cumFinalGuessDeletions / (testNum * strandLen) << endl;
	cout << "Avg. guess insertions:\t" << 1000 * (double)cumFinalGuessInsertions / (testNum * strandLen) << endl;
	cout << "Avg. guess edit dist:\t" << 1000 * (double)cumTotalFinalGuessEditDist / (testNum * strandLen) << endl;
}

void TestOriginalRetention(int testNum, int strandLen, int cloneNum, int delPatternLen, const int subPriority,
						   const int delPriority, const int insPriority)
{
	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);

	int cumTotalFinalGuessEditDist = 0, roundFinalGuessEditDist = 0;

	for (int i = 0; i < testNum; i++)
	{
		Cluster2 cluster(strandLen, cloneNum, 0.1, 0.1, 0.1, generator);

		cluster.TestOriginalRetention(delPatternLen, roundFinalGuessEditDist, subPriority, delPriority, insPriority,
									  generator);

		cumTotalFinalGuessEditDist += roundFinalGuessEditDist;
	}

	cout << "Avg. original after fix edit dist:\t\t"
		 << 1000 * (double)cumTotalFinalGuessEditDist / (testNum * strandLen) << endl;
}

void TestStats(int testNum, int strandLen, int cloneNum, int subPatternLen, int delPatternLen, double insThreshold,
			   const int subPriority, const int delPriority, const int insPriority)
{
	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);
	for (int i = 0; i < testNum; i++)
	{
		Cluster2 cluster(strandLen, cloneNum, 0.1, 0.1, 0.1, generator);
		//cluster.Stats(0, subPatternLen, delPatternLen, insThreshold, subPriority, delPriority, insPriority, generator);
	}
}

void GetAllCopies(ifstream &clusters, ifstream &refs, string &original, vector<string> &copies)
{
	original.clear();
	copies.clear();

	// read reference strand
	getline(refs, original);
	if (original.empty())
	{
		return;
	}

	// read associated cluster
	string line;
	getline(clusters, line); // first line "<num_of_copies>" dump
	int num = stoi(line);
	for(int i = 1; i <= num; i++) {
		getline(clusters, line);
		if (line.empty()) continue;
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

vector<string> SortCopiesByAbsLenDiff(const vector<string> &copies, const int correctSize)
{
	compare c;
	int copiesNum = copies.size();
	vector<pair<int, int> > copiesLen(copiesNum);
	for (int i = 0; i < copiesNum; i++)
	{
		int currentSize = copies[i].size();
		copiesLen[i].first = abs(currentSize - correctSize);
		copiesLen[i].second = i;
	}
	sort(copiesLen.begin(), copiesLen.end(), c);
	vector<string> sortedCopies;
	for (int i = 0; i < copiesNum; i++)
	{
		int oldIndex = copiesLen[i].second;
		sortedCopies.push_back(copies[oldIndex]);
	}
	return sortedCopies;
}

void GetCaseWithCopiesLimit(ifstream &input, string &original, vector<string> &copies, const int maxCopies)
{
	string line;
	original.clear();
	copies.clear();
	getline(input, line);
	if (line.empty())
	{
		return;
	}
	// first line is original
	original = line;

	// second line "*****" dump
	getline(input, line);

	// each line a copy until 2 empty lines.
	int endCase = 0;
	int copyIndex = 0;
	while (getline(input, line))
	{
		if (line.empty())
		{
			endCase++;
		}
		else
		{
			if (copyIndex < maxCopies)
			{
				copies.push_back(line);
				copyIndex++;
			}
		}
		if (endCase == 2)
		{
			break;
		}
	}
}

void TestFromFile(const string &inputFilename, int testNum, int strandLen, int maxCopies, int delPatternLen,
				  const int subPriority, const int delPriority, const int insPriority, const int maxReps)
{
	ifstream input;
	input.open(inputFilename.c_str());
	if (!input.is_open())
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
		GetCaseWithCopiesLimit(input, original, copies, maxCopies);
		//		set<int>::iterator it = filter.find(i);
		//		if (it != filter.end()) {
		//			cout << "filtered case #" << i << endl;
		//			countFiltered++;
		//			continue;
		//		}
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
			cout << i << endl;
			finalGuess = dividerBMA(copies, original.length());
		}
		else
		{
			finalGuess = cluster.TestBest(delPatternLen, roundFinalGuessEditDist, subPriority, delPriority, insPriority,
										  generator, maxReps);
		}
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
	map<int, int>::reverse_iterator rit = editDistanceHist.rbegin(); // points to last element in map
	int highestED = rit->first;
	int cumDist = 0;
	cout << "Edit distance hist:" << endl;
	for (int i = 0; i <= highestED; i++)
	{
		cumDist += editDistanceHist[i];
		cout << i << "\t" << cumDist << endl;
	}
	cout << "Number of filtered:\t" << countFiltered << endl;
	cout << "Avg. guess substitutions:\t"
		 << 1000 * (double)cumFinalGuessSubstitutions / ((testNum - countFiltered) * strandLen) << endl;
	cout << "Avg. guess deletions:\t"
		 << 1000 * (double)cumFinalGuessDeletions / ((testNum - countFiltered) * strandLen) << endl;
	cout << "Avg. guess insertions:\t"
		 << 1000 * (double)cumFinalGuessInsertions / ((testNum - countFiltered) * strandLen) << endl;
	cout << "Avg. guess edit dist:\t"
		 << 1000 * (double)cumTotalFinalGuessEditDist / ((testNum - countFiltered) * strandLen) << endl;

	input.close();
}

// counting from 1
// test [startCase,endCase]
void TestFromFileCaseRange(const string &clusterFileName, const string &refsFileName,
						   const string &outputFilename,
						   int testNum,
						   int strandLen, int maxCopies, int delPatternLen, const int subPriority, const int delPriority,
						   const int insPriority, const int maxReps)
{
	ifstream clusters;
	clusters.open(clusterFileName.c_str());
	if (!clusters.is_open())
	{
		cout << "Failed opening cluster file!" << endl;
		return;
	}

	ifstream refs;
	refs.open(refsFileName.c_str());
	if (!refs.is_open())
	{
		cout << "Failed opening refs file!" << endl;
		return;
	}

	ofstream output;
	output.open(outputFilename.c_str());
	if (not output.is_open())
	{
		cout << "Error opening output file!" << endl;
		return;
	}

	string original;
	vector<string> copies;


	unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator(sd);
	int roundFinalGuessEditDist = 0;
	double cumTotalFinalGuessEditDist = 0;
	double cumFinalGuessSubstitutions = 0, cumFinalGuessInsertions = 0, cumFinalGuessDeletions = 0;
	double error_rate = 0.0;
	map<int, int> editDistanceHist;
	int majoCounter = 0;
	for (int i = 1; i <= testNum; i++)
	{
		GetAllCopies(clusters, refs, original, copies);
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
		/*else if(copies.size() >= 20 || (avglen/original.length()>0.99 &&
            avglen/original.length()<1.01)){
            majoCounter++;
            finalGuess=dividerBMA(copies, original.length());
            //finalGuess=laMajority(copies, original.length());
        }*/
		else
		{
			finalGuess = cluster.TestBest(delPatternLen, roundFinalGuessEditDist, subPriority, delPriority, insPriority,
										  generator, maxReps);
		}

		roundFinalGuessEditDist = ComputeEditDistanceNum(cluster.Original(), finalGuess);
		editDistanceHist[roundFinalGuessEditDist]++;
		cumTotalFinalGuessEditDist += roundFinalGuessEditDist;

		output << finalGuess << endl;

		vector<LetterOps> result = ComputeEditDistancePriority(finalGuess, cluster.Original(), 0, generator);
		map<string, double> countOperations = CountOperations(result);
		assert(countOperations["I"] + countOperations["D"] + countOperations["R"] == roundFinalGuessEditDist);

		if (original.length())
		{
			cumTotalFinalGuessEditDist = ((i - 1) * cumTotalFinalGuessEditDist + double(roundFinalGuessEditDist) / (original.length())) / i;
			cumFinalGuessSubstitutions = ((i - 1) * cumFinalGuessSubstitutions + (countOperations["R"]) / (original.length())) / i;
			cumFinalGuessInsertions = ((i - 1) * cumFinalGuessInsertions + (countOperations["I"]) / (original.length())) / i;
			cumFinalGuessDeletions = ((i - 1) * cumFinalGuessDeletions + ((countOperations["D"]) / (original.length()))) / i;
			error_rate = ((i - 1) * error_rate + (double(roundFinalGuessEditDist) / (original.length()))) / i;
		}
	}

	cout << error_rate << endl;


	clusters.close();
	refs.close();
	output.close();
}

void CasesFromFile(const string &inputFilename, const int maxCopies)
{
	ifstream input;
	input.open(inputFilename.c_str());
	if (!input.is_open())
	{
		cout << "Failed opening input file!" << endl;
		return;
	}
	int countCases = 0;
	int minCopiesNum = 1000, maxCopiesNum = 0, currentCopiesNum, cumCopiesNum = 0;
	string original;
	vector<string> copies;
	while (1)
	{
		GetCaseWithCopiesLimit(input, original, copies, maxCopies);
		if (original.empty())
		{
			break;
		}
		countCases++;
		currentCopiesNum = copies.size();
		cumCopiesNum += currentCopiesNum;
		if (currentCopiesNum < minCopiesNum)
		{
			minCopiesNum = currentCopiesNum;
		}
		if (currentCopiesNum > maxCopiesNum)
		{
			maxCopiesNum = currentCopiesNum;
		}
		if (original.size() != 150)
		{
			cout << countCases << "\t" << currentCopiesNum << "\t" << original.size() << endl;
		}
	}
	cout << "min copies num: " << minCopiesNum << endl;
	cout << "max copies num: " << maxCopiesNum << endl;
	cout << "average copies num: " << (double)cumCopiesNum / countCases << endl;
	cout << "cases num: " << countCases << endl;
	input.close();
}

double CalcMed(vector<int> scores)
{
	size_t size = scores.size();

	if (size == 0)
	{
		return 0; // Undefined, really.
	}
	else
	{
		sort(scores.begin(), scores.end());
		if (size % 2 == 0)
		{
			return (scores[size / 2 - 1] + scores[size / 2]) / 2;
		}
		else
		{
			return scores[size / 2];
		}
	}
}

double AVGEditDistance(const string &original, const vector<string> &copies, int &minED, int &maxED, double &medED)
{
	double cumEditDistance = 0;
	int currentED;
	minED = 1000;
	maxED = 0;
	vector<int> eds;
	for (unsigned i = 0; i < copies.size(); i++)
	{
		currentED = ComputeEditDistanceNum(original, copies[i]);
		eds.push_back(currentED);
		cumEditDistance += currentED;
		if (currentED > maxED)
		{
			maxED = currentED;
		}
		if (currentED < minED)
		{
			minED = currentED;
		}
	}
	medED = CalcMed(eds);
	return cumEditDistance / copies.size();
}

void CasesStats(const string &inputFilename, const int maxCopies, const double EDThresh, const set<int> &filter)
{
	ifstream input;
	input.open(inputFilename.c_str());
	if (!input.is_open())
	{
		cout << "Failed opening input file!" << endl;
		return;
	}
	int countCases = 0;
	int currentCopiesNum, cumCopiesNum = 0;
	int minED = 0, maxED = 0;
	double medED;
	string original;
	vector<string> copies;
	while (1)
	{
		GetCaseWithCopiesLimit(input, original, copies, maxCopies);
		if (original.empty())
		{
			break;
		}
		countCases++;
		currentCopiesNum = copies.size();
		cumCopiesNum += currentCopiesNum;
		double currentAVGED = AVGEditDistance(original, copies, minED, maxED, medED);

		//		if (medED > EDThresh or (currentCopiesNum < 10 and currentAVGED > 10) or currentCopiesNum == 1) {
		//		if (currentAVGED > 1.5 or currentCopiesNum < 10) {
		set<int>::iterator it = filter.find(countCases);
		if (it != filter.end())
		{
			cout << countCases << "\t" << currentCopiesNum << "\t" << minED << "\t" << maxED << "\t" << currentAVGED
				 << "\t" << medED << endl;
		}
	}
	//	cout << "min copies num: " << minCopiesNum << endl;
	//	cout << "max copies num: " << maxCopiesNum << endl;
	//	cout << "average copies num: " << (double) cumCopiesNum / countCases << endl;
	//	cout << "cases num: " << countCases << endl;
	input.close();
}

int getTestNum(const string &refsFileName)
{
	ifstream refs;
	refs.open(refsFileName.c_str());
	if (!refs.is_open())
	{
		cout << "Failed opening input file!" << endl;
		return -1;
	}
	string line;
	int count = 0;
	while (!refs.eof()) {
		getline(refs, line);
		count++;
	}
	
	refs.close();

	return count;
}

// TODO:	decide fix order by copy len. too long -> prioritize fix inserts, too short -> prioritize fix deletions

int main(int argc, char *argv[])
{
	clock_t begin = clock();

	int maxCopies = 25;
	int delPatternLen = 3;

	int subPriority = 0;
	int delPriority = 0;
	int insPriority = 0;

	int maxReps = 2;


	string clusterFileName = "evyat.txt";
	string refsFileName = "evyat.txt";
	string outputFileName = "output.txt";
	int strandLen = 120;

	if (argc > 1)
		clusterFileName = argv[1];
	if (argc > 2)
		refsFileName = argv[2];
	if (argc > 3)
		strandLen = atoi(argv[3]);
	if (argc > 4)
		outputFileName = argv[4];
	
	cout << "Doing\n";
 
	int testNum = getTestNum(refsFileName);

	cout << "Done\n";


	TestFromFileCaseRange(clusterFileName, refsFileName, outputFileName, testNum, strandLen, maxCopies, delPatternLen,
						  subPriority, delPriority, insPriority, maxReps);


	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << endl;
	cout << "Time elapsed: " << (int)elapsed_secs << "\tseconds" << endl;
	return 0;
}

//	double delProb;
//
//	delProb = 0.005;
//	cout << "Test num:\t" << testNum << endl;
//	cout << "Strand len:\t" << strandLen << endl;
//	cout << "Cluster size:\t" << cloneNum << endl;
//	cout << "All prob:\t" << delProb << endl;
//	TestFixAll(testNum, strandLen, cloneNum, delPatternLen, subPriority, delPriority, insPriority, maxReps, delProb,
//			delProb, delProb);
//	cout << "*******************************" << endl;
//
//	for (int i = 1; i <= 10; i++) {
//		delProb = (double) i / 100;
//		cout << "Test num:\t" << testNum << endl;
//		cout << "Strand len:\t" << strandLen << endl;
//		cout << "Cluster size:\t" << cloneNum << endl;
//		cout << "All prob:\t" << delProb << endl;
//		TestFixAll(testNum, strandLen, cloneNum, delPatternLen, subPriority, delPriority, insPriority, maxReps, delProb,
//				delProb, delProb);
//		cout << "*******************************" << endl;
//	}

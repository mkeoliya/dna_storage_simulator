
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <map>
#include <set>

#include <string>
#include <iterator>
#include <time.h>
#include <math.h>
#include <random>
#include <chrono>
#include <cstdint>
#include <climits>
#include <fstream>
#include <stdlib.h>

using namespace std;
using namespace std::chrono;

#define DES_LEN 110
#define HISTO_LEN 1000
#define MAX 200
#define SAME 0
#define DEL 1
#define INS 2
#define SUB 3
#define NUMBER_OF_BP 4

int lookup[MAX][MAX];

string _correct_binary_indel(int n, int m, int a, string y);
string _correct_q_ary_indel(int n, int m, int a, int b, int q, string y);
string convert_y_to_alpha(string y);
string genRandomSeq(int len);
vector<int> compute_indices(char y1, string x);
int count_embeddings(string a, string b);
string maximum_likelihood_old(set<string> &candidates, const vector<string> &cluster);
string maximum_likelihood(set<string> &candidates, const vector<string> &cluster);
void createNoisyCluster(vector<string> &cluster, string original, int cluster_size, float del_prob);
string reconstruction_new(vector<string> &cluster, int des_len);
string reconstruction_new_alg4(vector<string> &cluster, int des_len);
void reomoveByLength(set<string> &candidates, int len);
void removeNonSuperseq_withLen(set<string> &candidates, vector<string> &cluster, int len);
void removeNonSuperseq(set<string> &candidates, vector<string> &cluster, int len);
int compute_syndrome_binary(int m, int a, string y);
string convert_y_to_alpha(string y);
int compute_syndrome_q_ary(int m, int a, int b, int q, string y, int &sec);
bool is_codeword(string &y, int q);
string decode_VT(string y, int q, int n, int m, int a, int b);

vector<string> LCS(string X, string Y, int m, int n);

//set<string> LCS(string X, string Y, int m, int n);
void LCSLength(string X, string Y, int m, int n);

int design_len = DES_LEN;
int total_dis_tests = 0;
float prob = 0.8;
int clus_size = 6;
int num_of_test = 100000;
int count = 0;
int q = 4;
int scs_size = 0;
int scs_after_filter_size = 0;
int scs_most_likely_size = 0;
int most_likelihood_count = 0;
int pairs_histo[HISTO_LEN] = {0};
int three_histo[HISTO_LEN] = {0};
int four_histo[HISTO_LEN] = {0};
int SR_histo[DES_LEN] = {0};
int SR_histo_BMA[DES_LEN] = {0};
int SR_histo_MATIKA[DES_LEN] = {0};
int SR_histo_D_R[DES_LEN] = {0};
int SR_histo_karin[DES_LEN] = {0};

double i_prob = 0.01;
double d_prob = 0.01;
double s_prob = 0.01;

set<string> ML;
map<int, int> scs_sizes;
map<int, int> scs_sizes_count;
map<int, int> mst_like;
map<int, int> mst_like_count;
//int table3[4][4]={{0}};
//int table4[4][4]={{0}};
//int table_l_3={0,1,3,4}=;
//int table_r_3=;
double error_rate_by_length[DES_LEN + 50] = {0};
double count_by_length[DES_LEN + 50] = {0};
map<string, int> table_1_rev;
map<string, int> table_2_rev;
vector<int> systematic_positions_step_1;
int lucky = 0;
int two = 0;
int three = 0;
int four = 0;
int nothing = 0;
ofstream myfile;


int DNAtoInt(char c)
{
    if (c == 'A')
        return 0;
    if (c == 'C')
        return 1;
    if (c == 'G')
        return 2;
    if (c == 'T')
        return 3;
    return -1;
}

char letToBp(int i)
{
    switch (i)
    {
    case 0:
        return 'A';
    case 1:
        return 'C';
    case 2:
        return 'G';
    case 3:
        return 'T';
    }
    return 'N';
}

int letToInt(char c)
{
    switch (c)
    {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    }
    return 4;
}

struct LetterOps
{
    std::string insert;  // string to insert before letter
    std::string CDR;     // Copy or Delete or Substitute letter
    std::string opParam; // Which letter to copy or which to delete or which to substitute with
};

unsigned int edit_distance(const std::string &s1, const std::string &s2)
{
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<std::vector<unsigned int> > d(len1 + 1, std::vector<unsigned int>(len2 + 1));

    d[0][0] = 0;
    for (unsigned int i = 1; i <= len1; ++i)
        d[i][0] = i;
    for (unsigned int i = 1; i <= len2; ++i)
        d[0][i] = i;

    for (unsigned int i = 1; i <= len1; ++i)
        for (unsigned int j = 1; j <= len2; ++j)
            // note that std::min({arg1, arg2, arg3}) works only in C++11,
            // for C++98 use std::min(std::min(arg1, arg2), arg3)
            d[i][j] = std::min({d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1)});
    return d[len1][len2];
}

void reverseStr(string &str)
{
    int n = str.length();

    // Swap character starting from two
    // corners
    for (int i = 0; i < n / 2; i++)
        swap(str[i], str[n - i - 1]);
}

int get_max(vector<int> &ptrs)
{
    int max_i = ptrs[0];
    for (int p : ptrs)
    {
        if (p > max_i)
        {
            max_i = p;
        }
    }
    return max_i;
}

char getMajorityByPointers(vector<string> &cluster, vector<int> &ptrs)
{
    int cl_size = cluster.size();
    int majo[NUMBER_OF_BP] = {0, 0, 0, 0};
    for (int i = 0; i < cl_size; i++)
    {
        if (ptrs[i] < cluster[i].length())
            majo[letToInt(cluster[i][ptrs[i]])]++;
    }
    int max_i = 0;
    int max_value = majo[max_i];
    for (int i = 1; i < NUMBER_OF_BP; i++)
    {
        if (majo[i] > max_value)
        {
            max_i = i;
            max_value = majo[i];
        }
    }

    return letToBp(max_i);
}


char getMajorityByPointers_window(vector<string> &cluster, vector<int> &ptrs, int w, vector<int> variant_reads)
{
    int cl_size = cluster.size();
    int majo[NUMBER_OF_BP + 1] = {0, 0, 0, 0, 0}; // the last position refer to the case of an emptystring
    for (int i = 0; i < cl_size; i++)
    {
        if (find(variant_reads.begin(), variant_reads.end(), i) != variant_reads.end())
        {
            continue;
        }
        if (ptrs[i] + w < cluster[i].length())
        {
            majo[letToInt(cluster[i][ptrs[i] + w])]++;
        }
        else
        {
            continue;
        }
    }
    int max_i = 0;
    int max_value = majo[max_i];
    for (int i = 1; i < NUMBER_OF_BP; i++)
    {
        if (majo[i] > max_value)
        {
            max_i = i;
            max_value = majo[i];
        }
    }
    if (max_i < NUMBER_OF_BP)
        return letToBp(max_i);
    else
        return 0;
}

string majority_window(vector<string> &cluster, int des_len, int w, vector<int> &ptrs, vector<int> variant_reads)
{
    string window_maj = "";
    for (int i = 0; i <= w; i++)
    {
        char c = getMajorityByPointers_window(cluster, ptrs, i, variant_reads);
        window_maj.push_back(c);
        if (c == 0)
            break;
    }
    return window_maj;
}

string LOOKAhead_Majority(vector<string> &cluster, int des_len, int w)
{
    string output = "";
    char base;
    char prev_base = 'N';
    vector<int> ptrs;
    for (int i = 0; i < cluster.size(); i++)
    {
        ptrs.push_back(0);
    }
    for (int i = 0; i < des_len; i++)
    {
        if (get_max(ptrs) > des_len)
        {
            break;
        }
        base = getMajorityByPointers(cluster, ptrs);
        vector<int> variant_reads;
        for (int j = 0; j < cluster.size(); j++)
        {
            if (base != cluster[j][ptrs[j]])
            {
                variant_reads.push_back(j);
            }
        }
        output.push_back(base);
        string window_maj = majority_window(cluster, des_len, w, ptrs, variant_reads);
        int new_w = window_maj.length();
        for (int j = 0; j < cluster.size(); j++)
        {
            if (ptrs[j] >= cluster[j].length())
            {
                continue;
            }
            if (cluster[j][ptrs[j]] == base)
            {
                ptrs[j] += 1;
                continue; // Case #1: No error in trace j
            }
            else
            {
                string current_window;
                int k = 0;
                for (k = 0; k <= w && cluster[j][ptrs[j] + k] != 0; k++)
                {
                    current_window.push_back(cluster[j][ptrs[j] + k]);
                }
                int max_w = min(k, new_w);
                if (max_w <= 1)
                {
                    ptrs[j] += 1;
                    continue;
                }
                /* case of substitution */
                if (current_window.compare(1, max_w - 1, window_maj, 1, max_w - 1) == 0)
                {
                    ptrs[j] += 1;
                    continue;
                }
                /* case of insertion */
                if (current_window.compare(1, max_w - 1, window_maj, 0, max_w - 1) == 0)
                {
                    ptrs[j] += 2;
                    continue;
                }

                /* case of deletion */
                if (current_window.compare(0, max_w - 1, window_maj, 1, max_w - 1) == 0)
                {
                    ptrs[j] += 0;
                    continue;
                }

                // Second Check by Edit Distance.
                /* case of substitution */
                if (edit_distance(current_window.substr(1, max_w - 1), window_maj.substr(1, max_w - 1)) < 2)
                {
                    ptrs[j] += 1;
                    continue;
                }
                /* case of insertion */
                if (edit_distance(current_window.substr(1, max_w - 1), window_maj.substr(0, max_w - 1)) < 2)
                {
                    ptrs[j] += 2;
                    continue;
                }

                /* case of deletion */
                if (edit_distance(current_window.substr(0, max_w - 1), window_maj.substr(1, max_w - 1)) < 2)
                {
                    ptrs[j] += 0;
                    continue;
                }

                ptrs[j] += 1;
            }
        }
        prev_base = base;
    }
    return output;
}


// main function
int main(int argc, char *argv[])
{
    auto start = high_resolution_clock::now();

    double succes_rate = 0.0;

    // default files
    string clusterFileName = "evyat.txt";
    string outputFileName = "output.txt";

    if (argc > 1)
        clusterFileName = argv[1];
    if (argc > 2)
        outputFileName = argv[2];
    
    try
    {
        int counterOfGoodLen = 0;
        string original;
        vector<string> cluster;

        double res_del;
        double res_ins;
        double res_sub;

        
            unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
            mt19937 generator(sd);
            int roundFinalGuessEditDist = 0;
            double cumTotalFinalGuessEditDist = 0;
            double cumFinalGuessSubstitutions = 0, cumFinalGuessInsertions = 0, cumFinalGuessDeletions = 0;
            double error_rate = 0.0;
            map<int, int> editDistanceHist;


            cluster.clear();

            // open the cluster file
            std::ifstream cluster_file;
            cluster_file.open(clusterFileName);

            // open the output file
            std::ofstream output_file;
            output_file.open(outputFileName);

            if (!cluster_file)
            {
                return -1;
            }

            int i = 0;

            while (!cluster_file.eof())
            {
                cluster.clear();

                if (cluster_file.eof())
                {
                    break;
                }

                // Read cluster
                string line;
                getline(cluster_file, line); // line containing number
                int num = stoi(line);
                for(int i = 0; i < num; i++) {
                    getline(cluster_file, line);
                    cluster.push_back(line);
                }


                // Reverse cluster
                vector<string> reverse_cluster;
                for (string cp : cluster)
                {
                    reverseStr(cp);
                    reverse_cluster.push_back(cp);
                }

                // Compute BMA reconstruction
                string karin_w3 = LOOKAhead_Majority(cluster, DES_LEN, 3);
                string karin_rev_w3 = LOOKAhead_Majority(reverse_cluster, DES_LEN, 3);

                string nkarin_w3 = "";
                for (int k = 0; k < DES_LEN; k++)
                {
                    if (k < DES_LEN / 2)
                        nkarin_w3.insert(nkarin_w3.begin(), karin_rev_w3[k]);
                }
                for (int k = DES_LEN / 2 - 1; k >= 0; k--)
                {
                    nkarin_w3.insert(nkarin_w3.begin(), karin_w3[k]);
                }

                karin_w3 = nkarin_w3;
                string recon = karin_w3;

                output_file << recon << endl;
                
                i++;
            }

            cout << "Total clusters: " << i << endl;
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
        cerr << e.what() << endl;
    }
    

    return 0;
}

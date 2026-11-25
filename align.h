#include<stdio.h>
#include<list>
#include<unordered_map>
#include<algorithm>
#include<string>
#include<string.h>
#include<vector>
#include "func.h"
#include<stdlib.h>
using namespace std;

std::tuple<std::vector<std::array<float, 2>>, std::vector<uint32_t>> cluster(vector<read> *split);
inline void check_concat(list<data>::iterator &iter, list<data> &be, list<data> &af, data tmp,int &is_find,int K);
void read_classify(read &r, vector<data> &res);
void max_chunk(list<data> &chunk, int K, vector<data> &res,read &r);
void push_chunk(list<data> &chunk,list<data> &now,int K);
void make_chunk(list<data> &chunk,char *query, int K, unordered_map<long long int,vector<int>> &kmer_ref,int l);
void align_query(int K, unordered_map<long long int,vector<int>> &kmer_ref, read &r,unordered_map<long long int,vector<int>> &kmer_ref2,int flag2);
void print_result(vector<data> &res, vector<cigar> &cigar_string, list<data> &chunk, read &r);

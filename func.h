#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <string.h>
#include <list>
#include <vector>
#include <unordered_map>
using namespace std;

#ifndef ST
#define ST


struct cigar
{
	char type;
	int len;
};
struct read
{
	char query_name[300];
	char query[300];
	char query_qual[300];
	vector<cigar> cigar_string;
	int map_f;  //=-1 un_mapped ;; 1 = perfect ;; 2 = split 
	int ref_st;
	int ref_ed;
	int overlap;
	int query_len;
	int map_q;
	int sam_f;
	int pos;
	int p_next;
	int pair_st;
	int tlen;
	char strand;
	int insert_size;
	int tg_miss;
	int tg_match;
	int del_size;
	int dist; // dist from clevege
	int ind; // mutual exclusive? 1 = no 0 = yes
	bool ins,del,miss;
	vector<string> tags;
	string qual;
//	string query_qual;
	string ref;
	
};

struct result
{
	int normal;
	int long_d;
	int t_align;
	int tg_align_0,tg_align_1,tg_align_2;
	int tg_ins,tg_ins1,tg_ins2;
	int tg_del,tg_del1,tg_del2;
	int tg_miss,tg_miss1,tg_miss2;
	int only;
	int large_d;
};
struct data
{
    int st,ed;
    int st2,ed2;
    int miss;
};

struct dataComparator{
bool operator ()(const data &w, const data &v){
			if(w.st==v.st)
			{
				if(w.ed==v.ed){return w.st2<v.st2;}
				return w.ed<v.ed;
			}
			return w.st<v.st;}
};
#endif

void parse_read(read &r);
void reverse_comp(char *query,int l);
char* num_to_char(char *query,int l);
void char_to_num(char *query,int l);
void make_sam(read &r1, read &r2,int flag,string gene_name);
int long_inner_check(read &r1,read &r2,vector<read> *split,result &res);
int make_ref(string refname,int K, unordered_map<long long int,vector<int>> &kmer_ref);
void print_read(read &r,FILE *f);
string rgb(string hex);

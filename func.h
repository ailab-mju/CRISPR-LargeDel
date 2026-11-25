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

/* Open a file that may be gzip/bzip2/xz compressed.  A compressed file is
   piped through the matching decompressor, so the rest of the code keeps
   using plain fgets().  *is_pipe is set to 1 for a pipe and must be handed
   back to ka_close().  Decompression then runs in its own process, which is
   also faster than decompressing inline. */
FILE *ka_open(const char *path,int *is_pipe);
void  ka_close(FILE *f,int is_pipe);


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
	bool sclip;    // terminal soft clip near the cut site -> insertion evidence
	int clip_len;  // length of that clip, 0 when sclip is 0
	char clip_side;// 'L' or 'R' (reference orientation), 0 when sclip is 0
	string pat;    // mutation pattern this pair was classified as; appended to
	               // the read name in every SAM this tool writes
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
	int tg_clip;   // pairs whose only evidence is a terminal soft clip
	int tg_nocov;  // pairs that straddle the cut but whose mates never cover it
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
int mate_gap_call(const read &r1,const read &r2);
int make_ref(string refname,int K, unordered_map<long long int,vector<int>> &kmer_ref);
std::string sam_record(read &r);   /* one SAM record, no trailing newline */
bool        read_aligned(const read &r);
string rgb(string hex);
/* clip rescue (see func.cpp): terminal-clip pairs are held back until the end
   of a sample and re-assigned to small indels seen there with full anchors */
struct SamOut;
void clip_rescue_observe(const read &r1,const read &r2);
int  clip_rescue_defer(const read &r1,const read &r2);
int  clip_rescue_waiting_full();
void clip_rescue_attach_full(const read *r1,const read *r2);
void clip_rescue_finish(result &res,SamOut *full);

#include<stack>
#include<stdio.h>
#include<limits>
#include<list>
#include<unordered_map>
#include<algorithm>
#include<string>
#include<string.h>
#include<vector>
#include<stdlib.h>
#include<functional>
#include "arg_parse.h"
#include "align.h"
#include "align_wrap.h"
#define MIN_match 13
#define MIN_skip 50
using namespace std;
extern char g_ref_seq[100000],view_q[1000];
extern Args arg;
extern int skip_read;
void side_align(vector<data> &res,read &r)
{
	if(res.size()==0){return;}
	data tmp;
	char t[5];
	t[0] = 'A',t[1] = 'C',t[2] = 'G',t[3] = 'T';
	int i,l=res.size(),last=-1,p,miss=0,min,mini,j,sw=0;
	if(arg.targetst+arg.cle-20 <= res[0].st && res[0].st<=arg.targetst+arg.cle+20)
	{
		for(i=res[0].st2;i<=res[0].ed2;i++)
		{
			if(i-res[0].st2>20){break;}
			p=res[0].st+(i-res[0].st2);
			if(g_ref_seq[p]!=t[r.query[i]])
			{
				miss++;
				last=i;
			}
		}
		min=miss;
		miss=0;
		for(i=-1;i>=-30;i--)
		{
			miss=0;
			for(j=0;j<=last;j++)
			{
				if(t[r.query[j]]!=g_ref_seq[res[0].st+i+j])
				{
					miss++;
				}
			}
			if(miss<min)
			{
				sw=1;
				min = miss;
				mini = i;
				tmp.st2 = 0;
				tmp.ed2 = last;
				tmp.miss=min;
				tmp.st = res[0].st+mini;
				tmp.ed = res[0].st+mini+last;
			}
		}
		if(sw){
		res[0].st = res[0].st+(tmp.ed2+1-res[0].st2);
		res[0].st2 = tmp.ed2+1;
		res.push_back(tmp);
		for(i=res.size()-1;i>=1;i--)
		{
			res[i] = res[i-1];
			}
		res[0] = tmp;
		}
//		printf("%d %d\n",min,mini);
//		printf("%d %d %d %d\n",tmp.st,tmp.ed,tmp.st2,tmp.ed2);

	}
	else if(arg.targetst+arg.cle-20 <= res[l-1].ed && res[l-1].ed<=arg.targetst+arg.cle+20)
	{
		for(i=res[l-1].ed2;i>=res[l-1].st2;i--)
		{
			if(res[l-1].ed2-i>20){break;}
			p=(res[l-1].ed-(res[l-1].ed2-i));
			if(g_ref_seq[p]!=t[r.query[i]])
			{
				miss++;
				last=i;
			}
		}
		min=miss;
		miss=0;
		for(i=1;i<=30;i++)
		{
			miss=0;
			for(j=last;j<=res[l-1].ed2;j++)
			{
				if(t[r.query[j]]!=g_ref_seq[res[l-1].st+i+j])
				{
					miss++;
				}
			}
			if(miss<min)
			{
				sw=1;
				min = miss;
				mini = i;
				tmp.st2 = last;
				tmp.ed2 = res[l-1].ed2;
				tmp.miss=min;
				tmp.st = res[l-1].ed+mini-(res[l-1].ed2-last);
				tmp.ed = res[l-1].ed+mini;
			}
		}
		if(sw){
		res[l-1].ed = res[l-1].ed-(res[l-1].ed2-last)-1;
		res[l-1].ed2 = last-1;
		res.push_back(tmp);}

//		printf("%d %d\n",min,mini);
//		printf("%d %d %d %d\n",tmp.st,tmp.ed,tmp.st2,tmp.ed2);

	}
	else
	{
		return;
	}
	if(last==-1){return;}

//	printf("%d\n",last);

}
void read_classify(read &r,vector<data> &res)
{
	int i;
	r.ins=r.del=0;
	int m=0,pos=0,dist=50000,skip=-1,n_del=0;
	if(r.cigar_string.empty() || res.empty()|| r.cigar_string[0].len==-1){r.map_f=-1;r.ref_st = 0;r.ref_ed = 0; return;}
	vector<int> match;
	cigar now;
	pos = res[0].st+1;
	if(pos<40){r.map_f=-1;r.ref_st=0;r.ref_ed=0;return;}
	r.ref_st = pos;
	for(i=0;i<r.cigar_string.size();i++)
	{
		dist = min(dist, abs(pos-arg.targetst-arg.cle));
		now = r.cigar_string[i];
		switch(now.type){
			case 'M':
				m+=now.len;
				pos+=now.len;
				match.push_back(now.len);
				break;
			case 'D':
				if(now.len>=MIN_skip){n_del++;}
				if(!(pos+now.len<arg.targetst+arg.cle-13|| arg.targetst+arg.cle+13<pos)){r.del=1;}
				skip = max(skip,now.len);
				pos+=now.len;
				break;
			case 'I':
				if(!(pos<arg.targetst+arg.cle-13|| arg.targetst+arg.cle+13<pos-now.len)){r.ins=1;}
				pos-=now.len;
				break;
			default:
				break;
			}

	}
	r.ref_ed = pos;
	dist = min(dist, abs(pos-arg.targetst));
	dist = min(dist, abs(pos-arg.targeted));
	r.dist = 50000;
	sort(match.begin(),match.end(), greater<int>() );
	if(n_del>=2){r.map_f=-1;return;}
	if(r.ref_ed < arg.targetst || arg.targeted < r.ref_st){
		r.overlap = 0;
	}
	else{
		if(dist<arg.dist_cut){
			r.overlap = 1;
		}
		else
		{
			r.overlap = 0;
		}
	}
	if(match.size()==1){r.map_f = 1;return;}
	if(skip>=MIN_skip && match[1] >= MIN_match && r.overlap ==1)
	{
		r.dist = dist;
		if(abs(arg.ref_len-r.ref_ed)<40)
		{
			r.map_f=-1;r.ref_st=0;r.ref_ed=0;return;
		}
		r.map_f=2;
		return;
	}
	else if(skip>=MIN_skip && match[1] >=MIN_match && r.overlap==0)
	{
		r.map_f=-1;
		return;
	}
	r.map_f=1;
}
int pre_chunk(char *query, unordered_map<long long int,vector<int>> &kmer_ref2,int l)
{
	int flag=0,i;
	int K=17;
	long long int h=0,last_bit=(1LL<<(2*(K-1))),h1=0,hh1=0,hh2=0;
	for(i=0;i<K;i++)
	{
		h*=4;
		h+=query[i];
	}
	for(i=l-K;i<l;i++)
	{
		h1*=4;
		h1+=query[i];
	}
	hh1 = h;
	hh2 = h1;
	if(kmer_ref2[h].size()!=0&&kmer_ref2[h][0]>arg.targeted+50){flag=kmer_ref2[h][0]; return flag;}
	if(kmer_ref2[h1].size()!=0&&kmer_ref2[h1][0]+K<arg.targetst-50){flag=kmer_ref2[h1][0]-l+K; return flag;}
	if(kmer_ref2[h].size()==0||kmer_ref2[h1].size()==0){return -1;}
	if(abs(l-(kmer_ref2[h1][0]-kmer_ref2[h][0]+K)) <= 3)
	{
        flag=-1;
        //flag=kmer_ref2[h][0];
		return flag;
	}
	else
	{
		flag=-1;
		return flag;
	}
}
void LIS(vector<data> &kmer_list,char *query, int K, unordered_map<long long int,vector<int>> &kmer_ref,int l,vector<int> &d, vector<int> &via)
{
	d.clear();
	via.clear();
	int i;
	long long int h=0,last_bit=(1LL<<(2*(K-1)));
	vector<data> res;
	data tmp;
	for(i=0;i<=l;i++)
	{
		if(i>=K)
		{
			for(int j=0;j<kmer_ref[h].size();j++)
			{
				tmp.st = kmer_ref[h][j],tmp.ed = kmer_ref[h][j]+K-1; //ref 위치
				tmp.st2 = i-K,tmp.ed2 = i-1; //query 위치
				kmer_list.push_back(tmp);
			}
			h-=last_bit*query[i-K];
		}
		h*=4;
		h+=query[i];
	}
	sort(kmer_list.begin(),kmer_list.end(),dataComparator());
	d.push_back(-1);
	for(i=0;i<kmer_list.size();i++)
	{
		if(d.back()<kmer_list[i].st2)
		{
			d.push_back(kmer_list[i].st2);
			via.push_back(d.size()-2);
		}
		else{
			auto it = lower_bound(d.begin(),d.end(),kmer_list[i].st2);
			*it = kmer_list[i].st2;
			via.push_back((int)(it-d.begin()));
		}
	}
}

void make_cigar(vector<data> &res, int K, vector<cigar> &cigar_string,int query_len)
{
	if(res.empty()){return;}
	cigar tmp;
	int match=0,r_match=0,skip=0;
	if(res.size()==0){tmp.len = -1; tmp.type='*'; cigar_string.push_back(tmp);return;}
	if(res[0].st2>0){ tmp.len = res[0].st2; tmp.type='S'; match = tmp.len;cigar_string.push_back(tmp);}
	tmp.len = res[0].ed2-res[0].st2+1;
	if(tmp.len>0){
	match+=tmp.len;
	tmp.type= 'M';
	r_match+=tmp.len;
	cigar_string.push_back(tmp);
	}
	for(int i=1;i<res.size();i++)
	{
		tmp.len=0;
		tmp.type='Z';

		if(res[i].st<=res[i-1].ed)
		{
			tmp.len = res[i-1].ed-res[i].st+1;
			tmp.type='I';
			res[i].st2+=tmp.len;
			res[i].st+=tmp.len;
			match+=tmp.len;
		}
		else if(res[i].st>res[i-1].ed+1)
		{
			tmp.len = res[i].st-(res[i-1].ed+1);
			tmp.type='D';
		}
		if(tmp.len!=0)
		{
			cigar_string.push_back(tmp);
		}
		tmp.len = res[i].ed2-res[i].st2+1;
		tmp.type='M';
		r_match+=tmp.len;
		match+=tmp.len;
		cigar_string.push_back(tmp);
	}
	if(query_len > match)
	{
		tmp.len = query_len-match;
		tmp.type = 'S';
		cigar_string.push_back(tmp);
	}
	if(r_match < 0.8*query_len)
	{
		cigar_string.clear();
		tmp.len = -1;
		tmp.type='*';
		cigar_string.push_back(tmp);
		res.clear();
	}

	return;
}
void chunk_align(vector<data> &lis_res,vector<data> &lis_res2)
{
	data tmp;
	lis_res2.clear();
	tmp.st=-1;tmp.ed=-1;
	tmp.miss=0;
	for(int i=0;i<lis_res.size();i++)
	{
		if(lis_res[i].st2 >= lis_res[i].ed2){continue;}
		if(tmp.st==-1)
		{
			tmp = lis_res[i];
		}
		else if(tmp.ed-lis_res[i].ed == tmp.ed2 - lis_res[i].ed2)// && tmp.miss<4)
		{
//			tmp.miss++;
			tmp.ed = lis_res[i].ed;
			tmp.ed2 = lis_res[i].ed2;
		}
		else if(tmp.st <= lis_res[i].st && lis_res[i].ed <= tmp.ed)
		{
			continue;
		}
		else
		{
			if(tmp.ed-tmp.st+1>=MIN_match&& tmp.st2<tmp.ed2){lis_res2.push_back(tmp);}
			tmp = lis_res[i];
		}
	}
	if(tmp.ed-tmp.st+1>=MIN_match&& tmp.st2<tmp.ed2){
	lis_res2.push_back(tmp);
	}
}

void align_query(int K, unordered_map<long long int,vector<int>> &kmer_ref, read &r,unordered_map<long long int,vector<int>> &kmer_ref2, int flag2)
{
	int i,pre_flag=-1;
	list<data> chunk;
	vector<data> res,kmer_list;
	vector<int> d,via;
	r.cigar_string.clear();
	r.miss=r.ins=r.del=0;
	pre_flag=pre_chunk(r.query,kmer_ref2,r.query_len);
    //pre_flag=-1;
    if(pre_flag!=-1 && flag2==0)
	{
		data tt;
		cigar tt2;
		tt.st = pre_flag;
		tt.st2=0;
		tt.ed2=r.query_len-1;
		tt.ed = pre_flag+r.query_len-1;
		res.push_back(tt);
		tt2.type='M';
		tt2.len=r.query_len;
		r.cigar_string.push_back(tt2);
		skip_read++;
		read_classify(r,res);
	}
	else if(K>0)// &&strcmp(r.query_name+1,"MN00350:100:000H2TV2H:1:13101:8940:17307")==0)
	{
/*		if(strcmp(r.query_name+1,"MN00254:225:000H2JMTW:1:22104:21400:10375")==0){
			printf("asdfasdfasdf");
		}*/
		LIS(kmer_list,r.query,K,kmer_ref,r.query_len,d,via);
		stack<int> st;
		int l=d.size()-2;
		for(i=via.size()-1;i>=0;i--)
		{
			if(l==0)
			{
				st.push(i);
				break;
			}
			if(l==via[i])
			{
				st.push(i);
				l--;
			}
		}
		data tmp;
		tmp.st = -1,tmp.ed=-1;
		vector<data> lis_res,lis_res2;
		while(!st.empty())
		{
			if(tmp.st==-1)
			{
				tmp.st = kmer_list[st.top()].st;
				tmp.ed = kmer_list[st.top()].ed;
				tmp.st2 = kmer_list[st.top()].st2;
				tmp.ed2 = kmer_list[st.top()].ed2;
				tmp.miss=0;
			}
			else
			{
				if(tmp.ed+1==kmer_list[st.top()].ed && tmp.ed2+1 == kmer_list[st.top()].ed2)
				{
					tmp.ed++;
					tmp.ed2++;
				}
				else
				{
					lis_res.push_back(tmp);
					tmp.st=tmp.ed=tmp.st2=tmp.ed2=-1;
				}
			}
			st.pop();
		}
		if(tmp.st!=-1){
		lis_res.push_back(tmp);
		}
		chunk_align(lis_res,lis_res2);
		lis_res = lis_res2;
		for(i=0;i<lis_res2.size();i++){
			if(lis_res2[i].st2<K-3)
			{
				lis_res2[i].st-=lis_res2[i].st2;
				lis_res2[i].st2=0;
			}
			if(r.query_len-lis_res2[i].ed2<K-3)
			{
				lis_res2[i].ed+=r.query_len-1-lis_res2[i].ed2;
				lis_res2[i].ed2=r.query_len-1;
			}
//			printf("%d %d %d %d|| ",lis_res2[i].st,lis_res2[i].ed,lis_res2[i].st2,lis_res2[i].ed2);
		}

		tmp.st2=-1;tmp.ed2=-1;
		for(i=0;i<lis_res2.size();i++)
		{
			if(tmp.ed2>=lis_res2[i].st2)
			{
				int gap = tmp.ed2-lis_res2[i].st2+1;
				lis_res2[i].st+=gap,lis_res2[i].st2+=gap;
			}
			tmp=lis_res2[i];
		}
		chunk_align(lis_res2,lis_res);
		char t[5];
	        t[0] = 'A',t[1] = 'C',t[2] = 'G',t[3] = 'T';
		for(i=0;i<lis_res.size();i++)
		{
			int miss=0;
			while(lis_res[i].st2>0)
			{
				if(t[r.query[lis_res[i].st2-1]] == g_ref_seq[lis_res[i].st-1])
				{
					lis_res[i].st2--;
					lis_res[i].st--;
				}
				else
				{
					if(miss<=0){
					lis_res[i].st2--;
					lis_res[i].st--;
					lis_res[i].miss++;
					miss++;
					}
					else
					{
						break;
					}
				}
			}
			miss=0;
			while(lis_res[i].ed2<r.query_len-1)
			{
				if(t[r.query[lis_res[i].ed2+1]] == g_ref_seq[lis_res[i].ed+1])
				{
					lis_res[i].ed2++;
					lis_res[i].ed++;
				}
				else
				{
					if(miss<=0){
					lis_res[i].ed2++;
					lis_res[i].ed++;
					miss++;

					lis_res[i].miss++;
					}
					else
					{
						break;
					}

				}
			}
		}
		tmp.ed2=-1;
		for(i=0;i<lis_res.size();i++)
		{
			if(tmp.ed2>=lis_res[i].st2)
			{
				int gap = tmp.ed2-lis_res[i].st2+1;
				lis_res[i].st+=gap,lis_res[i].st2+=gap;
			}
			tmp=lis_res[i];
		}
		chunk_align(lis_res,lis_res2);
		lis_res = lis_res2;
		for(i=1;i<lis_res.size();i++){
			 int gap = lis_res[i].st2-lis_res[i-1].ed2-1;
                        if(gap>K+K/2)
                        {
                                lis_res.clear();
                                return;
                        }
                        lis_res[i].st2-=gap;
                        lis_res[i].st-=gap;
		}
		chunk_align(lis_res,lis_res2);
		lis_res = lis_res2;
		int pos=0,miss=0,miss_tg=0,match_tg=0;
		for(i=0;i<lis_res.size();i++)
		{
			pos=0;
			for(int j=lis_res[i].st;j<=lis_res[i].ed;j++)
			{
				if(t[r.query[lis_res[i].st2+pos]]!=g_ref_seq[j])
				{
					miss++;
					if(arg.targetst+arg.cle-5<=j && j<=arg.targetst+arg.cle+5)
					{
						miss_tg++;
					}
				}
				else
				{
					if(arg.targetst+arg.cle-15<=j && j<=arg.targetst+arg.cle+15)
					{
						match_tg++;
					}
				}

				pos++;
			}
		}
		r.tg_miss= miss_tg;
		if(miss_tg>0){r.miss=1;}
		r.tg_match= match_tg;
		r.tags[1]= "MISS:"+(to_string(miss)+":"+to_string(miss_tg));
		//if(miss-miss_tg>=8){lis_res.clear();}//last_modify
		if(miss>=8){lis_res.clear();}//last_modify
		for(i=1;i<lis_res.size();i++)
		{
				while(1)
				{
					if((t[r.query[lis_res[i-1].ed2]] != g_ref_seq[lis_res[i-1].ed]) && (t[r.query[lis_res[i-1].ed2]] == g_ref_seq[lis_res[i].st-1]))
					{
						lis_res[i-1].ed2--;
						lis_res[i-1].ed--;
						lis_res[i].st--;
						lis_res[i].st2--;
					}
					else
					{
						break;
					}
				}
		}
		side_align(lis_res,r);
		make_cigar(lis_res,K,r.cigar_string,r.query_len);
		read_classify(r,lis_res);

	}
	return;
}



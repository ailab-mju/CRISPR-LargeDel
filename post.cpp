#include <vector>
#include <string>
#include <algorithm>
#include <list>
#include <cmath>
#include "post.h"
#include "dkm.hpp"
#include "func.h"
#include "samout.h"

using namespace std;
extern int mut_only,ctr_only,mut_common,ctr_common;
void read_cluster(vector<read> *split_ctr, vector<read> *split_mut,Args args)
{
	bool mutual1[20000],mutual2[20000];
	int count1[20000],count2[20000];
	std::tuple<std::vector<std::array<float, 2>>, std::vector<uint32_t>> cluster_data1 = cluster(split_ctr),cluster_data2 = cluster(split_mut);
	vector<std::array<float,2>> center1 = get<0>(cluster_data1),center2 = get<0>(cluster_data2);
	//printf("%d %d\n",center1.size(),center2.size());
	vector<uint32_t> label1 = get<1>(cluster_data1),label2 = get<1>(cluster_data2);
	for(int i=0;i<center1.size();i++){mutual1[i]=0;count1[i]=0;}
	for(int i=0;i<center2.size();i++){mutual2[i]=0;count2[i]=0;}
	for(int i=0;i<split_ctr[0].size();i++){count1[label1[i]]++;}
	for(int i=0;i<split_mut[0].size();i++){count2[label2[i]]++;}
	

	for(int i=0;i<center1.size();i++)
	{
		for(int j=0;j<center2.size();j++)
		{
			if(count1[i]<=3||count2[j]<=3){continue;}
			float interval = min(abs(center1[i][0]-center1[i][1]),abs(center2[j][0]-center2[j][1]));
			interval*=0.05;
			interval=max((float)interval,(float)50);
			if(abs(center1[i][0]-center2[j][0])<interval && abs(center1[i][1]-center2[j][1])<interval)
			{
				mutual1[i]=1;
				mutual2[j]=1;
			}
		}
	}
	for(int i=0;i<label1.size();i++)
	{
		split_ctr[0][i].ind = mutual1[label1[i]];
		split_ctr[1][i].ind = mutual1[label1[i]];
	}
	for(int i=0;i<label2.size();i++)
	{
		split_mut[0][i].ind = mutual2[label2[i]];
		split_mut[1][i].ind = mutual2[label2[i]];
	}
}
/* FIX: this used to stash its count in the global `mut_only` and the caller
   then copied it into `ctr_only`, which made the control count depend on call
   order and lost the "common" (PCR artefact) count entirely.  Both counts are
   now returned explicitly; the numbers themselves are unchanged. */
void write_read(vector<read> *split, string only,string common,int *n_only,int *n_common)
{
	SamOut *f[2];
	f[0] = samout_open(only);
	f[1] = samout_open(common);
	int c=0;
	for(int i=0;i<split[0].size();i++)
	{
		c+=(!split[0][i].ind);
		split[0][i].pat=split[1][i].pat=(split[0][i].ind?"LDEL_COMMON":"LDEL_ONLY");
		samout_write(f[split[0][i].ind],split[0][i]);
		samout_write(f[split[1][i].ind],split[1][i]);
	}
	if(n_only!=NULL){*n_only=c;}
	if(n_common!=NULL){*n_common=(int)split[0].size()-c;}
	samout_close(f[0]);
	samout_close(f[1]);	
}
void post_process(vector<read> *split_ctr, vector<read> *split_mut,Args args,long n_ctr,long n_mut)
{
	if(args.common_mode==1){read_cluster(split_ctr,split_mut,args);}
	else{local_common(split_ctr,split_mut,n_ctr,n_mut);}
	write_read(split_ctr,args.resctronly,args.resctrcom,&ctr_only,&ctr_common);
	write_read(split_mut,args.resmutonly,args.resmutcom,&mut_only,&mut_common);
}

static const char *k_colour[]={
        "#FF4A46", "#008941", "#006FA6", "#A30059", "#7A4900", "#0000A6", "#63FFAC", "#B79762",
        "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF",
        "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA",
        "#D16100", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F"};

/* ===================================================== local common test
   Which large deletions are "common" (present in the control as well, i.e. a
   PCR / alignment artefact) used to be decided by k-means: both samples were
   clustered on (start, end), k was picked from the grid 1,3,6,9,... by a
   penalised mean error, and two clusters were matched when their CENTROIDS
   were within max(50 bp, 5 % of the span).  Removing 1 % of the reads moved
   the only/common split of 2-9 kb deletions by up to 30 pp, because a slightly
   different point set picks a different k or different seeds, and one centroid
   can swallow real deletions together with an artefact cloud next to them.

   The test is now made per read pair, with no clustering at all:

     tol  = max(20 bp, 2 % of this pair's span)
     n_c  = control pairs with |start diff| < tol and |end diff| < tol
     n_m  = treated pairs in the same window (this pair included)
     common  <=>  n_c > 3                       (same floor as before)
             and  n_c / N_c  >=  0.5 * n_m / N_m

   N_c, N_m are the aligned pairs of each sample, so the second condition asks
   that the control carries at least half the treated sample's rate at these
   breakpoints - an artefact turns up at a similar rate in both, a real CRISPR
   deletion is enriched in the treated sample even with a few artefacts near
   it.  The centroid test used max(50 bp, 5 %); a single pair's breakpoints
   are far more precise than a centroid, and 5 % (+/-250 bp at 5 kb) swept
   up the alignment artefacts a WT control carries near real multi-kb
   deletions (HPRT1: 29 % of real 5 kb deletions removed, 11 % at 2 %).  The
   floor is 20 bp, not 50: at 50 a real 200 bp deletion 30 bp from an
   artefact shared with the control (stress_cluster) was removed with it.
   The result depends only on the set of breakpoints: it is independent
   of read order, of seeds and of k, and changes smoothly with the input.
   The control side is tested the same way against the treated sample.
   --common_mode kmeans restores the previous behaviour. */
#define COMMON_RATIO   0.5
#define COMMON_MIN     4
#define COMMON_TOL_MIN  20
#define COMMON_TOL_FRAC 0.02
static array<float,2> ldel_point(vector<read> *split,size_t i)
{
	const read &t = (split[0][i].map_f==2) ? split[0][i] : split[1][i];
	array<float,2> p={(float)t.ref_st,(float)t.ref_ed};
	return p;
}
static int count_near(const vector<array<float,2>> &sorted,float st,float ed,float tol)
{
	auto lo = lower_bound(sorted.begin(),sorted.end(),array<float,2>{st-tol,-1e9f});
	int n=0;
	for(auto it=lo; it!=sorted.end() && (*it)[0] < st+tol; ++it)
	{
		if((*it)[0] > st-tol && fabs((*it)[1]-ed) < tol){n++;}
	}
	return n;
}
/* deterministic colour groups for the viewer (YC tag): greedy, in sorted order */
static void colour_groups(vector<read> *split)
{
	size_t n=split[0].size();
	vector<pair<array<float,2>,size_t>> v;
	for(size_t i=0;i<n;i++){v.push_back(make_pair(ldel_point(split,i),i));}
	sort(v.begin(),v.end());
	vector<array<float,2>> seeds;
	vector<int> grp(n,0);
	for(auto &e : v)
	{
		float tol=max(50.0f,0.05f*fabs(e.first[1]-e.first[0]));
		int g=-1;
		for(int s=(int)seeds.size()-1; s>=0 && s>=(int)seeds.size()-64; s--)
		{
			if(fabs(seeds[s][0]-e.first[0])<tol && fabs(seeds[s][1]-e.first[1])<tol){g=s;break;}
		}
		if(g<0){seeds.push_back(e.first);g=(int)seeds.size()-1;}
		grp[e.second]=g;
	}
	for(size_t i=0;i<n;i++)
	{
		string c="YC:Z:"+rgb(k_colour[grp[i]%32]);
		split[0][i].tags.push_back(c);
		split[1][i].tags.push_back(c);
	}
}
void local_common(vector<read> *split_ctr, vector<read> *split_mut,long n_ctr,long n_mut)
{
	vector<array<float,2>> pc,pm;
	for(size_t i=0;i<split_ctr[0].size();i++){pc.push_back(ldel_point(split_ctr,i));}
	for(size_t i=0;i<split_mut[0].size();i++){pm.push_back(ldel_point(split_mut,i));}
	vector<array<float,2>> sc=pc, sm=pm;
	sort(sc.begin(),sc.end()); sort(sm.begin(),sm.end());
	double Nc = n_ctr>0 ? (double)n_ctr : 1.0, Nm = n_mut>0 ? (double)n_mut : 1.0;
	for(int side=0; side<2; side++)
	{
		vector<read> *sp      = side ? split_mut : split_ctr;
		vector<array<float,2>> &pts = side ? pm : pc;
		const vector<array<float,2>> &own   = side ? sm : sc;
		const vector<array<float,2>> &other = side ? sc : sm;
		double Nown = side ? Nm : Nc, Noth = side ? Nc : Nm;
		for(size_t i=0;i<pts.size();i++)
		{
			float st=pts[i][0], ed=pts[i][1];
			float tol=max((float)COMMON_TOL_MIN,(float)COMMON_TOL_FRAC*fabs(ed-st));
			int n_o = count_near(other,st,ed,tol);
			int n_s = count_near(own,st,ed,tol);
			bool common;
			if(side){ common = n_o>=COMMON_MIN && (n_o/Noth) >= COMMON_RATIO*(n_s/Nown); }
			else    { /* control pair: shared when the treated sample has it at
			             least at the rate that would make the treated pair common */
			          common = n_s>=COMMON_MIN && n_o>=COMMON_MIN && (n_s/Nown) >= COMMON_RATIO*(n_o/Noth); }
			sp[0][i].ind = common;
			sp[1][i].ind = common;
		}
	}
	colour_groups(split_ctr);
	colour_groups(split_mut);
}

std::tuple<std::vector<std::array<float, 2>>, std::vector<uint32_t>> cluster(vector<read> *split)
{
	if(split[0].size()==0)
	{
		std::tuple<std::vector<std::array<float, 2>>, std::vector<uint32_t>> a;
		return a;
	}
	
	string color[]={
        "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"
	};
	vector<array<float,2>> point;
	int i;
	for(i=0;i<split[0].size();i++)
	{
		read tmp;
		if(split[0][i].map_f==2)
		{
			tmp = split[0][i];
		}
		else
		{
			tmp = split[1][i];
		}
		point.push_back({(float)tmp.ref_st,(float)tmp.ref_ed});
	}
	int cl_k=1,tt=1;
	double min_error=numeric_limits<double>::max(),min_k,reg=2000,error;
	while(cl_k<500)
	{
		error=0;
		if(cl_k>point.size())
		{
			break;
		}
		auto cluster_data = dkm::kmeans_lloyd(point,cl_k);
		auto label = std::get<1>(cluster_data);
		auto center = std::get<0>(cluster_data);
		for(i=0;i<point.size();i++)
		{
			error+=((center[label[i]][0]-point[i][0])*(center[label[i]][0]-point[i][0]))+((center[label[i]][1]-point[i][1])*(center[label[i]][1]-point[i][1]));
		}
		error/=(double)point.size();
		error = error+reg*cl_k;
		if(error<min_error)
		{
			min_error = error;
			min_k = cl_k;
		}
		cl_k=3*tt;
		tt++;
	}
	auto cluster_data = dkm::kmeans_lloyd(point,min_k);
	auto label = std::get<1>(cluster_data);
	for(i=0;i<split[0].size();i++)
	{
		/* k can reach 498 but there are only 128 colours: the old code read
		   past the end of the array */
		split[0][i].tags.push_back("YC:Z:"+rgb(color[label[i]%128]));
		split[1][i].tags.push_back("YC:Z:"+rgb(color[label[i]%128]));
	}
	return cluster_data;
}

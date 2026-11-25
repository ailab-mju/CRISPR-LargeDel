#include <vector>
#include <string>
#include <algorithm>
#include <list>
#include "post.h"
#include "dkm.hpp"
#include "func.h"

using namespace std;
extern int mut_only,ctr_only;
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
void write_read(vector<read> *split, string only,string common)
{
	FILE *f[2];
	f[0] = fopen(only.c_str(),"w");
	f[1] = fopen(common.c_str(),"w");
	int c=0;
	for(int i=0;i<split[0].size();i++)
	{
		c+=(!split[0][i].ind);
		print_read(split[0][i],f[split[0][i].ind]);
		print_read(split[1][i],f[split[1][i].ind]);
	}
	mut_only=c;
	//printf("%d\n",c);
	fclose(f[0]);
	fclose(f[1]);	
}
void post_process(vector<read> *split_ctr, vector<read> *split_mut,Args args)
{
	read_cluster(split_ctr,split_mut,args);
	//printf("ctr only large del size :");
	write_read(split_ctr,args.resctronly,args.resctrcom);
	ctr_only=mut_only;
	//printf("mut only large del size :");
	write_read(split_mut,args.resmutonly,args.resmutcom);
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
		split[0][i].tags.push_back("YC:Z:"+rgb(color[label[i]]));
		split[1][i].tags.push_back("YC:Z:"+rgb(color[label[i]]));
	}
	return cluster_data;
}

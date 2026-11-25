/* FIX: modern libstdc++ no longer pulls these in transitively. */
#include <cstdint>
#include <array>
#include <tuple>
#include <string>
#include <vector>
#include "arg_parse.h"
#include "func.h"

using namespace std;
void read_cluster(vector<read> *split_ctr, vector<read> *split_mut,Args args);
void write_read(vector<read> *split, string only,string common,int *n_only,int *n_common);
void post_process(vector<read> *split_ctr, vector<read> *split_mut,Args args,long n_ctr=0,long n_mut=0);
void local_common(vector<read> *split_ctr, vector<read> *split_mut,long n_ctr,long n_mut);

std::tuple<std::vector<std::array<float, 2>>, std::vector<uint32_t>> cluster(vector<read> *split);

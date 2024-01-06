#include <string>
#include <vector>
#include "arg_parse.h"
#include "func.h"

using namespace std;
void read_cluster(vector<read> *split_ctr, vector<read> *split_mut,Args args);
void write_read(vector<read> *split, string only,string common);
void post_process(vector<read> *split_ctr, vector<read> *split_mut,Args args);

std::tuple<std::vector<std::array<float, 2>>, std::vector<uint32_t>> cluster(vector<read> *split);

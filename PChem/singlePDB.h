#include <iostream>
#include <fstream>

#include "boost/filesystem.hpp"
namespace bfs = boost::filesystem;

struct ParamData
{
	double bind_dist;

	int min_res;
	int max_res;
	double min_perc;
	double max_dist;
	double len_factor;
};
void SinglePDB (bfs::path input, bfs::path, ParamData input_params);


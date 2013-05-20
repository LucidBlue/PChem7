#include <iostream>
#include <fstream>
#include "atomList.h"
#include "singlePDB.h"

void SinglePDB (bfs::path input, bfs::path output, ParamData params)
{
//std::cout << params.max_dist
//	<< " " << params.min_res
//	<< " " << params.max_res
//	<< " " << params.min_perc
//	<< std::endl;

	ProteinComplex testComplex;

	// make this function (and maybe others) part of a separate header and namespace
	// include this namespace so the functions can be called easily

	std::string filename_clean = testComplex.ChangeFilename(input, "_c");
	std::string filename_dist = testComplex.ChangeFilename(input, "_dist");
	std::string filename_res = testComplex.ChangeFilename(input, "_res");

	bfs::path file_clean = output;
	file_clean /= filename_clean;

	bfs::path file_dist = output;
	file_dist /= filename_dist;

	bfs::path file_res = output;
	file_res /= filename_res;

	testComplex.FindDuplicates(input);

	//testComplex.CleanPDB(input, file_clean);
	testComplex.CleanPDB2(input, file_clean);

	//testComplex.LoadPDB(file_clean);
	testComplex.LoadPDB2(file_clean);

	testComplex.RemoveDuplicates();
	
	testComplex.AllAtomsDistCalc(params.bind_dist);		 
	//testComplex.PrintAtomDist(file_dist, params.bind_dist);
	
	testComplex.ExtractResidues();
	testComplex.PrintResidues(file_res);
	
	testComplex.ExtractLoops(params.max_dist, params.min_res, params.max_res, params.min_perc, params.len_factor);
	
	testComplex.PrintLoops(input, output);
}

// 
// atomList.h
// Created by Brad Sheneman
// 10/20/2012
//

#include <string>
#include <vector>

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
namespace bfs = boost::filesystem;

// this needs to become a class
// with a copy constructor
class AtomData
{
public:
	int atom_num;
	int residue_num;
	char chain_id;
	std::string atom_type;
	std::string atom_element;
	std::string residue_type;

	float x_cd;
	float y_cd;
	float z_cd;
	
	float occupancy;
	float temp_factor;
	
	int charge;
	bool interface_atom;
	
	AtomData();
	AtomData(const AtomData &source);
	AtomData& operator= (const AtomData &source);
};

struct AtomPairData
{
	AtomData atom1;
	AtomData atom2;
	
	double distance;
};

struct ResidueData
{
	AtomData aCarbon;
	bool interface_res;
	double distance_to_start;
};

class ProteinComplex
{
public:
	void FindDuplicates(bfs::path filename);
	bool IsDuplicate(char chain_ID_test);
	std::string ChangeFilename(bfs::path input_file, std::string append);
	void CleanPDB(bfs::path input, bfs::path output);
	void CleanPDB2(bfs::path input, bfs::path output);
	void LoadPDB(bfs::path filename);
	void LoadPDB2(bfs::path filename);
	void RemoveDuplicates();
	void InsertAtomData(AtomData atom);
	double AtomDistCalc(AtomData atom1, AtomData atom2);
	void TestCalc();
	void AllAtomsDistCalc(float bind_distance);
	void PrintAtomDist(bfs::path filename, float bind_distance);
	void ExtractResidues();
	void PrintResidues(bfs::path filename);
	void LoopFinder(int chain_ID, double max_dist, int min_res, int max_res, double min_perc, double len_factor);
	void ExtractLoops(double max_dist, int min_res, int max_res, double min_perc, double len_factor);
	void PrintLoops(bfs::path input_file, bfs::path output);
	
private:
	std::vector<std::vector<char> > ChainDuplicates;
	std::vector<std::vector<AtomData> > ComplexAtomData;
	std::vector<std::vector<ResidueData> > ComplexResidues;
	std::vector<std::vector<ResidueData> > ComplexLoops;
};
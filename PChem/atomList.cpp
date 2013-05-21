/*
 *  atomList.cpp
 *  PChemAnalysis
 *
 *  Created by Brad Sheneman
 *  10/20/12
 *
 */

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "atomList.h"
#include <cmath>
#include <locale>

// may need to change the defaults if this runs into errors
AtomData::AtomData()
{
	atom_num = 0;
	residue_num = 0;
	chain_id = ' ';
	atom_type = "";
	atom_element = "";
	residue_type = "";
	
	x_cd = 0;
	y_cd = 0;
	z_cd = 0;
	
	occupancy = 0;
	temp_factor = 0;
	
	charge = 0;
	interface_atom = false; 
}

AtomData::AtomData(const AtomData &source)
{
	atom_num = source.atom_num;
	residue_num = source.residue_num;
	chain_id = source.chain_id;
	atom_type = source.atom_type;
	atom_element = source.atom_element;
	residue_type = source.residue_type;
	
	x_cd = source.x_cd;
	y_cd = source.y_cd;
	z_cd = source.z_cd;
	
	occupancy = source.occupancy;
	temp_factor = source.temp_factor;
	
	charge = source.charge;
	interface_atom = source.interface_atom;
}

AtomData& AtomData::operator= (const AtomData &source)
{
	atom_num = source.atom_num;
	residue_num = source.residue_num;
	chain_id = source.chain_id;
	atom_type = source.atom_type;
	atom_element = source.atom_element;
	residue_type = source.residue_type;
	
	x_cd = source.x_cd;
	y_cd = source.y_cd;
	z_cd = source.z_cd;
	
	occupancy = source.occupancy;
	temp_factor = source.temp_factor;
	
	charge = source.charge;
	interface_atom = source.interface_atom;
	
	return *this;
}

// can only run on an un-edited PDB file (requires some specific header lines)
void ProteinComplex::FindDuplicates(bfs::path filename)
{
	bfs::ifstream infile(filename);
	std::string current_line;
	
	if(!infile.is_open())
	{
		std::cout << "File not found\n";
		return;
	}
	
	for (int chain_iter = 0; !infile.eof(); chain_iter++)
	{
 		std::getline(infile, current_line);
		
		if (current_line.find("SOURCE", 0) == 0)
			break;
		
		if (current_line.find("COMPND", 0) == 0)
		{
			if (current_line.find("CHAIN", 0) < 12)
			{
				std::vector<char> TempDuplicates;
				
				for (int iter = current_line.find("CHAIN", 0) + 6; iter < current_line.size(); iter++)
				{
					if (std::isalpha(current_line[iter], std::locale()))
						TempDuplicates.push_back(current_line[iter]);
				}
				
				ChainDuplicates.push_back(TempDuplicates);
			}
			
		}
	}
	infile.close();
}

std::string ProteinComplex::ChangeFilename(bfs::path input_file, std::string append)
{

	std::string file_name = input_file.filename().string();
	std::string file_name_new;

	if (file_name.find(".pdb", 0) || file_name.find(".PDB", 0))
			file_name.erase(file_name.end() - 4, file_name.end());
	
	file_name_new = file_name + append + ".pdb";

	return file_name_new;
}
/*
void ProteinComplex::CleanPDB(bfs::path input, bfs::path output)
{
	bfs::ifstream infile(input);
	bfs::ofstream outfile(output);
	
	if(!infile.is_open())
	{
		std::cout << "File not found\n";
		return;
	}
	
	for (int i = 0; !infile.eof(); i++)
	{
		std::string current_line;
		std::getline(infile, current_line);
	
		if (current_line.find("ATOM", 0) == 0)
			outfile << current_line << std::endl;
	}
	
	infile.close();
	outfile.close();
}
*/
bool ProteinComplex::IsDuplicate(char chain_ID_test)
{
	bool found = false;

	for (int i = 0; i < ChainDuplicates.size(); i++)
	{
		for (int j = 1; j < ChainDuplicates[i].size(); j++)
		{
			if (chain_ID_test == ChainDuplicates[i][j])
			{
				found = true;
			}
		}
	}
	return found;
}

void ProteinComplex::CleanPDB2(bfs::path input, bfs::path output)
{
	bfs::ifstream infile(input);
	bfs::ofstream outfile(output);
	
	std::vector<char> Visited_Chains;
	

	if(!infile.is_open())
	{
		std::cout << "File not found\n";
		return;
	}

	char current_chain = '*';

	for (int i = 0; !infile.eof(); i++)
	{
		std::string current_line;
		std::getline(infile, current_line);
		

		if (current_line.find("ATOM", 0) == 0)
		{
			bool visited = false;

			for (int a = 0; a < Visited_Chains.size(); a++)
			{
				if (current_line[21] == Visited_Chains[a])
					visited = true;
			}

			if (current_chain != current_line[21] && current_chain != '*')
				Visited_Chains.push_back(current_chain);

			if (visited == false && !IsDuplicate(current_line[21]))
				outfile << current_line << std::endl;

			current_chain = current_line[21];
		}
	}

	infile.close();
	outfile.close();
}

// as of now, loads all atom data from a cleaned PDB file
// will deal with how to add built-in cleaning method to this function
/*void ProteinComplex::LoadPDB(bfs::path filename)
{
	bfs::ifstream infile(filename);
	std::string carry_over;

	if(!infile.is_open())
	{
		std::cout << "File not found\n";
		return;
	}
	
	while (!infile.eof())
	{
		infile >> carry_over;
		if (!infile.eof())
		{
			if(carry_over == "ATOM")
				break;
		}
	}
	
	while (!infile.eof()) 
	{	
		if (carry_over == "ATOM")
		{
			AtomData current_atom;
			infile >> current_atom.atom_num;
			infile >> current_atom.atom_type;
			infile >> current_atom.residue_type;
			infile >> current_atom.chain_id;
			infile >> current_atom.residue_num;
			infile >> current_atom.x_cd;
			infile >> current_atom.y_cd;
			infile >> current_atom.z_cd;
			infile >> current_atom.occupancy;
			infile >> current_atom.temp_factor;
			infile >> current_atom.atom_element;
		
			InsertAtomData(current_atom);
			infile >> carry_over;
		}
		else
			infile >> carry_over;
	}
	
	infile.close();
}*/

void ProteinComplex::LoadPDB2(bfs::path filename)
{
	bfs::ifstream infile(filename);
	std::string current_line;

	if(!infile.is_open())
	{
		std::cout << "File not found\n";
		return;
	}
	
	while (!infile.eof())
	{
		std::getline(infile, current_line);
		if (!infile.eof())
		{
			if(current_line.substr(0, 4) == "ATOM")
				break;
		}
	}

	while (!infile.eof()) 
	{	
		if(current_line.substr(0, 4) == "ATOM")
		{
			AtomData current_atom;

			std::istringstream current_val(current_line.substr(6, 5));
			current_val >> current_atom.atom_num;
			current_val.clear();

			current_val.str(current_line.substr(12, 4));
			current_val >> current_atom.atom_type;
			current_val.clear();

			current_val.str(current_line.substr(17, 3));
			current_val >> current_atom.residue_type;
			current_val.clear();

			current_val.str(current_line.substr(21, 1));
			current_val >> current_atom.chain_id;
			current_val.clear();

			current_val.str(current_line.substr(22, 4));
			current_val >> current_atom.residue_num;
			current_val.clear();

			current_val.str(current_line.substr(30, 8));
			current_val >> current_atom.x_cd;
			current_val.clear();

			current_val.str(current_line.substr(38, 8));
			current_val >> current_atom.y_cd;
			current_val.clear();

			current_val.str(current_line.substr(46, 8));
			current_val >> current_atom.z_cd;
			current_val.clear();

			current_val.str(current_line.substr(54, 6));
			current_val >> current_atom.occupancy;
			current_val.clear();

			current_val.str(current_line.substr(60, 6));
			current_val >> current_atom.temp_factor;
			current_val.clear();

			current_val.str(current_line.substr(76, 2));
			current_val >> current_atom.atom_element;
			current_val.clear();
		
			InsertAtomData(current_atom);
			std::getline(infile, current_line);
		}
		else
			std::getline(infile, current_line);
	}
	infile.close();
}

// can run only once PDB file is loaded into ComplexAtomData
void ProteinComplex::RemoveDuplicates()
{	
	for (int i = 0; i < ChainDuplicates.size(); i++)
	{
		for (int j = 1; j < ChainDuplicates[i].size(); j++)
		{
			std::vector<std::vector<AtomData> >::iterator atom_iter = ComplexAtomData.begin();
			while (atom_iter < ComplexAtomData.end())
			{
				if ((*atom_iter)[0].chain_id == ChainDuplicates[i][j])
					atom_iter = ComplexAtomData.erase(atom_iter);
				else
					atom_iter++;
			}
		}
	}
}

void ProteinComplex::InsertAtomData(AtomData atom)
{
	char current_chain = ' ';
	
	for (int i = 0; i < ComplexAtomData.size(); i++)
	{	
		current_chain = ComplexAtomData[i][0].chain_id;
		if (current_chain == atom.chain_id)
		{
			ComplexAtomData[i].push_back(atom);
			break;
		}
	}
	
	if (current_chain != atom.chain_id)
	{
		std::vector<AtomData> Current_Chain;
		Current_Chain.push_back(atom);
		ComplexAtomData.push_back(Current_Chain);
	}
}

// calculates distances between a pair of AtomData objects and returns value as double
double ProteinComplex::AtomDistCalc(AtomData atom1, AtomData atom2)
{
	float x_1, x_2, y_1, y_2, z_1, z_2;
	x_1 = atom1.x_cd;
	x_2 = atom2.x_cd;
	y_1 = atom1.y_cd;
	y_2 = atom2.y_cd;
	z_1 = atom1.z_cd;
	z_2 = atom2.z_cd;
	
	return (sqrt(pow((x_2-x_1),2)+pow((y_2-y_1),2)+pow((z_2-z_1),2)));
}

// calls distance calc function for all atom pairs on separate chains
// and stores them in a vector list of atom pairs

// write function that performs this nested loop
// this traverses the vector of vectors, visiting each item once
// return AtomData?

void ProteinComplex::TestCalc()
{
	std::cout<< ComplexAtomData[0][0].atom_num << " ";
	std::cout<< ComplexAtomData[0][6].atom_num << " ";
	std::cout<< ComplexAtomData[0][14].atom_num << " ";
	std::cout<< ComplexAtomData[0][18].atom_num << " " << std::endl;
	std::cout<< AtomDistCalc(ComplexAtomData[0][0], ComplexAtomData[0][6]) << std::endl;
	std::cout<< AtomDistCalc(ComplexAtomData[0][6], ComplexAtomData[0][14]) << std::endl;
	std::cout<< AtomDistCalc(ComplexAtomData[0][14], ComplexAtomData[0][18]) << std::endl;
}

void ProteinComplex::AllAtomsDistCalc(float bind_distance)
{
	for (int i = 0; i < ComplexAtomData.size(); i++)
	{
		for (int j = 0; j < ComplexAtomData[i].size(); j++)
		{
			for (int a = i + 1; a < ComplexAtomData.size(); a++)
			{
				for (int b = 0; b < ComplexAtomData[a].size(); b++)
				{
					AtomPairData data;
					
					data.distance = AtomDistCalc(ComplexAtomData[i][j], ComplexAtomData[a][b]);
					if (data.distance < bind_distance)
					{
						ComplexAtomData[i][j].interface_atom = true;
						ComplexAtomData[a][b].interface_atom = true;
					}
					data.atom1 = ComplexAtomData[i][j];
					data.atom2 = ComplexAtomData[a][b];
				}
			}
		}
	}
}

// prints out all atom pair distance data into a file
// output formatting needs some work

/*
void ProteinComplex::PrintAtomDist(bfs::path filename, float bind_distance)
{
	bfs::ofstream outfile(filename);
	
	for (int i = 0; i < ComplexDist.size(); i++)
	{
		// redundant test to be sure both atoms being printed are actually interface
		if (ComplexDist[i].distance < bind_distance)
		{
			outfile << std::left << std::setw(10) << ComplexDist[i].atom1.atom_num;
			outfile << std::left << std::setw(10) << ComplexDist[i].atom1.atom_type;
			outfile << std::left << std::setw(10) << ComplexDist[i].atom1.chain_id;
			outfile << std::left << std::setw(10) << ComplexDist[i].atom1.residue_num;
			outfile << std::left << std::setw(10) << ComplexDist[i].atom1.interface_atom;
			outfile << std::left << std::setw(10) << ComplexDist[i].atom2.atom_num;
			outfile << std::left << std::setw(10) << ComplexDist[i].atom2.atom_type;
			outfile << std::left << std::setw(10) << ComplexDist[i].atom2.chain_id;
			outfile << std::left << std::setw(10) << ComplexDist[i].atom2.residue_num;
			outfile << std::left << std::setw(10) << ComplexDist[i].atom2.interface_atom;
			outfile << ComplexDist[i].distance << std::endl;
		}
	}
}
*/

// extracts all the residues and marks them as interface or not
// inserted with data corresponding to the a-carbon of that residue
void ProteinComplex::ExtractResidues()
{
	for (int i = 0; i < ComplexAtomData.size(); i++)
	{	
		std::vector<ResidueData> Current_Chain;
		AtomData test_res = ComplexAtomData[i][0];
		
		bool interface = false;
		AtomData current_ACarbon;
		
		for (int j = 0; j < ComplexAtomData[i].size(); j++)
		{
			AtomData current_res = ComplexAtomData[i][j];
			
			if (current_res.residue_num == test_res.residue_num)
			{
				if (current_res.atom_type == "CA")
					current_ACarbon = ComplexAtomData[i][j];
				// if any of the atoms for this residue are interface atoms
				// set the local bool for this residue to true
				if (current_res.interface_atom == true)
					interface = true;
			}
			
			else
			{
				// otherwise, we've moved on to the next residue
				// and we have to insert the data for the last residue
				ResidueData ResData;
				ResData.aCarbon = current_ACarbon;
				ResData.interface_res = interface;
				Current_Chain.push_back(ResData);
				
				// reset test res to the new residue;
				test_res = current_res;
				interface = false;
				AtomData blank_atom;
				current_ACarbon = blank_atom;
				
				// and inner loop needs to move back one? (in order to revisit this atom)
				j--;
			}

		}
		ComplexResidues.push_back(Current_Chain);
	}
}

void ProteinComplex::PrintResidues(bfs::path filename)
{
	bfs::ofstream outfile(filename);
	
	for (int i = 0; i < ComplexResidues.size(); i++)
	{
		for (int j = 0; j < ComplexResidues[i].size(); j++)
		{
			outfile << std::left << std::setw(10) << ComplexResidues[i][j].aCarbon.chain_id;
			outfile << std::left << std::setw(10) << ComplexResidues[i][j].aCarbon.residue_num;
			outfile << std::left << std::setw(10) << ComplexResidues[i][j].aCarbon.residue_type;
			outfile << std::left << std::setw(10) << ComplexResidues[i][j].interface_res << std::endl;
		}
	}
}

void ProteinComplex::LoopFinder(int chain_ID, double max_dist, int min_res, int max_res, double min_perc, double len_factor)
{
	int chain_size = ComplexResidues[chain_ID].size();
	
	// visit each residue in the chain
	for (int i = 0; i < chain_size; i++)
	{
		// for each residue, restart the count of total and interface residues
		std::vector<ResidueData> TempLoop;
		int num_res_interface = 0;
		int num_res_total = 0;
		int current_loop_length = 0;
		int first_res_num = ComplexResidues[chain_ID][i].aCarbon.residue_num;

		// prev_res_num and current_res_num are counting controls to make sure gaps in the protein chain are accounted for
		// start counting loops starting at that residue 
		// and ending at the maximum allowed loop size (or the end of the chain)
		for (int j = i; (current_loop_length < max_res) && (j < chain_size); j++)
		{
			double distance;
			double percent_interface;
			
			ResidueData current_res = ComplexResidues[chain_ID][j];
			int current_res_num = current_res.aCarbon.residue_num;
			current_loop_length = current_res_num - first_res_num + 1;

			TempLoop.push_back(current_res);
			num_res_total++;

			// average length of an amino acid is about 3.2 angstroms (this can be modified if necessary)
			// linker length shouldn't be more than about half of the loop length
			double max_linker_length = len_factor*(TempLoop.size());
			
			if (ComplexResidues[chain_ID][j].interface_res)
				num_res_interface++;
			// if the loop is large enough, measure the distance end to end
			if (j - i >= min_res - 1)
			{
				distance = AtomDistCalc(TempLoop.front().aCarbon, TempLoop.back().aCarbon);
				percent_interface = (double(num_res_interface)/double(num_res_total))*100;
				//std::cout << "nrt: " << num_res_total << " nri: " << num_res_interface << " pi: " << percent_interface << std::endl;
				

				// if distance and % interface criteria meet, add to 'loops' vector
				if (distance <= max_dist && distance <= max_linker_length && percent_interface >= min_perc)
				{
					TempLoop.back().distance_to_start = distance;
					//std::cout << percent_interface << "% " << distance << std::endl;
					ComplexLoops.push_back(TempLoop);
				}
			}

		}
	}
}

void ProteinComplex::ExtractLoops(double max_dist, int min_res, int max_res, double min_perc, double len_factor)
{
	for (int i = 0; i < ComplexResidues.size(); i++)
		LoopFinder(i, max_dist, min_res, max_res, min_perc, len_factor);
}

void ProteinComplex::PrintLoops(bfs::path input_file, bfs::path output)
{
	bfs::path output_file = output;
	output_file /= "Loops.txt";

	bfs::ofstream outfile(output_file, bfs::ofstream::app);

	outfile << input_file.filename().string() << std::endl << std::endl;
	
	for (int i = 0; i < ComplexLoops.size(); i++)
	{
		outfile << " linker distance: " << std::left << std::setw(9) << ComplexLoops[i].back().distance_to_start << " ";
		for (int j = 0; j < ComplexLoops[i].size(); j++)
		{
			ResidueData Current_Res = ComplexLoops[i][j];
			outfile << Current_Res.aCarbon.chain_id << " " 
			<< Current_Res.aCarbon.residue_num
			<< " " << Current_Res.aCarbon.residue_type;
			if (j < (ComplexLoops[i].size() -1))
				outfile << ", ";
		}
		
		outfile << std::endl;
	}
	outfile << std::endl;
	outfile.close();
}
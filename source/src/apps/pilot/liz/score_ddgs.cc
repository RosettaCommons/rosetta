// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Liz Kellogg ekellogg@u.washington.edu

// libRosetta headers

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/sasa.hh>
#include <core/io/raw_data/ScoreMap.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <devel/init.hh>

#include <basic/options/option_macros.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <protocols/ddg/ddGData.hh>
#include <utility/file/FileName.hh>
#include <utility/excn/Exceptions.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <basic/database/open.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>

#include <basic/Tracer.hh>

using namespace core;
using namespace scoring;
using namespace protocols;
using namespace protocols::ddG;
using namespace ObjexxFCL;

int min_index(utility::vector1<double>scores){
	double min=999;
	int min_index=-1;
	for(unsigned int i=1;i<=scores.size();i++){
		if(scores[i] < min) {
			min_index=i;
			min = scores[i];
		}
	}
	return min_index;
}

utility::vector1<std::string>
header(pose::Pose p,ScoreFunctionOP sfxn){
	std::map<std::string,core::Real> sm;
	(*sfxn)(p);
	/**
	ScoreMap::score_map_from_scored_pose(sm,p);
	std::map< std::string, Real >::const_iterator pair;
	Size width (8);
	for ( pair=sm.begin(); pair!=sm.end(); ++pair )
	{
		if ( pair->first.length() > 8 ) width = pair->first.length();
		out << ' ' << A( width, pair->first );
	}
	out << std::endl;
	**/
	utility::vector1<std::string> nonzero_weights;
	core::scoring::EnergyMap weights = p.energies().weights();
	nonzero_weights.push_back(name_from_score_type(core::scoring::total_score));
	for(unsigned int i = 1; i <= core::scoring::n_score_types;++i){
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);
		if(weights[ii] != 0) nonzero_weights.push_back(name_from_score_type(ii));
	}
	return nonzero_weights;
}

utility::vector1<double>
weights(pose::Pose p, ScoreFunctionOP sfxn){
	(*sfxn)(p);
	utility::vector1<double> nonzero_weights;
	core::scoring::EnergyMap weights = p.energies().weights();
	nonzero_weights.push_back(1); //total score gets a weight of 1
	for(unsigned int i = 1; i <= core::scoring::n_score_types;++i){
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);
		if(weights[ii] != 0) nonzero_weights.push_back(weights[ii]);
	}
	return nonzero_weights;
}

bool
sort_numerically_ascending(double a, double b){
	return a < b;
}

utility::vector1<double>
calc_ddg(utility::vector1<pose::Pose> mut, utility::vector1<pose::Pose> wt,
		ScoreFunctionOP sfxn,bool mean,bool min){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::ddg;

	utility::vector1<double> ddg;

	/*
		modification added in to take not the average of all structure energies
		but only the lowest 3 scoring structures.
		bool average_lowest_x turns this on and off
	*/
	bool average_lowest_x = false;
	int input_avg_lowest = 0;
	core::Size avg_lowest = 0;
	if(basic::options::option[OptionKeys::ddg::lowest_x_decoys].user()){
		average_lowest_x = true;
		input_avg_lowest = basic::options::option[OptionKeys::ddg::lowest_x_decoys]();
		if( input_avg_lowest  <= 0){
			std::cerr << "ERROR you need to specify a value > 0 for OPTION: ddg::lowest_x_decoys\nsetting avg_lowest to arbitrary value of 10" << std::endl;
			avg_lowest = 10;
		} else {
			avg_lowest = static_cast< Size > ( input_avg_lowest );
		}
	}

	utility::vector1<double> scores_mut;
	utility::vector1<double> scores_wt;


	for(core::Size i=1;i <= mut.size();i++){
		double sm = (*sfxn)(mut[i]);
		if(average_lowest_x)scores_mut.push_back(sm);
	}
	for(core::Size i=1;i <= wt.size();i++){
		double sw = (*sfxn)(wt[i]);
		if(average_lowest_x)scores_wt.push_back(sw);
	}

	utility::vector1<std::string> h = header(mut[1],sfxn);
	//	assert(env_wts.size() == h.size());
	if(mean){
		if(average_lowest_x){
			//get scores of all structures
			sort(scores_mut.begin(),scores_mut.end(),sort_numerically_ascending);
			sort(scores_wt.begin(),scores_wt.end(),sort_numerically_ascending);

			//std::cout << "checking output" << std::endl;
			//for(unsigned int i = 1; i <= scores_mut.size(); i++){
			//	std::cout << scores_mut[i] << " ";
			//}
			//std::cout << std::endl;
			//for(unsigned int i = 1; i <= scores_wt.size(); i++){
			//	std::cout << scores_wt[i] << " ";
			//}
			//std::cout << std::endl;

			//get lowest x and remove others from mut and wt vectors
			utility::vector1<pose::Pose> new_muts;
			utility::vector1<pose::Pose> new_wts;

			if(avg_lowest > mut.size() || avg_lowest > wt.size()){
				std::cerr << "cannot take lowest: " << avg_lowest << " decoys because there are only " << mut.size() << " " << wt.size() << " poses! " << std::endl;
				if( mut.size() < wt.size()){ avg_lowest = mut.size();}
				else{ avg_lowest = wt.size();}
			}
			for(unsigned int index = 1; index <= avg_lowest; index++){
				for(unsigned int i = 1;i <= mut.size(); i++){
					double s = (*sfxn)(mut[i]);
					if(s == scores_mut[index]){
						new_muts.push_back(mut[i]);
						break;
					}
				}
			}
			for(unsigned int index = 1; index <= avg_lowest; index++){
				for(unsigned int j = 1; j <= wt.size();j++){
					double s = (*sfxn)(wt[j]);
					if(s == scores_wt[index]){
						new_wts.push_back(wt[j]);
						break;
					}
				}
			}
			//replace old mut and wt with the new ones
			mut = new_muts;
			wt = new_wts;
		}
		//	std::cout << "dbl-checking muts size " << mut.size() << " and wt size " << wt.size() << std::endl;
		ddg.resize(h.size(),0.0);
		for(unsigned int ii = 1;ii <= h.size(); ii++){
			double avg_mut = 0.0;
			double avg_wt  = 0.0;

			for(unsigned int i=1; i<= mut.size(); i++){
				std::map<std::string, core::Real> scoremap_mut;
				ScoreMap::nonzero_energies(scoremap_mut,sfxn,mut[i]);
				avg_mut+=((scoremap_mut[h[ii]])*((double)1/(double)mut.size()));
				//std::cout << " avg mut sc is " << (scoremap_mut[h[ii]]) << " and modulated by " << ((double)1/(double)mut.size()) << std::endl;
			}

			for(unsigned int i=1; i<= wt.size(); i++){
				std::map<std::string, core::Real> scoremap_wt;
				ScoreMap::nonzero_energies(scoremap_wt,sfxn,wt[i]);
				avg_wt+=((scoremap_wt[h[ii]]))*((double)1/(double)wt.size());
			}
			//std::cout << "final check: " << avg_mut << " " << avg_wt << std::endl;
			ddg[ii] = (avg_mut - avg_wt);
		}
	}
	else if(min){
		std::map<std::string, core::Real> scoremap_mut;
		std::map<std::string, core::Real> scoremap_wt;
		utility::vector1<double>mut_scores;
		utility::vector1<double>wt_scores;
		for(unsigned int i=1;i <= mut.size();i++){
			mut_scores.push_back((*sfxn)(mut[i]));
		}
		for(unsigned int i=1;i <= wt.size();i++){
			wt_scores.push_back((*sfxn)(wt[i]));
		}
		pose::Pose min_mut = mut[min_index(mut_scores)];
		pose::Pose min_wt = wt[min_index(wt_scores)];

		ScoreMap::nonzero_energies(scoremap_mut,sfxn,min_mut);
		ScoreMap::nonzero_energies(scoremap_wt,sfxn,min_wt);
		ddg.resize(scoremap_mut.size(),0.0);
		for(unsigned int ii = 1;ii <= h.size(); ii++){
			ddg[ii]=((scoremap_mut[h[ii]])-(scoremap_wt[h[ii]]));
		}
	}
	return ddg;
}


std::string
mut_info(pose::Pose m,pose::Pose w){
	std::string info = "";
	assert(m.size() == w.size());
	for(unsigned int ii = 1;ii <= m.size(); ii++){
		if(m.residue(ii).name1() != w.residue(ii).name1()){
			std::ostringstream convert_to_string;
			convert_to_string << ii;
			info= w.residue(ii).name() + " "
				+ convert_to_string.str() + " "
				+ m.residue(ii).name();
			return info;
		}
	}
	return info;
}


void
print_header(pose::Pose p,ScoreFunctionOP sfxn,
						 std::ostream & out){
	utility::vector1<std::string> h = header(p,sfxn);
	for(unsigned int j = 1; j<= h.size(); j++){
		out << h[j] << " ";
	}
	out << std::endl;
}

void
print_ddgs(utility::vector1<double> ddg,
					 utility::vector1<pose::Pose> mut,
					 utility::vector1<pose::Pose> wt,
					 core::Real experimental_value,
					 std::ostream & out){
	//check to see that all sequences match across ensembles
	std::string sequence_check="";
	using namespace ObjexxFCL::format;
	//mjo commenting out 'width' and 'precision' because they are unused and cause warnings
	//Size width (8), precision (3);
	for(unsigned int i =1; i <= mut.size(); i++){
		if(i == 1){
			sequence_check = mut[i].sequence();
		}else{
			if((mut[i].sequence()).compare(sequence_check) != 0){
				std::cerr << "[ERROR] nth sequence " << i << " " << mut[i].sequence() << " doesn't match first sequence: " << sequence_check << std::endl;
				exit(1);
			}
		}
	}
	for(unsigned int j =1; j <= wt.size(); j++){
		if(j == 1){
			sequence_check = wt[j].sequence();
		}else{
			if((wt[j].sequence()).compare(sequence_check) != 0){
				std::cerr << "[ERROR] nth sequence " << j << " " << wt[j].sequence() << " doesn't match first sequence: " << sequence_check << std::endl;
				exit(1);
			}
		}
	}
	std::string curr_mut_info = mut_info(mut[1],wt[1]);
	out << curr_mut_info << " " ;
	for(unsigned int ii = 1; ii <= ddg.size(); ii++){
		//out << F(width, precision, ddg[ii]) << ' ';
		out << ddg[ii] << ' ';
	}
	out << ' ' << experimental_value << std::endl;
}

//calcs ddgs and outputs.
void
print_verbose_ddgs(utility::vector1<pose::Pose> mut,
									 utility::vector1<pose::Pose> wt,
									 ScoreFunctionOP sfxn,
									 bool mean,
									 bool min,
									 core::Real experimental_value,
									 std::ostream & outfile){
	//check to see that all sequences match across ensembles

	core::Size width (8), precision (3);
	{ //sequence check
	std::string sequence_check="";
	for(unsigned int i =1; i <= mut.size(); i++){
		if(i == 1){
			sequence_check = mut[i].sequence();
		}else{
			if((mut[i].sequence()).compare(sequence_check) != 0){
				std::cerr << "[ERROR] nth sequence " << i << " " << mut[i].sequence() << " doesn't match first sequence: " << sequence_check << std::endl;
				exit(1);
			}
		}
	}
	for(unsigned int j =1; j <= wt.size(); j++){
		if(j == 1){
			sequence_check = wt[j].sequence();
		}else{
			if((wt[j].sequence()).compare(sequence_check) != 0){
				std::cerr << "[ERROR] nth sequence " << j << " " << wt[j].sequence() << " doesn't match first sequence: " << sequence_check << std::endl;
				exit(1);
			}
		}
	}
	} //end sequence check
	std::string curr_mut_info = mut_info(mut[1],wt[1]);

	//non-zero weights
	utility::vector1<std::string> nonzero_weights = header(mut[1],sfxn);
	utility::vector1<double> wts = weights(mut[1],sfxn);
	assert(wts.size() == nonzero_weights.size());

	//output for each residue:
	//mutation info + residue number + energies

	utility::vector1<double> avg_res_total_E_diff(mut[1].size(),0.0);

	if(mean){
		//compute mean residue difference across all poses and residues
		FArray2D<double> avg_score_components(mut[1].size(),nonzero_weights.size(),0.0);

		for(unsigned int ii = 1; ii <= mut.size(); ii++){ //iterate over mutant poses
			//mut_i.accumulate_residue_total_energies();
			(*sfxn)(mut[ii]); //in order to update energies
			for(unsigned int resn = 1; resn <= mut[1].size(); resn++){
				EnergyMap mut_i_resn = mut[ii].energies().residue_total_energies(resn);
				avg_res_total_E_diff[resn]+=mut[ii].energies().residue_total_energy(resn)*((double)1/(double)mut.size());
				//		std::cout << "DEBUG: " << mut[ii].energies().residue_total_energy(resn)*((double)1/(double)mut.size()) << std::endl;
				for(unsigned int s = 1; s <= nonzero_weights.size(); s++){
					avg_score_components(resn,s)+=(mut_i_resn[core::scoring::score_type_from_name(nonzero_weights[s])]*wts[s]*((double)1/(double)mut.size()));
				}
			}
		}
		for(unsigned int ii = 1; ii <= wt.size(); ii++){//iterate over wt poses
			//wt_i.accumulate_residue_total_energies();
			(*sfxn)(wt[ii]);
			for(unsigned int resn = 1; resn <= wt[1].size(); resn++){
				EnergyMap wt_i_resn = wt[ii].energies().residue_total_energies(resn);
				avg_res_total_E_diff[resn]-=(wt[ii].energies().residue_total_energy(resn))*((double)1/(double)wt.size());
				for(unsigned int s = 1; s <= nonzero_weights.size(); s++){
					avg_score_components(resn,s)-=(wt_i_resn[core::scoring::score_type_from_name(nonzero_weights[s])]*wts[s]*((double)1/(double)wt.size()));
				}
			}
		}

		//output score components
		for( int i = 1; i <= avg_score_components.u1(); i++){
			outfile << curr_mut_info << ' ' << i << ' ' /**<< avg_res_total_E_diff[i] << ' ' **/;
			for( int j = 1; j <= avg_score_components.u2(); j++){
				outfile << ObjexxFCL::format::F(width,precision,avg_score_components(i,j)) << ' ';
			}
			outfile << '\n';
		}
	}
	else if(min){
		//compute min residue difference across all residues
		utility::vector1<double> mut_scores, wt_scores;
		for(unsigned int ii = 1; ii <= mut.size(); ii++){
			mut_scores.push_back((*sfxn)(mut[ii]));
		}
		for(unsigned int jj = 1; jj <= wt.size(); jj++){
			wt_scores.push_back((*sfxn)(wt[jj]));
		}
		//find lowest energy mut and wt
		int min_mut_ind = min_index(mut_scores);
		int min_wt_ind = min_index(wt_scores);

		//		std::cout << "min mut index " << min_mut_ind << " wt " << min_wt_ind << std::endl;

		//compute per-residue diff between mut and wt min energy structures
		//output information
		Energies me = mut[min_mut_ind].energies();
		Energies we = wt[min_wt_ind].energies();

		//total energy with hbond_lr_bb and hbond_sr_bb
		for(unsigned int resn = 1; resn <= mut[1].size(); resn++){
			EnergyMap m(mut[min_mut_ind].energies().residue_total_energies(resn));
			EnergyMap w(wt[min_wt_ind].energies().residue_total_energies(resn));
			outfile <<curr_mut_info << ' ' << resn << ' ';
			for(unsigned int s = 1; s <= nonzero_weights.size(); s++){
				outfile << ObjexxFCL::format::F(width, precision, (m[core::scoring::score_type_from_name(nonzero_weights[s])] -
																												w[core::scoring::score_type_from_name(nonzero_weights[s])])*wts[s]) << ' ' ;
			}
			outfile << '\n';
		}
	}else{

	}
	outfile << ' ' << experimental_value << std::endl;
}


utility::vector1<int>
num_nbrs(pose::Pose p){
	utility::vector1<int> nbrs;
	core::scoring::TenANeighborGraph const tenA_neighbor_graph( p.energies().tenA_neighbor_graph() );

	nbrs.resize( p.size() );
	for ( Size i=1; i<= p.size(); ++i ) {
		{
			nbrs[i] = 1;
			for ( utility::graph::Graph::EdgeListConstIter
							ir  = tenA_neighbor_graph.get_node( i )->const_edge_list_begin(),
							ire = tenA_neighbor_graph.get_node( i )->const_edge_list_end();
						ir != ire; ++ir ) {
				core::Size const neighbor_id( (*ir)->get_other_ind( i ) );
				chemical::ResidueType const & nbr_rsd( p.residue_type( neighbor_id ) );
				if ( nbr_rsd.is_protein() ) {
					nbrs[i] += 1;
				}
			}
		}
	}
	return nbrs;
}

/**
utility::vector1<core::Real>
sasa(pose::Pose p){
	const core::Real probe_radius = 1.2;
	id::AtomID_Map< Real > atom_sasa;
	utility::vector1< core::Real > rsd_sasa;
	core::scoring::calc_per_atom_sasa(p,atom_sasa,rsd_sasa,probe_radius);
	return rsd_sasa;
}

double
burial_weight(pose::Pose p, int position){
	//takes in a position in a pose and returns a weight based on some sort
	//of metric for burial
	utility::vector1<core::Real> sa;
	sa = sasa(p);
	if(true){
		if(sa[position] < 0.1){
			return 1.0;
		}else if(sa[position] > 0.4){
			return 0.5;
		}else{
			return 0.75;
		}
	}
}
**/
int
main(int argc, char* argv []){
    try {
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::ddg;

	OPT(ddg::mean);
	OPT(ddg::min);
	OPT(in::file::s);

  devel::init( argc, argv );

	utility::vector1<utility::file::FileName> in = basic::options::option[in::file::s]();

	bool mean, min;

	if(basic::options::option[OptionKeys::ddg::mean].user()){
		mean = basic::options::option[OptionKeys::ddg::mean]();
	}else{
		mean = false;
	}
	if(basic::options::option[OptionKeys::ddg::min].user()){
		min = basic::options::option[OptionKeys::ddg::min]();
	}else{
		min = false;
	}


	ScoreFunctionOP sfxn(get_score_function());

	protocols::ddG::ddGData dat(in[1].name());
	bool header_printed = false;

	std::ofstream outfile((basic::options::option[out::file::scorefile]()).c_str(),
												std::ios_base::app);

	while(!dat.end()){
		dat.get_next_filenames();
		//read in wt pdbs and mutant pdbs
		utility::vector1<pose::Pose> mut=dat.read_mut_data();
		utility::vector1<pose::Pose> wt=dat.read_wt_data();
		core::Real experimental_value = dat.read_exp_data();
		//		utility::vector1<double> env_wts((header(mut[1],sfxn)).size(),1.0);
		utility::vector1<double> ddgs = calc_ddg(mut, wt, sfxn,mean,min);//hard-coded only for debugging
		if(!header_printed){
			header_printed=true;
			print_header(mut[1],sfxn,outfile);
		}
		if(!basic::options::option[OptionKeys::ddg::print_per_res_diff]()){
			print_ddgs(ddgs,mut,wt,
								 experimental_value,
								 outfile);
		}else{
			print_verbose_ddgs(mut,wt,sfxn,mean,min,experimental_value,outfile);
		}
	}
    } catch (utility::excn::Exception const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}

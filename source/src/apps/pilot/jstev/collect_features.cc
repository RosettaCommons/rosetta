// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Collects features to determine tryptophan viability for indole complementation candidacy
/// @author Jason Stevens

//C++ Headers
#include <iostream>
#include <fstream>
#include <sstream>

//Boost Headers
#include <boost/lexical_cast.hpp>

//Protocol Headers
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NonlocalContactsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

//Core Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>

//Utility Headers
#include <utility/vector1.hh>
#include <utility/io/util.hh>
#include <utility/io/irstream.hh>
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;

class Features
{
private:
	Size index_;
	Size pdb_seq_pos_;
	Size bb_hbonds_;
	Size sc_hbonds_;
	Real sc_sasa_;
	Real packstat_;
	Size nlcontacts_;
	Size num_neighbors_;
	Size num_hyneighbors_;
	Real sc_temp_ave_;

public:
	void set_index(Size val) { index_ = val; }
	void set_pdb_seq_pos(Size val) { pdb_seq_pos_ = val; };
	void set_bb_hbonds(Size val) { bb_hbonds_ = val; }
	void set_sc_hbonds(Size val) { sc_hbonds_ = val; }
	void set_sc_sasa(Real val) { sc_sasa_ = val; }
	void set_packstat(Real val) { packstat_ = val; }
	void set_nlcontacts(Size val) { nlcontacts_ = val; }
	void set_num_neighbors(Size val) { num_neighbors_ = val; }
	void set_num_hyneighbors(Size val) { num_hyneighbors_ = val; }
	void set_sc_temp_ave(Real val) { sc_temp_ave_ = val; }
	void normalize_sc_temp_ave(Real val) { sc_temp_ave_ = sc_temp_ave_/val; }

	void print_features(std::ostream &output);

	Features();
};

Features::Features()
{
	index_ = 0;
	pdb_seq_pos_ = 0;
	bb_hbonds_ = 0;
	sc_hbonds_ = 0;
	sc_sasa_ = 0;
	packstat_ = 0;
	nlcontacts_ = 0;
	num_neighbors_ = 0;
	num_hyneighbors_ = 0;
	sc_temp_ave_ = 0;
}

void Features::print_features(std::ostream &output)
{
	output << index_ << "  " << pdb_seq_pos_ << "  " << bb_hbonds_ << "  " << sc_hbonds_ << "  " << sc_sasa_ << "  " << packstat_ << "  ";
	output << nlcontacts_ << "  " << num_neighbors_ << "  " << num_hyneighbors_ << "  " << sc_temp_ave_ << std::endl << std::endl;
}


int
main(int argc, char* argv[]){
    try {

	devel::init(argc, argv);

	pose::Pose input_pose;
	std::string input_pdb_name (basic::options::start_file());
	std::string input_pdb_atom = input_pdb_name; input_pdb_atom.append("_atom.pdb");
	core::import_pose::pose_from_file(input_pose, input_pdb_atom, core::import_pose::PDB_file);

	scoring::ScoreFunctionOP scorefxn(get_score_function());
	(*scorefxn)(input_pose);

	// mjo commented out to fix unused variable warning
	// Real const starting_total_score = input_pose.energies().total_energies()[total_score];
	//std::cout << "Total score at start is: " << starting_total_score << std::endl;

	//each trp residue has a vector of features
	utility::vector1<Features> trp_features;

	//indices vectors
	//utility::vector1<Size> trp_site_indices;
	utility::vector1<Size> trp_residue_indices;

	//total temperature
	Real total_temp = 0;

	//look at the SITE lines of the pdb file, store them in a string vector
	std::string input_pdb_site = input_pdb_name; input_pdb_site.append("_site.pdb");
	utility::vector1<std::string> pdb_site_strings;
	std::string line;
	std::ifstream in; in.open(input_pdb_site.c_str());

	while (in >> line) {

		pdb_site_strings.push_back(line);
	}

	//find the indices of the trp residues in the active site
	/*for (utility::vector1<std::string>::iterator iter = pdb_site_strings.begin(); iter != pdb_site_strings.end(); ++iter) {

		if (*iter == "TRP") {
			++iter; ++iter;
			trp_site_indices.push_back(utility::string2int(*iter));
		}
	}*/

	//loop through all residues
	for (Size i = 1; i <= input_pose.size(); ++i){

		//find tryptophans, and put their indices (NOT the actual sequence position) into a vector
		if (name_from_aa(input_pose.aa(i)) == "TRP") {

			trp_residue_indices.push_back(i);
			trp_features.push_back(Features());
		}
	}

	//set up a calculator for calculating the number of hydrogen bonds
	core::pose::metrics::PoseMetricCalculatorOP hb_calc = new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
	core::pose::metrics::CalculatorFactory::Instance().register_calculator("hbonds", hb_calc);

	//find the number of hydrogen bonds for each trp residue
	/*basic::MetricValue<utility::vector1<Size> > hb_metric;
	input_pose.metric("hbonds", "residue_Hbonds", hb_metric);
	utility::vector1<Size> total_hbonds = hb_metric.value();*/

	//find the number of hydrogen bonds for each trp residue
	basic::MetricValue<core::id::AtomID_Map<Size> > hb_metric;
	input_pose.metric("hbonds", "atom_Hbonds", hb_metric);
	core::id::AtomID_Map<Size> atom_hbonds = hb_metric.value();

	//set up a calculator for calculating the SASA for each trp residue
	core::pose::metrics::PoseMetricCalculatorOP sasa_calc = new core::pose::metrics::simple_calculators::SasaCalculatorLegacy();
	core::pose::metrics::CalculatorFactory::Instance().register_calculator("sasa", sasa_calc);

	//find the sasa value for each residue
	/*basic::MetricValue<utility::vector1<Real> > sasa_metric;
	input_pose.metric("sasa", "residue_sasa", sasa_metric);
	utility::vector1<Real> total_sasas = sasa_metric.value();*/

	//find the sasa value for each
	basic::MetricValue<core::id::AtomID_Map<Real> > sasa_metric;
	input_pose.metric("sasa", "atom_sasa", sasa_metric);
	core::id::AtomID_Map<Real> atom_sasas = sasa_metric.value();

	//set up a calculator for calculating the packstat for each trp residue
	core::pose::metrics::PoseMetricCalculatorOP packstat_calc = new protocols::toolbox::pose_metric_calculators::PackstatCalculator();
	core::pose::metrics::CalculatorFactory::Instance().register_calculator("packstat", packstat_calc);

	//find the packstat value for each residue
	basic::MetricValue<utility::vector1<Real> > packstat_metric;
	input_pose.metric("packstat", "residue_packstat", packstat_metric);
	utility::vector1<Real> total_packstats = packstat_metric.value();

	//set up a calculator for calculating the nlcontacts for each trp residue
	core::pose::metrics::PoseMetricCalculatorOP nlcontacts_calc = new protocols::toolbox::pose_metric_calculators::NonlocalContactsCalculator();
	core::pose::metrics::CalculatorFactory::Instance().register_calculator("nlcontacts", nlcontacts_calc);

	//find the packstat value for each residue
	basic::MetricValue<utility::vector1<Size> > nlcontacts_metric;
	input_pose.metric("nlcontacts", "residue_nlcontacts", nlcontacts_metric);
	utility::vector1<Size> total_nlcontacts = nlcontacts_metric.value();

	Size curr_index = 1;

	//loop through all of the tryptophan residues
	for (utility::vector1<Size>::iterator iter = trp_residue_indices.begin(); iter != trp_residue_indices.end(); ++iter){

		Size curr_trp_pos = input_pose.residue(*iter).seqpos();

		//create residue
		conformation::Residue current_trp_residue(input_pose.residue(*iter));

//---NON-FEATURE FEATURES---

		//add the residue sequence position to the feature vector
		trp_features[curr_index].set_pdb_seq_pos(input_pose.pdb_info()->number(*iter));

		//add the residue index to the feature vector
		trp_features[curr_index].set_index(curr_trp_pos);


//---FEATURE 1: Number of hydrogen bonds for each trp residue (backbone and sidechain)---

		//loop through the AtomID_Map of atom hydrogen bonds and separate backbone from sidechain
		//********NOTE: make sure that seqpos and ii are the right ways to index************
		Size bb_hbonds = 0, sc_hbonds = 0;

		for(Size ii = 1; ii <= current_trp_residue.nheavyatoms(); ++ii) {

			if (current_trp_residue.atom_is_backbone(ii))
				bb_hbonds += atom_hbonds(curr_trp_pos, ii);
			else //can assume is sidechain?
				sc_hbonds += atom_hbonds(curr_trp_pos, ii);
		}

		trp_features[curr_index].set_bb_hbonds(bb_hbonds);
		trp_features[curr_index].set_sc_hbonds(sc_hbonds);


		//---(DEPRECATED) FEATURE 2: Binary indication of active site inclusion---

		/*bool active = false;

		//if the current trp residue is in the active site set its boolean state to TRUE
		for (utility::vector1<std::string>::iterator iter1 = pdb_site_strings.begin(); iter1 != pdb_site_strings.end(); ++iter1) {

			if (utility::string2int(*iter1) == input_pose.pdb_info()->number(*iter)) active = true;
		}

		trp_features[i].push_back(active);*/


		//---(DEPRECATED) FEATURE 3: SASA value for each trp residue---

		//get the sasa value for the current trp and add it to the feature vector
		//trp_features[i].set_total_sasa(total_sasas[*iter]);


//---FEATURE 3.5: Sidechain SASA values---

		//loop through the AtomID_Map of atom sasa bonds and separate backbone from sidechain
		//********NOTE: make sure that seqpos and ii are the right ways to index************
		Real sc_sasa = 0;

		for(Size ii = 1; ii <= current_trp_residue.nheavyatoms(); ++ii) {

			if (!current_trp_residue.atom_is_backbone(ii)) //can assume is sidechain?
				sc_sasa += atom_sasas(curr_trp_pos, ii);
		}

		trp_features[curr_index].set_sc_sasa(sc_sasa);


//---FEATURE 4: Packstat value for each trp residue---

		//get the packstat value for the current trp and add it to the feature vector
		trp_features[curr_index].set_packstat(total_packstats[*iter]);


//---FEATURE 5: Nonlocal contacts for each trp residue---

		//get the packstat value for the current trp and add it to the feature vector
		trp_features[curr_index].set_nlcontacts(total_nlcontacts[*iter]);


//---FEATURE 6: Neighbors by distance for each trp residue---

		//set up a calculator for calculating the neighbors by distance for each trp residue
		//******NOTE: I'm not sure if the central_residue value is meant to be the index or the real seqpos (I'm doing the latter)*******
		core::pose::metrics::PoseMetricCalculatorOP nbydist_calc =
			new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator(input_pose.pdb_info()->number(*iter), 5);
		std::string calc_name = "nbydist";
		calc_name.append(boost::lexical_cast<std::string>(*iter));
		core::pose::metrics::CalculatorFactory::Instance().register_calculator(calc_name, nbydist_calc);

		//find the number of neighbors for the current residue
		basic::MetricValue<Size> nbydist_metric;
		input_pose.metric(calc_name, "num_neighbors", nbydist_metric);
		Size total_neighbors = nbydist_metric.value();

		//get the neighbors by distance value for each trp and add it to the feature vector
		trp_features[curr_index].set_num_neighbors(total_neighbors);


//---FEATURE 7: Number of Hydrophobic Neighbors---

		//get the aa of the neighbors (may need to convert to 1-20)
		basic::MetricValue<std::set<Size> > neighbors_metric;
		input_pose.metric(calc_name, "neighbors", neighbors_metric);
		std::set<Size> neighbors = neighbors_metric.value();

		Size num_hyd = 0;

		for(std::set<Size>::iterator ssiter = neighbors.begin(); ssiter != neighbors.end(); ++ssiter){

			//convert pdb sequence position to pose number.
			//********NOTE: that i'm assuming the neighbors are on the same chain as the current trp residue************
			Size num = input_pose.pdb_info()->pdb2pose(input_pose.pdb_info()->chain(*iter), *ssiter, ' ');

			//********NOTE: i'm assuming non-charged, not polar residues are non-polar*************
			if ((input_pose.residue(num).is_charged() == false) && (input_pose.residue(num).is_polar() == false))
				++num_hyd;
		}

		trp_features[curr_index].set_num_hyneighbors(num_hyd);

		std::cout << std::endl << std::endl;


//---FEATURE 8: Average, normalized b-factor/temp for each trp sidechain---

		Real temp = 0, ave_res_temp = 0;

		//loop through all the atoms and average the temperatures; increase the total
		for (Size iii = 1; iii != current_trp_residue.nheavyatoms(); ++iii) {

			if (!current_trp_residue.atom_is_backbone(iii)) { //can assume is sidechain?

				temp = input_pose.pdb_info()->temperature(curr_trp_pos, iii);
				total_temp += temp;
				ave_res_temp += temp;
			}
		}

		trp_features[curr_index].set_sc_temp_ave(ave_res_temp/current_trp_residue.nheavyatoms());

		++curr_index;
	}

	std::cout << "ind  sp  bbh  sch  scsa  pkst  nlc  neig  hydr  temp\n\n";


	//normalize the sidechain temps for each residue and print the features for each residue
	for (utility::vector1<Features>::iterator iter = trp_features.begin(); iter != trp_features.end(); ++iter){

		iter->normalize_sc_temp_ave(total_temp);
		iter->print_features(std::cout);
	}

    } catch (utility::excn::Exception const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
    }
    return 0;
}

/*

current features list:

-number of hydrogen bonds -- sidechain and backbone
-sidechain sasa
-packstat for each residue
-non-local contacts for each residue
-number of neighbors within 5 A
-number of hydrophobic neighbors within 5 A
-normalized (divided by entire structure's temp) temperature (b-factor) average for sidechain


other, non-rosetta features list (all done in pymol):

-binary indication of whether it touches the substrate - pymol
-overall stability of the protein (non-rosetta: literature or temp the organism grows at)
-shortest distance to substrate (closest two) - pymol


training (can use more data for training)
1) if i make a mutation will i lose function - functional => no rescue potential
2) predict aggregation => no rescue potential
3) if passes 1 and 2, test for rescue


NOTES:

make a program to chop up the pdb files into atom and site files
(or maybe you can find some other way through Rosetta)


shitcode:

*/

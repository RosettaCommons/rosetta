// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file calc_pair_stats.cc
/// @brief
/// @author Robert Vernon

#include <core/types.hh>
#include <devel/init.hh>

#include <core/chemical/AA.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/pose/Pose.hh>
// Auto-header: duplicate removed #include <devel/init.hh>

#include <core/kinematics/RT.hh>
#include <core/io/pdb/pdb_writer.hh>


#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/robert.OptionKeys.gen.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>

#include <numeric/random/random.hh>


#include <numeric/xyzMatrix.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <protocols/evaluation/RmsdEvaluator.hh>

// Unit headers
#include <protocols/scoring/methods/pcs/PseudocontactShiftEnergy.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftEnergyCreator.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/jumping/StrandPairing.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>


///////////////////////////////////////////////////////////////////////////////

using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using ObjexxFCL::FArray2A_float;
using ObjexxFCL::FArray2D_float;


using namespace core;
using namespace pose;
using namespace conformation;


void invert_exclude_residues( Size nres, utility::vector1<int> const& exclude_list, protocols::simple_filters::ResidueSelection& residue_selection ) {

	residue_selection.clear();

	for( Size ir = 1; ir <= nres; ++ir ) {
		bool exclude_residue = false;
		for( Size ex = 1; ex <= exclude_list.size(); ex ++ ){
			if( int(exclude_list[ex]) == int(ir) ) {
				exclude_residue = true;
				break;
			}
		}

		if ( !exclude_residue ) {
			residue_selection.push_back( ir );
		}
	} // for ( Size ir = 1; ir <= native_pose.size(); ++ir )
}


int
main( int argc, char* argv [] )
{

	try {

	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;

	using namespace core::scoring::methods::pcs;

	// options, random initialization
	devel::init( argc, argv );

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	// configure score function
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::scoring::EnergyMap emap;
	emap.zero();

	core::io::silent::SilentFileData sfd;

	std::string infile	= *(option[ in::file::silent ]().begin());

	if ( option[ in::file::silent ].user() ) {
		sfd.read_file( infile );
	}

	utility::io::ozstream outfile(infile+".30.chunks");

	core::Real const grid_edge = 16;
	core::Real const grid_step = 1.5;
	core::Real const grid_small_cutoff = 4;
	core::Real const grid_large_cutoff = 8;
	core::Real const grid_cone_angle_cutoff = 91;
	std::string const grid_atom_name_1 = "CB";
	std::string const grid_atom_name_2 = "CA";
	core::SSize const grid_residue_num_1 = 12;
	core::SSize const grid_residue_num_2 = 12;
	core::Real const grid_k_vector = 0;
	bool const minimize_best_tensor = true;
	core::Real const pcs_weight = 10;

	PCS_Energy_parameters_manager::get_instance()->set_grid_param(grid_edge,
																																grid_step,
																																grid_small_cutoff,
																																grid_large_cutoff,
																																grid_cone_angle_cutoff,
																																grid_atom_name_1,
																																grid_atom_name_2,
																																grid_residue_num_1,
																																grid_residue_num_2,
																																grid_k_vector,
																																minimize_best_tensor,
																																pcs_weight
																																);

	utility::vector1<std::string> vec_filename;
	vec_filename.push_back("/work/rvernon/christophe/targets_v4/epsilon/PCS_eps_Dy_CNH.npc.4rescore");
	vec_filename.push_back("/work/rvernon/christophe/targets_v4/epsilon/PCS_eps_Er_CNH.npc.4rescore");
	vec_filename.push_back("/work/rvernon/christophe/targets_v4/epsilon/PCS_eps_Tb_CNH.npc.4rescore");

	utility::vector1<core::Real> vec_individual_weight;
	vec_individual_weight.push_back(1.0);
	vec_individual_weight.push_back(1.0);
	vec_individual_weight.push_back(1.0);

	PCS_Energy_parameters_manager::get_instance()->set_vector_name_and_weight(vec_filename,
																																						vec_individual_weight);


// 	if (option[ in::file::native_exclude_res ].user()) {

// 		utility::vector1< Size > vec_exclude;
// 		vec_exclude = option[ in::file::native_exclude_res ]();
// 		PCS_Energy_parameters_manager::get_instance()->set_vector_exclude_residues(vec_exclude);


// 		utility::vector1< bool > vec_exclude_bool;

// 		vec_exclude_bool = PCS_Energy_parameters_manager::get_instance()->get_vector_exclude_residues();


// 	}

	core::pose::Pose native_pose, pose;

	if ( option[ in::file::native ].user() ) {
		// read in pdb and constraints if necessary
		core::import_pose::pose_from_file(
			native_pose, *rsd_set, option[ in::file::native ]()
		);
	}

	Size const window_size(30);

	Size count = 0;
	for ( SilentFileData::iterator iter = sfd.begin(), end = sfd.end();
				iter != end; ++iter
	) {
		count++;

		iter->fill_pose( pose, *rsd_set );

		for (Size i = 1; i <= pose.size() - window_size; ++i) {

			Size r_start(i), r_end(i+window_size);
			utility::vector1< Size > vec_exclude;

			for (Size o = 1; o < r_start; ++o) {
				vec_exclude.push_back(o);
			}

			for (Size o = (r_end + 1); o <= pose.size(); ++o) {
				vec_exclude.push_back(o);
			}

			PCS_Energy_parameters_manager::get_instance()->set_vector_exclude_residues(vec_exclude);

			PCS_Energy pcs_energy;

			core::Real pcs = pcs_energy.calculate_pcs_score(pose, false);

			outfile << "WINDOW T_" << count << " " << i << " " << pcs << " ";

			if ( option[ in::file::native ].user() ) {
				protocols::simple_filters::ResidueSelection residues;
				invert_exclude_residues( native_pose.size(), vec_exclude, residues );
				core::Real rmsd = core::scoring::CA_rmsd( pose, native_pose, residues );
				core::Real gdtmm = core::scoring::CA_gdtmm( pose, native_pose, residues );

				outfile << rmsd << " " << gdtmm;
			}

			outfile << std::endl;

		}

		std::cout << "DONE #" << count << std::endl;

	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

} // int main( int argc, char * argv [] )

//}
//}

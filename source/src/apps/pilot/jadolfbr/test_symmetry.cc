// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/jadolfbr/test_neighborhood_selector.cc
/// @brief Test neighborhood selector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>


// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

#include <utility/string_util.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/SymmetricalResidueSelector.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>

static basic::Tracer TR("test_symmetry");


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;



	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );

}

core::Size count_subset( utility::vector1< bool > subset ){
	core::Size counts = 0;
	for ( bool const & b : subset ) {
		if ( b ) {
			counts+=1;
		}
	}
	return counts;
}

std::string
get_attachment_point_string( core::pose::Pose const & pose, core::Size resnum){
	using utility::to_string;
	using namespace core::chemical::carbohydrates;

	CarbohydrateInfoCOP info = pose.residue(resnum).carbohydrate_info();
	std::string outstring = "";
	std::string attach = "_->";

	if ( info->mainchain_glycosidic_bond_acceptor() ) {
		outstring = attach + to_string(info->mainchain_glycosidic_bond_acceptor());
	}

	for ( uint i = 1; i <= info->n_branches(); ++i ) {
		outstring = outstring + "," +attach + to_string( info->branch_point( i ));
	}
	return outstring;
}

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::select::residue_selector;
		using namespace protocols::simple_moves::symmetry;
		devel::init( argc, argv );
		register_options();


		core::pose::PoseOP pose = core::import_pose::pose_from_file( "/Users/jadolfbr/2ciw.pdb");

		std::cout << "Total Residues: " << pose->total_residue() << std::endl;


		NeighborhoodResidueSelector nbr_selector = NeighborhoodResidueSelector();
		GlycanResidueSelectorOP glycan_selector = GlycanResidueSelectorOP( new GlycanResidueSelector() );
		SymmetricalResidueSelector symm_selector = SymmetricalResidueSelector();
		core::scoring::ScoreFunctionOP score  = core::scoring::get_score_function();

		score->score(*pose);

		glycan_selector->set_select_from_branch_residue( pose->pdb_info()->pdb2pose('A', 93) );
		glycan_selector->set_include_root( true );
		utility::vector1< bool > glycan_residues = glycan_selector->apply( *pose );

		std::cout << "Total Glycans: " << pose->glycan_tree_set()->size() << std::endl;
		std::cout << "Chains: " << pose->num_chains() << std::endl;

		std::cout << "Glycan Selection: " << count_subset( glycan_residues ) << std::endl;
		nbr_selector.set_focus( glycan_residues );
		nbr_selector.set_distance( 12 );
		nbr_selector.set_include_focus_in_subset( false );

		utility::vector1< bool > nbr_residues = nbr_selector.apply( *pose );
		std::cout << "NBR Selection: " << count_subset( nbr_residues ) << std::endl;

		core::Size protein_branches = 0;
		core::Size carbohydrate_residues = 0;
		for ( core::Size resnum = 1; resnum <= pose->size(); ++resnum ) {
			if ( pose->residue( resnum ).is_carbohydrate() ) {
				std::string attachment_points = get_attachment_point_string( *pose, resnum);
				core::Size parent_res = pose->glycan_tree_set()->get_parent( resnum );
				bool bp = pose->residue( resnum ).is_branch_point();


				std::cout << "Carbohydrate: "<< resnum  <<" "<< pose->pdb_info()->pose2pdb(resnum) << " Parent: " << parent_res << " BP: "<<bp <<" "<< pose->pdb_info()->pose2pdb(resnum) << " " << " CON: " << " DIS: " << pose->glycan_tree_set()->get_distance_to_start( resnum )
					<< utility::pad_right( attachment_points, 10) << pose->residue( resnum ).carbohydrate_info()->short_name() << std::endl;

				carbohydrate_residues += 1;

			} else if ( pose->residue( resnum ).is_branch_point() ) {
				std::cout << "Branch Point: " << pose->residue( resnum ).name3()<<" "<< resnum <<" " <<pose->pdb_info()->pose2pdb(resnum) << std::endl;
				protein_branches += 1;

			}
		}
		std::cout << "Glycan Residues: " << carbohydrate_residues <<std::endl;
		std::cout << "Protein BPs: " << protein_branches << std::endl;



		std::cout << "Symmetrizing Pose. " << std::endl;

		///Now we print info on the symmetrical version of the pose.
		SetupForSymmetryMover symm = SetupForSymmetryMover();
		symm.apply(*pose);
		score->score(*pose);

		std::cout << "Total Residues: " << pose->total_residue() << std::endl;

		std::cout << "Total Glycans: " << pose->glycan_tree_set()->size() << std::endl;
		std::cout << "Chains: " << pose->num_chains() << std::endl;


		symm_selector.set_selector( glycan_selector );
		glycan_residues = symm_selector.apply( *pose );

		std::cout << "Glycan Selection: " << count_subset( glycan_residues ) << std::endl;
		nbr_selector.set_focus( glycan_residues );
		nbr_selector.set_distance( 6 );
		nbr_selector.set_include_focus_in_subset( false );

		nbr_residues = nbr_selector.apply( *pose );
		std::cout << "NBR Selection: " << count_subset( nbr_residues ) << std::endl;

		for ( core::Size resnum = 1; resnum <= pose->size(); ++resnum ) {
			if ( pose->residue( resnum ).is_carbohydrate() ) {
				std::string attachment_points = get_attachment_point_string( *pose, resnum);
				core::Size parent_res = pose->glycan_tree_set()->get_parent( resnum );
				bool bp = pose->residue( resnum ).is_branch_point();


				std::cout << "Carbohydrate: "<< resnum  <<" "<< pose->pdb_info()->pose2pdb(resnum) << " Parent: " << parent_res << " BP: "<<bp <<" "<< pose->pdb_info()->pose2pdb(resnum) << " " << " CON: " << " DIS: " << pose->glycan_tree_set()->get_distance_to_start( resnum )
					<< utility::pad_right( attachment_points, 10) << pose->residue( resnum ).carbohydrate_info()->short_name() << std::endl;

				carbohydrate_residues += 1;

			} else if ( pose->residue( resnum ).is_branch_point() ) {
				std::cout << "Branch Point: " << pose->residue( resnum ).name3()<<" "<< resnum <<" " <<pose->pdb_info()->pose2pdb(resnum) << std::endl;
				protein_branches += 1;

			}
		}
		std::cout << "Glycan Residues: " << carbohydrate_residues <<std::endl;
		std::cout << "Protein BPs: " << protein_branches << std::endl;
		pose->dump_pdb("/Users/jadolfbr/2ciw_symm.pdb");


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

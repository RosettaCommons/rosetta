// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/match_enzdes_util/ReactionState.hh
/// @brief  class to model one state of a reaction
/// @author Florian Richter, flosopher@gmail.com, feb/mar 2012

// Unit headers
#include <protocols/toolbox/match_enzdes_util/InvrotTree.hh>

// package headers
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTarget.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>


// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/rotamer_set/bb_independent_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <basic/Tracer.hh>

//utility headers
#include <utility>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// c++ headers
#include <fstream>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static basic::Tracer tr( "protocols.toolbox.match_enzdes_util.InvrotTree" );


InvrotTree::InvrotTree()
: ReferenceCount()
{
	invrot_targets_.clear();
	invrot_tree_constraints_.clear();
}

InvrotTree::~InvrotTree() = default;

core::scoring::constraints::ConstraintCOP
InvrotTree::get_constraint_for_target_state( Size target_state ) const {
	return invrot_tree_constraints_[target_state];
}

void
InvrotTree::generate_inverse_rotamer_constraints(
	core::pose::Pose const & pose,
	AllowedSeqposForGeomCstCOP geomcst_seqpos
)
{
	if ( invrot_targets_.size() == 0 ) utility_exit_with_message("InvrotTree is asked to generate constraints even though no target states exist. Something is wrong somewhere.");

	invrot_tree_constraints_.clear();

	for ( Size i = 1; i <= invrot_targets_.size(); ++i ) {
		invrot_tree_constraints_.push_back( invrot_targets_[i]->generate_constraints( pose, geomcst_seqpos ) );
	}

}

utility::vector1< InvrotCollectorCOP >
InvrotTree::collect_all_inverse_rotamers( ) const
{
	utility::vector1< InvrotCollectorOP > invrot_collectors;
	for ( Size i =1; i <= invrot_targets_.size(); ++i ) invrot_targets_[i]->collect_all_inverse_rotamers( invrot_collectors );
	return invrot_collectors;
}

void
InvrotTree::dump_invrots_tree_as_multimodel_pdbs( std::string filename_base) const
{
	tr << "Writing InvrotTree to file(s).... " << std::endl;
	//1. collect all the invrots
	utility::vector1< InvrotCollectorCOP > invrot_collectors( this->collect_all_inverse_rotamers() );
	if ( invrot_collectors.size() == 0 ) {
		tr <<"Error when trying to dump inverse rotamer tree to file. No inverse rotamers found in tree.";
		return;
	}

	//2. write them
	Size files_to_write( invrot_collectors.size() );
	tr <<"A total of " << files_to_write << " unique definitions of the invrot tree exist." << std::endl;

	for ( Size i =1; i <= files_to_write; ++i ) {
		std::string filename( filename_base + "_" + utility::to_string( i ) + ".pdb" );
		std::vector< std::list< core::conformation::ResidueCOP > > const & invrot_lists( invrot_collectors[i]->invrots() );
		Size num_rotamer_lists( invrot_lists.size() );
		std::vector< std::list< core::conformation::ResidueCOP >::const_iterator > res_iterators( num_rotamer_lists );

		tr << "Writing definition " << i << " to file " << filename << "... " << std::endl;
		for ( Size j =0; j < num_rotamer_lists; ++j ) {
			res_iterators[ j ] = invrot_lists[j].begin();
			tr << invrot_lists[j].size() << " invrots for list " << j << ", ";
		}
		tr << std::endl;

		core::Size atomcounter(1), modelcount(1);
		bool all_res_iterators_at_end( false );
		std::ofstream file_out( filename.c_str() );
		file_out << "MODEL   1\n";

		while ( !all_res_iterators_at_end ) {
			all_res_iterators_at_end = true;
			for ( core::Size j =0; j < num_rotamer_lists; ++j ) {
				if ( res_iterators[ j ] != invrot_lists[j].end() ) {
					all_res_iterators_at_end = false;
					core::io::pdb::dump_pdb_residue( **(res_iterators[ j ]), atomcounter, file_out );
					res_iterators[ j ]++;
				}
			}
			if ( !all_res_iterators_at_end ) {
				file_out << "ENDMDL \n";
				file_out << "MODEL    "+utility::to_string( ++modelcount )+"\n";
			}
		} // while( !all_res_iterators_at_end )
		file_out.close();
	} //loop over files to write
}

void
InvrotTree::clear_target_states()
{
	invrot_targets_.clear();
}

void
InvrotTree::add_to_targets( InvrotTargetOP invrot_target )
{
	invrot_targets_.push_back( invrot_target );
}


TheozymeInvrotTree::TheozymeInvrotTree( EnzConstraintIOCOP enzcst_io )
: InvrotTree(), enzcst_io_(std::move(enzcst_io))
{}

TheozymeInvrotTree::~TheozymeInvrotTree()= default;


bool
TheozymeInvrotTree::check_pose_tree_compatibility(
	core::pose::Pose & ) const
{
	utility_exit_with_message("stubbed out");
	return false;  // required for compilation on Windows
}

void
TheozymeInvrotTree::generate_targets_and_inverse_rotamers()
{
	this->clear_target_states();

	//1. make a residue from the first block ligand
	utility::vector1< core::chemical::ResidueTypeCOP > ds_restypes( enzcst_io_->mcfi_list( 1 )->mcfi( 1 )->allowed_restypes( enzcst_io_->mcfi_list( 1 )->mcfi( 1 )->downstream_res() ) );
	core::conformation::ResidueCOP ligres( core::conformation::ResidueFactory::create_residue( *(ds_restypes[1]) ) );

	utility::vector1< core::conformation::ResidueCOP > all_rots;
	//amino acid target?
	if ( ligres->is_protein() ) all_rots = core::pack::rotamer_set::bb_independent_rotamers( ligres->type_ptr() );

	//ligand?
	else if ( ligres->is_ligand() ) {
		//for now
		//all_rots.push_back( ligres );
		all_rots = core::pack::rotamer_set::bb_independent_rotamers( ligres->type_ptr() );
	} else {
		//one residue only
		all_rots.push_back( ligres );
	}
	//2. determine redundancy, same way as it's determined in matcher task
	// ..... no code yet .. will probably have to refactor redundancy determination a bit
	// easiest will be to put the redundancy determining code in LigandConformerBuilder
	//into a utility function in LigandConformer, and move the relevant atom determining
	//code from the matcher task to EnzConstraintIO
	//
	//hack for now
	utility::vector1<  utility::vector1< core::conformation::ResidueCOP > > conformer_groups;
	conformer_groups.push_back( all_rots );

	//but at afer this, we'll generate a bunch of ligand targets
	for ( Size i =1; i <= conformer_groups.size(); ++i ) {

		SingleResidueInvrotTargetOP invtarg( new SingleResidueInvrotTarget( conformer_groups[i] ) );
		if ( !invtarg->initialize_tree_nodes_from_enzcst_io( enzcst_io_ ) ) tr << "Target from conformer group " << i << " failed to initialize, cstfile geometry (clashes?) is bad somewhere." << std::endl;
		else this->add_to_targets( invtarg );
	}
	if ( this->num_target_states() == 0 ) utility_exit_with_message("No targets could be initialized by TheozymeInvrotTree. check cstfile");
}

} //match_enzdes_util
} //toolbox
} //protocols

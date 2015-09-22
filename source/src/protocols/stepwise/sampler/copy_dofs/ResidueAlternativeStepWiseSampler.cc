// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSampler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSampler.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeSet.hh>
#include <protocols/simple_moves/CopyDofMover.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.copy_dofs.ResidueAlternativeStepWiseSampler" );

using namespace core;
using namespace core::conformation;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace copy_dofs {

//constructor
ResidueAlternativeStepWiseSampler::ResidueAlternativeStepWiseSampler( ResidueAlternativeSet const & residue_alternative_set,
	core::pose::Pose const & starting_pose ):
	CopyDofStepWiseSampler( residue_alternative_set.pose_list(),
	residue_alternative_set.res_map(),
	starting_pose ),
	representative_seqpos_( residue_alternative_set.representative_seqpos() )
{
	initialize_residues();
}

//constructor
ResidueAlternativeStepWiseSampler::ResidueAlternativeStepWiseSampler( ResidueAlternativeSet const & residue_alternative_set ):
	CopyDofStepWiseSampler( residue_alternative_set.pose_list(),
	residue_alternative_set.res_map() ),
	representative_seqpos_( residue_alternative_set.representative_seqpos() )
{
	initialize_residues();
}

//Constructor
ResidueAlternativeStepWiseSampler::ResidueAlternativeStepWiseSampler( utility::vector1< core::pose::PoseOP > const & pose_list,
	std::map< Size, Size > const & res_map,
	Size const representative_seqpos,
	core::pose::Pose const & starting_pose ):
	CopyDofStepWiseSampler( pose_list, res_map, starting_pose ),
	representative_seqpos_( representative_seqpos )
{
	initialize_residues();
}

//Constructor
ResidueAlternativeStepWiseSampler::ResidueAlternativeStepWiseSampler( utility::vector1< pose::PoseOP > const & pose_list,
	std::map< Size, Size > const & res_map,
	Size const representative_seqpos ):
	CopyDofStepWiseSampler( pose_list, res_map ),
	representative_seqpos_( representative_seqpos )
{
	initialize_residues();
}

//Constructor
ResidueAlternativeStepWiseSampler::ResidueAlternativeStepWiseSampler( utility::vector1< pose::PoseOP > const & pose_list,
	Size const representative_seqpos ):
	CopyDofStepWiseSampler( pose_list, simple_res_map( representative_seqpos ) ),
	representative_seqpos_( representative_seqpos )
{
	initialize_residues();
}

//Destructor
ResidueAlternativeStepWiseSampler::~ResidueAlternativeStepWiseSampler()
{}

///////////////////////////////////////////////////////
std::map< Size, Size >
ResidueAlternativeStepWiseSampler::simple_res_map( Size const i ){
	std::map< Size, Size > res_map;
	res_map[ i ] = i;
	return res_map;
}

///////////////////////////////////////////////////////
void
ResidueAlternativeStepWiseSampler::initialize_residues(){
	original_type_ = pose_list_[ 1 ]->residue( representative_seqpos_ ).name(); // could be a hash if speed is needed
	utility::vector1< ResidueOP > residues;
	for ( Size n = 1; n <= pose_list_.size(); n++ ) {
		pose::PoseOP pose = pose_list_[ n ];
		runtime_assert( pose->residue( representative_seqpos_ ).name() == original_type_ );
		residues.push_back( pose->residue( representative_seqpos_ ).clone() );
	}
	residues_for_each_type_[ original_type_ ] = residues;
}

///////////////////////////////////////////////////////
void
ResidueAlternativeStepWiseSampler::initialize_residues_for_type( Residue const & rsd_in ){
	utility::vector1< ResidueOP > residues;
	for ( Size n = 1; n <= pose_list_.size(); n++ ) {
		pose::PoseOP pose = pose_list_[ n ];
		Residue const & original_rsd = pose->residue( representative_seqpos_ );
		residues.push_back( ResidueOP( new Residue( rsd_in.type(), original_rsd, pose->conformation() /*pray this works*/ ) ) );
	}
	std::string type_in = rsd_in.name(); // could be a hash.
	residues_for_each_type_[ type_in ] = residues;
}

///////////////////////////////////////////////////////
Residue const &
ResidueAlternativeStepWiseSampler::get_residue_at_origin(){
	return *residues_for_each_type_[ original_type_ ][ id_ ];
}

///////////////////////////////////////////////////////
Residue const &
ResidueAlternativeStepWiseSampler::get_residue_at_origin_with_matching_type( Residue const & rsd_in ){
	std::string type_in = rsd_in.name(); // could be a hash.
	if (  residues_for_each_type_.find( type_in ) == residues_for_each_type_.end() ) {
		initialize_residues_for_type( rsd_in );
	}
	return *residues_for_each_type_[ type_in ][ id_ ];
}

//////////////////////////////////////////////////////////////////////////
/// @brief Name of the class
std::string
ResidueAlternativeStepWiseSampler::get_name() const {
	return "ResidueAlternative residue:" + utility::to_string( representative_seqpos_ );
}

} //copy_dofs
} //sampler
} //stepwise
} //protocols

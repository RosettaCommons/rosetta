// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/constraints/InverseRotamersRCG.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2010

#include <protocols/forge/constraints/InverseRotamersRCG.hh>
#include <protocols/forge/build/Interval.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/id/SequenceMapping.hh>
#include <basic/Tracer.hh>

static basic::Tracer tr( "protocols.forge.constraints.InverseRotamersRCG" );

namespace protocols{
namespace forge{
namespace constraints{

InverseRotamersRCG::InverseRotamersRCG(
	core::Size lstart,
	core::Size lstop,
	std::list< core::conformation::ResidueCOP > const & inverse_rotamers
	) : constraint_func_(NULL), func_sd_(0.4)
{
	intervals_.clear();
	inverse_rotamers_.clear();
	intervals_.push_back( forge::build::Interval( lstart, lstop ) );
	for( std::list< core::conformation::ResidueCOP >::const_iterator rot_it( inverse_rotamers.begin() ), rot_end( inverse_rotamers.end() );
			 rot_it != rot_end; ++rot_it ){
		inverse_rotamers_.push_back( *rot_it );
	}
}

InverseRotamersRCG::~InverseRotamersRCG(){}

void
InverseRotamersRCG::generate_remodel_constraints(
	core::pose::Pose const & pose )
{
	using namespace core::scoring::constraints;
	//safeguard against bad user input
	if( inverse_rotamers_.size() == 0 ){
		std::cerr << "WARNING: InverseRotamersRCG is asked to produce constraints but was not given any inverse rotamers. Something's probably wrong somewhere." << std::endl;
		return;
	}

	//if no constraint func has been set, we'll create a default one
	if( !constraint_func_ ){
		constraint_func_ = new BoundFunc( 0, 0.05, func_sd_, "invrot");
	}

	//see the comment in protocols/ligand_docking/LigandBaseProtocol.cc::restrain_protein_Calphas
	core::id::AtomID fixed_pt( pose.atom_tree().root()->atom_id() );

	utility::vector1< ConstraintCOP > all_res_invrot_csts;
	core::Size totrescount(0);

	for( core::Size i(1); i <= intervals_.size(); ++i ){

		//eventually remap intervals according to vlb seqmap
		if( this->seqmap() ){
			intervals_[i].left = (*(this->seqmap() ))[ intervals_[i].left ];
			intervals_[i].right = (*(this->seqmap() ))[ intervals_[i].right ];
		}

		for( core::Size remres( intervals_[i].left ); remres <= intervals_[i].right; ++remres ){
			core::conformation::ResidueCOP cur_remodel_res( &pose.residue( remres ) );
			if( cur_remodel_res->name3() == "GLY" ) continue;
			totrescount++;

			core::id::AtomID rem_CA( cur_remodel_res->type().atom_index("CA"), remres );
			core::id::AtomID rem_CB( cur_remodel_res->type().atom_index("CB"), remres );
			core::id::AtomID rem_N( cur_remodel_res->type().atom_index("N"), remres );
			//core::id::AtomID rem_C( cur_remodel_res->type().atom_index("C"), remres );

			for( core::Size invrot(1); invrot <= inverse_rotamers_.size(); ++invrot ){

				utility::vector1< ConstraintCOP > cur_res_invrot_csts;
				cur_res_invrot_csts.push_back( new BackboneStubConstraint( pose, remres, fixed_pt, *(inverse_rotamers_[invrot]), -20.0, 0.8) );
				//old style: coordinate constraints for all atoms, backbone stub csts
				// might be working better

				//utility::vector1< ConstraintCOP > cur_res_invrot_csts;
				cur_res_invrot_csts.push_back( new CoordinateConstraint( rem_CA, fixed_pt, inverse_rotamers_[invrot]->xyz("CA"), constraint_func_ ) );
				cur_res_invrot_csts.push_back( new CoordinateConstraint( rem_CB, fixed_pt, inverse_rotamers_[invrot]->xyz("CB"), constraint_func_ ) );
				cur_res_invrot_csts.push_back( new CoordinateConstraint( rem_N, fixed_pt, inverse_rotamers_[invrot]->xyz("N"), constraint_func_ ) );
				//cur_res_invrot_csts.push_back( new CoordinateConstraint( rem_C, fixed_pt, inverse_rotamers_[invrot]->xyz("C"), constraint_func_ ) );

				all_res_invrot_csts.push_back( new MultiConstraint( cur_res_invrot_csts ) );


			} //loop over inverse rotamers
		} // loop over all remodel residues in this interval
	} //loop over all intervals
	tr << "Created a total of " << all_res_invrot_csts.size() << " constraints between " << inverse_rotamers_.size() << " inverse rotamers and " << totrescount << " residues in " << intervals_.size() << " intervals." << std::endl;

	this->add_constraint( new AmbiguousConstraint( all_res_invrot_csts ) );

	//we can probably delete the inverse rotamers now, to save some memory
	this->clear_inverse_rotamers();
}

void
InverseRotamersRCG::set_constraint_func(
	core::scoring::constraints::FuncOP constraint_func ){
	constraint_func_ = constraint_func;
}

void
InverseRotamersRCG::clear_inverse_rotamers()
{
	inverse_rotamers_.clear();
}

} //namespace remodel
} //namespace forge
} //namespace protocols

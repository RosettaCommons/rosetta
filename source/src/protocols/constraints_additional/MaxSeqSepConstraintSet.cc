// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MaxSeqSepConstraintSet.cc
/// @brief
/// @details
/// @author Oliver Lange


// Unit Headers
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>

// Package Headers
//#include <protocols/jumping/JumpSetup.fwd.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
/*
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/EnergyMap.fwd.hh>
*/
/*
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh> // SequenceMover
#include <protocols/moves/TrialMover.hh>
*/

// ObjexxFCL Headers

// Utility headers
#include <utility>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/numeric.functions.hh>

#include <basic/Tracer.hh>


//// C++ headers
#include <cstdlib>
#include <string>

#include <core/id/SequenceMapping.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.constraints_additional.MaxSeqSepConstraintSet", basic::t_info );

using core::scoring::constraints::ConstraintSet;
using core::scoring::constraints::ConstraintSetOP;
using core::kinematics::ShortestPathInFoldTree;
using core::kinematics::ShortestPathInFoldTreeOP;
using core::kinematics::ShortestPathInFoldTreeCOP;

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace constraints_additional {

using namespace core;

/// @brief a ConstraintsSet whose constraints can be switched off, according to sequence separation in residues
/// between residue pair constraints.
MaxSeqSepConstraintSet::~MaxSeqSepConstraintSet() = default;
MaxSeqSepConstraintSet::MaxSeqSepConstraintSet( ConstraintSet const & other, core::kinematics::FoldTree const & f ) :
	ConstraintSet( other )
{
	tr.Trace << f << std::endl;
	shortest_path_ = core::kinematics::ShortestPathInFoldTreeOP( new ShortestPathInFoldTree( f ) );
}

/// @copy constructor. Performs a shallow copy of the constraints and the ShortestPathInFoldTree object
MaxSeqSepConstraintSet::MaxSeqSepConstraintSet( MaxSeqSepConstraintSet const & ) = default;

MaxSeqSepConstraintSet::MaxSeqSepConstraintSet( ConstraintSet const &other, ShortestPathInFoldTreeCOP sp ) :
	ConstraintSet( other ),
	shortest_path_(std::move( sp ))
{}

ConstraintSet &
MaxSeqSepConstraintSet::operator = ( ConstraintSet const & rhs )
{
	if ( this != &rhs ) {
		MaxSeqSepConstraintSet const * msscs_rhs = dynamic_cast< MaxSeqSepConstraintSet const * > ( & rhs );
		if ( ! msscs_rhs ) {
			throw utility::excn::EXCN_Msg_Exception( "MaxSeqSepConstraintSet handed a non MaxSeqSepConstraintSet in operator =" );
		}
		ConstraintSet::operator = ( rhs );

		if ( max_seq_sep_ != msscs_rhs->max_seq_sep_ ) {
			max_seq_sep_ = msscs_rhs->max_seq_sep_;
			mark_revision_id_expired();
		}
		if ( shortest_path_ != msscs_rhs->shortest_path_ ) {
			shortest_path_ = msscs_rhs->shortest_path_;
			mark_revision_id_expired();
		}

		// mark_revision_id_expired();
	}
	return *this;
}

ConstraintSetOP
MaxSeqSepConstraintSet::clone() const {
	return ConstraintSetOP( new MaxSeqSepConstraintSet( *this ) );
}

void MaxSeqSepConstraintSet::detached_copy( ConstraintSet const & src ) {
	MaxSeqSepConstraintSet const * msscs_src = dynamic_cast< MaxSeqSepConstraintSet const * > ( & src );
	if ( ! msscs_src ) {
		throw utility::excn::EXCN_Msg_Exception( "MaxSeqSepConstraintSet handed a non MaxSeqSepConstraintSet in detatched_copy" );
	}
	deep_copy( src );
	max_seq_sep_ = msscs_src->max_seq_sep_;
	shortest_path_ = core::kinematics::ShortestPathInFoldTreeOP( new core::kinematics::ShortestPathInFoldTree( *msscs_src->shortest_path_ ));
}


ConstraintSetOP
MaxSeqSepConstraintSet::detached_clone() const
{
	MaxSeqSepConstraintSetOP newset( new MaxSeqSepConstraintSet );
	newset->detached_copy( *this );
	return newset;
}

bool
MaxSeqSepConstraintSet::same_type_as_me( ConstraintSet const & other, bool recurse /* = true */ ) const
{
	MaxSeqSepConstraintSet const * msscs_other = dynamic_cast< MaxSeqSepConstraintSet const * > ( & other );
	if ( ! msscs_other ) return false;
	if ( recurse ) {
		return other.same_type_as_me( *this, false );
	}
	return true;
}


ConstraintSetOP MaxSeqSepConstraintSet::remapped_clone(
	core::pose::Pose const& src,
	core::pose::Pose const& dest,
	core::id::SequenceMappingCOP smap
) const {
	MaxSeqSepConstraintSetOP clone_ptr( new MaxSeqSepConstraintSet( *Parent::remapped_clone( src, dest, smap ), shortest_path_ ) );
	clone_ptr->set_max_seq_sep( max_seq_sep() );
	return clone_ptr;
}


void
MaxSeqSepConstraintSet::residue_pair_energy(
	Residue const & rsd1,
	Residue const & rsd2,
	Pose const & pose,
	core::scoring::ScoreFunction const & scorefxn,
	core::scoring::EnergyMap & emap
) const
{
	int const pos1( rsd1.seqpos() ), pos2 ( rsd2.seqpos() );
	// if ( tr.Trace.visible() ) tr.Trace << "max_seq_sep(): " << max_seq_sep() << " check seqsep for residues " << rsd1.seqpos() << "  and  " << rsd2.seqpos() << "....";
	if ( too_far( pos1, pos2 ) ) {
		//if ( tr.Trace.visible() ) tr.Trace << "\n";
		return; //cast avoids warning
	}
	// if ( tr.Trace.visible() ) tr.Trace << "evaluated\n";
	ConstraintSet::residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap );
}

void
MaxSeqSepConstraintSet::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & minmap,
	core::scoring::ResSingleMinimizationData const & res1_data_cache,
	core::scoring::ResSingleMinimizationData const & res2_data_cache,
	core::scoring::ResPairMinimizationData & respair_data_cache
) const
{
	if ( too_far( rsd1.seqpos(), rsd2.seqpos() ) ) return;
	Parent::setup_for_minimizing_for_residue_pair( rsd1, rsd2, pose, sfxn, minmap, res1_data_cache, res2_data_cache, respair_data_cache );
}


Size
MaxSeqSepConstraintSet::show_violations( std::ostream& out, pose::Pose& pose, Size verbose_level, Real threshold ) {
	out << " total constr: " << get_all_constraints().size() << "  ";
	out << " max_seq_sep: " << max_seq_sep() << " ";
	return Parent::show_violations( out, pose, verbose_level, threshold );
}

bool
MaxSeqSepConstraintSet::too_far( int const pos1, int const pos2 ) const {
	if ( shortest_path_ ) return max_seq_sep() < shortest_path_->dist( pos1, pos2 );
	return max_seq_sep() < (core::Size) std::abs( pos1 - pos2 );
}

Size
MaxSeqSepConstraintSet::largest_possible_sequence_sep( core::pose::Pose const& pose ) const {
	if ( shortest_path_ ) return shortest_path_->max_dist();
	return pose.size();
}

/// Does *NOT* zero the emap values, just adds the additional contribution to the
/// existing emap energies (so can be called inside finalize_total_energies)
void
MaxSeqSepConstraintSet::eval_non_residue_pair_energy(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::scoring::EnergyMap & emap
) const
{
	//non_residue_pair_constraints_.conformation_energy( pose.conformation(), sfxn.weights(), emap );
	core::scoring::func::ConformationXYZ const xyz_func( pose.conformation() );
	core::scoring::constraints::ConstraintCOPs const& csts( non_residue_pair_constraints().constraints() );
	for ( auto const & it : csts ) {
		core::scoring::constraints::Constraint const & cst( *it );
		if ( cst.effective_sequence_separation( shortest_path() ) < max_seq_sep() ) cst.score( xyz_func, sfxn.weights(), emap );
	}
}

MaxSeqSepConstraintSet::MaxSeqSepConstraintSet() {}

} //abinitio
} //protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::constraints_additional::MaxSeqSepConstraintSet::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::constraints::ConstraintSet >( this ) );
	arc( CEREAL_NVP( max_seq_sep_ ) ); // Size
	arc( CEREAL_NVP( shortest_path_ ) ); // core::kinematics::ShortestPathInFoldTreeCOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::constraints_additional::MaxSeqSepConstraintSet::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::constraints::ConstraintSet >( this ) );
	arc( max_seq_sep_ ); // Size
	std::shared_ptr< core::kinematics::ShortestPathInFoldTree > local_shortest_path;
	arc( local_shortest_path ); // core::kinematics::ShortestPathInFoldTreeCOP
	shortest_path_ = local_shortest_path; // copy the non-const pointer(s) into the const pointer(s)
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::constraints_additional::MaxSeqSepConstraintSet );
CEREAL_REGISTER_TYPE( protocols::constraints_additional::MaxSeqSepConstraintSet )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_constraints_additional_MaxSeqSepConstraintSet )
#endif // SERIALIZATION

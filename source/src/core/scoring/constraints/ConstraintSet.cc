// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// Unit headers
#include <core/scoring/constraints/ConstraintSet.hh>

// Package headers
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/DOF_Constraint.hh> //Special constraints.
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/constraints/CstMinimizationData.hh>

// Project headers
#include <core/conformation/Conformation.hh> //for attaching/detaching
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>  //need to pass in weights
#include <basic/Tracer.hh>
#include <basic/prof.hh>
#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>

// C++ Headers
#include <set>

#include <core/id/SequenceMapping.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/LengthEvent.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace constraints {

/// @details Auto-generated virtual destructor
ResidueConstraints::~ResidueConstraints() {}

#ifdef    SERIALIZATION
template < class Archive >
void
ResidueConstraints::save( Archive & arc ) const
{
	arc( map_ );
}

template < class Archive >
void
ResidueConstraints::load(
	Archive & arc
)
{
	arc( map_ );
}

SAVE_AND_LOAD_SERIALIZABLE( ResidueConstraints );

#endif // SERIALIZATION

static THREAD_LOCAL basic::Tracer tr( "core.scoring.ConstraintSet" );

ConstraintSet::ConstraintSet()
:
	total_residue_( 0 ),
	sequence_constraints_(),
	revision_id_( 0 ),
	revision_id_current_( false ),
	conformation_pt_( /* NULL */ )
{}

ConstraintSet::ConstraintSet( ConstraintSet const & other )
:
	ReferenceCount(),
	total_residue_( other.total_residue_ ),
	sequence_constraints_( other.sequence_constraints_.size() ),
	residue_pair_constraints_( other.residue_pair_constraints_.size() ),
	revision_id_( 0 ),
	revision_id_current_( false ),
	conformation_pt_( /* NULL */ )
{
	shallow_copy( other, 1, other.total_residue_ );
}

ConstraintSet::ConstraintSet(
	ConstraintSet const & other,
	Size start_residue,
	Size end_residue
) :
	ReferenceCount(),
	total_residue_( end_residue ),
	sequence_constraints_( other.sequence_constraints_.size() ),
	residue_pair_constraints_( end_residue ),
	revision_id_( 0 ),
	revision_id_current_( false ),
	conformation_pt_( /* NULL */ )
{
	shallow_copy( other, start_residue, end_residue );
}

/// @brief Destructor must detach from conformation
ConstraintSet::~ConstraintSet() { this->detach_from_conformation(); }

ConstraintSetOP ConstraintSet::clone() const {
	return ConstraintSetOP( new ConstraintSet( *this ) );
}

/// @details This can be called by derived classes to make sure
/// that all of the base class data is efficiently copied.
ConstraintSet &
ConstraintSet::operator = ( ConstraintSet const & rhs ) {

	if ( this == &rhs ) return *this;

	if ( residue_pair_constraints_.size() != rhs.residue_pair_constraints_.size() ) {
		residue_pair_constraints_.resize( rhs.residue_pair_constraints_.size() );
		mark_revision_id_expired();
	}

	for ( Size ii = 1; ii <= residue_pair_constraints_.size(); ++ii ) {
		if ( ! rhs.residue_pair_constraints_[ ii ] && ! residue_pair_constraints_[ ii ] ) continue;

		if ( ! rhs.residue_pair_constraints_[ ii ] ) {
			residue_pair_constraints_[ ii ].reset();
			mark_revision_id_expired();
			continue;
		} else if ( ! residue_pair_constraints_[ ii ] ) {
			residue_pair_constraints_[ ii ] = ResidueConstraintsOP( new ResidueConstraints );
			mark_revision_id_expired();
		}

		ResiduePairConstraintsIterator rhs_iter = rhs.residue_pair_constraints_[ ii ]->begin();
		ResiduePairConstraintsIterator this_iter = residue_pair_constraints_[ ii ]->begin();
		ResiduePairConstraintsIterator rhs_iter_end = rhs.residue_pair_constraints_[ ii ]->end();
		ResiduePairConstraintsIterator this_iter_end = residue_pair_constraints_[ ii ]->end();
		while ( rhs_iter != rhs_iter_end && this_iter != this_iter_end ) {
			if ( rhs_iter->first == this_iter->first ) {
				// assignment operator for the Constraints class will perform pointer comparisons
				// to avoid copying data it doesn't have to
				(*this_iter->second) = (*rhs_iter->second );
				++this_iter;
				++rhs_iter;
			} else if ( this_iter->first < rhs_iter->first ) {
				mark_revision_id_expired();
				ResiduePairConstraintsIterator next_this_iter( this_iter );
				++next_this_iter;
				residue_pair_constraints_[ ii ]->erase( this_iter->first );
				this_iter = next_this_iter;
			} else {
				mark_revision_id_expired();
				residue_pair_constraints_[ ii ]->insert( rhs_iter->first, rhs_iter->second->clone() );
				++rhs_iter;
			}
		}

		// handle extra constraints in this; i.e. there are some constraints that must
		// be deleted because they aren't in the rhs ConstrainSet
		while ( this_iter != this_iter_end ) {
			mark_revision_id_expired();
			ResiduePairConstraintsIterator next_this_iter( this_iter );
			++next_this_iter;
			residue_pair_constraints_[ ii ]->erase( this_iter->first );
			this_iter = next_this_iter;
		}
		// handle extra constraints in rhs; i.e. there are some constraints in the rhs
		// ConstraintSet that should be copied over into this
		while ( rhs_iter != rhs_iter_end ) {
			mark_revision_id_expired();
			residue_pair_constraints_[ ii ]->insert( rhs_iter->first, rhs_iter->second->clone() );
			++rhs_iter;
		}
	}

	ResiduePairConstraintsIterator rhs_intra_iter = rhs.intra_residue_constraints_.begin();
	ResiduePairConstraintsIterator this_intra_iter = intra_residue_constraints_.begin();
	ResiduePairConstraintsIterator rhs_intra_iter_end = rhs.intra_residue_constraints_.end();
	ResiduePairConstraintsIterator this_intra_iter_end = intra_residue_constraints_.end();

	while ( rhs_intra_iter != rhs_intra_iter_end && this_intra_iter != this_intra_iter_end ) {
		if ( rhs_intra_iter->first == this_intra_iter->first ) {
			// assignment operator for the Constraints class will perform pointer comparisons
			// to avoid copying data it doesn't have to
			(*this_intra_iter->second) = (*rhs_intra_iter->second );
			++this_intra_iter;
			++rhs_intra_iter;
		} else if ( this_intra_iter->first < rhs_intra_iter->first ) {
			mark_revision_id_expired();
			ResiduePairConstraintsIterator next_this_intra_iter( this_intra_iter );
			++next_this_intra_iter;
			intra_residue_constraints_.erase( this_intra_iter->first );
			this_intra_iter = next_this_intra_iter;
		} else {
			mark_revision_id_expired();
			intra_residue_constraints_.insert( rhs_intra_iter->first, rhs_intra_iter->second->clone() );
			++rhs_intra_iter;
		}
	}

	// handle extra constraints in this
	while ( this_intra_iter != this_intra_iter_end ) {
		ResiduePairConstraintsIterator next_this_intra_iter( this_intra_iter );
		++next_this_intra_iter;
		intra_residue_constraints_.erase( this_intra_iter->first );
		this_intra_iter = next_this_intra_iter;
	}
	// handle extra constraints in rhs
	while ( rhs_intra_iter != rhs_intra_iter_end ) {
		intra_residue_constraints_.insert( rhs_intra_iter->first, rhs_intra_iter->second->clone() );
		++rhs_intra_iter;
	}

	non_residue_pair_constraints_ = rhs.non_residue_pair_constraints_;
	dof_constraints_ = rhs.dof_constraints_;

	total_residue_ = rhs.total_residue_;

	// shallow copy of the sequence constraints.
	sequence_constraints_.resize( rhs.sequence_constraints_.size() );
	for ( Size ii = 1; ii <= rhs.sequence_constraints_.size(); ++ii ) {
		if ( sequence_constraints_[ ii ] != rhs.sequence_constraints_[ ii ] ) {
			// avoid smart-pointer copy if unnecessary.
			sequence_constraints_[ ii ] = rhs.sequence_constraints_[ ii ];
		}
	}

	return *this;
}

void ConstraintSet::detached_copy( ConstraintSet const & src ) {
	deep_copy( src );
}

ConstraintSetOP
ConstraintSet::detached_clone() const
{
	ConstraintSetOP newset( new ConstraintSet );
	newset->detached_copy( *this );
	return newset;
}

bool
ConstraintSet::same_type_as_me( ConstraintSet const & other, bool recurse /* = true */ ) const {
	if ( recurse ) {
		return other.same_type_as_me( *this, false );
	}
	return true;
}


/// @brief Copies the data from this ConstraintSet into a new object and returns
/// an OP atoms are mapped to atoms with the same name in dest pose ( e.g. for
/// switch from centroid to fullatom ) if a sequence_mapping is present it is
/// used to map residue numbers .. NULL = identity mapping to the new object.
/// This will really clone all constraints since they have to change their
/// atom-numbers and residue-numbers
ConstraintSetOP ConstraintSet::remapped_clone(
	pose::Pose const & src,
	pose::Pose const & dest,
	id::SequenceMappingCOP smap /*default NULL*/
) const {
	ConstraintCOPs all_cst = get_all_constraints();
	ConstraintSetOP new_set( new ConstraintSet );
	for ( ConstraintCOPs::const_iterator it = all_cst.begin(), eit = all_cst.end(); it!=eit; ++it ) {
		ConstraintCOP new_cst = (*it)->remapped_clone( src, dest, smap );
		if ( new_cst ) new_set->add_constraint( new_cst );
	}
	return new_set;
}

/// @brief
/// same as remapped clone, but it also steals coordinates from src pose
ConstraintSetOP ConstraintSet::steal_def_clone(
	pose::Pose const & src,
	pose::Pose const & dest,
	id::SequenceMappingCOP smap /*default NULL*/
) const {
	ConstraintCOPs all_cst = get_all_constraints();
	ConstraintSetOP new_set( new ConstraintSet );
	for ( ConstraintCOPs::const_iterator it = all_cst.begin(), eit = all_cst.end(); it!=eit; ++it ) {
		ConstraintOP new_cst = (*it)->clone();
		new_cst->steal_def( src );
		new_cst = new_cst->remapped_clone( src, dest, smap );
		if ( new_cst ) new_set->add_constraint( new_cst );
	}
	return new_set;
}


void
ConstraintSet::remap_residue_positions(
	id::SequenceMapping const & smap
) {

	if ( ! this->has_constraints() ) return; //in this case we don't have to worry about anything

	ConstraintCOPs all_cst = get_all_constraints();

	debug_assert( all_cst.size() != 0 );

	//nuke the current constraints
	clear();

	for ( ConstraintCOPs::const_iterator it = all_cst.begin(), eit = all_cst.end(); it!=eit; ++it ) {
		ConstraintCOP new_cst = (*it)->remap_resid( smap );
		if ( new_cst ) this->add_constraint( new_cst );
		else tr.Debug << "when remapping the constraint set, one constraint could not be remapped. :( "<< std::endl;
	}
}

void
ConstraintSet::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData & res_data_cache
) const
{
	ResidueConstraints::const_iterator iter = intra_residue_constraints_.find( rsd.seqpos() );
	if ( iter != intra_residue_constraints_.end() ) {
		res_data_cache.set_data( cst_res_data, basic::datacache::CacheableDataOP( new CstMinimizationData( iter->second ) ) );
	}
}


void
ConstraintSet::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData & respair_data_cache
) const
{
	using basic::datacache::CacheableDataOP;
	Size resno1 = rsd1.seqpos();
	Size resno2 = rsd2.seqpos();
	if ( ! residue_pair_constraints_[ resno1 ] ) return;
	ResidueConstraints::const_iterator iter = residue_pair_constraints_[ resno1 ]->find( resno2 );
	if ( iter != residue_pair_constraints_[ resno1 ]->end() ) {
		respair_data_cache.set_data( cst_respair_data, CacheableDataOP( new CstMinimizationData( iter->second ) )  );
	}
}

void
ConstraintSet::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const & scfxn
) const {
	// iterate through all constraints and call setup_for_scoring
	func::ConformationXYZ confxyz( pose.conformation() );
	for ( int i=1; i <= (int) residue_pair_constraints_.size() ; ++i ) {
		if ( !residue_pair_constraints_[i] ) continue;

		ResidueConstraints const & seqpos_constraints( *residue_pair_constraints_[ i ] );
		for ( ResidueConstraints::const_iterator
				it= seqpos_constraints.begin(),
				ite= seqpos_constraints.end();
				it != ite; ++it ) {
			it->second->setup_for_scoring( confxyz, scfxn );
		}
	}

	for ( ResidueConstraints::const_iterator it=intra_residue_constraints_.begin(), it_end=intra_residue_constraints_.end(); it != it_end;
			++it ) {
		it->second->setup_for_scoring( confxyz, scfxn );
	}

	non_residue_pair_constraints_.setup_for_scoring( confxyz, scfxn );
}

void
ConstraintSet::setup_for_derivatives( pose::Pose &pose, ScoreFunction const &scfxn ) const {
	// iterate through all constraints and call setup_for_scoring
	func::ConformationXYZ confxyz( pose.conformation() );
	for ( int i=1; i <= (int) residue_pair_constraints_.size() ; ++i ) {
		if ( !residue_pair_constraints_[i] ) continue;

		ResidueConstraints const & seqpos_constraints( *residue_pair_constraints_[ i ] );
		for ( ResidueConstraints::const_iterator it= seqpos_constraints.begin(),
				ite= seqpos_constraints.end(); it != ite; ++it ) {
			it->second->setup_for_derivatives( confxyz, scfxn );
		}
	}

	for ( ResidueConstraints::const_iterator it=intra_residue_constraints_.begin(), it_end=intra_residue_constraints_.end(); it != it_end;
			++it ) {
		it->second->setup_for_derivatives( confxyz, scfxn );
	}

	non_residue_pair_constraints_.setup_for_derivatives( confxyz, scfxn );
}


void
ConstraintSet::residue_pair_energy(
	Residue const & rsd1,
	Residue const & rsd2,
	Pose const &, // pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const
{
	Size const pos1( rsd1.seqpos() ), pos2( rsd2.seqpos() );
	if ( residue_pair_constraints_.size() < pos1 || !residue_pair_constraints_[pos1] ) return;
	ResidueConstraints const & pos1_constraints( *residue_pair_constraints_[ pos1 ] );
	ResidueConstraints::const_iterator it( pos1_constraints.find( pos2 ) );
	if ( it == pos1_constraints.end() ) return;
	Constraints const & pair_constraints( *(it->second) );
	pair_constraints.residue_pair_energy( rsd1, rsd2, scorefxn.weights(), emap );
}

/// This could be made much more efficient by a mapping from AtomID's to constraints
void
ConstraintSet::deprecated_eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	ScoreFunction const & scfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	// Derivatives for intraresidue and two-body constraints are evaluated using the MinimizationGraph

	Size const seqpos( atom_id.rsd() );
	deprecated_eval_atom_derivative_for_residue_pairs( atom_id, pose, scfxn, weights, F1, F2 );

	{ // intraresidue constraints
		ResidueConstraints::const_iterator it( intra_residue_constraints_.find( seqpos ) );
		if ( it != intra_residue_constraints_.end() ) {
			it->second->eval_ws_atom_derivative( atom_id, pose.conformation(), weights, F1, F2 );
		}
	}

	{ // nonpair constraints
		non_residue_pair_constraints_.eval_ws_atom_derivative( atom_id, pose.conformation(), weights, F1, F2 );
	}

}

void
ConstraintSet::eval_multibody_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	non_residue_pair_constraints_.eval_ws_atom_derivative( atom_id, pose.conformation(), weights, F1, F2 );
}

void
ConstraintSet::deprecated_eval_atom_derivative_for_residue_pairs(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	// residue pair constraints:
	Size const seqpos( atom_id.rsd() );
	if ( residue_pair_constraints_.size() >= seqpos && residue_pair_constraints_[ seqpos ] ) {
		ResidueConstraints const & seqpos_constraints( *residue_pair_constraints_[ seqpos ] );
		for ( ResidueConstraints::const_iterator it= seqpos_constraints.begin(), ite= seqpos_constraints.end(); it != ite;
				++it ) {
			it->second->eval_ws_atom_derivative( atom_id, pose.conformation(), weights, F1, F2 );
		}
	}
}


void
ConstraintSet::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const &, // pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	ResidueConstraints::const_iterator it( intra_residue_constraints_.find( rsd.seqpos() ) );
	if ( it != intra_residue_constraints_.end() ) {
		it->second->intra_residue_energy( rsd, sfxn.weights(), emap );
	}
}


void
ConstraintSet::eval_intrares_energy(
	conformation::Residue const & rsd,
	EnergyMap & emap
) const
{
	ResidueConstraints::const_iterator it( intra_residue_constraints_.find( rsd.seqpos() ) );
	if ( it != intra_residue_constraints_.end() ) {
		it->second->intra_residue_energy( rsd, emap /*dummy -- not actually used in this function? */, emap );
	}
}


/// Does *NOT* zero the emap values, just adds the additional contribution to the
/// existing emap energies (so can be called inside finalize_total_energies)
void
ConstraintSet::eval_non_residue_pair_energy(
	Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	non_residue_pair_constraints_.conformation_energy( pose.conformation(), sfxn.weights(), emap );

	// Real dof_score(0.0);
	// for ( DOF_ConstraintOPs::const_iterator it=dof_constraints_.begin(), ite = dof_constraints_.end(); it != ite; ++it ) {
	//  DOF_ConstraintOP const & dof_constraint( *it );
	//  dof_score = dof_constraint->func( pose.dof( dof_constraint->dof_id() ) );
	//  emap[ dof_constraint->score_type() ] += dof_score;
	// }

}

/*Real
ConstraintSet::eval_dof_derivative(
id::DOF_ID const & id,
id::TorsionID const &, // tor,
pose::Pose const & pose,
ScoreFunction const &, // scorefxn,
EnergyMap const & weights
) const
{
if ( dof_constraints_.empty() ) return 0.0;

//  DOF_ConstraintOPs::const_iterator it( dof_constraints_.find( id ) ); // will this be too slow??
//  Real deriv(0.0);
//  if ( it != dof_constraints_.end() ) {
//   DOF_ConstraintOP const & dof_constraint( it->second );
//   deriv = weights[ dof_constraint->score_type() ] * dof_constraint->dfunc( pose.dof( it->first ) );
//  }
Real deriv( 0.0 );
for ( DOF_ConstraintOPs::const_iterator it=dof_constraints_.begin(), ite = dof_constraints_.end(); it != ite; ++it ) {
DOF_ConstraintOP const & dof_constraint( *it );
if ( dof_constraint->dof_id() != id ) continue;
deriv += weights[ dof_constraint->score_type() ] * dof_constraint->dfunc( pose.dof( dof_constraint->dof_id() ) );
}

return deriv;
}*/


/// helper, static
void
add_constraint_to_residue_constraints(
	int const seqpos,
	ConstraintCOP cst,
	ResidueConstraints & residue_constraints
)
{
	if ( !residue_constraints.has( seqpos ) ) {
		residue_constraints.insert( seqpos, ConstraintsOP( new Constraints() ) );
		//residue_constraints.insert( std::make_pair( seqpos, ConstraintsOP( new Constraints() ) ) );
	}
	residue_constraints.find( seqpos )->second->add_constraint( cst );
}

void
ConstraintSet::add_residue_pair_constraint( Size const pos1, Size const pos2, ConstraintCOP cst )
{
	if (  total_residue_ < pos1 ) total_residue_ = pos1;
	if (  total_residue_ < pos2 ) total_residue_ = pos2;
	if (  residue_pair_constraints_.size() < pos1 ) residue_pair_constraints_.resize( pos1, 0 );
	if ( !residue_pair_constraints_[ pos1 ] ) residue_pair_constraints_[ pos1 ] = ResidueConstraintsOP( new ResidueConstraints() );

	add_constraint_to_residue_constraints( pos2, cst, *(residue_pair_constraints_[ pos1 ] ) );
}


void
ConstraintSet::add_constraints( ConstraintCOPs cst_list ) {
	for ( ConstraintCOPs::iterator it = cst_list.begin(), end = cst_list.end();
			it != end; ++it ) {
		add_constraint( *it );
	}
}

/// @details copy another constraint set into this one
void
ConstraintSet::add_constraints( ConstraintSetCOP const cst_set ) {
	add_constraints( cst_set->get_all_constraints() );
	return;
}


void
ConstraintSet::add_constraint( ConstraintCOP cst )
{
	mark_revision_id_expired();

	aa_composition_energy::SequenceConstraintCOP seq_cst( utility::pointer::dynamic_pointer_cast< aa_composition_energy::SequenceConstraint const > (cst) ); //see whether this is a SequenceConstraint
	if ( seq_cst ) { // If it is a sequence constraint, store it as such
		add_sequence_constraint( seq_cst );
	} else {

		// figure out if it's inter-res, residue_pair, or 3+body
		utility::vector1< int > pos_list( cst->residues() );

		if ( pos_list.size() == 1 ) {
			// intra-res
			//  tr.Trace << "add intra-res constraint " << std::endl;
			if ( total_residue_ < Size(pos_list[1]) ) total_residue_ = pos_list[1];
			add_constraint_to_residue_constraints( pos_list[1], cst, intra_residue_constraints_ );
		} else if ( pos_list.size() == 2 ) {
			// rsd-pai
			//  tr.Trace << "add res constraint " << std::endl;
			add_residue_pair_constraint( pos_list[1], pos_list[2], cst );
			add_residue_pair_constraint( pos_list[2], pos_list[1], cst );
		} else {
			// 3+ body
			//  tr.Trace << "add 3+body constraint " << std::endl;
			for ( core::Size ii = 1; ii <= pos_list.size(); ++ii ) {
				if ( total_residue_ < Size(pos_list[ ii ]) ) total_residue_ = pos_list[ii];
			}
			non_residue_pair_constraints_.add_constraint( cst );
		}
	}
	return;
}

/// helper, static
bool
remove_constraint_from_residue_constraints(
	int const seqpos,
	ConstraintCOP cst,
	ResidueConstraints & residue_constraints,
	bool object_comparison
)
{
	if ( !residue_constraints.has( seqpos ) ) return false;

	ResidueConstraints::iterator csts_it = residue_constraints.find( seqpos );

	//if( residue_constraints.find( seqpos )->second->remove_constraint( cst ) ) {
	if ( csts_it->second->remove_constraint( cst, object_comparison ) ) {

		if ( csts_it->second->size() == 0 ) residue_constraints.erase( seqpos );

		return true;
	}

	return false;
}


/// private
// this function is called twice with interchanged pos1/pos2 parameters ---> all constraints are removed symmetrically
bool
ConstraintSet::remove_residue_pair_constraint(
	Size const pos1,
	Size const pos2,
	ConstraintCOP cst,
	bool object_comparison  )
{
	if (  residue_pair_constraints_.size() < pos1 ) return false;
	if ( !residue_pair_constraints_[ pos1 ] ) return false;

	bool return_val = remove_constraint_from_residue_constraints( pos2, cst, *(residue_pair_constraints_[ pos1 ] ), object_comparison );

	//don't forget to resize residue_pair_constraints_ if we removed the last constraint from the back
	if ( ( (*residue_pair_constraints_[ pos1 ]).size() == 0 ) && residue_pair_constraints_.size() == pos1 ) {

		for ( core::Size ii = residue_pair_constraints_.size(); ii >= 1; --ii ) {

			if ( residue_pair_constraints_[ ii ] ) {
				if ( (*residue_pair_constraints_[ ii ]).size() != 0 ) break;
			}
			residue_pair_constraints_.pop_back();
		}
	}

	return return_val;
}


void
ConstraintSet::add_sequence_constraint(
	core::scoring::aa_composition_energy::SequenceConstraintCOP cst
)
{
	sequence_constraints_.push_back( cst );
}

bool
ConstraintSet::remove_constraints(
	ConstraintCOPs cst_list,
	bool object_comparison  )
{
	bool success = false;
	for ( ConstraintCOPs::iterator it = cst_list.begin(), end = cst_list.end();
			it != end; ++it ) {
		success = remove_constraint( *it, object_comparison );
		if ( success == false ) return false;
	}

	return success;
} //remove_constraints function


bool
ConstraintSet::remove_constraint(
	ConstraintCOP cst,
	bool object_comparison
)
{
	mark_revision_id_expired();

	// figure out if it's inter-res, residue_pair, or 3+body
	utility::vector1< int > pos_list( cst->residues() );

	bool success = false;
	if ( pos_list.size() == 1 ) {
		// intra-res
		tr.Trace << "remove intra-res constraint " << std::endl;
		success = remove_constraint_from_residue_constraints( pos_list[1], cst, intra_residue_constraints_, object_comparison );
	} else if ( pos_list.size() == 2 ) {
		// rsd-pai
		tr.Trace << "remove res constraint " << std::endl;
		success = remove_residue_pair_constraint( pos_list[1], pos_list[2], cst, object_comparison )
			&& remove_residue_pair_constraint( pos_list[2], pos_list[1], cst, object_comparison );
	} else {
		// 3+ body
		tr.Trace << "remove 3+body constraint " << std::endl;
		success = non_residue_pair_constraints_.remove_constraint( cst, object_comparison );
	}
	return success;
}


void
ConstraintSet::add_dof_constraint( DOF_ID const & id, func::FuncOP func, ScoreType const & t )
{
	mark_revision_id_expired();
	dof_constraints_.push_back( DOF_ConstraintOP( new DOF_Constraint( id, func, t ) ) );
}


/// @brief Returns all constraints in the set as a flat list, regardless of type.
/// @details This will be fairly inefficient if there are many, many constraints.
utility::vector1< ConstraintCOP >
ConstraintSet::get_all_constraints() const
{
	// NOTE: the old implementation which relied on Constraint::show to create a
	// string describing each constraint and then put that constraint in a
	// map from strings to ConstraintOPs had a pretty big bug: if two constraints
	// produced the same output in their show methods, then only one of them
	// would be returned by this function.
	//
	// This implementation avoids sorting constraints by their address, which can
	// produce "instabilities" if the address a constraint is allocated in varies
	// from one run to the next.  The use of std::set< ConstraintOP > would produce
	// this effect.
	utility::vector1< ConstraintCOP > all;
	for ( ResidueConstraints::const_iterator j = intra_residue_constraints_.begin(), i_end = intra_residue_constraints_.end(); j != i_end; ++j ) {
		for ( Constraints::const_iterator i = j->second->begin(), i_end = j->second->end(); i != i_end; ++i ) {
			all.push_back( *i );
		}
	}

	// since the constraint between residues i and j is held twice in the ConstraintSet,
	// (once in residue i's residue pair constraints and once in residue j's residue pair
	// constraints), we cannot blithely insert all residue pair constraints into the list,
	// but instead, should insert a constraint between i and j only when i < j (because later,
	// the constraint will appear again in the opposite order -- i.e. as between j and i).
	for ( core::Size ii = 1; ii <= residue_pair_constraints_.size(); ++ii ) {
		if ( ! residue_pair_constraints_[ ii ] ) continue; // some entries may be null
		for ( ResidueConstraints::const_iterator ij_csts = residue_pair_constraints_[ii]->begin(),
				ij_csts_end = residue_pair_constraints_[ii]->end(); ij_csts != ij_csts_end; ++ij_csts ) {
			if ( ij_csts->first < ii ) continue; // "j" is ij_csts->first; compare this against the index for i.
			for ( Constraints::const_iterator cst_iter = ij_csts->second->begin(),
					cst_iter_end = ij_csts->second->end(); cst_iter != cst_iter_end; ++cst_iter ) {
				all.push_back( *cst_iter );
			}
		}
	}

	for ( Constraints::const_iterator i = non_residue_pair_constraints_.begin(), i_end = non_residue_pair_constraints_.end(); i != i_end; ++i ) {
		all.push_back( *i );
	}

	for ( Size ii = 1; ii <= sequence_constraints_.size(); ++ii ) {
		all.push_back( sequence_constraints_[ ii ] );
	}

	return all;
}

ConstraintSet::ResiduePairConstraintsIterator
ConstraintSet::residue_pair_constraints_begin( Size resid ) const
{
	if ( residue_pair_constraints_.size() < resid || !residue_pair_constraints_[ resid ] ) return empty_rsdcst_.begin();
	return residue_pair_constraints_[ resid ]->begin();
}

ConstraintSet::ResiduePairConstraintsIterator
ConstraintSet::residue_pair_constraints_end( Size resid ) const
{
	if ( residue_pair_constraints_.size() < resid || !residue_pair_constraints_[ resid ]  ) return empty_rsdcst_.end();
	return residue_pair_constraints_[ resid ]->end();
}


void
ConstraintSet::on_length_change( conformation::signals::LengthEvent const & event )
{
	if ( ! utility::pointer::equal(conformation_pt_, event.conformation) ) {
		std::cerr << "HUH?!? weird stuff is going on. ConstraintSet is hearing length voices that it shouldn't: " << conformation_pt_.lock() << " != " << event.conformation << std::endl;
		return;
	}

	//if the signal is invalidate, let's just detach from the conformation
	if ( event.tag == conformation::signals::LengthEvent::INVALIDATE ) {

		this->detach_from_conformation();
		return;
	}


	id::SequenceMapping smap( event );

	this->remap_residue_positions( smap );
}


void
ConstraintSet::on_connection_change( core::conformation::signals::ConnectionEvent const & event )
{
	using core::conformation::signals::ConnectionEvent;

	switch ( event.tag ) {

	case ConnectionEvent::DISCONNECT :
		if ( utility::pointer::equal(conformation_pt_, event.conformation) ) {
			this->detach_from_conformation();
		} else {
			if ( !conformation_pt_.expired() ) {
				tr.Error << "ERROR: HUH?!? weird stuff is going on. ConstraintSet is hearing disconnection voices that it shouldn't" << std::endl;
			}
		}
		break;

	case ConnectionEvent::TRANSFER :
		// Disconnect -- ConstraintSet does not honor TRANSFER tag.
		break;

	default : // do nothing
		break;

	}
}

void
ConstraintSet::attach_to_conformation( core::conformation::ConformationCAP conformation ) {

	if ( !conformation_pt_.expired() ) this->detach_from_conformation();

	conformation_pt_ = conformation;

	core::conformation::ConformationCOP conformation_pt( conformation_pt_ );
	conformation_pt->attach_length_obs( &ConstraintSet::on_length_change, this );
	conformation_pt->attach_connection_obs( &ConstraintSet::on_connection_change, this );

}

void
ConstraintSet::detach_from_conformation() {

	if ( conformation_pt_.expired() ) return;

	core::conformation::ConformationCOP conformation_pt( conformation_pt_ );
	conformation_pt->detach_length_obs( &ConstraintSet::on_length_change, this );
	conformation_pt->detach_connection_obs( &ConstraintSet::on_connection_change, this );
	conformation_pt_.reset();
}


Size
ConstraintSet::revision_id() const
{
	if ( ! revision_id_current_ ) {
		++revision_id_;
		revision_id_current_ = true;
	}
	return revision_id_;
}

void
ConstraintSet::mark_revision_id_expired()
{
	revision_id_current_ = false;
}


// only prints out pair ResiduePairConstraints at the moment
void
ConstraintSet::show(
	std::ostream& out
) const {
	using namespace core::scoring::constraints;
	out << "ResiduePairConstraints: total: " << residue_pair_constraints_.size() << "   plotting active..." << std::endl;
	for ( Size ii = 1; ii <= residue_pair_constraints_.size(); ++ii ) {

		if ( ! residue_pair_constraints_exists( ii ) ) continue;

		for ( ResiduePairConstraintsIterator
				iter = residue_pair_constraints_[ ii ]->begin(),
				iter_end = residue_pair_constraints_[ ii ]->end();
				iter  != iter_end; ++iter ) {

			// iter->first is the other seqpos, iter->second is the ConstraintCOP
			if ( residue_pair_constraint_exists( ii, iter->first ) ) {
				out << "ResiduePairConstraints (" << ii << "," << iter->first << ")" << std::endl;
				iter->second->show( out );
				out << std::endl;
			}
			// print out the residue pair constraints for residue ii
		}
	} // for ( Size ii = 1; ii <= other.residue_pair_constraints_.size(); ++ii )

}

void
ConstraintSet::show_definition(
	std::ostream& out,
	pose::Pose const& pose
) const {
	using namespace core::scoring::constraints;

	// Intra-Residue
	for ( ResidueConstraints::const_iterator it = intra_residue_constraints_.begin(), eit = intra_residue_constraints_.end();
			it != eit; ++it ) {
		it->second->show_definition( out, pose );
	}

	// Residue-Pairs
	for ( Size ii = 1; ii <= residue_pair_constraints_.size(); ++ii ) {
		if ( ! residue_pair_constraints_exists( ii ) ) continue;

		for ( ResiduePairConstraintsIterator
				iter = residue_pair_constraints_[ ii ]->begin(),
				iter_end = residue_pair_constraints_[ ii ]->end();
				iter  != iter_end; ++iter ) {

			// iter->first is the other seqpos, iter->second is the ConstraintCOP
			if ( residue_pair_constraint_exists( ii, iter->first ) ) {
				if ( ii < iter->first ) {
					//out << "ResiduePairConstraints (" << ii << "," << iter->first << ")" << std::endl;
					iter->second->show_definition( out, pose );
					//  out << std::endl;
				}
			}
			// print out the residue pair constraints for residue ii
		}
	} // for ( Size ii = 1; ii <= other.residue_pair_constraints_.size(); ++ii )

	// 3+ body Constraints
	for ( Constraints::const_iterator it = non_residue_pair_constraints_.begin(),
			eit = non_residue_pair_constraints_.end(); it != eit; ++it ) {
		(*it)->show_def( out, pose );
	}
}

void
ConstraintSet::show_numbers(
	std::ostream& out
) const {
	using namespace core::scoring::constraints;

	core::Size intercount=0;
	for ( Size ii = 1; ii <= residue_pair_constraints_.size(); ++ii ) {
		if ( ! residue_pair_constraints_exists( ii ) ) continue;
		for ( ResiduePairConstraintsIterator
				iter = residue_pair_constraints_[ ii ]->begin(),
				iter_end = residue_pair_constraints_[ ii ]->end();
				iter  != iter_end; ++iter ) {

			// iter->first is the other seqpos, iter->second is the ConstraintCOP
			if ( residue_pair_constraint_exists( ii, iter->first ) ) {
				if ( ii > iter->first ) {
					intercount++;
				}
			}
		}
	} // for ( Size ii = 1; ii <= other.residue_pair_constraints_.size(); ++ii )
	out <<  "IntraRes: " << intra_residue_constraints_.size()
		<< " InterRes: " << intercount
		<<   " NonRes: " << non_residue_pair_constraints_.size()
		<< std::endl;
}


// only prints out pair ResiduePairConstraints at the moment
Size
ConstraintSet::show_violations(
	std::ostream& out,
	pose::Pose& pose,
	Size verbose_level,
	Real threshold
) const {
	using namespace core::scoring::constraints;
	Size total_viol=0;
	core::scoring::ScoreFunction empty_scorefxn;
	setup_for_scoring( pose, empty_scorefxn ); //make sure the constraints are in good shape. (eg. named->numbers)
	// Intra-Residue
	if ( intra_residue_constraints_.size() ) {
		if ( verbose_level>0 ) out << "IntraResidueConstraints: ... ";
		if ( verbose_level>50 ) out << std::endl;
		for ( ResidueConstraints::const_iterator it = intra_residue_constraints_.begin(), eit = intra_residue_constraints_.end();
				it != eit; ++it ) {
			if ( verbose_level>50 ) out << "IntraResidueConstraints ( " << it->first << " ) ";
			Size viol = it->second->show_violations( out, pose, verbose_level, threshold );
			if ( verbose_level>50 ) out << " " << viol << " violated" << std::endl;
			total_viol += viol;
		}
	}

	if ( verbose_level>0 ) out << "ResiduePairConstraints: ... ";
	if ( verbose_level>50 ) out << std::endl;
	for ( Size ii = 1; ii <= residue_pair_constraints_.size(); ++ii ) {
		if ( ! residue_pair_constraints_exists( ii ) ) continue;
		for ( ResiduePairConstraintsIterator
				iter = residue_pair_constraints_[ ii ]->begin(),
				iter_end = residue_pair_constraints_[ ii ]->end();
				iter  != iter_end; ++iter ) {

			// iter->first is the other seqpos, iter->second is the ConstraintCOP
			if ( residue_pair_constraint_exists( ii, iter->first ) ) {
				Size viol ( 0 );
				if ( ii > iter->first ) {
					if ( verbose_level>50 ) out << "ResiduePairConstraints ( " << ii << " , " << iter->first << " ) ";
					if ( verbose_level>80 ) out << std::endl;
					viol=iter->second->show_violations( out, pose, verbose_level, threshold );
					if ( verbose_level>50 ) out << " " << viol << " violated" << std::endl;
					total_viol+=viol;
				}
			}
			// print out the residue pair constraints for residue ii
		}
	} // for ( Size ii = 1; ii <= other.residue_pair_constraints_.size(); ++ii )

	// 3+ body Constraints
	if ( non_residue_pair_constraints_.size() ) {
		if ( verbose_level>0 ) out << "MultiConstraints: ... ";
		if ( verbose_level>50 ) out << std::endl;
		for ( Constraints::const_iterator it = non_residue_pair_constraints_.begin(),
				eit = non_residue_pair_constraints_.end(); it != eit; ++it ) {
			Size viol( 0 );
			if ( verbose_level>80 ) out << std::endl;
			viol=(*it)->show_violations( out, pose, verbose_level, threshold );
			if ( verbose_level>50 ) out << " " << viol << " violated" << std::endl;
			total_viol+=viol;
		}
	}
	if ( verbose_level>0 ) out << "total violations: "<<  total_viol << std::endl;
	return total_viol;
}

void
ConstraintSet::clear()
{
	sequence_constraints_.clear();
	intra_residue_constraints_.clear();
	residue_pair_constraints_.clear();
	non_residue_pair_constraints_.clear();
	dof_constraints_.clear();
	total_residue_ = 0;

	mark_revision_id_expired();

}

void
ConstraintSet::clear_sequence_constraints() {
	sequence_constraints_.clear();
	return;
}


core::Size
ConstraintSet::n_sequence_constraints() const { return sequence_constraints_.size(); }

core::scoring::aa_composition_energy::SequenceConstraintCOP
ConstraintSet::sequence_constraint( core::Size const index ) const {
	runtime_assert_string_msg( index > 0 && index <= sequence_constraints_.size(), "Error in core::scoring::constraints::ConstraintSet::sequence_constraint(): Index is out of range.  Must be greater than zero and less than or equal to the number of sequence constraints." );
	return sequence_constraints_[index];
}


bool
ConstraintSet::is_empty() const
{
	bool empty = true;
	if ( intra_residue_constraints_.size() ||
			residue_pair_constraints_.size() ||
			non_residue_pair_constraints_.size() ||
			dof_constraints_.size() ) empty = false;

	return empty;
}

#ifdef    SERIALIZATION

/// @brief Serialize this object
template < class Archive >
void
ConstraintSet::save( Archive & arc ) const
{
	arc( total_residue_ );
	arc( sequence_constraints_ );
	arc( residue_pair_constraints_ );
	arc( intra_residue_constraints_ );
	arc( non_residue_pair_constraints_ );
	arc( empty_rsdcst_ );
	arc( dof_constraints_ );
	arc( revision_id_ );
	arc( revision_id_current_ );
	// When a Pose is deserialized, it must re-establish the observation
	// of the constraint set on the conformation; this pointer cannot
	// be serialized.
	// EXEMPT conformation_pt_
}

/// @brief Deserialize this object
template < class Archive >
void
ConstraintSet::load( Archive & arc )
{
	arc( total_residue_ );
	arc( sequence_constraints_ );
	arc( residue_pair_constraints_ );
	arc( intra_residue_constraints_ );
	arc( non_residue_pair_constraints_ );
	arc( empty_rsdcst_ );
	arc( dof_constraints_ );
	arc( revision_id_ );
	arc( revision_id_current_ );

	// EXEMPT conformation_pt_
}

SAVE_AND_LOAD_SERIALIZABLE( ConstraintSet );

#endif // SERIALIZATION


std::ostream & operator << (std::ostream & os, ConstraintSet const & set)
{
	set.show(os);
	return os;
}

/// @details Copy all of the Constraint pointers from other into this so that the
/// two sets points to the same Constraint objects, but different Constraints. (Reminder,
/// a Constraint is a single constraint and a Constraints is a container of
/// Constraint objects). All of the Constraints objects must be cloned -- cloning
/// a Constraints object creates a shallow copy.  Only constraints that have at least
/// one residue between start_residue and end_residue will be copied.
void ConstraintSet::shallow_copy(
	ConstraintSet const & other,
	Size start_residue,
	Size end_residue
)
{
	basic::ProfileThis doit( basic::CONSTRAINT_SET_COPY );

	residue_pair_constraints_.clear();
	residue_pair_constraints_.resize( other.residue_pair_constraints_.size() );

	for ( Size ii = 1; ii <= other.residue_pair_constraints_.size(); ++ii ) {

		if ( ! other.residue_pair_constraints_exists( ii ) ) continue;

		bool first_insert ( true );
		for ( ResiduePairConstraintsIterator
				iter = other.residue_pair_constraints_[ ii ]->begin(),
				iter_end = other.residue_pair_constraints_[ ii ]->end();
				iter  != iter_end; ++iter ) {

			if ( ((ii < start_residue) || (ii > end_residue)) &&
					(( iter->first < start_residue) || (iter->first > end_residue )) ) {
				continue; // do not insert unless eihter end of constraint is in range
			}

			if ( first_insert ) {
				residue_pair_constraints_[ ii ] = ResidueConstraintsOP( new ResidueConstraints );
				first_insert = false;
			}
			// Shallow copy of the Constraint objects all held in a single Constraints object
			// which must be cloned
			residue_pair_constraints_[ ii ]->insert( iter->first, iter->second->clone() );
		}
	}

	intra_residue_constraints_.clear();
	for ( ResiduePairConstraintsIterator
			iter = other.intra_residue_constraints_.begin(),
			iter_end = other.intra_residue_constraints_.end();
			iter  != iter_end; ++iter ) {
		// All Constraints are now immutable, and so do not ever need to be cloned.
		//intra_residue_constraints_.insert( iter->first, iter->second->clone() );
		if ( (( iter->first < start_residue) || (iter->first > end_residue )) ) {
			continue; // do not insert unless eihter end of constraint is in range
		}
		intra_residue_constraints_.insert( iter->first, iter->second->clone() );
	}

	// Constraints assignment operator performs a shallow copy
	non_residue_pair_constraints_ = other.non_residue_pair_constraints_;

	// DOF_Constraints
	dof_constraints_ = other.dof_constraints_;

	// Shallow copy of the sequence constraints:
	sequence_constraints_.resize( other.sequence_constraints_.size() );
	for ( core::Size ii = 1, iimax = other.sequence_constraints_.size(); ii <= iimax; ++ii ) {
		if ( sequence_constraints_[ ii ] != other.sequence_constraints_[ii] ) {
			// pointer comparison -- avoid pointer copy if the two constraint sets are already
			// pointing at the same constraint object
			sequence_constraints_[ii] = other.sequence_constraints_[ ii ];
		}
	}

}

/// @details Clone all of the constraints from other into this so that the
/// two sets do not point at any shared data
void ConstraintSet::deep_copy(
	ConstraintSet const & other
)
{
	clear();
	ConstraintCOPs all_other_csts = other.get_all_constraints();
	for ( ConstraintCOPs::const_iterator iter = all_other_csts.begin(),
			iter_end = all_other_csts.end(); iter != iter_end; ++iter ) {
		add_constraint( (*iter)->clone() );
	}
}


} // constraints
} // scoring
} // core

#ifdef SERIALIZATION
CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_ConstraintSet )
#endif

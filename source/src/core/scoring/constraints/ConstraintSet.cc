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


// C++ Headers
#include <set>

#include <core/id/SequenceMapping.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/LengthEvent.hh>





namespace core {
namespace scoring {
namespace constraints {

/// @details Auto-generated virtual destructor
ResidueConstraints::~ResidueConstraints() {}

static thread_local basic::Tracer tr( "core.scoring.ConstraintSet" );

ConstraintSet::ConstraintSet()
:
	revision_id_( 0 ),
	revision_id_current_( false ),
	conformation_pt_( /* NULL */ )
{}

ConstraintSet::ConstraintSet( ConstraintSet const & other )
:
	ReferenceCount(),
	residue_pair_constraints_( other.residue_pair_constraints_.size() ),
	revision_id_( 0 ),
	revision_id_current_( false ),
	conformation_pt_( /* NULL */ )
{
	basic::ProfileThis doit( basic::CONSTRAINT_SET_COPY );

	// Loop over residue 1
	for ( Size ii = 1; ii <= other.residue_pair_constraints_.size(); ++ii ) {

		if ( ! other.residue_pair_constraints_exists( ii ) ) continue;

		bool first_insert ( true );

		// Loop over residue 2
		for ( ResiduePairConstraintsIterator
				iter = other.residue_pair_constraints_[ ii ]->begin(),
				iter_end = other.residue_pair_constraints_[ ii ]->end();
				iter  != iter_end; ++iter ) {
			if ( first_insert ) {
				residue_pair_constraints_[ ii ] = ResidueConstraintsOP( new ResidueConstraints );
				first_insert = false;
			}
			residue_pair_constraints_[ ii ]->insert( iter->first, iter->second->clone() );
		}
	}

 	for ( ResiduePairConstraintsIterator
			iter = other.intra_residue_constraints_.begin(),
			iter_end = other.intra_residue_constraints_.end();
			iter  != iter_end; ++iter ) {
		intra_residue_constraints_.insert( iter->first, iter->second->clone() );
	}

	non_residue_pair_constraints_ = other.non_residue_pair_constraints_;

	dof_constraints_ = other.dof_constraints_;

}

ConstraintSet::ConstraintSet( ConstraintSet const & other,
															Size start_residue,
															Size end_residue )
:
	ReferenceCount(),
	residue_pair_constraints_( other.residue_pair_constraints_.size() ),
	revision_id_( 0 ),
	revision_id_current_( false ),
	conformation_pt_( /* NULL */ )
{
	basic::ProfileThis doit( basic::CONSTRAINT_SET_COPY );

	for ( Size ii = 1; ii <= other.residue_pair_constraints_.size(); ++ii ) {

		if ( ! other.residue_pair_constraints_exists( ii ) ) continue;

		bool first_insert ( true );
		for ( ResiduePairConstraintsIterator
				iter = other.residue_pair_constraints_[ ii ]->begin(),
				iter_end = other.residue_pair_constraints_[ ii ]->end();
				iter  != iter_end; ++iter ) {

			if(   ((ii < start_residue) || (ii > end_residue)) &&
			      (( iter->first < start_residue) || (iter->first > end_residue )) ) {
				continue; // do not insert unless eihter end of constraint is in range
			}


			if ( first_insert ) {
				residue_pair_constraints_[ ii ] = ResidueConstraintsOP( new ResidueConstraints );
				first_insert = false;
			}
			// All Constraints are now immutable, and so do not ever need to be cloned.
			//residue_pair_constraints_[ ii ]->insert( iter->first, iter->second->clone() );
			residue_pair_constraints_[ ii ]->insert( iter->first, iter->second );
		}
	}

 	for ( ResiduePairConstraintsIterator
			iter = other.intra_residue_constraints_.begin(),
			iter_end = other.intra_residue_constraints_.end();
			iter  != iter_end; ++iter ) {
		// All Constraints are now immutable, and so do not ever need to be cloned.
		//intra_residue_constraints_.insert( iter->first, iter->second->clone() );
		if( (( iter->first < start_residue) || (iter->first > end_residue )) ) {
				continue; // do not insert unless eihter end of constraint is in range
		}
		intra_residue_constraints_.insert( iter->first, iter->second );
	}

	non_residue_pair_constraints_ = other.non_residue_pair_constraints_;

	dof_constraints_ = other.dof_constraints_;

}

ConstraintSetOP ConstraintSet::clone() const {
	return ConstraintSetOP( new ConstraintSet( *this ) );
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

	if( ! this->has_constraints() ) return; //in this case we don't have to worry about anything

	ConstraintCOPs all_cst = get_all_constraints();

	assert( all_cst.size() != 0 );

	//nuke the current constraints
	clear();

	for ( ConstraintCOPs::const_iterator it = all_cst.begin(), eit = all_cst.end(); it!=eit; ++it ) {

		ConstraintCOP new_cst = (*it)->remap_resid( smap );

		if( new_cst ) this->add_constraint( new_cst );

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
		if (!residue_pair_constraints_[i]) continue;

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
		if (!residue_pair_constraints_[i]) continue;

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

///
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


///
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

///
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

	//	Real dof_score(0.0);
	//	for ( DOF_ConstraintOPs::const_iterator it=dof_constraints_.begin(), ite = dof_constraints_.end(); it != ite; ++it ) {
	//		DOF_ConstraintOP const & dof_constraint( *it );
	//		dof_score = dof_constraint->func( pose.dof( dof_constraint->dof_id() ) );
	//		emap[ dof_constraint->score_type() ] += dof_score;
	//	}

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

// 	DOF_ConstraintOPs::const_iterator it( dof_constraints_.find( id ) ); // will this be too slow??
// 	Real deriv(0.0);
// 	if ( it != dof_constraints_.end() ) {
// 		DOF_ConstraintOP const & dof_constraint( it->second );
// 		deriv = weights[ dof_constraint->score_type() ] * dof_constraint->dfunc( pose.dof( it->first ) );
// 	}
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

///
/// private
// this function is called twice with interchanged pos1/pos2 parameters ---> all constraints are added symmetrically
void
ConstraintSet::add_residue_pair_constraint( Size const pos1, Size const pos2, ConstraintCOP cst )
{
	if (  residue_pair_constraints_.size() < pos1 ) residue_pair_constraints_.resize( pos1, 0 );
	if ( !residue_pair_constraints_[ pos1 ] ) residue_pair_constraints_[ pos1 ] = ResidueConstraintsOP( new ResidueConstraints() );

	add_constraint_to_residue_constraints( pos2, cst, *(residue_pair_constraints_[ pos1 ] ) );
}


void
ConstraintSet::add_constraints( ConstraintCOPs cst_list ) {
	for( ConstraintCOPs::iterator it = cst_list.begin();
	     it != cst_list.end();
			 it++) {
		add_constraint(*it);
	}
}

///@details copy another constraint set into this one
void
ConstraintSet::add_constraints( ConstraintSetCOP const cst_set ) {
	add_constraints(cst_set->get_all_constraints());
	return;
}

///
void
ConstraintSet::add_constraint( ConstraintCOP cst )
{
	mark_revision_id_expired();

	// figure out if it's inter-res, residue_pair, or 3+body
	utility::vector1< int > pos_list( cst->residues() );

	if ( pos_list.size() == 1 ) {
		// intra-res
		//		tr.Trace << "add intra-res constraint " << std::endl;
		add_constraint_to_residue_constraints( pos_list[1], cst, intra_residue_constraints_ );
	} else if ( pos_list.size() == 2 ) {
		// rsd-pai
		//		tr.Trace << "add res constraint " << std::endl;
		add_residue_pair_constraint( pos_list[1], pos_list[2], cst );
		add_residue_pair_constraint( pos_list[2], pos_list[1], cst );
	} else {
		// 3+ body
		//		tr.Trace << "add 3+body constraint " << std::endl;
		non_residue_pair_constraints_.add_constraint( cst );
	}
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
	if( csts_it->second->remove_constraint( cst, object_comparison ) ) {

		if( csts_it->second->size() == 0 ) residue_constraints.erase( seqpos );

		return true;
	}

	return false;
}

///
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
	if( ( (*residue_pair_constraints_[ pos1 ]).size() == 0 ) && residue_pair_constraints_.size() == pos1 ) {

		for( core::Size ii = residue_pair_constraints_.size(); ii >= 1; --ii) {

			if( residue_pair_constraints_[ ii ] ) {
				if( (*residue_pair_constraints_[ ii ]).size() != 0)	break;
			}
			residue_pair_constraints_.pop_back();
		}
	}

	return return_val;
}



bool
ConstraintSet::remove_constraints(
	ConstraintCOPs cst_list,
	bool object_comparison  )
{

	bool success = false;
	for( ConstraintCOPs::iterator it = cst_list.begin();
	     it != cst_list.end();
			 it++) {
		success = remove_constraint(*it, object_comparison);
		if( success == false ) return false;
	}

	return success;
} //remove_constraints function


///
bool
ConstraintSet::remove_constraint(
	ConstraintCOP cst,
	bool object_comparison  )
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

///
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
	// OPs implement operator< so they're OK to use in sets and maps.
	// Set takes care of duplicate insertions (from residue pair constraints).
	std::map<std::string, ConstraintCOP> all_constr;
	std::stringstream constraintString;
	for(ResidueConstraints::const_iterator j = intra_residue_constraints_.begin(), i_end = intra_residue_constraints_.end(); j != i_end; ++j) {
		for(Constraints::const_iterator i = j->second->begin(), i_end = j->second->end(); i != i_end; ++i) {
			constraintString.str("");
			(*i)->show(constraintString);
			all_constr.insert(std::make_pair(constraintString.str(),*i));
		}
	}
	for( ResiduePairConstraints::const_iterator k = residue_pair_constraints_.begin(), j_end = residue_pair_constraints_.end(); k != j_end; ++k ) {
		if( ! *k ) continue; // some entries may be null
		for(ResidueConstraints::const_iterator j = (**k).begin(), i_end = (**k).end(); j != i_end; ++j) {
			for(Constraints::const_iterator i = j->second->begin(), i_end = j->second->end(); i != i_end; ++i) {
				constraintString.str("");
				(*i)->show(constraintString);
				all_constr.insert(std::make_pair(constraintString.str(),*i));
			}
		}
	}
	for(Constraints::const_iterator i = non_residue_pair_constraints_.begin(), i_end = non_residue_pair_constraints_.end(); i != i_end; ++i) {
		constraintString.str("");
		(*i)->show(constraintString);
		all_constr.insert(std::make_pair(constraintString.str(),*i));
	}
	// Copy final set contents into a list to return to user...
	utility::vector1< ConstraintCOP > all;
	for(std::map<std::string, ConstraintCOP>::iterator i = all_constr.begin(), i_end = all_constr.end(); i != i_end; ++i) {
		all.push_back(i->second);
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
ConstraintSet::on_length_change( conformation::signals::LengthEvent const & event ) {

	if( ! utility::pointer::equal(conformation_pt_, event.conformation) ) {
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
ConstraintSet::on_connection_change( core::conformation::signals::ConnectionEvent const & event ) {
	using core::conformation::signals::ConnectionEvent;

	switch ( event.tag ) {

		case ConnectionEvent::DISCONNECT:
			if( utility::pointer::equal(conformation_pt_, event.conformation) ) {
				this->detach_from_conformation();
			} else {
				tr.Error << "ERROR: HUH?!? weird stuff is going on. ConstraintSet is hearing disconnection voices that it shouldn't" << std::endl;
			}
			break;

		case ConnectionEvent::TRANSFER:
			// Disconnect -- ConstraintSet does not honor TRANSFER tag.
			break;

		default: // do nothing
			break;

	}
}

void
ConstraintSet::attach_to_conformation( core::conformation::ConformationCAP conformation ) {

	if( !conformation_pt_.expired() ) this->detach_from_conformation();

	conformation_pt_ = conformation;

	core::conformation::ConformationCOP conformation_pt( conformation_pt_ );
	conformation_pt->attach_length_obs( &ConstraintSet::on_length_change, this );
	conformation_pt->attach_connection_obs( &ConstraintSet::on_connection_change, this );

}

void
ConstraintSet::detach_from_conformation() {

	if( conformation_pt_.expired() ) return;

#ifdef PTR_MODERN
	core::conformation::ConformationCOP conformation_pt( conformation_pt_ );
	conformation_pt->detach_length_obs( &ConstraintSet::on_length_change, this );
	conformation_pt->detach_connection_obs( &ConstraintSet::on_connection_change, this );
	conformation_pt_.reset();
#else
	// With ReferenceCount, this gets interesting:
	// This function is called from Conformation::~Conformation via
	// notify_connection_obs() and then on_connection_change().
	// The conformation_pt_->count_ is 0 at this point.
	// Putting conformation_pt_ into an OP (as above) increments the
	// counter again, which is decremented at the end of this function
	// when conformation_pt expires and Conformation::~Conformation
	// is called again, resulting in a call loop and stack overflow.
	conformation_pt_->detach_length_obs( &ConstraintSet::on_length_change, this );
	conformation_pt_->detach_connection_obs( &ConstraintSet::on_connection_change, this );
	conformation_pt_ = NULL;
#endif

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
					//		out << std::endl;
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
		if ( verbose_level>0) out << "IntraResidueConstraints: ... ";
		if ( verbose_level>50 ) out << std::endl;
		for ( ResidueConstraints::const_iterator it = intra_residue_constraints_.begin(), eit = intra_residue_constraints_.end();
					it != eit; ++it ) {
			if ( verbose_level>50 ) out << "IntraResidueConstraints ( " << it->first << " ) ";
			Size viol = it->second->show_violations( out, pose, verbose_level, threshold );
			if ( verbose_level>50 ) out << " " << viol << " violated" << std::endl;
			total_viol += viol;
		}
	}

	if ( verbose_level>0) out << "ResiduePairConstraints: ... ";
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
		if ( verbose_level>0) out << "MultiConstraints: ... ";
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
	intra_residue_constraints_.clear();
	residue_pair_constraints_.clear();
	non_residue_pair_constraints_.clear();
	dof_constraints_.clear();

	mark_revision_id_expired();

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


std::ostream & operator << (std::ostream & os, ConstraintSet const & set)
{
	set.show(os);
    return os;
}


} // constraints
} // scoring
} // core

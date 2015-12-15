// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file
/// @brief contains declarations for a type of constraint that holds a number of constraints
/// @brief where only the lowest N are evaluated
/// @author Florian Richter (floric@u.washington.edu, march 2008)


#include <protocols/constraints_additional/AmbiguousMultiConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <utility/sort_predicates.hh>

//C++ headers
#include <list>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace constraints_additional {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
AmbiguousMultiConstraint::AmbiguousMultiConstraint(
	core::Size num_act_csts
):
	AmbiguousConstraint(),
	num_active_constraints_(num_act_csts)
{
	active_constraints_.clear();
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
AmbiguousMultiConstraint::AmbiguousMultiConstraint(
	core::Size num_act_csts,
	core::scoring::constraints::ConstraintCOPs & cst_in ):
	AmbiguousConstraint( cst_in ),
	num_active_constraints_(num_act_csts)
{
	active_constraints_.clear();
}

bool AmbiguousMultiConstraint::operator == ( core::scoring::constraints::Constraint const & other ) const
{
	if ( ! core::scoring::constraints::AmbiguousConstraint::operator == ( other ) ) return false;

	AmbiguousMultiConstraint const & other_downcast( static_cast< AmbiguousMultiConstraint const & > ( other ) );

	return num_active_constraints_ == other_downcast.num_active_constraints_;
}

bool AmbiguousMultiConstraint::same_type_as_me( core::scoring::constraints::Constraint const & other ) const
{
	return dynamic_cast< AmbiguousMultiConstraint const * > (&other);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief ScoreFunction, scores all member constraints but only reports the lowest N ones
/// @brief note: this could potentially be made faster if the cur_emap isn't copied, but instead
/// @brief pushed back right away into score_cst_pairs and then score_cst_pairs.last() is handed
/// @brief to the member_it->score() function
void
AmbiguousMultiConstraint::score(
	core::scoring::func::XYZ_Func const & xyz_func,
	core::scoring::EnergyMap const & weights,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core::scoring::constraints;
	using namespace core::scoring;

	runtime_assert( member_constraints().size() >= num_active_constraints_ );

	typedef std::pair< ConstraintCOP, EnergyMap > cst_emap_pair;
	typedef std::pair< core::Real, cst_emap_pair > score_cst_pair;
	//std::cout << "scoring ambiguous constraint..." << std::endl;

	std::list< score_cst_pair > score_cst_pairs;

	//low_total_cst_score_ = 1000000;
	//low_EMap_.zero( cst_score_types_ );

	// first, we score every constraint and save the result in a list
	for ( ConstraintCOPs::const_iterator
			member_it = member_constraints().begin(), end = member_constraints().end(); member_it != end; ++member_it ) {
		EnergyMap cur_emap;
		(*member_it)->score(xyz_func, weights, cur_emap);
		core::Real cur_score = calculate_total_cst_score( weights, cur_emap);
		score_cst_pairs.push_back( score_cst_pair( cur_score, cst_emap_pair( *member_it, cur_emap ) ) );
	}

	//then we sort the list by score
	score_cst_pairs.sort( utility::SortFirst< core::Real, cst_emap_pair >() );

	active_constraints_.clear();
	//and go through and add the information of the relevant constraints to the emap,
	//as well as note the active constraints
	core::Size counter(1);
	for ( std::list< score_cst_pair >::iterator cst_it = score_cst_pairs.begin() ;
			counter <= num_active_constraints_; ++cst_it ) {


		emap[constant_constraint] += cst_it->second.second[constant_constraint];
		emap[coordinate_constraint] += cst_it->second.second[coordinate_constraint];
		emap[atom_pair_constraint] += cst_it->second.second[atom_pair_constraint];
		emap[angle_constraint] += cst_it->second.second[angle_constraint];
		emap[dihedral_constraint] += cst_it->second.second[dihedral_constraint];
		emap[backbone_stub_constraint] += cst_it->second.second[backbone_stub_constraint];

		active_constraints_.push_back( cst_it->second.first );

		//don't forget to increment the counter
		counter++;
	}

} //score


core::scoring::constraints::ConstraintOP
AmbiguousMultiConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	using namespace core::scoring::constraints;

	ConstraintCOPs new_csts;
	for ( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); ++cst_it ) {
		ConstraintOP new_cst = (*cst_it)->remap_resid( seqmap );
		if ( new_cst ) new_csts.push_back( new_cst );
	}

	if ( new_csts.size() > 0 ) {
		return ConstraintOP( new AmbiguousMultiConstraint( num_active_constraints_, new_csts ) );
	} else return NULL;

}


/// @brief function to minimize the N lowest scoring member constraints
void
AmbiguousMultiConstraint::fill_f1_f2(
	core::id::AtomID const & atom,
	core::scoring::func::XYZ_Func const & xyz,
	core::Vector & F1,
	core::Vector & F2,
	core::scoring::EnergyMap const & weights
) const
{
	using namespace core::scoring::constraints;

	runtime_assert( active_constraints_.size() == num_active_constraints_ );
	for ( ConstraintCOPs::const_iterator cst_it = active_constraints_.begin(); cst_it != active_constraints_.end(); ++cst_it ) {
		assert( *cst_it ) ;
		(*cst_it)->fill_f1_f2(atom, xyz, F1, F2, weights);
	}
}

void
AmbiguousMultiConstraint::show( std::ostream& out) const
{
	using namespace core::scoring::constraints;
	out << "AmbiguousMultiConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
	for ( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(), end = member_constraints().end(); cst_it != end; ++cst_it ) {
		(*cst_it)->show(out);
	}

	out << " ...all member constraints of this AmbiguousMultiConstraint shown." << std::endl;
}

core::Size
AmbiguousMultiConstraint::show_violations( std::ostream& out, core::pose::Pose const& pose, core::Size verbose_level, core::Real threshold ) const
{
	using namespace core::scoring::constraints;

	if ( verbose_level >= 70 ) {
		Size total_viol( 0 );
		if ( verbose_level >=80 ) out << type() << " " << member_constraints().size() << " ";
		for ( ConstraintCOPs::const_iterator cst_it = active_constraints_.begin(), end = active_constraints_.end(); cst_it != end; ++cst_it ) {

			//out << "active ";
			total_viol += (*cst_it)->show_violations( out, pose, verbose_level,   threshold );
		}
		return total_viol;
	}
	return 0;
}

}
}

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::constraints_additional::AmbiguousMultiConstraint::AmbiguousMultiConstraint() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::constraints_additional::AmbiguousMultiConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::constraints::AmbiguousConstraint >( this ) );
	arc( CEREAL_NVP( num_active_constraints_ ) ); // core::Size
	arc( CEREAL_NVP( active_constraints_ ) ); // core::scoring::constraints::ConstraintCOPs
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::constraints_additional::AmbiguousMultiConstraint::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::constraints::AmbiguousConstraint >( this ) );
	arc( num_active_constraints_ ); // core::Size
	utility::vector1< std::shared_ptr< core::scoring::constraints::Constraint > > local_active_constraints;
	arc( local_active_constraints ); // core::scoring::constraints::ConstraintCOPs
	active_constraints_ = local_active_constraints; // copy the non-const pointer(s) into the const pointer(s)
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::constraints_additional::AmbiguousMultiConstraint );
CEREAL_REGISTER_TYPE( protocols::constraints_additional::AmbiguousMultiConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_constraints_additional_AmbiguousMultiConstraint )
#endif // SERIALIZATION

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief contains declarations for a type of constraint that holds multiple other constrains that belong to each other and are all evaluate at once
/// @author Florian Richter (floric@u.washington.edu, march 2008)

#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/id/SequenceMapping.hh>
#include <numeric/random/random.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer tr( "core.io.constraints" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/// @brief default Constructor
MultiConstraint::MultiConstraint( ScoreType const & t /* = dof_constraint */ ) :
	Constraint( t ),
	report_this_as_effective_sequence_separation_( 0 )
{}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
MultiConstraint::MultiConstraint( ConstraintCOPs const & cst_in, ScoreType const & t ):
	Constraint( t ),  //this is temporary, multi constraint shouldn't have a score type
	report_this_as_effective_sequence_separation_( 0 )
{
	for ( ConstraintCOPs::const_iterator it = cst_in.begin(); it != cst_in.end(); ++it ) {
		add_individual_constraint( *it );
	} //loop over all input csts that make up this multi constraint
}// constructor



/// @details Clone all of the member constraints into a new vector and then invoke the
/// (shallow-copy) constructor that accepts a vector of ConstriantCOPs -- if this
/// %MultiConstraint doesn't have any member constraints, then return a default-constructed
/// one.
ConstraintOP
MultiConstraint::clone() const {
	return ConstraintOP( new MultiConstraint( *this ) );
}


MultiConstraintOP
MultiConstraint::empty_clone() const {
	return MultiConstraintOP( new MultiConstraint );
}

/// @brief number of atoms involved in this MultiConstraint container
Size
MultiConstraint::natoms() const
{
	return member_atoms_.size();
}

/// @brief number of constraints data
Size
MultiConstraint::size() const { return member_constraints_.size(); }

std::string
MultiConstraint::type() const {
	return "MultiConstraint";
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details read definition of a multiconstraint. Since a MultiConstraint is essentially a vector of ?????? (Please someone finish this)
void
MultiConstraint::read_def(
	std::istream& data,
	core::pose::Pose const& pose,
	func::FuncFactory const & func_factory
)
{
	Size ct ( 0 );
	ConstraintCOP cst_op = ConstraintIO::get_instance()->read_individual_constraint_new( data, pose, func_factory );
	if ( cst_op ) ++ct;
	while ( cst_op ) {
		add_individual_constraint( cst_op );
		cst_op = ConstraintIO::get_instance()->read_individual_constraint_new( data, pose, func_factory );
		if ( cst_op ) ++ct;
	}
	if ( member_constraints_.size() > 0 ) {
		if ( tr.Debug.visible() ) {
			show_def(tr.Debug,pose);
			tr.Debug << std::endl;
		}
	} else {
		tr.Error << type() << " read_def: no constraints defined! " << ct << " " << size() << std::endl;
	}
}

/// @detail this function only checks whether the member_constraints_
/// are identical. should mean by inference that the member_residues_, member_atoms_
/// and atomid_to_csts_ are also identical
bool
MultiConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( !           same_type_as_me( other_cst ) ) return false;
	if ( ! other_cst.same_type_as_me(     *this ) ) return false;

	MultiConstraint const & other( static_cast< MultiConstraint const & > (other_cst) );

	if ( report_this_as_effective_sequence_separation_ != other.report_this_as_effective_sequence_separation_ ) return false;
	if ( member_constraints_.size() != other.member_constraints_.size() ) return false;
	for ( core::Size i =1; i <= member_constraints_.size(); ++i ) {
		if ( *(member_constraints_[i]) != *(other.member_constraints_[i]) ) return false;
	}
	if ( this->score_type() != other.score_type() ) return false;

	return true;
}

bool MultiConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< MultiConstraint const * > ( &other );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultiConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{
	for ( auto const & member_constraint : member_constraints_ ) {
		member_constraint->score( xyz_func, weights, emap );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultiConstraint::add_individual_constraint( ConstraintCOP cst_in )
{
	member_constraints_.push_back( cst_in );

	//examine this constraint with respect to whether it contains
	//not yet seen atoms
	for ( Size i = 1; i <= cst_in->natoms(); ++i ) {
		AtomID cur_atomid = cst_in->atom(i);

		if ( std::find( member_residues_.begin(), member_residues_.end(), cur_atomid.rsd() ) == member_residues_.end() ) {
			member_residues_.push_back( cur_atomid.rsd() );
		}

		auto map_it = atomid_to_csts_.find(cur_atomid);

		//if it does, add the atom id to the atom_members and the atomid/vect pair to the map
		if ( map_it == atomid_to_csts_.end() ) {
			member_atoms_.push_back(cur_atomid);
			ConstraintCOPs cst_vect;
			cst_vect.push_back( cst_in );

			atomid_to_csts_.insert( std::pair< AtomID, ConstraintCOPs > (cur_atomid, cst_vect) );
		} else { //if it doesn't we have to add this constraint to the right list in the map
			map_it->second.push_back( cst_in );
		}
	}
}

//@brief translates the atom-names into numbers
void MultiConstraint::setup_for_scoring( func::XYZ_Func const & xyz, ScoreFunction const & scfxn) const {
	for ( auto const & cst : member_constraints_ ) {
		cst->setup_for_scoring( xyz, scfxn );
	}
}


ConstraintOP
MultiConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	ConstraintCOPs new_csts;
	for ( auto const & cst : member_constraints_ ) {
		ConstraintOP new_cst = cst->remap_resid( seqmap );
		if ( new_cst ) new_csts.push_back( new_cst );
	}
	if ( new_csts.size() > 0 ) {
		return ConstraintOP( new MultiConstraint( new_csts ) );
	} else return NULL;
}

ConstraintOP MultiConstraint::remapped_clone(
	pose::Pose const& src,
	pose::Pose const& dest,
	id::SequenceMappingCOP map ) const
{
	MultiConstraintOP new_multi = empty_clone();
	for ( auto const & cst : member_constraints_ ) {
		new_multi->add_individual_constraint( cst->remapped_clone( src, dest, map ) );
	}
	return new_multi;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief function that figures out what constraints this atom is part of
/// and calculates the derivative for those
void
MultiConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{

	//std::cout << "we are in a multicst containing the following atoms (test AtomID is " << atom.rsd() << " " << atom.atomno() << ": ";
	//  for(Size i = 1; i <= member_atoms_.size(); ++i){
	//std::cout << member_atoms_[i].rsd() << " " << member_atoms_[i].atomno() << " ";
	// }
	//std::cout << ", the atom map has " << atomid_to_csts_.size() << " atoms." << std::endl;

	auto map_it = atomid_to_csts_.find(atom);

	if ( map_it != atomid_to_csts_.end() ) {
		ConstraintCOPs cur_csts = map_it->second;
		//std::cout << "now taking deriv of atom " << map_it->first.atomno() << " , vector 1 element before is " << F1[1];

		for ( ConstraintCOPs::const_iterator cst_it = cur_csts.begin(), end = cur_csts.end(); cst_it != end; ++cst_it ) {
			//std::cout << ", type of cst is " << (*cst_it)->score_type() << ",  ";
			(*cst_it)->fill_f1_f2(atom, xyz, F1, F2, weights);
		}
		//std::cout << " and after " << F1[1] << std::endl;

	} else return;
}

void
MultiConstraint::show( std::ostream& out) const
{
	out << "MultiConstraint containing the following " << member_constraints_.size() << " constraints: "
		<< std::endl;
	for ( ConstraintCOPs::const_iterator cst_it = member_constraints_.begin(), end = member_constraints_.end(); cst_it != end; ++cst_it ) {
		(*cst_it)->show(out);
	}

	out << " ...all member constraints of this MultiConstraint shown." << std::endl;
}

void MultiConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	out << type() << std::endl;
	for ( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(), end = member_constraints().end(); cst_it != end; ++cst_it ) {
		(*cst_it)->show_def( out, pose );
	}
	out << "End_"<< type() << std::endl;
}

Size
MultiConstraint::show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold ) const
{
	if ( verbose_level > 80 ) out << "Violations for MultiConstraint: " << std::endl;

	core::Size biggest_violation(0);
	for ( ConstraintCOPs::const_iterator cst_it = member_constraints_.begin(), end = member_constraints_.end(); cst_it != end; ++cst_it ) {
		core::Size cur_viol = (*cst_it)->show_violations( out, pose, verbose_level, threshold);
		if ( cur_viol > biggest_violation ) biggest_violation = cur_viol;
	}
	return biggest_violation;
}

core::Size
MultiConstraint::choose_effective_sequence_separation(
	core::kinematics::ShortestPathInFoldTree const& sp,
	numeric::random::RandomGenerator& RG
) {
	utility::vector1< core::Size > collected_seq_separations;
	ConstraintCOPs const& cst_list( member_constraints() );
	for ( ConstraintCOPs::const_iterator cst_it = cst_list.begin(), end = cst_list.end(); cst_it != end; ++cst_it ) {
		Size seq_sep( 0 );
		Constraint& cst = const_cast< Constraint& >( **cst_it );
		seq_sep = cst.choose_effective_sequence_separation( sp, RG );
		collected_seq_separations.push_back( seq_sep );
	}
	Size index = static_cast< Size >( collected_seq_separations.size()*RG.uniform() ) + 1;
	report_this_as_effective_sequence_separation_ = collected_seq_separations[ index ];
	return report_this_as_effective_sequence_separation_;
}

MultiConstraint::MultiConstraint( MultiConstraint const & src ) :
	Constraint( src ),
	report_this_as_effective_sequence_separation_( src.report_this_as_effective_sequence_separation_ )
{
	for ( ConstraintCOPs::const_iterator it = src.member_constraints_.begin(); it != src.member_constraints_.end(); ++it ) {
		add_individual_constraint( (*it)->clone() );
	}
}

/// @brief Return a vector of Constraints that are clones of the member constraints.
ConstraintCOPs
MultiConstraint::cloned_member_constraints() const
{
	ConstraintCOPs member_constraint_clones( member_constraints_.size() );
	for ( Size ii = 1; ii <= member_constraints_.size(); ++ii ) {
		member_constraint_clones[ ii ] = ConstraintCOP( member_constraints_[ii]->clone() );
	}
	return member_constraint_clones;
}

} //constraints
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::MultiConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( member_constraints_ ) ); // ConstraintCOPs
	arc( CEREAL_NVP( member_residues_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( member_atoms_ ) ); // utility::vector1<AtomID>
	arc( CEREAL_NVP( atomid_to_csts_ ) ); // std::map<AtomID, ConstraintCOPs>
	arc( CEREAL_NVP( report_this_as_effective_sequence_separation_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::MultiConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	utility::vector1< std::shared_ptr< core::scoring::constraints::Constraint > > local_member_constraints;
	arc( local_member_constraints ); // ConstraintCOPs
	member_constraints_ = local_member_constraints; // copy the non-const pointer(s) into the const pointer(s)
	arc( member_residues_ ); // utility::vector1<core::Size>
	arc( member_atoms_ ); // utility::vector1<AtomID>

	std::map< AtomID, core::scoring::constraints::ConstraintOPs > local_atomid_to_csts;
	arc( local_atomid_to_csts );
	for ( std::map< AtomID, ConstraintOPs >::const_iterator
			iter = local_atomid_to_csts.begin(), iter_end = local_atomid_to_csts.end();
			iter != iter_end; ++iter ) {
		atomid_to_csts_[ iter->first ] = iter->second;
	}
	// atomid_to_csts_ = local_atomid_to_csts; // copy the non-const pointer(s) into the const pointer(s)

	arc( report_this_as_effective_sequence_separation_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::MultiConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::MultiConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_MultiConstraint )
#endif // SERIALIZATION

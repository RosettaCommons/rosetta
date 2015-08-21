// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

static thread_local basic::Tracer tr( "core.io.constraints" );

namespace core {
namespace scoring {
namespace constraints {

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
	if ( cst_op ) ct++;
	while ( cst_op ) {
		add_individual_constraint( cst_op );
		cst_op = ConstraintIO::get_instance()->read_individual_constraint_new( data, pose, func_factory );
		if ( cst_op ) ct++;
	}
	if ( member_constraints_.size() > 0 ) {
		if ( tr.Debug.visible() ) {
			show_def(tr.Debug,pose);
			tr.Debug << std::endl;
		}
	} else {
		tr.Error << "ERROR: " << type() << " read_def: no constraints defined! " << ct << " " << size() << std::endl;
	}
}

/// @detail this function only checks whether the member_constraints_
/// are identical. should mean by inference that the member_residues_, member_atoms_
/// and AtomID_to_Csts_ are also identical
bool
MultiConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( !dynamic_cast< MultiConstraint const * > ( &other_cst ) ) return false;

	MultiConstraint const & other( static_cast< MultiConstraint const & > (other_cst) );

	if ( report_this_as_effective_sequence_separation_ != other.report_this_as_effective_sequence_separation_ ) return false;

	core::Size num_csts = member_constraints_.size();
	if ( num_csts != other.member_constraints_.size() ) return false;
	for ( core::Size i =1; i <= num_csts; ++i ) {
		if ( *(member_constraints_[i]) != *(other.member_constraints_[i]) ) return false;
	}
	if ( this->score_type() != other.score_type() ) return false;

	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultiConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{
	for ( ConstraintCOPs::const_iterator member_it = member_constraints_.begin(), end = member_constraints_.end(); member_it != end; ++member_it ) {
		(*member_it)->score( xyz_func, weights, emap );
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
MultiConstraint::MultiConstraint( const ConstraintCOPs & cst_in, ScoreType const & t ):
	Constraint( t ),  //this is temporary, multi constraint shouldn't have a score type
	report_this_as_effective_sequence_separation_( 0 )
{
	for ( ConstraintCOPs::const_iterator it = cst_in.begin(), end = cst_in.end(); it != end; ++it ) {
		add_individual_constraint( *it );
	} //loop over all input csts that make up this multi constraint
}// constructor
////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultiConstraint::add_individual_constraint( ConstraintCOP cst_in )
{
	member_constraints_.push_back( cst_in );

	//examine this constraint with respect to whether it contains
	//not yet seen atoms
	for ( Size i = 1; i <= cst_in->natoms(); i++ ) {
		AtomID cur_atomid = cst_in->atom(i);

		if ( std::find( member_residues_.begin(), member_residues_.end(), cur_atomid.rsd() ) == member_residues_.end() ) {
			member_residues_.push_back( cur_atomid.rsd() );
		}

		std::map< AtomID, ConstraintCOPs >::iterator map_it = AtomID_to_Csts_.find(cur_atomid);

		//if it does, add the atom id to the atom_members and the atomid/vect pair to the map
		if ( map_it == AtomID_to_Csts_.end() ) {
			member_atoms_.push_back(cur_atomid);
			ConstraintCOPs cst_vect;
			cst_vect.push_back( cst_in );

			AtomID_to_Csts_.insert( std::pair< AtomID, ConstraintCOPs > (cur_atomid, cst_vect) );
		} else { //if it doesn't we have to add this constraint to the right list in the map
			map_it->second.push_back( cst_in );
		}
	}
}

//@brief translates the atom-names into numbers
void MultiConstraint::setup_for_scoring( func::XYZ_Func const & xyz, ScoreFunction const & scfxn) const {
	for ( ConstraintCOPs::const_iterator cst_it = member_constraints_.begin(); cst_it != member_constraints_.end(); ++cst_it ) {
		(*cst_it)->setup_for_scoring( xyz, scfxn );
	}
}


ConstraintOP
MultiConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	ConstraintCOPs new_csts;
	for ( ConstraintCOPs::const_iterator cst_it = member_constraints_.begin(); cst_it != member_constraints_.end(); ++cst_it ) {
		ConstraintOP new_cst = (*cst_it)->remap_resid( seqmap );
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
	for ( ConstraintCOPs::const_iterator it = member_constraints_.begin();
			it != member_constraints_.end(); ++it ) {
		new_multi->add_individual_constraint( (*it)->remapped_clone( src, dest, map ) );
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
	//  for(Size i = 1; i <= member_atoms_.size(); i++){
	//std::cout << member_atoms_[i].rsd() << " " << member_atoms_[i].atomno() << " ";
	// }
	//std::cout << ", the atom map has " << AtomID_to_Csts_.size() << " atoms." << std::endl;

	std::map< AtomID, ConstraintCOPs >::const_iterator map_it = AtomID_to_Csts_.find(atom);

	if ( map_it != AtomID_to_Csts_.end() ) {
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


} //constraints
} //scoring
} //core

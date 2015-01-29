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
/// @brief where only the lowest one is evaluated
/// @author Florian Richter (floric@u.washington.edu, march 2008)

// Unit headers
#include <core/scoring/constraints/AmbiguousConstraint.hh>

// Package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
// AUTO-REMOVED #include <core/scoring/func/XYZ_Func.hh>

// AUTO-REMOVED #include <core/scoring/Energies.hh>

// Project headers
// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
AmbiguousConstraint::AmbiguousConstraint():
	MultiConstraint()
{
	active_constraint_ = NULL;

	init_cst_score_types();

debug_assert ( member_constraints().size() == 0 );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
AmbiguousConstraint::AmbiguousConstraint( ConstraintCOPs & cst_in ):
	MultiConstraint( cst_in )
{
	active_constraint_ = NULL;

	init_cst_score_types();

debug_assert ( member_constraints().size() > 0 );

}
////////////////////////////////////////////////////////////////////////////////////////////////////
void
AmbiguousConstraint::init_cst_score_types()
{
	cst_score_types_.clear();
	cst_score_types_.push_back(constant_constraint);
	cst_score_types_.push_back(coordinate_constraint);
	cst_score_types_.push_back(atom_pair_constraint);
	cst_score_types_.push_back(angle_constraint);
	cst_score_types_.push_back(dihedral_constraint);
	cst_score_types_.push_back(backbone_stub_constraint);
	cst_score_types_.push_back(backbone_stub_linear_constraint);
}


/// @detail this function only checks whether the passed in cst
/// is the same type and then hands off to the multi cst. probably
/// overkill to also demand that the active cst is identical
bool
AmbiguousConstraint::operator == ( Constraint const & other_cst ) const
{
	if( !dynamic_cast< AmbiguousConstraint const * > ( &other_cst ) ) return false;

	return MultiConstraint::operator==( other_cst);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief ScoreFunction, scores all member constraints but only reports the lowest one
void
AmbiguousConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{

	Real low_total_cst_score = 0;
	low_EMap_.zero( cst_score_types_ );

	core::Real cur_score = 0;
	bool first_pass = true;

	for( ConstraintCOPs::const_iterator
			member_it = member_constraints().begin();
			member_it != member_constraints().end(); member_it++){

		temp_EMap_.zero(cst_score_types_ );

		(*member_it)->score(xyz_func, weights, temp_EMap_ );

		cur_score = calculate_total_cst_score( weights, temp_EMap_ );

		if ( first_pass || (cur_score < low_total_cst_score ) ){
			first_pass = false;

			low_total_cst_score = cur_score;
			active_constraint_ = *member_it;

			low_EMap_[constant_constraint]      = temp_EMap_[constant_constraint];
			low_EMap_[coordinate_constraint]    = temp_EMap_[coordinate_constraint];
			low_EMap_[atom_pair_constraint]     = temp_EMap_[atom_pair_constraint];
			low_EMap_[angle_constraint]         = temp_EMap_[angle_constraint];
			low_EMap_[dihedral_constraint]      = temp_EMap_[dihedral_constraint];
			low_EMap_[backbone_stub_constraint] = temp_EMap_[backbone_stub_constraint];
			low_EMap_[backbone_stub_linear_constraint] = temp_EMap_[backbone_stub_linear_constraint];

		}

	}
	//add lowest score to the actual emap
	emap[constant_constraint]      += low_EMap_[constant_constraint];
	emap[coordinate_constraint]    += low_EMap_[coordinate_constraint];
	emap[atom_pair_constraint]     += low_EMap_[atom_pair_constraint];
	emap[angle_constraint]         += low_EMap_[angle_constraint];
	emap[dihedral_constraint]      += low_EMap_[dihedral_constraint];
	emap[backbone_stub_constraint] += low_EMap_[backbone_stub_constraint];
	emap[backbone_stub_linear_constraint] += low_EMap_[backbone_stub_linear_constraint];


} //score


/// @brief helper function to accumulate all constraint scores into one number
core::Real
AmbiguousConstraint::calculate_total_cst_score( EnergyMap const & weights, EnergyMap & emap ) const
{

	core::Real total_score = emap[constant_constraint] * weights[constant_constraint] + emap[coordinate_constraint] * weights[coordinate_constraint] + emap[atom_pair_constraint] * weights[atom_pair_constraint] + emap[angle_constraint] * weights[angle_constraint] + emap[dihedral_constraint] * weights[dihedral_constraint] + emap[backbone_stub_constraint] * weights[backbone_stub_constraint] + emap[backbone_stub_linear_constraint] * weights[backbone_stub_linear_constraint];

	return total_score;
}


ConstraintOP
AmbiguousConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{

	ConstraintCOPs new_csts;

	for( ConstraintCOPs::const_iterator cst_it = member_constraints_.begin(); cst_it != member_constraints_.end(); ++cst_it ){

		ConstraintOP new_cst = (*cst_it)->remap_resid( seqmap );

		if( new_cst ) new_csts.push_back( new_cst );

	}

	if( new_csts.size() > 0 ){
		return ConstraintOP( new AmbiguousConstraint( new_csts ) );
	}
	else return NULL;

}


/// @brief function to minimize lowest scoring member constraint
void
AmbiguousConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
//debug_assert( active_constraint_ );
	//std::cerr << "Error: attempted minimization of an ambiguous constraint before the lowest scoring member constraint has been determined. Was the pose not scored?" << std::endl;
	//utility::exit( EXIT_FAILURE, __FILE__, __LINE__);

	// fpd
	if (member_constraints_.size() == 0)
		return;

	// fpd
	if ( !active_constraint_ ) {
		EnergyMap dummy;
		score( xyz, weights, dummy);
	}

	active_constraint_->fill_f1_f2( atom, xyz, F1, F2, weights );
}

void
AmbiguousConstraint::show( std::ostream& out) const
{
	/// APL -- you cannot show the active constraint in the absence of a Pose.
	//out << "AmbiguousConstraint Active constraint:" << std::endl;
	//active_constraint()->show(out);
	out << "AmbiguousConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
	for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
		(*cst_it)->show(out);
	}

	out << " ...all member constraints of this AmbiguousConstraint shown." << std::endl;
}

ConstraintCOP
AmbiguousConstraint::active_constraint(
//	func::XYZ_Func const & xyz, // these arguments would allow active-constraint determination in a purley const fashion
//	EnergyMap const & weights
) const
{
	return active_constraint_;
	/* APL Code below is for when we remove non-constness to the Ambiguous constraints
	ConstraintCOP active;

	Real low_total_cst_score = 0;
	EnergyMap temp_EMap;
	core::Real cur_score = 0;
	bool first_pass = true;

	for( ConstraintCOPs::const_iterator
			member_it = member_constraints().begin();
			member_it != member_constraints().end(); member_it++){
		temp_EMap.zero(cst_score_types_ );
		(*member_it)->score( xyz, weights, temp_EMap );
		cur_score = calculate_total_cst_score( weights, temp_EMap );
		if ( first_pass || (cur_score < low_total_cst_score ) ){
			first_pass = false;
			low_total_cst_score = cur_score;
			active = *member_it;
		}
	}

	return active;*/
}

Size
AmbiguousConstraint::show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold ) const
{
	if ( verbose_level >= 70 ) {
		Size total_viol( 0 );
		if ( verbose_level >=80 ) { out << type() << " " << member_constraints().size() << " "; }
		for( ConstraintCOPs::const_iterator
				cst_it = member_constraints().begin();
				cst_it != member_constraints().end(); cst_it++) {
			if ( active_constraint_ ) {
				if ( (*cst_it).get() == active_constraint_.get() ) {
					// figure out if it's inter-res, residue_pair, or 3+body
					utility::vector1< int > pos_list( (*cst_it)->residues() );
					if ( pos_list.size() == 2 ) {
						out << "ResiduePairConstraints ( " << pos_list[ 1 ] << " , " << pos_list[ 2 ] << " ) ";
						if ( verbose_level>80 ) out << std::endl;
						//out << "active ";
						total_viol += (*cst_it)->show_violations( out, pose, verbose_level, threshold );
						out << " of total: " << 1 << " ";
					}
				} else {
					//					out << "       ";
					//don't count the violations from this one
					//					(*cst_it)->show_violations( out, pose, verbose_level, threshold );
				}
			}
		}
		return total_viol;
	} else {
		if( active_constraint_ ) return active_constraint_->show_violations( out, pose, verbose_level );
		else out << "WARNING: requested to show violations of an ambiguous constraint before"
			<<" the pose was scored and thus, the active constraint selected." << std::endl;

		return 0;
	}
}
	/*if ( verbose_level >= 70 ) {
		Size total_viol( 0 );
		ConstraintCOP active( active_constraint( ConformationXYZ( pose.conformation()), pose.energies().weights() ) );
		if ( verbose_level >=80 )	{
			// out << type() << " " << member_constraints().size() << " ";
			//for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++) {
				//if ( active_constraint_ ) {
				//	if ( (*cst_it).get() == active_constraint_.get() ) {
						// figure out if it's inter-res, residue_pair, or 3+body
						utility::vector1< int > pos_list;
						// generate list of all involved residues
						for ( Size i=1; i<= active->natoms(); ++i ) {
							int const seqpos( active->atom(i).rsd() );
							// seqpos already in list?
							if ( std::find( pos_list.begin(), pos_list.end(), seqpos )== pos_list.end() ) {
								pos_list.push_back( seqpos );
							}
						}
						if ( pos_list.size() == 2 ) {
							out << "ResiduePairConstraints ( " << pos_list[ 1 ] << " , " << pos_list[ 2 ] << " ) ";
							if ( verbose_level>80 ) out << std::endl;
							//out << "active ";
							total_viol += active->show_violations( out, pose, verbose_level, threshold );
							out << " of total: " << 1 << " ";
						}
				//	} else {
						//					out << "       ";
						//don't count the violations from this one
						//					(*cst_it)->show_violations( out, pose, verbose_level, threshold );
				//	}
				//}
			//}
			return total_viol;
		} else {
			return active->show_violations( out, pose, verbose_level );
			//else out << "WARNING: requested to show violations of an ambiguous constraint before"
			//				 <<" the pose was scored and thus, the active constraint selected." << std::endl;
			//
			//return 0;
		}
	}
	return 0;*/
//}

void
AmbiguousConstraint::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	func::FuncFactory const & func_factory
) {
	ConstraintOP constr;
	while( ( constr = ConstraintIO::read_individual_constraint_new( data, pose, func_factory ) ) != 0) {
		add_individual_constraint(constr);
	}
}

} //constraints
} //scoring
} //core

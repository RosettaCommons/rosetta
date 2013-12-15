// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief meta constraint where N constraints declared
/// @brief only the lowest K are evaluated
/// @author


#include <core/scoring/constraints/KofNConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <core/scoring/constraints/ConstraintIO.hh>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/conformation/Conformation.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/SequenceMapping.fwd.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintCreator.fwd.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/KofNConstraint.fwd.hh>
#include <core/scoring/constraints/MultiConstraint.fwd.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/random/random.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>


namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
KofNConstraint::KofNConstraint( core::Size K /*=0*/ ) : MultiConstraint() {
	K_ = K;
	init_cst_score_types();
	assert ( member_constraints().size() == 0 );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
KofNConstraint::KofNConstraint(  ConstraintCOPs & cst_in , core::Size K /*=0*/  ) : MultiConstraint( cst_in ) {
	K_ = K;
	init_cst_score_types();
	assert ( member_constraints().size() > 0 );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
KofNConstraint::init_cst_score_types()
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

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief ScoreFunction, scores all member constraints; reports the lowest k
void
KofNConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const {
	//std::cout << "scoring K of N constraint..." << std::endl;

	if (K_ == 0) {
		return; // ? warning msg here?
	}
	runtime_assert( K_ <= member_constraints().size() );  // bomb out if K<N

  cutoff_cst_score_ = 1000000;
  active_constraints_.clear();

  utility::vector1<core::Real> all_scores;
	utility::vector1<EnergyMap> tmp_EMaps;

	// step 1 score
	for( ConstraintCOPs::const_iterator member_it = member_constraints().begin(); member_it != member_constraints().end(); member_it++) {
    EnergyMap emap_i;
		(*member_it)->score(xyz_func, weights, emap_i);

		tmp_EMaps.push_back( emap_i );
    all_scores.push_back( calculate_total_cst_score( weights, emap_i ) );
	}

	// step 2 choose cutoff
	//   what if there is a tie?
	//   allow more than K_ to be active <-- this will mess up scores but derivatives will be ok
  utility::vector1<core::Real> sort_scores = all_scores;
	std::sort( sort_scores.begin(), sort_scores.end() );
	cutoff_cst_score_ = sort_scores[ K_ ];

	for (core::Size i=1; i<=all_scores.size(); ++i) {
		if( all_scores[i] <= cutoff_cst_score_ ){
			active_constraints_.push_back( member_constraints()[i] );

			emap[constant_constraint]      += tmp_EMaps[i][constant_constraint];
			emap[coordinate_constraint]    += tmp_EMaps[i][coordinate_constraint];
			emap[atom_pair_constraint]     += tmp_EMaps[i][atom_pair_constraint];
			emap[angle_constraint]         += tmp_EMaps[i][angle_constraint];
			emap[dihedral_constraint]      += tmp_EMaps[i][dihedral_constraint];
			emap[backbone_stub_constraint] += tmp_EMaps[i][backbone_stub_constraint];
			emap[backbone_stub_linear_constraint] += tmp_EMaps[i][backbone_stub_linear_constraint];
		}
	}
} //score


/// @brief helper function to accumulate all constraint scores into one number
core::Real
KofNConstraint::calculate_total_cst_score( EnergyMap const & weights, EnergyMap & emap ) const
{
	core::Real total_score =
			emap[constant_constraint] * weights[constant_constraint] +
			emap[coordinate_constraint] * weights[coordinate_constraint] +
			emap[atom_pair_constraint] * weights[atom_pair_constraint] +
			emap[angle_constraint] * weights[angle_constraint] +
			emap[dihedral_constraint] * weights[dihedral_constraint] +
			emap[backbone_stub_constraint] * weights[backbone_stub_constraint] +  // + was a semicolon! ~Labonte
			emap[backbone_stub_linear_constraint] * weights[backbone_stub_linear_constraint];

	return total_score;
}


ConstraintOP
KofNConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const {
	ConstraintCOPs new_csts;
	for( ConstraintCOPs::const_iterator cst_it = member_constraints_.begin(); cst_it != member_constraints_.end(); ++cst_it ){
		ConstraintOP new_cst = (*cst_it)->remap_resid( seqmap );
		if( new_cst ) new_csts.push_back( new_cst );
	}
	if( new_csts.size() > 0 ){
		return ConstraintOP( new KofNConstraint( new_csts ) );
	}
	else return NULL;
}


/// @brief function to minimize lowest scoring member constraint
void
KofNConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const {
	if (K_ == 0) {
		return; // warning msg?
	}
	runtime_assert( active_constraints_.size() != 0 );
	for (core::Size i=1; i<=active_constraints_.size(); ++i) {
		active_constraints_[i]->fill_f1_f2(atom, xyz, F1, F2, weights);
	}
}

void
KofNConstraint::show( std::ostream& out) const
{
	out << "KofNConstraint active constraints (K=" << K_ << "):" << std::endl;
	for (core::Size i=1; i<=active_constraints_.size(); ++i) {
		active_constraints_[i]->show(out);
	}
  out << "KofNConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
  for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
    (*cst_it)->show(out);
  }
  out << " ...all member constraints of this KofNConstraint shown." << std::endl;
}

utility::vector1<ConstraintCOP>
KofNConstraint::active_constraints() const {
	return active_constraints_;
}

Size
KofNConstraint::show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold ) const {
	if (K_ == 0) { return 0; }

	Size total_viol( 0 );
	if ( verbose_level >=80 )
		out << type() << " " << K_ << " of " << member_constraints().size() << " ";

	for (core::Size i=1; i<=active_constraints_.size(); ++i) {
		total_viol += active_constraints_[i]->show_violations( out, pose, verbose_level, threshold );
	}

	return total_viol;
}


void
KofNConstraint::read_def( std::istream& data, core::pose::Pose const& pose,func::FuncFactory const& func_factory ) {
	data >> K_;
	ConstraintOP constr;
	while( ( constr = ConstraintIO::read_individual_constraint_new( data, pose, func_factory ) ) != 0) {
		add_individual_constraint(constr);
	}
	std::cout << "Read K of N constraints! K = " << K_ << " N = " << member_constraints().size() << std::endl;
}

// void
// KofNConstraint::read_def(
// 	std::istream& data,
// 	core::pose::Pose const& pose,
// 	FuncFactory const& func_factory
// ) {
// 	Size res1, res2;
// 	std::string name1, name2;
// 	std::string func_type;
// 	std::string type;

// 	data
// 		>> name1 >> res1
// 		>> name2 >> res2
// 		>> func_type;

// 	tr.Debug << "read: " << name1 << " " << name2 << " " << res1 << " " << res2 << " func: " << func_type << std::endl;
// 	if ( res1 > pose.total_residue() || res2 > pose.total_residue() ) {
// 		tr.Warning 	<< "ignored constraint (no such atom in pose!)"
// 			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
// 		data.setstate( std::ios_base::failbit );
// 		return;
// 	}

// 	atom1_ = id::AtomID( id::NamedAtomID( name1, res1 ), pose );
// 	atom2_ = id::AtomID( id::NamedAtomID( name2, res2 ), pose );
// 	if ( atom1_.atomno() == 0 || atom2_.atomno() == 0 ) {
// 		tr.Warning << "Error reading atoms: read in atom names("
// 			<< name1 << "," << name2 << "), "
// 			<< "and found AtomIDs (" << atom1_ << "," << atom2_ << ")" << std::endl;
// 			data.setstate( std::ios_base::failbit );
// 			return;
// 	}

// 	func_ = func_factory.new_func( func_type );
// 	func_->read_data( data );

// 	if ( tr.Debug.visible() ) {
// 		func_->show_definition( std::cout );
// 		std::cout << std::endl;
// 	}
// } // parse_ambigous_constraint




} //constraints
} //scoring
} //core

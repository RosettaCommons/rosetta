// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ResidueCouplingConstraint.cc
///
/// @brief
/// @author Moritz Ertelt


#include <core/scoring/constraints/ResidueCouplingConstraint.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility>
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

namespace core {
namespace scoring {
namespace constraints {

static basic::Tracer TR( "core.scoring.constraints.ResidueCouplingConstraint" );

ResidueCouplingConstraint::ResidueCouplingConstraint():
	Constraint( core::scoring::res_type_linking_constraint )

{}

ResidueCouplingConstraint::ResidueCouplingConstraint(
	core::pose::Pose const &, //pose,
	Size const& seqpos1,
	Size const& seqpos2,
	Size const& tensor_index1,
	Size const& tensor_index2,
	CouplingTensorCOP const& tensor_const,
	Real const& strength,
	std::string const& alphabet
):
	Constraint( core::scoring::res_type_linking_constraint ),

	seqpos1_( seqpos1 ),
	seqpos2_( seqpos2 ),
	tensor_index1_ ( tensor_index1 ),
	tensor_index2_ ( tensor_index2 ),
	tensor_const_(tensor_const),
	strength_(strength),
	alphabet_(alphabet)
{
	// TR << "Using alphabet for TENSOR lookup: '" << get_alphabet() << "'" << std::endl;
	for ( size_t i = 0; i < get_alphabet().size(); ++i ) {
		char aa = get_alphabet()[i];
		amino_acids_[aa] = i;
	}
}

ResidueCouplingConstraint::~ResidueCouplingConstraint() = default;

ConstraintOP
ResidueCouplingConstraint::clone() const
{
	return ConstraintOP( utility::pointer::make_shared< ResidueCouplingConstraint >( *this ) );
}

Size ResidueCouplingConstraint::natoms() const {
	return 0;
}

id::AtomID const &
ResidueCouplingConstraint::atom( Size const ) const {
	utility_exit_with_message("ResidueCouplingConstraint is not atom-based!.");
	return core::id::GLOBAL_BOGUS_ATOM_ID;  // required for compilation on Windows
}


utility::vector1< core::Size >
ResidueCouplingConstraint::residues() const {
	utility::vector1< core::Size > pos_list;
	pos_list.push_back(seqpos1_);
	pos_list.push_back(seqpos2_);
	return pos_list;
}

std::string
ResidueCouplingConstraint::get_alphabet() const {
	return alphabet_;
}

void
ResidueCouplingConstraint::set_alphabet(std::string const& alphabet) {
	alphabet_ = alphabet;
}

void
ResidueCouplingConstraint::show( std::ostream & out ) const {
	out << "ResidueCouplingConstraint; ";
	out << "seqpos1: " << seqpos1_;
	out << "seqpos2: " << seqpos2_;
}

bool
ResidueCouplingConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( !           same_type_as_me( other_cst ) ) return false;
	if ( ! other_cst.same_type_as_me( *this ) ) return false;

	auto const & other( static_cast< ResidueCouplingConstraint const & > (other_cst) );

	if ( seqpos1_ != other.seqpos1_ ) return false;
	if ( seqpos2_ != other.seqpos2_ ) return false;

	return true;
}

bool ResidueCouplingConstraint::same_type_as_me( Constraint const & other ) const {
	return dynamic_cast< ResidueCouplingConstraint const * > (&other);
}

// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
ResidueCouplingConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{
	Real const weight(weights[ this->score_type() ]);

	if ( weight == 0 ) return;

	conformation::Residue const & rsd1( xyz_func.residue(seqpos1_) );
	conformation::Residue const & rsd2( xyz_func.residue(seqpos2_) );

	// looking up the coupling strength in the tensor (!Tensor index starts at ZERO!)
	utility::fixedsizearray1< numeric::Size, 4 > look_up {
		(tensor_index1_),
		amino_acids_.find(rsd1.type().name1())->second,
		(tensor_index2_),
		amino_acids_.find(rsd2.type().name1())->second};

	double coupling_strength((*tensor_const_)(look_up));

	// adding the coupling_strength as bonus to the emap
	emap[ this->score_type() ] -= coupling_strength  * strength_;
}

core::Real
ResidueCouplingConstraint::dist( core::scoring::func::XYZ_Func const & xyz_func ) const {
	return xyz_func.residue(seqpos1_).aa() != xyz_func.residue(seqpos2_).aa();
}

void
ResidueCouplingConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	func::XYZ_Func const &,
	Vector & ,//F1,
	Vector & ,//F2,
	EnergyMap const & //weights
) const
{
	// Do nothing.
	// Derivative of this restraint is effectively zero
	// so we just "add zero" to F1 and F2.
}

} // namespace constraints
} // namespace scoring
} // namespace core

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/constraints/FabConstraint.hh
/// @brief This class is specific to antibodies and penalizes presence of residues flanking
///        antibody cdr residues at Antigen-Antibody interfaces (ported from Fab constraint
///        in rosetta++ which uses a constant constraint score of 0.5/flanking residue)
/// @author Krishna Kilambi (kkpraneeth@jhu.edu, April 2012)

#ifndef INCLUDED_core_scoring_constraints_FabConstraint_hh
#define INCLUDED_core_scoring_constraints_FabConstraint_hh

// Unit header

#include <core/scoring/constraints/FabConstraint.fwd.hh>
#include <core/scoring/constraints/MultiConstraint.hh>




#include <core/pose/Pose.fwd.hh>


//Utility Headers

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


class FabConstraint : public MultiConstraint {
public:

	/// @brief Constructor
	FabConstraint();

	/// @brief Constructor
	FabConstraint(ConstraintCOPs const & cst_in) ;

	///
	ConstraintOP clone() const override;

	std::string type() const override;

	void
	show(std::ostream& out) const override;

	/// @brief read in constraint definition
	void
	read_def(std::istream& data, pose::Pose const& pose,func::FuncFactory const& func_factory) override;

	Size
	pose_res_no(core::pose::Pose const & pose, std::string const & res_designation);

	utility::vector1<Real>
	calc_penalty_vec(Size start_res, Size stop_res, utility::vector1<Size> res1, utility::vector1<Size> res2);

	void
	setup_csts(core::pose::Pose const & pose, utility::vector1<Size> const & res1, utility::vector1<Size> const & res2, utility::vector1<std::string> const & antchains);

	bool operator==( Constraint const & rhs ) const override;
	bool same_type_as_me( Constraint const & other ) const override;

private:

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //FabConstraint

} //constraints
} //scoring
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_FabConstraint )
#endif // SERIALIZATION


#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
///
/// @file   protocols/cyclic_peptide/CreateTorsionConstraint.hh
/// @brief  Add torsion constraints to the current pose conformation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_cyclic_peptide_CreateTorsionConstraint_hh
#define INCLUDED_protocols_cyclic_peptide_CreateTorsionConstraint_hh

#include <protocols/moves/Mover.hh>

namespace protocols {
namespace cyclic_peptide {

class CreateTorsionConstraint : public moves::Mover
{
public:
	CreateTorsionConstraint();
	virtual ~CreateTorsionConstraint();

	virtual void apply( core::pose::Pose &pose );
	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, core::pose::Pose const & );

private:
	utility::vector1<Size> res1_;
	utility::vector1<std::string> atom1_;
	utility::vector1<Size> res2_;
	utility::vector1<std::string> atom2_;
	utility::vector1<Size> res3_;
	utility::vector1<std::string> atom3_;
	utility::vector1<Size> res4_;
	utility::vector1<std::string> atom4_;
	utility::vector1<std::string> cst_func_;
};

} // moves
} // protocols

#endif //INCLUDED_protocols_cyclic_peptide_CreateTorsionConstraint_hh

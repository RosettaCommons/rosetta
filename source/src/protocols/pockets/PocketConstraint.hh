// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/constraints/PocketConstraint.hh
///
/// @brief
/// @author David Johnson


#ifndef INCLUDED_protocols_pockets_PocketConstraint_hh
#define INCLUDED_protocols_pockets_PocketConstraint_hh

#include <core/scoring/constraints/Constraint.hh>
#include <protocols/pockets/PocketConstraint.fwd.hh>
// AUTO-REMOVED #include <protocols/pockets/PocketGrid.hh>

// AUTO-REMOVED #include <math.h>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <list>

#include <core/conformation/Residue.fwd.hh>
#include <protocols/pockets/PocketGrid.fwd.hh>


namespace protocols {
namespace pockets {

///@brief This constraint favors creating a pocket suitable for a small-molecule
///
class PocketConstraint : public core::scoring::constraints::Constraint
{
public:
	virtual std::string type() const {
		return "Pocket";
	}

	PocketConstraint( );
	PocketConstraint( core::pose::Pose const & pose );
	PocketConstraint( const PocketConstraint& old );

	void init(core::pose::Pose const & pose);

	core::Size target_res() const { return seqpos_; }
	virtual ~PocketConstraint();

	virtual core::Size natoms() const { return atom_ids_.size(); };

	virtual core::id::AtomID const & atom( core::Size const index ) const { return atom_ids_[index]; };

	void show_def( std::ostream& out, core::pose::Pose const & pose ) const;
	void read_def( std::istream& in, core::pose::Pose const & pose, core::scoring::func::FuncFactory const & func_factory );

	virtual
	void score( core::scoring::func::XYZ_Func const & xyz_func, core::scoring::EnergyMap const & weights, core::scoring::EnergyMap & emap ) const;

	virtual
	void
	fill_f1_f2(
		core::id::AtomID const & ,
		core::scoring::func::XYZ_Func const & ,
		core::Vector & ,
		core::Vector & ,
		core::scoring::EnergyMap const & weights
	) const;

	virtual
	core::scoring::constraints::ConstraintOP clone() const;

	void set_target_res( core::pose::Pose const & pose, core::Size new_seqpos );
	void set_target_res_pdb(core::pose::Pose const & pose, std::string resid );

private:

	core::Size seqpos_;
	core::Size totalres_;
	core::Size angles_;
	core::Real weight_;
	mutable protocols::pockets::PocketGridOP pocketgrid_;
	utility::vector1< AtomID > atom_ids_;
	bool dumppdb_;
	std::vector< core::conformation::ResidueOP > residues_;

}; // PocketConstraint


} // namespace constraints_additional
} // namespace protocols


#endif // INCLUDED_protocols_pockets_PocketConstraint_HH

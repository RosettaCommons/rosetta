// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pockets/constraints/PocketConstraint.hh
///
/// @brief
/// @author David Johnson


#ifndef INCLUDED_protocols_pockets_PocketConstraint_hh
#define INCLUDED_protocols_pockets_PocketConstraint_hh

#include <core/scoring/constraints/Constraint.hh>
#include <protocols/pockets/PocketConstraint.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/conformation/Residue.fwd.hh>
#include <protocols/pockets/PocketGrid.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace pockets {

/// @brief This constraint favors creating a pocket suitable for a small-molecule
///
class PocketConstraint : public core::scoring::constraints::Constraint
{
public:
	std::string type() const override {
		return "Pocket";
	}

	PocketConstraint( );
	PocketConstraint( core::pose::Pose const & pose );
	PocketConstraint( const PocketConstraint& old );

	void init(core::pose::Pose const & pose);

	core::Size target_res() const { return seqpos_; }
	~PocketConstraint() override;

	core::Size natoms() const override { return atom_ids_.size(); };

	core::id::AtomID const & atom( core::Size const index ) const override { return atom_ids_[index]; };

	void show_def( std::ostream& out, core::pose::Pose const & pose ) const override;
	void read_def( std::istream& in, core::pose::Pose const & pose, core::scoring::func::FuncFactory const & func_factory ) override;

	
	void score( core::scoring::func::XYZ_Func const & xyz_func, core::scoring::EnergyMap const & weights, core::scoring::EnergyMap & emap ) const override;

	
	void
	fill_f1_f2(
		core::id::AtomID const & ,
		core::scoring::func::XYZ_Func const & ,
		core::Vector & ,
		core::Vector & ,
		core::scoring::EnergyMap const & weights
	) const override;

	
	core::scoring::constraints::ConstraintOP clone() const override;

	bool operator == ( core::scoring::constraints::Constraint const & other ) const override;
	bool same_type_as_me ( core::scoring::constraints::Constraint const & other ) const override;

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
	std::vector< core::conformation::ResidueCOP > residues_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // PocketConstraint


} // namespace constraints_additional
} // namespace protocols


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_pockets_PocketConstraint )
#endif // SERIALIZATION


#endif // INCLUDED_protocols_pockets_PocketConstraint_HH

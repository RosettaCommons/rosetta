// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/idealize/IdealizeMover.hh
/// @brief protocols for idealizing a Pose
/// @author


#ifndef INCLUDED_protocols_idealize_IdealizeMover_hh
#define INCLUDED_protocols_idealize_IdealizeMover_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/idealize/IdealizeMover.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
//Auto Headers
#include <utility/vector1.hh>
namespace protocols {
namespace idealize {


/// @brief Mover class for transforming a Pose to ideal bonds of the given Pose.
/// The idea is that this Mover stochastically picks a move-able position, forces
/// that position into ideal geometry, and tries to use minimization to bring
/// the coordinates back to very near their starting points.

class IdealizeMover : public moves::Mover {

public:
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

public:
	IdealizeMover():
		Mover("IdealizeMover"),
		atom_pair_constraint_weight_( 0.0 ),
		coordinate_constraint_weight_( 0.01 ),
		cis_omega_( false ),
		fast_( false ),
		chainbreaks_( false ),
		report_CA_rmsd_(true),
		impose_constraints_( true ),
		constraints_only_( false )
	{
		ignore_residues_in_csts_.clear();
	}

	/// @brief clone has to be overridden only if clone invocation is expected.
	virtual moves::MoverOP clone() const {
		return moves::MoverOP( new IdealizeMover( *this ) );
	}

	virtual moves::MoverOP fresh_instance() const {
		return moves::MoverOP( new IdealizeMover );
	}

	IdealizeMover &
	atom_pair_constraint_weight( core::Real const setting )
	{
		atom_pair_constraint_weight_ = setting;
		return *this;
	}


	IdealizeMover &
	coordinate_constraint_weight( core::Real const setting )
	{
		coordinate_constraint_weight_ = setting;
		return *this;
	}

	IdealizeMover &
	fast( bool const setting )
	{
		fast_ = setting;
		return *this;
	}

	IdealizeMover &
	chainbreaks( bool const setting ){
		chainbreaks_ = setting;
		return *this;
	}

	IdealizeMover &
	cis_omega( bool const setting ){
		cis_omega_ = setting;
		return *this;
	}

	IdealizeMover &
	report_CA_rmsd( bool const setting )
	{
		report_CA_rmsd_ = setting;
		return *this;
	}

	void
	apply( Pose & pose );

	virtual std::string get_name() const;

	/// @brief Sets the list of residue positions to idealize.
	void set_pos_list( utility::vector1< Size > pos_list ) {
		pos_list_ = pos_list;
	}

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	utility::vector1< core::Size > ignore_residues_in_csts() const;
	void ignore_residues_in_csts( utility::vector1< core::Size > const & i );
	void impose_constraints( bool const i ){ impose_constraints_ = i; }
	bool impose_constraints() const{ return( impose_constraints_ ); }
	bool constraints_only() const{ return constraints_only_; }
	void constraints_only( bool const c ){ constraints_only_ = c; }

private:
	// methods
	void
	setup_idealize_constraints( Pose & pose );

private:
	// data
	utility::vector1< Size > pos_list_;

	Real atom_pair_constraint_weight_;
	Real coordinate_constraint_weight_;
	utility::vector1< core::Size > ignore_residues_in_csts_; // dflt empty; residues that don't carry constraints, e.g., residues that were inserted without closing loops etc.
	bool cis_omega_;
	bool fast_;
	bool chainbreaks_;
	bool report_CA_rmsd_;
	bool impose_constraints_; //dflt true; should the mover impose its own constraints or respect the incoming ones;
	bool constraints_only_; //dflt false; jump out after imposing the idealize constraints?
};

} // idealize
} // protocols

#endif

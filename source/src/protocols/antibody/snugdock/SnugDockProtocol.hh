// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/SnugDockProtocol.hh
/// @brief Dock and antigen to an antibody while optimizing the rigid body orientation of the VH and VL chains and performing CDR loop minimization.
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )


#ifndef INCLUDED_protocols_antibody_SnugDockProtocol_HH
#define INCLUDED_protocols_antibody_SnugDockProtocol_HH

// Unit headers
#include <protocols/antibody/snugdock/SnugDockProtocol.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/RefineOneCDRLoop.fwd.hh>

// Project headers
#include <protocols/docking/DockingProtocol.fwd.hh>

// C++ headers
#include <iostream>

using namespace core;
namespace protocols {
namespace antibody {

class SnugDockProtocol: public moves::Mover {
public: // boiler plate / virtuals
	// default constructor
	SnugDockProtocol();

	// copy constructor
	SnugDockProtocol( SnugDockProtocol const & rhs );

	// assignment operator
	SnugDockProtocol & operator=( SnugDockProtocol const & rhs );

	// destructor
	virtual ~SnugDockProtocol();

	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	virtual bool reinitialize_for_new_input() const;

	/// @brief Associates relevant options with the SnugDockProtocol class
	static void register_options();

	// Accessors for auto kink-constraint options
	bool auto_generate_kink_constraint() const { return auto_generate_kink_constraint_; }
	void auto_generate_kink_constraint( bool const setting ) { auto_generate_kink_constraint_ = setting; }

	bool high_res_kink_constraint() const { return high_res_kink_constraint_; }
	void high_res_kink_constraint( bool const setting ) { high_res_kink_constraint_ = setting; }


public:
	void show( std::ostream & out=std::cout ) const;
	friend std::ostream & operator<<(std::ostream& out, SnugDockProtocol const & snugdockprotocol );

private: // methods
	void setup_objects( Pose const & pose );
	void setup_loop_refinement_movers();
	void init();
	void init_for_equal_operator_and_copy_constructor( SnugDockProtocol & lhs, SnugDockProtocol const & rhs);
	void init_from_options();
	void set_default();

	docking::DockingProtocolOP docking() const;

private: // data
	AntibodyInfoOP antibody_info_;

	// Movers
	RefineOneCDRLoopOP low_res_refine_cdr_h2_;
	RefineOneCDRLoopOP low_res_refine_cdr_h3_;
	mutable docking::DockingProtocolOP docking_;

	std::string loop_refinement_method_;

	// H3 filter options
	bool h3_filter_;
	Size h3_filter_tolerance_;

	// auto kink-constraint options
	bool auto_generate_kink_constraint_;
	bool high_res_kink_constraint_;


}; // class SnugDockProtocol

} // namespace antibody
} // namespace protocols

#endif // INCLUDED_protocols_antibody_SnugDockProtocol_HH

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Assigns a ResidueTypeConstraint to a pose.
/// @author Doo Nam Kim (doonam.kim@gmail.com)

#ifndef INCLUDED_protocols_constraint_movers_ResidueTypeConstraintMover_hh
#define INCLUDED_protocols_constraint_movers_ResidueTypeConstraintMover_hh

#include <protocols/constraint_movers/ResidueTypeConstraintMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/constraints/ResidueTypeConstraint.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace constraint_movers {

// TODO:
// 1. hold a std::set of AA_name3s, and error if the values provided in the string setter
// or in parse_my_tag have duplicates, because that's probably a user mistake with
// semantic consequences.
// 2. make the string setter safe to use if it's a CSV itself.
// 3. ensure that the name3s provided are real name3s

class ResidueTypeConstraintMover : public protocols::moves::Mover {

public:
	typedef core::scoring::constraints::ResidueTypeConstraint ResidueTypeConstraint;
	//  typedef core::scoring::constraints::ResidueTypeConstraintOP ResidueTypeConstraintOP;
	//  typedef core::scoring::constraints::ResidueTypeConstraintCOP ResidueTypeConstraintCOP;

public:
	ResidueTypeConstraintMover();
	~ResidueTypeConstraintMover() override;
	ResidueTypeConstraintMover( std::string const & );

	//  void constraint_file( std::string const & );
	//
	//  void constraint_set( ResidueTypeConstraintCOP );
	//  ResidueTypeConstraintOP constraint_set();
	//  ResidueTypeConstraintCOP constraint_set() const;

	void apply( Pose & ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	void set_AA_name3( utility::vector1< std::string > const & setting ) { AA_name3s_ = setting; }
	void set_AA_name3( std::string const & setting ) { AA_name3s_ = utility::vector1< std::string >( 1, setting ); }
	void set_favor_bonus( core::Real const setting ) { favor_bonus_ = setting; }

private:
	utility::vector1< std::string > AA_name3s_;
	core::Real favor_bonus_;
};

} // moves
} // protocols

#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/MgMonteCarlo.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_MgMonteCarlo_HH
#define INCLUDED_protocols_magnesium_MgMonteCarlo_HH

#include <protocols/moves/Mover.hh>
#include <protocols/magnesium/MgMonteCarlo.fwd.hh>

namespace protocols {
namespace magnesium {

class MgMonteCarlo: public moves::Mover {

public:

	//constructor
	MgMonteCarlo();

	//destructor
	~MgMonteCarlo();

	void set_temperature( core::Real const & setting ){ temperature_ = setting; }
	core::Real temperature() const { return temperature_; }

	void set_add_delete_frequency( core::Real const & setting ){ add_delete_frequency_ = setting; }
	core::Real add_delete_frequency() const { return add_delete_frequency_; }

	void set_cycles( core::Size const & setting ){ cycles_ = setting; }
	core::Size cycles() const { return cycles_; }

	void set_output_pdb( bool const & setting ){ output_pdb_ = setting; }
	bool output_pdb() const { return output_pdb_; }

public:

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const{ return "MgMonteCarlo"; }

private:

	void
	setup_mg_water_fold_tree( core::pose::Pose & pose ) const;

private:

	core::Size cycles_;
	core::Real temperature_;
	core::Real add_delete_frequency_;
	bool output_pdb_;

};

} //magnesium
} //protocols

#endif

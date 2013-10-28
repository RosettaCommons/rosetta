// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/SpinMover.hh
/// author Eva-Maria Strauch ( evas01@u.washington.edu )

#ifndef INCLUDED_protocols_protein_interface_design_movers_SpinMover_HH
#define INCLUDED_protocols_protein_interface_design_movers_SpinMover_HH

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class SpinMover : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
//	typedef core::conformation::Residue Residue;

public:
	SpinMover();
	SpinMover( core::Size jump_num );
	/// @param jump_num The jump number of the interface. 0 for no interface
	/// @note Pass everything through the final filter (True Filter)

	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	//virtual std::string protocols::moves::Mover::get_name();
	void parse_my_tag( utility::tag::TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	protocols::moves::MoverOP clone() const { return( protocols::moves::MoverOP( new SpinMover( *this ) ) ); }
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new SpinMover ); }
	virtual ~SpinMover() {};
private:
	core::Size jump_num_;
};


} //movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_SpinMover_HH*/


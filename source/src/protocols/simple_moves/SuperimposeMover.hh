// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SuperimposeMover.hh
/// @brief
/// @author Ingemar Andre

#ifndef INCLUDED_protocols_simple_moves_SuperimposeMover_hh
#define INCLUDED_protocols_simple_moves_SuperimposeMover_hh

// Unit headers
#include <protocols/simple_moves/SuperimposeMover.fwd.hh>
#include <protocols/moves/Mover.hh> // we need to store a pose
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

///@brief Superimpose pose onto reference pose. 
///  Default CA only.  All residues.  Optionally set to superimpose regions or backbone only (N, C, CA, O)
class SuperimposeMover : public moves::Mover {

public:
	/// @brief
	/// 	empty constructor
	SuperimposeMover();

	SuperimposeMover( core::pose::Pose const & ref_pose );

	SuperimposeMover( core::pose::Pose const & ref_pose , core::Size ref_start, core::Size ref_end, core::Size target_start, core::Size target_end, bool CA_only);

	~SuperimposeMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void set_reference_pose( core::pose::Pose const & pose, Size start=1, Size end=0);
	// Undefined, commenting out to fix PyRosetta build  void set_target_pose( Size start=1, Size end=-1);
	void set_target_range( Size start, Size end );
	
	void set_ca_only(bool setting);
	
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;

private:
	core::Real superimpose( core::pose::Pose & mod_pose, core::pose::Pose const & ref_pose, Size ref_start, Size ref_end, Size target_start, Size target_end);
	core::Real superimposebb( core::pose::Pose & mod_pose, core::pose::Pose const & ref_pose, Size ref_start, Size ref_end, Size target_start, Size target_end);

private:

	core::pose::PoseOP ref_pose_;
	core::Size ref_start_, ref_end_;
	core::Size target_start_, target_end_;
        bool CA_only_;

};

} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_SuperimposeMover_HH

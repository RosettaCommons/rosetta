// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/StartFrom.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_StartFrom_hh
#define INCLUDED_protocols_ligand_docking_StartFrom_hh

#include <protocols/ligand_docking/StartFrom.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class StartFrom: public protocols::moves::Mover
{
public:
	StartFrom();
	virtual ~StartFrom();
	StartFrom(StartFrom const & that);

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual std::string get_name() const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	void coords(core::Vector const & coords,std::string const & pdb_tag);

	void chain(std::string const & chain);

	void apply(core::pose::Pose & pose);

	// Undefined, commenting out to make PyRosetta compile
	//StartFrom(core::pose::Pose & pose);

	void parse_startfrom_file(std::string filename);

private:
	std::string chain_;
	std::map< std::string, utility::vector1<core::Vector> > starting_points_;
	std::map<std::string,core::Vector > potential_starting_positions_;
	bool use_file_name_;
};

void move_ligand_to_desired_centroid(
	core::Size const jump_id,
	core::Vector const desired_centroid,
	core::pose::Pose & pose
);

} //namespace ligand_docking
} //namespace protocols

#endif

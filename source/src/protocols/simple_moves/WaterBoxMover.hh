// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file WaterBoxMover.hh

#ifndef INCLUDED_protocols_moves_WaterBoxMover_hh
#define INCLUDED_protocols_moves_WaterBoxMover_hh

#include <protocols/simple_moves/WaterBoxMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// C++ Headers
#include <string>

namespace protocols {
namespace simple_moves {

class WaterBoxMover : public protocols::moves::Mover {
public:
	WaterBoxMover();

	static std::string mover_name();

	virtual void apply( Pose & pose );

	void add_water( Pose & pose, core::Vector O, core::Size resnum );

	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & pose );

	void del_virt( bool dv ) { del_virt_ = dv; }
	bool del_virt() const { return del_virt_; }

	void point_water( bool pw ) { point_waters_ = pw; }
	bool point_water() const { return point_waters_; }

	void
	delete_virtual_waters( core::pose::Pose  & pose );

	void
	bb_sol( Pose & pose, core::pack::task::PackerTaskCOP task );

	void
	bb_sol( core::pose::Pose & pose );

	void
	sc_sol( Pose & pose, core::pack::task::PackerTaskCOP task );

	void
	sc_sol( core::pose::Pose & pose );

	void
	lig_sol( core::pose::Pose & pose );

	static
	utility::tag::XMLSchemaComplexTypeGeneratorOP
	define_water_box_mover_schema();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	core::pack::task::PackerTaskCOP task_;
	core::pack::task::TaskFactoryCOP task_factory_;

	bool point_waters_;
	bool make_rotatable_, make_point_, /*make_point_1_,*/ del_virt_, remove_all_;
	bool bb_sol_, sc_sol_, lig_sol_;
};

} // moves
} // protocols

#endif

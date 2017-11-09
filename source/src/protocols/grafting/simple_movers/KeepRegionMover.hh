// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/grafting/simple_movers/KeepRegionMover.hh
/// @brief Keep a region of a pose, delete the rest
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_grafting_simple_movers_KeepRegionMover_hh
#define INCLUDED_protocols_grafting_simple_movers_KeepRegionMover_hh

#include <protocols/grafting/simple_movers/KeepRegionMover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace grafting {
namespace simple_movers {

/// @brief Keep a contiguous region of a pose, delete the rest.
/// Re-detect disulfides
class KeepRegionMover : public protocols::moves::Mover {
public:

	KeepRegionMover();
	KeepRegionMover(std::string const & res_start, std::string const & res_end);

	KeepRegionMover(KeepRegionMover const & src);

	virtual ~KeepRegionMover();

	void
	apply(core::pose::Pose & pose) override;


public:

	/// @brief Set the region of the pose where we keep the residues
	void
	region(std::string const & res_start, std::string const & res_end);

	std::pair<std::string, std::string>
	region() const;

	/// @brief Set the first residue that we will keep
	void
	start(std::string const & res_start);

	void
	start( core::Size res_start );

	std::string const &
	start() const;

	/// @brief Set the last residue that we will keep
	void
	end(std::string const & res_end);

	void
	end( core::Size res_end );

	std::string const &
	end() const;



public:


	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		Pose const & pose) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string start_;
	std::string end_;
	core::Size nter_overhang_;
	core::Size cter_overhang_;
};


}
}
}


#endif



// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/PatchdockTransform.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PatchdockTransform_hh
#define INCLUDED_protocols_protein_interface_design_movers_PatchdockTransform_hh
#include <protocols/protein_interface_design/movers/PatchdockTransform.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/protein_interface_design/PatchdockReader.fwd.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief wrapper around protocols::protein_interface_design::PatchdockReader class. That class is derived from JobInputter and handles input situations that involve patchdock output files. Here, we provide an entry point for a patchdock transformation within the RosettaScripts trajectory.
class PatchdockTransform : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
	typedef protocols::protein_interface_design::PatchdockReaderOP PatchdockReaderOP;
public:
	PatchdockTransform();
	void apply( Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new PatchdockTransform ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual ~PatchdockTransform();
	PatchdockReaderOP pd_reader() const;
	void pd_reader( PatchdockReaderOP p );
private:
	PatchdockReaderOP pd_reader_;
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_PatchdockTransform_HH*/

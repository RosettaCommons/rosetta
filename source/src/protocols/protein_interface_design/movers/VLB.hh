// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/VLB.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_VLB_hh
#define INCLUDED_protocols_protein_interface_design_movers_VLB_hh

//#include <protocols/protein_interface_design/movers/DesignRepackMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/forge/build/BuildManager.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief user interface for YAB's Variable Length Build.
class VLB : public protocols::moves::Mover
{
public:
	VLB(); // default ctor//design_ = true;
	VLB( VLB const & init ); // copy ctor
	~VLB(); // dtor
	VLB & operator=( VLB const & init );

	VLB( protocols::forge::build::BuildManagerCOP manager, core::scoring::ScoreFunctionCOP scorefxn );
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose );

private:
	protocols::forge::build::BuildManagerOP manager_;
	core::scoring::ScoreFunctionOP scorefxn_;
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_VLB_HH*/

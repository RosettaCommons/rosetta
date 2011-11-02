// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/DomainAssembly.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_DomainAssembly_hh
#define INCLUDED_protocols_protein_interface_design_movers_DomainAssembly_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief insert fragments in a linker region. Very similar to what's in looprelax_main
class DomainAssembly : public protocols::moves::Mover
{
typedef core::fragment::FragSetOP FragSetOP;
typedef core::fragment::FragSetCOP FragSetCOP;
public:
	DomainAssembly();
	DomainAssembly(
		core::Size const linker_start,
		core::Size const linker_end,
		FragSetOP fragset_large,
		FragSetOP fragset_small
	);
	virtual ~DomainAssembly();
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);
private:
	core::Size linker_start_, linker_end_;
	core::fragment::FragSetOP fragset_large_, fragset_small_;
	bool fragments_set_;
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_DomainAssembly_HH*/

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/CompoundTranslate.hh
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_CompoundTranslate_hh
#define INCLUDED_protocols_ligand_docking_CompoundTranslate_hh

// Unit Headers
#include <protocols/ligand_docking/Translate.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

class CompoundTranslate : public protocols::moves::Mover
{
public:
	CompoundTranslate();
	virtual ~CompoundTranslate();
	CompoundTranslate(CompoundTranslate const & that);

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual std::string get_name() const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	void apply(core::pose::Pose & pose);

private:
	TranslateOPs translates_;
	bool randomize_order_;
	bool allow_overlap_;
};

} //namespace ligand_docking
} //namespace protocols

#endif

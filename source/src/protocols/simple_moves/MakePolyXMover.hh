// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./protocols/simple_moves/MakePolyXMover.hh
/// @brief  header file of MakePolyXMover.cc
/// @author Nobuyasu Koga ( nobuyasau@uw.edu )

#ifndef INCLUDED_protocols_simple_moves_MakePolyXMover_HH
#define INCLUDED_protocols_simple_moves_MakePolyXMover_HH

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MakePolyXMover.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {


class MakePolyXMover : public protocols::moves::Mover {
public:

	typedef protocols::moves::MoverOP MoverOP;
	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;

public:


	// @brief default constructor
	MakePolyXMover();

	// @brief value constructor
	MakePolyXMover( std::string aa, bool keep_pro, bool keep_gly, bool keep_disulfide_cys );

	// @brief destructor
	~MakePolyXMover() override;

	/// @brief clone this object
	protocols::moves::MoverOP clone() const override;

	/// @brief create this type of object
	protocols::moves::MoverOP fresh_instance() const override;

	// @brief virtual main operation
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;


private:

	/// @brief using amino acid for converting pose to poly XXX
	std::string aa_;

	/// @brief if true, proline, proline are not converted
	bool keep_pro_;

	/// @brief if true, proline, glycine are not converted
	bool keep_gly_;

	/// @brief if true, proline, cystein are not converted
	bool keep_disulfide_cys_;


};


} // simple_moves
} // protocols


#endif //INCLUDED_protocols_simple_moves_MakePolyXMover_HH

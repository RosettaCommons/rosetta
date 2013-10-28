// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file protocols/protein_interface_design/movers/DisulfideMover.hh
/// @brief
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date 4/30/2009

#ifndef INCLUDED_protocols_protein_interface_design_movers_DisulfideMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_DisulfideMover_hh

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers

#include <core/scoring/disulfides/CentroidDisulfidePotential.fwd.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>



// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace movers {


class DisulfideMover : public protocols::simple_moves::DesignRepackMover
{
private:
	typedef protocols::simple_moves::DesignRepackMover parent;
public:
	///@brief default ctor
	DisulfideMover();
	///@brief copy ctor
	DisulfideMover(DisulfideMover const& dm);
	///@brief Constructor with a single target residue
	DisulfideMover( core::Size targetResidue );
	///@brief Constructor with multiple target residues
	DisulfideMover( utility::vector1<core::Size> const& targetResidues );
	virtual ~DisulfideMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const {
		return (protocols::moves::MoverOP( new DisulfideMover( *this ) ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new DisulfideMover );
	}

	void parse_my_tag( utility::tag::TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
public:
	/// @brief Find all residues which could disulfide bond to a target
	/// @return pairs of residues (target, host) from the target protein and the
	///   docking protein.
	/// @note This is implemented as a static method so that DisulfidedFilter
	///  can share code.
	static void disulfide_list( core::pose::Pose const & pose,
		utility::vector1< core::Size > const& targets, Size rb_jump,
		utility::vector1< std::pair<core::Size,core::Size> > & disulfides);
private:
	/// @brief Modify the pose to define a disulfide bond between the two specified
	///   residues.
	/// @details Does not do the repacking & minimization required to place the
	///   disulfide correctly.
	static void form_disulfide(core::pose::Pose & pose, core::Size lower_res, core::Size upper_res);

private:
	static const core::scoring::disulfides::CentroidDisulfidePotential potential_;
private:
	Size rb_jump_;
};

} // movers
} // protein_interface_design
} // protocols

#endif //INCLUDED_protocols_protein_interface_design_movers_DisulfideMover_HH_


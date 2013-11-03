// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/AddHydrogens.hh
///
/// @brief
/// @author Gordon Lemmon


#ifndef INCLUDED_protocols_ligand_docking_AddHydrogens_hh
#define INCLUDED_protocols_ligand_docking_AddHydrogens_hh

#include <protocols/ligand_docking/AddHydrogens.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>


#include <protocols/moves/Mover.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

///@brief
///
///@details
///
class AddHydrogens : public protocols::moves::Mover{

public:
	AddHydrogens();
	virtual ~AddHydrogens();
	AddHydrogens(AddHydrogens const & that);
	virtual void apply( core::pose::Pose & pose );

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


private:
	std::string chain_;

}; // class AddHydrogens

} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_AddHydrogens_HH

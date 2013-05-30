// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/MapHotspot.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_MapHotspot_hh
#define INCLUDED_protocols_protein_interface_design_movers_MapHotspot_hh

#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <map>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/DataMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class MapHotspot : public protocols::moves::Mover
{
public:
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::pose::Pose Pose;
	typedef core::conformation::Residue Residue;
	typedef core::pack::rotamer_set::RotamerSetOP RotamerSetOP;
	typedef std::map< core::Size, protocols::filters::FilterCOP > SizeFilter_map;

public:
	MapHotspot();
	// for direct access
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	/// @brief this is the recursive function where the functionality takes place
	void GenerateMap( core::pose::Pose const & start_pose, core::pose::Pose & curr_pose, core::Size const jump_number );
	/// @brief minimizes rb and sc dofs for all of the hotspots
	void MinimizeHotspots( core::pose::Pose & pose );
	void parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	protocols::moves::MoverOP clone() const { return( protocols::moves::MoverOP( new MapHotspot( *this ) ) ); }
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new MapHotspot ); }
	void output_pose( core::pose::Pose const & pose ) const;
	RotamerSetOP create_rotamer_set( core::pose::Pose const &, core::Size const hotspot_resnum, core::Size const explosion ) const;
	virtual ~MapHotspot();
private:
	bool clash_check_;
	std::map< core::Size, core::Size > explosion_; // rotamer explosion
	SizeFilter_map jump_filters_; // filter. Defaults to TrueFilter
	std::map< core::Size, std::string > allowed_aas_per_jump_; //defaults to "ADEIKLMNQRSTVWY"
	std::map< core::Size, core::scoring::ScoreFunctionCOP > minimization_scorefxns_;
	std::map< core::Size, protocols::moves::MoverOP > jump_movers_; //defaults to null mover
	std::string file_name_prefix_;
};


} //movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_TryRotamers_HH*/


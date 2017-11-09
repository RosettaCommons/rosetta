// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/TryRotamers.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_TryRotamers_hh
#define INCLUDED_protocols_protein_interface_design_movers_TryRotamers_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.hh> /// EVIL

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class TryRotamers : public protocols::moves::Mover
{
public:
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::pose::Pose Pose;
	typedef core::conformation::Residue Residue;
	typedef core::pack::rotamer_set::RotamerSetOP RotamerSetOP;

public:
	TryRotamers();
	TryRotamers( std::string const & resnum,
		core::scoring::ScoreFunction const& scorefxn,
		protocols::filters::Filter const& final_filter,
		core::Size explosion = 0, // rotamer explosion
		core::Size jump_num = 1,
		bool clash_check = false,
		bool solo_res = false,
		bool include_current = true
	);

	/// @param jump_num The jump number of the interface. 0 for no interface
	/// @note Pass everything through the final filter (True Filter)
	TryRotamers( std::string const & resnum,
		core::scoring::ScoreFunction const& scorefxn,
		core::Size explosion = 0, // rotamer explosion
		core::Size jump_num = 1,
		bool clash_check = false,
		bool solo_res = false,
		bool include_current = true
	);

	// for direct access
	void set_resnum( std::string const & r ) { resnum_ = r; }

	/// EVIL
	void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ){ scorefxn_ = scorefxn; }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;
	core::pack::rotamer_set::Rotamers::const_iterator begin() const { return rotset_->begin(); }
	core::pack::rotamer_set::Rotamers::const_iterator end() const { return rotset_->end(); }

	protocols::moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new TryRotamers( *this ) ) ); }
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new TryRotamers ); }
	virtual ~TryRotamers();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

protected:

	void setup_rotamer_set(core::pose::Pose & pose, core::Size resnum );

	/// EVIL
	core::pack::rotamer_set::RotamerSetOP rotamer_set() { return rotset_ ; }

	/// EVIL
	void get_rotamer_set( core::pack::rotamer_set::RotamerSetOP const & rs ){
		rotset_ = rs;
		rotamer_it_ = begin();
	}

private:
	core::scoring::ScoreFunctionCOP scorefxn_;
	core::pack::rotamer_set::Rotamers::const_iterator rotamer_it_;
	core::pack::rotamer_set::RotamerSetOP rotset_;
	std::string resnum_;
	core::Size jump_num_;
	bool clash_check_;
	bool solo_res_;
	bool include_current_;
	bool automatic_connection_; // should TryRotamers decide on the foldtree on its own? default true
	core::Size explosion_; // rotamer explosion
	protocols::filters::FilterOP final_filter_; // filter. Defaults to TrueFilter
	core::select::residue_selector::ResidueSelectorCOP shove_residues_; // residues for which to use the shove_bb atom type, so that backbone atoms might clash.
};


} //movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_TryRotamers_HH*/

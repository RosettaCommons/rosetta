// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/InterfaceHbondsFilter.hh
/// @brief Counts the number of cross interface hydrogen bonds
/// @author Longxing Cao (longxing@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_InterfaceHbondsFilter_hh
#define INCLUDED_protocols_simple_filters_InterfaceHbondsFilter_hh

#include <protocols/simple_filters/InterfaceHbondsFilter.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreType.hh>


// c++ headers
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>

namespace protocols {
namespace simple_filters {

class InterfaceHbondsFilter : public filters::Filter
{
public:

	InterfaceHbondsFilter();
	InterfaceHbondsFilter(
		core::Size hbonds_num_threshold_in,
		core::Real hbonds_energy_cutoff_in = 0.0,
		core::Size jump_in = 1,
		core::scoring::ScoreFunctionOP scorefxn_in = nullptr,
		bool salt_bridge_mode_in = false,
		bool include_His_chain_terminus_in = false,
		core::Real salt_bridge_distance_cutoff_in = 4.0
	);



	bool apply( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override {
		return utility::pointer::make_shared< InterfaceHbondsFilter >( *this );
	}
	filters::FilterOP fresh_instance() const override{
		return utility::pointer::make_shared< InterfaceHbondsFilter >();
	}
	InterfaceHbondsFilter &
	operator=( InterfaceHbondsFilter const & ot ) {
		if ( this != & ot ) {
			hbonds_num_threshold_ = ot.hbonds_num_threshold_;
			hbonds_energy_cutoff_ = ot.hbonds_energy_cutoff_;
			jump_ = ot.jump_;
			scorefxn_ = ot.scorefxn_;
			salt_bridge_mode_ = ot.salt_bridge_mode_;
			include_His_chain_terminus_ = ot.include_His_chain_terminus_;
			salt_bridge_distance_cutoff_ = ot.salt_bridge_distance_cutoff_;
			charged_res_ = ot.charged_res_;
		}
		return *this;
	}
	InterfaceHbondsFilter( InterfaceHbondsFilter const & init ) :
		Filter( init ) {
		*this = init;
	}
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Size compute( core::pose::Pose const & pose ) const;
	~InterfaceHbondsFilter() override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void set_jump( core::Size const jump_in ) {
		jump_ = jump_in;
	}
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn_in ) {
		scorefxn_ = scorefxn_in;
	}
	void set_hbonds_num_threshold( core::Size hbonds_num_threshold_in ) {
		hbonds_num_threshold_ = hbonds_num_threshold_in;
	}
	void set_hbonds_energy_cutoff( core::Real hbonds_energy_cutoff_in ) {
		hbonds_energy_cutoff_ = hbonds_energy_cutoff_in;
	}
	void set_salt_bridge_mode( bool salt_bridge_mode_in ) {
		salt_bridge_mode_ = salt_bridge_mode_in;
		if ( salt_bridge_mode_ ) {
			initialize_charged_residue_map();
		}
	}
	void set_include_His_chain_terminus( bool include_His_chain_terminus_in ) {
		include_His_chain_terminus_ = include_His_chain_terminus_in;
		if ( salt_bridge_mode_ ) {
			initialize_charged_residue_map();
		}
	}
	void set_salt_bridge_distance_cutoff( bool salt_bridge_distance_cutoff_in ) {
		salt_bridge_distance_cutoff_ = salt_bridge_distance_cutoff_in;
	}

private:
	// the pose energy graph must be updated, so set this as private functions
	// and the compute function guantees this.
	core::Size compute_hbonds( core::pose::Pose const & pose ) const;
	core::Size compute_salt_bridges( core::pose::Pose const & pose ) const;
	void       initialize_charged_residue_map();
private:
	core::Size hbonds_num_threshold_;
	core::Real hbonds_energy_cutoff_;
	core::Size jump_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool salt_bridge_mode_;
	bool include_His_chain_terminus_;
	core::Real salt_bridge_distance_cutoff_;
	// data structure is bad? but this makes the code cleaner.
	enum CHARGE {
		POSITIVE,
		NEGATIVE
	};
	typedef std::pair<std::string, CHARGE> Charged_Group;
	std::unordered_map<std::string, std::vector<Charged_Group > > charged_res_;
};

}
}

#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/seeded_abinitio/SegmentHybridizer.hh
/// @brief repurposing logic and some functions from CartesianHybridze protocols for segment insertions and chain closure
/// @author Eva-Maria Strauch

#ifndef INCLUDED_protocols_seeded_abinitio_SegmentHybridizer_HH
#define INCLUDED_protocols_seeded_abinitio_SegmentHybridizer_HH

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <boost/unordered/unordered_map.hpp>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.hh>

namespace protocols {
namespace seeded_abinitio {

class SegmentHybridizer : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;

public:
	SegmentHybridizer();
	~SegmentHybridizer() override;

	void apply( core::pose::Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	// XRW TEMP  std::string get_name() const override;
	void parse_my_tag(  utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;
	void init();
	void set_scorefunction(core::scoring::ScoreFunctionOP scorefxn_in);
	void hybridize( core::pose::Pose & pose , core::Size insert_pos_start, core::Size insert_pos_stop);
	void apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame );
	void check_and_create_fragments( core::pose::Pose & pose, core::Size insert_start, core::Size insert_stop );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	// to output and score full atom again
	core::scoring::ScoreFunctionOP highres_scorefxn_;
	core::scoring::ScoreFunctionOP lowres_scorefxn_;
	core::scoring::ScoreFunctionOP min_scorefxn_;
	core::scoring::ScoreFunctionOP bonds_scorefxn_;
	core::scoring::ScoreFunctionOP nocst_scorefxn_;

	// for span/segment declarations
	utility::vector1< std::pair < std::string,std::string > > seg_vector_;

	// movemap for minimization
	core::kinematics::MoveMapOP mm_;
	core::kinematics::MoveMapOP extended_mm_;

	// for cartesian alignment
	core::Size cartfrag_overlap_;

	// how much outside of the replaced segment should be remodeled
	core::Size extend_outside_;
	core::Size extend_inside_;
	bool auto_mm_;

	/// fragments parts
	core::Real rms_;
	core::Size nfrags_ ;
	core::Size big_;
	core::Size small_;
	core::fragment::FragSetOP fragments_big_;
	core::fragment::FragSetOP fragments_small_;
	boost::unordered_map<core::Size, core::fragment::Frame> library_;
	bool use_seq_;
	int tries_;
	core::Size mc_cycles_;
	core::Real temp_;
	bool use_frags_;
	core::Size min_cycles_;
	bool all_movable_;
	bool extra_min_;
};


} // seeded_abinitio
} // protocols

#endif /*INCLUDED_protocols_seeded_abinitio_movers_SegmentHybridizer_HH*/

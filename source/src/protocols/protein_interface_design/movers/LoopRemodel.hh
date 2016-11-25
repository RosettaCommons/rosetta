// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/LoopRemodel.hh
/// @brief Header for parseable class to run loop perturbation or refinement between a given loop between start/end (inclusive)
/// @author Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_LoopRemodel_hh
#define INCLUDED_protocols_protein_interface_design_movers_LoopRemodel_hh

#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <protocols/protein_interface_design/movers/LoopRemodel.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class LoopRemodel : public simple_moves::DesignRepackMover
{
public:
	LoopRemodel();
	LoopRemodel(
		std::string const & protocol,
		core::Size const loop_start,
		core::Size const loop_end,
		core::Size const cycles,
		bool const auto_loops,
		//bool const design,
		bool const perturb,
		bool const refine,
		bool const hurry,
		core::scoring::ScoreFunctionOP hires_score,
		core::scoring::ScoreFunctionOP lores_score,
		protocols::loops::LoopsCOP loops,
		core::fragment::FragSetOP frag1,
		core::fragment::FragSetOP frag3,
		core::fragment::FragSetOP frag9
	);
	// various setters and getters
	bool perturb(){ return perturb_; }
	void perturb( bool const setting ) { perturb_ = setting; }
	bool refine(){ return refine_; }
	void refine( bool const setting ) { refine_ = setting; }
	bool hurry(){ return hurry_; }
	void hurry( bool const setting ) { hurry_ = setting; }


	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new LoopRemodel ); }
	void apply( core::pose::Pose & pose ) override;
	// XRW TEMP  virtual std::string get_name() const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	virtual ~LoopRemodel();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string protocol_;
	core::scoring::ScoreFunctionOP hires_score_, lores_score_;
	core::Size loop_start_, loop_end_, cycles_;
	bool auto_loops_, perturb_, refine_, hurry_; // design_,
	protocols::loops::LoopsOP loops_;
	core::fragment::FragSetOP frag1_;
	core::fragment::FragSetOP frag3_;
	core::fragment::FragSetOP frag9_;

	bool pick_loop_frags( protocols::loops::LoopsCOP loops, std::string const & full_seqeuence, std::string const & full_ss );
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_LoopRemodel_HH*/

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/movers/LoopMoverFromCommandLine.hh
/// @brief Header for parseable class to build loops from command line
/// @author Jordan Willis(Jordan.r.willis@vanderbilt.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_LoopMoverFromCommandLine_HH
#define INCLUDED_protocols_protein_interface_design_movers_LoopMoverFromCommandLine_HH

#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/types.hh>
#include <protocols/protein_interface_design/movers/LoopMoverFromCommandLine.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class LoopMoverFromCommandLine : public protocols::simple_moves::DesignRepackMover
{
public:
		LoopMoverFromCommandLine();
		LoopMoverFromCommandLine(
					std::string const protocol,
					bool const perturb,
					bool const refine,
					core::scoring::ScoreFunctionOP & hires_score,
					core::scoring::ScoreFunctionOP & lores_score,
					std::string const loop_file_name,
					protocols::loops::LoopsCOP loops
					);
		// various setters and getters
		bool perturb() const { return perturb_;}
		void perturb( bool const setting ) { perturb_ = setting; }
		bool refine() const { return refine_; }
		void refine( bool const setting ) { refine_ = setting; }
		void refine( std::string const setting ){ string_refine_ = setting; }
//		void string_refine( std::string const setting ){ string_refine_ = setting;}
		std::string string_refine() const{ return string_refine_; }
		void intermedrelax( std::string const setting ){ intermedrelax_ = setting;}
		std::string intermedrelax() const{ return intermedrelax_; }
		void relax( std::string const setting ){ relax_ = setting; }
		std::string relax() const{ return relax_; }
	  void remodel( std::string const setting ){ remodel_ = setting; }
		std::string remodel() const{ return remodel_; }
		protocols::moves::MoverOP clone() const;
		protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new LoopMoverFromCommandLine ); }
		void apply( core::pose::Pose & pose );
		virtual std::string get_name() const;
		void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
		virtual ~LoopMoverFromCommandLine();

private:
		std::string protocol_;
		core::scoring::ScoreFunctionOP hires_score_, lores_score_;
		bool perturb_, refine_;
		std::string loop_file_name_;
		protocols::loops::LoopsOP loops_;
		std::string intermedrelax_, remodel_, relax_, string_refine_;
};

} //moves
} //protocols
}

#endif

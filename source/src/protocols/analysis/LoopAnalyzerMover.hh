// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/analysis/LoopAnalyzerMover.hh
/// @brief LoopAnalyzerMover examines loop structures and packages extra scores into a Job object
/// @author Steven Lewis

#ifndef INCLUDED_protocols_analysis_LoopAnalyzerMover_hh
#define INCLUDED_protocols_analysis_LoopAnalyzerMover_hh

// Unit Headers
#include <protocols/analysis/LoopAnalyzerMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace analysis {

class LoopAnalyzerMover : public protocols::moves::Mover {

public:

	LoopAnalyzerMover( protocols::loops::Loops const & loops, bool const tracer = false );

	~LoopAnalyzerMover() override;

	LoopAnalyzerMover();

	LoopAnalyzerMover( LoopAnalyzerMover const & rhs );

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	LoopAnalyzerMover & operator = ( LoopAnalyzerMover const & rhs );

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	bool reinitialize_for_new_input() const override {return false;} //reset allows this; true here would break loops

public:
	/// @brief apply function will calculate data about the input pose.  It is not intended to modify the pose itself (conformation and energies objects) although it may toss data into the DataCache or a Job object.
	void apply( core::pose::Pose & input_pose ) override;

	// XRW TEMP  std::string get_name() const override;

public: ///////////////////getters, setters/////////////
	/// @brief set loops object, because public setters/getters are a rule
	void set_loops( protocols::loops::LoopsCOP loops );

	/// @brief get loops object, because public setters/getters are a rule
	protocols::loops::LoopsCOP const & get_loops( void ) const ;

	/// @brief set tracer bool, because public setters/getters are a rule
	inline void set_use_tracer( bool tracer ) { tracer_ = tracer; }

	/// @brief get tracer bool, because public setters/getters are a rule
	inline bool get_use_tracer( void ) const { return tracer_; }

	/// @brief Return the total score found from the last apply call
	core::Real get_total_score() const;
	core::Real get_max_rama() const;
	core::Real get_max_omega() const;
	core::Real get_max_pbond() const;
	core::Real get_max_chainbreak() const;

	/// @brief Return the vector of chainbreak scores
	utility::vector1<core::Real>
	get_chainbreak_scores();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	/// @brief reset stored data
	void reset();

	/// @brief places cutpoints in the loops, scores chainbreak, removes cutpoints
	void calculate_all_chainbreaks( core::pose::Pose & pose );

	/// @brief ctor helper: create scorefunction
	void set_sf();

	/// @brief convert loops into positions (must wait until pose, thus not in ctor)
	void find_positions( core::pose::Pose const & pose );

private:
	/// @brief used to store a copy of the input loops
	protocols::loops::LoopsCOP loops_;

	/// @brief output to tracer or PDB/silent file
	bool tracer_;

	/// @brief used to calculate positions to examine - loops +- 1 position are interesting, but vary w/termini, etc
	utility::vector1< core::Size > positions_;

	/// @brief scorefunction used to apply multiple individual terms at once, not as a cohesive unit
	core::scoring::ScoreFunctionOP sf_;

	/// @brief scorefunction for chainbreak score
	core::scoring::ScoreFunctionOP chbreak_sf_;

	///brief remember chainbreak scores
	utility::vector1< core::Real > scores_;

	core::Real total_score_;
	core::Real max_rama_;
	core::Real max_chainbreak_;
	core::Real max_omega_;
	core::Real max_pbond_;

}; //class LoopAnalyzerMover

}//analysis
}//protocols

#endif //INCLUDED_protocols_analysis_LoopAnalyzerMover_HH

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rosetta_scripts/ParsedProtocol.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_rosetta_scripts_ParsedProtocol_HH
#define INCLUDED_protocols_rosetta_scripts_ParsedProtocol_HH

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobOutputterObserver.hh>

#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/exit.hh>
#include <protocols/moves/ResId.hh>

// C++ headers
#include <string>

// Unit headers
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace rosetta_scripts {

class ParsedProtocol :
		public protocols::moves::Mover,
		public protocols::moves::ResId,
		public protocols::jd2::JobOutputterObserver
{
public:
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef std::pair< std::pair< protocols::moves::MoverOP, std::string >, protocols::filters::FilterOP > mover_filter_pair;
	typedef utility::vector1< mover_filter_pair > MoverFilterVector;
	typedef MoverFilterVector::iterator iterator;
	typedef MoverFilterVector::const_iterator const_iterator;
public:
	ParsedProtocol();
	virtual ~ParsedProtocol();
	virtual void apply( Pose & pose );
	virtual core::pose::PoseOP get_additional_output( );
	virtual std::string get_name() const;
	/// Tricky! movers are cloned into the protocol b/c their apply functions (which are nonconst) could accumulate
	/// state information. Filters are safe and are therefore merely registered.
	/// Under this state of affairs, a mover or filter may be called many times in the protocol, and it will be
	/// guaranteed to have no state accumulation.
	void add_mover( protocols::moves::MoverCOP mover, std::string const mover_name, protocols::filters::FilterOP filter );
	void final_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );
	core::scoring::ScoreFunctionCOP final_scorefxn() const;
	void final_score(core::pose::Pose & pose) const;
	void report_all( Pose const & pose ) const; // cycles over all filter->report methods to output their values to a common stream.
	void report_filters_to_job( Pose const & pose ) const;  // as above but reports to job object
	//as above but is called directly from JobOutputter via Observer pattern
	virtual
	void add_values_to_job( Pose const & pose, protocols::jd2::JobOP ) const;

	void report_all_sm( std::map< std::string, core::Real > & score_map, Pose const & pose ) const; // ditto, but outputs filter values into score_map object
	protocols::moves::MoverCOP get_mover( core::Size const mover_number ) const {
		runtime_assert( movers_.size() >= mover_number && mover_number > 0 );
		return( movers_[ mover_number ].first.first );
	}
	void set_resid( core::Size const resid );
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new ParsedProtocol ); }
	virtual void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ); // this is defined as public here, b/c I need to circumvent the name-check, since this is called both by the Movers section (as ParsedProtocol) and the PROTOCOLS section.
	void clear() { movers_.clear(); }
	std::string mode() const{ return mode_; }
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	void apply_probability( utility::vector1< core::Real > const a );
	utility::vector1< core::Real > apply_probability();
	core::Size size() { return movers_.size(); }
	core::Size last_attempted_mover_idx() { return last_attempted_mover_idx_; }
	void last_attempted_mover_idx( core::Size const s ){ last_attempted_mover_idx_ = s;}
	bool report_call_order() const { return report_call_order_; }
	void report_call_order( bool const c ) { report_call_order_ = c; }
	std::string call_order() const{ return call_order_; }
private:
	void finish_protocol(Pose & pose);
	bool apply_mover_filter_pair(Pose & pose, mover_filter_pair const & mover_pair);
	bool apply_filter(Pose & pose, mover_filter_pair const & mover_pair);

	void sequence_protocol(Pose & pose, utility::vector1< mover_filter_pair >::const_iterator mover_it_in);
	void random_order_protocol(Pose & pose);
	void random_single_protocol(Pose & pose);
private:

	MoverFilterVector movers_;
	core::scoring::ScoreFunctionCOP final_scorefxn_;
	std::string mode_;
	utility::vector1< core::Real > apply_probability_; // if mode_="single_random", assigns a probability of execution to each mover/filter pair. Defaults to equal probabilities to all.
	core::Size last_attempted_mover_idx_; //index to last attempted mover; useful for adaptive monte carlo
	bool report_call_order_; //dflt false; At the end of the run, write to out the sequence of mover/filter calls (good for stochastic application
	std::string call_order_; // saved call order, not writeable
	protocols::moves::MoverOP last_mover_;
	bool resume_support_;

};

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_ParsedProtocol_HH

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/FilterScan.hh
/// @brief Scans a task factory, mutates to each designable residue and evaluates a filter. Mutations that pass the filter are output in a resfile format
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_FilterScan_hh
#define INCLUDED_protocols_protein_interface_design_filters_FilterScan_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/protein_interface_design/filters/FilterScan.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1.hh>
#include <core/chemical/AA.hh>
#include <protocols/simple_filters/DeltaFilter.fwd.hh>

// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace filters {

class FilterScanFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	FilterScanFilter();
	virtual ~FilterScanFilter();

	/// @brief Constructor with a single target residue
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	/// Undefined, commenting out to fix PyRosetta build  core::Real compute( core::pose::Pose const & pose ) const;
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	void resfile_name( std::string const & resfile_name);
	void score_log_file( std::string const & score_log_file);
	std::string resfile_name() const;
	std::string score_log_file() const;
	protocols::filters::FilterOP triage_filter() const;
	void triage_filter( protocols::filters::FilterOP filter );

	protocols::filters::FilterOP filter() const;
	void filter( protocols::filters::FilterOP filter );
	std::string resfile_general_property() const;
	void resfile_general_property( std::string const & );
	void relax_mover( protocols::moves::MoverOP  mover );
	protocols::moves::MoverOP relax_mover() const;
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	bool unbound() const;
	void unbound( bool const u );
	core::Size jump() const;
	void jump( core::Size const j );
	bool delta() const;
	void delta( bool const d );
	bool report_all() const;
	void report_all( bool const ra );
	void dump_pdb( bool const d );
	bool dump_pdb() const;
	bool rtmin() const{ return rtmin_; }
	void rtmin( bool const r ){ rtmin_ = r; }
	void single_substitution( core::pose::Pose & pose, core::Size const resi, core::chemical::AA const target_aa ) const;
	utility::vector1< protocols::simple_filters::DeltaFilterOP > delta_filters() const;
	void delta_filters( utility::vector1< protocols::simple_filters::DeltaFilterOP > const & d );

	std::string dump_pdb_name() const{ return dump_pdb_name_; }
	void dump_pdb_name( std::string const & s ){ dump_pdb_name_ = s; }
	bool keep_native() const{ return keep_native_; }
	void keep_native( bool const k ){ keep_native_ = k; }

	utility::vector1< core::Real > delta_filter_thresholds() const;
	void delta_filter_thresholds( utility::vector1< core::Real > const & );
private:
	core::pack::task::TaskFactoryOP task_factory_;
	protocols::filters::FilterOP triage_filter_;//dflt null; mutations that are categorically rejected, no matter what
	protocols::filters::FilterOP filter_;//dflt null; a filter to use for its report functionality to probe the pose's state
	std::string score_log_file_;
	std::string resfile_name_;
	std::string resfile_general_property_; //dflt nataa; what to write in the resfile above the 'start' line
	protocols::moves::MoverOP relax_mover_; //dflt nullmover; what to do after mutation
	core::scoring::ScoreFunctionOP scorefxn_; // which scorefxn to use during packing
	bool delta_; // dflt false; compute as delta? If true, all values are reported relative to the baseline input pose's filter evaluation.
	bool unbound_;
	bool report_all_;
	core::Size jump_;
	void unbind( core::pose::Pose & ) const; //utility function for unbinding the pose
	bool dump_pdb_; // dflt false; dump a pdb for each substitution (with extensions signifying the substitution).
	bool rtmin_; //dflt false; shall we do rtmin after each mutation (and at baseline)?
	utility::vector1< protocols::simple_filters::DeltaFilterOP > delta_filters_;
	std::string dump_pdb_name_; //dflt "" ; give a user-defined name to the dumped-pdbs
	bool keep_native_; // dflt false ; always write the native residue into the resfile
	utility::vector1< core::Real > delta_filter_thresholds_; //dflt empty; in case you want to write resfiles for a variety of energy cutoffs
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_FilterScanFilter_HH_

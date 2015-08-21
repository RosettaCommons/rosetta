// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  A filter that, for each dimer in a pose, outputs the peptide which contributes most to the interface.

/// @author Nir London
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @date   Jul. 30, 2014

#ifndef INCLUDED_protocols_analysis_PeptideDeriverFilter_hh
#define INCLUDED_protocols_analysis_PeptideDeriverFilter_hh

// Unit headers
#include <protocols/analysis/PeptideDeriverFilter.fwd.hh>

// Project headers
#include <basic/options/keys/OptionKeys.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/filters/Filter.hh>
#include <utility/io/orstream.hh>

// RosettaScripts header
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ headers
#include <iosfwd>
#include <string>

namespace protocols {
namespace analysis {

/// @brief an enum specifying the type of peptide entry, when calling PeptideDeriverOuputter::peptide_entry()
enum PeptideDeriverEntryType {
	/// a general peptide from the sliding window run
	ET_GENERAL,
	/// the peptide with the best interface score (as calculated for the linear peptide)
	ET_BEST_LINEAR,
	/// the peptide with the best interface score (as calculated for the cyclic peptide)
	ET_BEST_CYCLIC
};

/// @brief the output format that PeptideDeriverFilter reports in
enum PeptideDeriverReportFormat {
	/// a pretty, readable but verbose format (Markdown)
	PRF_MARKDOWN,
	/// a stripped-down, easily parsable format
	PRF_BASIC
};

/// @brief an abstract base class for handling calculation outputs from
///        PeptideDeriverFilter
/// Since PeptideDeriverFilter might have a set of outputs for each residue,
/// for each chain in each chain-pair, outputting is quite elaborate. This is
/// an attempt to decouple the calculation from the representation of results.
/// Representation of results is delegated to implementors of this class.
class PeptideDeriverOutputter {
public:
	/// @brief called by PeptideDeriverFilter when processing of a strucuture (possibly multi-chain) starts
	virtual void begin_structure(core::pose::Pose const &, std::string const &) = 0;

	/// @brief called by PeptideDeriverFilter when a chain-pair pose is prepared (minimized and disulfide-bridge resolved)
	virtual void chain_pair_pose_prepared(core::pose::Pose const & pose) = 0;

	/// @brief called by PeptideDeriverFilter when calculation commences on a receptor-partner pair
	virtual void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const total_isc,
		std::string const & options_string) = 0;

	/// @brief called by PeptideDeriverFilter when calculation commences for
	///        the specified peptide length
	/// this is called for every peptide length specified to
	/// PeptideDeriverFilter, per each call to
	/// PeptideDeriverOutputter::begin_receptor_partner_pair()
	virtual void peptide_length(core::Size const pep_length) = 0;

	/// @brief called by PeptideDeriverFilter for each 'peptide' entry that should be output.
	/// besides the peptides for each sliding window, this is being called for the 'best' peptides (@see PeptideDeriverEntryType)
	virtual void peptide_entry(core::pose::Pose const & pose, PeptideDeriverEntryType const entry_type, core::Size const pep_start,
		core::Real const linear_isc, std::string const & disulfide_info, bool const was_cyclic_pep_modeled,
		core::pose::Pose const & cyclic_pose, core::Real const cyclic_isc) = 0;

	/// @brief called by PeptideDeriverFilter when calculation concludes for a receptor-partner pair (for all the different peptide lengths)
	virtual void end_receptor_partner_pair() = 0;

	/// @brief called by PeptideDeriverFilter when processing of a strucuture (all chain-pairs considered) ends
	virtual void end_structure() = 0;
};

/// @brief a container class which holds a list of PeptideDeriverOutputter instances and delegates calls to those.
class PeptideDeriverOutputterContainer : public PeptideDeriverOutputter {
public:
	/// @brief add an outputter to the list.
	void push_back(PeptideDeriverOutputterOP const item);

	/// @brief clear all outputters in the list.
	void clear();

	virtual void begin_structure(core::pose::Pose const & pose, std::string const &);

	virtual void chain_pair_pose_prepared(core::pose::Pose const & pose);

	virtual void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const total_isc,
		std::string const & options_string);

	virtual void peptide_length(core::Size const pep_length);

	virtual void peptide_entry(core::pose::Pose const & pose,
		PeptideDeriverEntryType const entry_type, core::Size const pep_start,
		core::Real const linear_isc, std::string const & disulfide_info,
		bool const was_cyclic_pep_modeled, core::pose::Pose const & cyclic_pose,
		core::Real const cyclic_isc);

	virtual void end_receptor_partner_pair();

	virtual void end_structure();

private:
	utility::vector1<PeptideDeriverOutputterOP> list_;

}; // PeptideDeriverOutputterContainer

// TODO : find a nice way to make a base class for this and the Markdown outputter
/// @brief outputs a Peptiderive report to a stream in a basic, easily parsable format
class PeptideDeriverBasicStreamOutputter : public PeptideDeriverOutputter {
public:
	/// @param out    an output stream to which output will be diverted
	///               should belong to the thread on which the outputter
	///               methods are called in, and outlive the outputter object.
	/// @param prefix a string prefix to prepend to each output line
	PeptideDeriverBasicStreamOutputter(utility::io::orstream & out, std::string prefix);

	virtual void begin_structure(core::pose::Pose const &, std::string const &);

	virtual void chain_pair_pose_prepared(core::pose::Pose const & ) {
		// do nothing
	}
	virtual void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const total_isc,
		std::string const & options_string);

	virtual void peptide_length(core::Size const pep_length);

	// TODO : consider using a PeptideDeriverEntry struct to make this signature less verbose and more tolerant to changes in the protocol
	//        e.g.,  virtual void peptide_entry(core::pose::Pose const & , PeptideDeriverEntry const & peptide_deriver_entry) --yuvals
	virtual void peptide_entry(core::pose::Pose const & , PeptideDeriverEntryType const entry_type, core::Size const pep_start,
		core::Real const linear_isc, std::string const & disulfide_info, bool const was_cyclic_pep_modeled,
		core::pose::Pose const &, core::Real const cyclic_isc);

	virtual void end_receptor_partner_pair();

	virtual void end_structure();

private:
	utility::io::orstream * out_p_;
	std::string prefix_;

}; // PeptideDeriverBasicStreamOutputter

// TODO : add an option to allow output in this format
// TODO : allow easy switching with reST format
// TODO : conside creating a separate class(es) for simple Markdown/reST/HTML constructs (headers, tables, lists)
/// @brief outputs a Markdown formatted Peptiderive report to a stream
class PeptideDeriverMarkdownStreamOutputter : public PeptideDeriverOutputter {
public:

	/// @param out    an output stream to which output will be diverted
	///               should belong to the thread on which the outputter
	///               methods are called in, and outlive the outputter object.
	/// @param prefix a string prefix to prepend to each output line
	PeptideDeriverMarkdownStreamOutputter(utility::io::orstream & out, std::string prefix);

	virtual void begin_structure(core::pose::Pose const & pose, std::string const &);

	// TODO : do we want to output anything here? Probably not.
	virtual void chain_pair_pose_prepared(core::pose::Pose const &) {
		// do nothing
	}

	virtual void end_structure();

	virtual void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const total_isc,
		std::string const & options_string);

	virtual void peptide_length(core::Size const pep_length);

	virtual void peptide_entry(core::pose::Pose const & pose,
		PeptideDeriverEntryType const entry_type, core::Size const pep_start,
		core::Real const linear_isc, std::string const & disulfide_info,
		bool const was_cyclic_pep_modeled, core::pose::Pose const &,
		core::Real const cyclic_isc );

	virtual void end_receptor_partner_pair();

private:
	/// clear accumulating member variables
	void clear_buffers();

	utility::io::orstream * out_p_;
	std::string prefix_;

	core::Size current_pep_length_;
	core::Real current_total_isc_;
	char current_receptor_chain_letter_;
	char current_partner_chain_letter_;
	std::string current_sequence_;

	std::ostringstream header_;
	std::ostringstream best_linear_peptides_;
	std::ostringstream best_cyclic_peptides_;
	std::ostringstream all_peptides_;
	std::ostringstream footer_;

}; // PeptideDeriverMarkdownStreamOutputter

/// @brief outputs poses at different points of the Peptiderive protocol, according to the given set of options
class PeptideDeriverPoseOutputter : public PeptideDeriverOutputter {
public:
	/// @param is_dump_best_peptide_pose flag to output the best (in terms of linea/cyclic isc) poses. Max. 4 poses per chain per.
	/// @param is_dump_prepared_pose     flag to output each chain-pair pose as Peptiderive sees it
	/// @param is_dump_cyclic_poses      flag to output each modeled cyclic peptide (@see basic::options::OptionKeys::peptide_deriver::optimize_cyclic_threshold)
	/// @param scorefxn                  score function to score the pose with before outputting
	PeptideDeriverPoseOutputter( bool const is_dump_best_peptide_pose,
		bool const is_dump_prepared_pose, bool const is_dump_cyclic_poses,
		core::scoring::ScoreFunctionCOP scorefxn);

	virtual void begin_structure(core::pose::Pose const &, std::string const &) {
		// do nothing
	}

	virtual void chain_pair_pose_prepared(core::pose::Pose const & pose);

	virtual void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const,
		std::string const &);

	virtual void peptide_length(core::Size const pep_length);

	virtual void end_receptor_partner_pair() {
		// do nothing
	}

	virtual void end_structure() {
		// do nothing
	}

	virtual void peptide_entry(core::pose::Pose const & pose,
		PeptideDeriverEntryType const entry_type, core::Size const,
		core::Real const, std::string const & disulfide_info,
		bool const was_cyclic_pep_modeled, core::pose::Pose const & cyclic_pose,
		core::Real const);

private:
	/// output a post to the given file
	void output_pose( core::pose::Pose & pose, std::string const & pose_name );

	bool is_dump_best_peptide_pose_;
	bool is_dump_prepared_pose_;
	bool is_dump_cyclic_poses_;

	char current_receptor_chain_letter_;
	char current_partner_chain_letter_;
	core::Size current_peptide_length_;
	core::pose::PoseOP current_chain_pair_pose_;
	bool is_chain_pair_new_;
	core::scoring::ScoreFunctionCOP scorefxn_;

}; // PeptideDeriverPoseOutputter


/// @brief implementation of the Peptiderive protocol - given a pose, it
///        evaluates the binding energy of linear stretches, of given lengths,
///        and indicates which of these contributes more than others to the
///        interaction energy. It also reports on and optionally models
///        cyclizable peptides.
///
/// Data is either output to a log (using report()) or passed to the given
/// PeptideDeriverOutputter (using derive_peptide()).
///
/// When calling report(), the output format is determined by a list of flags:
/// - is_dump_peptide_pose
/// - is_dump_cyclic_poses
/// - is_dump_prepared_pose
/// - is_dump_report_file
/// - is_report_gzip
/// See each flag's documentation for their effect. These flags are settable
/// via command-line or RosettaScripts options, and may also be set
/// programmatically. These only affect the output coming from report(), since
/// report() uses these to determine which PeptideDeriverOutputter instances
/// it passes to derive_peptide().
///
/// When calling derive_peptide(), however, it is up to the consumer to provide
/// the PeptideDeriverOutputter instance of choice. Naturally, it may use any
/// of the subclasses provided in this file:
/// - PeptideDeriverOutputterContainer
/// - PeptideDeriverBasicStreamOutputter
/// - PeptideDeriverMarkdownStreamOutputter
/// - PeptideDeriverPoseOutputter
///
/// Note: the PeptideDeriverFilter doesn't actually filter anything, so it
/// might have been better to suffix it with -Reporter, but wasn't so not to
/// confuse it with the FeatureReporter class hierarchy.
///
class PeptideDeriverFilter : public protocols::filters::Filter {
public:
	/// @brief default constructor
	PeptideDeriverFilter();

	/// @brief copy constructor
	PeptideDeriverFilter(PeptideDeriverFilter const &);

	/// @brief assignment operator
	PeptideDeriverFilter & operator=( PeptideDeriverFilter const & rhs );

	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;

	// inline overrides
	virtual protocols::filters::FilterOP clone() const {
		return protocols::filters::FilterOP( new PeptideDeriverFilter( *this ) );
	}
	virtual protocols::filters::FilterOP fresh_instance() const {
		return protocols::filters::FilterOP( new PeptideDeriverFilter );
	}

	// NOTE : I have no idea why this virtual method exists together with the get_type() method, which
	//        returns the type_ member initialized upon construction (see constructor).
	virtual std::string name() const { return "PeptideDeriverFilter"; }

	// accessors
	utility::vector1<core::Size> get_pep_lengths() const { return pep_lengths_; }
	void set_pep_lengths(utility::vector1<core::Size> const & values) { pep_lengths_ = values; }

	utility::vector1<char> get_restrict_receptors_to_chains() const { return restrict_receptors_to_chains_; }
	void set_restrict_receptors_to_chains(utility::vector1<char> const & values) { restrict_receptors_to_chains_ = values; }

	utility::vector1<char> get_restrict_partners_to_chains() const { return restrict_partners_to_chains_; }
	void set_restrict_partners_to_chains(utility::vector1<char> const & values) { restrict_partners_to_chains_ = values; }

	PeptideDeriverReportFormat get_report_format() const { return report_format_; }
	void set_report_format(PeptideDeriverReportFormat const value) { report_format_ = value; }
	static PeptideDeriverReportFormat parse_report_format_string(std::string const value);

	bool get_is_skip_zero_isc() const { return is_skip_zero_isc_; }
	void set_is_skip_zero_isc(bool const value) { is_skip_zero_isc_ = value; }

	bool get_is_dump_peptide_pose() const { return is_dump_peptide_pose_; }
	void set_is_dump_peptide_pose(bool const value) { is_dump_peptide_pose_ = value; }

	bool get_is_dump_cyclic_poses() const { return is_dump_cyclic_poses_; }
	void set_is_dump_cyclic_poses(bool const value) { is_dump_cyclic_poses_ = value; }

	bool get_is_dump_prepared_pose() const { return is_dump_prepared_pose_; }
	void set_is_dump_prepared_pose(bool const value) { is_dump_prepared_pose_ = value; }

	bool get_is_dump_report_file() const { return is_dump_report_file_; }
	void set_is_dump_report_file(bool const value) { is_dump_report_file_ = value; }

	bool get_is_report_gzip() const { return is_report_gzip_; }
	void set_is_report_gzip(bool const value) { is_report_gzip_ = value; }

	bool get_is_do_minimize() const { return is_do_minimize_; }
	void set_is_do_minimize(bool const value) { is_do_minimize_ = value; }

	core::Real get_optimize_cyclic_threshold() const {return optimize_cyclic_threshold_; }
	void set_optimize_cyclic_threshold(core::Real const value) { optimize_cyclic_threshold_ = value; }

	core::scoring::ScoreFunctionOP get_scorefxn_deriver() const { return scorefxn_deriver_; }
	void set_scorefxn_deriver(core::scoring::ScoreFunctionOP value) { scorefxn_deriver_ = value; }

	// NOTE : scorefxn_minimizer_ has no accessor because that seems less of something
	//        we'd like to allow chaging, as we impose constraints on CAs

	// RosettaScripts implementation
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	// Filter-specific methods

	/// @brief Prepare pose before derivation, by minimizing it and finding disulfides
	void prepare_pose(
		PeptideDeriverOutputter & output,
		core::pose::Pose & pose) const;

	/// @brief Calculate the interface score, i.e. addition to energy
	///        resulting from the complex coming together.
	core::Real calculate_interface_score(core::pose::Pose const & pose, core::Size const jump_id) const;

	// @brief Given a two-monomer pose, calculates each the interface score (delta between energy of complex and of the seperated monomers) is calculated for each residue.
	void calculate_per_residue_interface_score(core::pose::Pose & chain_pair_pose, std::set<core::Size> interface_residues, std::ostream & report_out, core::Size const jump_id) const;

	// @brief given two protein chains, build one pose containing the first pose and a peptide cut from the second pose
	core::pose::PoseOP build_receptor_peptide_pose (core::pose::Pose const & receptor_pose, core::pose::Pose const & partner_pose, core::Size peptide_start, core::Size peptide_end, core::Size & jump_id) const;


	/// @brief overload of PeptideDeriverFilter::derive_peptide() which accepts
	///        a many-chain pose and chain indices, rather than single-chain poses.
	void derive_peptide(
		PeptideDeriverOutputter & output,
		core::pose::Pose const & pose,
		core::Size const first_chain_index,
		core::Size const second_chain_index,
		bool const both_ways ) const;

	/// @brief The core of the filter: calculate and output a table of the
	///        energetic contribution of each linear segment of the
	///        pep_chain_pose chain (going through it with a sliding window),
	///        together with other interesting data.
	void derive_peptide(
		PeptideDeriverOutputter & output,
		core::pose::Pose const & receptor_pose,
		core::pose::Pose const & partner_pose,
		core::Real const total_isc ) const;

private:
	/// @brief Minimization step, run on an input once before derivation.
	void minimize(core::pose::Pose &) const;

	/// @brief common implementation of cctor and operator=
	static void assign(PeptideDeriverFilter &, PeptideDeriverFilter const &);

	/// @brief helper function, creates a vector containing the range of residues indices that span the given chain
	static void push_chain_residue_indices(core::pose::Pose const & many_chain_pose, core::Size const chain_index, utility::vector1<core::Size> & residue_indices );

	/// @brief helper function, returns a list of chain indices. If restrict_to_chains is specified, only indices of chains with these letters
	///        are retuned. Otherwise, indices of all chains (i.e., 1..num_chains) are returned.
	static utility::vector1<core::Size> get_chain_indices( core::pose::Pose const & pose, utility::vector1<char> const & restrict_to_chains );

	/// @brief read options from option system
	void parse_options();

	// @brief For a two-monomer pose, returns a set of interface residues
	static std::set<core::Size> find_interface_residues (core::pose::Pose const & chain_pair_pose);

	/// @brief return zero when isc is close enough to zero; used to prevent outputs such as alternating -0 and 0.
	static core::Real normalize_isc(core::Real const isc);

	// private fields

	/// @brief the window size(s) of peptides to derive
	utility::vector1< core::Size > pep_lengths_;

	/// @brief whether to optimize by jumping over the windows in which the total Isc = 0
	bool is_skip_zero_isc_;

	/// @brief whether to output the pose of the each (best) derived peptide. Only affects the way report() does output.
	bool is_dump_peptide_pose_;

	/// @brief whether to output each of the modelled cyclic poses (only cyclizable poses with enough contribution, i.e. beyond optimize_cyclic_threshold_, are modelled). Only affects the way report() does output.
	bool is_dump_cyclic_poses_;

	/// @brief whether to output the pose of the model after minimization. Only affects the way report() does output.
	bool is_dump_prepared_pose_;

	/// @brief the format in which PeptideDeriver outputs. See PeptideDeriverReportFormat for options. Only affects the way report() does output.
	PeptideDeriverReportFormat report_format_;

	/// @brief whether to output the peptide information report to a separate file (rather than the log). Only affects the way report() does output.
	bool is_dump_report_file_;

	/// @brief whether to gzip report file; only relevant when is_dump_report_file is true. Only affects the way report() does output.
	bool is_report_gzip_;

	/// @brief whether to minimize the structure before starting the algorithm
	bool is_do_minimize_;

	/// @brief set a threshold from which a peptide contributing such percent of total binding energy will be optimized
	core::Real optimize_cyclic_threshold_;

	/// @brief list of chains to use as receptors. When empty, all chains are included.
	///        Peptides won't be derived from these chains (unless they are also specified as partners).
	/// @see restrict_partners_to_chains_
	utility::vector1<char> restrict_receptors_to_chains_;

	/// @brief list of chain to use as partners. When empty, all chains are included.
	///        We only go over the chains listed as partners to derive the peptides.
	/// @see restrict_receptors_to_chains_
	utility::vector1<char> restrict_partners_to_chains_;

	/// @brief The score function used to score derived peptides
	core::scoring::ScoreFunctionOP scorefxn_deriver_;

	/// @brief The score function used to minimize the initial pose
	core::scoring::ScoreFunctionOP scorefxn_minimizer_;

};


} // namespace analysis
} // namespace protocols

#endif

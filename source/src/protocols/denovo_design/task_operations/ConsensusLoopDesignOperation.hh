// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/task_operations/ConsensusLoopDesignOperation.hh
/// @brief Restrict loop residue identities to those commonly found in nature
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_task_operations_ConsensusLoopDesignOperation_hh
#define INCLUDED_protocols_denovo_design_task_operations_ConsensusLoopDesignOperation_hh

// unit headers
#include <protocols/denovo_design/task_operations/ConsensusLoopDesignOperation.fwd.hh>

// protocol headers

// project headers

// core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/types.hh>

// utility Headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/SingletonBase.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace task_operations {

struct SurroundingSS {
	SurroundingSS() :
		before( char( 0 ) ),
		after( char( 0 ) ) {};

	SurroundingSS( char const bef, char const aft ) :
		before( bef ),
		after( aft ) {};

	char before;
	char after;

	bool
	operator<( SurroundingSS const & ab2 ) const;

	friend std::ostream &
	operator<<( std::ostream & os, SurroundingSS const & ss );
};

class AAFrequency : public utility::pointer::ReferenceCount {
public:
	AAFrequency( core::Real const frequency, core::Real const enrichment );

	core::Real
	frequency() const;

	core::Real
	enrichment() const;

private:
	char aa_;
	core::Real frequency_;
	core::Real enrichment_;

private:
	AAFrequency();
};
typedef utility::pointer::shared_ptr< AAFrequency const > AAFrequencyOP;
typedef utility::pointer::shared_ptr< AAFrequency const > AAFrequencyCOP;

class AAFrequencies : public utility::pointer::ReferenceCount {
public:
	typedef std::map< char, AAFrequencyCOP > AAFrequencyMap;

public:
	AAFrequencies();

	AAFrequency const &
	frequency( char const aa ) const;

	bool
	has_frequency( char const aa ) const;

	void
	set_frequency( char const aa, AAFrequencyCOP frequency );

	AAFrequencyMap::const_iterator
	begin() const;

	AAFrequencyMap::const_iterator
	end() const;

private:
	AAFrequencyMap aa_freq_;
};
typedef utility::vector1< AAFrequencies > LoopAAs;

class ConsensusLoopDatabase : public utility::SingletonBase< ConsensusLoopDatabase > {
public:
	typedef std::string Abego;
	typedef std::map< Abego, LoopAAs > AbegoToLoopAAsMap;
	typedef std::map< SurroundingSS, AbegoToLoopAAsMap > SurroundingSSMap;

public:
	ConsensusLoopDatabase();

	virtual ~ConsensusLoopDatabase();

	static ConsensusLoopDatabase *
	create_singleton_instance();

public:
	void
	set_frequency(
		SurroundingSS const & surrounding_ss,
		Abego const & abego,
		core::Size const loop_resid,
		char const aa,
		AAFrequencyCOP const frequency );

	AAFrequencies const &
	frequencies(
		SurroundingSS const & surrounding_ss,
		Abego const & abego,
		core::Size const loop_resid ) const;

	AAFrequency const &
	frequency(
		SurroundingSS const & surrounding_ss,
		Abego const & abego,
		core::Size const loop_resid,
		char const aa ) const;

	bool
	has_frequencies(
		SurroundingSS const & surrounding_ss,
		Abego const & abego ) const;

	bool
	has_frequencies(
		SurroundingSS const & surrounding_ss,
		Abego const & abego,
		core::Size const loop_resid ) const;

	bool
	has_frequency(
		SurroundingSS const & surrounding_ss,
		Abego const & abego,
		core::Size const loop_resid,
		char const aa ) const;

private:
	void
	read_db();

private:
	SurroundingSSMap ss_map_;
};

struct LoopInfo {
	SurroundingSS ss_around;
	core::Size startres;
	std::string abego;
};
std::ostream & operator<<( std::ostream & os, LoopInfo const & info );
typedef utility::vector1< LoopInfo > LoopInfoVec;

class ConsensusLoopDesignOperation : public core::pack::task::operation::TaskOperation {
public:
	typedef std::string AAs;

public:
	/// @brief default constructor
	ConsensusLoopDesignOperation();

	/// @brief destructor
	virtual ~ConsensusLoopDesignOperation();

	static std::string
	class_name();

	/// @brief make clone
	virtual core::pack::task::operation::TaskOperationOP
	clone() const;

	/// @brief apply
	virtual void
	apply(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask & task ) const;

	/// @brief Returns the name of the class
	virtual std::string
	get_name() const;

public:
	// mutators

	/// @brief Sets the secondary structure to be used for selection of loops.
	///        If a residue selector is set, this is ignored.
	/// @param[in] secstruct  Secondary structure to be used. Must match pose length
	void
	set_secstruct( std::string const & secstruct );

	/// @brief Sets residue selector. If set, only selected residues will
	///        be operated upon. (default = select by secondary structure type L)
	/// @param[in] selector_val  Residue selector to be used.  Cloned by this function and
	///                          set as selector_
	void
	set_residue_selector( core::select::residue_selector::ResidueSelector const & selector_val );

public:
	virtual void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map );

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	AAFrequencies const &
	aa_frequencies(
		LoopInfo const & info,
		core::Size const resid ) const;

	AAs
	forbidden_aas( AAFrequencies const & frequencies ) const;

	void
	disallow_aas(
		core::pack::task::PackerTask & task,
		LoopInfo const & info ) const;

	LoopInfoVec
	get_loop_info( core::pose::Pose const & pose ) const;

	LoopInfoVec
	loop_info_from_subset(
		core::pose::Pose const & pose,
		std::string const & ss,
		core::select::residue_selector::ResidueSubset const & subset ) const;

	/// @brief if true, residues adjacent to loops will be restricted. Otherwise, just the loop. (default=false)
	void
	set_include_adjacent_residues( bool const include_res );

	/// @brief Reads/parses a blueprint file and calls set_secstruct() on the
	///        resulting secondary structure.  Length of blueprint file must
	///        match the length of the input pose.
	void
	set_secstruct_from_blueprint( std::string const & bp_file );

	void
	set_enrichment_threshold( core::Real const threshold );

private:
	/// @brief Gets secondary structure to be used for determining what is surrounding
	///        the loops
	/// @param[in] pose input pose
	/// @returns Secondary structure string according to the following rules:
	///          1. If secstruct_ is set, return that.
	///          2. If use_dssp_ is set, return string computed from DSSP.
	///          3. return pose.secstruct()
	std::string
	get_secstruct( core::pose::Pose const & pose ) const;

	/// @brief Returns a residue selector to be used to choose loops
	/// @param[in] secstruct  Secondary structure to be used, if default selector is used
	/// @details If selector_ is provided, simply return that.  If not,
	///          create a default secondary structure selector using the given secstruct
	core::select::residue_selector::ResidueSelectorCOP
	residue_selector( std::string const & secstruct ) const;

	AAs
	compute_aas_after_disallowing(
		AAs const & aas,
		core::pack::task::ResidueLevelTask const & task ) const;

	AAs
	compute_best_allowed_aas(
		AAFrequencies const & aa_freqs,
		core::pack::task::ResidueLevelTask const & task ) const;

private:
	std::string secstruct_;
	bool include_adjacent_residues_;
	bool use_dssp_;
	core::Real enrichment_threshold_;
	core::select::residue_selector::ResidueSelectorCOP selector_;
};

} // task_operations
} // denovo_design
} // protocols

#endif

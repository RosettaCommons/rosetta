// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/denovo_design/task_operations/ConsensusLoopDesign.hh
/// @brief Restrict loop residue identities to those commonly found in nature
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_task_operations_ConsensusLoopDesign_hh
#define INCLUDED_protocols_denovo_design_task_operations_ConsensusLoopDesign_hh

// unit headers
#include <protocols/denovo_design/task_operations/ConsensusLoopDesign.fwd.hh>

// protocol headers

// project headers

// core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/types.hh>

// utility Headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

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
	bool operator<( SurroundingSS const & ab2 ) const;
};
typedef std::list< std::string > LoopAAs;
typedef std::map< std::string, LoopAAs > AbegoToSequenceMap;
typedef std::map< SurroundingSS, AbegoToSequenceMap > ConsensusSequenceTable;

struct LoopInfo {
	SurroundingSS ss_around;
	core::Size startres;
	std::string abego;
};
std::ostream & operator<<( std::ostream & os, LoopInfo const & info );
typedef utility::vector1< LoopInfo > LoopInfoVec;

class ConsensusLoopDesignOperation : public core::pack::task::operation::TaskOperation {
public:
	/// @brief default constructor
	ConsensusLoopDesignOperation();

	/// @brief destructor
	virtual ~ConsensusLoopDesignOperation();

	/// @brief make clone
	virtual core::pack::task::operation::TaskOperationOP clone() const;

	/// @brief apply
	virtual void apply(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask & task ) const;

	/// @brief Returns the name of the class
	virtual std::string get_name() const;

public:
	// mutators
	void set_selector( core::pack::task::residue_selector::ResidueSelectorCOP selector_val );

public:
	void parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map );

	void read_db();

	std::list< std::string > const & allowed_aas(
		SurroundingSS const & surrounding,
		std::string const & loop_abego ) const;

	void set_allowed_aas(
		SurroundingSS const & surrounding,
		std::string const & loop_abego,
		std::list< std::string > const & allowed_aa );

	void disallow_aas(
		core::pack::task::PackerTask & task,
		LoopInfo const & info ) const;

	LoopInfoVec get_loop_info( core::pose::Pose const & pose ) const;

	LoopInfoVec loop_info_from_subset(
		core::pose::Pose const & pose,
		std::string const & ss,
		core::pack::task::residue_selector::ResidueSubset const & subset ) const;

	/// @brief if true, residues adjacent to loops will be restricted. Otherwise, just the loop. (default=false)
	void set_include_adjacent_residues( bool const include_res );

	void set_secstruct( std::string const & secstruct );

	void set_secstruct_from_blueprint( std::string const & bp_file );

private:
	std::string secstruct_;
	bool include_adjacent_residues_;
	core::pack::task::residue_selector::ResidueSelectorCOP selector_;
	ConsensusSequenceTable seqtable_;
};

} // task_operations
} // denovo_design
} // protocols

#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/farna/RNAThreadAndMinimizeMover.hh
/// @brief Thread a new sequence over a given RNA scaffold and do a little optimization
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_farna_RNAThreadAndMinimizeMover_hh
#define INCLUDED_protocols_farna_RNAThreadAndMinimizeMover_hh

// Unit headers
#include <protocols/farna/RNAThreadAndMinimizeMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#ifdef WIN32
#include <core/chemical/ResidueType.hh>
#else
#include <core/chemical/ResidueType.fwd.hh>
#endif
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/sequence/Sequence.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>

namespace protocols {
namespace farna {

///@brief Thread a new sequence over a given RNA scaffold and do a little optimization
class RNAThreadAndMinimizeMover : public protocols::moves::Mover {

public:

	RNAThreadAndMinimizeMover( std::string const & seq,
		std::string const & template_sequence_from_alignment,
		bool long_strategy = false,
		std::string const & input_sequence_type = "",
		std::string const & mutation_list = "",
		std::string const & insertion_list = ""
	);
	
	RNAThreadAndMinimizeMover():
		protocols::moves::Mover( RNAThreadAndMinimizeMover::class_name() ),
		seq_( "" ),
		template_sequence_from_alignment_( "" ),
		long_strategy_( false ),
		input_sequence_type_( "" ),
		scorefxn_( nullptr ),
		target_sequence_( "" ),
		working_res_( utility::vector1< Size >() ),
		mutation_list_( "" ),
		insertion_list_( "" )
	{}
	
	// copy constructor (not needed unless you need deep copies)
	//RNAThreadAndMinimize( RNAThreadAndMinimize const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~RNAThreadAndMinimizeMover() {};

	static std::string
	class_name();

public:

	void accomodate_length_change( Pose & pose, Size const insertion_begin, Size const seqpos );
	void process_insertions( core::pose::Pose & pose );
	void process_deletions( core::pose::Pose & pose );

	void set_up_target_sequence( core::pose::Pose & pose );	
	
	void
	obtain_rtypes_for_target_sequence();
	
	void
	long_mutate_strategy( core::pose::Pose & pose );
	
	void
	mutate_all_at_once( core::pose::Pose & pose );
	
	core::kinematics::MoveMapOP mm_from_residues( core::pose::Pose const & pose, utility::vector1< Size > const & changed_pos, bool add_nearby=false );
	
	////////////////////////////////////////////////////////////////////////////
	void
	setup_alignment_map( std::map< Size, Size > & mapping,
						std::string const & sequence_from_alignment ){
		Size count( 0 );
		for ( Size i = 1; i <= sequence_from_alignment.size(); i++ ) {
			if ( sequence_from_alignment[ i-1 ] != '-' ) {
				count++;
				mapping[ i ] = count;
			} else {
				mapping[ i ] = 0;
			}
		}
	}
	
	// mover virtual API
	virtual void
	apply( core::pose::Pose & pose );

	virtual void
	show( std::ostream & output = std::cout ) const;

	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//RNAThreadAndMinimizeMover & operator=( RNAThreadAndMinimizeMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

private:
	// The input sequence
	std::string seq_;
	
	// In prior versions, the sequence from an input alignment. May return
	// to that soon.
	std::string template_sequence_from_alignment_;

	// Use the slower method - mutations one at a time
	bool long_strategy_;
	
	// MODOMICS, IUPAC, or anything else = annotated
	std::string input_sequence_type_;
	
	// Scorefunction - set from create_score_function
	core::scoring::ScoreFunctionOP scorefxn_;
	
	// The target ANNOTATED sequence - set after parsing
	std::string  target_sequence_;
	
	// Manages mapping seqpos
	utility::vector1< Size > working_res_;
	
	// Target vector of residue types
	core::chemical::ResidueTypeCOPs rtypes_;
	
	// Set of mutations
	std::string mutation_list_;
	
	// Set of insertions
	std::string insertion_list_;
};

std::ostream &
operator<<( std::ostream & os, RNAThreadAndMinimizeMover const & mover );

} //protocols
} //farna

#endif //protocols/farna_RNAThreadAndMinimizeMover_hh

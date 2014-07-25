// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/MembraneDDGMover.cc
///
/// @brief      Compute ddG Scores for comparison with experimental ddGs in Membrane Proteins
/// @details	Initialize a membrane pose, compute an initial membrane position,
///				compute a native score, make mutation, repack sidechains within radius, and
///				score new structure. Uses a quick mutation object interface instead of resfiles
///				for easily breaking down tasks.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/20/14)

// Unit Headers
#include <protocols/membrane/ddG/MembraneDDGMover.hh>
#include <protocols/membrane/ddG/MembraneDDGMoverCreator.hh>

// Project Headers
#include <protocols/membrane/ddG/Mutation.hh>

#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/InitialMembranePositionMover.hh>

// Package headers
#include <protocols/moves/Mover.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/chemical/AA.hh>

#include <core/conformation/Conformation.hh> 
#include <core/conformation/Residue.hh>

#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/PointGraphData.hh>

#include <core/graph/UpperEdgeGraph.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh> 

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>
#include <set>

#include <iterator>
#include <algorithm>
#include <iostream>

static basic::Tracer TR( "protocols.membrane.ddG.MembraneDDGMover" );

namespace protocols {
namespace membrane {
namespace ddG {

using namespace core;
using namespace core::scoring;
using namespace protocols::moves;

////////////////////
/// Constructors ///
////////////////////

/// @brief Default Constructor (private)
/// @details Construct a defauilt version of this mover: pack_radius = 0.0, membrane
/// env smooth sfxn, and no mutations
MembraneDDGMover::MembraneDDGMover() :
	Mover(),
	pack_radius_( 0.0 )
{
	using namespace core::scoring;
	sfxn_ = ScoreFunctionFactory::create_score_function( "fa_menv_smooth_2014" );
}

/// @brief Custom Constructor
/// @details Create a ddG task given a pack radius, score funciton, and list of mutations
MembraneDDGMover::MembraneDDGMover(
	Real pack_radius_in,
	ScoreFunctionOP sfxn_in,
	utility::vector1< MutationOP > mutations_in
	) :
	Mover(),
	pack_radius_( pack_radius_in ),
	sfxn_( sfxn_in ),
	mutations_( mutations_in )
{}

/// @brief Copy Constructor
/// @details Create a deep copy of this object
MembraneDDGMover::MembraneDDGMover( MembraneDDGMover const & src ) :
	Mover( src ),
	pack_radius_( src.pack_radius_ ),
	sfxn_( src.sfxn_ )
{}

/// @brief Assignment Operator
/// @details Create a deep copy of this object overloading the assignment operator
MembraneDDGMover &
MembraneDDGMover::operator=( MembraneDDGMover const & src ) {
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new MembraneDDGMover( *this ) );
	
}

/// @brief Destructor
MembraneDDGMover::~MembraneDDGMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MembraneDDGMover::clone() const {
	return new MembraneDDGMover( *this );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MembraneDDGMover::fresh_instance() const {
	return new MembraneDDGMover();
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MembraneDDGMover::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
	) {
	// TODO: Implement this guy!
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
MembraneDDGMoverCreator::create_mover() const {
	return new MembraneDDGMover;
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MembraneDDGMoverCreator::keyname() const {
	return MembraneDDGMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MembraneDDGMoverCreator::mover_name() {
	return "MembraneDDGMover";
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (MembraneDDGMover)
std::string
MembraneDDGMover::get_name() const {
	return "MembraneDDGMover";
}

/// @brief Compute ddG Scores for Membranes
void
MembraneDDGMover::apply( Pose & pose ) {
	
	using namespace basic::options::OptionKeys;
	using namespace core::pack;
	using namespace protocols::simple_moves;
	using namespace protocols::membrane;
	
	// Add the membrane framework (membrane residue, topology, etc)
	// to the pose
	AddMembraneMoverOP add_memb = new AddMembraneMover();
	add_memb->apply( pose );
	
	// Compute the initial position of the membrane based on its xyz coords
	// and transmembrane spans
	InitialMembranePositionMoverOP init_position = new InitialMembranePositionMover();
	init_position->apply( pose );
	
	// Compute the score of the starting structure
	core::Real native_score = (*sfxn_)(pose);

	// Compute ddG score for each mutation and print result
	for ( Size i = 1; i <= mutations_.size(); ++i ) {
		
		core::Real ddG = compute_ddG_score(pose, mutations_[i]->aa(), mutations_[i]->position(), native_score );
		TR << mutations_[i]->aa() << " " << ddG << std::endl;
	}
	
}

/// @brief Compute ddG Score
/// @details Compute ddG from a copy of the native pose, native score
/// and provided resfile. Doing this PyRosetta Style - thanks to Evan Baugh's Mutate.py for
/// instructions
core::Real
MembraneDDGMover::compute_ddG_score(
	Pose & pose,
	AA aa,
	Size position,
	core::Real native_score
	) {

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pack::task;
	using namespace protocols::simple_moves;

	// Check the pose is still a fullatom pose
	if (! pose.is_fullatom() ) {
		utility_exit_with_message( "Cannot make mutation and repack rotamers on a non-fullatom pose" );
	}

	// Create a deep copy of the input pose
	PoseOP copy_pose = new Pose( pose );

	// Create a new packer task for making a mutation in the pose
	PackerTaskOP repack_task = TaskFactory::create_packer_task( *copy_pose );
	
	// Get residue neighbors within pack radius
	std::set< Size > neighbors = get_residue_neighbors( pose, position, pack_radius_ );
	
	// For residues within the pack radii, restrict to repacking (no design). For all others,
	// prevent repacking & design.
	for ( Size i = 1; i <= copy_pose->total_residue(); ++i ) {
		
		// If we are not looking at the desired mutated posiiton, only allow
		// repacking
		if ( i != position ) {
			repack_task->nonconst_residue_task( i ).restrict_to_repacking();
		}
		
		// If we ar enot looking at the desired mutated posiiton and not within pack
		// radius, prevent repacking all together
		if ( i != position && ( std::find( neighbors.begin(), neighbors.end(), i ) == neighbors.end() ) ) {
			repack_task->nonconst_residue_task( i ).prevent_repacking();
		}
	}
	
	
	// Restrict designed residues to the amino-acod of interest (aa => mutant_aa)
	SSize const canonical_aas = 20;
	utility::vector1< bool > allowed_residues;
	for ( SSize i = 1; i <= canonical_aas; ++i ) {
		allowed_residues.push_back( i == aa );
	}
	repack_task->nonconst_residue_task( position ).restrict_absent_canonical_aas( allowed_residues );
	
	// Setup a new pack rotamers mover and apply the packer task
	PackRotamersMoverOP packer = new PackRotamersMover( sfxn_, repack_task );
	packer->apply( *copy_pose );
		
	// Score the pose
	core::Real mutant_score = (*sfxn_)(*copy_pose);
				
	// Compute ddG
	core::Real ddG_score = mutant_score - native_score;
	
	// Return a new score!
	return ddG_score;
				
}

/// @brief Access Neighbors within Radius
/// @details Determine the number of neighbors within
/// the user provided radius
std::set< Size >
MembraneDDGMover::get_residue_neighbors( Pose & pose, core::Size position, core::Real radius ) {

	using namespace core::conformation;
	
	// Create a new residue point graph from the current pose conformaiton
	PointGraphOP pg( new PointGraph );
	residue_point_graph_from_conformation( pose.conformation(), *pg );
	find_neighbors<PointGraphVertexData,PointGraphEdgeData>( pg, radius ); //create edges
	
	std::set< core::Size > neighbors;
	// Iterate through the point graph and pick out residues
	for ( core::Size r(1); r <= position; ++r){
		for ( PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex(r).upper_edge_list_begin(),
			 edge_end_iter = pg->get_vertex(r).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			 
			  core::Size const other = edge_iter->upper_vertex();
			 
			if (other == position ) {
			
				TR << other << " " << position << " " << r << std::endl;
				neighbors.insert(r);
			}
			else if ( r == position ){
			TR << other << " " << position << " " << r << std::endl;
				neighbors.insert(other);
			}
		}
	}
	
	return neighbors;
}


/// @brief List of AAs to test
/// @details Return a vector 1 of characters representing the 20 AAs to substitute. This is
/// possibly too specific to my current task as well
utility::vector1< char >
MembraneDDGMover::designed_amino_acids() {

	// Listing the 20 canonical amino acids in a character array of size 20
	utility::vector1< char > aa;
	aa.push_back( 'A' );
	aa.push_back( 'C' );
	aa.push_back( 'D' );
	aa.push_back( 'E' );
	aa.push_back( 'F' );
	aa.push_back( 'G' );
	aa.push_back( 'H' );
	aa.push_back( 'I' );
	aa.push_back( 'K' );
	aa.push_back( 'L' );
	aa.push_back( 'M' );
	aa.push_back( 'N' );
	aa.push_back( 'P' );
	aa.push_back( 'Q' );
	aa.push_back( 'R' );
	aa.push_back( 'S' );
	aa.push_back( 'T' );
	aa.push_back( 'V' );
	aa.push_back( 'W' );
	aa.push_back( 'Y' );
	
	return aa;
		
}

} // ddG
} // membrane
} // protocols


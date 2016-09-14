// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/HelixFromSequence.cc
/// @brief Generates a (transmembrane) helix from a fasta file containing the sequence
/// @author Julia Koehler Leman (julia.koehler1982@gmail.com)

// Project Headers
#include <protocols/membrane/HelixFromSequence.hh>
#include <protocols/membrane/HelixFromSequenceCreator.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/Sequence.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/OptimizeProteinEmbeddingMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/relax/membrane/MPRangeRelaxMover.hh>
#include <protocols/relax/RangeRelaxMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>

// Package Headers
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/types.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <core/sequence/util.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>
#include <string>

static THREAD_LOCAL basic::Tracer TR( "protocols.membrane.HelixFromSequence" );

namespace protocols {
namespace membrane {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
HelixFromSequence::HelixFromSequence():
	protocols::moves::Mover( "HelixFromSequence" )
{
	register_options();
	set_defaults();
	init_from_cmd();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
HelixFromSequence::HelixFromSequence( HelixFromSequence const & )= default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
HelixFromSequence::~HelixFromSequence()= default;

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
HelixFromSequence::apply( core::pose::Pose & pose ){

	using namespace core::scoring;
	using namespace core::conformation::membrane;
	using namespace protocols::relax;
	using namespace protocols::relax::membrane;
	using namespace protocols::simple_moves;

	// create ideal helix pose from the sequence:
	// 1. create pose from sequence
	make_pose_from_sequence( pose, seq_, core::chemical::FA_STANDARD );
	TR << "pose:total_residue: " << pose.size() << std::endl;

	// 2. need to set the PDBInfo object in the pose, because
	// make_pose_from_sequence doesn't take care of that!
	PDBInfoOP pdbinfo( new PDBInfo( pose ) );
	pose.pdb_info( pdbinfo );

	// 3. create ideal helices from SSEs
	for ( Size i = 1; i <= pose.size(); ++i ) {
		pose.set_phi(   i, -62 );
		pose.set_psi(   i, -41 );
		pose.set_omega( i, 180 );
	}

	// if a membrane protein: transform helix into membrane:
	if ( mem_ == true ) {

		// create topology object from first to last residue
		SpanningTopologyOP topo( new SpanningTopology() );
		topo->add_span( 1, pose.size() );

		// in case of optimizing the membrane embedding
		if ( opt_mem_ == true ) {

			// run AddMembraneMover with topology
			AddMembraneMoverOP addmem( new AddMembraneMover( topo, 1, 0 ));
			addmem->apply( pose );

			topo->show();

			// transform into membrane and optimize embedding
			// runs TransformIntoMembrane underneath
			OptimizeProteinEmbeddingMoverOP opt( new OptimizeProteinEmbeddingMover() );
			opt->apply( pose );

			// run MPRangeRelax with default values (nres and 0.1 angle max)
			if ( ! skip_rlx_ ) {
				MPRangeRelaxMoverOP relax( new MPRangeRelaxMover() );
				relax->apply( pose );
			}

			// membrane protein but not optimizing embedding
		} else {

			// run AddMembraneMover with topology
			AddMembraneMoverOP addmem( new AddMembraneMover( topo, 1, 0 ));
			addmem->apply( pose );

			// transform into membrane
			TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover() );
			transform->apply( pose );

			// run MPRangeRelax with default values (nres and 0.1 angle max)
			if ( ! skip_rlx_ ) {
				MPRangeRelaxMoverOP relax( new MPRangeRelaxMover() );
				relax->apply( pose );
			}
		}
	} else {
		// if soluble protein
		// run RangeRelax with default values (nres and 0.1 angle max)
		if ( ! skip_rlx_ ) {
			RangeRelaxMoverOP relax( new RangeRelaxMover() );
			relax->apply( pose );
		}
	}

	// score pose for scorefile
	if ( mem_ ) {
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );
		ScoreMoverOP score( new ScoreMover( sfxn ) );
		score->apply( pose );
	} else {
		ScoreMoverOP score( new ScoreMover() );
		score->apply( pose );
	}

} // apply

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
HelixFromSequence::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

/// @brief Get the name of the Mover
std::string
HelixFromSequence::get_name() const {
	return "HelixFromSequence";
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
HelixFromSequence::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
HelixFromSequence::fresh_instance() const
{
	return protocols::moves::MoverOP( new HelixFromSequence );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
HelixFromSequence::clone() const{
	return protocols::moves::MoverOP( new HelixFromSequence( *this ) );
}

/*
HelixFromSequence & HelixFromSequenceoperator=( HelixFromSequence const & src){
return HelixFromSequence( src );
}
*/

std::ostream &operator<< (std::ostream &os, HelixFromSequence const &mover)
{
	mover.show(os);
	return os;
}

////////////////////////////////////////////////////////////////////////////////
/// Creator ///
///////////////

protocols::moves::MoverOP
HelixFromSequenceCreator::create_mover() const {
	return protocols::moves::MoverOP( new HelixFromSequence );
}

std::string
HelixFromSequenceCreator::keyname() const {
	return HelixFromSequenceCreator::mover_name();
}

std::string
HelixFromSequenceCreator::mover_name(){
	return "HelixFromSequence";
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

/// @brief register options
void HelixFromSequence::register_options(){

	using namespace basic::options;

	option.add_relevant( OptionKeys::in::file::fasta );
	option.add_relevant( OptionKeys::mp::setup::transform_into_membrane );
	option.add_relevant( OptionKeys::mp::transform::optimize_embedding );
	option.add_relevant( OptionKeys::relax::range::skip_relax );

} // register options

////////////////////////////////////////////////////////////////////////////////

/// @brief set defaults
void HelixFromSequence::set_defaults() {

	// fasta sequence
	seq_ = "";

	// membrane protein?
	mem_ = false;

	// optimize membrane position?
	opt_mem_ = false;

	// skip relax?
	skip_rlx_ = false;

} // set defaults

////////////////////////////////////////////////////////////////////////////////

/// @brief init from commandline
void HelixFromSequence::init_from_cmd() {

	using namespace basic::options;

	// read sequence from fasta
	if ( option[OptionKeys::in::file::fasta].user() ) {

		seq_ = core::sequence::read_fasta_file( option[ OptionKeys::in::file::fasta ]()[1] )[1]->sequence();
		TR << "Read in fasta file" << option[OptionKeys::in::file::fasta]()[1] << std::endl;
	} else {
		utility_exit_with_message("Please provide fasta file!");
	}

	// transform into membrane?
	if ( option[OptionKeys::mp::setup::transform_into_membrane].user() ) {
		mem_ = option[OptionKeys::mp::setup::transform_into_membrane]();

		if ( mem_ == true ) {
			TR << "Pose is a membrane protein and will be transformed into the membrane" <<  std::endl;
		}
	}

	// optimize membrane embedding?
	if ( option[OptionKeys::mp::transform::optimize_embedding].user() ) {
		opt_mem_ = option[OptionKeys::mp::transform::optimize_embedding]();

		if ( opt_mem_ == true ) {
			TR << "Protein embedding in the membrane will be optimized" <<  std::endl;
		}
	}

	// skip relax step after creating the helix?
	if ( option[OptionKeys::relax::range::skip_relax].user() ) {
		skip_rlx_ = option[OptionKeys::relax::range::skip_relax]();

		if ( skip_rlx_ == true ) {
			TR << "Skipping refinement step..." << std::endl;
		}
	}

} // init from cmd




} //protocols
} //membrane



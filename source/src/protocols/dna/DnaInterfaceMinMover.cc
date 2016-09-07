// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DnaInterfaceMinMover.cc
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/dna/DnaInterfaceMinMover.hh>
#include <protocols/dna/DnaInterfaceMinMoverCreator.hh>

// Package headers
#include <basic/datacache/DataMap.hh>
#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/filters/Filter.fwd.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <basic/options/option.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

//// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace dna {

using utility::vector1;
using namespace core;
using namespace kinematics;
using namespace optimization;
using namespace basic::options;
using namespace scoring;
using namespace moves;

using basic::t_warning;
using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static THREAD_LOCAL basic::Tracer TR( "protocols.dna.DnaInterfaceMinMover", t_info );

std::string
DnaInterfaceMinMoverCreator::keyname() const
{
	return DnaInterfaceMinMoverCreator::mover_name();
}

protocols::moves::MoverOP
DnaInterfaceMinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DnaInterfaceMinMover );
}

std::string
DnaInterfaceMinMoverCreator::mover_name()
{
	return "DnaInterfaceMinMover";
}

DnaInterfaceMinMover::DnaInterfaceMinMover()
: protocols::simple_moves::MinMover("DnaInterfaceMinMover"),
	interface_(/* 0 */),
	chi_(true),
	bb_(false)
{
	min_type( option[ OptionKeys::run::min_type ]() );
	tolerance( option[ OptionKeys::run::min_tolerance ]() );
}

DnaInterfaceMinMover::DnaInterfaceMinMover( DnaInterfaceMinMover const & src ) :
	//utility::pointer::ReferenceCount(),
	protocols::simple_moves::MinMover( src )
{
	*this = src;
}

DnaInterfaceMinMover &
DnaInterfaceMinMover::operator = ( DnaInterfaceMinMover const & src )
{
	if ( src.interface_ ) interface_ = src.interface_->clone();
	else interface_.reset();
	chi_ = src.chi_;
	bb_ = src.bb_;
	return *this;
}

DnaInterfaceMinMover::~DnaInterfaceMinMover()= default;

DnaInterfaceMinMover::DnaInterfaceMinMover( DnaInterfaceFinderOP interface )
: protocols::simple_moves::MinMover("DnaInterfaceMinMover"),
	interface_( interface ),
	chi_(true),
	bb_(false)
{
	runtime_assert( interface != nullptr );
	reset_from_interface();
	min_type( option[ OptionKeys::run::min_type ]() );
	tolerance( option[ OptionKeys::run::min_tolerance ]() );
}

void
DnaInterfaceMinMover::use_interface( DnaInterfaceFinderOP interface )
{
	runtime_assert( interface != nullptr );
	interface_ = interface;
	reset_from_interface();
}

void
DnaInterfaceMinMover::reset_from_interface()
{
	runtime_assert( interface_ != nullptr );
	MoveMapOP new_movemap( new MoveMap );
	// sidechain minimization for all amino acids in the vicinity of nucleotide bases
	// (according to 'interface', which need not represent the entire protein-DNA interface)
	for (const auto & itr : interface_->protein_neighbors()) {
		if ( itr.second.close() ) {
			// avoid adding pointless 'false' values to std::map if chi_ or bb_ are false
			if ( chi_ ) new_movemap->set_chi( itr.first, true );
			if ( bb_ ) new_movemap->set_bb( itr.first, true );
		}
	}
	movemap( new_movemap );
}

void
DnaInterfaceMinMover::apply( pose::Pose & pose )
{
	if ( ! interface_ ) {
		interface_ = DnaInterfaceFinderOP( new DnaInterfaceFinder );
		interface_->determine_protein_interface( pose );
		reset_from_interface();
	}
	if ( ! movemap() ) reset_from_interface(); // just in case
	if ( ! score_function() ) {
		score_function( get_score_function() );
	}

	AtomTreeMinimizer minimizer;
	// (*score_function())(pose_);
	minimizer.run( pose, *movemap(), *score_function(), *min_options() );
}

std::string
DnaInterfaceMinMover::get_name() const {
	return DnaInterfaceMinMoverCreator::mover_name();
}

/// @brief parse XML (in the context of the parser/scripting scheme)
void
DnaInterfaceMinMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
)
{
	if ( tag->getName() != "DnaInterfaceMinMover" ) {
		TR(t_warning) << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( tag->hasOption("scorefxn") ) {
		std::string const key( tag->getOption<std::string>("scorefxn") );
		if ( datamap.has( "scorefxns", key ) ) {
			score_function( datamap.get_ptr< ScoreFunction >( "scorefxns", key ) );
		}
	}
	if ( tag->hasOption("dna_interface_finder") ) {
		std::string const key( tag->getOption<std::string>("dna_interface_finder") );
		if ( datamap.has( "dna_interface_finder", key ) ) {
			use_interface( datamap.get_ptr< DnaInterfaceFinder >( "dna_interface_finder", key ) );
		}
	}
	if ( tag->hasOption("min_type") ) {
		min_type( tag->getOption<std::string>("min_type") );
	}
	if ( tag->hasOption("tolerance") ) {
		tolerance( tag->getOption<Real>("tolerance") );
	}
}

/// @brief required in the context of the parser/scripting scheme
MoverOP
DnaInterfaceMinMover::fresh_instance() const
{
	return MoverOP( new DnaInterfaceMinMover );
}

/// @brief required in the context of the parser/scripting scheme
MoverOP
DnaInterfaceMinMover::clone() const
{
	return MoverOP( new DnaInterfaceMinMover( *this ) );
}

} // namespace dna
} // namespace protocols

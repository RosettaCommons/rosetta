// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <protocols/abinitio/TemplateJumpSetup.hh>

// Package Headers
#include <protocols/abinitio/PairingStatistics.hh>
#include <protocols/abinitio/Template.hh>
#include <protocols/abinitio/Templates.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <protocols/jumping/JumpSample.hh>

//numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

// Utility headers
#include <utility>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/io/util.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

#include <utility/vector0.hh>


static basic::Tracer tr( "protocols.abinitio.Templates" );

namespace protocols {
namespace abinitio {


using namespace core;
using namespace fragment;
using namespace jumping;

TemplateJumpSetup::~TemplateJumpSetup() = default;
FixTemplateJumpSetup::~FixTemplateJumpSetup() = default;

TemplateJumpSetup::TemplateJumpSetup(
	TemplatesCOP templates,
	core::fragment::SecondaryStructureCOP secstruct,
	PairingStatisticsCOP strand_stats,
	core::scoring::dssp::PairingsList const& helix_pairings
) : templates_(std::move( templates )),
	secstruct_(std::move( secstruct )),
	strand_stats_(std::move( strand_stats )),
	helix_pairings_( helix_pairings )
{}

JumpSample
TemplateJumpSetup::create_jump_sample() const {
	if ( templates_ ) {
		tr.Debug << "create JumpSample from Templates: nres = " << templates_->target_total_residue() << "\n" << templates_->pairings() << std::endl;
	}
	runtime_assert( strand_stats_ != nullptr );

	if ( !strand_stats_->nr_models() ) {
		core::scoring::dssp::PairingList no_jumps; //create empty pairing list
		return JumpSample( secstruct_->total_residue(), no_jumps, *secstruct_ );
	}

	// pick a topology
	std::string topol_id;
	core::scoring::dssp::StrandPairingSet const& target_strands( strand_stats_->suggest_topology( topol_id ) );
	tr.Debug << "topology selected! " << target_strands << std::endl;
	if ( !target_strands.size() ) {
		core::scoring::dssp::PairingList no_jumps; //create empty pairing list
		return JumpSample( secstruct_->total_residue(), no_jumps, *secstruct_ );
	}

	// generate random sequence from 1 .. N
	utility::vector1< Size > strand_ids;
	Real total_weight(0.0); //avoids nan
	for ( Size i = 1; i <= target_strands.size(); i ++ ) {
		total_weight += strand_stats_->strand_weight( target_strands.strand_pairing( i ) );
		strand_ids.push_back( i );
	}

	numeric::random::random_permutation( strand_ids, numeric::random::rg() );
	numeric::random::random_permutation( strand_ids, numeric::random::rg() ); //want it really random
	tr.Debug << "total weight: " << total_weight << " size: " << target_strands.size() << "strand_ids.size(): "
		<< strand_ids.size()
		<< " permutated sequence: ";
	utility::io::write_vector( tr.Debug, strand_ids );
	tr.Debug << std::endl;

	// take first 2..nr_strands strand_ids
	Size const rg_nrstrands( std::min( static_cast< int >( numeric::random::rg().uniform() * (target_strands.size()) ) + 2 , (int) target_strands.size() ) );
	tr.Debug << "random choice:  aim for selection of " << rg_nrstrands << std::endl;
	core::scoring::dssp::PairingList target_pairings;
	Size nstrand_selected( 0 );
	for ( Size i = 1; i<= target_strands.size(); i++ ) {
		//decide if we pair this strand:
		Real weight( strand_stats_->strand_weight( target_strands.strand_pairing( strand_ids[ i ] ) ) );
		if ( total_weight > 0 && weight/total_weight < 0.005 ) continue;
		nstrand_selected++;
		tr.Debug << "pre-selected " << target_strands.strand_pairing( strand_ids[ i ] ) << std::endl;
		core::scoring::dssp::PairingList possible_pairings;
		target_strands.strand_pairing( strand_ids[ i ] ).get_beta_pairs( possible_pairings );
		tr.Debug << "pre-selected " << target_strands.strand_pairing( strand_ids[ i ] )
			<< "provides: " << possible_pairings.size() << " pairings"<< std::endl;
		// bias selection of individual pairing to center of strand
		// we do that now simply by adding the center pairings twice
		if ( possible_pairings.size() > 2 ) {
			Size const ori_size = possible_pairings.size();
			for ( Size i=2; i<=ori_size - 1; i++ ) {
				possible_pairings.push_back( possible_pairings[ i ] );
			}
		}
		Size const rg_pairing( static_cast< int > ( numeric::random::rg().uniform() * possible_pairings.size() ) + 1 );
		tr.Debug << "RANDOM: choose pairing " << rg_pairing << std::endl;
		//check if pairing is in templates...
		bool found( false );
		if ( templates_ ) {
			for ( auto it=templates_->begin(),
					eit = templates_->end(); it!=eit && !found; ++it ) {
				//check if at least 1 template has pairing .. we know all its pairings already aligned in targe-sequence...
				// its in the strand_stats_
				found = strand_stats_->topology( it->first ).has_pairing( possible_pairings[ rg_pairing ]);
			}
		}
		if ( found || !templates_ ) target_pairings.push_back( possible_pairings[ rg_pairing ]);
	}

	// go throu selected pairings and throw out stupid stuff...
	// 1.) don't want to have two pairings with same register in close vicinity...
	// 2.) don't want to apply very local pairings very often ... use only 20% of time
	core::scoring::dssp::PairingList final_selection;
	for ( auto & target_pairing : target_pairings ) {

		bool ignore( false );

		// is a very short sequence separation?
		if ( std::abs( (int) target_pairing.Pos2() - (int) target_pairing.Pos1() ) < 15 ) {
			tr.Debug << "pairing " << target_pairing << " is close in sequence..." << std::endl;
			if ( numeric::random::rg().uniform() < 0.2 ) {
				tr.Debug << "selected with 20% chance " << std::endl;
				final_selection.push_back( target_pairing );
			} else ignore = true;
		}

		// are already pairings selected with same register and not far away ?
		Size const my_reg = target_pairing.get_register();

		for ( auto fit = final_selection.begin(), efit = final_selection.end();
				fit!=efit && !ignore; ++fit ) {
			if ( my_reg == fit->get_register() ) {
				if ( std::abs( (int) target_pairing.Pos1() - (int) fit->Pos1() ) < 15
						|| std::abs( (int) target_pairing.Pos1() - (int) fit->Pos2() ) < 15
						|| std::abs( (int) target_pairing.Pos2() - (int) fit->Pos1() ) < 15
						|| std::abs( (int) target_pairing.Pos2() - (int) fit->Pos2() ) < 15 ) {
					ignore = true;
				}
			}
		}
		if ( !ignore ) final_selection.push_back( target_pairing );
		else tr.Debug << "pairing " << target_pairing << " ignored because it is redundant " << std::endl;
	}

	if ( helix_pairings_.size() ) {
		Size const nr_helix_jumps(
			std::min( static_cast< int >( (numeric::random::rg().uniform()+0.3) * ( helix_pairings_.size() ) + 1 ), (int) helix_pairings_.size() ) );
		core::scoring::dssp::PairingList perm = helix_pairings_;
		numeric::random::random_permutation( perm , numeric::random::rg() );
		for ( Size i = 1; i<= nr_helix_jumps; i++ ) {
			final_selection.push_back( helix_pairings_[ i ] );
		}
	}
	tr.Debug << "selected jumps: " << final_selection << std::endl;
	// TODO in distant future: take care of bulges if there are any:
	return JumpSample( secstruct_->total_residue(), final_selection, *secstruct_ );
}

jumping::JumpSample
TemplateJumpSetup::clean_jumps( JumpSample const& target_jumps ) const {
	runtime_assert( strand_stats_ != nullptr );

	if ( !strand_stats_->nr_models() ) {
		core::scoring::dssp::PairingList no_jumps; //create empty pairing list
		return jumping::JumpSample( secstruct_->total_residue(), no_jumps, *secstruct_ );
	}

	core::scoring::dssp::PairingList filtered_jumps;

	core::fragment::FrameList jump_frames;
	kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_jump( true );

	target_jumps.generate_jump_frames( jump_frames,mm );

	// treat each jump individually ---> want to select frags only from templates that have compatible pairings
	for ( auto & jump_frame : jump_frames ) {
		core::scoring::dssp::Pairing target_pairing( target_jumps.get_pairing( jump_frame->start(), jump_frame->stop() ) );
		bool found( false );
		for ( PairingStatistics::const_iterator it = strand_stats_->begin(); !found && it != strand_stats_->end(); ++it ) {
			found =  it->second.pairing().has_pairing( target_pairing );
		}
		if ( found ) {
			filtered_jumps.push_back( target_pairing );//add
		}
	}
	return jumping::JumpSample( secstruct_->total_residue(), filtered_jumps, *secstruct_ );
}

bool TemplateJumpSetup::is_helix_jump( core::scoring::dssp::Pairing const& p ) const {
	auto iter =  find( helix_pairings_.begin(), helix_pairings_.end(), p );
	if ( iter != helix_pairings_.end() ) return true;
	core::scoring::dssp::Pairing rev = p;
	rev.reverse();
	iter =  find( helix_pairings_.begin(), helix_pairings_.end(), rev );
	return ( iter != helix_pairings_.end() );
}

/// @brief returns an ordered FragSet that is compatible with the JumpSample
/// default: generate jumps from ss-library according to JumpSample
core::fragment::FragSetOP
TemplateJumpSetup::generate_jump_frags( JumpSample const& target_jumps, kinematics::MoveMap const& mm ) const {
	tr.Debug <<" generate jump fragments... " << std::endl;
	FragSetOP jump_frags( new OrderedFragSet );
	//generate fragments from all homologes
	core::fragment::FrameList jump_frames;
	core::scoring::dssp::PairingsList library_pairings;

	target_jumps.generate_jump_frames( jump_frames,mm );

	// treat each jump individually ---> want to select frags only from templates that have compatible pairings
	for ( auto & jump_frame : jump_frames ) {
		core::scoring::dssp::Pairing target_pairing( target_jumps.get_pairing( jump_frame->start(), jump_frame->stop() ) );
		Size nr_frags( 0 );
		tr.Debug << "get frags for pairing " << target_pairing << std::endl;
		if ( templates_ ) {
			if ( !is_helix_jump( target_pairing ) ) {
				for ( auto const & it : *templates_ ) {

					//check if template has pairing .. we know all its pairings already aligned in targe-sequence...
					// its in the strand_stats_
					std::string const& model_id( it.first );
					if ( strand_stats_->topology( model_id ).has_pairing( target_pairing ) ) {
						nr_frags++;
						FrameList aFrame; aFrame.push_back( jump_frame );
						it.second->steal_frags( aFrame, *jump_frags  );
					}
				}
			} else { //get here if pairing is helix jump
				tr.Debug << "has been found in helix-list blindly collect all jump-geometries from models with an H" << std::endl;
				FrameList aFrame; aFrame.push_back( jump_frame );
				for ( auto const & it : templates_->helixjump_picks() ) {
					nr_frags++;
					it->steal_frags( aFrame, *jump_frags );
				}
			}
		}
		if ( nr_frags == 0 && templates_ ) {
			utility_exit_with_message("selected a pairing for which no pairing could be found in templates... shouldn't ever happen");
		}
		if ( nr_frags == 0 ) {
			//use standard ss-library
			library_pairings.push_back( target_pairing );
		}
	}

	jumping::StandardPairingLibrary::get_instance()->generate_jump_frags(
		library_pairings,
		mm,
		true /*with Torsion*/,
		*jump_frags
	);

	return jump_frags; //it is better to generate each time, since JumpingFoldConstraints might add other fragments ( steal native etc )
}

FixTemplateJumpSetup::FixTemplateJumpSetup(
	TemplatesCOP templates,
	core::fragment::SecondaryStructureCOP secstruct,
	PairingStatisticsCOP strand_stats,
	core::scoring::dssp::PairingsList const& helix_pairings,
	BaseJumpSetupOP jump_def
) : TemplateJumpSetup(  templates , secstruct, strand_stats, helix_pairings ),
	jump_def_(std::move( jump_def ))
{}

FixTemplateJumpSetup::FixTemplateJumpSetup(
	TemplateJumpSetup const& templ,
	BaseJumpSetupOP jump_def
) : TemplateJumpSetup( templ ),
	jump_def_(std::move( jump_def ))
{}

JumpSample
FixTemplateJumpSetup::create_jump_sample() const {
	return jump_def_->create_jump_sample();
}

} //abinitio
} //protocols

#if 0
// take first 2..nr_strands strand_ids
Size const rg_nrstrands( std::min( static_cast< int >( numeric::random::rg().uniform() * (target_strands.size()) ) + 2 , (int) target_strands.size() ) );
tr.Debug << "random choice:  aim for selection of " << rg_nrstrands << std::endl;
core::scoring::dssp::PairingList target_pairings;
Size nstrand_selected( 0 );
for ( Size i = 1; i<= target_strands.size(); i++ ) {
	//decide if we pair this strand:
	Real weight( strand_stats_->strand_weight( target_strands.strand_pairing( strand_ids[ i ] ) ) );
	if ( total_weight > 0 && weight*rg_nrstrands / total_weight < numeric::random::rg().uniform() ) continue;
	nstrand_selected++;
	tr.Debug << "pre-selected " << target_strands.strand_pairing( strand_ids[ i ] ) << std::endl;
	core::scoring::dssp::PairingList possible_pairings;
	target_strands.strand_pairing( strand_ids[ i ] ).get_beta_pairs( possible_pairings );
	tr.Debug << "pre-selected " << target_strands.strand_pairing( strand_ids[ i ] )
		<< "provides: " << possible_pairings.size() << " pairings"<< std::endl;
	// bias selection of individual pairing to center of strand
	// we do that now simply by adding the center pairings twice
	if ( possible_pairings.size() > 2 ) {
		Size const ori_size = possible_pairings.size();
		for ( Size i=2; i<=ori_size - 1; i++ ) {
			possible_pairings.push_back( possible_pairings[ i ] );
		}
	}
	Size const rg_pairing( static_cast< int > ( numeric::random::rg().uniform() * possible_pairings.size() ) + 1 );
	tr.Debug << "RANDOM: choose pairing " << rg_pairing << std::endl;
	//check if pairing is in templates...
	bool found( false );
	if ( templates_ ) {
		for ( Templates::const_iterator it=templates_->begin(),
				eit = templates_->end(); it!=eit && !found; ++it ) {
			//check if at least 1 template has pairing .. we know all its pairings already aligned in targe-sequence...
			// its in the strand_stats_
			found = strand_stats_->topology( it->first ).has_pairing( possible_pairings[ rg_pairing ]);
		}
	}

	if ( found || !templates_ ) target_pairings.push_back( possible_pairings[ rg_pairing ]);
}

// go throu selected pairings and throw out stupid stuff...
// 1.) don't want to have two pairings with same register in close vicinity...
// 2.) don't want to apply very local pairings very often ... use only 20% of time
#endif

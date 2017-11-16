// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/AbscriptMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/AbscriptMover.hh>

// Package headers
#include <protocols/abinitio/abscript/AbscriptMoverCreator.hh>
#include <protocols/abinitio/abscript/FragmentCM.hh>
#include <protocols/abinitio/abscript/FragmentJumpCM.hh>
#include <protocols/abinitio/abscript/AbscriptStageMover.hh>
#include <protocols/abinitio/abscript/StagePreparer.hh>

#include <protocols/environment/claims/EnvClaim.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>

#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/SecondaryStructure.hh>

#include <core/pose/Pose.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <protocols/moves/MonteCarlo.hh>

#include <protocols/simple_moves/SmoothFragmentMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/GunnCost.hh>

#include <protocols/canonical_sampling/mc_convergence_checks/util.hh>

#include <protocols/jd2/util.hh>

#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/JumpSample.hh>

#include <protocols/checkpoint/CheckPointer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>

#include <utility/excn/Exceptions.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

//option includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/fold_cst.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

#ifdef WIN32
#include <basic/datacache/WriteableCacheableMap.hh>
#endif

static basic::Tracer tr( "protocols.abinitio.abscript.AbscriptMover", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace protocols::environment;
using namespace protocols::environment::claims;

// If the temperature needs to be altered, the obvious way to do it is through
// a command line option, which could/will replace this. For some reason, though
// in the FragmentSampler it's not set that way...
core::Real const DEFAULT_TEMP = 2.0;

// I would really like to define this in StageID.hh, but that gives a linker error
// and I'm not sure how to fix it without adding a .cc file specifically for this.

StageID& increment_stageid( StageID& id ) {
	id = static_cast< StageID >( static_cast<int>( id ) + 1 );
	debug_assert( !( id > END ) );
	return id;
}


class AbscriptMover::StageTracker {

public:
	StageTracker( moves::MonteCarloOP mc );

	void begin_stage( std::string const& stagename, core::pose::Pose& );

	void do_step( core::pose::Pose&, AbscriptStageMover&, core::Real const& progress );

	void end_stage( core::pose::Pose& pose );

private:
	clock_t starttime_;
	std::string stagename_;
	moves::MonteCarloOP mc_;
	core::Size iterations_;
	checkpoint::CheckPointer checkpointer_;
};


// creator
// XRW TEMP std::string
// XRW TEMP AbscriptMoverCreator::keyname() const {
// XRW TEMP  return AbscriptMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AbscriptMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AbscriptMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AbscriptMover::mover_name() {
// XRW TEMP  return "AbscriptMover";
// XRW TEMP }

core::scoring::ScoreFunctionOP setup_score( std::string const& scorename,
	basic::options::FileVectorOptionKey optkey ) {
	using namespace core::scoring;
	using namespace basic::options;

	core::scoring::ScoreFunctionOP scorefn;

	if ( option[ optkey ].user() ) {
		scorefn = ScoreFunctionFactory::create_score_function( scorename, option[ optkey ]() );
	} else {
		scorefn = ScoreFunctionFactory::create_score_function( scorename );
	}

	scorefn->set_weight( core::scoring::atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ] );

	debug_assert( scorefn );
	return scorefn;
}

AbscriptMover::AbscriptMover():
	ClientMover(),
	id_map_(),
	stage_movers_(),
	mc_()
{
	using namespace basic::options::OptionKeys;
	using namespace basic::options;

	//I don't love this. This should be const static data.
	id_map_["I"] = I; id_map_["II"] = II;
	id_map_["IIIa"] = IIIa; id_map_["IIIb"] = IIIb;
	id_map_["IVa"] = IVa; id_map_["IVb"] = IVb;

	mc_ = moves::MonteCarloOP( new moves::MonteCarlo( *setup_score( "score3", OptionKeys::abinitio::stage4_patch ), DEFAULT_TEMP ) );
	canonical_sampling::mc_convergence_checks::setup_convergence_checks_from_cmdline( *mc_ );

	core::scoring::ScoreFunctionOP tmp_score = setup_score( "score0", OptionKeys::abinitio::stage1_patch );
	tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage1 ]() );
	stage_movers_[ I ] = AbscriptStageMoverOP( new AbscriptStageMover( I, mc_, tmp_score, 2000 ) );

	tmp_score = setup_score( "score1", OptionKeys::abinitio::stage2_patch );
	tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage2 ]() );
	stage_movers_[ II ] = AbscriptStageMoverOP( new AbscriptStageMover( II, mc_, tmp_score, 2000 ) );

	tmp_score = setup_score( "score2", OptionKeys::abinitio::stage3a_patch );
	tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage3 ]() );
	stage_movers_[ IIIa ] = AbscriptStageMoverOP( new AbscriptStageMover( IIIa, mc_, tmp_score, 2000 ) );

	tmp_score = setup_score( "score5", OptionKeys::abinitio::stage3b_patch );
	tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage3 ]() );
	stage_movers_[ IIIb ] = AbscriptStageMoverOP( new AbscriptStageMover( IIIb, mc_, tmp_score, 2000 ) );

	tmp_score = setup_score( "score3", OptionKeys::abinitio::stage4_patch );
	tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage4 ]() );
	stage_movers_[ IVa ] = AbscriptStageMoverOP( new AbscriptStageMover( IVa, mc_, tmp_score, 4000 ) );

	tmp_score = setup_score( "score3", OptionKeys::abinitio::stage4_patch );
	tmp_score->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage4 ]() );
	stage_movers_[ IVb ] = AbscriptStageMoverOP( new AbscriptStageMover( IVb, mc_, tmp_score, 4000 ) );

	core::Real seqsep_stage1 = 0.15;
	core::Real seqsep_stage3 = 1.0; //have basically all but the furthest chainbreaks active at end of stage3.
	core::Real seqsep_stage4 = 1.2; //make sure that all chainbreaks are switched on in the end

	if ( option[ OptionKeys::fold_cst::seq_sep_stages ].user() ) {
		if ( option[ OptionKeys::fold_cst::seq_sep_stages ]().size() != 3 ) {
			throw utility::excn::EXCN_BadInput("Option seq_sep_stages requires exactly 3 values.");
		}
		seqsep_stage1 = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 1 ]; //default 15
		seqsep_stage3 = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 2 ]; //default 15
		seqsep_stage4 = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 3 ]; //default 15
	}

	stage_movers_[ I    ]->set_seq_sep_ramping( 0.0, seqsep_stage1 );
	stage_movers_[ II   ]->set_seq_sep_ramping( 0.0, seqsep_stage1 );
	stage_movers_[ IIIa ]->set_seq_sep_ramping( seqsep_stage3 - seqsep_stage1 , seqsep_stage1 );
	stage_movers_[ IIIb ]->set_seq_sep_ramping( seqsep_stage3 - seqsep_stage1 , seqsep_stage1 );
	stage_movers_[ IVa  ]->set_seq_sep_ramping( seqsep_stage4 - seqsep_stage3 , seqsep_stage3 );
	stage_movers_[ IVb  ]->set_seq_sep_ramping( seqsep_stage4 - seqsep_stage3 , seqsep_stage3 );

	stage_movers_[ I    ]->set_chainbreak_ramping( 0.0       , 0.00       );
	stage_movers_[ II   ]->set_chainbreak_ramping( 0.0       , 0.25 / 3.0 );
	stage_movers_[ IIIa ]->set_chainbreak_ramping( 2.5 / 3.0 , 0.00       );
	stage_movers_[ IIIb ]->set_chainbreak_ramping( 0.5 / 3.0 , 0.00       );
	stage_movers_[ IVa  ]->set_chainbreak_ramping( 1.5 / 3.0 , 2.50 / 3.0 );
	stage_movers_[ IVb  ]->set_chainbreak_ramping( 1.5 / 3.0 , 2.50 / 3.0 );
}

AbscriptMover::AbscriptMover( AbscriptMover const& src ):
	ClientMover( src ),
	id_map_( src.id_map_ ),
	stage_movers_( src.stage_movers_ ),
	mc_( src.mc_ )
{}

moves::MoverOP AbscriptMover::fresh_instance() const {
	return moves::MoverOP( new AbscriptMover() );
}

moves::MoverOP AbscriptMover::clone() const {
	return moves::MoverOP( new AbscriptMover( *this ) );
}

void AbscriptMover::apply( core::pose::Pose& pose ){

	using namespace basic::options;

	if ( ! pose.is_centroid() ) {
		throw utility::excn::EXCN_BadInput( "AbscriptMover recieved a non-centroid pose. Only centroid poses are accepted because of the score functions used in this algorithm." );
	}
	StageTracker tracker( mc_ );

	std::map< core::Size, core::Size > ITERATIONS = calculate_iterations( pose );

	// Stage I ---------------------------------------------------------------------
	mc_->clear_poses();
	mc_->reset( pose );

	tracker.begin_stage( "I", pose );
	for ( core::Size i = 1; i <= ITERATIONS[1]; ++i ) {
		core::Real progress = core::Real( i - 1 )/ITERATIONS[1];
		tracker.do_step( pose,
			*stage_movers_[I],
			progress );
	}
	tracker.end_stage( pose );

	// Stage II --------------------------------------------------------------------
	tracker.begin_stage( "II", pose );
	tracker.do_step( pose, *stage_movers_[II], 0.0 );
	tracker.end_stage( pose );

	// Stage III -------------------------------------------------------------------
	tracker.begin_stage( "III", pose );
	for ( core::Size i = 1; i <= ITERATIONS[3]; ++i ) {
		// sequence is IIIb, IIIa, IIIb, IIIa, IIIb, IIIa, IIIb, IIIa, IIIa, IIIa
		core::Real progress = core::Real( i - 1 )/ITERATIONS[3];
		tracker.do_step( pose,
			*stage_movers_[ ( ( i % 2 ) == 0 || ( i == 9 ) ) ? IIIa : IIIb ],
			progress );
	}
	tracker.end_stage( pose );

	// Stage IV ---------------------------------------------------------------------
	tracker.begin_stage( "IV", pose );
	for ( core::Size i = 1; i <= ITERATIONS[4]; ++i ) {
		// sequence is IVa, IVb, IVb
		core::Real progress = core::Real( i - 1 )/ITERATIONS[4];
		tracker.do_step( pose,
			*stage_movers_[ ( i == 1 ) ? IVa : IVb ],
			progress );
	}
	tracker.end_stage( pose );

	mc_->recover_low( pose );
}

void AbscriptMover::yield_submovers( MoverSet& movers_out ) const {
	for ( std::map< StageID, AbscriptStageMoverOP >::const_iterator it = stage_movers_.begin();
			it != stage_movers_.end(); ++it ) {
		it->second->yield_submovers( movers_out );
	}
}

void
AbscriptMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & ) {

	using namespace basic::options;

	//Add abscript-configured scorefunctions to DataMap?
	utility::vector1< std::string > id_strs( (Size) IVb );
	id_strs[ I ] = "I"; id_strs[ II ] = "II"; id_strs[ IIIa ] = "IIIa";
	id_strs[ IIIb ] = "IIIb"; id_strs[IVa] = "IVa"; id_strs[ IVb ] = "IVb";

	for ( StageID id = I; id <= IVb; increment_stageid(id) ) {
		std::string const scorefxn_name = tag->getOption<std::string>( "name" )+"_stage"+id_strs[ id ];

		if ( datamap.has( "scorefxns", scorefxn_name ) ) {
			//This error condition would be pretty easy to remove by just changing the way
			throw utility::excn::EXCN_BadInput( "Automatic score function " + scorefxn_name +
				" was already registered." );
		}

		datamap.add( "scorefxns", scorefxn_name,
			stage_movers_[id]->scorefxn()->clone() );
	}

	//cycles control default
	core::Real cycles_prefactor = 1.0;

	if ( tag->hasOption( "cycles" ) ) {
		cycles_prefactor = tag->getOption< core::Real >( "cycles", 1.0 );
		if ( option[ OptionKeys::abinitio::increase_cycles ].user() ) {
			throw utility::excn::EXCN_RosettaScriptsOption("Flag increase_cycles and AbscriptMover 'cycles' option are incompatible.");
		}
	} else if ( option[ OptionKeys::abinitio::increase_cycles ].user() ) {
		cycles_prefactor = option[ OptionKeys::abinitio::increase_cycles ];
	}

	if ( option[ OptionKeys::run::test_cycles ].user() ) {
		// This has to go here rather than in the above if block so that we don't get die for unaccessed tags problems.
		tr.Info << "Found -run:test_cycles; overriding cycles behavior to test mode." << std::endl;
		cycles_prefactor = 0.01;
	}

	for ( std::map< StageID, AbscriptStageMoverOP >::const_iterator it = stage_movers_.begin();
			it != stage_movers_.end(); ++it ) {
		it->second->set_cycles_adjust( cycles_prefactor );
	}

	utility::vector1< StageID > skipped_stages;
	if ( tag->hasOption("skip_stages") ) {
		if ( option[ OptionKeys::abinitio::skip_stages ].user() ) {
			throw utility::excn::EXCN_RosettaScriptsOption("Flag abinitio::skip_stages and AbscriptMover 'skip_stages' option are incompatible.");
		} else {
			utility::vector1< std::string > const skipped_str( utility::string_split( tag->getOption< std::string >( "skip_stages" ), ',' ) );
			for ( utility::vector1< std::string >::const_iterator stage_it = skipped_str.begin();
					stage_it != skipped_str.end(); ++stage_it ) {
				if ( id_map_.find( *stage_it ) != id_map_.end() ) {
					skipped_stages.push_back( id_map_[*stage_it] );
				} else {
					Size stage_int = 0;
					try {
						stage_int = utility::string2int( *stage_it );
					} catch (...) {
						throw utility::excn::EXCN_RosettaScriptsOption( "Stage '"+*stage_it+"' is not recognized." );
					}
					if ( stage_int > 0 && stage_int <= 4 ) {
						skipped_stages.push_back( StageID( stage_int ) );
					} else {
						tr.Error << "Skipped stage " << stage_int << " must be between 1 and 4." << std::endl;
						utility::excn::EXCN_RosettaScriptsOption( "Stage number "+*stage_it+" is invalid." );
					}
				}
			}
		}
	} else if ( option[ OptionKeys::abinitio::skip_stages ].user() ) {
		for ( IntegerVectorOption::const_iterator it = option[ OptionKeys::abinitio::skip_stages ]().begin();
				it != option[ OptionKeys::abinitio::skip_stages ]().end(); ++it ) {
			if ( *it < 1 || *it > 4 ) {
				throw utility::excn::EXCN_BadInput( "The option abinitio::skip_stages specified value "+
					utility::to_string( *it )+", which is not a valid stage." );
			}
			skipped_stages.push_back( StageID(*it) );
		}
	}

	for ( utility::vector1< StageID >::const_iterator stage = skipped_stages.begin();
			stage != skipped_stages.end(); ++stage ) {
		stage_movers_[ *stage ]->set_cycles_adjust( 0.0 );
		tr.Info << "AbscriptMover found abintio:skip_stages for stage " << *stage
			<< ". This stage will not be run." << std::endl;
	}

	//Load submovers using subtags.
	typedef utility::vector0< utility::tag::TagCOP > TagVector;
	TagVector const& subtags = tag->getTags();

	for ( TagVector::const_iterator subtag_it = subtags.begin();
			subtag_it != subtags.end(); ++subtag_it ) {
		TagCOP stagetag = *subtag_it;

		if ( stagetag->getName() == "Stage" ) {
			std::string const& id_str = stagetag->getOption<std::string>("ids");
			tr.Debug << "Adding movers for stages: " << id_str << std::endl;
			StageIDs stages = parse_stage_id( id_str );

			TagVector const& movertags = stagetag->getTags();

			for ( TagVector::const_iterator movertag_it = movertags.begin();
					movertag_it != movertags.end(); ++movertag_it ) {
				TagCOP movertag = *movertag_it;

				if ( movertag->getName() == "Mover" ||
						movertag->getName() == "Preparer" ) {
					std::string name = movertag->getOption<std::string>( "name", "null" );
					if ( movers.find( name ) == movers.end() ) {
						tr.Error << "Mover not found for XML tag:\n" << movertag << std::endl;
						throw utility::excn::EXCN_RosettaScriptsOption("Mover '"+name+"' not found.");
					}

					if ( movertag->getName() == "Mover" ) {
						core::Real weight = movertag->getOption< core::Real >( "weight", 1.0 );
						register_submover( movers.find( name )->second, stages, weight );
					} else if ( movertag->getName() == "Preparer" ) {
						register_preparer( movers.find( name )->second, stages );
					} else {
						throw utility::excn::EXCN_RosettaScriptsOption( "Stage subtag "+movertag->getName()+
							" is not valid. Only 'Mover' and 'Preparer' are accepted." );
					}
				}
			}
		} else if ( stagetag->getName() == "Fragments" ||
				stagetag->getName() == "fragments" ) {
			core::select::residue_selector::ResidueSelectorCOP selector( NULL );
			if ( stagetag->hasOption( "selector" ) ) {
				selector = datamap.get_ptr<core::select::residue_selector::ResidueSelector const>( "ResidueSelector", stagetag->getOption<std::string>( "selector" ) );
			}

			add_frags( stagetag->getOption< std::string >( "small_frags" ),
				stagetag->getOption< std::string >( "large_frags" ),
				selector,
				stagetag->getOption< bool >( "initialize", true ));
		} else {
			tr.Error << "AbscriptMover recieved illegal tag ('Stage' and 'Fragments' are acceptable): '" << stagetag << "'" << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption("Illegal AbscriptMover subtag.");
		}
	}

	if ( tag->hasOption( "std_frags" ) ) {
		// The die for unread tags will cause a crash, but we want to elaborate a bit on why.
		tr.Error << "The AbscriptMover tag 'std_frags' is deprecated. Use 'Fragments' instead. For example:\n"
			<< "                 <Fragments small_frags='frag3.txt' large_frags='frag9.txt' selector='ChainA' />"
			<< std::endl;
	}

	if ( !tag->hasTag( "Fragments" ) && !tag->hasTag( "fragments" ) ) {
		if ( !tag->hasTag( "no_fragments" ) &&
				option[ OptionKeys::in::file::frag3 ].active() &&
				option[ OptionKeys::in::file::frag9 ].active() ) {
			add_frags( tag->getOption<std::string>("small_frags", option[ OptionKeys::in::file::frag3 ]() ),
				tag->getOption<std::string>("large_frags", option[ OptionKeys::in::file::frag9 ]() ) );
		} else {
			tr.Warning << this->get_name() << " being defined without any fragments." << std::endl;
		}
	}
}

void AbscriptMover::add_frags(
	std::string const& small_fragfile,
	std::string const& large_fragfile,
	core::select::residue_selector::ResidueSelectorCOP selector,
	bool initialize ){

	using namespace core::fragment;
	using namespace basic::options;
	using namespace protocols::simple_moves;

	FragmentIO small_frag_io( option[ OptionKeys::abinitio::number_3mer_frags ](), 1,
		option[ OptionKeys::frags::annotate ]() );

	FragmentIO big_frag_io( option[ OptionKeys::abinitio::number_9mer_frags ](), 1,
		option[ OptionKeys::frags::annotate ]() );

	FragSetOP frags_small = small_frag_io.read_data( small_fragfile );
	FragSetOP frags_large = big_frag_io.read_data( large_fragfile );

	//embed standard movers in movers that will claim and unlock for them.
	FragmentCMOP claim_large(
		new FragmentCM( simple_moves::FragmentMoverOP( new ClassicFragmentMover( frags_large ) ) ) );
	FragmentCMOP claim_small(
		new FragmentCM( simple_moves::FragmentMoverOP( new ClassicFragmentMover( frags_small ) ) ) );
	FragmentCMOP claim_smooth(
		new FragmentCM( simple_moves::FragmentMoverOP( new SmoothFragmentMover( frags_small, simple_moves::FragmentCostOP( new GunnCost() ) ) ) ) );


	// Initialize controls if fragments are inserted along the
	// length of the protein as the first move.
	claim_large->initialize(initialize);
	claim_small->initialize(false);
	claim_smooth->initialize(false);

	if ( selector ) {
		claim_large->set_selector( selector );
		claim_small->set_selector( selector );
		claim_smooth->set_selector( selector );
	}

	claim_small->yield_cut_bias( true );

	//apply to appropriate stages
	for ( StageID id = I; id <= IIIb; increment_stageid(id) ) {
		stage_movers_[id]->add_submover( claim_large, 1.0 );
	}

	stage_movers_[ IVa ]->add_submover( claim_small, 1.0 );
	stage_movers_[ IVb ]->add_submover( claim_smooth, 1.0 );
}

// XRW TEMP std::string AbscriptMover::get_name() const {
// XRW TEMP  return "AbscriptMover";
// XRW TEMP }

void verify_stage_ID( std::map< std::string, StageID > const& id_map, std::string const& id ){
	if ( id_map.find( id ) == id_map.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Stage id '" + id + "' is not recognized." );
	}
}

StageIDs AbscriptMover::parse_stage_id( std::string const& id_str ) const {
	typedef utility::vector0< std::string > Strings;
	Strings comma_delin = utility::string_split( id_str, ',' );

	StageIDs ids;

	for ( Strings::const_iterator range = comma_delin.begin();
			range != comma_delin.end(); ++range ) {
		Strings hyphen_delin = utility::string_split( *range, '-' );
		if ( hyphen_delin.size() == 1 ) {
			verify_stage_ID( id_map_, hyphen_delin[0] );
			ids.push_back( id_map_.find( hyphen_delin[0] )->second );
		} else if ( hyphen_delin.size() == 2 ) {
			verify_stage_ID( id_map_, hyphen_delin[0] );
			verify_stage_ID( id_map_, hyphen_delin[1] );

			for ( StageID id = id_map_.find( hyphen_delin[0] )->second;
					id <= id_map_.find( hyphen_delin[1] )->second; increment_stageid(id) ) {
				ids.push_back( id );
			}
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("Stage ID range syntax error for '" + *range + "'" );
		}
	}

	return ids;
}

void AbscriptMover::register_submover( protocols::moves::MoverOP mover_in,
	StageIDs const& ids,
	core::Real weight ){
	ClientMoverOP mover = utility::pointer::dynamic_pointer_cast< protocols::environment::ClientMover > ( mover_in );
	if ( !mover ) {
		tr.Error << "The mover " << mover_in->get_name()
			<< " cannot be attached to the Abscript mover because it is not a ClientMover "
			<< std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("Abscript mover takes only ClientMovers.");
	}

	for ( StageIDs::const_iterator id = ids.begin(); id != ids.end(); ++id ) {
		stage_movers_[ *id ]->add_submover( mover, weight );
	}
}

std::map< core::Size, core::Size > AbscriptMover::calculate_iterations( core::pose::Pose const& pose ){
	std::map< core::Size, core::Size > ITERATIONS;

	//Stage1 has FT-dependent seqsep ramping ---------------------------------------
	core::kinematics::ShortestPathInFoldTree shortestpath( pose.fold_tree() );
	debug_assert( shortestpath.max_dist() <= pose.size() );

	core::Size end_seqsep = (core::Size)(stage_movers_[ I ]->seq_sep_intercept() * shortestpath.max_dist());
	core::Size begin_seqsep = std::min( end_seqsep, Size(3) );

	core::Real end_seqsep_factor   = core::Real( end_seqsep   ) / shortestpath.max_dist();
	core::Real begin_seqsep_factor = core::Real( begin_seqsep ) / shortestpath.max_dist();

	stage_movers_[ I ]->set_seq_sep_ramping( 0.0, begin_seqsep_factor );

	// Calculate number of interations for each phase ------------------------------
	ITERATIONS[ 1 ] = 1;
	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		stage_movers_[ I ]->set_seq_sep_ramping( end_seqsep_factor - begin_seqsep_factor, begin_seqsep_factor );
		ITERATIONS[ 1 ] += (end_seqsep - begin_seqsep) / 2;
	}

	ITERATIONS[ 2 ] = 1;
	ITERATIONS[ 3 ] = 10;
	ITERATIONS[ 4 ] = 3;

	return ITERATIONS;
}

EnvClaims AbscriptMover::yield_claims( core::pose::Pose const&,
	basic::datacache::WriteableCacheableMapOP ) {
	return claims_;
}

void AbscriptMover::register_preparer( protocols::moves::MoverOP mover, StageIDs const& ids ){
	StagePreparerOP preparer = utility::pointer::dynamic_pointer_cast< protocols::abinitio::abscript::StagePreparer > ( mover );
	if ( !preparer ) {
		tr.Error << "The mover '" << mover->get_name() << "' is not a preparer." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("Non-preparer mover added as preparer to Abscript.");
	}

	for ( StageIDs::const_iterator id = ids.begin(); id != ids.end(); ++id ) {
		stage_movers_[ *id ]->add_preparer( preparer );
	}

}

std::string AbscriptMover::get_name() const {
	return mover_name();
}

std::string AbscriptMover::mover_name() {
	return "AbscriptMover";
}

std::string AbscriptMover::stage_complex_namer( std::string tag_name ){
	return "abscript_stage_" + tag_name + "_complex_type";
}

void AbscriptMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"name", xs_string,
		"Prefix for the name of the scorefunction. Scorefunctions should have been previously defined with names prefix + '_stage' + I, II, IIIa, IIIb, IVa, and IVb.");
	attlist + XMLSchemaAttribute(
		"cycles", xsct_real,
		"Equivalent to the -run:increase_cycles in standard ab initio; increases the number of ab initio cycles by that factor.");
	attlist + XMLSchemaAttribute(
		"skip_stages", xsct_positive_integer_cslist,
		"Comma-separated list of ab initio stages (1, 2, 3, or 4) that you wish to skip");

	XMLSchemaSimpleSubelementList subelements;
	//XMLSchemaSimpleSubelementList subelements_fragments, subelements_stages;
	XMLSchemaSimpleSubelementList stage_subelements;
	AttributeList stage_attlist;
	AttributeList mover_attlist;
	AttributeList preparer_attlist;
	stage_attlist + XMLSchemaAttribute::required_attribute(
		"ids", xs_string,
		"Stage of ab intio to which this subtag refers.");
	mover_attlist
		+ XMLSchemaAttribute::attribute_w_default( "name", xs_string, "Name of a mover that this stage should use with the specified weight relative to other specified movers", "null" )
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "How heavily should this mover be weighted relative to other movers in this stage, i.e. how frequently should it be used relative to others?", "1.0" );
	preparer_attlist
		+ XMLSchemaAttribute::attribute_w_default( "name", xs_string, "Name of the mover that this stage should use to prepare the structure", "null" );
	stage_subelements
		.add_simple_subelement( "Mover", mover_attlist, "Tag specifying a mover that should be used in this stage of the ab initio protocol" )
		.add_simple_subelement( "Preparer", preparer_attlist, "Tag specifying a StagePreparer that should be run before other movers are sampled during this stage" );

	AttributeList frag_attlist;

	frag_attlist + XMLSchemaAttribute(
		"selector", xs_string,
		"Name of previously defined residue selector specifying where fragments should be inserted")
		+ XMLSchemaAttribute(
		"small_frags", xs_string,
		"Fragments file containing small fragments (i.e. 3mers) ")
		+ XMLSchemaAttribute(
		"large_frags", xs_string,
		"Fragments file containing large fragments (i.e. 9mers)")
		+ XMLSchemaAttribute::attribute_w_default(
		"initialize", xsct_rosetta_bool,
		"Use small fragments insert a fragment at all positions before starting",
		"true");

	XMLSchemaComplexTypeGenerator stage_ctgen;
	stage_ctgen.element_name( "Stage" )
		.set_subelements_repeatable( stage_subelements )
		.add_attributes( stage_attlist )
		.complex_type_naming_func( & stage_complex_namer )
		.description( "Tag specifying behavior for a specific stage of the ab initio protocol" )
		.write_complex_type_to_schema( xsd );
	subelements.add_already_defined_subelement( "Stage", & stage_complex_namer );
	//subelements_stages.add_already_defined_subelement( "Stage", & stage_complex_namer );

	subelements.add_simple_subelement( "Fragments", frag_attlist, "Macro used to add the appropriate ClassicFragmentMovers for insertion of large fragments, normal insertion of small fragments, and smooth insertion of small fragments" );
	//subelements_fragments.add_simple_subelement("Fragments", frag_attlist, "Macro used to add the appropriate ClassicFragmentMovers for insertion of large fragments, normal insertion of small fragments, and smooth insertion of small fragments" ); //, 0, 1);

	//XMLSchemaComplexTypeGenerator whole_shebang;
	//whole_shebang.complex_type_naming_func( & moves::complex_type_name_for_mover )
	// .element_name( mover_name() )
	// .description( "A special mover used to replicate the state of ab initio in early 2014" )
	// .add_attributes( attlist )
	// //.add_optional_name_attribute() -- Name attribute has already been added as required.
	// //.set_subelements_repeatable( subelements )
	// .add_ordered_subelement_set_as_optional( subelements_fragments )
	// .add_ordered_subelement_set_as_repeatable( subelements_stages )
	// .write_complex_type_to_schema( xsd );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"A special mover used to replicate the state of ab initio in early 2014",
		attlist, subelements);
}

std::string AbscriptMoverCreator::keyname() const {
	return AbscriptMover::mover_name();
}

protocols::moves::MoverOP
AbscriptMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AbscriptMover );
}

void AbscriptMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AbscriptMover::provide_xml_schema( xsd );
}


AbscriptMover::StageTracker::StageTracker( moves::MonteCarloOP mc ):
	mc_( mc ),
	checkpointer_( "Abscript" )
{}

void AbscriptMover::StageTracker::begin_stage( std::string const& stagename,
	core::pose::Pose& pose ){
	stagename_ = stagename;
	iterations_ = 0;

	if ( stagename_ == "I" &&
			basic::options::option[ basic::options::OptionKeys::run::intermediate_structures ]() ) {
		protocols::jd2::output_intermediate_pose( pose, "stage0" );
	}

	// part X ----------------------------------------
	tr.Info <<  "\n===================================================================\n";
	tr.Info <<  "   Stage " << stagename << "\n";
	tr.Info <<  "--------------------------------------------------------------------" << std::endl;
	// PROF_START( stagename() );
	starttime_ = clock();
}

void AbscriptMover::StageTracker::do_step( core::pose::Pose& pose,
	AbscriptStageMover& mover,
	core::Real const& progress ){
	iterations_ += 1;

	//std::string checkpoint_id = "stage_"+stagename_+"_it_"+utility::to_string( iterations_ );

	//  Checkpointing goes here, but doesn't work with ProtectedConformation yet, and so isn't implemented.
	//  if ( !checkpointer_.recover_checkpoint( pose, protocols::jd2::current_output_name(), checkpoint_id,
	//                                          false /* fullatom */, true /* fold_tree */ ) ) {
	if ( mover.setup_stage( pose, progress ) ) {
		mover.apply( pose );
	} else {
		tr.Info << "Stage " << stagename_ << " step at " << progress << "% skipped." << std::endl;
	}
	//    checkpointer_.checkpoint( pose, protocols::jd2::current_output_name(), checkpoint_id, true /*fold tree */ );
	//  }
	//  checkpointer_.debug( protocols::jd2::current_output_name(), checkpoint_id, pose.energies().total_energy() );
}

void AbscriptMover::StageTracker::end_stage( core::pose::Pose& pose ){

	if ( tr.Info.visible() ) {
		tr.Info << "\n";
		mc_->score_function().show( tr.Info, pose );
	}
	tr.Info << std::endl;

	mc_->show_counters();
	mc_->reset_counters();

	clock_t endtime = clock();
	// PROF_STOP( stagename() );
	if ( basic::options::option[ basic::options::OptionKeys::run::profile ] ) basic::prof_show();
	if ( basic::options::option[ basic::options::OptionKeys::abinitio::debug ]() ) {
		tr.Info << "Time/step: " << ( double(endtime) - starttime_ )/( CLOCKS_PER_SEC ) << std::endl;
	}

	if ( basic::options::option[ basic::options::OptionKeys::run::intermediate_structures ]() &&
			iterations_ > 0 ) {
		protocols::jd2::output_intermediate_pose( pose, "stage"+stagename_ );
	}
}

} // abscript
} // abinitio
} // protocols

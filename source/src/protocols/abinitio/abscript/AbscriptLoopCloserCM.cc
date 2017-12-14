// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/AbscriptLoopCloserCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/AbscriptLoopCloserCM.hh>
#include <protocols/abinitio/abscript/AbscriptLoopCloserCMCreator.hh>

// Package headers
#include <protocols/environment/claims/TorsionClaim.hh>

#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/EnvExcn.hh>

#include <core/environment/DofPassport.hh>

// Project headers
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/select/residue_selector/TrueResidueSelector.hh>

#include <protocols/loops/loop_closure/ccd/WidthFirstSlidingWindowLoopClosure.hh>
#include <protocols/loops/Exceptions.hh>

#include <protocols/checkpoint/CheckPointer.hh>

#include <protocols/jumping/util.hh>

#include <protocols/idealize/IdealizeMover.hh>

#include <protocols/rosetta_scripts/util.hh>

#ifdef WIN32
#include <basic/datacache/WriteableCacheableMap.hh>
#endif

//Utility Headers

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>

#include <utility>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

#include <boost/functional/hash.hpp>
#include <boost/foreach.hpp>

// tracer
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr( "protocols.abinitio.abscript.AbscriptLoopCloserCM", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace protocols::environment;

// creator
// XRW TEMP std::string
// XRW TEMP AbscriptLoopCloserCMCreator::keyname() const {
// XRW TEMP  return AbscriptLoopCloserCM::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AbscriptLoopCloserCMCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AbscriptLoopCloserCM );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AbscriptLoopCloserCM::mover_name() {
// XRW TEMP  return "AbscriptLoopCloserCM";
// XRW TEMP }

AbscriptLoopCloserCM::AbscriptLoopCloserCM():
	Parent(),
	fragset_(),
	scorefxn_(),
	selector_( new core::select::residue_selector::TrueResidueSelector() )
{}

AbscriptLoopCloserCM::AbscriptLoopCloserCM( core::fragment::FragSetCOP fragset,
	core::scoring::ScoreFunctionOP scorefxn ):
	Parent(),
	fragset_(std::move( fragset )),
	scorefxn_(std::move( scorefxn )),
	selector_( new core::select::residue_selector::TrueResidueSelector() )
{}

claims::EnvClaims AbscriptLoopCloserCM::yield_claims( core::pose::Pose const&,
	basic::datacache::WriteableCacheableMapOP ){
	claims::EnvClaims claims;

	// We want to control everything that will be relevant to the output pose (which should be the same as the input.
	claims::TorsionClaimOP claim( new claims::TorsionClaim( utility::pointer::static_pointer_cast< ClientMover > ( get_self_ptr() ), selector() ) );
	claim->strength( claims::CAN_CONTROL, claims::DOES_NOT_CONTROL );

	claims.push_back( claim );

	return claims;
}

core::Real angle_diff( core::Real const& a, core::Real const& b ){
	//use arcsin to get back degrees, but need an f(0)=0, so sin not cos.
	return std::asin( std::abs( std::cos( a ) - std::cos ( b ) ) );
}

bool angle_cpy( core::pose::Pose& target, core::pose::Pose const& source, core::id::TorsionID t_id ) {
	// Empirically, a tolerance of 2 has worked for me. It's in degrees, so it's small, but not small enough
	// for my taste. Ideally, I would run some kind of minimization procedure to try to pull target back to its
	// original configuration, but I'm not sure how to do that. Typically, these structures will be all-atom
	// relaxed after this anyway, so that helps.
	core::Real const TOLERANCE = 2;
	std::string angle;
	if ( t_id.torsion() == core::id::psi_torsion ) {
		angle = "PSI";
	} else if ( t_id.torsion() == core::id::phi_torsion ) {
		angle = "PHI";
	} else if ( t_id.torsion() == core::id::omega_torsion ) {
		angle = "OMEGA";
	} else {
		angle = "OTHER";
	}

	if ( angle_diff( target.torsion( t_id ), source.torsion( t_id ) ) > TOLERANCE ) {
		try {
			target.set_torsion( t_id, source.torsion( t_id ) );
		} catch ( EXCN_Env_Security_Exception& ){
			std::ostringstream ss;
			ss << "[ERROR] Loop closure tried to write to residue " << t_id.rsd() << " " << angle << " angle ("
				<< t_id << "). Target angle was " << target.torsion( t_id ) << " and source angle was "
				<< source.torsion( t_id ) << " (delta=" << angle_diff( target.torsion( t_id ), source.torsion( t_id ) )
				<< "; tolerance is " << TOLERANCE << ")." << std::endl;
			throw CREATE_EXCEPTION(protocols::loops::EXCN_Loop_not_closed,  ss.str() );
		}
	} else {
		tr.Debug << "  ignoring modified " << angle << " @ " << t_id.torsion() << " (input=" << source.torsion( t_id ) << ", output=" << target.torsion( t_id ) << std::endl;
		return false;
	}

	return true;
}

void AbscriptLoopCloserCM::apply( core::pose::Pose& in_pose ){

	debug_assert( passport() );
	debug_assert( final_ft_ );

	//Produce unprotected pose
	core::pose::Pose pose( in_pose );
	core::conformation::ConformationOP deprotected( new core::conformation::Conformation( in_pose.conformation() ) );
	pose.set_new_conformation( deprotected );

	tr.Debug << "Closing loops with current fold tree: " << in_pose.fold_tree();
	tr.Debug << "Using final fold tree: " << *final_ft_ << std::endl;

	// Close Loops
	bool success;

	attempt_idealize( pose );

	success = attempt_ccd( pose );

	if ( success ) {
		attempt_idealize( pose );
	} else {
		// We could implement "don't fail unclosed" here, but I don't think we need to?
		return;
	}

	// Copy result into protected conformation in in_pose
	DofUnlock unlock( in_pose.conformation(), passport() );

	Size warn = 0;

	try {
		for ( Size i = 1; i <= in_pose.size(); ++i ) {
			tr.Debug << "Movemap for " << in_pose.residue( i ).name3() << i
				<< ": bb(" << ( movemap_->get_bb( i ) ? "T" : "F" )
				<< ") chi(" << ( movemap_->get_chi( i ) ? "T" : "F" ) << ")" << std::endl;

			if ( angle_cpy( in_pose, pose, core::id::TorsionID( i, core::id::BB, core::id::omega_torsion ) ) ) {
				warn += 1;
			}
			if ( angle_cpy( in_pose, pose, core::id::TorsionID( i, core::id::BB, core::id::phi_torsion ) ) ) {
				warn += 1;
			}
			if ( angle_cpy( in_pose, pose, core::id::TorsionID( i, core::id::BB, core::id::psi_torsion ) ) ) {
				warn += 1;
			}
		}
	} catch( utility::excn::Exception& e ) {
		tr.Warning << this->get_name() << " failed to close loops without violating the ClientMover contract. Report: "
			<< e.msg() << std::endl;
		this->set_last_move_status( moves::FAIL_DO_NOT_RETRY );
	}
	if ( warn ) {
		tr.Warning << "AbscriptLoopMoverCM ignoring " << warn
			<< " values that differ between idealized and unidealized poses. Use '-out:levels protocols.abinitio.abscript.AbscriptLoopCloserCM:debug' to see more info." << std::endl;
	}
}

bool AbscriptLoopCloserCM::attempt_ccd( core::pose::Pose& pose ){
	using namespace loops::loop_closure::ccd;

	try {
		checkpoint::CheckPointer checkpointer( "AbscriptLoopCloserCM" );
		SlidingWindowLoopClosureOP closing_protocol;
		closing_protocol = SlidingWindowLoopClosureOP( new WidthFirstSlidingWindowLoopClosure( fragset_,
			scorefxn_,
			movemap_ ) );

		jumping::close_chainbreaks( closing_protocol,
			pose,
			checkpointer,
			get_current_tag(),
			*final_ft_ );
	} catch (loops::EXCN_Loop_not_closed& excn ) {
		tr.Warning << this->get_name() << " failed to close a loop. Report: " << excn << std::endl;
		set_last_move_status( moves::FAIL_DO_NOT_RETRY );
		set_current_tag( "C_"+get_current_tag().substr(std::min(2,(int)get_current_tag().size())) );
		return false;
	}

	return true;
}

void AbscriptLoopCloserCM::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const&,
	protocols::moves::Movers_map const&,
	core::pose::Pose const& ) {

	using namespace basic::options::OptionKeys;
	using namespace basic::options;

	if ( tag->hasOption( "selector" ) ) {
		set_selector( datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", tag->getOption<std::string>( "selector" ) ) );
	} else {
		set_selector( core::select::residue_selector::ResidueSelectorCOP( new core::select::residue_selector::TrueResidueSelector() ) );
	}

	std::string fragfile;
	if ( tag->hasOption("fragments") ) {
		fragfile = tag->getOption<std::string>("fragments");
	} else if ( !option[OptionKeys::in::file::frag3].user() ) {
		tr.Error <<"enter frags with either fragments= or -in:file:frag3"  << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "fragment file wrong type" );
	} else {
		fragfile = option[ OptionKeys::in::file::frag3 ]();
	}

	core::fragment::FragmentIO frag_io( option[ OptionKeys::abinitio::number_3mer_frags ](), 1,
		option[ OptionKeys::frags::annotate ]() );

	fragset_ = frag_io.read_data( fragfile );

	if ( tag->hasOption( "scorefxn" ) ) {
		try {
			scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );
		} catch ( ... ) {
			tr.Error << "AbscriptLoopCloserCM failed to find the score function '"
				<< tag->getOption< std::string >( "scorefxn" ) << std::endl;
			throw;
		}
	} else {
		tr.Warning << "Configuring AbscriptLoopCloserCM '" << tag->getName() << "' with default abinitio loop closure score function." << std::endl;

		option[ OptionKeys::abinitio::stage4_patch ].activate();

		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage4_patch ]() );
		scorefxn_->set_weight( core::scoring::linear_chainbreak, option[ jumps::chainbreak_weight_stage4 ]() );
		scorefxn_->set_weight( core::scoring::atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ] );
	}
}

void AbscriptLoopCloserCM::attempt_idealize( core::pose::Pose& pose ) {
	idealize::IdealizeMoverOP idealizer( new idealize::IdealizeMover );
	idealizer->fast( false );

	utility::vector1< core::Size > pos_list;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( movemap_->get_bb( i ) ) {
			pos_list.push_back( i );
		}
	}

	if ( tr.Trace.visible() ) {
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			tr.Trace << ( movemap_->get_bb(i) ? "T" : "F" );
		}
		tr.Trace << std::endl;
		tr.Trace << "movable positions: " << pos_list << std::endl;
	}


	if ( pos_list.size() == 0 &&
			pose.fold_tree().num_cutpoint() - final_ft_->num_cutpoint() > 0 ) {
		std::ostringstream ss;
		ss << this->get_name() << " tried to start closing the fold tree by idealization ("
			<< pose.fold_tree().num_cutpoint() - final_ft_->num_cutpoint()
			<< " cuts are to be closed) but was not given any movable positions." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
	}

	idealizer->set_pos_list( pos_list );

	core::kinematics::FoldTree fold_tree( pose.fold_tree() );
	pose.fold_tree( *final_ft_ );
	if ( pos_list.size() > 0 ) { // req'd because an empty pos_list is treated as all
		idealizer->apply( pose );
	}
	pose.fold_tree( fold_tree );
	(*scorefxn_)( pose );
	//jd2::output_intermediate_pose( pose, "loops_closed_preprocessed" );
}

void AbscriptLoopCloserCM::passport_updated() {
	if ( has_passport() ) {
		movemap_ = passport()->render_movemap();
	} else {
		movemap_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );
		movemap_->set_bb( true );
	}
}


void AbscriptLoopCloserCM::broking_finished( EnvClaimBroker::BrokerResult const& broker ) {
	//set the target fold tree for the closer.
	final_ft_ = broker.closer_ft;
	debug_assert( final_ft_ );
}

// XRW TEMP std::string AbscriptLoopCloserCM::get_name() const {
// XRW TEMP  return "AbscriptLoopCloserCM";
// XRW TEMP }

moves::MoverOP AbscriptLoopCloserCM::clone() const {
	return moves::MoverOP( new AbscriptLoopCloserCM( *this ) );
}

std::string AbscriptLoopCloserCM::get_name() const {
	return mover_name();
}

std::string AbscriptLoopCloserCM::mover_name() {
	return "AbscriptLoopCloserCM";
}

void AbscriptLoopCloserCM::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"selector", xs_string,
		"Residue selector specifying the region where fragments will be inserted");
	attlist + XMLSchemaAttribute(
		"fragments", xs_string,
		"Path to fragment file containing fragments that will be used to close the loop");

	rosetta_scripts::attributes_for_parse_score_function(attlist);
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Uses the WidthFirstSlidingWindowLoopCloser to fix loops using fragment insertion in the specified region",
		attlist );
}

std::string AbscriptLoopCloserCMCreator::keyname() const {
	return AbscriptLoopCloserCM::mover_name();
}

protocols::moves::MoverOP
AbscriptLoopCloserCMCreator::create_mover() const {
	return protocols::moves::MoverOP( new AbscriptLoopCloserCM );
}

void AbscriptLoopCloserCMCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AbscriptLoopCloserCM::provide_xml_schema( xsd );
}


} // abscript
} // abinitio
} // protocols

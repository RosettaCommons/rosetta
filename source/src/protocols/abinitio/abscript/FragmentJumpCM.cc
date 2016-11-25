// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/FragmentJumpCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/FragmentJumpCM.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <core/environment/LocalPosition.hh>
#include <core/environment/LocalPosition.hh>

#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/TorsionClaim.hh>

#include <protocols/abinitio/abscript/JumpSampleData.hh>
#include <protocols/abinitio/abscript/FragmentJumpCMCreator.hh>

// Project headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/Frame.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/abinitio/PairingStatistics.hh>
#include <protocols/abinitio/TemplateJumpSetup.hh>
#include <protocols/abinitio/Templates.hh>

#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/JumpSample.hh>
#include <protocols/jumping/RandomSheetBuilder.hh>

#include <core/scoring/dssp/PairingsList.hh>

//Utility Headers
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>

#include <boost/functional/hash.hpp>

// tracer
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer tr( "protocols.abinitio.abscript.FragmentJumpCM", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace protocols::environment;

// creator
// XRW TEMP std::string
// XRW TEMP FragmentJumpCMCreator::keyname() const {
// XRW TEMP  return FragmentJumpCM::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FragmentJumpCMCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FragmentJumpCM );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FragmentJumpCM::mover_name() {
// XRW TEMP  return "FragmentJumpCM";
// XRW TEMP }

FragmentJumpCM::FragmentJumpCM():
	Parent(),
	moverkey_( "" )
{}

FragmentJumpCM::FragmentJumpCM( std::string const& topol_filename,
	core::select::residue_selector::ResidueSelectorCOP selector,
	std::string const& moverkey ) :
	moverkey_( moverkey )
{
	set_topology( topol_filename );
	Parent( mover(), selector );
}

void FragmentJumpCM::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const&,
	protocols::moves::Movers_map const&,
	core::pose::Pose const& ) {
	if ( tag->hasOption( "topol_file" ) ) {
		std::string const& topol_filename = tag->getOption< std::string >( "topol_file",
			"[VALUE_UNSET]" );
		set_topology( topol_filename );
	} else if ( tag->hasOption( "ss_info" ) &&
			tag->hasOption( "pairing_file" ) &&
			tag->hasOption( "n_sheets" ) ) {
		std::string const& ss_info_file = tag->getOption< std::string >( "ss_info",
			"[VALUE_UNSET]" );
		std::string const& pairing_file = tag->getOption< std::string >( "pairing_file",
			"[VALUE_UNSET]" );
		core::Size const& n_sheets = tag->getOption< core::Size >( "n_sheets", 0 );
		bool bRandomSheets = tag->getOption< bool >( "random_sheets", true );

		set_topology( ss_info_file, pairing_file, n_sheets, bRandomSheets );
	} else if ( !tag->getOption< bool >( "restart_only", false ) ) {
		tr.Error << "[ERROR] Parsing problem in FragmentJumpCM: "
			<< " option set not compatible. Valid sets are 'topol_file' "
			<< " or 'ss_info', 'pairing_file', and 'n_sheets'." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( "FragmentJumpCM incompatible options." );
	}

	initialize( tag->getOption< bool >( "initialize", true ) );

	if ( tag->hasOption( "selector" ) ) {
		set_selector( datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", tag->getOption<std::string>( "selector" ) ) );
	} else {
		set_selector( NULL );
	}

	set_moverkey( tag->getOption< std::string >( "name" , "" ) );

	tr.Debug << "Initialzed " << get_name();
	if ( selector() ) {
		tr.Debug << " with ResidueSelector";
	} else {
		tr.Debug << "without ResidueSelector";
	}
	tr.Debug << " and moverkey " << moverkey() << std::endl;

}

claims::EnvClaims FragmentJumpCM::yield_claims( core::pose::Pose const& pose,
	basic::datacache::WriteableCacheableMapOP map ){
	using namespace core::pose::datacache;
	using namespace basic::datacache;

	utility::vector1< bool > selection;
	if ( selector() ) {
		selection = selector()->apply( pose );
	} else {
		selection = utility::vector1< bool >( pose.size(), true );
	}

	if ( map->find( "JumpSampleData") != map->end() ) {
		std::set< WriteableCacheableDataOP > const& data_set = map->find( "JumpSampleData" )->second;

		for ( WriteableCacheableDataOP data_ptr : data_set ) {
			JumpSampleDataOP jumpdata_ptr = utility::pointer::dynamic_pointer_cast< protocols::abinitio::abscript::JumpSampleData > ( data_ptr );
			assert( jumpdata_ptr );

			if ( jumpdata_ptr->moverkey() == moverkey() ) {
				tr.Debug << "Found JumpSampleOP in cache with key '" << moverkey() << "'" << std::endl;
				if ( jump_def_ ) {
					tr.Warning << "Found JumpSampleOP found in FragmentJumpCM's input pose. Overwriting input"
						<< " topologies with input structure's information." << std::endl;
					jump_def_ = NULL;
				}
				return build_claims( selection, jumpdata_ptr->jump_sample() );
			} else if ( tr.Debug.visible() ) {
				tr.Debug << "Rejected JumpSampleDataOP with key " << jumpdata_ptr->moverkey()
					<< ", which conflicts with present key " << moverkey() << std::endl;
			}
		}
	} else {
		tr.Debug << get_name() << " found no JumpSampleData objects in input pose cache." << std::endl;
	}

	//If we don't hit the return statement in the above block of code, we have to calculate a new jump set.
	tr.Debug << get_name() << " calculating new jumps based on input topology file" << std::endl;

	jumping::JumpSample jump_sample;
	try {
		// called for each to fold tree get random jumps from top file.
		jump_sample = calculate_jump_sample();
	} catch ( ... ) {
		tr.Error << "[ERROR] " << get_name() << " failed to generate a jump sample from fragments. "
			<< "Were you expecting JumpSampleData in the input file? No fitting data was found."
			<< std::endl;
		throw;
	}

// cache the result as a WriteableCacheable in the pose DataCache for later retrieval.
// this is important for the restarting feature in Abinitio.
	JumpSampleDataOP data( new JumpSampleData( moverkey(), jump_sample ) );
	(*map)[ data->datatype() ].insert( data );

	return build_claims( selection, jump_sample );
}

claims::EnvClaims FragmentJumpCM::build_claims( utility::vector1< bool > const& residue_selection,
	jumping::JumpSample const& jump_sample ) {
	claims::EnvClaims claim_list;

	//Now that we've committed to a particular sample, configure fragments in the mover.
	setup_fragments( jump_sample );
	tr.Debug <<"mover.configured" << mover()->get_name() << std::endl;
	int shift = residue_selection.index( true )-1;
	ClientMoverOP this_ptr = utility::pointer::static_pointer_cast< ClientMover > ( get_self_ptr() );

	for ( int i = 1; i <= (int) jump_sample.size(); i++ ) {
		Size const up = jump_sample.jumps()( 1, i ) + shift;
		Size const dn = jump_sample.jumps()( 2, i ) + shift;

		std::string jump_name = get_name() + "Jump" + utility::to_string( i );

		claims::JumpClaimOP jclaim( new claims::JumpClaim( this_ptr,
			jump_name,
			LocalPosition( "BASE", up ),
			LocalPosition( "BASE", dn ) ) );
		if ( initialize() ) {
			jclaim->strength( claims::MUST_CONTROL, claims::CAN_CONTROL );
		} else {
			jclaim->strength( claims::MUST_CONTROL, claims::DOES_NOT_CONTROL );
		}

		jclaim->set_atoms( "", "" ); //allow jump atoms to be determined by ft at broker-time.
		jclaim->physical( false ); //jumps are not physical, and should be scored as chainbreaks

		jclaim->stubs_intra_residue( true );

		claim_list.push_back( jclaim );

		// Jump Fragments make the mover into a bit of an access primadonna because jump fragments
		// include torsions from the residues at takeoff and landing
		claims::TorsionClaimOP tclaim_up( new claims::TorsionClaim( this_ptr, LocalPosition( "BASE", up ) ) );
		claims::TorsionClaimOP tclaim_dn( new claims::TorsionClaim( this_ptr, LocalPosition( "BASE", dn ) ) );

		tclaim_up->strength( claims::CAN_CONTROL, claims::CAN_CONTROL );
		tclaim_dn->strength( claims::CAN_CONTROL, claims::CAN_CONTROL );

		claim_list.push_back( tclaim_up );
		claim_list.push_back( tclaim_dn );
	}

	return claim_list;
}

// XRW TEMP std::string FragmentJumpCM::get_name() const {
// XRW TEMP  return "FragmentJumpCM";
// XRW TEMP }

void FragmentJumpCM::set_topology( std::string const& ss_info_file,
	std::string const& pairing_file,
	core::Size const& n_sheets,
	bool bRandomSheets ) {
	if ( mover() ) {
		tr.Warning << "[WARNING] internal ClassicFragmentMover overwritten during FragmentJumpCM::process_topology_file call." << std::endl;
	}

	jumping::SheetBuilder::SheetTopology sheets;
	sheets.push_back( n_sheets );

	core::fragment::SecondaryStructureOP ss_def( new core::fragment::SecondaryStructure );
	ss_def->read_psipred_ss2( ss_info_file );

	core::scoring::dssp::PairingList pairlist;
	core::scoring::dssp::read_pairing_list( pairing_file, pairlist );

	if ( bRandomSheets ) {
		jump_def_ = jumping::BaseJumpSetupOP( new jumping::RandomSheetBuilder( ss_def, pairlist, sheets ) );
	} else {
		jump_def_ = jumping::BaseJumpSetupOP( new jumping::SheetBuilder( ss_def, pairlist, sheets ) );
	}

	// Fail faster with better information on bad topology files.
	jumping::JumpSample jump_sample;
	try{
		jump_sample = calculate_jump_sample();
	} catch( utility::excn::EXCN_BadInput ){
		std::stringstream ss;
		ss << "Was not able to construct a valid jump sample in 10 attempts using ss_info file "
			<< ss_info_file << ", pairing file " << pairing_file << " and " << n_sheets
			<< "( sheets are random = " << bRandomSheets << ".";

		tr.Error << ss.str() << std::endl;

		throw utility::excn::EXCN_BadInput( ss.str() );
	}
}

void FragmentJumpCM::set_topology( std::string const& topol_filename ){

	if ( mover() ) {
		tr.Warning << "[WARNING] internal ClassicFragmentMover overwritten during FragmentJumpCM::process_topology_file call."
			<< std::endl;
	}

	utility::io::izstream is( topol_filename );
	if ( !is.good() ) {
		throw utility::excn::EXCN_FileNotFound(" Topology file '" + topol_filename +"' not found." );
	}

	abinitio::PairingStatisticsOP ps( new abinitio::PairingStatistics );
	is >> *ps;
	tr.Info << *ps << std::endl;
	core::fragment::SecondaryStructureOP ss_def( new core::fragment::SecondaryStructure );
	ss_def->extend( 10000 ); //Set number of residues to unreasonably large.
	core::scoring::dssp::PairingList helix_pairings; // helix pairings not used, required by BaseJumpSetup.

	jump_def_ = jumping::BaseJumpSetupOP( new abinitio::TemplateJumpSetup( NULL, ss_def, ps, helix_pairings ) );

	// Fail faster with better information on bad topology files.
	jumping::JumpSample jump_sample;
	try{
		jump_sample = calculate_jump_sample();
	} catch( utility::excn::EXCN_BadInput ){
		throw utility::excn::EXCN_BadInput( "Was not able to construct a valid jump sample in 10 attempts using topology file "
			+ topol_filename + "." );
	}
}

void FragmentJumpCM::setup_fragments( jumping::JumpSample const& jump_sample ) {

	core::kinematics::MoveMapOP dummy_mm( new core::kinematics::MoveMap );
	dummy_mm->set_bb( true );
	dummy_mm->set_jump( true );

	core::fragment::FragSetOP jump_frags;
	if ( jump_def_ ) {
		jump_frags = jump_def_->generate_jump_frags( jump_sample, *dummy_mm );
	} else {
		jump_frags = core::fragment::FragSetOP( new core::fragment::OrderedFragSet );
		core::fragment::FrameList jump_frames;
		jump_sample.generate_jump_frames( jump_frames, *dummy_mm );
		jump_frags->add( jump_frames );
	}


	if ( !mover() ) {
		simple_moves::ClassicFragmentMoverOP mover( new simple_moves::ClassicFragmentMover( jump_frags, dummy_mm ) );
		mover->set_check_ss( false ); //for some reason it's all 'L' notated, which causes a rejection of fragments...
		mover->enable_end_bias_check( false );
		set_mover( mover );
	} else {
		mover()->set_fragments( jump_frags );
	}
}

jumping::JumpSample FragmentJumpCM::calculate_jump_sample() const {
	jumping::JumpSample jump_sample;

	if ( !jump_def_ ) {
		tr.Error << "[ERROR]" << get_name() << " tried to make jumps but couldn't because no "
			<< "(appropriate) JumpSampleData was cached in the pose and no topology information was "
			<< "provided through parse_my_tag." << std::endl;
		throw utility::excn::EXCN_BadInput( "FragmentJumpCM requires, but did not have, a JumpSample." );
	}

	core::Size attempts = 10;
	while ( !jump_sample.is_valid() ) {
		attempts -= 1;
		if ( attempts <= 0 ) {
			throw utility::excn::EXCN_BadInput( "Was not able to construct a valid jump sample in 10 attempts." );
		}

		jump_sample = jump_def_->create_jump_sample();
	}

	return jump_sample;
}

moves::MoverOP FragmentJumpCM::fresh_instance() const {
	return moves::MoverOP( new FragmentJumpCM() );
}

moves::MoverOP FragmentJumpCM::clone() const{
	return moves::MoverOP( new FragmentJumpCM( *this ) );
}

std::string FragmentJumpCM::get_name() const {
	return mover_name();
}

std::string FragmentJumpCM::mover_name() {
	return "FragmentJumpCM";
}

void FragmentJumpCM::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"topol_file", xs_string,
		"File specifying topology for the pose. If not provided, must specify ss_info, pairing_file, and n_sheets");
	attlist + XMLSchemaAttribute(
		"ss_info", xs_string,
		"File containing secondary structure information for the pose. Incompatible with topol_file.");
	attlist + XMLSchemaAttribute(
		"pairing_file", xs_string,
		"File containing info on secondary structure pairing for the pose.");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"n_sheets", xsct_non_negative_integer,
		"Number of beta sheets", "0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"random_sheets", xsct_rosetta_bool,
		"Should beta sheets be built randomly?", "true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"restart_only", xsct_rosetta_bool,
		"Do not reset the topology for this pose", "false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"initialize", xsct_rosetta_bool,
		"Apply FragmentJumpCM after brokering is complete?", "true");
	attlist + XMLSchemaAttribute(
		"selector", xs_string,
		"Residue selector specifying region where this mover will be applied");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"name", xs_string,
		"Unique name for this client mover", "");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Inserts strand-strand rigid body translations into jumps between predicted adjacent beta strands.",
		attlist );
}

std::string FragmentJumpCMCreator::keyname() const {
	return FragmentJumpCM::mover_name();
}

protocols::moves::MoverOP
FragmentJumpCMCreator::create_mover() const {
	return protocols::moves::MoverOP( new FragmentJumpCM );
}

void FragmentJumpCMCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FragmentJumpCM::provide_xml_schema( xsd );
}


} // abscript
} // abinitio
} // protocols

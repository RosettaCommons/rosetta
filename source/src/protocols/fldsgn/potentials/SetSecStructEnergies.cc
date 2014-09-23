// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/fldsgn/potentials/SetSecStructEnergy.cc
/// @brief mover for setting centroid score of secondary structure through parser
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// unit headers
#include <protocols/fldsgn/potentials/SetSecStructEnergies.hh>
#include <protocols/fldsgn/potentials/SetSecStructEnergiesCreator.hh>

// package headers
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
//#include <core/conformation/symmetry/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasSecondaryStructureEnergy.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasHelicesSheetPotential.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <utility/tag/Tag.hh>
// C++ headers
#include <utility>

// boost
// AUTO-REMOVED #include <boost/lexical_cast.hpp>

#include <protocols/fldsgn/potentials/sspot/NatbiasStrandPairPotential.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>



namespace protocols {
namespace fldsgn {
namespace potentials {

static thread_local basic::Tracer TR( "protocols.fldsgn.SetSecStructEnergies" );

std::string
SetSecStructEnergiesCreator::keyname() const
{
	return SetSecStructEnergiesCreator::mover_name();
}

protocols::moves::MoverOP
SetSecStructEnergiesCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetSecStructEnergies );
}

std::string
SetSecStructEnergiesCreator::mover_name()
{
	return "SetSecStructEnergies";
}


/// @brief default constructor
SetSecStructEnergies::SetSecStructEnergies() :
	Super( "SetSecStructEnergies" ),
	loaded_( false ),
	blueprint_( /* NULL */ ),
	ss_from_blueprint_( true ),
	sfx_( /* NULL */ ),
	hh_weight_( 1.0 ),
	hs_weight_( 1.0 ),
	ss_weight_( 1.0 ),
	stwist_weight_( 1.0 ),
	hs_pair_weight_( 0.0 ),
	ss_pair_weight_( 0.0 ),
	rsigma_weight_( 0.0 ),
	hpairpot_( NatbiasHelixPairPotentialOP( new NatbiasHelixPairPotential ) ),
	hspot_( NatbiasHelicesSheetPotentialOP( new NatbiasHelicesSheetPotential ) )
{}


/// @brief value constructor
SetSecStructEnergies::SetSecStructEnergies( ScoreFunctionOP const sfx, String const & filename, bool const ss_from_blueprint ) :
	Super( "SetSecStructEnergies" ),
	loaded_( false ),
	ss_from_blueprint_( ss_from_blueprint ),
	sfx_( sfx ),
	hh_weight_( 1.0 ),
	hs_weight_( 1.0 ),
	ss_weight_( 1.0 ),
	stwist_weight_( 1.0 ),
	hs_pair_weight_( 0.0 ),
	ss_pair_weight_( 0.0 ),
	rsigma_weight_( 0.0 ),
	hpairpot_( NatbiasHelixPairPotentialOP( new NatbiasHelixPairPotential ) ),
	hspot_( NatbiasHelicesSheetPotentialOP( new NatbiasHelicesSheetPotential ) )
{
	set_blueprint( filename );
}

/// @brief value constructor
SetSecStructEnergies::SetSecStructEnergies( ScoreFunctionOP const sfx, BluePrintOP const blueprintOP, bool const ss_from_blueprint ) :
	Super( "SetSecStructEnergies" ),
	loaded_( false ),
	blueprint_( blueprintOP ),
	ss_from_blueprint_( ss_from_blueprint ),
	sfx_( sfx ),
	hh_weight_( 1.0 ),
	hs_weight_( 1.0 ),
	ss_weight_( 1.0 ),
	stwist_weight_( 1.0 ),
	hs_pair_weight_( 0.0 ),
	ss_pair_weight_( 0.0 ),
	rsigma_weight_( 0.0 ),
	hpairpot_( NatbiasHelixPairPotentialOP( new NatbiasHelixPairPotential ) ),
	hspot_( NatbiasHelicesSheetPotentialOP( new NatbiasHelicesSheetPotential ) )
{}

/// @Brief copy constructor
SetSecStructEnergies::SetSecStructEnergies( SetSecStructEnergies const & rval ) :
	Super( rval ),
	loaded_( rval.loaded_ ),
	blueprint_( rval.blueprint_ ),
	ss_from_blueprint_( rval.ss_from_blueprint_ ),
	sfx_( rval.sfx_ ),
	hh_weight_( rval.hh_weight_ ),
	hs_weight_( rval.hs_weight_ ),
	ss_weight_( rval.ss_weight_ ),
	stwist_weight_( rval.stwist_weight_ ),
	hs_pair_weight_( rval.hs_pair_weight_ ),
	ss_pair_weight_( rval.ss_pair_weight_ ),
	rsigma_weight_( rval.rsigma_weight_ ),
	hpairpot_( rval.hpairpot_ ),
	hspot_( rval.hspot_ )
{}

/// @brief default destructor
SetSecStructEnergies::~SetSecStructEnergies() {}

/// @brief clone this object
SetSecStructEnergies::MoverOP
SetSecStructEnergies::clone() const
{
	return SetSecStructEnergies::MoverOP( new SetSecStructEnergies( *this ) );
}

/// @brief create this type of object
SetSecStructEnergies::MoverOP
SetSecStructEnergies::fresh_instance() const
{
	return SetSecStructEnergies::MoverOP( new SetSecStructEnergies() );
}

/// @brief set the centroid level score function
void
SetSecStructEnergies::scorefunction( ScoreFunction const & sfx )
{
	sfx_ = sfx.clone();
}

/// @brief set the centroid level score function
void
SetSecStructEnergies::scorefunction( ScoreFunctionOP sfx )
{
	sfx_ = sfx->clone();
}

/// @brief use blueprint
void
SetSecStructEnergies::set_blueprint( String const & filename )
{
	blueprint_ = BluePrintOP( new BluePrint( filename ) );
}

/// @brief use blueprint
void
SetSecStructEnergies::set_blueprint( BluePrintOP const blp )
{
	blueprint_ = blp;
}

/// @brief make symmetric secstruct
SetSecStructEnergies::String
SetSecStructEnergies::symmetric_secstruct( SymmetryInfoOP const syminfo, String const & ss )
{
	runtime_assert( ss.length() == syminfo->num_independent_residues() );

	String secstruct("");
	for( Size i=1; i<= syminfo->subunits(); i++ ) {
		secstruct += ss;
	}
	for( Size i=1; i<=syminfo->num_virtuals(); i++ ) {
		secstruct += "L";
	}

	return secstruct;
}

/// @brief apply defined moves to given Pose
void SetSecStructEnergies::apply( Pose & pose )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using core::scoring::dssp::Dssp;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;

	using core::scoring::ScoreType;
	using protocols::fldsgn::potentials::sspot::NatbiasSecondaryStructureEnergy;

	using protocols::fldsgn::potentials::sspot::NatbiasStrandPairPotential;
	using protocols::fldsgn::potentials::sspot::NatbiasStrandPairPotentialOP;

	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::StrandPairingSet;
	using protocols::fldsgn::topology::StrandPairingSetOP;
	using protocols::fldsgn::topology::HelixPairingSet;
	using protocols::fldsgn::topology::HelixPairingSetOP;
	using protocols::fldsgn::topology::HSSTripletSet;
	using protocols::fldsgn::topology::HSSTripletSetOP;

	using core::conformation::symmetry::SymmetricConformation;
	using core::pose::symmetry::is_symmetric;
	using protocols::simple_moves::symmetry::SetupForSymmetryMover;
	using protocols::simple_moves::symmetry::SetupForSymmetryMoverOP;

	if( loaded_ ) return;

	runtime_assert( blueprint_ != 0 );
	runtime_assert( sfx_ != 0 );

	// assign secondary structure
	if( !ss_from_blueprint_ ){
		Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );
	}else{
		blueprint_->insert_ss_into_pose( pose );
	}

	Pose scratch( pose );
	String ss("");
	if( option[ basic::options::OptionKeys::symmetry::symmetry_definition ].active() ) {

		if ( !is_symmetric( scratch ) ) {
			SetupForSymmetryMoverOP symm_setup_mover( new SetupForSymmetryMover );
			symm_setup_mover->apply( scratch );
		}

		SymmetricConformation & symm_conf ( dynamic_cast<SymmetricConformation & > ( scratch.conformation()) );
		SymmetryInfoOP syminfo( symm_conf.Symmetry_Info() );
		ss = symmetric_secstruct( syminfo, blueprint_->secstruct() );

	} else {
		ss = blueprint_->secstruct();
	}

	// set NatbiasSecondaryStructure energy
	SS_Info2_OP ssinfo( new SS_Info2( ss ) );
	StrandPairingSetOP spairset( new StrandPairingSet( blueprint_->strand_pairings(), ssinfo ) );
	HelixPairingSetOP  hpairset( new HelixPairingSet ( blueprint_->helix_pairings() ) );
	HSSTripletSetOP     hss3set( new HSSTripletSet   ( blueprint_->hss_triplets() ) );

	NatbiasSecondaryStructureEnergy sspot;
	NatbiasStrandPairPotentialOP   spairpot( new NatbiasStrandPairPotential( spairset ) );
	hpairpot_->set_hpairset( hpairset );
	hpairpot_->show_params();

	hspot_->hss_triplet_set( hss3set );
	hspot_->show_params();

	sspot.native_secstruct( ss );
	if( ss_weight_ != 0.0 ) sspot.set_natbias_spairpot( spairpot );
	if( hs_weight_ != 0.0 ) sspot.set_natbias_helices_sheet_pot( hspot_ );
	if( hh_weight_ != 0.0 ) sspot.set_natbias_hpairpot( hpairpot_ );

	// turn Original SecondaryStructureEnergy off
	std::map< ScoreType, Real > new_weights;
	new_weights.insert( std::map< ScoreType, Real >::value_type( core::scoring::natbias_hs, hs_weight_ ) );
	new_weights.insert( std::map< ScoreType, Real >::value_type( core::scoring::natbias_ss, ss_weight_ ) );
	new_weights.insert( std::map< ScoreType, Real >::value_type( core::scoring::natbias_hh, hh_weight_ ) );
	new_weights.insert( std::map< ScoreType, Real >::value_type( core::scoring::natbias_stwist, stwist_weight_ ) );
	sfx_->add_extra_method( new_weights, sspot );

	sfx_->set_weight( core::scoring::hs_pair, hs_pair_weight_ );
	if( hs_pair_weight_ <= 0.0 ) {
		TR << "original hs_pair term was turned off." << std::endl;
	}

	sfx_->set_weight( core::scoring::ss_pair, ss_pair_weight_ );
	if( ss_pair_weight_ <= 0.0 ) {
		TR << "original ss_pair term was turned off." << std::endl;
	}

	sfx_->set_weight( core::scoring::rsigma, rsigma_weight_ );
	if( rsigma_weight_ <= 0.0 ) {
		TR << "original rsigma term was turned off." << std::endl;
	}

	// this mover have to be called only once
	loaded_ = true;

}

/// @brief
std::string
SetSecStructEnergies::get_name() const {
	return SetSecStructEnergiesCreator::mover_name();
}


/// @brief parse xml
void
SetSecStructEnergies::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	std::string const blueprint( tag->getOption<std::string>( "blueprint", "" ) );
	if( blueprint == "" ){
		TR << "No input of blueprint file ! " << std::endl;
		runtime_assert( false );
	}
	set_blueprint( blueprint );

	bool const ss_from_blueprint( tag->getOption<bool>( "ss_from_blueprint", true ) );
	ss_from_blueprint_ = ss_from_blueprint;

	// set scorefxn
	String const sfxn ( tag->getOption<String>( "scorefxn", "" ) );
	if( sfxn == "" ) {
		TR << "[ERROR] No input of scorefxn !"  << std::endl;
		runtime_assert( false );
	}

	sfx_ = data.get_ptr<ScoreFunction>( "scorefxns", sfxn );
	TR << "score function, " << sfxn << ", is used. " << std::endl;

	///
	if( tag->hasOption( "natbias_ss" ) ) ss_weight_ = tag->getOption< Real >( "natbias_ss" );
	if( tag->hasOption( "natbias_hh" ) ) hh_weight_ = tag->getOption< Real >( "natbias_hh" );
	if( tag->hasOption( "natbias_hs" ) ) hs_weight_ = tag->getOption< Real >( "natbias_hs" );
	if( tag->hasOption( "natbias_stwist" ) ) stwist_weight_ = tag->getOption< Real >( "natbias_stwist" );

	/// original secondary structure potential except for sheet potential
	if( tag->hasOption( "hs_pair" ) ) hs_pair_weight_ = tag->getOption< Real >( "hs_pair" );
	if( tag->hasOption( "ss_pair" ) ) ss_pair_weight_ = tag->getOption< Real >( "ss_pair" );
	if( tag->hasOption( "rsigma"  ) ) rsigma_weight_  = tag->getOption< Real >( "rsigma" );

	// params for NatbiasHelixPairPotential
	if( tag->hasOption( "hh_dist_wts" ) ) {
		hpairpot_->set_dist_wts( tag->getOption< Real >( "hh_dist_wts" ) );
	}
	if( tag->hasOption( "hh_dist" ) ) {
		hpairpot_->set_dist( tag->getOption< Real >( "hh_dist" ) );
	}
	if( tag->hasOption( "hh_dist_s2" ) ) {
		hpairpot_->set_dist_sigma2( tag->getOption< Real >( "hh_dist_s2" ) );
	}

	if( tag->hasOption( "hh_cross_angle_wts" ) ) {
		hpairpot_->set_angle_wts( tag->getOption< Real >( "hh_cross_angle_wts" ) );
	}
	if( tag->hasOption( "hh_cross_angle" ) ) {
		hpairpot_->set_angle( tag->getOption< Real >( "hh_cross_angle" ) );
	}
	if( tag->hasOption( "hh_cross_angle_s2" ) ) {
		hpairpot_->set_angle_sigma2( tag->getOption< Real >( "hh_cross_angle_s2" ) );
	}

	// params for NatbiasHelicesSheetPotential
	if( tag->hasOption( "hs_atr_dist_wts" ) )	{
		hspot_->set_hs_atr_dist_wts( tag->getOption< Real >( "hs_atr_dist_wts" ) );
	}
	if( tag->hasOption( "hs_atr_dist" ) ) {
		hspot_->set_hs_atr_dist( tag->getOption< Real >( "hs_atr_dist" ) );
	}
	if( tag->hasOption( "hs_atr_dist_s2" ) ) {
		hspot_->set_hs_atr_dist_sigma2( tag->getOption< Real >( "hs_atr_dist_s2" ) );
	}

	if( tag->hasOption( "hs_angle_wts" ) ) {
		hspot_->set_hs_angle_wts( tag->getOption< Real >( "hs_angle_wts" ) );
	}
	if( tag->hasOption( "hs_angle" ) ) {
		hspot_->set_hs_angle( tag->getOption< Real >( "hs_angle" ) );
	}
	if( tag->hasOption( "hs_angle_s2") ) {
		hspot_->set_hs_angle_sigma2( tag->getOption< Real >( "hs_angle_s2" ) );
	}

	if( tag->hasOption( "hsheet_repl_dist" ) ) {
		hspot_->set_hsheet_repl_dist( tag->getOption< Real >( "hsheet_repl_dist" ) );
	}

	if( tag->hasOption( "hh_align_angle_wts" ) ) {
		hspot_->set_hh_angle_wts( tag->getOption< Real >( "hh_align_angle_wts" ) );
	}
	if( tag->hasOption( "hh_align_angle" ) ) {
		hspot_->set_hh_angle( tag->getOption< Real >( "hh_align_angle" ) );
	}
	if( tag->hasOption( "hh_align_angle_s2" ) ) {
		hspot_->set_hh_angle_sigma2( tag->getOption< Real >( "hh_align_angle_s2" ) );
	}

}

} // Namespace potentials
} // namespace fldsgn
} // namespace protocols

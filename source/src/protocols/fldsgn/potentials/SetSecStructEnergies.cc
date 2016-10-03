// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

#include <protocols/fldsgn/potentials/sspot/NatbiasStrandPairPotential.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.SetSecStructEnergies" );

namespace protocols {
namespace fldsgn {
namespace potentials {

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
	secstruct_( "" ),
	use_dssp_( true ),
	hh_pair_( "" ),
	ss_pair_( "" ),
	hss_triplet_( "" ),
	sfx_orig_(),
	sfx_( /* NULL */ ),
	hh_weight_( 1.0 ),
	hs_weight_( 1.0 ),
	ss_weight_( 1.0 ),
	stwist_weight_( 1.0 ),
	hs_pair_weight_( 0.0 ),
	ss_pair_weight_( 0.0 ),
	rsigma_weight_( 0.0 ),
	add_symmetry_( false ),
	hpairpot_( NatbiasHelixPairPotentialOP( new NatbiasHelixPairPotential ) ),
	hspot_( NatbiasHelicesSheetPotentialOP( new NatbiasHelicesSheetPotential ) )
{}

/// @brief value constructor
SetSecStructEnergies::SetSecStructEnergies( ScoreFunctionOP const sfx, String const & filename, bool const ss_from_blueprint ) :
	Super( "SetSecStructEnergies" ),
	secstruct_( "" ),
	use_dssp_( true ),
	hh_pair_( "" ),
	ss_pair_( "" ),
	hss_triplet_( "" ),
	sfx_orig_( sfx->clone() ),
	sfx_( sfx ),
	hh_weight_( 1.0 ),
	hs_weight_( 1.0 ),
	ss_weight_( 1.0 ),
	stwist_weight_( 1.0 ),
	hs_pair_weight_( 0.0 ),
	ss_pair_weight_( 0.0 ),
	rsigma_weight_( 0.0 ),
	add_symmetry_( false ),
	hpairpot_( NatbiasHelixPairPotentialOP( new NatbiasHelixPairPotential ) ),
	hspot_( NatbiasHelicesSheetPotentialOP( new NatbiasHelicesSheetPotential ) )
{
	protocols::jd2::parser::BluePrint bp( filename );
	init_from_blueprint( bp, ss_from_blueprint );
}

/// @brief value constructor
SetSecStructEnergies::SetSecStructEnergies( ScoreFunctionOP const sfx, BluePrintOP const blueprintOP, bool const ss_from_blueprint ) :
	Super( "SetSecStructEnergies" ),
	secstruct_( "" ),
	use_dssp_( true ),
	hh_pair_( "" ),
	ss_pair_( "" ),
	hss_triplet_( "" ),
	sfx_orig_( sfx->clone() ),
	sfx_( sfx ),
	hh_weight_( 1.0 ),
	hs_weight_( 1.0 ),
	ss_weight_( 1.0 ),
	stwist_weight_( 1.0 ),
	hs_pair_weight_( 0.0 ),
	ss_pair_weight_( 0.0 ),
	rsigma_weight_( 0.0 ),
	add_symmetry_( false ),
	hpairpot_( NatbiasHelixPairPotentialOP( new NatbiasHelixPairPotential ) ),
	hspot_( NatbiasHelicesSheetPotentialOP( new NatbiasHelicesSheetPotential ) )
{
	init_from_blueprint( *blueprintOP, ss_from_blueprint );
}

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

/// @brief access ptr to the modified score function
SetSecStructEnergies::ScoreFunctionOP
SetSecStructEnergies::scorefunction_ptr() const
{
	return sfx_;
}

/// @brief set the centroid level score function
void
SetSecStructEnergies::set_scorefunction_ptr( ScoreFunctionOP sfx )
{
	/// TL: we don't want to clone here. The given scorefunction is modified and we assume
	///     they have access to the object
	sfx_ = sfx;
	if ( sfx ) {
		sfx_orig_ = sfx->clone();
	} else {
		sfx_orig_ = ScoreFunctionOP();
	}
}

/// @brief sets the secondary structure to be used for computation
/// @param[in] secstruct Secondary structure to be forced on the pose
/// @details If this is non-empty, it will be used as the pose secondary structure to
///          determine secondary structure elements in the pose
void
SetSecStructEnergies::set_secstruct( String const & secstruct )
{
	secstruct_ = secstruct;
}

/// @brief make symmetric secstruct
SetSecStructEnergies::String
SetSecStructEnergies::symmetric_secstruct( SymmetryInfoCOP const syminfo, String const & ss ) const
{
	runtime_assert( ss.length() == syminfo->num_independent_residues() );

	String secstruct("");
	for ( Size i=1; i<= syminfo->subunits(); i++ ) {
		secstruct += ss;
	}
	for ( Size i=1; i<=syminfo->num_virtuals(); i++ ) {
		secstruct += "L";
	}

	return secstruct;
}

/// @brief apply defined moves to given Pose
void SetSecStructEnergies::apply( Pose & pose )
{
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

	// TL: Movers should be repeatedly callable for different poses.
	//     In this case, I might be building a different structure, in which case
	//     the potentials need to be re-created.
	//if ( loaded_ ) return;

	runtime_assert( sfx_ != 0 );
	runtime_assert( sfx_orig_ );

	// reset to original so this can be called repeatedly
	sfx_->reset();
	sfx_->assign( *sfx_orig_ );

	std::string const ss = get_secstruct( pose );

	TR.Debug << "Using secstruct " << ss << std::endl;
	TR.Debug << "Using strand pairs " << ss_pair_ << std::endl;
	TR.Debug << "Using helix pairs " << hh_pair_ << std::endl;
	TR.Debug << "Using hss-triplets" << hss_triplet_ << std::endl;

	// set NatbiasSecondaryStructure energy
	SS_Info2_OP ssinfo( new SS_Info2( ss ) );

	StrandPairingSetOP spairset( new StrandPairingSet( ss_pair_, ssinfo ) );
	spairset->finalize();

	HelixPairingSetOP  hpairset( new HelixPairingSet ( hh_pair_ ) );
	HSSTripletSetOP     hss3set( new HSSTripletSet   ( hss_triplet_ ) );

	NatbiasSecondaryStructureEnergy sspot;
	NatbiasStrandPairPotentialOP   spairpot( new NatbiasStrandPairPotential( spairset ) );
	hpairpot_->set_hpairset( hpairset );
	hpairpot_->show_params();

	hspot_->hss_triplet_set( hss3set );
	hspot_->show_params();

	sspot.native_secstruct( ss );
	if ( ss_weight_ != 0.0 ) sspot.set_natbias_spairpot( spairpot );
	if ( hs_weight_ != 0.0 ) sspot.set_natbias_helices_sheet_pot( hspot_ );
	if ( hh_weight_ != 0.0 ) sspot.set_natbias_hpairpot( hpairpot_ );

	// turn Original SecondaryStructureEnergy off
	std::map< ScoreType, Real > new_weights;
	new_weights.insert( std::map< ScoreType, Real >::value_type( core::scoring::natbias_hs, hs_weight_ ) );
	new_weights.insert( std::map< ScoreType, Real >::value_type( core::scoring::natbias_ss, ss_weight_ ) );
	new_weights.insert( std::map< ScoreType, Real >::value_type( core::scoring::natbias_hh, hh_weight_ ) );
	new_weights.insert( std::map< ScoreType, Real >::value_type( core::scoring::natbias_stwist, stwist_weight_ ) );
	sfx_->add_extra_method( new_weights, sspot );

	sfx_->set_weight( core::scoring::hs_pair, hs_pair_weight_ );
	if ( hs_pair_weight_ <= 0.0 ) {
		TR << "original hs_pair term was turned off." << std::endl;
	}

	sfx_->set_weight( core::scoring::ss_pair, ss_pair_weight_ );
	if ( ss_pair_weight_ <= 0.0 ) {
		TR << "original ss_pair term was turned off." << std::endl;
	}

	sfx_->set_weight( core::scoring::rsigma, rsigma_weight_ );
	if ( rsigma_weight_ <= 0.0 ) {
		TR << "original rsigma term was turned off." << std::endl;
	}

	// this mover have to be called only once
	//loaded_ = true;

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
	std::string const blueprintfile( tag->getOption<std::string>( "blueprint", "" ) );
	if ( !blueprintfile.empty() ) {
		protocols::jd2::parser::BluePrint const bp( blueprintfile );
		bool const ss_from_blueprint( tag->getOption<bool>( "ss_from_blueprint", true ) );
		init_from_blueprint( bp, ss_from_blueprint );
	}

	secstruct_ = tag->getOption< std::string >( "secstruct", secstruct_ );
	use_dssp_ = tag->getOption< bool >( "use_dssp", use_dssp_ );
	hh_pair_ = tag->getOption< std::string >( "hh_pair", hh_pair_ );
	ss_pair_ = tag->getOption< std::string >( "ss_pair", ss_pair_ );
	hss_triplet_ = tag->getOption< std::string >( "hss_triplets", hss_triplet_ );

	// set scorefxn
	String const sfxn ( tag->getOption<String>( "scorefxn", "" ) );
	if ( sfxn == "" ) {
		TR << "[ERROR] No input of scorefxn !"  << std::endl;
		runtime_assert( false );
	}

	add_symmetry_ = tag->getOption< bool >( "add_symmetry", add_symmetry_ );
	if ( !tag->hasOption( "add_symmetry" ) ) {
		add_symmetry_ =
			basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].active();
	}

	set_scorefunction_ptr( data.get_ptr< ScoreFunction >( "scorefxns", sfxn ) );
	TR << "score function, " << sfxn << ", is used. " << std::endl;


	if ( tag->hasOption( "natbias_ss" ) ) ss_weight_ = tag->getOption< Real >( "natbias_ss" );
	if ( tag->hasOption( "natbias_hh" ) ) hh_weight_ = tag->getOption< Real >( "natbias_hh" );
	if ( tag->hasOption( "natbias_hs" ) ) hs_weight_ = tag->getOption< Real >( "natbias_hs" );
	if ( tag->hasOption( "natbias_stwist" ) ) stwist_weight_ = tag->getOption< Real >( "natbias_stwist" );

	/// original secondary structure potential except for sheet potential
	if ( tag->hasOption( "hs_pair" ) ) hs_pair_weight_ = tag->getOption< Real >( "hs_pair" );
	if ( tag->hasOption( "ss_pair" ) ) ss_pair_weight_ = tag->getOption< Real >( "ss_pair" );
	if ( tag->hasOption( "rsigma"  ) ) rsigma_weight_  = tag->getOption< Real >( "rsigma" );

	// params for NatbiasHelixPairPotential
	if ( tag->hasOption( "hh_dist_wts" ) ) {
		hpairpot_->set_dist_wts( tag->getOption< Real >( "hh_dist_wts" ) );
	}
	if ( tag->hasOption( "hh_dist" ) ) {
		hpairpot_->set_dist( tag->getOption< Real >( "hh_dist" ) );
	}
	if ( tag->hasOption( "hh_dist_s2" ) ) {
		hpairpot_->set_dist_sigma2( tag->getOption< Real >( "hh_dist_s2" ) );
	}

	if ( tag->hasOption( "hh_cross_angle_wts" ) ) {
		hpairpot_->set_angle_wts( tag->getOption< Real >( "hh_cross_angle_wts" ) );
	}
	if ( tag->hasOption( "hh_cross_angle" ) ) {
		hpairpot_->set_angle( tag->getOption< Real >( "hh_cross_angle" ) );
	}
	if ( tag->hasOption( "hh_cross_angle_s2" ) ) {
		hpairpot_->set_angle_sigma2( tag->getOption< Real >( "hh_cross_angle_s2" ) );
	}

	// params for NatbiasHelicesSheetPotential
	if ( tag->hasOption( "hs_atr_dist_wts" ) ) {
		hspot_->set_hs_atr_dist_wts( tag->getOption< Real >( "hs_atr_dist_wts" ) );
	}
	if ( tag->hasOption( "hs_atr_dist" ) ) {
		hspot_->set_hs_atr_dist( tag->getOption< Real >( "hs_atr_dist" ) );
	}
	if ( tag->hasOption( "hs_atr_dist_s2" ) ) {
		hspot_->set_hs_atr_dist_sigma2( tag->getOption< Real >( "hs_atr_dist_s2" ) );
	}

	if ( tag->hasOption( "hs_angle_wts" ) ) {
		hspot_->set_hs_angle_wts( tag->getOption< Real >( "hs_angle_wts" ) );
	}
	if ( tag->hasOption( "hs_angle" ) ) {
		hspot_->set_hs_angle( tag->getOption< Real >( "hs_angle" ) );
	}
	if ( tag->hasOption( "hs_angle_s2") ) {
		hspot_->set_hs_angle_sigma2( tag->getOption< Real >( "hs_angle_s2" ) );
	}

	if ( tag->hasOption( "hsheet_repl_dist" ) ) {
		hspot_->set_hsheet_repl_dist( tag->getOption< Real >( "hsheet_repl_dist" ) );
	}

	if ( tag->hasOption( "hh_align_angle_wts" ) ) {
		hspot_->set_hh_angle_wts( tag->getOption< Real >( "hh_align_angle_wts" ) );
	}
	if ( tag->hasOption( "hh_align_angle" ) ) {
		hspot_->set_hh_angle( tag->getOption< Real >( "hh_align_angle" ) );
	}
	if ( tag->hasOption( "hh_align_angle_s2" ) ) {
		hspot_->set_hh_angle_sigma2( tag->getOption< Real >( "hh_align_angle_s2" ) );
	}

}

/// @brief chooses and return the secondary structure to be used in computation
/// @param[in] pose Input pose
/// @details  If secstruct_ is set, returns that.
///           If use_dssp_ is set, returns DSSP secstruct string
///           Otherwise, returns pose secondary structure
std::string
SetSecStructEnergies::get_secstruct( core::pose::Pose const & pose ) const
{
	using core::conformation::symmetry::SymmetricConformation;
	using core::pose::symmetry::is_symmetric;
	using protocols::simple_moves::symmetry::SetupForSymmetryMover;
	using protocols::simple_moves::symmetry::SetupForSymmetryMoverOP;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string secstruct;
	if ( !secstruct_.empty() ) {
		secstruct = secstruct_;
	} else if ( use_dssp_ ) {
		core::scoring::dssp::Dssp dssp( pose );
		secstruct = dssp.get_dssp_secstruct();
	} else {
		secstruct = pose.secstruct();
	}

	if ( add_symmetry_ ) {
		Pose scratch( pose );
		if ( !is_symmetric( scratch ) ) {
			SetupForSymmetryMoverOP symm_setup_mover( new SetupForSymmetryMover );
			symm_setup_mover->apply( scratch );
		}

		SymmetricConformation & symm_conf ( dynamic_cast<SymmetricConformation & > ( scratch.conformation()) );
		SymmetryInfoOP syminfo( symm_conf.Symmetry_Info() );
		return symmetric_secstruct( syminfo, secstruct );
	}

	return secstruct;
}

/// @brief initializes pairings and secondary structure from blueprint file
/// @param[in] bp                BluePrint object
/// @param[in] ss_from_blueprint If true, secstruct_ will be set from the blueprint.  If false,
///                              secstruct_ will be unchanged, and only the pairings will be set.
void
SetSecStructEnergies::init_from_blueprint( protocols::jd2::parser::BluePrint const & bp, bool const ss_from_blueprint )
{
	if ( ss_from_blueprint ) secstruct_ = bp.secstruct();
	hh_pair_ = bp.helix_pairings();
	ss_pair_ = bp.strand_pairings();
	hss_triplet_ = bp.hss_triplets();
}


} // Namespace potentials
} // namespace fldsgn
} // namespace protocols

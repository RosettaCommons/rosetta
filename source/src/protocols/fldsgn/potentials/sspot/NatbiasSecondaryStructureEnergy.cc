// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/potentials/sspot/NatbiasSecondaryStructureEnergy.cc
/// @brief native biased centroid score for secondary structures
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit headers
#include <protocols/fldsgn/potentials/sspot/NatbiasSecondaryStructureEnergy.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasSecondaryStructureEnergyCreator.hh>

// Package headers
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasHelicesSheetPotential.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasStrandPairPotential.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethod.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.fldsgn.potentials.sspot.NatbiasSecondaryStructureEnergy", basic::t_info );

namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

//////////////////////////////////////////////////////////////////////////////////////////////////////
/// SecodaryStructureEnergy2Creator
/// @details This must return a fresh instance of the SecondaryStructureEnergy class,
/// never an instance already in use
NatbiasSecondaryStructureEnergyCreator::EnergyMethodOP
NatbiasSecondaryStructureEnergyCreator::create_energy_method(	EnergyMethodOptions const & ) const
{
	return NatbiasSecondaryStructureEnergyCreator::EnergyMethodOP( new NatbiasSecondaryStructureEnergy );
}

NatbiasSecondaryStructureEnergyCreator::ScoreTypes
NatbiasSecondaryStructureEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( core::scoring::natbias_hs );
	sts.push_back( core::scoring::natbias_ss );
	sts.push_back( core::scoring::natbias_hh );
	sts.push_back( core::scoring::natbias_stwist );
	return sts;
}


/// @brief default constructor
NatbiasSecondaryStructureEnergy::NatbiasSecondaryStructureEnergy() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new NatbiasSecondaryStructureEnergyCreator ) ),
	native_secstruct_( "" ),
	use_sspot_( false ),
	use_hhpot_( false ),
	use_hspot_( false ),
	use_nobias_( false ),
	sspot_( /* NULL */ ),
	hhpot_( /* NULL */ ),
	hspot_( /* NULL */ )
{}


/// @brief copy constructor
NatbiasSecondaryStructureEnergy::NatbiasSecondaryStructureEnergy( NatbiasSecondaryStructureEnergy const & src ):
	parent( src ),
	native_secstruct_( src.native_secstruct_ ),
	use_sspot_( src.use_sspot_ ),
	use_hhpot_( src.use_hhpot_ ),
	use_hspot_( src.use_hspot_ ),
	use_nobias_( src.use_nobias_ ),
	sspot_( src.sspot_ ),
	hhpot_( src.hhpot_ ),
	hspot_( src.hspot_ )
{}


/// @brief clone
NatbiasSecondaryStructureEnergy::EnergyMethodOP
NatbiasSecondaryStructureEnergy::clone() const
{
	return NatbiasSecondaryStructureEnergy::EnergyMethodOP( new NatbiasSecondaryStructureEnergy( *this ) );
}

/// @brief set native secondary structure
void
NatbiasSecondaryStructureEnergy::native_secstruct( String const & secstruct )
{
	native_secstruct_ = secstruct;
}

/// @brief set native NatbiasStrandPairPotential
void
NatbiasSecondaryStructureEnergy::set_natbias_spairpot( StrandPairingSetOP const spairset )
{
	using protocols::fldsgn::potentials::sspot::NatbiasStrandPairPotential;
	use_sspot_ = true;
	sspot_ = NatbiasStrandPairPotentialOP( new NatbiasStrandPairPotential( spairset ) );
}


/// @brief set NatbiasHelixPairPotential
void
NatbiasSecondaryStructureEnergy::set_natbias_hpairpot( HelixPairingSetOP const hpairset )
{
	using protocols::fldsgn::potentials::sspot::NatbiasHelixPairPotential;
	use_hhpot_ = true;
	hhpot_ = NatbiasHelixPairPotentialOP( new NatbiasHelixPairPotential( hpairset ) );
}

/// @brief set HelicesSheetPotential
void
NatbiasSecondaryStructureEnergy::set_natbias_helices_sheet_pot( HSSTripletSetOP const hss3set )
{
	using protocols::fldsgn::potentials::sspot::NatbiasHelicesSheetPotential;
	use_hspot_ = true;
	hspot_ = NatbiasHelicesSheetPotentialOP( new NatbiasHelicesSheetPotential( hss3set ) );
}


/// @brief set native NatbiasStrandPairPotential
void
NatbiasSecondaryStructureEnergy::set_natbias_spairpot( NatbiasStrandPairPotentialOP const sspot )
{
	use_sspot_ = true;
	sspot_ = sspot;
}

/// @brief set NatbiasHelixPairPotential
void
NatbiasSecondaryStructureEnergy::set_natbias_hpairpot( NatbiasHelixPairPotentialOP const hhpot )
{
	use_hhpot_ = true;
	hhpot_ = hhpot;
}

/// @brief set NatbiasHelicesSheetPotential
void
NatbiasSecondaryStructureEnergy::set_natbias_helices_sheet_pot( NatbiasHelicesSheetPotentialOP const hspot )
{
	use_hspot_ = true;
	hspot_ = hspot;
}

/// @brief set up for scoring
void
NatbiasSecondaryStructureEnergy::setup_for_scoring( Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

/// @brief all scoring happens here
void
NatbiasSecondaryStructureEnergy::finalize_total_energy(
	Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;

	// this runtime_assert causes rosetta to abort if there is a ligand in the pose
	// However, a ligand-aware version of it occurs in SS_Info2, so it is not necessary to do it here
	//runtime_assert( pose.total_residue() == native_secstruct_.size() );

	Real ss_score( 0.0 ), hh_score( 0.0 ), hs_score( 0.0 );

	SS_Info2_OP ssinfo( new SS_Info2( pose, native_secstruct_ ) );

	if( use_sspot_ ) sspot_->score( pose, *ssinfo, ss_score );

	if( use_hhpot_ ) hhpot_->score( ssinfo, hh_score );

	if( use_hspot_ ) hspot_->score( ssinfo, hh_score, hs_score );

	// set calculated score
	totals[ core::scoring::natbias_hs ] = hs_score;
	totals[ core::scoring::natbias_ss ] = ss_score;
	totals[ core::scoring::natbias_hh ] = hh_score;

	// TR << "hs_score: " << hs_score << " hh_score: " << hh_score << " ss_score: " << ss_score << std::endl;

}


/// @brief SecondaryStructureEnergy distance cutoff
core::Distance
NatbiasSecondaryStructureEnergy::atomic_interaction_cutoff() const
{
	return 6.0; /// now subtracted off 6.0 from cutoffs in centroid params files
// 	return 0.0; /// since all the cutoffs for centroid mode are rolled into the cendist check
}

/// @brief SecondaryStructureEnergy
void
NatbiasSecondaryStructureEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{}

core::Size
NatbiasSecondaryStructureEnergy::version() const
{
	return 1;
}

} // sspot
} // potentials
}	// fldsgn
}	// protocols

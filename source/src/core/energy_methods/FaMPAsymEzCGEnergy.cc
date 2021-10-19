// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/FaMPAsymEzCG.hh
///
/// @brief  Fullatom implementation of asymetric EZ potential from
/// @details See Schramm et al 2012 at 10.1016/j.str.2012.03.016 for specific details.
///    Implemented in bins of 1A along the Z-axis; depth of 0 is middle of membrane. Positions
///    more than 30A from center are assigned a score of 0. Assigned scores are based on residue identity
///    and on Z-coordinate of the CG atom.
///    Last Modified: 7/3/18
///
/// @author  Meghan Franklin (meghanwfranklin@gmail.com)

// Unit headers
#include <core/energy_methods/FaMPAsymEzCGEnergy.hh>
#include <core/energy_methods/FaMPAsymEzCGEnergyCreator.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/conformation/Atom.hh>
#include <core/id/AtomID.hh>

#include <core/scoring/Energies.hh>


#include <core/scoring/membrane/MembraneData.hh>

#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

static basic::Tracer TR( "core.energy_methods.FaMPAsymEzCGEnergy" );


using namespace core::scoring::methods;

namespace core {
namespace energy_methods {

// Membrane distance limits
static const core::Real MAX_Z_POSITION = 30.5;
static const core::Real Z_BIN_SHIFT = 31;
// e table dimensions
static const core::Real MAX_AA = 20;
static const core::Real ASYMEZ_TABLE_BINS = 61;

// Creator Methods //////////////////////////////////////////

/// @brief Return a fresh instance of the energy method
core::scoring::methods::EnergyMethodOP
FaMPAsymEzCGEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< FaMPAsymEzCGEnergy >();
}

/// @brief Return relevant score types
scoring::ScoreTypes
FaMPAsymEzCGEnergyCreator::score_types_for_method() const {
	scoring::ScoreTypes sts;
	sts.push_back( scoring::FaMPAsymEzCG );
	return sts;
}

// Constructors /////////////////////////////////////////////

FaMPAsymEzCGEnergy::FaMPAsymEzCGEnergy() :
	parent( utility::pointer::make_shared< FaMPAsymEzCGEnergyCreator >() )
{
	Size const max_aa( MAX_AA );
	Size const asymEZ_table_bins( ASYMEZ_TABLE_BINS );

	std::string tag,line;
	chemical::AA aa;
	asymEZ_CG_.dimension( max_aa, asymEZ_table_bins);

	// Read in updated table
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/MembranePotential/AsymEZ_CG.txt" );

	Size i=1;
	while ( stream  ) {
		getline( stream, line );
		if ( line[0] == '#' ) continue;
		std::istringstream l(line);
		l >> aa;
		for ( Size j=1; j<=asymEZ_table_bins; ++j ) {
			if ( l.fail() ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/AsymEZ_CG.txt (FaMPAsymEzCGEnergy)");
			l >> asymEZ_CG_(aa,j);
		}
		if ( ++i > max_aa ) break;
	}
}

/// @brief Create a clone of this energy method
scoring::methods::EnergyMethodOP
FaMPAsymEzCGEnergy::clone() const
{
	return utility::pointer::make_shared< FaMPAsymEzCGEnergy >( *this );
}

// Scoring Methods ////////////////////////////////////////////////

/// @details looks up score for depth of each CG in membrane.
void
FaMPAsymEzCGEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	scoring::EnergyMap & emap
) const
{
	// asymEZ was only developed from/for proteins
	if ( ! ( rsd.is_protein() &&  rsd.type().is_canonical_aa() ) ) return;


	Size const atomindex_i = rsd.atom_index( representative_atom_name( rsd.aa() ));

	Real score = 0;

	Real const z_position( pose.conformation().membrane_info()->atom_z_position( pose.conformation(), rsd.seqpos(), atomindex_i ) );


	if ( (z_position <= (-1*MAX_Z_POSITION)) || (z_position >= MAX_Z_POSITION) ) {
		score = 0; //outside the membrane
	} else {
		//z-bin is the rounded z-coord + 30 (because Rosetta is indexed from 1, not -30)
		Size z_bin = static_cast< Size > ( round(z_position) + Z_BIN_SHIFT);
		score = asymEZ_CG_(rsd.aa(), z_bin);
		//std::cout << rsd.aa() << rsd.seqpos() << "\tz = " << z_position << "\tz_bin = " << z_bin << "\tFaMPAsymEzCB = " << score << std::endl;
	}

	if ( rsd.aa() == core::chemical::aa_ile || rsd.aa() == core::chemical::aa_val ) {
		Size const atomindex_i = rsd.atom_index( "CG2" );
		Real const z_position( pose.conformation().membrane_info()->atom_z_position( pose.conformation(), rsd.seqpos(), atomindex_i ) );
		Size z_bin = static_cast< Size > ( round(z_position) + Z_BIN_SHIFT);
		if ( (z_position <= (-1*MAX_Z_POSITION)) || (z_position >= MAX_Z_POSITION) ) {
			score += 0;
		} else {
			// only average if both atoms were within range
			if ( 0 == score ) {
				score += asymEZ_CG_(rsd.aa(), z_bin);
			} else {
				score += asymEZ_CG_(rsd.aa(), z_bin);
				score /= 2;
			}
			//std::cout << rsd.aa() << rsd.seqpos() << "\tz = " << z_position << "\tz_bin = " << z_bin << "\tFaMPAsymEzCB = " << score << std::endl;
		}
	}

	emap[ scoring::FaMPAsymEzCG ] += score;
}


/// @details Only about 2/3 of residues have a specifically labeled CG atom; therefore, this
///    deals with the special cases of the remainder of residues (M, C, A, G, I, S,T, V)
std::string const &
FaMPAsymEzCGEnergy::representative_atom_name( chemical::AA const aa ) const
{
	debug_assert( aa >= 1 && aa <= chemical::num_canonical_aas );

	static std::string const cbeta_string(  "CB"  );
	static std::string const calpha_string( "CA"  );
	static std::string const cgamma_string( "CG"  );
	static std::string const sulf_string( "SG"  );
	static std::string const oxy_string( "OG"  );
	static std::string const cgamma1_string( "CG1"  );
	static std::string const cgamma2_string( "CG2"  );

	switch ( aa ) {
	case ( chemical::aa_ala ) : return cbeta_string;
	case ( chemical::aa_cys ) : return sulf_string;
	case ( chemical::aa_asp ) : return cgamma_string;
	case ( chemical::aa_glu ) : return cgamma_string;
	case ( chemical::aa_phe ) : return cgamma_string;
	case ( chemical::aa_gly ) : return calpha_string;
	case ( chemical::aa_his ) : return cgamma_string;
	case ( chemical::aa_ile ) : return cgamma1_string;
	case ( chemical::aa_lys ) : return cgamma_string;
	case ( chemical::aa_leu ) : return cgamma_string;
	case ( chemical::aa_met ) : return cgamma_string;
	case ( chemical::aa_asn ) : return cgamma_string;
	case ( chemical::aa_pro ) : return cgamma_string;
	case ( chemical::aa_gln ) : return cgamma_string;
	case ( chemical::aa_arg ) : return cgamma_string;
	case ( chemical::aa_ser ) : return oxy_string;
	case ( chemical::aa_thr ) : return cgamma2_string;
	case ( chemical::aa_val ) : return cgamma1_string;
	case ( chemical::aa_trp ) : return cgamma_string;
	case ( chemical::aa_tyr ) : return cgamma_string;
	default :
		utility_exit_with_message( "ERROR: Failed to find amino acid " + chemical::name_from_aa( aa ) + " in FAMPAsymEzCG::representative_atom_name" );
		break;
	}
	// unreachable
	return calpha_string;
}

void
FaMPAsymEzCGEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}

core::Size
FaMPAsymEzCGEnergy::version() const
{
	return 1; // Initial versioning
}

} // energy_methods
} // core

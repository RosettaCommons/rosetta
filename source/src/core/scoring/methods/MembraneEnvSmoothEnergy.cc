// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MembraneEnvSmoothEnergy.cc
/// @brief  Statistically derived smooth membrane environment potential class implementation
/// @brief  and based on EnvSmoothEnergy.cc developed by Andrew Leaver-Fay, Mike Tyka and Phil Bradley
/// @author Vladimir Yarov-Yarovoy

// Unit headers
#include <core/scoring/methods/MembraneEnvSmoothEnergy.hh>
#include <core/scoring/methods/MembraneEnvSmoothEnergyCreator.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/MembranePotential.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

// C++
static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.MembraneEnvSmoothEnergy", basic::t_info );

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the MembraneEnvSmoothEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MembraneEnvSmoothEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MembraneEnvSmoothEnergy );
}

ScoreTypes
MembraneEnvSmoothEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( Menv_smooth );
	return sts;
}

Distance const start_sig = 9.8;
Distance const end_sig   = 10.2;

DistanceSquared const start_sig2 = start_sig*start_sig;
DistanceSquared const end_sig2   = end_sig*end_sig;

////////////////////////////////////////////////////////////////////////////

MembraneEnvSmoothEnergy::MembraneEnvSmoothEnergy() :
	parent( methods::EnergyMethodCreatorOP( new MembraneEnvSmoothEnergyCreator ) )
{
	Size const max_aa( 20 );
	Size const env_log_table_cen10_bins( 40 );
	Size const max_mem_layers( 3 );
	std::string tag,line;
	chemical::AA aa;
	mem_env_log10_.dimension( max_aa, max_mem_layers, env_log_table_cen10_bins );
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/MembranePotential/CEN10_Menv_smooth_log.txt" );
	for ( Size i=1; i<= max_aa; ++i ) {
		getline( stream, line );
		std::istringstream l(line);
		l >> tag >> aa;
		if ( l.fail() || tag != "MENV_SMOOTH_LOG_CEN10:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/CEN10_Menv_smooth_log.txt (Menv_smooth)");
		for ( Size j=1; j<=max_mem_layers; ++j ) {
			for ( Size k=1; k<=env_log_table_cen10_bins; ++k ) {
				l >> mem_env_log10_(aa,j,k);
			}
		}
	}
}

/// clone
EnergyMethodOP
MembraneEnvSmoothEnergy::clone() const
{
	return EnergyMethodOP( new MembraneEnvSmoothEnergy( *this ) );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

inline Real sqr ( Real x ) {
	return x*x;
}


/// @details stores dScore/dNumNeighbors so that when neighbor atoms on adjacent
/// residues move, their influence on the score of the surrounding residues is
/// rapidly computed.
void
MembraneEnvSmoothEnergy::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const &
) const
{
	pose.update_residue_neighbors();
	Size const nres( pose.total_residue() );

	residue_N_.clear();
	residue_E_.clear();
	residue_dEdN_.clear();


	// iterate over all the residues in the protein and count their neighbours
	// and save values of E, N, and dEdN
	for ( Size i = 1; i <= nres; ++i ) {

		// get the appropriate residue from the pose.
		conformation::Residue const & rsd( pose.residue(i) );
		// currently this is only for protein residues
		if ( ! rsd.is_protein() ) continue; //return;

		Size const atomindex_i = rsd.atom_index( representative_atom_name( rsd.aa() ));

		core::conformation::Atom const & atom_i = rsd.atom(atomindex_i);

		const Energies & energies( pose.energies() );
		const TwelveANeighborGraph & graph ( energies.twelveA_neighbor_graph() );

		Real countN = 0.0;
		// iterate across neighbors within 12 angstroms
		for ( graph::Graph::EdgeListConstIter
				ir  = graph.get_node(i)->const_edge_list_begin(),
				ire = graph.get_node(i)->const_edge_list_end();
				ir != ire; ++ir ) {
			Size const j( (*ir)->get_other_ind( i ) );
			conformation::Residue const & rsd_j( pose.residue(j) );
			Size atomindex_j( rsd_j.type().nbr_atom() );

			core::conformation::Atom const & atom_j = rsd_j.atom(atomindex_j);

			Real sqdist = atom_i.xyz().distance_squared(atom_j.xyz());
			countN += sigmoidish_neighbor( sqdist );
		}

		Real score = 0;
		Real dscoredN = 0;

		calc_energy( rsd, pose, countN, rsd.aa(), score, dscoredN );

		residue_N_.push_back( countN );
		residue_E_.push_back( score );
		residue_dEdN_.push_back( dscoredN );

		//std::cout << "ENV:  " << i << "  " << score << std::endl;
	}
}

void
MembraneEnvSmoothEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const
{
	pose.update_residue_neighbors();
}

/// @details counts the number of nbr atoms within a given radius of the for the input
/// residue.  Because the representative atom on the input residue may be in a different
/// location than the representative atom on the same residue when scoring_begin() is called,
/// these neighbor counts cannot be reused; therefore, scoring_begin does not keep
/// neighbor counts.
void
MembraneEnvSmoothEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	// currently this is only for protein residues
	if ( ! rsd.is_protein() ) return;

	TwelveANeighborGraph const & graph ( pose.energies().twelveA_neighbor_graph() );
	Size const atomindex_i = rsd.atom_index( representative_atom_name( rsd.aa() ));

	core::conformation::Atom const & atom_i = rsd.atom(atomindex_i);

	Real countN    =  0.0;
	// iterate across neighbors within 12 angstroms
	using namespace ObjexxFCL::format;
	//TR << I(4,rsd.seqpos()) << std::endl;
	for ( graph::Graph::EdgeListConstIter
			ir  = graph.get_node( rsd.seqpos() )->const_edge_list_begin(),
			ire = graph.get_node( rsd.seqpos() )->const_edge_list_end();
			ir != ire; ++ir ) {
		Size const j( (*ir)->get_other_ind( rsd.seqpos() ) );
		conformation::Residue const & rsd_j( pose.residue(j) );

		// if virtual residue, don't score
		if ( rsd_j.aa() == core::chemical::aa_vrt ) continue;

		Size atomindex_j( rsd_j.nbr_atom() );

		core::conformation::Atom const & atom_j = rsd_j.atom(atomindex_j);

		Real sqdist = atom_i.xyz().distance_squared( atom_j.xyz() );
		countN += sigmoidish_neighbor( sqdist );
		//TR << I(4,j) << F(8,3,sqdist) << F(8,3,countN) << std::endl;
	}

	Real score = 0;
	Real dscoredN = 0;

	calc_energy( rsd, pose, countN, rsd.aa(), score, dscoredN );

	emap[ Menv_smooth ] += score;
}


/// @details Special cases handled for when an atom is both the representative
/// atom for an amino acid, and its nbr_atom.
void
MembraneEnvSmoothEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & ,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	conformation::Residue const & rsd = pose.residue( atom_id.rsd() );

	if ( ! rsd.is_protein() ) return;

	Size const i = rsd.seqpos();
	Size const i_nbr_atom = rsd.type().nbr_atom();
	Size const i_rep_atom = rsd.atom_index( representative_atom_name( rsd.aa() ));

	// forces act only on the nbr atom (CB or CA) or the representative atom
	if ( i_nbr_atom != (Size) atom_id.atomno() && i_rep_atom != (Size) atom_id.atomno() ) return;

	core::conformation::Atom const & atom_i = rsd.atom( atom_id.atomno() );

	TwelveANeighborGraph const & graph ( pose.energies().twelveA_neighbor_graph() );

	// its possible both of these are true
	bool const input_atom_is_nbr( i_nbr_atom == Size (atom_id.atomno()) );
	bool const input_atom_is_rep( i_rep_atom == Size ( atom_id.atomno() ));

	Vector f1(0.0), f2(0.0);

	for ( graph::Graph::EdgeListConstIter
			ir  = graph.get_node(i)->const_edge_list_begin(),
			ire = graph.get_node(i)->const_edge_list_end();
			ir != ire; ++ir ) {
		Size const j( (*ir)->get_other_ind( i ) );
		conformation::Residue const & rsd_j( pose.residue(j) );

		// if virtual residue, don't score
		if ( rsd_j.aa() == core::chemical::aa_vrt ) continue;

		if ( input_atom_is_nbr && input_atom_is_rep && rsd_j.is_protein() ) {
			Size const resj_rep_atom = rsd_j.atom_index( representative_atom_name( rsd_j.aa() ));
			Size const resj_nbr_atom = rsd_j.nbr_atom();
			if ( resj_rep_atom == resj_nbr_atom ) {
				/// two birds, one stone
				increment_f1_f2_for_atom_pair(
					atom_i, rsd_j.atom( resj_rep_atom ),
					weights[ Menv_smooth ] * ( residue_dEdN_[ j ] + residue_dEdN_[ i ] ),
					F1, F2 );
			} else {
				increment_f1_f2_for_atom_pair(
					atom_i, rsd_j.atom( resj_rep_atom ),
					weights[ Menv_smooth ] * ( residue_dEdN_[ j ] ),
					F1, F2 );

				increment_f1_f2_for_atom_pair(
					atom_i, rsd_j.atom( resj_nbr_atom ),
					weights[ Menv_smooth ] * ( residue_dEdN_[ i ] ),
					F1, F2 );
			}
		} else if ( input_atom_is_nbr && rsd_j.is_protein() ) {
			Size const resj_rep_atom = rsd_j.atom_index( representative_atom_name( rsd_j.aa() ));

			increment_f1_f2_for_atom_pair(
				atom_i, rsd_j.atom( resj_rep_atom ),
				weights[ Menv_smooth ] * ( residue_dEdN_[ j ] ),
				F1, F2 );

		} else {
			Size const resj_nbr_atom = rsd_j.nbr_atom();
			increment_f1_f2_for_atom_pair(
				atom_i, rsd_j.atom( resj_nbr_atom ),
				weights[ Menv_smooth ] * ( residue_dEdN_[ i ] ),
				F1, F2 );

		}
	}
}

/// @details returns const & to static data members to avoid expense
/// of string allocation and destruction.  Do not call this function
/// on non-canonical aas
std::string const &
MembraneEnvSmoothEnergy::representative_atom_name( chemical::AA const aa ) const
{
	debug_assert( aa >= 1 && aa <= chemical::num_canonical_aas );

	static std::string const cbeta_string(  "CB"  );
	static std::string const sgamma_string( "SG"  );
	static std::string const cgamma_string( "CG"  );
	static std::string const cdelta_string( "CD"  );
	static std::string const czeta_string(  "CZ"  );
	static std::string const calpha_string( "CA"  );
	static std::string const ceps_1_string( "CE1" );
	static std::string const cdel_1_string( "CD1" );
	static std::string const ceps_2_string( "CE2" );
	static std::string const sdelta_string( "SD"  );

	switch ( aa ) {
	case ( chemical::aa_ala ) : return cbeta_string;
	case ( chemical::aa_cys ) : return sgamma_string;
	case ( chemical::aa_asp ) : return cgamma_string;
	case ( chemical::aa_glu ) : return cdelta_string;
	case ( chemical::aa_phe ) : return czeta_string;
	case ( chemical::aa_gly ) : return calpha_string;
	case ( chemical::aa_his ) : return ceps_1_string;
	case ( chemical::aa_ile ) : return cdel_1_string;
	case ( chemical::aa_lys ) : return cdelta_string;
	case ( chemical::aa_leu ) : return cgamma_string;
	case ( chemical::aa_met ) : return sdelta_string;
	case ( chemical::aa_asn ) : return cgamma_string;
	case ( chemical::aa_pro ) : return cgamma_string;
	case ( chemical::aa_gln ) : return cdelta_string;
	case ( chemical::aa_arg ) : return czeta_string;
	case ( chemical::aa_ser ) : return cbeta_string;
	case ( chemical::aa_thr ) : return cbeta_string;
	case ( chemical::aa_val ) : return cbeta_string;
	case ( chemical::aa_trp ) : return ceps_2_string;
	case ( chemical::aa_tyr ) : return czeta_string;
	default :
		utility_exit_with_message( "ERROR: Failed to find amino acid " + chemical::name_from_aa( aa ) + " in EnvSmooth::representative_atom_name" );
		break;
	}

	// unreachable
	return calpha_string;
}

/// @brief MembraneEnvSmoothEnergy distance cutoff
Distance
MembraneEnvSmoothEnergy::atomic_interaction_cutoff() const
{
	return 0.0;
}

/// @brief MembraneEnvSmoothEnergy
void
MembraneEnvSmoothEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ twelve_A_neighbor_graph ] = true;
}

void
MembraneEnvSmoothEnergy::calc_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	Real const neighbor_count,
	chemical::AA const aa,
	Real & score,
	Real & dscore_dneighbor_count
) const
{
	Size low_bin = static_cast< Size > ( floor(neighbor_count));
	Size high_bin = static_cast< Size > ( ceil(neighbor_count));
	Real inter = neighbor_count - low_bin;
	if ( neighbor_count <= 1 ) {
		low_bin = high_bin = 1; inter = 0;
	}

	Real const MembraneDepth (MembraneEmbed_from_pose( pose ).depth( rsd.seqpos() ) );

	Real thickness = 2.0;
	int  slope = 14;
	int  layer1,layer2,layer;
	Real f,z,zn,layer_edge;

	Real const env10_weight=1.0;

	if ( high_bin < 41 ) {
		if ( ( MembraneDepth < 11.0 ) || ( MembraneDepth > 49.0 ) ) {
			//pure water layer
			layer = 3;

			score = env10_weight * mem_env_log10_( aa, layer, low_bin ) * (1.0-inter) +
				env10_weight * mem_env_log10_( aa, layer, high_bin ) * (inter);

			dscore_dneighbor_count = env10_weight * mem_env_log10_( aa, layer, high_bin ) -
				env10_weight * mem_env_log10_( aa, layer, low_bin );

		} else if ( ( MembraneDepth >= 11.0 && MembraneDepth <= 13.0 ) || ( MembraneDepth >= 47.0 && MembraneDepth <= 49.0 ) ) {
			//interpolate between water and interface phases
			layer1 = 2; //interface layer
			layer2 = 3; //water layer

			if ( MembraneDepth <= 13.0 ) {
				layer_edge = 13.0;
			} else {
				layer_edge = 47.0;
			}
			z = 2*std::abs( (MembraneDepth - layer_edge ) ) / thickness;
			zn = std::pow( z, slope );
			f = zn/(1 + zn);

			score = f * ( env10_weight * mem_env_log10_( aa, layer2, low_bin ) * (1.0-inter) +
				env10_weight * mem_env_log10_( aa, layer2, high_bin ) * (inter) ) +
				( 1 - f ) * ( env10_weight * mem_env_log10_( aa, layer1, low_bin ) * (1.0-inter) +
				env10_weight * mem_env_log10_( aa, layer1, high_bin ) * (inter) );

			dscore_dneighbor_count = f * ( env10_weight * mem_env_log10_( aa, layer2, high_bin ) -
				env10_weight * mem_env_log10_( aa, layer2, low_bin ) ) +
				( 1 - f ) * ( env10_weight * mem_env_log10_( aa, layer1, high_bin ) -
				env10_weight * mem_env_log10_( aa, layer1, low_bin ) );

			if ( MembraneDepth <= 12.0 || MembraneDepth >= 48.0 ) {
				layer = 2;
			} else {
				layer = 3;
			}

		} else if ( ( MembraneDepth > 13.0 && MembraneDepth < 17.0 ) || ( MembraneDepth > 43.0 && MembraneDepth < 47.0 ) ) {
			//pure interface phase
			layer = 2;

			score = env10_weight * mem_env_log10_( aa, layer, low_bin ) * (1.0-inter) +
				env10_weight * mem_env_log10_( aa, layer, high_bin ) * (inter);

			dscore_dneighbor_count = env10_weight * mem_env_log10_( aa, layer, high_bin ) -
				env10_weight * mem_env_log10_( aa, layer, low_bin );

		} else if ( ( MembraneDepth >= 17.0 && MembraneDepth <= 19.0 ) || ( MembraneDepth >= 41.0 && MembraneDepth <= 43.0 ) ) {
			//interpolate between interface and hydrophobic phases
			layer1 = 1; //hydrophobic layer
			layer2 = 2; //interface layer

			if ( MembraneDepth <= 19.0 ) {
				layer_edge = 19.0;
			} else {
				layer_edge = 41.0;
			}
			z = 2*std::abs( (MembraneDepth - layer_edge ) ) / thickness;
			zn = std::pow( z, slope );
			f = zn/(1 + zn);

			score = f * ( env10_weight * mem_env_log10_( aa, layer2, low_bin ) * (1.0-inter) +
				env10_weight * mem_env_log10_( aa, layer2, high_bin ) * (inter) ) +
				( 1 - f ) * ( env10_weight * mem_env_log10_( aa, layer1, low_bin ) * (1.0-inter) +
				env10_weight * mem_env_log10_( aa, layer1, high_bin ) * (inter) );

			dscore_dneighbor_count = f * ( env10_weight * mem_env_log10_( aa, layer2, high_bin ) -
				env10_weight * mem_env_log10_( aa, layer2, low_bin ) ) +
				( 1 - f ) * ( env10_weight * mem_env_log10_( aa, layer1, high_bin ) -
				env10_weight * mem_env_log10_( aa, layer1, low_bin ) );

		} else {
			//pure hydrophobic phase
			layer = 1;

			score = env10_weight * mem_env_log10_( aa, layer, low_bin ) * (1.0-inter) +
				env10_weight * mem_env_log10_( aa, layer, high_bin ) * (inter);

			dscore_dneighbor_count = env10_weight * mem_env_log10_( aa, layer, high_bin ) -
				env10_weight * mem_env_log10_( aa, layer, low_bin );
		}

	} else { //high_bin >= 40 neighbors - does it ever happen?

		dscore_dneighbor_count = 0;

		if ( ( MembraneDepth < 11.0 ) || ( MembraneDepth > 49.0 ) ) {
			//pure water layer
			layer = 3;

			score = env10_weight * mem_env_log10_( aa, layer, 40 );

		} else if ( ( MembraneDepth >= 11.0 && MembraneDepth <= 13.0 ) || ( MembraneDepth >= 47.0 && MembraneDepth <= 49.0 ) ) {
			//interpolate between water and interface phases
			layer1 = 2; //interface layer
			layer2 = 3; //water layer

			if ( MembraneDepth <= 13.0 ) {
				layer_edge = 13.0;
			} else {
				layer_edge = 47.0;
			}
			z = 2*std::abs( (MembraneDepth - layer_edge ) ) / thickness;
			zn = std::pow( z, slope );
			f = zn/(1 + zn);

			score = f * ( env10_weight * mem_env_log10_( aa, layer2, 40 ) ) +
				( 1 - f ) * ( env10_weight * mem_env_log10_( aa, layer1, 40 ) );

			if ( MembraneDepth <= 12.0 || MembraneDepth >= 48.0 ) {
				layer = 2;
			} else {
				layer = 3;
			}

		} else if ( ( MembraneDepth > 13.0 && MembraneDepth < 17.0 ) || ( MembraneDepth > 43.0 && MembraneDepth < 47.0 ) ) {
			//pure interface phase
			layer = 2;

			score = env10_weight * mem_env_log10_( aa, layer, 40 );

		} else if ( ( MembraneDepth >= 17.0 && MembraneDepth <= 19.0 ) || ( MembraneDepth >= 41.0 && MembraneDepth <= 43.0 ) ) {
			//interpolate between interface and hydrophobic phases
			layer1 = 1; //hydrophobic layer
			layer2 = 2; //interface layer

			if ( MembraneDepth <= 19.0 ) {
				layer_edge = 19.0;
			} else {
				layer_edge = 41.0;
			}
			z = 2*std::abs( (MembraneDepth - layer_edge ) ) / thickness;
			zn = std::pow( z, slope );
			f = zn/(1 + zn);

			score = f * ( env10_weight * mem_env_log10_( aa, layer2, 40 ) ) +
				( 1 - f ) * ( env10_weight * mem_env_log10_( aa, layer1, 40 ) );

		} else {
			//pure hydrophobic phase
			layer = 1;

			score = env10_weight * mem_env_log10_( aa, layer, 40 );
		}
	}

	score *= 2.019; // this factor is from rosetta++
	dscore_dneighbor_count *= 2.019;
}


Real
MembraneEnvSmoothEnergy::sigmoidish_neighbor( DistanceSquared const sqdist ) const
{
	if ( sqdist > end_sig2 ) {
		return 0.0;
	} else if ( sqdist < start_sig2 ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0  - sqr( (dist - start_sig) / (end_sig - start_sig) ) );
	}
}


void
MembraneEnvSmoothEnergy::increment_f1_f2_for_atom_pair(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real weighted_dScore_dN,
	Vector & F1,
	Vector & F2
) const
{
	DistanceSquared dist2 = atom1.xyz().distance_squared(atom2.xyz());
	Distance dist( 0.0 ); // only used if start_sig2 <= dist2 <= end_sig2

	Real dNdd = 0;
	if ( dist2 > end_sig2 ) {
		dNdd = 0;
	} else if ( dist2 < start_sig2 ) {
		dNdd = 0.0;
	} else {
		dist = sqrt( dist2 );
		Real x = (dist - start_sig)/ (end_sig - start_sig);
		dNdd = 4*x*(-1 + x*x) / (end_sig - start_sig);
	}

	Real dscoredd = ( weighted_dScore_dN ) * dNdd;
	if ( dscoredd != 0 ) {

		Vector const f1( cross( atom1.xyz(), atom2.xyz() ));
		Vector const f2( atom1.xyz() - atom2.xyz() );

		dscoredd /= dist;
		F1 += dscoredd * f1;
		F2 += dscoredd * f2;
	}
}
core::Size
MembraneEnvSmoothEnergy::version() const
{
	return 1; // Initial versioning
}
}
}
}

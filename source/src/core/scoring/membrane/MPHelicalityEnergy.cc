// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPHelicalityEnergy.cc
///
/// @brief  Fullatom and centroid level smooth membrane non-helicality penalty
/// @details
/// @FlesihmanLab
/// @Last Modified: 20/2/17
///
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)
/// @author Assaf Elazar
/// @author Sarel Fleishman

// Unit headers
#include <core/scoring/membrane/MPHelicalityEnergy.hh>
#include <core/scoring/membrane/MPHelicalityEnergyCreator.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>

#include <core/conformation/Atom.hh>
#include <core/id/AtomID.hh>

#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>

#include <utility/graph/Graph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/membrane/MembraneData.hh>

#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <numeric/interpolation/spline/CubicSpline.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathVector.srlz.hh>
#include <boost/lexical_cast.hpp>
#include <map>
#include <math.h>
#include <set>


// Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
// Auto using namespaces end
using namespace core::scoring;
using namespace core::scoring::methods;

static basic::Tracer TR( "core.scoring.membrane.MPHelicalityEnergy" );

namespace core {
namespace scoring {
namespace membrane {

// sigmoid constants
Real const slope_6 = 0.15; //0.5; //0.05;//0.197; //0.05231;
Real const offset_6 = 20.0; //16.0; //30.0;

Real const slope_12 = 0.5; //0.02;//0.197; //0.05231;
Real const offset_12 = 475.0; //350.0;

Real const sixA_slope = 0.15;
Real const sixA_offset = 20;

Real const twelveA_slope = 5;
Real const twelveA_offset = 220.0;
// Creator Methods //////////////////////////////////////////

/// @brief Return a fresh instance of the energy method
methods::EnergyMethodOP
MPHelicalityEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MPHelicalityEnergy );
}

/// @brief Return relevant score types
ScoreTypes
MPHelicalityEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPHelicality );
	return sts;
}

// Constructors /////////////////////////////////////////////

MPHelicalityEnergy::MPHelicalityEnergy() :
	parent( EnergyMethodCreatorOP( new MPHelicalityEnergyCreator ) )
{
}

/// @brief Create a clone of this energy method
EnergyMethodOP
MPHelicalityEnergy::clone() const
{
	return EnergyMethodOP( new MPHelicalityEnergy( *this ) );
}

/// @details stores dScore/dNumNeighbors so that when neighbor atoms on adjacent
/// residues move, their influence on the score of the surrounding residues is
/// rapidly computed.
void
MPHelicalityEnergy::setup_for_derivatives(
	pose::Pose &,
	ScoreFunction const &
) const
{
}

void
MPHelicalityEnergy::setup_for_scoring(
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
MPHelicalityEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	// currently this is only for protein residues
	if ( ! rsd.is_protein() ) return;

	if ( oneletter_code_from_aa( rsd.aa() ) == 'P' ) return;
	if ( rsd.is_lower_terminus() || rsd.is_upper_terminus() ) return;
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Error: Cannot use mpframework energy term using a non membrane pose!");
	}
	Real score = 0;
	// Skip Membrane Residues
	//core::Size nres = pose.size()-1;
	core::Size nres = core::pose::nres_protein( pose );
	if ( rsd.seqpos() > nres ) return;

	// Skip Cases
	if ( rsd.seqpos() == 0 ) return;
	if ( rsd.aa() == core::chemical::aa_vrt ) return;
	Real const burial_sig_6x12 = calc_residue_burial( pose, rsd );


	calc_energy( rsd, pose, score );

	emap[ MPHelicality ] += score * burial_sig_6x12;
}


/// @brief approximate residue burial for either centroid mode: using tenA or for
/// full-atom using counts of atoms 6 and 12A away. either way a sigmoid is used
/// to approximate burial from the counts
core::Real
MPHelicalityEnergy::calc_residue_burial( pose::Pose const & pose,
	conformation::Residue const & rsd) const
{

	if ( !pose.is_centroid() ) {
		utility::vector1 < core::Size > const atom_counts = neighboring_atoms(pose, rsd,  6.0, 12.0);

		core::Size atom_count_1 = atom_counts[1];
		core::Size atom_count_2 = atom_counts[2];
		Real const burial_sig_6A = burial_sigmoid( atom_count_1, slope_6, offset_6 );
		Real const burial_sig_12A = burial_sigmoid( atom_count_2, slope_12, offset_12 );

		Real burial_sig_6x12 = burial_sig_6A * burial_sig_12A;
		if ( burial_sig_6x12 < 0 ) burial_sig_6x12 = 0.0;
		if ( burial_sig_6x12 > 1 ) burial_sig_6x12 = 1.0;

		return( burial_sig_6x12 );
	} else {
		//int tenA_neighbors = pose.energies().tenA_neighbor_graph().get_node( rsd.seqpos() )->num_neighbors_counting_self();
		utility::vector1< core::Size > centroid_nbrs = centroid_neighbors(pose, rsd);
		//core::Real tenA_sig = burial_sigmoid(tenA_neighbors, slope_tenA, offset_tenA);
		core::Size sixA_n = centroid_nbrs[1];
		core::Size twelveA_n = centroid_nbrs[2];

		core::Real sixA_sig = burial_sigmoid(sixA_n, sixA_slope, sixA_offset); //slope_tenA, offset_tenA);
		core::Real twelveA_sig = burial_sigmoid(twelveA_n, twelveA_slope, twelveA_offset); //slope_twelveA, offset_twelveA);
		//int twelveA_neighbors = pose.energies().tenA_neighbor_graph().get_node( rsd.seqpos()  )->num_neighbors_counting_self();
		//core::Real twelveA_sig = burial_sigmoid(twelveA_neighbors, slope_twelveA, offset_twelveA);
		core::Real cen_bur = sixA_sig * twelveA_sig;
		if ( cen_bur > 1 ) cen_bur = 1.0;
		if ( cen_bur < 0 ) cen_bur = 0.0;

		return( cen_bur );
	}
}

utility::vector1 < core::Size >
MPHelicalityEnergy::centroid_neighbors(
	pose::Pose const & pose,
	conformation::Residue const & rsd) const
{
	std::set < core::Size > n_12A;
	std::set < core::Size > n_6A;
	n_12A.clear();
	n_6A.clear();

	core::Size atom_id = 0;

	//int twelveA_neighbors = pose.energies().tenA_neighbor_graph().get_node( rsd.seqpos() )->num_neighbors_counting_self();
	//const TwelveANeighborGraph & graph ( pose.energies.twelveA_neighbor_graph() );
	core::scoring::TwelveANeighborGraph const & graph ( pose.energies().twelveA_neighbor_graph()  );
	utility::graph::Node const * node = graph.get_node( rsd.seqpos()  );
	//core::graph::Node const * node = pose.energies().twelveA_neighbor_graph().get_node( rsd.seqpos() );

	Real distance = 0.0;
	for ( utility::graph::EdgeListConstIterator edge = node->const_edge_list_begin(); edge != node->const_edge_list_end(); ++edge ) {
		core::Real edge_res_i = (*edge)->get_other_ind( rsd.seqpos() );
		if ( edge_res_i >= rsd.seqpos()-4 && edge_res_i <= rsd.seqpos()+4 ) {
			continue;
		}

		for ( core::Size atom_i=1; atom_i <= pose.residue( edge_res_i ).natoms(); atom_i++ ) {
			distance = pose.residue( edge_res_i ).xyz( atom_i ).distance( rsd.xyz( rsd.nbr_atom() ) );
			atom_id = (edge_res_i * 1000) + atom_i;
			if ( distance <= 12.0 ) {
				n_12A.insert( atom_id );
				if ( distance <= 6.0 ) {
					n_6A.insert( atom_id );
				}
			}
		}
	}

	utility::vector1 < core::Size > result;
	result.empty();
	result.push_back( static_cast < core::Size > ( n_6A.size() ) );
	result.push_back( static_cast < core::Size > ( n_12A.size() ) );
	return( result );
}

/// @details Special cases handled for when an atom is both the representative
/// atom for an amino acid, and its nbr_atom.
void
MPHelicalityEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & ,
	ScoreFunction const &,
	EnergyMap const &,
	Vector &,
	Vector &
) const
{
	conformation::Residue const & rsd = pose.residue( atom_id.rsd() );

	if ( ! rsd.is_protein() ) return;
	// leaving the derivatives as 0.0, because residue solvation is relatively
	// flat, so the derivative is almost 0...
	Vector f1(0.0), f2(0.0);
}

/// @brief MPHelicalityEnergy distance cutoff
Distance
MPHelicalityEnergy::atomic_interaction_cutoff() const
{
	return 0.0;
}

/// @brief MPHelicalityEnergy
void
MPHelicalityEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	// context_graphs_required[ twelve_A_neighbor_graph ] = true;
	//context_graphs_required[ ten_A_neighbor_graph ] = true;
	context_graphs_required[ twelve_A_neighbor_graph  ] = true;
}

void
MPHelicalityEnergy::calc_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	Real & score
) const
{
	Real const z_position( pose.conformation().membrane_info()->residue_z_position( pose.conformation(), rsd.seqpos() ) );
	//bool full_atom = !pose.is_centroid();

	Real res_energy = 0.0;
	core::Real phi = pose.phi( rsd.seqpos() );
	core::Real psi = pose.psi( rsd.seqpos() );

	// values for ideal membrane helix taken from Page, R. C., Kim, S., & Cross, T. A. (2008).
	// Transmembrane Helix Uniformity Examined by Spectral Mapping of Torsion Angles.
	// Structure, 16(5), 787–797. http://doi.org/10.1016/j.str.2008.02.018
	Real phi_score = std::pow( phi + 60, 2 );
	Real psi_score = std::pow( psi + 45, 2 );

	Real za = std::pow( z_position / 10, 4 );
	Real f_z = za / ( 1 + za );
	//TR << rsd.seqpos() << " " << rsd.aa() << " " << phi << " " << psi << " " << phi_score << " " << psi_score << "  " << z_position << " " << (1 / std::pow( 25, 4  ))* std::pow( phi_score + psi_score , 2  ) * f_z << std::endl;
	res_energy += (1 / std::pow( 25, 4 ))* std::pow( phi_score + psi_score , 2 ) * f_z; // a parabola around (0, 0) with (35, 1)
	// where 35 is the maximal distance from the center that can still be considered an a-helix. sqrt()^4 = ^2
	score += res_energy;
}


bool
MPHelicalityEnergy::minimize_in_whole_structure_context( pose::Pose const &   ) const { return false; }


// returns the number of atoms closer to rsd than both cutoffs
utility::vector1 < core::Size >
MPHelicalityEnergy::neighboring_atoms(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Real const & cutoff_1,
	Real const & cutoff_2
) const
{
	// Tg << "NERIGHBORING ATOMS" << std::endl;
	Real residue_distance_cutoff = 16.0;
	utility::vector1 < core::Size > result;
	result.empty();

	std::set < core::Size > atoms_seen_1;
	std::set < core::Size > atoms_seen_2;
	atoms_seen_1.clear();
	atoms_seen_2.clear();

	core::Size atom_id = 0;
	Real distance = 0.0;
	// has to be Real for negative numbers...
	core::Real rsd_seqpos = rsd.seqpos();

	// iter all residues, if they are closer to rsd than residue_distance_cutoff
	// than check all their atoms
	for ( core::Size res_i=1; res_i <= pose.size(); res_i++ ) {

		if ( res_i >= rsd_seqpos-4 && res_i <= rsd_seqpos+4 ) {
			continue;
		} else { }

		core::conformation::Residue const & residue_i = pose.residue(res_i);
		if ( residue_i.xyz(residue_i.nbr_atom()).distance(rsd.xyz(rsd.nbr_atom())) <= residue_distance_cutoff ) {

			for ( core::Size atom_i=1; atom_i <= residue_i.natoms(); atom_i++ ) {

				// count any atom closer than cutoff and push to the seen vector
				distance = residue_i.xyz(atom_i).distance(rsd.xyz(rsd.nbr_atom()));
				if ( distance <= cutoff_2 ) {
					atom_id = (res_i * 1000) + atom_i;
					atoms_seen_2.insert( atom_id );
					if ( distance <= cutoff_1 ) {
						atoms_seen_1.insert( atom_id );
					}
				}
			}
		}
	}
	result.push_back( static_cast < core::Size > (atoms_seen_1.size()) );
	result.push_back( static_cast < core::Size > (atoms_seen_2.size()) );
	return result;
}

Real
MPHelicalityEnergy::burial_sigmoid(
	core::Size const n_atoms,
	Real const slope,
	Real const offset
) const
{
	return 1.0 / ( 1.0 + exp( ( n_atoms - offset ) * slope ));
}

core::Size
MPHelicalityEnergy::version() const
{
	return 1; // Initial versioning
}

} // membrane
} // scoring
} // core

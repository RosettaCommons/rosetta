// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPResidueLipophilicityEnergy.cc
///
/// @brief  Fullatom Smoothed Membrane Environment Energy
/// @details residue speicific enrgy by membrane depth, according to the Elazar
/// hydrophobicity scale
///    @FlesihmanLab.
///    Last Modified: 4/4/16
///
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)
/// @author Assaf Elazar
/// @author Sarel Fleishman

// Unit headers
#include <core/scoring/membrane/MPResidueLipophilicityEnergy.hh>
#include <core/scoring/membrane/MPResidueLipophilicityEnergyCreator.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>

#include <core/conformation/Atom.hh>
#include <core/id/AtomID.hh>

#include <core/scoring/Energies.hh>

#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>

#include <utility/graph/Graph.hh>
#include <core/scoring/Energies.hh>
//#include <core/scoring/TenANeighborGraph.hh>
//#include <core/graph/Graph.hh>
#include <core/scoring/membrane/MembraneData.hh>

#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <numeric/MathMatrix.hh>
#include <numeric/MathVector.srlz.hh>
#include <boost/lexical_cast.hpp>
#include <map>
#include <math.h>
#include <set>


// this is stupid...we really shouldn't have to do things like this anymore...
// Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
// Auto using namespaces end
using namespace core::scoring;
using namespace core::scoring::methods;

static basic::Tracer TR( "core.scoring.membrane.MPResidueLipophilicityEnergy" );

namespace core {
namespace scoring {
namespace membrane {

// Distance Constants

// sigmoid constants
Real const slope_6 = 0.15;
Real const offset_6 = 20.0;

Real const slope_12 = 0.5;
Real const offset_12 = 475.0;

Real const sixA_slope = 0.15;
Real const sixA_offset = 20;

Real const twelveA_slope = 5;
Real const twelveA_offset = 220.0;
// Creator Methods //////////////////////////////////////////

/// @brief Return a fresh instance of the energy method
methods::EnergyMethodOP
MPResidueLipophilicityEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MPResidueLipophilicityEnergy );
}

/// @brief Return relevant score types
ScoreTypes
MPResidueLipophilicityEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPResidueLipophilicity );
	return sts;
}

// Constructors /////////////////////////////////////////////

MPResidueLipophilicityEnergy::MPResidueLipophilicityEnergy() :
	parent( EnergyMethodCreatorOP( new MPResidueLipophilicityEnergyCreator ) )
{
	//parse_elazar_polynom_table();
	parse_elazar_spline_table();
}

/// @brief Create a clone of this energy method
EnergyMethodOP
MPResidueLipophilicityEnergy::clone() const
{
	return EnergyMethodOP( new MPResidueLipophilicityEnergy( *this ) );
}

/// @details stores dScore/dNumNeighbors so that when neighbor atoms on adjacent
/// residues move, their influence on the score of the surrounding residues is
/// rapidly computed.
void
MPResidueLipophilicityEnergy::setup_for_derivatives(
	pose::Pose &,
	ScoreFunction const &
) const
{
}

void
MPResidueLipophilicityEnergy::setup_for_scoring(
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
MPResidueLipophilicityEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	// currently this is only for protein residues
	if ( ! rsd.is_protein() ) return;

	Real score = 0;

	Real const burial_sig_6x12 = calc_residue_burial( pose, rsd );


	calc_energy( rsd, pose, score );

	emap[ MPResidueLipophilicity ] += score * burial_sig_6x12;
}

/// @brief called by the ResidueLipophilicityFilter to report the measurements made by ResidueLipophilicity
/// mostly useful for debugging
core::Real
MPResidueLipophilicityEnergy::report_ressolv(
	std::ostream & out,
	pose::Pose const & pose,
	bool print_splines
) const
{
	core::Real total_score = 0;
	out << "residue\tZ\tna1\tna2\tsig1\tsig2\tsigx\tcen6N\tcen6sig\tcen12N\tcen12sig\tcen10x12\tcen_score\tscorex\tscore" << std::endl;
	for ( core::Size res_i=1; res_i <= pose.size(); res_i++ ) {
		core::conformation::Residue const & residue_i = pose.residue(res_i);
		if ( ! residue_i.is_protein()  ) continue;

		Real const z_position( pose.conformation().membrane_info()->residue_z_position( pose.conformation(), res_i ));

		utility::vector1 < core::Size > const atom_counts = neighboring_atoms(pose, residue_i,  6.0, 12.0);

		core::Size atom_count_1 = atom_counts[1];
		core::Size atom_count_2 = atom_counts[2];
		Real const burial_sig_6A = burial_sigmoid( atom_count_1, slope_6, offset_6 );
		Real const burial_sig_12A = burial_sigmoid( atom_count_2, slope_12, offset_12 );

		Real burial_sig_6x12 = ( burial_sig_6A * burial_sig_12A ); // - 0.38) / ( 0.82 - 0.38 );
		if ( burial_sig_6x12 < 0 ) burial_sig_6x12 = 0.0;
		if ( burial_sig_6x12 > 1 ) burial_sig_6x12 = 1.0;

		utility::vector1< core::Size > centroid_nbrs = centroid_neighbors(pose, residue_i);
		core::Size sixA_n = centroid_nbrs[1];
		core::Size twelveA_n = centroid_nbrs[2];
		core::Real sixA_sig = burial_sigmoid(sixA_n, sixA_slope, sixA_offset) ; //slope_tenA, offset_tenA);
		core::Real twelveA_sig = burial_sigmoid(twelveA_n, twelveA_slope, twelveA_offset); //slope_twelveA, offset_twelveA);
		Real cen_burial = sixA_sig * twelveA_sig;
		if ( cen_burial < 0  ) cen_burial = 0.0;
		if ( cen_burial > 1  ) cen_burial = 1.0;

		Real score = 0;
		calc_energy( residue_i, pose, score );
		Real score_burial = score * burial_sig_6x12;

		out << residue_i.seqpos() << residue_i.aa() << "\t" << z_position;
		out << "\t" << atom_count_1 << "\t" << atom_count_2;
		out << "\t" << burial_sig_6A << "\t" << burial_sig_12A;
		out << "\t" << burial_sig_6x12 << "\t" << sixA_n;
		out << "\t" << sixA_sig << "\t" << twelveA_n;
		out << "\t" << twelveA_sig << "\t" << cen_burial;
		out << "\t"<< score * cen_burial;
		out << "\t" << score_burial << "\t" << score << std::endl;

		total_score += score_burial;
	}
	out << "total_mp_res_lipo\t" << total_score << std::endl;

	if ( print_splines ) {
		bool full_atom = !pose.is_centroid();
		out << "ResidueLipophilicity splines see:" << std::endl;
		std::string AAs = "ACDEFGHIKLMNPQRSTVWY";
		for ( std::string::iterator it=AAs.begin(); it!=AAs.end(); ++it ) {
			out << *it << "\t";
			for ( core::Real z=-134.5/2; z<=134.5/2; ++z ) {
				out << spline_by_z( *it, z, full_atom ) << "\t";
			}
			out << std::endl;
		}
	}
	return( total_score );
}

/// @brief approximate residue burial for either centroid mode: using tenA or for
/// full-atom using counts of atoms 6 and 12A away. either way a sigmoid is used
/// to approximate burial from the counts
core::Real
MPResidueLipophilicityEnergy::calc_residue_burial( pose::Pose const & pose,
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
		utility::vector1< core::Size > centroid_nbrs = centroid_neighbors(pose, rsd);
		core::Size sixA_n = centroid_nbrs[1];
		core::Size twelveA_n = centroid_nbrs[2];

		core::Real sixA_sig = burial_sigmoid(sixA_n, sixA_slope, sixA_offset); //slope_tenA, offset_tenA);
		core::Real twelveA_sig = burial_sigmoid(twelveA_n, twelveA_slope, twelveA_offset); //slope_twelveA, offset_twelveA);
		core::Real cen_bur = sixA_sig * twelveA_sig;
		if ( cen_bur > 1 ) cen_bur = 1.0;
		if ( cen_bur < 0 ) cen_bur = 0.0;

		return( cen_bur );
	}
}

utility::vector1 < core::Size >
MPResidueLipophilicityEnergy::centroid_neighbors(
	pose::Pose const & pose,
	conformation::Residue const & rsd) const
{
	std::set < core::Size > n_12A;
	std::set < core::Size > n_6A;
	n_12A.clear();
	n_6A.clear();

	core::Size atom_id = 0;

	const Energies & energies( pose.energies() );
	const core::scoring::TwelveANeighborGraph & graph ( energies.twelveA_neighbor_graph()  );

	Real distance = 0.0;
	for ( utility::graph::Graph::EdgeListConstIter ir = graph.get_node( rsd.seqpos() )->const_edge_list_begin(),
			ire = graph.get_node( rsd.seqpos() )->const_edge_list_end();
			ir != ire; ++ir ) {
		core::Real edge_res_i = (*ir)->get_other_ind( rsd.seqpos() );
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
MPResidueLipophilicityEnergy::eval_atom_derivative(
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

/// @brief MPResidueLipophilicityEnergy distance cutoff
Distance
MPResidueLipophilicityEnergy::atomic_interaction_cutoff() const
{
	return 0.0;
}

/// @brief MPResidueLipophilicityEnergy
void
MPResidueLipophilicityEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	// context_graphs_required[ twelve_A_neighbor_graph ] = true;
	//context_graphs_required[ ten_A_neighbor_graph ] = true;
	context_graphs_required[ twelve_A_neighbor_graph ] = true;
}

void
MPResidueLipophilicityEnergy::calc_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	Real & score
) const
{
	Real const z_position( pose.conformation().membrane_info()->residue_z_position( pose.conformation(), rsd.seqpos() ) );
	bool full_atom = !pose.is_centroid();

	char const res_a = static_cast< char > (oneletter_code_from_aa(rsd.aa()));
	Real res_energy = 0.0;
	res_energy = spline_by_z( res_a, z_position, full_atom );

	// before 9Oct216 penalised termini by D/K respective to terminus type. on 9Oct2016 moved to a solid +50...:
	if ( z_position >= -10 && z_position <= 10 ) {
		if ( rsd.is_upper_terminus() ) {
			res_energy += 50;
		} else if ( rsd.is_lower_terminus() ) {
			res_energy += 50;
		}
	}

	score += res_energy;
}

void
MPResidueLipophilicityEnergy::parse_elazar_polynom_table()
{
	// Read in updated table
	utility::io::izstream stream_1;
	basic::database::open( stream_1, "scoring/score_functions/MembranePotential/ELazaridis_polynom_table.txt" );
	std::string line_1;

	while ( stream_1  ) {
		getline( stream_1, line_1 );
		std::istringstream l_1(line_1);

		if ( line_1.size() == 0 ) break;

		res_polynom_map_.insert( std::pair< char, utility::vector1< Real > > (line_1[0], split_line_to_polyval( line_1 )));
	}

	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/MembranePotential/ELazaridis_cen_polynom_table.txt" );
	std::string line;

	while ( stream  ) {
		getline( stream, line );
		std::istringstream l(line);

		if ( line.size() == 0 ) break;

		cen_res_polynom_map_.insert( std::pair< char, utility::vector1< Real > > ( line[0], split_line_to_polyval( line ) ) );
	}

}

bool
MPResidueLipophilicityEnergy::minimize_in_whole_structure_context( pose::Pose const &   ) const { return false; }

void
MPResidueLipophilicityEnergy::parse_elazar_spline_table()
{
	TR << "parse_elazar_spline_table()" << std::endl;

	utility::io::izstream stream_1;
	basic::database::open( stream_1, "scoring/score_functions/MembranePotential/spline_test_fa.txt" );
	std::string line_1;

	TR << "parsing for full atom" << std::endl;
	while ( stream_1  ) {
		getline( stream_1, line_1 );
		if ( line_1[0] == '#' ) continue;
		std::istringstream l_1(line_1);
		if ( line_1.size() == 0 ) break;
		res_spline_map_.insert( std::pair< char, numeric::interpolation::spline::CubicSpline > (line_1[0], split_line_to_spline( line_1 )));
	}

	utility::io::izstream stream_2;
	basic::database::open( stream_2, "scoring/score_functions/MembranePotential/spline_test_cen.txt" );
	std::string line_2;

	TR << "parsing for cen" << std::endl;
	while ( stream_2  ) {
		getline( stream_2, line_2 );
		if ( line_2[0] == '#' ) continue;
		std::istringstream l_2(line_2);
		if ( line_2.size() == 0 ) break;
		cen_res_spline_map_.insert( std::pair< char, numeric::interpolation::spline::CubicSpline > (line_2[0], split_line_to_spline( line_2 )));
	}

}

numeric::interpolation::spline::CubicSpline
MPResidueLipophilicityEnergy::split_line_to_spline(
	std::string line_
) const
{
	numeric::interpolation::spline::CubicSpline spline;
	utility::vector1< Real > result;
	utility::vector1< std::string > spl;

	spl = utility::string_split( line_, ' ');

	numeric::Size matrix_size = 86;
	numeric::MathVector<numeric::Real> input_values(matrix_size, core::Real( 0.0 ));
	for ( core::Size i = 2; i <= spl.size(); ++i ) {
		input_values( i-2 ) = boost::lexical_cast< Real >( spl[ i ] );
	}

	using namespace numeric::interpolation::spline;
	BorderFlag behavior = e_Natural;
	numeric::Real start = -67.25;
	numeric::Real delta = 134.5/85.0;
	const std::pair<numeric::Real, numeric::Real> first_be = std::make_pair( 0, 0);

	spline.train(behavior, start, delta, input_values, first_be);

	return spline;
}

utility::vector1< Real >
MPResidueLipophilicityEnergy::split_line_to_polyval(
	std::string line_
) const
{
	// return every line as vector1 of 5 numbers, coefficients of highest power to lowest
	utility::vector1< Real > result;
	utility::vector1< std::string > spl;
	spl = utility::string_split( line_, ' ');
	result.push_back(boost::lexical_cast< Real >(spl[2]));
	result.push_back(boost::lexical_cast< Real >(spl[3]));
	result.push_back(boost::lexical_cast< Real >(spl[4]));
	result.push_back(boost::lexical_cast< Real >(spl[5]));
	result.push_back(boost::lexical_cast< Real >(spl[6]));

	return result;
}

// returns the result of the polynom of res at z
Real
MPResidueLipophilicityEnergy::spline_by_z(
	char const & res,
	Real const & z,
	bool & full_atom
) const
{
	Real x = 0.0;
	if ( z <= 50 && z >= -50 ) {
		if ( full_atom ) {
			x = res_spline_map_.at( res ).F(z);
		} else {
			x = cen_res_spline_map_.at( res ).F(z);
		}
	} else {
		return 0.0;
	}
	return x;
}

// returns the result of the polynom of res at z
Real
MPResidueLipophilicityEnergy::polynom_by_z(
	char const & res,
	Real const & z,
	bool & full_atom
) const
{
	Real x = 0.0;
	if ( z <= 15 && z >= -15 ) {
		if ( full_atom ) {
			x += res_polynom_map_.at(res)[1] * pow(z, 4.0);
			x += res_polynom_map_.at(res)[2] * pow(z, 3.0);
			x += res_polynom_map_.at(res)[3] * pow(z, 2.0);
			x += res_polynom_map_.at(res)[4] * z;
			x += res_polynom_map_.at(res)[5];
		} else {
			x += cen_res_polynom_map_.at(res)[1] * pow(z, 4.0);
			x += cen_res_polynom_map_.at(res)[2] * pow(z, 3.0);
			x += cen_res_polynom_map_.at(res)[3] * pow(z, 2.0);
			x += cen_res_polynom_map_.at(res)[4] * z;
			x += cen_res_polynom_map_.at(res)[5];
		}
	} else {
		return 0.0;
	}
	return x;
}

// returns the number of atoms closer to rsd than both cutoffs
utility::vector1 < core::Size >
MPResidueLipophilicityEnergy::neighboring_atoms(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Real const & cutoff_1,
	Real const & cutoff_2
) const
{
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
MPResidueLipophilicityEnergy::burial_sigmoid(
	core::Size const n_atoms,
	Real const slope,
	Real const offset
) const
{
	Real result = 1.0 / ( 1.0 + exp( ( Real( n_atoms ) - offset ) * slope ));
	if ( result < 0 ) result = 0.0;
	if ( result > 1 ) result = 1.0;

	return result;
}

core::Size
MPResidueLipophilicityEnergy::version() const
{
	return 1; // Initial versioning
}

} // membrane
} // scoring
} // core

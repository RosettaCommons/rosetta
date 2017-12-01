// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane/MPSpanInsertionEnergy.cc
/// @brief  penalize spans with unnatural dG of insertion
/// @details the distribution of dG of insertion to the membrane is quite narrow. rosetta will easily
/// design TMs with much more more negative dG of insertion. this energy term penalizes the pose for
/// TMs with a dG of insertion very differnet than the natural distribution. as calculated by the dsTbL
/// hydrophobicity scale
/// @author Jonathan Weinstein


//Unit headers
#include <core/scoring/membrane/MPSpanInsertionEnergy.hh>
#include <core/scoring/membrane/MPSpanInsertionEnergyCreator.hh>

//Package headers
#include <core/pose/Pose.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/interpolation/spline/CubicSpline.hh>
#include <boost/lexical_cast.hpp>
#include <map>

//utility headers
#include <utility/exit.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <basic/Tracer.hh>

//C++ headers
#include <iostream>
#include <utility/io/izstream.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>
#include <basic/database/open.hh>
#include <core/scoring/dssp/Dssp.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

static basic::Tracer TR( "core.scoring.membrane.MPSpanInsertionEnergy"  );

using namespace core::scoring;
using namespace core::scoring::methods;

namespace core {
namespace scoring {
namespace membrane {


/// @details This must return a fresh instance of the MPSpanInsertionEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MPSpanInsertionEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MPSpanInsertionEnergy );
}

ScoreTypes
MPSpanInsertionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( span_ins );
	return sts;
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
MPSpanInsertionEnergy::MPSpanInsertionEnergy() :
	parent( methods::EnergyMethodCreatorOP( new MPSpanInsertionEnergyCreator ) )
{
	parse_elazar_spline_table();
	// these values were calculted from the Rost dataset of PDBTM/OPM defiend topologies.
	// the dG of indertion of all TMs in the dataset was calulated (from sequence only, z by seq position)
	// using the same splines as are used by this energy term.
	// the script that makes the calculation is called retrive_natural_TMs_scores.py
	avg_dg_span_single_ = -5.71;
	std_dg_span_single_ = 5.03;

	avg_dg_span_multi_ = -3.58;
	std_dg_span_multi_ = 4.52;
}



//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
void
MPSpanInsertionEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	core::Real result = compute( const_cast< const pose::Pose &>( pose ) );
	totals[ span_ins ] = result;
}

/// @brief Create a clone of this energy method
EnergyMethodOP
MPSpanInsertionEnergy::clone() const
{
	return EnergyMethodOP( new MPSpanInsertionEnergy( *this ) );
}

void
MPSpanInsertionEnergy::setup_for_derivatives(
	pose::Pose &,
	ScoreFunction const &
) const
{
}


void
MPSpanInsertionEnergy::eval_atom_derivative(
	id::AtomID const & aid,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const &,
	Vector &,
	Vector &
) const {
	conformation::Residue const & rsd = pose.residue( aid.rsd() );
	if ( ! rsd.is_protein() ) return;
	Vector f1(0.0), f2(0.0);
}

core::Real
MPSpanInsertionEnergy::compute( pose::Pose const & pose ) const
{
	utility::vector1< core::conformation::membrane::Span > spans = create_updated_span( pose );

	utility::vector1< core::Real > span_scores;
	core::Real span_scores_sum = 0.0;
	span_scores.empty();
	core::Real temp( 0.0 ); core::Real p( 0.0 );

	core::Real _avg( 0.0 ); core::Real _std( 0.0 );
	if ( spans.size() == 1 ) {
		_avg = avg_dg_span_single_;
		_std = std_dg_span_single_;
	} else {
		_avg = avg_dg_span_multi_;
		_std = std_dg_span_multi_;
	}
	for ( core::Size span_i = 1; span_i <= spans.size(); ++span_i  ) {
		temp = calc_span_score( pose, spans[ span_i  ].start(), spans[ span_i  ].end() );
		p = exp( pow( temp - _avg, 2 ) / ( 2 * _std * _std ) );
		TR.Debug << "for span " << spans[span_i].start() << "-" << spans[span_i].end() << " = " << temp << " " << p << std::endl;
		span_scores.push_back( p );
		span_scores_sum += p;
	}
	TR.Debug << "the sum for all spans is " << span_scores_sum << std::endl;
	return( span_scores_sum );
}

//////////////////////////////////////////////////////
//@brief find all residues that are in the membrane (within membrane_thickness from the membrane mid-plane)
//////////////////////////////////////////////////////
utility::vector1< core::conformation::membrane::Span >
MPSpanInsertionEnergy::create_updated_span( pose::Pose const & pose ) const
{
	core::Real thickness( pose.conformation().membrane_info()->membrane_thickness() );
	utility::vector1< bool > span_vec;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( !pose.residue( i ).is_protein() ) continue;
		if ( pose.residue( i ).nbr_atom_xyz().z() >= -thickness &&
				pose.residue( i ).nbr_atom_xyz().z() <= thickness ) {
			span_vec.push_back( true );
		} else {
			span_vec.push_back( false );
		}
	}
	bool open_tm( false );
	core::Size chain( pose.chain( 1 ) );
	core::Size start( 0 );
	core::conformation::membrane::Orientation orient( core::conformation::membrane::out );

	utility::vector1< core::conformation::membrane::Span > all_spans;

	for ( core::Size i = 1; i <= span_vec.size(); ++i ) {
		if ( pose.chain( i ) != chain ) { // moved chains
			if ( open_tm ) {
				if ( pose.residue( start ).nbr_atom_xyz().z() < 0 && pose.residue( i - 1 ).nbr_atom_xyz().z() > 0 ) {
					orient = core::conformation::membrane::in;
				} else if ( pose.residue( start ).nbr_atom_xyz().z() > 0 && pose.residue( i - 1 ).nbr_atom_xyz().z() < 0 ) {
					orient = core::conformation::membrane::out;
				} else {
					TR.Debug << "for span " << start << "-" << i-1 << " orientation not clear, skipping..." << std::endl;
					open_tm = false;
					continue;
					//orient = core::conformation::membrane::in;
					//throw utility::excn::EXCN_Msg_Exception( "check the pose" );
				}
				all_spans.push_back( core::conformation::membrane::Span( start, i-1, orient ) );
				open_tm = false;
			}
			chain = pose.chain( i );
		}// else {
		if ( span_vec[ i ] && ! open_tm ) { // found new TM
			start = i;
			open_tm = true;
		} else if ( !span_vec[ i ] && open_tm ) { // reached end of TM
			if ( pose.residue( start ).nbr_atom_xyz().z() < 0 && pose.residue( i - 1 ).nbr_atom_xyz().z() > 0 ) {
				orient = core::conformation::membrane::in;
			} else if ( pose.residue( start ).nbr_atom_xyz().z() > 0 && pose.residue( i - 1 ).nbr_atom_xyz().z() < 0 ) {
				orient = core::conformation::membrane::out;
			} else {
				TR.Debug << "for span " << start << "-" << i-1 << " orientation not clear, skipping..." << std::endl;
				open_tm = false;
				continue;
				//orient = core::conformation::membrane::in;
				//throw utility::excn::EXCN_Msg_Exception( "check the pose"  );
			}
			all_spans.push_back( core::conformation::membrane::Span( start, i-1, orient ) );
			open_tm = false;
		}
		//}
	}
	bool add_tm( false );
	if ( open_tm ) {
		core::Size i = span_vec.size()+1;
		if ( pose.residue( start ).nbr_atom_xyz().z() < 0 && pose.residue( i - 1 ).nbr_atom_xyz().z() > 0 ) {
			orient = core::conformation::membrane::in;
			add_tm = true;
		} else if ( pose.residue( start ).nbr_atom_xyz().z() > 0 && pose.residue( i - 1 ).nbr_atom_xyz().z() < 0 ) {
			orient = core::conformation::membrane::out;
			add_tm = true;
		} else {
			TR.Debug << "for span " << start << "-" << i-1 << " orientation not clear, skipping..." << std::endl;
		}
		if ( add_tm ) all_spans.push_back( core::conformation::membrane::Span( start, i-1, orient ) );
	}
	TR.Debug << "the updated span list is:" << std::endl;
	for ( core::Size i = 1; i <= all_spans.size(); ++i ) {
		TR.Debug << "span #" << i << " " << all_spans[ i ].start() << "-" << all_spans[ i ].end() << " ~ " << all_spans[ i ].orientation() << std::endl;
	}
	return( all_spans );
}

core::Real
MPSpanInsertionEnergy::calc_span_score( pose::Pose const & pose, core::Size start, core::Size end ) const
{
	core::Real dg_span = 0.0;
	for ( core::Size i = start; i <= end; ++i ) {
		core::conformation::Residue rsd( pose.residue( i ) );
		dg_span += res_spline_map_.at( oneletter_code_from_aa(rsd.aa()) ).F( pose.conformation().membrane_info()->residue_z_position( pose.conformation(), rsd.seqpos() ) );
	}
	return( dg_span );
}

void
MPSpanInsertionEnergy::parse_elazar_spline_table()
{
	utility::io::izstream stream_1;
	basic::database::open( stream_1, "scoring/score_functions/MembranePotential/elazar_spline_mp_span_ins_fa.txt" );
	std::string line_1;

	while ( stream_1  ) {
		getline( stream_1, line_1 );
		if ( line_1[0] == '#' ) continue;
		std::istringstream l_1(line_1);
		if ( line_1.size() == 0 ) break;
		res_spline_map_.insert( std::pair< char, numeric::interpolation::spline::CubicSpline > (line_1[0], split_line_to_spline( line_1 )));
	}
}

numeric::interpolation::spline::CubicSpline
MPSpanInsertionEnergy::split_line_to_spline(
	std::string line
) const
{
	numeric::interpolation::spline::CubicSpline spline;
	utility::vector1< Real > result;
	result.empty();
	utility::vector1< std::string > spl;
	spl.empty();

	spl = utility::string_split( line, ' ');

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

core::Real
MPSpanInsertionEnergy::spline_by_z(
	char const & res,
	core::Real const & z
) const
{
	core::Real x = 0.0;
	if ( z <= 50 && z >= -50 ) {
		x = res_spline_map_.at( res ).F(z);
	} else {
		return 0.0;
	}
	return x;
}

core::Size
MPSpanInsertionEnergy::version() const
{
	return 1; // Initial versioning
}


} // membrane
} // scoring
} // core

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane/MPSpanAngleEnergy.cc
/// @brief  penalize spans with unnatural dG of insertion
/// @author Jonathan Weinstein


//Unit headers
#include <core/scoring/membrane/MPSpanAngleEnergy.hh>
#include <core/scoring/membrane/MPSpanAngleEnergyCreator.hh>

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

// term specific headers
#include <numeric/conversions.hh>
#include <basic/svd/SVD_Solver.hh>
#include <numeric/PCA.hh>

static basic::Tracer TR( "core.scoring.membrane.MPSpanAngleEnergy"  );

using namespace core::scoring;
using namespace core::scoring::methods;

namespace core {
namespace scoring {
namespace membrane {


/// @details This must return a fresh instance of the MPSpanAngleEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MPSpanAngleEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MPSpanAngleEnergy );
}

ScoreTypes
MPSpanAngleEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( mp_span_ang );
	return sts;
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
MPSpanAngleEnergy::MPSpanAngleEnergy() :
	parent( methods::EnergyMethodCreatorOP( new MPSpanAngleEnergyCreator ) )
{
	core::scoring::membrane::MPSpanInsertionEnergy mp_span_ins = core::scoring::membrane::MPSpanInsertionEnergy();
}



//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
void
MPSpanAngleEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	utility::vector1< core::Real > result = compute( const_cast< const pose::Pose &>( pose ), TR, false );
	core::Real total = 0;
	for ( auto const & sc : result ) total += sc;
	totals[ mp_span_ang ] = total;
}

EnergyMethodOP
MPSpanAngleEnergy::clone() const
{
	return EnergyMethodOP( new MPSpanAngleEnergy( *this ) );
}

void
MPSpanAngleEnergy::setup_for_derivatives(
	pose::Pose &,
	ScoreFunction const &
) const
{
}


void
MPSpanAngleEnergy::setup_for_scoring(
	pose::Pose &,
	ScoreFunction const &
) const
{
}


void
MPSpanAngleEnergy::eval_atom_derivative(
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


utility::vector1< core::Real >
MPSpanAngleEnergy::compute( pose::Pose const & pose, std::ostream & out, bool report ) const
{

	utility::vector1< core::conformation::membrane::Span > spans = mp_span_ins_.create_updated_span( pose );

	utility::vector1< core::Real > span_scores;
	utility::vector1< core::Real > span_scores_sum;
	span_scores.empty();

	numeric::xyzVector < core::Real > memb_normal( 0, 0, 1 );

	for ( core::Size span_i = 1; span_i <= spans.size(); ++span_i  ) {
		utility::vector1< numeric::xyzVector < core::Real > > span_vec( find_helix_vector( pose, spans[ span_i ].start(), spans[ span_i ].end() ) );
		numeric::xyzVector < core::Real > span_vec_norm( span_vec[2] - span_vec[1] );
		span_vec_norm.normalize_or_zero();
		core::Real angle = numeric::conversions::degrees( angle_of( span_vec_norm, memb_normal ) );

		if ( angle >= 90 ) angle = std::abs( angle - 180 );

		span_scores_sum.push_back( calc_ang_score( angle ) );
		if ( report ) {
			out << "span " << span_i << " " << spans[ span_i ].start() << " " << spans[ span_i ].end() << std::endl;
			out << "vec0 " << span_vec[1].x() << " " << span_vec[1].y() << " " << span_vec[1].z() << std::endl;
			out << "vec1 " << span_vec[2].x() << " " << span_vec[2].y() << " " << span_vec[2].z() << std::endl;
			out << "nrm " << span_vec_norm.x() << " " << span_vec_norm.y() << " " << span_vec_norm.z() << std::endl;
			out << "found angle " << angle << " scored " << calc_ang_score( angle ) << std::endl;
			span_scores_sum.push_back( angle );
		}
	}
	return( span_scores_sum );
}


utility::vector1< numeric::xyzVector < core::Real > >
MPSpanAngleEnergy::find_helix_vector( core::pose::Pose const & pose, core::Size start, core::Size end ) const
{
	utility::vector1< numeric::xyzVector <core::Real > > vec;

	utility::vector1< numeric::xyzVector< core::Real> > bb_coords; bb_coords.clear();

	for ( core::Size i=start; i <= end; i++ ) {
		core::conformation::Residue rsd( pose.residue( i ) );
		for ( core::Size j=1; j <= rsd.nheavyatoms(); ++j ) {
			if ( rsd.atom_is_backbone( j ) ) {
				core::Vector const & bbvec( rsd.xyz( j ) );
				bb_coords.push_back( bbvec );
			}
		}
	}
	numeric::xyzVector<core::Real> com = center_of_mass( bb_coords );
	numeric::xyzVector<core::Real> first_principal_component = numeric::first_principal_component( bb_coords );

	numeric::xyzVector<core::Real> com_principal_component = first_principal_component += com;

	//p0 is point on principal component vector closest to N-term
	numeric::xyzVector<core::Real> p0 = numeric::closest_point_on_line( com, com_principal_component, bb_coords[1]);

	//p1 is point on principal component vector closest to C-term
	numeric::xyzVector<core::Real> p1 = numeric::closest_point_on_line( com, com_principal_component, bb_coords.back());
	//TR << "com " << com.x() << " " << com.y() << " " << com.z() << std::endl;
	//TR << "p0 " << p0.x() << " " << p0.y() << " " << p0.z() << std::endl;
	//TR << "p1 " << p1.x() << " " << p1.y() << " " << p1.z() << std::endl;
	vec.push_back( p0 );
	vec.push_back( p1 );
	return( vec );
}

core::Real
MPSpanAngleEnergy::calc_ang_score( core::Real ang ) const
{
	return ( 0.00015 * std::pow(ang, 3) ) - ( 0.0088 * std::pow(ang, 2) ) + ( 0.18 * ang ) - 0.5;
}


core::Size
MPSpanAngleEnergy::version() const
{
	return 1; // Initial versioning
}


} // membrane
} // scoring
} // core

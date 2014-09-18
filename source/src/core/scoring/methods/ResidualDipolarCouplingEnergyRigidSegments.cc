// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ResidualDipolarCouplingEnergyRigidSegments.cc
/// @brief  RDC energy - comparing experimental RDC values to calculated values
/// @author Srivatsan Raman


//Unit headers
#include <core/scoring/methods/ResidualDipolarCouplingEnergyRigidSegments.hh>
#include <core/scoring/methods/ResidualDipolarCouplingEnergyRigidSegmentsCreator.hh>
#include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/scoring/ResidualDipolarCoupling.fwd.hh>
#include <core/scoring/ScoreType.hh>
//Package headers

#include <core/conformation/Residue.hh>
//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/datacache/CacheableDataType.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
// AUTO-REMOVED #include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

#include <core/id/NamedAtomID.hh>

//utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

//Objexx headers
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/Fmath.hh>

#include <utility/io/ozstream.hh> //for dump_weights

#include <basic/options/option.hh>
#include <basic/options/keys/rdc.OptionKeys.gen.hh>

//C++ headers
#include <iostream>

//Auto Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>


static thread_local basic::Tracer tr( "core.scoring.ResidualDipolarCoupling" );

namespace core {
namespace scoring {
namespace methods {

using namespace ObjexxFCL::format;

/// @details This must return a fresh instance of the ResidualDipolarCouplingEnergyRigidSegments class,
/// never an instance already in use
methods::EnergyMethodOP
ResidualDipolarCouplingEnergyRigidSegmentsCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new ResidualDipolarCouplingEnergyRigidSegments;
}

ScoreTypes
ResidualDipolarCouplingEnergyRigidSegmentsCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rdc_segments );
	return sts;
}



//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
ResidualDipolarCouplingEnergyRigidSegments::ResidualDipolarCouplingEnergyRigidSegments() :
	parent( new ResidualDipolarCouplingEnergyRigidSegmentsCreator )
{}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
EnergyMethodOP
ResidualDipolarCouplingEnergyRigidSegments::clone() const
{

  return new ResidualDipolarCouplingEnergyRigidSegments();

}

void ResidualDipolarCouplingEnergyRigidSegments::setup_for_scoring(
  pose::Pose & pose,
  ScoreFunction const &
) const
{
	//ResidualDipolarCoupling& rdc_data( rdc_from_pose( pose ) );
	dip_score_ = eval_dipolar( pose );
}

void ResidualDipolarCouplingEnergyRigidSegments::finalize_total_energy(
  pose::Pose &,
  ScoreFunction const &,
  EnergyMap & totals
) const
{
	totals[ rdc ] = dip_score_;
}

void ResidualDipolarCouplingEnergyRigidSegments::setup_for_minimizing(
  pose::Pose & pose,
  ScoreFunction const &,
	optimization::MinimizerMap const &
) const
{
	ResidualDipolarCoupling const& rdc_data( * retrieve_RDC_segments_from_pose( pose ) );
	residual_dipolar_coupling::RDC_lines const& All_RDC_lines( rdc_data.get_RDC_data() );
	residual_dipolar_coupling::RDC_lines::const_iterator it;
	Size ct = 0;
	for( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it) {
		id::AtomID atom1( id::NamedAtomID( it->atom1(), it->res1() ), pose );
		id::AtomID atom2( id::NamedAtomID( it->atom2(), it->res2() ), pose );
		tr.Trace << "insert in atom-map " << atom1 << " " << atom2 << std::endl;
		++ct;
		atom2rdc_map_.set( atom1, ct );
		atom2rdc_map_.set( atom2, ct );
	}
}



//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
ResidualDipolarCoupling &
ResidualDipolarCouplingEnergyRigidSegments::rdc_from_pose(
	pose::Pose & pose
) const
{
// 	//using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA;

// 	if( pose.data().has( RESIDUAL_DIPOLAR_COUPLING_DATA ) )
// 		return *( static_cast< ResidualDipolarCoupling const * >( pose.data().get_const_ptr( RESIDUAL_DIPOLAR_COUPLING_DATA )() ) );

// 	ResidualDipolarCouplingOP rdc_info = new ResidualDipolarCoupling;
// 	pose.data().set( RESIDUAL_DIPOLAR_COUPLING_DATA, rdc_info );
 	ResidualDipolarCouplingOP rdc_info( retrieve_RDC_segments_from_pose( pose ) );
	if ( !rdc_info ) {
		rdc_info = new ResidualDipolarCouplingRigidSegments;
		store_RDC_segments_in_pose( rdc_info, pose );
	}
	return *rdc_info;
}

//////////////////////////////////////////////////////
//@brief main computation routine for RDC energy... everything is happening here right now.
// this has to be spread out over different routines to make this energy yield derivatives
//////////////////////////////////////////////////////
Real ResidualDipolarCouplingEnergyRigidSegments::eval_dipolar(
  pose::Pose & pose
) const
{

	ResidualDipolarCoupling& rdc_data( rdc_from_pose( pose ) );
	utility::vector1< core::scoring::RDC > const& All_RDC_lines( rdc_data.get_RDC_data() );
	Real score;
	//Size const nrow( All_RDC_lines.size() ); //number of experimental couplins
	if ( basic::options::option[ basic::options::OptionKeys::rdc::iterate_weights ].user() ) {
		Real const sigma2( basic::options::option[ basic::options::OptionKeys::rdc::iterate_weights ] );
		Real const tol( basic::options::option[ basic::options::OptionKeys::rdc::iterate_tol ] );
		bool const reset( basic::options::option[ basic::options::OptionKeys::rdc::iterate_reset ] );
		score = rdc_data.iterate_tensor_weights( pose, sigma2, tol, reset );

		if ( basic::options::option[ basic::options::OptionKeys::rdc::dump_weight_trajectory ].user() ) {
			std::string const filename( basic::options::option[ basic::options::OptionKeys::rdc::dump_weight_trajectory ]() );
			utility::io::ozstream out( filename , std::ios_base::out | std::ios_base::app );
			utility::vector1< core::scoring::RDC >::const_iterator it;
			for ( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it) {
				out << RJ( 4, it->weight()) << " ";
			}
			out << std::endl;
		} //dump_weights
	} else {
		score = rdc_data.compute_dipscore( pose );
	}
	return score;
}


void
ResidualDipolarCouplingEnergyRigidSegments::eval_atom_derivative(
			 id::AtomID const & aid,
			 pose::Pose const & pose,
    	 kinematics::DomainMap const &,
	   	 ScoreFunction const &,
		   EnergyMap const & score_weights,
		   Vector & F1,
		   Vector & F2
) const {

	if ( !atom2rdc_map_.has( aid ) ) return; //damn this "has" isn't correct at all
	Size const rdc_nr( atom2rdc_map_[ aid ] );
	if ( rdc_nr == 0 ) {
		//		tr.Trace << "no RDC entry for " << aid << " skipping.. "<< std::endl;
		return;
	}
	ResidualDipolarCoupling const& rdc_cache( *retrieve_RDC_segments_from_pose( pose ) );
	utility::vector1< core::scoring::RDC > All_RDC_lines( rdc_cache.get_RDC_data() );
	runtime_assert( rdc_nr <= All_RDC_lines.size() );
	RDC const& rdc_data( All_RDC_lines[ rdc_nr ] );
	conformation::Residue const& rsd1( pose.residue( rdc_data.res1() ) );
	conformation::Residue const& rsd2( pose.residue( rdc_data.res2() ) );
	Vector fij;
	if ( aid.rsd() == rdc_data.res1() && utility::trimmed_compare( rsd1.atom_name( aid.atomno() ), rdc_data.atom1() ) ) {
		fij = rdc_data.fij();
	} else if ( aid.rsd() == rdc_data.res2() && utility::trimmed_compare( rsd2.atom_name( aid.atomno() ), rdc_data.atom2() ) ){
		fij = -rdc_data.fij();
	} else return;

	//	Real fij;
	//thanks to Will Sheffler:
	numeric::xyzVector<core::Real> atom_x = pose.xyz(aid);
	numeric::xyzVector<core::Real> const f2( -fij );
	numeric::xyzVector<core::Real> const atom_y = atom_x - f2;   // a	"fake" atom in the direcion of the gradient
	numeric::xyzVector<core::Real> const f1( atom_x.cross( atom_y ) );

	F1 += score_weights[ rdc ] * f1;
	F2 += score_weights[ rdc ] * f2;

}

} // methods
} // scoring
} // core

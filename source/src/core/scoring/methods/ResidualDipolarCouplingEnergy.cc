// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ResidualDipolarCouplingEnergy.cc
/// @brief  RDC energy - comparing experimental RDC values to calculated values
/// @author Srivatsan Raman


//Unit headers
#include <core/scoring/methods/ResidualDipolarCouplingEnergy.hh>
#include <core/scoring/methods/ResidualDipolarCouplingEnergyCreator.hh>
#include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/scoring/ResidualDipolarCoupling.fwd.hh>
#include <core/scoring/ScoreType.hh>
//Package headers

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzVector.hh>


//utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

//Objexx headers
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/io/ozstream.hh> //for dump_weights

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rdc.OptionKeys.gen.hh>

//C++ headers
#include <iostream>

#include <core/scoring/EnergyMap.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>


static THREAD_LOCAL basic::Tracer tr( "core.scoring.ResidualDipolarCoupling" );

namespace core {
namespace scoring {
namespace methods {

using namespace ObjexxFCL::format;

/// @details This must return a fresh instance of the ResidualDipolarCouplingEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ResidualDipolarCouplingEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new ResidualDipolarCouplingEnergy );
}

ScoreTypes
ResidualDipolarCouplingEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rdc );
	return sts;
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
ResidualDipolarCouplingEnergy::ResidualDipolarCouplingEnergy() :
	parent( methods::EnergyMethodCreatorOP( new ResidualDipolarCouplingEnergyCreator ) )
{}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
EnergyMethodOP
ResidualDipolarCouplingEnergy::clone() const
{

	return EnergyMethodOP( new ResidualDipolarCouplingEnergy() );

}

void ResidualDipolarCouplingEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const
{
	//ResidualDipolarCoupling& rdc_data( rdc_from_pose( pose ) );
	dip_score_ = eval_dipolar( pose );
}

void ResidualDipolarCouplingEnergy::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	totals[ rdc ] = dip_score_;
}

void ResidualDipolarCouplingEnergy::setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &
) const
{
	ResidualDipolarCoupling const& rdc_data( * retrieve_RDC_from_pose( pose ) );
	ResidualDipolarCoupling::RDC_lines const& All_RDC_lines( rdc_data.get_RDC_data() );
	ResidualDipolarCoupling::RDC_lines::const_iterator it;
	Size ct = 0;
	for ( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it ) {
		id::AtomID atom1( pose.residue(it->res1()).atom_index(it->atom1()), it->res1());
		id::AtomID atom2( pose.residue(it->res2()).atom_index(it->atom2()), it->res2());
		tr.Trace << "insert in atom-map " << atom1 << " " << atom2 << std::endl;
		++ct;
		utility::vector1< core::Size > atm1_map = atom2rdc_map_.get( atom1 );
		utility::vector1< core::Size > atm2_map = atom2rdc_map_.get( atom2 );
		atm1_map.push_back( ct );
		atm2_map.push_back( ct );
		atom2rdc_map_.set( atom1, atm1_map );
		atom2rdc_map_.set( atom2, atm2_map );
	}
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
ResidualDipolarCoupling &
ResidualDipolarCouplingEnergy::rdc_from_pose(
	pose::Pose & pose
) const
{
	ResidualDipolarCouplingOP rdc_info( retrieve_RDC_from_pose( pose ) );
	if ( !rdc_info ) {
		rdc_info = ResidualDipolarCouplingOP( new ResidualDipolarCoupling );
		store_RDC_in_pose( rdc_info, pose );
	}
	return *rdc_info;
}

//////////////////////////////////////////////////////
//@brief main computation routine for RDC energy... everything is happening here right now.
// this has to be spread out over different routines to make this energy yield derivatives
//////////////////////////////////////////////////////
Real ResidualDipolarCouplingEnergy::eval_dipolar(
	pose::Pose& pose
) const
{
	ResidualDipolarCoupling& rdc_data( rdc_from_pose( pose ) );
	return eval_dipolar( pose, rdc_data );
}

Real ResidualDipolarCouplingEnergy::eval_dipolar(
	pose::Pose const& pose,
	ResidualDipolarCoupling& rdc_data
) const
{
	Real score;
	//Size const nrow( All_RDC_lines.size() ); //number of experimental couplins
	if ( basic::options::option[ basic::options::OptionKeys::rdc::iterate_weights ].user() ) {
		utility::vector1< core::scoring::RDC > const& All_RDC_lines( rdc_data.get_RDC_data() );
		Real const sigma2( basic::options::option[ basic::options::OptionKeys::rdc::iterate_weights ] );
		Real const tol( basic::options::option[ basic::options::OptionKeys::rdc::iterate_tol ] );
		bool const reset( basic::options::option[ basic::options::OptionKeys::rdc::iterate_reset ] );
		score = rdc_data.iterate_tensor_weights( pose, sigma2, tol, reset );

		if ( basic::options::option[ basic::options::OptionKeys::rdc::dump_weight_trajectory ].user() ) {
			std::string const filename( basic::options::option[ basic::options::OptionKeys::rdc::dump_weight_trajectory ]() );
			utility::io::ozstream out( filename , std::ios_base::out | std::ios_base::app );
			utility::vector1< core::scoring::RDC >::const_iterator it;
			for ( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it ) {
				out << RJ( 4, it->weight()) << " ";
			}
			out << std::endl;
		} //dump_weights
	} else {
		static std::string const fit_method(
			basic::options::option[ basic::options::OptionKeys::rdc::fit_method ]()
		);
		if ( fit_method == "svd" ) {
			tr.Trace << "residual-energy method chosen: 'svd' " << std::endl;
			score = rdc_data.compute_dipscore( pose );
		} else {
			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			if ( option[ OptionKeys::rdc::fixDa].user() ) {
				if ( option[ OptionKeys::rdc::fixR].user() ) {
					//Real const tensorDa( option[ OptionKeys::rdc::fixDa] );
					//Real const tensorR( option[ OptionKeys::rdc::fixR] );
					utility::vector1<Real> const tensorDa = option[ OptionKeys::rdc::fixDa]();
					utility::vector1<Real> const tensorR = option[ OptionKeys::rdc::fixR]();

					//make sure R is between 0 and 2/3
					for ( core::Size i = 1; i <= option[ OptionKeys::rdc::fixR ]().size(); ++i ) {
						if ( (tensorR[i] < 0) || (tensorR[i] > 2.0/3.0) ) {
							utility_exit_with_message("0=< R <=2/3");
						}
					}

					//make sure user provide the same number Da and R
					if (  ( option[ OptionKeys::rdc::fixDa ]().size() != option[ OptionKeys::in::file::rdc ]().size() )
							|| ( option[ OptionKeys::rdc::fixR ]().size() !=  option[ OptionKeys::in::file::rdc ]().size() )   )  {
						utility_exit_with_message("Number of Da and R must be the same as in number of experiment");
					}

					score = rdc_data.compute_dipscore_nlsDaR( pose, tensorDa , tensorR );
				} else { //end of Da and R
					//Real const tensorDa( option[ OptionKeys::rdc::fixDa] );
					utility::vector1<Real> const tensorDa = option[ OptionKeys::rdc::fixDa]();
					//tr.Trace << "nls and fixDa: " << tensorDa  << std::endl;
					score = rdc_data.compute_dipscore_nlsDa( pose, tensorDa);
				} //end of Da and noR
			} else { //end of fixDa
				if ( option[ OptionKeys::rdc::fixR].user() ) {
					//Real const tensorR( option[ OptionKeys::rdc::fixR] );
					utility::vector1<Real> const tensorR = option[ OptionKeys::rdc::fixR]();
					//tr.Trace << "nls and fixR: " << tensorR << std::endl;
					for ( core::Size i = 1; i <= option[ OptionKeys::rdc::fixR ]().size(); ++i ) {
						if ( (tensorR[i] < 0) || (tensorR[i] > 2.0/3.0) ) {
							utility_exit_with_message("0=< R <=2/3");
						}
					}
					score = rdc_data.compute_dipscore_nlsR( pose, tensorR);
				} else { //end of noDa and R
					//tr.Trace << "nls" << std::endl;
					score = rdc_data.compute_dipscore_nls( pose );
				}//end of noDa and noR
			}//end of no fixDa
		}//end of nls
	}//end of else dump
	return score;
}


void
ResidualDipolarCouplingEnergy::eval_atom_derivative(
	id::AtomID const & aid,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & score_weights,
	Vector & F1,
	Vector & F2
) const {

	if ( !atom2rdc_map_.has( aid ) ) return; //damn this "has" isn't correct at all
	utility::vector1< Size > const rdc_nrs( atom2rdc_map_[ aid ] );
	if ( rdc_nrs.size() == 0 ) {
		//  tr.Trace << "no RDC entry for " << aid << " skipping.. "<< std::endl;
		return;
	}
	Vector fij(0,0,0);
	for ( core::Size ii=1; ii<=rdc_nrs.size(); ++ii ) {
		core::Size rdc_nr = rdc_nrs[ ii ];
		ResidualDipolarCoupling const& rdc_cache( *retrieve_RDC_from_pose( pose ) );
		utility::vector1< core::scoring::RDC > All_RDC_lines( rdc_cache.get_RDC_data() );
		runtime_assert( rdc_nr <= All_RDC_lines.size() );
		RDC const& rdc_data( All_RDC_lines[ rdc_nr ] );
		conformation::Residue const& rsd1( pose.residue( rdc_data.res1() ) );
		conformation::Residue const& rsd2( pose.residue( rdc_data.res2() ) );
		if ( aid.rsd() == rdc_data.res1() && utility::trimmed_compare( rsd1.atom_name( aid.atomno() ), rdc_data.atom1() ) ) {
			fij += rdc_data.fij();
		} else if ( aid.rsd() == rdc_data.res2() && utility::trimmed_compare( rsd2.atom_name( aid.atomno() ), rdc_data.atom2() ) ) {
			fij -= rdc_data.fij();
		} else return;
	}

	// Real fij;
	//thanks to Will Sheffler:
	numeric::xyzVector<core::Real> atom_x = pose.xyz(aid);
	numeric::xyzVector<core::Real> const f2( -fij );
	numeric::xyzVector<core::Real> const atom_y = atom_x - f2;   // a "fake" atom in the direcion of the gradient
	numeric::xyzVector<core::Real> const f1( atom_x.cross( atom_y ) );

	F1 += score_weights[ rdc ] * f1;
	F2 += score_weights[ rdc ] * f2;
}

core::Size
ResidualDipolarCouplingEnergy::version() const
{
	return 1; // Initial versioning
}

} // methods
} // scoring
} // core

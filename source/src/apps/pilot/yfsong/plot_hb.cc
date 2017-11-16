// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/options/option.hh>

#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <ObjexxFCL/format.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>

#include <sstream>
#include <string>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


namespace plot_hb {
// mjo the new HBEvalTypes breaks this reporting mechanism!
basic::options::IntegerOptionKey hb_type("plot_hb:hb_type");

basic::options::RealOptionKey   dist("plot_hb:dist");
basic::options::RealOptionKey   m_cos_theta("plot_hb:m_cos_theta");
basic::options::RealOptionKey   m_cos_psi("plot_hb:m_cos_psi");
basic::options::BooleanOptionKey plot_dist("plot_hb:plot_dist");
basic::options::BooleanOptionKey plot_theta("plot_hb:plot_theta");
basic::options::BooleanOptionKey plot_psi("plot_hb:plot_psi");
// basic::options::BooleanOptionKey show_poly("plot_hb:show_poly");
}

static basic::Tracer tr( "apps.pilot.yfsong.plot_hb" );

int
main (int argc, char *argv[]){
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// what was SP3SC is now hbe_GENERIC_SP3SCSC_LR
		option.add( plot_hb::hb_type,     "hbond type" ).def(224);
		option.add( plot_hb::dist,        "hbond distance" ).def(1.9);
		option.add( plot_hb::m_cos_theta, "-cos(theta)" ).def(0.95);
		option.add( plot_hb::m_cos_psi,   "-cos(psi)" ).def(0.95);
		option.add( plot_hb::plot_dist,   "plot E vs distance" ).def(false);
		option.add( plot_hb::plot_theta,  "plot E vs -cos(theta)" ).def(false);
		option.add( plot_hb::plot_psi,    "plot E vs -cos(psi)" ).def(false);
		// option.add( plot_hb::show_poly,    "plot polynomials" ).def(false);
		devel::init( argc, argv );

		core::scoring::hbonds::HBEvalTuple hbe;

		if ( option[plot_hb::hb_type]() <= core::scoring::hbonds::hbe_MAX ) {
			hbe = core::scoring::hbonds::HBEvalType( option[plot_hb::hb_type]());
		} else {
			std::cerr << "Error! can't recognize hbond type: " << option[ plot_hb::hb_type ]() << std::endl;
			return -1;
		}

		core::Real dist = option[ plot_hb::dist ]();
		core::Real m_cos_theta = option[ plot_hb::m_cos_theta ]();
		core::Real m_cos_psi   = option[ plot_hb::m_cos_psi ]();
		core::Real dummy_chi = 0.0;
		core::Real energy = 0.0;
		core::Real dE_dr;
		core::Real dE_dxD;
		core::Real dE_dxH;
		HBondDatabaseCOP hb_database( HBondDatabase::get_database( "score12_params" ) );
		HBondOptions hboptions;

		if ( option[ plot_hb::plot_dist ]() ) {
			for ( dist = 0.0; dist < 3.01; dist += 0.1 ) {
				core::scoring::hbonds::hbond_compute_energy(*hb_database, hboptions, hbe, dist, m_cos_theta, m_cos_psi, dummy_chi, energy, dE_dr, dE_dxD, dE_dxH);
				tr << F(9,3,dist) << " " << F(9,3,m_cos_theta) << " " << F(9,3,m_cos_psi) << " " << F(9,3,energy) << std::endl;
			}
		} else if ( option[ plot_hb::plot_theta ]() ) {
			for ( m_cos_theta = -0.95; m_cos_theta < 0.99; m_cos_theta += 0.05 ) {
				core::scoring::hbonds::hbond_compute_energy(*hb_database, hboptions, hbe, dist, m_cos_theta, m_cos_psi, dummy_chi, energy, dE_dr, dE_dxD, dE_dxH);
				tr << F(9,3,dist) << " " << F(9,3,m_cos_theta) << " " << F(9,3,m_cos_psi) << " " << F(9,3,energy) << std::endl;
			}
		} else if ( option[ plot_hb::plot_psi ]() ) {
			for ( m_cos_psi = -0.95; m_cos_psi < 0.99; m_cos_psi += 0.05 ) {
				core::scoring::hbonds::hbond_compute_energy(*hb_database, hboptions, hbe, dist, m_cos_theta, m_cos_psi, dummy_chi, energy, dE_dr, dE_dxD, dE_dxH);
				tr << F(9,3,dist) << " " << F(9,3,m_cos_theta) << " " << F(9,3,m_cos_psi) << " " << F(9,3,energy) << std::endl;
			}
		} else {
			core::scoring::hbonds::hbond_compute_energy(*hb_database, hboptions, hbe, dist, m_cos_theta, m_cos_psi, dummy_chi, energy, dE_dr, dE_dxD, dE_dxH);
			tr << F(9,3,dist) << " " << F(9,3,m_cos_theta) << " " << F(9,3,m_cos_psi) << " " << F(9,3,energy) << std::endl;
		}

		// if ( option[ plot_hb::show_poly ]() ) core::scoring::hbonds::show_poly();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

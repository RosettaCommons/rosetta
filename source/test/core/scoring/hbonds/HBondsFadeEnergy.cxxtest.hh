// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBonds.cxxtest.hh
/// @brief  Test evaluation of hbond potentials for each HBEvalType across parameter ranges.
/// @author Matthew O'Meara

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/core/init_util.hh>

// Package headers
#include <core/scoring/hbonds/hbonds_geom.hh>

static basic::Tracer TR("core.scoring.hbonds.HBondsFadeEnergy");

class HBondsFadeEnergy : public CxxTest::TestSuite {

public:
  void setUp() {
		core_init();
	}

	void tearDown(){}


	void test_fade_energy() {
		using namespace core;
		using namespace scoring;
    using namespace hbonds;

		for(Real raw_energy=-.2; raw_energy < .2; raw_energy +=.05){
			Real energy(raw_energy);
			fade_energy(energy);
			TR
				<< raw_energy << "\t" << energy << std::endl;
		}

		//library(ggplot2)
		//z <- read.table("/tmp/data.csv")
		//names(z) <- c("RawEnergy", "FadedEnergy")
		//p <- ggplot(data=z) + theme_bw() +
		//  geom_line(aes(x=RawEnergy, y=RawEnergy)) +
		//  geom_line(aes(x=RawEnergy, y=FadedEnergy)) +
		//  ggtitle("Raw vs Faded Energy as computed by unit test.")
		//ggsave("raw_faded_energy.pdf")
	}
	void test_fade_energy_derivs() {
		using namespace core;
		using namespace scoring;
    using namespace hbonds;

		Real dE_dr(1.0), dE_dxD(1.0), dE_dxH(1.0), dE_dBAH(1.0), dE_dchi(1.0);
		Real delta( 1e-5 );
		for ( Real raw_energy = -.20004; raw_energy < .2; raw_energy +=.05 ) {
			Real energy(raw_energy);
			dE_dr = dE_dxD = dE_dxH = dE_dBAH = dE_dchi = 1.0;
			fade_energy(energy, dE_dr, dE_dxD, dE_dxH, dE_dBAH, dE_dchi);
			Real elow  = raw_energy - delta; fade_energy( elow );
			Real ehigh = raw_energy + delta; fade_energy( ehigh );
			Real numeric_deriv = (ehigh-elow)/(2*delta);
			//std::cout << "fade energy: " << raw_energy << " -> " << energy << "; deriv: " << numeric_deriv << " " << dE_dr << std::endl;
			TS_ASSERT_DELTA( numeric_deriv, dE_dr,   1e-6 );
			TS_ASSERT_DELTA( numeric_deriv, dE_dxD,  1e-6 );
			TS_ASSERT_DELTA( numeric_deriv, dE_dxH,  1e-6 );
			TS_ASSERT_DELTA( numeric_deriv, dE_dBAH, 1e-6 );
			TS_ASSERT_DELTA( numeric_deriv, dE_dchi, 1e-6 );
		}
	}
};

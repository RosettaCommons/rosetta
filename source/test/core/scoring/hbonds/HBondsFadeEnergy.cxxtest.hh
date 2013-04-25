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
	
};

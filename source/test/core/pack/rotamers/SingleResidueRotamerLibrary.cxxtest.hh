// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleResidueRotamerLibrary.cxxtest.hh
/// @brief
/// @author Rocco Moretti

// Test headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/pack/task/ResidueLevelTask_.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("core.pack.rotamers.SingleResidueRotamerLibrary.cxxtest");

using namespace core;
using namespace core::pack;
using namespace core::pack::task;

class DummySRRL: public core::pack::rotamers::SingleResidueRotamerLibrary {

public:
	virtual
	core::Real
	rotamer_energy_deriv(
		conformation::Residue const &,
		dunbrack::RotamerLibraryScratchSpace &
	) const {
		utility_exit_with_message("UNIMPLEMENTED");
	}

	virtual
	core::Real
	rotamer_energy(
		conformation::Residue const &,
		dunbrack::RotamerLibraryScratchSpace &
	) const { utility_exit_with_message("UNIMPLEMENTED"); }

	virtual
	core::Real
	best_rotamer_energy(
		conformation::Residue const &,
		bool,
		dunbrack::RotamerLibraryScratchSpace &
	) const { utility_exit_with_message("UNIMPLEMENTED"); }

	virtual
	void
	assign_random_rotamer_with_bias(
		conformation::Residue const & ,
		pose::Pose const &,
		dunbrack::RotamerLibraryScratchSpace &,
		numeric::random::RandomGenerator &,
		dunbrack::ChiVector &,
		bool
	) const { utility_exit_with_message("UNIMPLEMENTED"); }

	virtual
	void
	fill_rotamer_vector(
		pose::Pose const &,
		scoring::ScoreFunction const &,
		pack::task::PackerTask const &,
		utility::graph::GraphCOP,
		chemical::ResidueTypeCOP,
		conformation::Residue const&,
		utility::vector1< utility::vector1< Real > > const &,
		bool,
		core::pack::rotamers::RotamerVector &
	) const { utility_exit_with_message("UNIMPLEMENTED"); }

	virtual
	void
	write_to_file( utility::io::ozstream & ) const {
		utility_exit_with_message("UNIMPLEMENTED");
	}

};

class DummyRLT: public core::pack::task::ResidueLevelTask_ {

public:
	DummyRLT(int chi) { chi_ = chi; }

	void chi(int setting) { chi_ = setting; }

	virtual
	ExtraRotSample
	extrachi_sample_level(
		bool buried,
		int chi,
		chemical::ResidueType const &
	) const {
		if ( buried && chi <= chi_ ) {
			//std::cout << " Buried ex " << chi << std::endl;
			return EX_ONE_STDDEV;
		} else {
			//std::cout << " No Exp " << chi << std::endl;
			return NO_EXTRA_CHI_SAMPLES;
		}
	}

private:
	int chi_;

};

class SingleResidueRotamerLibraryTests : public CxxTest::TestSuite
{
public:

	core::chemical::ResidueTypeCOP hydroxyls_;

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-extra_res_fa core/pack/rotamers/6SA.params" );
		hydroxyls_ = core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)->name_map("6SA").get_self_ptr();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_extra_chi_expansion() {
		DummySRRL rotlib;
		DummyRLT rlt( 10 );

		utility::vector1< utility::vector1< core::Real > > chi_samples;

		// With no extra expansion.
		chi_samples = rotlib.compute_proton_chi_samplings( *hydroxyls_, rlt, false );

		TS_ASSERT_EQUALS( chi_samples.size(), hydroxyls_->n_proton_chi() );
		for ( core::Size ii(1); ii <= hydroxyls_->n_proton_chi(); ++ii ) {
			TS_ASSERT_EQUALS( chi_samples[ii].size(), 3 ); // Should be 60, -60, 180
		}

		// With extra expansion.
		chi_samples = rotlib.compute_proton_chi_samplings( *hydroxyls_, rlt, true );

		TS_ASSERT_EQUALS( chi_samples.size(), hydroxyls_->n_proton_chi() );
		for ( core::Size ii(1); ii <= hydroxyls_->n_proton_chi(); ++ii ) {
			//std::cout << " CHI " << ii << " samples " << chi_samples[ii].size() << std::endl;
			if ( ii <= 10 ) { // The first 19 chis in the params file are proton chis, so we don't have to worry about expansion
				TS_ASSERT_EQUALS( chi_samples[ii].size(), 9 ); // Should be 60, 40, 80, -60, -60, -40, -80, 180, 160, 200
			} else {
				TS_ASSERT_EQUALS( chi_samples[ii].size(), 3 ); // Should be 60, -60, 180
			}
		}

	}

	void test_expand_proton_chis() {
		DummySRRL rotlib;

		utility::vector1< core::Real > vect;
		vect.push_back( 60.0 );

		TS_ASSERT_EQUALS( hydroxyls_->n_proton_chi(), 19 );

		utility::vector1< utility::vector1< core::Real > > sampling( hydroxyls_->n_proton_chi(), vect );

		utility::vector1< dunbrack::ChiSetOP > chisets;

		chisets = rotlib.expand_proton_chis( sampling, *hydroxyls_ );
		TS_ASSERT_EQUALS( chisets.size(), 1 );

		for ( core::Size ii(1); ii <= sampling.size(); ii += 2 ) { // Should be 10
			sampling[ii].push_back( -60.0 );
		}

		chisets = rotlib.expand_proton_chis( sampling, *hydroxyls_, 50000 );
		TS_ASSERT_EQUALS( chisets.size(), 1024 );

		chisets = rotlib.expand_proton_chis( sampling, *hydroxyls_, 500 );
		TS_ASSERT_EQUALS( chisets.size(), 500 );

		// // Test to see if we bog down the machinery with the full expansion.
		// // This is turned off by default, as it takes a while
		// DummyRLT rlt( 0 );
		// chisets = rotlib.expand_proton_chis( rotlib.compute_proton_chi_samplings( *hydroxyls_, rlt, false ), *hydroxyls_, 50000 );

		// TS_ASSERT_EQUALS( chisets.size(), 50000 );
	}

};


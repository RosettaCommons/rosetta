// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/AmbiguousNMRDistanceConstraint.cxxtest.hh
/// @brief  test suite for AmbiguousNMRDistanceConstraint
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Clay Tydings (claiborne.w.tydings@vanderbilt.edu)

//Andrew did the serialization.
//Clay generalized the constraints and wrote unit tests to show they get the same result.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh> // DO NOT AUTO-REMOVE

//protocol headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>

//core headers
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>

// basic headers
#include <algorithm>
#include <cctype>

// unit headers

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class AmbiguousNMRDistanceConstraintTests : public CxxTest::TestSuite
{
private:
	std::string const test_label_ = "PEPTIDE_STUB_MOVER_TEST_LABEL";
public:

	void setUp(){
		core_init();
	}

	void tearDown(){

	}

	void test_ala() {
		//typedef utility::vector1< core::id::AtomID > Atoms;
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ALA", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		//utility::vector1< core::id::AtomID > atoms_old;
		//utility::vector1< core::id::AtomID > atoms_new;
		//Atoms atoms_new;
		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;


		core::Size res (1);
		//core::chemical::ResidueType type;
		//type(pose.residue_type(res));
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		core::scoring::constraints::parse_NMR_name("HA", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HA", res, pose, atoms_new);

		core::scoring::constraints::parse_NMR_name("H", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("H", res, pose, atoms_new);

		core::scoring::constraints::parse_NMR_name("1HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("1HB", res, pose, atoms_new);

		core::scoring::constraints::parse_NMR_name("2HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("2HB", res, pose, atoms_new);

		core::scoring::constraints::parse_NMR_name("3HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("3HB", res, pose, atoms_new);

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_gly() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "GLY", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HA", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HA", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QA", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QA", res, pose, atoms_new);

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_ile() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ILE", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HD1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HD1", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD1", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QQG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QQG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG1", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG1", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_leu() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "LEU", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HD1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HD1", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD1", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HD2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HD2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD2", res, pose, atoms_new);
		//The old function does not work for HD, just returns HD
		//core::scoring::constraints::parse_NMR_name("HD", res, aa, atoms_old);
		//core::scoring::constraints::parse_NMR_name_general("HD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QQD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QQD", res, pose, atoms_new);

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_val() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "VAL", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HG1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG1", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG1", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QQG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QQG", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_pro() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "PRO", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_cys() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "CYS", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_met() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "MET", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_phe() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "PHE", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_tyr() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "TYR", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_trp() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "TRP", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HZ", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HZ", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QZ", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QZ", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_his() {
		//typedef utility::vector1< core::id::AtomID > Atoms;
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "HIS", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos  && atm_name != "HD" && atm_name != "HE" ) {
				//for exclusion, see:
				//parse_NMR_name( std::string name, core::Size res, AmbiguousNMRDistanceConstraint::Atoms& atoms, core::pose::Pose const& pose )
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_arg() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ARG", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HH1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HH1", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QH1", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QH1", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HH2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HH2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QH2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QH2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HH", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HH", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QQH", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QQH", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_lys() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "LYS", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HZ", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HZ", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QZ", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QZ", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_glu() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "GLU", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_asp() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ASP", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_gln() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "GLN", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HE", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HE", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HE2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HE2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QE2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QE2", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_asn() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "ASN", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("HD2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HD2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QD2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QD2", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_ser() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "SER", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HB", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QB", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QB", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_thr() {
		core::pose::Pose pose( (core::pose::Pose()) );

		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(false);
		stubmover.add_residue( "Append", "THR", 0, false, "", 0, pose.total_residue(), nullptr, "" );
		stubmover.set_residue_label( test_label_ );
		stubmover.apply(pose);

		core::scoring::constraints::NamedAtoms atoms_old, atoms_new;
		core::Size res (1);
		core::chemical::AA const aa( pose.residue_type( res ).aa() );

		for ( core::Size num = 1; num <= pose.residue_type(res).natoms(); ++num ) {
			auto atm_name = pose.residue_type(res).atom_name(num);
			atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
			if ( atm_name.find('H') != std::string::npos ) {
				core::scoring::constraints::parse_NMR_name(atm_name, res, aa, atoms_old);
				core::scoring::constraints::parse_NMR_name_general(atm_name, res, pose, atoms_new);
			}
		}

		//        for (auto atom : atoms_old) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "Old: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }
		//        for (auto atom : atoms_new) {
		//            auto atm_name = atom.atom();
		//            //atm_name.erase(std::remove_if(atm_name.begin(), atm_name.end(), ::isspace), atm_name.end());
		//            std::cout << "New: " << atom.to_string() << " Atom:" << atm_name << "." << std::endl;
		//        }

		TS_ASSERT( atoms_new == atoms_old );

		//Ambiguouous constraints
		core::scoring::constraints::parse_NMR_name("HG2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("HG2", res, pose, atoms_new);
		core::scoring::constraints::parse_NMR_name("QG2", res, aa, atoms_old);
		core::scoring::constraints::parse_NMR_name_general("QG2", res, pose, atoms_new);

		std::sort(atoms_old.begin(), atoms_old.end());
		std::sort(atoms_new.begin(), atoms_new.end());

		TS_ASSERT( atoms_new == atoms_old );
	}

	void test_serialize_AmbiguousNMRDistanceConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		AmbiguousNMRDistanceConstraint::Atoms atoms1, atoms2;
		atoms1.push_back( AtomID( 1, 2 )); atoms1.push_back( AtomID( 1, 10 ) );
		atoms2.push_back( AtomID( 2, 3 )); atoms1.push_back( AtomID( 4, 16 ) );

		ConstraintOP instance( new AmbiguousNMRDistanceConstraint( atoms1, atoms2, some_func ) );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		ConstraintOP instance2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		TS_ASSERT( utility::pointer::dynamic_pointer_cast< AmbiguousNMRDistanceConstraint > ( instance2 ) );
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};


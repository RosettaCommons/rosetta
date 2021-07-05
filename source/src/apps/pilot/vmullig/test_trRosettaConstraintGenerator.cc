// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_trRosettaConstraintGenerator.cc
/// @brief A unit test (albeit one in the integration test suite) for the
/// trRosettaConstraintGenerator.  This compares the constraints generated
/// by the C++ ConstraintGenerator to known outputs from the original
/// Python version of trRosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// core headers

// protocol headers

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

#ifdef USE_TENSORFLOW
#include <protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGenerator.hh>
#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/SplineFunc.hh>
#include <basic/options/option_macros.hh>
#include <basic/citation_manager/CitationManager.hh>
#else
#include <basic/tensorflow_manager/util.hh>
#endif //ifndef USE_TENSORFLOW


static basic::Tracer TR("apps.pilot.vmullig.test_trRosettaConstraintGenerator");

#define FLOAT_COMPARISON_THRESHOLD 0.002

/// @brief Indicate which commandline flags are relevant to this application.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( in::file::fasta );

}

#ifdef USE_TENSORFLOW

OPT_KEY( String, csts_file )
OPT_KEY( String, msa_file )

/// @brief Print out a description of two constraints.
std::string
summarize_csts (
	core::scoring::constraints::ConstraintCOP const &cst1,
	core::scoring::constraints::ConstraintCOP const &cst2
) {
	using namespace core::scoring::constraints;
	std::ostringstream ss;
	for ( core::Size i(1); i<=2; ++i ) {
		Constraint const & cst( i==1 ? *cst1 : *cst2 );
		if ( i == 1 ) {
			ss << "\n1:\t";
		} else {
			ss << "\n2:\t";
		}
		cst.show( ss );
	}

	return ss.str();
}

/// @brief Are two distance constraints approximately equal?
bool
dist_csts_approx_equal(
	core::scoring::constraints::AtomPairConstraint const & cst1,
	core::scoring::constraints::AtomPairConstraint const & cst2
) {
	if ( cst1.atom1() != cst2.atom1() ) return false;
	if ( cst1.atom2() != cst2.atom2() ) return false;
	core::scoring::func::Func const & func1( cst1.get_func() );
	core::scoring::func::Func const & func2( cst2.get_func() );
	core::scoring::func::SplineFunc const * splinefunc1( dynamic_cast< core::scoring::func::SplineFunc const * >( &func1 ) );
	core::scoring::func::SplineFunc const * splinefunc2( dynamic_cast< core::scoring::func::SplineFunc const * >( &func2 ) );
	if ( splinefunc1 == nullptr || splinefunc2 == nullptr ) return false;
	return splinefunc1->is_approximately_equal( *splinefunc2, FLOAT_COMPARISON_THRESHOLD );
}

/// @brief Are two angle constraints approximately equal?
bool
ang_csts_approx_equal(
	core::scoring::constraints::AngleConstraint const & cst1,
	core::scoring::constraints::AngleConstraint const & cst2
) {
	if ( cst1.atom1() != cst2.atom1() ) return false;
	if ( cst1.atom2() != cst2.atom2() ) return false;
	if ( cst1.atom3() != cst2.atom3() ) return false;
	core::scoring::func::Func const & func1( cst1.get_func() );
	core::scoring::func::Func const & func2( cst2.get_func() );
	core::scoring::func::SplineFunc const * splinefunc1( dynamic_cast< core::scoring::func::SplineFunc const * >( &func1 ) );
	core::scoring::func::SplineFunc const * splinefunc2( dynamic_cast< core::scoring::func::SplineFunc const * >( &func2 ) );
	if ( splinefunc1 == nullptr || splinefunc2 == nullptr ) return false;
	return splinefunc1->is_approximately_equal( *splinefunc2, FLOAT_COMPARISON_THRESHOLD );
}

/// @brief Are two dihedral constraints approximately equal?
bool
dihed_csts_approx_equal(
	core::scoring::constraints::DihedralConstraint const & cst1,
	core::scoring::constraints::DihedralConstraint const & cst2
) {
	if ( cst1.atom(1) != cst2.atom(1) ) return false;
	if ( cst1.atom(2) != cst2.atom(2) ) return false;
	if ( cst1.atom(3) != cst2.atom(3) ) return false;
	if ( cst1.atom(4) != cst2.atom(4) ) return false;
	core::scoring::func::Func const & func1( cst1.get_func() );
	core::scoring::func::Func const & func2( cst2.get_func() );
	core::scoring::func::SplineFunc const * splinefunc1( dynamic_cast< core::scoring::func::SplineFunc const * >( &func1 ) );
	core::scoring::func::SplineFunc const * splinefunc2( dynamic_cast< core::scoring::func::SplineFunc const * >( &func2 ) );
	if ( splinefunc1 == nullptr || splinefunc2 == nullptr ) return false;
	return splinefunc1->is_approximately_equal( *splinefunc2, FLOAT_COMPARISON_THRESHOLD );
}

/// @brief Carry out the test.
/// @returns True for success, false for failure.
bool
do_test(
	std::string const & cst_file,
	std::string const & fasta_file,
	std::string const & msa_file
) {
	bool success(true);
	std::string const errmsg( "Error in test_trRosettaConstraintGenerator app do_test() function: " );

	//Build the pose:
	core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
	{
		core::chemical::ResidueTypeSetCOP restypeset(
			core::chemical::ChemicalManager::get_instance()->residue_type_set(
			core::chemical::CENTROID_t
			)
		);
		utility::vector1< core::sequence::SequenceOP > const seqs( core::sequence::read_fasta_file( fasta_file ) );
		runtime_assert_string_msg( seqs.size() == 1, errmsg + "Exactly one FASTA sequence must be contained in " + fasta_file + "." );
		core::pose::make_pose_from_sequence( *pose, seqs[1]->ungapped_sequence(), restypeset, true );
	}

	//Mutate the glycines to alanine:
	for ( core::Size ir(1), irmax(pose->total_residue()); ir<=irmax; ++ir ) {
		if ( pose->residue_type(ir).aa() == core::chemical::aa_gly ) {
			protocols::simple_moves::MutateResidue mutres(ir, "ALA");
			mutres.apply(*pose);
		}
	}
	pose->update_residue_neighbors();

	//Generate constraints:
	utility::vector1< core::scoring::constraints::ConstraintCOP > csts1;
	{
		protocols::trRosetta_protocols::constraint_generators::trRosettaConstraintGenerator cst_gen;
		cst_gen.set_msa_file( msa_file );
		csts1 = cst_gen.apply(*pose);
	}

	//Read in Python trRosetta constraints and apply these to the Pose:
	{
		protocols::constraint_movers::ConstraintSetMover cstmover;
		cstmover.constraint_file( cst_file );
		cstmover.apply(*pose);
	}

	//Compare Python and C++ trRosetta constraints:
	{
		using namespace core::scoring::constraints;

		utility::vector1< ConstraintCOP > csts2( pose->constraint_set()->get_all_constraints() );
		if ( csts2.empty() ) {
			TR.Error << errmsg << "No constraints were read in from the constraint file from the Python version of trRosetta!" << std::endl;
			success = false;
		}
		TR << "Found " << csts2.size() << " constraints in the constraint file from the Python version of trRosetta." << std::endl;
		std::map< std::pair< core::Size, core::Size >, AtomPairConstraintCOP > dist_csts;
		std::map< std::pair< core::Size, core::Size >, AngleConstraintCOP > phi_csts;
		std::map< std::pair< core::Size, core::Size >, DihedralConstraintCOP > theta_csts;
		std::map< std::pair< core::Size, core::Size >, DihedralConstraintCOP > omega_csts;
		for ( auto const & cst : csts2 ) {
			//Determine the type of the constraint:
			{ //Atom pair constraints:
				AtomPairConstraintCOP distcst( utility::pointer::dynamic_pointer_cast< AtomPairConstraint const >( cst ) );
				if ( distcst != nullptr ) {
					std::pair< core::Size, core::Size > const respair( std::make_pair( distcst->atom1().rsd(), distcst->atom2().rsd() ) );
					runtime_assert( dist_csts.count( respair ) == 0 );
					dist_csts[respair] = distcst;
					continue;
				}
			}
			{ //Angle constraints:
				AngleConstraintCOP angcst( utility::pointer::dynamic_pointer_cast< AngleConstraint const >( cst ) );
				if ( angcst != nullptr ) {
					std::pair< core::Size, core::Size > const respair( std::make_pair( angcst->atom1().rsd(), angcst->atom3().rsd() ) );
					runtime_assert( phi_csts.count( respair ) == 0 );
					phi_csts[respair] = angcst;
					continue;
				}
			}
			{ //Angle constraints:
				DihedralConstraintCOP dihedcst( utility::pointer::dynamic_pointer_cast< DihedralConstraint const >( cst ) );
				if ( dihedcst != nullptr ) {
					std::pair< core::Size, core::Size > const respair( std::make_pair( dihedcst->atom(1).rsd(), dihedcst->atom(4).rsd() ) );
					if (
							utility::strip( pose->residue_type(respair.first).atom_name( dihedcst->atom(1).atomno() ) ) == "CA" &&
							utility::strip( pose->residue_type(respair.second).atom_name( dihedcst->atom(4).atomno() ) ) == "CA"
							) { //This is an omega constraint
						runtime_assert( omega_csts.count( respair ) == 0 );
						omega_csts[respair] = dihedcst;
						continue;
					} else if (
							utility::strip( pose->residue_type(respair.first).atom_name( dihedcst->atom(1).atomno() ) ) == "N" &&
							utility::strip( pose->residue_type(respair.second).atom_name( dihedcst->atom(4).atomno() ) ) == "CB"
							) { //This is a theta constraint
						runtime_assert( theta_csts.count( respair ) == 0 );
						theta_csts[respair] = dihedcst;
						continue;
					}
				}
			}
			TR.Error << errmsg << "A constraint was unrecognized!" << std::endl;
			success = false;
		}
		if ( dist_csts.empty() ) {
			TR.Error << errmsg << "No distance constraints were found in the constraint file!" << std::endl;
			success = false;
		} else {
			TR << "Found " << dist_csts.size() << " distance constraints in the constraint file." << std::endl;
		}
		if ( phi_csts.empty() ) {
			TR.Error << errmsg << "No phi constraints were found in the constraint file!" << std::endl;
			success = false;
		} else {
			TR << "Found " << phi_csts.size() << " phi constraints in the constraint file." << std::endl;
		}
		if ( theta_csts.empty() ) {
			TR.Error << errmsg << "No theta constraints were found in the constraint file!" << std::endl;
			success = false;
		} else {
			TR << "Found " << theta_csts.size() << " theta constraints in the constraint file." << std::endl;
		}
		if ( omega_csts.empty() ) {
			TR.Error << errmsg << "No omega constraints were found in the constraint file!" << std::endl;
			success = false;
		} else {
			TR << "Found " << omega_csts.size() << " omega constraints in the constraint file." << std::endl;
		}

		//Iterate through all trRosetta C++ constraints and compare.
		core::Size distcounter(0), phicounter(0), thetacounter(0), omegacounter(0);
		for ( auto const & cst : csts1 ) {
			{ //Is this a distance constraint?
				AtomPairConstraintCOP distcst( utility::pointer::dynamic_pointer_cast< AtomPairConstraint const >(cst) );
				if ( distcst != nullptr ) {
					std::pair< core::Size, core::Size > const respair( distcst->atom1().rsd(), distcst->atom2().rsd() );
					if ( dist_csts.count(respair) != 0 ) {
						if ( !dist_csts_approx_equal( *distcst, *(dist_csts.at(respair) ) ) ) {
							TR.Error << errmsg << "Comparison failed between distance constraints!" << std::endl;
							success = false;
						}
						TR.Debug << summarize_csts( distcst, dist_csts.at(respair) ) << std::endl;
						++distcounter;
					}
					continue;
				}
			}
			{ //Is this an angle constraint?
				AngleConstraintCOP angcst( utility::pointer::dynamic_pointer_cast< AngleConstraint const >(cst) );
				if ( angcst != nullptr ) {
					std::pair< core::Size, core::Size > const respair( angcst->atom1().rsd(), angcst->atom3().rsd() );
					if ( phi_csts.count(respair) != 0 ) {
						if ( !ang_csts_approx_equal(*angcst, *(phi_csts.at(respair)) ) ) {
							TR.Error << errmsg << "Comparison failed between angle constraints!" << std::endl;
							success = false;
						}
						TR.Debug << summarize_csts( angcst, phi_csts.at(respair) ) << std::endl;
						++phicounter;
					}
					continue;
				}
			}
			{ //Is this a dihedral constraint?
				DihedralConstraintCOP dihedcst( utility::pointer::dynamic_pointer_cast< DihedralConstraint const >(cst) );
				if ( dihedcst != nullptr ) {
					std::pair< core::Size, core::Size > const respair( dihedcst->atom(1).rsd(), dihedcst->atom(4).rsd() );
					if ( utility::strip( pose->residue_type(respair.first).atom_name( dihedcst->atom(1).atomno() ) ) == "CA" ) {
						if ( omega_csts.count(respair) != 0 ) {
							if ( !dihed_csts_approx_equal( *dihedcst, *( omega_csts.at(respair) ) ) ) {
								TR.Error << errmsg << "Comparison failed between omega dihedral constraints!" << std::endl;
								success = false;
							}
							TR.Debug << summarize_csts( dihedcst, omega_csts.at(respair) ) << std::endl;
							++omegacounter;
						}
					} else {
						if ( theta_csts.count(respair) != 0 ) {
							if ( !dihed_csts_approx_equal( *dihedcst, *(theta_csts.at(respair)) ) ) {
								TR.Error << errmsg << "Comparison failed between theta dihedral constraints!" << std::endl;
								success = false;
							}
							TR.Debug << summarize_csts( dihedcst, theta_csts.at(respair) ) << std::endl;
							++thetacounter;
						}
					}
					continue;
				}
			}
		}

		if ( distcounter != dist_csts.size() ) {
			TR.Error << errmsg << "Only compared " << distcounter << " of " << dist_csts.size() << " atom pair constraints!" << std::endl;
			success = false;
		} else {
			TR << "Compared " << distcounter << " of " << dist_csts.size() << " atom pair constraints." << std::endl;
		}
		if ( phicounter != phi_csts.size() ) {
			TR.Error << errmsg << "Only compared " << phicounter << " of " << phi_csts.size() << " phi angle constraints!" << std::endl;
			success = false;
		} else {
			TR << "Compared " << phicounter << " of " << phi_csts.size() << " phi angle constraints." << std::endl;
		}
		if ( thetacounter != theta_csts.size() ) {
			TR.Error << errmsg << "Only compared " << thetacounter << " of " << theta_csts.size() << " theta dihedral constraints!" << std::endl;
			success = false;
		} else {
			TR << "Compared " << thetacounter << " of " << theta_csts.size() << " theta dihedral constraints constraints." << std::endl;
		}
		if ( omegacounter != omega_csts.size() ) {
			TR.Error << errmsg << "Only compared " << omegacounter << " of " << omega_csts.size() << " omega dihedral constraints!" << std::endl;
			success = false;
		} else {
			TR << "Compared " << omegacounter << " of " << omega_csts.size() << " omega dihedral constraints." << std::endl;
		}
	}

	return success;
}

#endif //USE_TENSORFLOW

/// @brief Program entry point.
int
main(
	int argc,
	char * argv []
) {
	try {
#ifdef USE_TENSORFLOW
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		NEW_OPT( csts_file, "The input constraints file.  This should be the temporary file generated by the Python version of trRosetta.  All of the temporary files that it needs to define splines must also be present.", "" );
		NEW_OPT( msa_file, "The input multiple sequence alignment file.", "" );

		register_options();
		devel::init( argc, argv );

		TR << "Starting test_trRosettaConstraintGenerator." << std::endl;
		TR << "This is a unit test (albeit one in the integration test suite) for the "
			"trRosettaConstraintGenerator.  This compares the constraints generated "
			"by the C++ ConstraintGenerator to known outputs from the original "
			"Python version of trRosetta." << std::endl;
		TR << "Test created 24 February 2021 by Vikram K. Mulligan, Flatiron "
			"Institute (vmulligan@flatironinstitute.org)." << std::endl;

		runtime_assert_string_msg(
			option[ csts_file ].user() && !option[ csts_file ]().empty(),
			"A constraints file must be provided with the -csts_file flag."
		);
		runtime_assert_string_msg(
			option[ msa_file ].user() && !option[ msa_file ]().empty(),
			"A multiple sequence alignment must be provided with the -msa_file flag."
		);
		runtime_assert_string_msg(
			option[ in::file::fasta ].user() && (option[ in::file::fasta ]().size() == 1),
				"A FASTA file must be provided with the -in:file:fasta flag."
		);

		bool const success( do_test( option[csts_file](), option[in::file::fasta]()[1], option[msa_file]() ) );
		TR << "Test " << (success ? "PASSED" : "FAILED" ) << "!" << std::endl;

		basic::citation_manager::CitationManager::get_instance()->write_all_citations_and_unpublished_author_info();

		if( !success ) {
			utility_exit(); //Exit with failure status if we don't succeed.
		}

#else
		devel::init( argc, argv );
		utility_exit_with_message(
			"The test_trRosettaConstraintGenerator application must be compiled with -extras=tensorflow or -extras=tensorflow_gpu!\n\n"
			+ basic::tensorflow_manager::get_tensorflow_compilation_instructions("test_trRosettaConstraintGenerator application")
		);
#endif

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/vmullig/simple_cycpep_predict.cc
/// @brief Predicts structure of N-to-C cyclic peptides from amino acid sequence.
/// @details
/// @author Vikram K. Mulligan (vmullig@uw.edu)

//General includes
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <numeric/conversions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <numeric/constants.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/protein_interface_design/filters/HbondsToResidueFilter.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/Energies.hh>

//Constraints
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

//Application-specific includes
#include <numeric/crick_equations/HelixParams.hh>
#include <protocols/helical_bundle/FitSimpleHelix.hh>
#include <protocols/helical_bundle/util.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <stdio.h>

#define PEPBOND_LENGTH 1.328685
#define PEPBOND_C_ANGLE 2.02807246864
#define PEPBOND_N_ANGLE 2.12406564732

//Tracer:
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.vmullig.simple_cycpep_predict" );

//Options (ugh -- global variables):
OPT_KEY (String, sequence_file)
OPT_KEY (Integer, genkic_closure_attempts)
OPT_KEY (Integer, genkic_min_solution_count)
OPT_KEY (Boolean, cyclic_permutations)
OPT_KEY (Boolean, use_rama_filter)
OPT_KEY (Real, rama_cutoff)
OPT_KEY (Real, high_hbond_weight_multiplier)
OPT_KEY (Real, min_genkic_hbonds)
OPT_KEY (Real, min_final_hbonds)
OPT_KEY (Real, hbond_energy_cutoff)
OPT_KEY (Integer, fast_relax_rounds)
OPT_KEY (Boolean, count_sc_hbonds)

/// @brief Set up the options for this pilot app.
void register_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( sequence_file, "Filename of a file specfying the sequence, as a series of whitespace-separated full residue names (e.g. \"ALA LYS DARG DPRO HYP\").  Required input.", "");
	NEW_OPT( genkic_closure_attempts, "How many closure attempts should we make for each job?  Default 10,000.", 10000);
	NEW_OPT( genkic_min_solution_count, "How many solutions should genKIC find before picking one?  Default 1.", 1);
	NEW_OPT( cyclic_permutations, "Should cyclic permutations of the sequence be considered when setting up the kinematic closure?  Default true.", true);
	NEW_OPT( use_rama_filter, "Should GenKIC solutions be filtered based on rama score?  True by default.", true);
	NEW_OPT( rama_cutoff, "The maximum rama score value that's permitted in the accepted GenKIC solutions.  Default 0.3.", 0.3 );
	NEW_OPT( high_hbond_weight_multiplier, "In parts of the protocol involving upweighting of the backbone hbond terms, by what factor should backbone hbond energy be upweighted?  Default 10.0.", 10.0 );
	NEW_OPT( min_genkic_hbonds, "The minimum number of backbone hbonds for a solution to pass during GenKIC closure.  Default 3.", 3.0);
	NEW_OPT( min_final_hbonds, "The minimum number of backbone hbonds for a solution to pass after final relaxtion.  Default 0 (report only).", 0.0);
	NEW_OPT( hbond_energy_cutoff, "The mainchain hbond energy threshold for something to be counted as a hydrogen bond.  Default -0.25.", -0.25);
	NEW_OPT( fast_relax_rounds, "The number of rounds of FastRelax to perform at each FastRelax step.  Note that there are two such steps: a high-hbond initial FastRelax applied to all GenKIC solutions, and a regular scorefunction final FastRelax applied to the best GenKIC solution.  Default 3.", 3);
	NEW_OPT( count_sc_hbonds, "Should sidechain-backbone and sidechain-sidechain hydrogen bonds be counted in the total hydrogen bond count?  Default false.", false);
	return;
}

/// @brief Actually build the geometry that we'll be working with.
///
void build_polymer(
	core::pose::PoseOP pose,
	utility::vector1<std::string> const &restypes
) {
	using namespace protocols::cyclic_peptide;
	core::Size const nres( restypes.size() );
	runtime_assert(restypes.size() >=4 );

	TR << "Building sequence ";
	for ( core::Size i=1; i<=nres; ++i ) {
		TR << restypes[i];
		if ( i<nres ) TR << " ";
	}
	TR << "." << std::endl;

	PeptideStubMover stubmover;

	stubmover.set_reset_mode(true);
	stubmover.reset_mover_data();
	for ( core::Size i=1; i<=nres; ++i ) {
		stubmover.add_residue( "Append", restypes[i], 0, false, "", 1, 0, "" );
	}

	stubmover.apply(*pose);

	TR << "Build successful." << std::endl;

	return;
}

/// @brief Read a sequence (as a series of full names, separated by whitespace) and store
/// it in a string vector.
void
read_sequence (
	std::string const &seqfile,
	utility::vector1 < std::string > &resnames
) {
	using namespace utility::io;
	resnames.clear();

	izstream infile;
	infile.open( seqfile );
	runtime_assert_string_msg( infile.good(), "Error in read_sequence() in app simple_cycpep_predict:  Unable to open sequence file for read!" );

	std::string curline(""); //Buffer for current line.
	utility::vector1< std::string > lines; //Storing all lines

	//Read the file:
	while ( getline(infile, curline) ) {
		if ( curline.size() < 1 ) continue; //Ignore blank lines.
		lines.push_back( curline );
	}
	infile.close();

	//Parse the lines:
	for ( core::Size i=1, imax=lines.size(); i<=imax; ++i ) { //Loop through all lines
		std::istringstream curline(lines[i]);
		std::string oneword("");
		while ( !curline.eof() ) {
			curline >> oneword;
			resnames.push_back( oneword );
		}
	}

	if ( TR.visible() ) {
		TR << "Parsed the following sequence:" << std::endl;
		for ( core::Size i=1, imax=resnames.size(); i<=imax; ++i ) {
			TR << resnames[i];
			if ( i<imax ) TR << ", ";
		}
		TR << "." << std::endl;
	}

	runtime_assert_string_msg( resnames.size() >= 4, "Error in simple_cycpcp_predict app read_sequence() function!  The minimum number of residues for a cyclic peptide is 4.  (GenKIC requires three residues, plus a fourth to serve as an anchor)." );

	return;
}

/// @brief Set up the DeclareBond mover used to connect the termini.
///
void
set_up_termini_mover (
	protocols::cyclic_peptide::DeclareBondOP termini,
	core::pose::PoseCOP pose,
	bool const native=false
) {
	core::Size const nres(pose->n_residue());

	runtime_assert_string_msg(pose->residue(1).has_lower_connect(), "Error in simple_cycpep_predict app set_up_termini_mover() function: residue 1 does not have a LOWER_CONNECT.");
	runtime_assert_string_msg(pose->residue(nres).has_upper_connect(), "Error in simple_cycpep_predict app set_up_termini_mover() function: the final residue does not have an UPPER_CONNECT.");
	std::string firstatom( pose->residue(1).atom_name( pose->residue(1).lower_connect_atom() ) );
	std::string lastatom( pose->residue(nres).atom_name( pose->residue(nres).upper_connect_atom() ) );

	if ( native ) {
		TR << "Setting up terminal bond for the native pose between residue 1, atom " << firstatom << " and residue " << nres << ", atom " << lastatom << "." << std::endl;
	} else {
		TR << "Setting up terminal bond between residue 1, atom " << firstatom << " and residue " << nres << ", atom " << lastatom << "." << std::endl;
	}

	termini->set( nres, lastatom, 1, firstatom, false, false, 0, 0, false  );

	return;
}

/// @brief Takes a vector of residue names, chooses a random number for cyclic offset, and
/// does a cyclic permutation.
/// @details Returns the offset and stores the new string vector in resnames_copy.
core::Size
do_cyclic_permutation ( utility::vector1 <std::string> const &resnames, utility::vector1 <std::string> &resnames_copy ) {
	core::Size const nname( resnames.size() );//Number of residue names
	core::Size const offset( static_cast<core::Size>(numeric::random::rg().random_range(0,resnames.size()-1)) );

	resnames_copy.clear();
	resnames_copy.resize(nname, "");
	core::Size counter(0);
	for ( core::Size i=offset+1; i<=nname; ++i ) {
		++counter;
		resnames_copy[counter] = resnames[i];
	}
	for ( core::Size i=1; i<=offset; ++i ) {
		++counter;
		resnames_copy[counter] = resnames[i];
	}

	TR << "Circularly shifted residue list by " << offset << ".  New list is: ";
	for ( core::Size i=1; i<=nname; ++i ) {
		TR << resnames_copy[i];
		if ( i<nname ) TR << ", ";
	}
	TR << std::endl;
	TR.flush();

	return offset;
}

/// @brief Imports the native pose and sets up a terminial peptide bond.
///
void
import_and_set_up_native (
	std::string const &native_file,
	core::pose::PoseOP native_pose,
	core::Size const expected_residue_count
) {
	core::import_pose::pose_from_pdb(*native_pose, native_file);
	TR << "Improrting native structure from " << native_file << "." << std::endl;
	runtime_assert_string_msg( native_pose->n_residue() == expected_residue_count, "Error in simple_cycpep_predict app!  The imported native pose has a different number of residues than the sequence provided." );

	TR << "Stripping termini from native structure." << std::endl;
	core::pose::remove_lower_terminus_type_from_pose_residue(*native_pose, 1);
	core::pose::remove_upper_terminus_type_from_pose_residue(*native_pose, expected_residue_count);

	//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
	protocols::cyclic_peptide::DeclareBondOP termini( new protocols::cyclic_peptide::DeclareBond );
	set_up_termini_mover( termini, native_pose, true );
	termini->apply(*native_pose);

	return;
}

//Function to add cyclic constraints to a pose:
void add_cyclic_constraints (core::pose::PoseOP pose) {
	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::id;

	TR << "Setting up cyclic constraints." << std::endl;

	core::Size const nres(pose->n_residue());

	//The four atoms defining the peptide bond:
	AtomID const atom_a( pose->residue(nres).type().icoor(pose->residue(nres).upper_connect_atom()).stub_atom1().atomno(), nres );
	AtomID const atom_b( pose->residue(nres).upper_connect_atom(), nres );
	AtomID const atom_c( pose->residue(1).lower_connect_atom(), 1 );
	core::Size atom_d_index(0);
	for ( core::Size i=1, imax=pose->residue(1).n_mainchain_atoms(); i<=imax; ++i ) { //Find the atom index of the first mainchain atom with the lower_connect atom as a parent.
		if ( i == atom_c.atomno() ) continue;
		if ( pose->residue(1).type().icoor(i).stub_atom1().atomno() == atom_c.atomno() ) {
			atom_d_index=i;
			break;
		}
	}
	AtomID const atom_d( atom_d_index, 1 );

	TR << "The following four atoms define the terminal bond:" << std::endl;
	TR << "1.\tRes=" << atom_a.rsd() << "\tAtom=" << pose->residue(atom_a.rsd()).atom_name(atom_a.atomno()) << std::endl;
	TR << "2.\tRes=" << atom_b.rsd() << "\tAtom=" << pose->residue(atom_b.rsd()).atom_name(atom_b.atomno()) << std::endl;
	TR << "3.\tRes=" << atom_c.rsd() << "\tAtom=" << pose->residue(atom_c.rsd()).atom_name(atom_c.atomno()) << std::endl;
	TR << "4.\tRes=" << atom_d.rsd() << "\tAtom=" << pose->residue(atom_d.rsd()).atom_name(atom_d.atomno()) << std::endl;

	{//Peptide bond length constraint:
		FuncOP harmfunc1( new HarmonicFunc( PEPBOND_LENGTH, 0.01) );
		ConstraintCOP distconst1( new AtomPairConstraint ( atom_b, atom_c, harmfunc1 ) );
		pose->add_constraint (distconst1);
	}

	{ //Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		FuncOP circharmfunc1( new CircularHarmonicFunc( numeric::constants::d::pi, 0.02) );
		ConstraintCOP dihedconst1( new DihedralConstraint ( atom_a, atom_b, atom_c, atom_d, circharmfunc1) );
		pose->add_constraint (dihedconst1);
	}

	{ //Peptide bond angle constraints:
		FuncOP circharmfunc2a( new CircularHarmonicFunc( PEPBOND_C_ANGLE, 0.02) );
		FuncOP circharmfunc2b( new CircularHarmonicFunc( PEPBOND_N_ANGLE, 0.02) );
		ConstraintCOP angleconst1( new AngleConstraint ( atom_a, atom_b, atom_c, circharmfunc2a) );
		ConstraintCOP angleconst2( new AngleConstraint ( atom_b, atom_c, atom_d, circharmfunc2b) );
		pose->add_constraint (angleconst1);
		pose->add_constraint (angleconst2);
	}

	TR << "Finished setting up constraints." << std::endl;

	return;
}

/// @brief Sets all omega values to 180, and randomizes mainchain torsions.
/// @details For alpha-amino acids, mainchain torsions are randomized by the Ramachandran plot.
/// For other residue types, just randomizes mainchain torsions other than peptide bonds.
void
set_mainchain_torsions (core::pose::PoseOP pose) {
	TR << "Randomizing mainchain torsions." << std::endl;
	core::Size const nres(pose->n_residue());
	for ( core::Size i=1; i<=nres; ++i ) { //Loop through all residues
		if ( pose->residue(i).type().is_alpha_aa() ) {
			core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran(); //Get the Rama scoring function
			core::Real phi(0.0), psi(0.0);
			//TR << "aa" << i << "=" << pose->residue_type(i).aa() << std::endl; //DELETE ME
			rama.random_phipsi_from_rama( pose->residue_type(i).aa(), phi, psi); //TODO -- use backbone_aa
			pose->set_phi(i,phi);
			pose->set_psi(i,psi);
			if ( i!=nres ) pose->set_omega(i, 180.0);
		} else { //If this is not an alpha-amino acid:
			for ( core::Size j=1, jmax=pose->residue(i).mainchain_torsions().size(); j<=jmax; ++j ) { //Loop through all mainchain torsions.
				if ( i==nres && j==jmax ) continue; //Skip the last mainchain torsion (not a DOF).
				core::Real setting(180.0);
				if ( j!=jmax ) {
					setting = numeric::random::rg().uniform()*360.0 - 180.0;
				}
				pose->set_torsion( core::id::TorsionID(i, core::id::BB, j), setting );
			}
		}
	}
	return;
}

void
set_up_hbond_filter(
	protocols::filters::CombinedFilterOP total_hbond,
	core::Size const nres,
	core::scoring::ScoreFunctionOP sfxn,
	core::Real const &min_hbonds
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	total_hbond->set_threshold( -1.0 * min_hbonds );
	for ( core::Size i=1; i<=nres; ++i ) { //Loop through all residues and add hbond counters
		protocols::protein_interface_design::filters::HbondsToResidueFilterOP hbondfilt( new protocols::protein_interface_design::filters::HbondsToResidueFilter );
		hbondfilt->set_resnum(i);
		hbondfilt->set_sidechain( option[count_sc_hbonds]() );
		hbondfilt->set_energy_cutoff( static_cast<core::Real>(option[hbond_energy_cutoff]()) );
		hbondfilt->set_partners(0);
		hbondfilt->set_scorefxn( sfxn );
		total_hbond->add_filter( hbondfilt, -0.5, false );
	}
	return;
}

/// @brief Use GeneralizedKIC to close the pose.
///
bool
genkic_close(
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP sfxn_highhbond,
	protocols::filters::CombinedFilterOP total_hbond
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::generalized_kinematic_closure;
	using namespace protocols::generalized_kinematic_closure::selector;

	TR << "Performing GeneralizedKIC closure of loop." << std::endl;

	//Number of residues in the pose:
	core::Size const nres( pose->n_residue() );
	runtime_assert( nres >= 4 ); //Already checked at sequence load time, so should be true, but let's make sure.

	//Randomly pick one of the middle residues to be the anchor residue:
	core::Size const anchor_res( numeric::random::rg().random_range(2, nres-1) );
	core::Size const first_loop_res( anchor_res + 1 );
	core::Size const last_loop_res( anchor_res - 1 );

	//Randomly pick a residue to be the middle pivot residue.  Can't be first in loop, last in loop, or anchor res.
	core::Size middle_loop_res( numeric::random::rg().random_range(1, nres-3 ) );
	if ( middle_loop_res == last_loop_res ) { middle_loop_res += 3; }
	else if ( middle_loop_res == anchor_res ) { middle_loop_res +=2; }
	else if ( middle_loop_res == first_loop_res ) { middle_loop_res +=1; }
	if ( middle_loop_res > nres ) { middle_loop_res -= nres; }

	//Create the pre-selection mover and set options.  TODO.
	protocols::rosetta_scripts::ParsedProtocolOP pp( new protocols::rosetta_scripts::ParsedProtocol );
	//protocols::filters::CombinedFilterOP total_hbond( new protocols::filters::CombinedFilter );
	//set_up_hbond_filter( total_hbond, nres, sfxn_default, static_cast<core::Real>( option[min_genkic_hbonds]() ) );
	pp->add_mover_filter_pair( NULL, "Total_Hbonds", total_hbond );
	protocols::relax::FastRelaxOP frlx( new protocols::relax::FastRelax(sfxn_highhbond, 1) );
	for ( core::Size i=1, imax=static_cast<core::Size>(option[fast_relax_rounds]()); i<=imax; ++i ) {
		pp->add_mover_filter_pair( frlx, "High_Hbond_FastRelax", NULL );
	}
	//CONTINUE HERE

	//Create the mover and set options:
	GeneralizedKICOP genkic( new GeneralizedKIC );
	genkic->set_selector_type( lowest_energy_selector );
	genkic->set_closure_attempts( static_cast<core::Size>( option[genkic_closure_attempts].value() ) );
	genkic->set_min_solution_count( static_cast<core::Size>( option[genkic_min_solution_count].value() ) );
	genkic->set_selector_scorefunction( sfxn_highhbond );
	genkic->set_preselection_mover(pp);

	//Define the loop residues:
	for ( core::Size i=first_loop_res; i<=nres; ++i ) { genkic->add_loop_residue(i); }
	for ( core::Size i=1; i<=last_loop_res; ++i ) { genkic->add_loop_residue(i); }

	//Set pivots:
	std::string at1(""), at2(""), at3("");
	if ( pose->residue(first_loop_res).type().is_alpha_aa() ) { at1="CA"; }
	else if ( pose->residue(first_loop_res).type().is_beta_aa() ) { at1="CM"; }
	else if ( pose->residue(first_loop_res).type().is_gamma_aa() ) { at1="C3"; }
	else { utility_exit_with_message( "Unrecognized residue type at loop start.  Currently, this app only works with alpha, beta, and gamma amino acids." ); }
	if ( pose->residue(middle_loop_res).type().is_alpha_aa() ) { at2="CA"; }
	else if ( pose->residue(middle_loop_res).type().is_beta_aa() ) { at2="CM"; }
	else if ( pose->residue(middle_loop_res).type().is_gamma_aa() ) { at2="C3"; }
	else { utility_exit_with_message( "Unrecognized residue type at loop midpoint.  Currently, this app only works with alpha, beta, and gamma amino acids." ); }
	if ( pose->residue(last_loop_res).type().is_alpha_aa() ) { at3="CA"; }
	else if ( pose->residue(last_loop_res).type().is_beta_aa() ) { at3="CM"; }
	else if ( pose->residue(last_loop_res).type().is_gamma_aa() ) { at3="C3"; }
	else { utility_exit_with_message( "Unrecognized residue type at loop midpoint.  Currently, this app only works with alpha, beta, and gamma amino acids." ); }
	genkic->set_pivot_atoms( first_loop_res, at1, middle_loop_res, at2, last_loop_res, at3 );

	//Close the bond:
	std::string const firstatom( pose->residue(1).atom_name( pose->residue(1).lower_connect_atom() ) );
	std::string const lastatom( pose->residue(nres).atom_name( pose->residue(nres).upper_connect_atom() ) );
	genkic->close_bond( nres, lastatom, 1, firstatom, 0, "", 0, "", PEPBOND_LENGTH, PEPBOND_C_ANGLE/numeric::constants::d::pi*180.0, PEPBOND_N_ANGLE/numeric::constants::d::pi*180.0, 180.0, false, false );

	//Add perturbers:
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( i==anchor_res ) continue; //Can't perturb the anchor residue.
		if ( pose->residue(i).type().is_alpha_aa() ) {
			genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_alpha_backbone_by_rama );
			genkic->add_residue_to_perturber_residue_list(i);
		} else {
			//TODO Randomize mainchain torsions here for beta- and gamma-amino acids.
			utility_exit_with_message( "Handling of beta- and gamma-amino acids in setup of the genKIC perturber in the simple_cycpep_predict app has not yet been written.  TODO." );
		}
	}

	//Add rama check filters:
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( i!=first_loop_res && i!=middle_loop_res && i!=last_loop_res ) continue; //Just filter the pivots.
		if ( pose->residue(i).type().is_alpha_aa() ) {
			genkic->add_filter( protocols::generalized_kinematic_closure::filter::alpha_aa_rama_check );
			genkic->set_filter_resnum(i);
			genkic->set_filter_rama_cutoff_energy( static_cast<core::Real>(option[rama_cutoff]()) );
		}
	}

	//Add bump check filter:
	genkic->add_filter( protocols::generalized_kinematic_closure::filter::loop_bump_check );

	//Apply the mover:
	genkic->apply( *pose );

	return genkic->last_run_successful();
}

/// @brief Given a pose that has undergone an N-residue cyclic permutation, restore
/// the original pose, without the permutation.
void
depermute (
	core::pose::PoseOP pose,
	core::Size const offset
) {

	/// 1 2 3 4 5 6 7 8
	/// 2 3 4 5 6 7 8 1
	/// 3 4 5 6 7 8 1 2
	/// 4 5 6 7 8 1 2 3

	if ( offset==0 ) return; //Do nothing if the pose was not offset.

	core::Size const nres(pose->n_residue());
	debug_assert(nres > offset);
	core::Size const old_first_res_index( nres-offset+1 );

	//TR << "nres=" << nres << " offset=" << offset << " old_first_res_index=" << old_first_res_index << std::endl; //DELETE ME

	core::pose::PoseOP newpose( new core::pose::Pose );

	for ( core::Size ir=old_first_res_index; ir<=nres; ++ir ) {
		if ( ir == old_first_res_index ) {
			newpose->append_residue_by_jump( *(pose->residue(ir).clone()), 0, "", "", true );
		} else {
			newpose->append_residue_by_bond( *(pose->residue(ir).clone()), false, 0, 0, 0, false, false );
		}
	}

	for ( core::Size ir=1; ir<old_first_res_index; ++ir ) {
		newpose->append_residue_by_bond( *(pose->residue(ir).clone()), false, 0, 0, 0, false, false );
	}

	//I don't bother to set up cyclic constraints, since we won't be doing any more minimization after calling this function.

	//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
	protocols::cyclic_peptide::DeclareBondOP termini( new protocols::cyclic_peptide::DeclareBond );
	set_up_termini_mover( termini, newpose );
	termini->apply(*newpose);

	(*pose) = (*newpose);

	return;
}

/// @brief Align pose to native_pose, and return the RMSD between the two poses.
/// @details Assumes that the pose has already been de-permuted (i.e. the native and the pose line up).
core::Real
align_and_calculate_rmsd(
	core::pose::PoseOP pose,
	core::pose::PoseCOP native_pose
) {
	core::Size const nres( pose->n_residue() );
	debug_assert( native_pose->n_residue() == nres ); //Should be true.

	core::id::AtomID_Map< core::id::AtomID > amap;
	core::pose::initialize_atomid_map(amap, *pose, core::id::BOGUS_ATOM_ID);
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		for ( core::Size ia=1, iamax=native_pose->residue(ir).type().first_sidechain_atom(); ia<iamax; ++ia ) { //Loop through all mainchain heavyatoms (including atoms coming off mainchain that are not sidechain atoms, like peptide "O").
			if ( native_pose->residue(ir).type().atom_is_hydrogen(ia) ) continue;
			amap[ core::id::AtomID(ia,ir) ] = core::id::AtomID(ia,ir);
			//TR << "Adding ia=" << ia << " ir=" << ir << " to map." << std::endl; //DELETE ME
		}
	}
	return core::scoring::superimpose_pose( *pose, *native_pose, amap ); //Superimpose the pose and return the RMSD.
}

int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try {
		register_options();
		devel::init(argc, argv);

		if ( TR.visible() ) {
			TR << "Starting simple_cycpep_predict.cc" << std::endl;
			TR << "Pilot app created 16 September 2015 by Vikram K. Mulligan, Ph.D., Baker laboratory." << std::endl;
			TR << "For questions, contact vmullig@uw.edu." << std::endl << std::endl;
			TR.flush();
		}

		//THE FOLLOWING THINGS ARE DONE ONCE PER PROGRAM EXECUTION:
		//Initial checks:
		runtime_assert_string_msg( option[genkic_closure_attempts]() >= 0, "Error in simple_cycpep_predict app: the number of GeneralizedKIC closure attempts (\"-genkic_closure_attempts\" flag) cannot be negative.  (Note also that setting this to zero is risky, since GenKIC will continue to seek solutions until the minimum number of solutions is reached.)" );
		runtime_assert_string_msg( option[genkic_min_solution_count]() >= 0, "Error in simple_cycpep_predict app: the minimum number of GenKIC solutions (\"-genkic_min_solution_count\" flag) cannot be negative.  (Note also that setting this to zero means no minimum.)" );
		runtime_assert_string_msg( !(option[genkic_closure_attempts]() == 0 && option[genkic_min_solution_count]() == 0), "Error in simple_cycpep_predict app: both the \"-genkic_closure_attempts\" and \"-genkic_min_solution_count\" flags were set to zero.  This would result in GenKIC looping infinitely." );
		runtime_assert_string_msg( option[min_genkic_hbonds]() >= 0.0, "Error in simple_cycpep_predict app: the minimum number of hbonds during GenKIC steps (\"-min_genkic_hbonds\" flag) can be zero, but cannot be negative." );
		runtime_assert_string_msg( option[min_final_hbonds]() >= 0.0, "Error in simple_cycpep_predict app: the minimum number of hbonds after relaxation steps (\"-min_final_hbonds\" flag) can be zero, but cannot be negative." );
		runtime_assert_string_msg( option[fast_relax_rounds]() > 0, "Error in simple_cycpep_predict app: the number of FastRelax rounds (\"-fast_relax_rounds\" flag) must be greater than zero." );
		runtime_assert_string_msg( !( option[out::file::silent].user() && option[out::file::o].user() ), "Error in simple_cycpep_predict app: either silent file output (\"-out:file:silent\" flag) or PDB output (\"-out:file:o\") output may be used, but not both." );

		//Set up output file names:
		std::string out_filename("S_");
		std::string out_scorefilename("default.sc");
		bool silent_out(false);
		if ( option[out::file::silent].user() ) {
			out_filename=option[out::file::silent]();
			silent_out=true;
			option[ basic::options::OptionKeys::out::file::silent_struct_type ].def( "binary"); //Force binary file output.
		} else if ( option[out::file::o].user() ) {
			out_filename=option[out::file::o]();
			silent_out=false;
		}
		if ( option[out::file::scorefile].user() ) { out_scorefilename=option[out::file::scorefile](); }

		//Get the scorefunction:
		core::scoring::ScoreFunctionOP sfxn_default( core::scoring::get_score_function() );
		//Create a scorefunction variant with upweighted backbone hbond terms:
		core::scoring::ScoreFunctionOP sfxn_highhbond( sfxn_default->clone() );
		core::Real const hbond_multiplier( static_cast<core::Real>(option[high_hbond_weight_multiplier]()) );
		sfxn_highhbond->set_weight( core::scoring::hbond_lr_bb, hbond_multiplier * sfxn_default->get_weight(core::scoring::hbond_lr_bb) ); //Upweight the long-range backbone hbonds
		sfxn_highhbond->set_weight( core::scoring::hbond_sr_bb, hbond_multiplier * sfxn_default->get_weight(core::scoring::hbond_sr_bb) ); //Upweight the short-range backbone hbonds
		//Create variants of the above two scorefunctions with constraint weights turned on:
		core::scoring::ScoreFunctionOP sfxn_default_cst( sfxn_default->clone() );
		core::scoring::ScoreFunctionOP sfxn_highhbond_cst( sfxn_highhbond->clone() );
		if ( sfxn_default->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::atom_pair_constraint, 1.0); }
		if ( sfxn_default->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::angle_constraint, 1.0); }
		if ( sfxn_default->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::dihedral_constraint, 1.0); }
		if ( sfxn_highhbond->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::atom_pair_constraint, 1.0); }
		if ( sfxn_highhbond->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::angle_constraint, 1.0); }
		if ( sfxn_highhbond->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::dihedral_constraint, 1.0); }


		//Get the sequence that we're considering:
		runtime_assert_string_msg( option[sequence_file].user(), "Error in simple_cycpep_predict app!  A sequence file must be provided (specified with the \"-sequence_file <filename>\" flag." );
		std::string const seqfile( option[sequence_file]() );
		utility::vector1 < std::string > resnames;
		read_sequence( seqfile, resnames );

		//Get the native sequence that we will compare to.
		core::pose::PoseOP native_pose;
		std::string native_file("");
		if ( option[in::file::native].user() ) {
			native_file=option[in::file::native]();
			native_pose=core::pose::PoseOP(new core::pose::Pose);
			import_and_set_up_native ( native_file, native_pose, resnames.size() );
			//native_pose->dump_pdb( "native.pdb" ); //DELETE ME
		} else {
			TR << "No native structure specified by the user.  No RMSD values will be calculated." << std::endl;
		}

		//Set up a filter for total number of hbonds:
		protocols::filters::CombinedFilterOP total_hbond( new protocols::filters::CombinedFilter );
		set_up_hbond_filter( total_hbond, resnames.size(), sfxn_default, static_cast<core::Real>( option[min_genkic_hbonds]() ) );

		//EVERYTHING ABOVE THIS POINT IS DONE ONCE PER PROGRAM EXECUTION.
		core::Size success_count(0);

		for ( core::Size irepeat=1, irepeat_max=static_cast<core::Size>(option[out::nstruct]()); irepeat<=irepeat_max; ++irepeat ) { //Loop nstruct times

			//Cyclic permutation of sequence.
			core::Size cyclic_offset(0);
			utility::vector1 < std::string > resnames_copy;
			if ( option[cyclic_permutations].value() ) {
				cyclic_offset = do_cyclic_permutation( resnames, resnames_copy );
			}
			runtime_assert(cyclic_offset < resnames_copy.size() ); //Should be true.

			//Create the pose:
			core::pose::PoseOP pose( new core::pose::Pose );
			build_polymer(pose, resnames_copy);

			//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
			protocols::cyclic_peptide::DeclareBondOP termini( new protocols::cyclic_peptide::DeclareBond );
			set_up_termini_mover( termini, pose );
			termini->apply(*pose);

			//Add cyclic constraints:
			add_cyclic_constraints(pose);

			//Set all omega values to 180 and randomize mainchain torsions:
			set_mainchain_torsions(pose);

			//Do the kinematic closure:
			bool const success( genkic_close(pose, sfxn_highhbond_cst, total_hbond) );

			if ( !success ) {
				TR << "Closure failed.";
				if ( irepeat < irepeat_max ) TR << "  Continuing to next job.";
				TR << std::endl;
				TR.flush();
				continue;
			}

			//If we reach here, then closure was successful.  Time to relax the pose.

			TR << "Closure successful." << std::endl;
			protocols::relax::FastRelaxOP frlx( new protocols::relax::FastRelax(sfxn_default_cst, 1) );
			for ( core::Size i=1, imax=static_cast<core::Size>(option[fast_relax_rounds]()); i<=imax; ++i ) {
				TR << "Applying final FastRelax, round " << i << "." << std::endl;
				frlx->apply( *pose );
			}

			//Undo the cyclic permutation in anticipation of re-aligning to the native:
			depermute( pose, cyclic_offset );

			core::Real native_rmsd(0.0);

			//Score the pose before output:
			(*sfxn_default)(*pose);

			//Re-filter based on number of Hbonds (using option[min_final_hbonds]()):
			core::Real const final_hbonds( total_hbond->compute( *pose ) );
			if ( final_hbonds > -1.0*static_cast<core::Real>( option[min_final_hbonds]() ) ) {
				TR << "Final hbond count is " << -1.0*final_hbonds << ", which is less than the minimum.  Failing job." << std::endl;
				continue;
			}

			++success_count; //Increment the count of number of successes.

			if ( native_pose ) {
				native_rmsd = align_and_calculate_rmsd(pose, native_pose);
			}
			TR << "Result\tRMSD\tEnergy\tHbonds" << std::endl;
			TR << irepeat << "\t";
			if ( native_pose ) { TR << native_rmsd; }
			else { TR << "--"; }
			TR << "\t" << pose->energies().total_energy() << "\t" << -1.0*final_hbonds << std::endl;

			if ( silent_out ) {
				core::io::silent::SilentFileDataOP silent_file (new io::silent::SilentFileData );
				silent_file->set_filename( out_filename );
				core::io::silent::SilentStructOP ss( io::silent::SilentStructFactory::get_instance()->get_silent_struct_out() );
				char tag[512];
				sprintf(tag, "result_%04lu", static_cast<unsigned long>(irepeat) );
				ss->fill_struct( *pose, std::string(tag) );
				silent_file->write_silent_struct( *ss, out_filename );
			} else { //if pdb output
				char outstring[512];
				sprintf(outstring, "%s%04lu.pdb", out_filename.c_str(), static_cast<unsigned long>(irepeat) );
				pose->dump_scored_pdb( std::string(outstring), *sfxn_default );
			}

		} //Looping through nstruct

		TR << option[out::nstruct]() << " jobs attempted.  " << success_count << " jobs returned solutions." << std::endl;

	} catch ( utility::excn::EXCN_Base& excn ) {
		TR.Error << "Exception caught: " << std::endl;
		excn.show( TR.Error );
		TR.Error.flush();
		return -1;
	}

	if ( TR.visible() ) {
		TR << "Finished simple_cycpep_predict.cc.  Exiting." << std::endl;
		TR.flush();
	}

	return 0;
}

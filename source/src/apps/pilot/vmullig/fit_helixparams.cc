// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file fit_helixparams.cc
/// @brief Constructs a piece of secondary structure with repeating mainchain torsion values,
/// and fits the Crick parameters for a helix to the structure.
/// @details Updated 29 May 2015 to permit the repeating unit to be more than one residue.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

//General includes
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
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

//Application-specific includes
#include <numeric/crick_equations/HelixParams.hh>
#include <protocols/helical_bundle/FitSimpleHelix.hh>
#include <protocols/helical_bundle/util.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//Tracer:
static basic::Tracer TR( "apps.pilot.vmullig.fit_helixparams" );

//Options (ugh -- global variables):
OPT_KEY (StringVector, residue_type)
OPT_KEY (Integer, repeats)
OPT_KEY (RealVector, mainchain_torsions)
OPT_KEY (String, min_type)
OPT_KEY (Real, min_tolerance)
OPT_KEY (Real, rms_offset)
OPT_KEY (Real, r1_guess)
OPT_KEY (Real, omega1_guess)
OPT_KEY (Real, dz1_guess)
OPT_KEY (String, reference_atom)
OPT_KEY (Integer, reference_residue)

OPT_KEY (RealVector, r1_guesses)
OPT_KEY (RealVector, delta_omega1_guesses)
OPT_KEY (RealVector, delta_z1_guesses)

OPT_KEY (StringVector, nonideal_angles)
OPT_KEY (StringVector, nonideal_bondlengths)

/// @brief Set up the options for this pilot app.
void register_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//The ideal phi, psi, and omega angles of a standard left-handed alpha helix:
	utility::vector1<core::Real> helix_phipsiomega;
	helix_phipsiomega.push_back(-64.8); //phi
	helix_phipsiomega.push_back(-41.0); //psi
	helix_phipsiomega.push_back(180.0); //omega

	//An empty real vector.
	utility::vector1 <core::Real> emptyreal;

	//An empty string vector.
	utility::vector1 <std::string> emptystring;

	//The default vector of residue types
	utility::vector1<std::string> alavect;
	alavect.push_back("ALA");

	NEW_OPT( residue_type, "What type of residue (full name) should the polymer be made of?  Default ALA.  Multiple residues can be specified if the repeating unit in this helix is more than one residue.", alavect);
	NEW_OPT( repeats, "How many repeating units long should the polymer used for fitting be?  Default 21.  Note that two additional residues will be added at the start and end, which will not be used in the fit.  Note also that this is the number of repeat units.  That is, if there is more than one residue in the repeating unit, the number of residues in total will be the number of residues per repeating unit times this number (plus two at the end).", 21);
	NEW_OPT( mainchain_torsions, "A list of mainchian dihedral values that will repeat for every residue in the polymer.  Default (-64.8, -41.0, 180.0), corresponding to phi, psi, and omega of an ideal left-handed alpha helix.  If specified, a value must be specified for every mainchain dihedral angle in the residue type from which the polymer will be built.", helix_phipsiomega);
	NEW_OPT( min_type, "The minimization type that will be used (default \"dfpmin\").", "dfpmin");
	NEW_OPT( min_tolerance, "The minimization tolerance that will be used (the default value, 1E-7, is a good place to start).", 0.0000001);
	NEW_OPT( rms_offset, "A small value used to offset certain numbers in the RMS calculation to avoid zero determinants.  Default 1.0e-5.  Raise this if you get bad geometry from a fit.", 1.0e-5);
	NEW_OPT( r1_guess, "Initial guess for the value of r1 of the reference atom (the helix radius, in Angstroms).  Default 1.5.", 1.5);
	NEW_OPT( omega1_guess, "Initial guess for the value of omega1 (the turn per residue, in radians).  Default 1.7.", 1.7);
	NEW_OPT( dz1_guess, "Initial guess for the value of dz1 (the rise per residue, in Angstroms).  Default 1.5.", 1.5);
	NEW_OPT( reference_atom, "The mainchain atom that will be fit first, and used as the reference against wish delta-z and delta-omega are computed.  (Default \"CA\").", "CA");
	NEW_OPT( reference_residue, "If there is more than one residue per repeat unit in the helix, this is the residue containing the mainchain atom that will be fit first, and used as the reference against wish delta-z and delta-omega are computed.  (Default 1).", 1);

	NEW_OPT( r1_guesses, "Initial guesses for the value of r1 of all atoms (the helix radius, in Angstroms).  The guess for the reference atom will be disregarded (r1_guess is used instead).  Unused if not specified, but if specified, one guess must be provided for each atom.", emptyreal);
	NEW_OPT( delta_omega1_guesses, "Initial guesses for the value of delta_omega1 of all atoms (the rotational offset, in Angstroms).  The guess for the reference atom will be disregarded, since delta_omega1 is always zero for the reference atom.  Unused if not specified, but if specified, one guess must be provided for each atom.", emptyreal);
	NEW_OPT( delta_z1_guesses, "Initial guesses for the value of delta_z1 of all atoms (the axial offset, in Angstroms).  The guess for the reference atom will be disregarded, since delta_z1 is always zero for the reference atom.  Unused if not specified, but if specified, one guess must be provided for each atom.", emptyreal);
	NEW_OPT( nonideal_angles, "Mainchain bond angles that are nonideal.  These must be specified as a list of the form \"<resnum1> <first_atomname> <resnum2> <second_atomname> <resnum3> <third_atomname> <angle>\".  For example, if we were working with alpha-amino acids and we wanted to specify that the N-CA-C bond angle for the second residue in the repeating unit was 110 degrees, we would use \"-nonideal_angles 2 N 2 CA 2 C 110.0\".  Not used if not specified.", emptystring);
	NEW_OPT( nonideal_bondlengths, "Mainchain bond lengths that are nonideal.  These must be specified as a list of the form \"<resnum1> <preceding_atomname> <resnum2> <following_atomname> <length>\".  For example, if we were working with alpha-amino acids and we wanted to specify that the CA-C bond length for the second residue in the repeating unit was 1.3 Angstroms, we would use \"-nonideal_bondlengths 2 CA 2 C 1.3\".  Not used if not specified.",  emptystring);

	return;
}


/// @brief Actually build the geometry that we'll be working with.
///
void build_polymer(
	core::pose::Pose &pose,
	utility::vector1<std::string> const &restypes,
	core::Size const residues_per_repeat,
	core::Size const repeats
)
{
	using namespace protocols::cyclic_peptide;

	TR << "Building " << repeats << " repeat polymer (" << residues_per_repeat*repeats << " residues) from ";
	for ( core::Size i=1; i<=residues_per_repeat; ++i ) {
		TR << restypes[i] << (i<residues_per_repeat ? ", " : " ");
	}
	TR << "building blocks." << std::endl;

	PeptideStubMover stubmover;

	stubmover.set_reset_mode(true);
	stubmover.reset_mover_data();
	stubmover.add_residue ("Append", "GLY", 0, false, "", 1, 0, "");
	for ( core::Size i=1; i<=repeats; ++i ) {
		for ( core::Size j=1; j<=residues_per_repeat; ++j ) {
			stubmover.add_residue( "Append", restypes[j], 0, false, "", 1, 0, "" );
		}
	}

	stubmover.add_residue ("Append", "GLY", 0, false, "", 1, 0, "");

	stubmover.apply(pose);

	return;
}


/// @brief Sets the polymer conformation to the repeating conformation (helix) specified with
/// the option[mainchain_torsions]() vector.
void set_pose_conformation(
	core::pose::Pose &pose,
	core::Size const total_repeats,
	core::Size const residues_per_repeat,
	utility::vector1 <core::Real> const &mainchain_torsions,
	utility::vector1< std::pair< utility::vector1< core::id::AtomID >, core::Real > > const &nonstandard_bondangle_list,
	utility::vector1< std::pair< utility::vector1< core::id::AtomID >, core::Real > > const &nonstandard_bondlength_list
)
{
	using namespace core;
	using namespace core::id;

	core::Size const itorsmax=mainchain_torsions.size();
	core::Size ir(1); //Current residue.  Start on residue 1, because we've got an extra residue at each terminus.

	{ //Scope 1: Check that the right number of mainchain torsions have been specified.
		core::Size counter(0);
		for ( core::Size i=1; i<=residues_per_repeat; ++i ) {
			++ir;
			counter += pose.residue(ir).mainchain_torsions().size();
		}
		runtime_assert_string_msg( counter==itorsmax, "User input error: the number of mainchain torsions specified does not equal the number in the repeating unit." );
	}

	ir=1;

	for ( core::Size irepeat=1; irepeat<=total_repeats; ++irepeat ) {
		++ir;

		//Set nonstandard bond angles and bond lengths:
		for ( core::Size j=1, jmax=nonstandard_bondangle_list.size(); j<=jmax; ++j ) {
			runtime_assert_string_msg( nonstandard_bondangle_list[j].first.size() == 3, "Internal program error: nonstandard_bondangle_list vectors are the wrong size." );
			core::id::AtomID at1( nonstandard_bondangle_list[j].first[1].atomno(), nonstandard_bondangle_list[j].first[1].rsd() + (irepeat - 1) * residues_per_repeat + 1 );
			core::id::AtomID at2( nonstandard_bondangle_list[j].first[2].atomno(), nonstandard_bondangle_list[j].first[2].rsd() + (irepeat - 1) * residues_per_repeat + 1 );
			core::id::AtomID at3( nonstandard_bondangle_list[j].first[3].atomno(), nonstandard_bondangle_list[j].first[3].rsd() + (irepeat - 1) * residues_per_repeat + 1 );
			pose.conformation().set_bond_angle( at1, at2, at3, nonstandard_bondangle_list[j].second );
		}
		pose.update_residue_neighbors();
		for ( core::Size j=1, jmax=nonstandard_bondlength_list.size(); j<=jmax; ++j ) {
			runtime_assert_string_msg( nonstandard_bondlength_list[j].first.size() == 2, "Internal program error: nonstandard_bondlength_list vectors are the wrong size." );
			core::id::AtomID at1( nonstandard_bondlength_list[j].first[1].atomno(), nonstandard_bondlength_list[j].first[1].rsd() + (irepeat - 1) * residues_per_repeat + 1 );
			core::id::AtomID at2( nonstandard_bondlength_list[j].first[2].atomno(), nonstandard_bondlength_list[j].first[2].rsd() + (irepeat - 1) * residues_per_repeat + 1 );
			pose.conformation().set_bond_length( at1, at2, nonstandard_bondlength_list[j].second );
		}
		pose.update_residue_neighbors();

		//Set torsions:
		core::Size itors_list(0); //Index of the current torsion in the torsion list.
		core::Size itors_pose(0); //Index of the current torsion in the pose.
		while ( itors_list < itorsmax ) {
			++itors_list;
			++itors_pose;
			if ( itors_pose > pose.residue(ir).mainchain_torsions().size() ) { //We've set all mainchain torsions in the current resiude, so go on to next.
				++ir;
				itors_pose=0;
				--itors_list; //Because this will be incremented twice, otherwise.
				continue;
			}
			pose.set_torsion( TorsionID( ir, id::BB, itors_pose ), mainchain_torsions[ itors_list ] );
		}
	}
	pose.update_residue_neighbors();

	//pose.dump_pdb( "temp.pdb" ); //DELETE ME

	return;
}


/// @brief Add chains of Cu (a convenient single-atom residue) to the pose in the ideal positions of each
/// mainchain atom, based on the fitted Crick parameters.
void add_Cu_chains(
	core::pose::Pose &pose,
	utility::vector1 < core::Real > const &r1_vals,
	core::Real const &omega1_val,
	core::Real const &z1_val,
	utility::vector1 < core::Real > const &delta_omega1_vals,
	utility::vector1 < core::Real > const &delta_z1_vals,
	core::Size const rescount, //Actually the repeat count.
	core::Size const res_per_repeat
) {
	using namespace protocols::cyclic_peptide;

	//A pose that will just be a bunch of Cu residues
	core::pose::Pose cupose;

	//Count the number of mainchain atoms that we'll be placing coppers over:
	core::Size nres(pose.n_residue() - 2);
	core::Size natoms(0);
	for ( core::Size ir=2, irmax=res_per_repeat+1; ir<=irmax; ++ir ) {
		natoms += pose.residue(ir).n_mainchain_atoms();
	}

	//Add the Cu residues:
	PeptideStubMover stubmover;
	stubmover.set_reset_mode(true);
	stubmover.reset_mover_data();
	for ( core::Size ir=1; ir<=rescount; ++ir ) {
		stubmover.add_residue (
			"Append",
			"CU",
			0,
			true,
			"",
			1,
			0,
			""
		);
	}
	stubmover.apply(cupose);

	core::Size atomno(0);
	for ( core::Size ir=2, irmax=res_per_repeat+1; ir<=irmax; ++ir ) {
		for ( core::Size ia=1, iamax=pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia ) {
			++atomno;
			core::Real t = -1.0*static_cast<core::Real>(nres)/2.0 + ir-2;
			for ( core::Size irepeat=1; irepeat<=rescount; ++irepeat ) {
				cupose.set_xyz( core::id::AtomID( 1, irepeat), numeric::crick_equations::xyz( r1_vals[atomno], omega1_val, t, z1_val, delta_omega1_vals[atomno], delta_z1_vals[atomno] ) );
				t+=res_per_repeat;
			}
			pose.append_pose_by_jump(cupose, 1);
		}
	}

	return;
}

/// @brief Parse the user flags for nonstandard bond angles, and convert this into a data object
/// that can be passed to the pose conformation setup function.
void parse_nonstandard_angles(
	utility::vector1< std::pair< utility::vector1< core::id::AtomID >, core::Real > > &nonstandard_angle_list,
	core::pose::Pose const &pose,
	core::Size const residues_per_repeat
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	nonstandard_angle_list.clear();

	std::stringstream ss("");
	core::Size const noptions(option[nonideal_angles]().size());

	for ( core::Size i=1; i<=noptions; ++i ) {
		ss << option[nonideal_angles]()[i] << " ";
	}

	for ( core::Size i=1; i<=noptions; i+=7 ) {
		core::Size res1;
		ss >> res1;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing angle list!" );
		runtime_assert_string_msg( res1 > 0 && res1 <= residues_per_repeat, "Error in parsing angle list: residue indices must lie in the repeating unit." );

		std::string at1str("");
		ss >> at1str;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing angle list!" );

		core::Size res2;
		ss >> res2;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing angle list!" );
		runtime_assert_string_msg( res2 > 0 && res2 <= residues_per_repeat, "Error in parsing angle list: residue indices must lie in the repeating unit." );

		std::string at2str("");
		ss >> at2str;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing angle list!" );

		core::Size res3;
		ss >> res3;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing angle list!" );
		runtime_assert_string_msg( res3 > 0 && res3 <= residues_per_repeat, "Error in parsing angle list: residue indices must lie in the repeating unit." );

		std::string at3str("");
		ss >> at3str;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing angle list!" );

		core::Real angleval(0.0);
		ss >> angleval;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing angle list!" );

		core::id::AtomID id1( pose.residue(res1+1).type().atom_index( at1str ), res1 );
		core::id::AtomID id2( pose.residue(res2+1).type().atom_index( at2str ), res2 );
		core::id::AtomID id3( pose.residue(res3+1).type().atom_index( at3str ), res3 );
		utility::vector1< core::id::AtomID > idvect;
		idvect.push_back(id1);
		idvect.push_back(id2);
		idvect.push_back(id3);
		core::Real const angleval_radians( numeric::conversions::radians(angleval) );

		nonstandard_angle_list.push_back( std::pair< utility::vector1<core::id::AtomID>, core::Real >(idvect, angleval_radians) );
	}

	return;
}

/// @brief Parse the user flags for nonstandard bond lengths, and convert this into a data object
/// that can be passed to the pose conformation setup function.
void parse_nonstandard_bondlengths(
	utility::vector1< std::pair< utility::vector1< core::id::AtomID >, core::Real > > &nonstandard_bondlength_list,
	core::pose::Pose const &pose,
	core::Size const residues_per_repeat
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	nonstandard_bondlength_list.clear();

	std::stringstream ss("");
	core::Size const noptions(option[nonideal_bondlengths]().size());

	for ( core::Size i=1; i<=noptions; ++i ) {
		ss << option[nonideal_bondlengths]()[i] << " ";
	}

	for ( core::Size i=1; i<=noptions; i+=5 ) {
		core::Size res1;
		ss >> res1;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing bondlength list!" );
		runtime_assert_string_msg( res1 > 0 && res1 <= residues_per_repeat, "Error in parsing bondlength list: residue indices must lie in the repeating unit." );

		std::string at1str("");
		ss >> at1str;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing bondlength list!" );

		core::Size res2;
		ss >> res2;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing bondlength list!" );
		runtime_assert_string_msg( res2 > 0 && res2 <= residues_per_repeat, "Error in parsing bondlength list: residue indices must lie in the repeating unit." );

		std::string at2str("");
		ss >> at2str;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing angle list!" );

		core::Real lengthval(0.0);
		ss >> lengthval;
		runtime_assert_string_msg( !ss.fail(), "Error in parsing angle list!" );

		core::id::AtomID id1( pose.residue(res1+1).type().atom_index( at1str ), res1 );
		core::id::AtomID id2( pose.residue(res2+1).type().atom_index( at2str ), res2 );
		utility::vector1< core::id::AtomID > idvect;
		idvect.push_back(id1);
		idvect.push_back(id2);

		nonstandard_bondlength_list.push_back( std::pair< utility::vector1<core::id::AtomID>, core::Real >(idvect, lengthval) );
	}

	return;
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
			TR << "Starting fit_helixparams.cc" << std::endl;
			TR << "Pilot app created 23 October 2014 by Vikram K. Mulligan, Ph.D., Baker laboratory." << std::endl;
			TR << "App updated 29 May 2015 to support multiple residues per repeating unit." << std::endl;
			TR << "For questions, contact vmullig@uw.edu." << std::endl << std::endl;
		}

		//Copy user inputs from the options system:
		utility::vector1<std::string> const restypes( option[residue_type]() );
		core::Size const residues_per_repeat( restypes.size() ); //How many residues are there per repeating unit?
		runtime_assert_string_msg( option[repeats]()>0, "Error in user input from command-line flags.  The number of repeats must be greater than zero." );
		core::Size const nrepeats( static_cast<core::Size>(option[repeats]()) );
		utility::vector1<core::Real> const mainchaintorsions( option[mainchain_torsions]() );

		utility::vector1<core::Real> const r1guesses( option[r1_guesses]() );
		utility::vector1<core::Real> const deltaomega1_guesses( option[delta_omega1_guesses]() );
		utility::vector1<core::Real> const deltaz1_guesses( option[delta_z1_guesses]() );

		utility::vector1< std::pair< utility::vector1< core::id::AtomID >, core::Real > > nonstandard_angle_list;
		utility::vector1< std::pair< utility::vector1< core::id::AtomID >, core::Real > > nonstandard_bondlength_list;

		core::pose::Pose pose; //Make the empty pose

		build_polymer(pose, restypes, residues_per_repeat, nrepeats); //Build the repeating secondary structure element

		//Set up the AtomIDs for the nonstandard bond angles and bond lengths:
		parse_nonstandard_angles( nonstandard_angle_list, pose, residues_per_repeat );
		parse_nonstandard_bondlengths( nonstandard_bondlength_list, pose, residues_per_repeat );

		set_pose_conformation(pose, nrepeats, residues_per_repeat, mainchaintorsions, nonstandard_angle_list, nonstandard_bondlength_list); //Set the pose conformation.

		//Count atoms:
		core::Size atomcount(0);
		for ( core::Size ir=2, irmax=residues_per_repeat+1; ir<=irmax; ++ir ) {
			atomcount += pose.residue(ir).n_mainchain_atoms();
		}
		runtime_assert_string_msg( r1guesses.size()==0 || r1guesses.size()==atomcount, "Error in user input from command-line flags.  If specified, the number of r1 guesses must equal the number of atoms." );
		runtime_assert_string_msg( deltaomega1_guesses.size()==0 || deltaomega1_guesses.size()==atomcount, "Error in user input from command-line flags.  If specified, the number of delta_omega1 guesses must equal the number of atoms." );
		runtime_assert_string_msg( deltaz1_guesses.size()==0 || deltaz1_guesses.size()==atomcount, "Error in user input from command-line flags.  If specified, the number of delta_z1 guesses must equal the number of atoms." );

		protocols::helical_bundle::FitSimpleHelix fitter;
		fitter.set_range(2,  pose.n_residue() - 1);
		fitter.set_min_type( option[min_type]() );
		fitter.set_min_tolerance( option[min_tolerance]());
		fitter.set_rms_offset( option[rms_offset]() );
		fitter.set_initial_guesses( option[r1_guess](), option[omega1_guess](), option[dz1_guess]() );
		fitter.set_reference_atom( option[reference_atom]() );
		fitter.set_residues_per_repeat( residues_per_repeat );
		fitter.set_reference_residue( static_cast<core::Size>(option[reference_residue]()) );

		if ( r1guesses.size() > 0 ) fitter.set_r1_guesses( r1guesses );
		if ( deltaomega1_guesses.size() > 0 ) fitter.set_delta_omega1_guesses( deltaomega1_guesses );
		if ( deltaz1_guesses.size() > 0 ) fitter.set_delta_z1_guesses( deltaz1_guesses );

		fitter.apply(pose); //Will shift the pose to match the ideal helix

		//Storage spots for the results of the fit:
		utility::vector1 < core::Real > r1_vals; //Radius for each mainchain atom
		utility::vector1 < core::Real > delta_z1_vals; //Z-offset for each mainchain atom
		utility::vector1 < core::Real > delta_omega1_vals; //Omega offset for each mainchain atom
		core::Real z1_val = 0.0; //Rise per residue, Angstroms
		core::Real omega1_val = 0.0; //Turn per residue, radians

		//Retrieve the results of the fit:
		fitter.get_crick_parameters( r1_vals, omega1_val, z1_val, delta_omega1_vals, delta_z1_vals );

		//Add chains B, C, D, etc. -- chains of Cu (a convenient single-atom residue) overlaying on each of the mainchain atoms.
		add_Cu_chains(pose, r1_vals, omega1_val, z1_val, delta_omega1_vals, delta_z1_vals, nrepeats, residues_per_repeat );

		if ( TR.visible() ) TR << "Writing result.pdb." << std::endl;
		pose.dump_pdb("result.pdb"); //dump out a pose.

		if ( TR.visible() ) TR << "Writing result.crick_params" << std::endl;
		if ( residues_per_repeat==1 ) protocols::helical_bundle::write_minor_helix_params("result.crick_params", r1_vals, omega1_val, z1_val, delta_omega1_vals, delta_z1_vals);
		else {
			utility::vector1 <core::Size> atoms_per_residue;
			core::Size ir=1;
			for ( core::Size irepeat=1; irepeat<=residues_per_repeat; ++irepeat ) {
				++ir;
				atoms_per_residue.push_back( pose.residue(ir).n_mainchain_atoms() );
			}
			protocols::helical_bundle::write_minor_helix_params("result.crick_params", residues_per_repeat, atoms_per_residue, r1_vals, omega1_val, z1_val, delta_omega1_vals, delta_z1_vals);
		}

		//Write out the Crick parameters
		if ( TR.visible() ) {
			TR << std::endl;
			TR << "Results from the fitter:" << std::endl;
			TR.width(20);
			TR.precision(8);
			TR << "Res\tID\tAtom\tr1\tomega1\tz1\tdelta_omega1\tdelta_z1\t" << std::endl;
			core::Size completed_residues(0);
			core::Size counter(0);
			for ( core::Size ir=2; completed_residues<residues_per_repeat; ++ir ) {
				for ( core::Size ia=1, iamax=pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia ) {
					++counter;
					TR << ir-1 << "\t" << ia << "\t" << pose.residue(ir).atom_name(ia);
					if ( (ir-1)==static_cast<core::Size>(option[reference_residue]()) && ia==pose.residue(ir).atom_index(option[reference_atom]()) ) TR << "*";
					TR << "\t" << r1_vals[counter] << "\t" << omega1_val << "\t" << z1_val << "\t" << delta_omega1_vals[counter] << "\t" << delta_z1_vals[counter] << std::endl;
				}
				++completed_residues;
			}
			TR <<"(*Reference atom.)" << std::endl;
			TR << std::endl;
		}

		if ( TR.visible() ) {
			TR << "Finished fit_helixparams.cc.  Exiting." << std::endl;
			TR.flush();
		}

	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return 0;
}

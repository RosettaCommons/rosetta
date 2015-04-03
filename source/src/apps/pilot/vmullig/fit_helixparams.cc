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
/// @author Vikram K. Mulligan (vmullig@uw.edu)

//General includes
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>

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
OPT_KEY (String, residue_type)
OPT_KEY (Integer, residue_repeats)
OPT_KEY (RealVector, mainchain_torsions)
OPT_KEY (String, min_type)
OPT_KEY (Real, min_tolerance)
OPT_KEY (Real, r1_guess)
OPT_KEY (Real, omega1_guess)
OPT_KEY (Real, dz1_guess)
OPT_KEY (String, reference_atom)


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

	NEW_OPT( residue_type, "What type of residue (full name) should the polymer be made of?  Default ALA.", "ALA");
	NEW_OPT( residue_repeats, "How many subunits long should the polymer used for fitting be?  Default 21.  Note that two additional residues will be added at the start and end, which will not be used in the fit.", 21);
	NEW_OPT( mainchain_torsions, "A list of mainchian dihedral values that will repeat for every residue in the polymer.  Default (-64.8, -41.0, 180.0), corresponding to phi, psi, and omega of an ideal left-handed alpha helix.  If specified, a value must be specified for every mainchain dihedral angle in the residue type from which the polymer will be built.", helix_phipsiomega);
	NEW_OPT( min_type, "The minimization type that will be used (default \"dfpmin\").", "dfpmin");
	NEW_OPT( min_tolerance, "The minimization tolerance that will be used (the default value, 1E-7, is a good place to start).", 0.0000001);
	NEW_OPT( r1_guess, "Initial guess for the value of r1 (the helix radius, in Angstroms).  Default 1.5.", 1.5);
	NEW_OPT( omega1_guess, "Initial guess for the value of omega1 (the turn per residue, in radians).  Default 1.7.", 1.7);
	NEW_OPT( dz1_guess, "Initial guess for the value of dz1 (the rise per residue, in Angstroms).  Default 1.5.", 1.5);
	NEW_OPT( reference_atom, "The mainchain atom that will be fit first, and used as the reference against wish delta-z and delta-omega are computed.  (Default \"CA\").", "CA");

	return;
}


/// @brief Actually build the geometry that we'll be working with.
void build_polymer( core::pose::Pose &pose )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::cyclic_peptide;

	TR << "Building " << option[residue_repeats]() << " residue polymer from " << option[residue_type]() << " building blocks." << std::endl; 

	PeptideStubMover stubmover;

	stubmover.set_reset_mode(true);
	stubmover.reset_mover_data();
	stubmover.add_residue ("Append", "GLY", 0, false, "", 1, 0, "");
	stubmover.add_residue (
		"Append",
		option[residue_type](),
		0,
		false,
		"",
		static_cast<core::Size>( option[residue_repeats]() ),
		0,
		""		
	);
	stubmover.add_residue ("Append", "GLY", 0, false, "", 1, 0, "");

	stubmover.apply(pose);

	return;
}


/// @brief Sets the polymer conformation to the repeating conformation (helix) specified with
/// the option[mainchain_torsions]() vector.
void set_pose_conformation( core::pose::Pose &pose)
{
	using namespace core;
	using namespace core::id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Size const itorsmax=option[mainchain_torsions]().size();

	for(core::Size ir=2,irmax=pose.n_residue()-1; ir<=irmax; ++ir) { //Loop through all pose residues
		runtime_assert_string_msg( itorsmax == pose.residue(ir).mainchain_torsions().size(), "In set_pose_conformation: the number of mainchain torsions must match the number of values provided with the -mainchain_torsions flag." );

		for(core::Size itors=1; itors<=itorsmax; ++itors) {
			pose.set_torsion( TorsionID( ir, id::BB, itors ), static_cast<core::Real>( option[mainchain_torsions]()[itors] ) );
		}

	}

	pose.update_residue_neighbors();

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
	std::string const &/*ref_atom*/,
	core::Size const rescount
) {
	using namespace protocols::cyclic_peptide;

	//A pose that will just be a bunch of Cu residues
	core::pose::Pose cupose;

	//Add the Cu residues:
	PeptideStubMover stubmover;
	stubmover.set_reset_mode(true);
	stubmover.reset_mover_data();
	for(core::Size ir=1; ir<=rescount; ++ir) {
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
	
	for(core::Size ia=1, iamax=pose.residue(2).n_mainchain_atoms(); ia<=iamax; ++ia) {
		core::Real t = -1.0*static_cast<core::Real>(rescount)/2.0;
		for(core::Size ir=1; ir<=rescount; ++ir) {
			cupose.set_xyz( core::id::AtomID( 1, ir), numeric::crick_equations::xyz( r1_vals[ia], omega1_val, t, z1_val, delta_omega1_vals[ia], delta_z1_vals[ia] ) );
			t+=1.0;
		}
		pose.append_pose_by_jump(cupose, 1);
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

		if(TR.visible()) {
			TR << "Starting fit_helixparams.cc" << std::endl;
			TR << "Pilot app created 23 October 2014 by Vikram K. Mulligan, Ph.D., Baker laboratory." << std::endl;
			TR << "For questions, contact vmullig@uw.edu." << std::endl << std::endl;
		}

		core::pose::Pose pose; //Make the empty pose

		build_polymer(pose); //Build the repeating secondary structure element
		set_pose_conformation(pose); //Set the pose conformation.

		protocols::helical_bundle::FitSimpleHelix fitter;
		fitter.set_range(2,  static_cast<core::Size>(option[residue_repeats]())+1);
		fitter.set_min_type( option[min_type]() );
		fitter.set_min_tolerance( option[min_tolerance]());
		fitter.set_initial_guesses( option[r1_guess](), option[omega1_guess](), option[dz1_guess]() );
		fitter.set_reference_atom( option[reference_atom]() );
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
		add_Cu_chains(pose, r1_vals, omega1_val, z1_val, delta_omega1_vals, delta_z1_vals, option[reference_atom](), static_cast<core::Size>(option[residue_repeats]()) );

		if(TR.visible()) TR << "Writing result.pdb." << std::endl;
		pose.dump_pdb("result.pdb"); //dump out a pose.

		if(TR.visible()) TR << "Writing result.crick_params" << std::endl;
		protocols::helical_bundle::write_minor_helix_params("result.crick_params", r1_vals, omega1_val, z1_val, delta_omega1_vals, delta_z1_vals);

		//Write out the Crick parameters
		if(TR.visible()) {
			TR << std::endl;
			TR << "Results from the fitter:" << std::endl;
			TR.width(20);
			TR.precision(8);
			TR << "ID\tAtom\tr1\tomega1\tz1\tdelta_omega1\tdelta_z1\t" << std::endl;
			for(core::Size ia=1, iamax=pose.residue(2).n_mainchain_atoms(); ia<=iamax; ++ia) {
				TR << ia << "\t" << pose.residue(2).atom_name(ia);
				if( ia==pose.residue(2).atom_index(option[reference_atom]()) ) TR << "*";
				TR << "\t" << r1_vals[ia] << "\t" << omega1_val << "\t" << z1_val << "\t" << delta_omega1_vals[ia] << "\t" << delta_z1_vals[ia] << std::endl;
			}
			TR <<"(*Reference atom.)" << std::endl;
			TR << std::endl;
		}

		if(TR.visible()) {
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

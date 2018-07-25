// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/util.cc
/// @brief  Utility functions for helical bundle construction.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// C++ headers
#include <cstdio>

// Unit Headers
#include <protocols/helical_bundle/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

#include <core/id/TorsionID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/util.tmpl.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/conversions.hh>
#include <numeric/crick_equations/BundleParams.hh>
#include <protocols/rosetta_scripts/util.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


using basic::Error;
using basic::Warning;

namespace protocols {
namespace helical_bundle {

static basic::Tracer TR("protocols.helical_bundle.util");

using namespace utility::tag;

void add_attributes_for_make_bundle_symmetry( AttributeList & attlist ) {
	attlist + XMLSchemaAttribute::attribute_w_default( "symmetry", xsct_non_negative_integer, "Symmetry setting (n-fold; 0 or 1 === no symmetry", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "symmetry_copies", xsct_non_negative_integer, "How many symmetry copies will be generated? 'All' if zero, only the first one if 1, but you can ask for any other number", "0" );
}

void add_attributes_for_make_bundle_dofs( AttributeList & attlist ) {
	attlist + XMLSchemaAttribute::attribute_w_default( "set_bondlengths", xsct_rosetta_bool, "Should bond lengths be set (true) or fixed at ideal values (false)?", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "set_bondangles", xsct_rosetta_bool, "Should bond angles be set (true) or fixed at ideal values (false)?", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "set_dihedrals", xsct_rosetta_bool, "Should dihedrals be set (true) or fixed at ideal values (false)?", "true" );
}

void add_attributes_for_make_bundle_minorhelix_defaults( AttributeList & attlist ) {
	if ( ! attribute_w_name_in_attribute_list( "crick_params_file", attlist ) ) {
		attlist + XMLSchemaAttribute( "crick_params_file", xs_string, "File name of a file containing Crick parameters for the secondary structure type desired." );
	}
}

void add_attributes_for_make_bundle_other_defaults( AttributeList & attlist ) {
	attlist + XMLSchemaAttribute( "residue_name", xs_string, "Residue, indicated by name, from which to build the helical bundle." )
		+ XMLSchemaAttribute::attribute_w_default( "repeating_unit_offset", xsct_non_negative_integer, "Offset between repeating units.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "invert", xsct_rosetta_bool, "Should this helix be flipped?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "helix_length", xsct_non_negative_integer, "Length, in residues, for this helix.", "0" );
}

void add_attributes_for_helix_params( AttributeList & attlist ) {
	// This may already exist if add_attributes_for_minor_helix_params has been called
	if ( ! attribute_w_name_in_attribute_list( "crick_params_file", attlist ) ) {
		attlist + XMLSchemaAttribute( "crick_params_file", xs_string, "File name of a file containing Crick parameters for the secondary structure type desired." );
	}
	// Per-helix control of bond lengths, angles, and dihedrals
	add_attributes_for_make_bundle_dofs( attlist );
}

void add_attributes_for_minor_helix_params( AttributeList & attlist ) {
	if ( ! attribute_w_name_in_attribute_list( "crick_params_file", attlist ) ) {
		attlist + XMLSchemaAttribute( "crick_params_file", xs_string, "File name of a file containing Crick parameters for the secondary structure type desired." );
	}
	attlist + XMLSchemaAttribute::attribute_w_default( "omega1", xsct_real, "Minor helical turn per residue, for a specific helix", "0" );
	if ( ! attribute_w_name_in_attribute_list( "delta_omega1", attlist ) ) {
		attlist + XMLSchemaAttribute::attribute_w_default( "delta_omega1", xsct_real, "Minor helical twist per residue, for a specific helix", "0" );
	}
	if ( ! attribute_w_name_in_attribute_list( "z1", attlist ) ) {
		attlist + XMLSchemaAttribute::attribute_w_default( "z1", xsct_real, "Helix rise per residue, for a specific helix", "0" );
	}
}

void add_attributes_for_other_helix_params( AttributeList & attlist ) {
	attlist + XMLSchemaAttribute( "residue_name", xs_string, "For a specific helix, residue, indicated by name, from which to build the helical bundle." )
		+ XMLSchemaAttribute::attribute_w_default( "repeating_unit_offset", xsct_non_negative_integer, "For a specific helix, Offset between repeating units", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "invert", xsct_rosetta_bool, "For a specific helix, should this helix be flipped?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "helix_length", xsct_non_negative_integer, "For a specific helix, length, in residues, for this helix", "0" );
}

/// @brief Actual write of the crick_params file data.
/// @details Called by both write_minor_helix_params variants.
/// The outfile ozstream must already be opened, and will not be closed
/// by this function.
void write_crick_params_file_data (
	utility::io::ozstream &outfile,
	utility::vector1 < core::Real > const &r1,
	core::Real const &omega1,
	core::Real const &z1,
	utility::vector1 < core::Real > const &delta_omega1,
	utility::vector1 < core::Real > const &delta_z1
) {
	using namespace utility::io;

	char linebuffer[256];
	std::string line_out;

	for ( core::Size i=1, imax=r1.size(); i<=imax; ++i ) {
		sprintf(linebuffer, "r1\t%.12f\n", r1[i]);
		line_out=linebuffer;
		outfile.write( line_out, line_out.length() );
	}

	sprintf(linebuffer, "omega1\t%.12f\nz1\t%.12f\n", omega1, z1);
	line_out=linebuffer;
	outfile.write(line_out, line_out.length() );

	for ( core::Size i=1, imax=delta_omega1.size(); i<=imax; ++i ) {
		sprintf(linebuffer, "delta_omega1\t%.12f\n", delta_omega1[i]);
		line_out=linebuffer;
		outfile.write( line_out, line_out.length() );
	}

	for ( core::Size i=1, imax=delta_z1.size(); i<=imax; ++i ) {
		sprintf(linebuffer, "delta_z1\t%.12f\n", delta_z1[i]);
		line_out=linebuffer;
		outfile.write( line_out, line_out.length() );
	}
}

/// @brief Write out a crick_params file.
/// @details Variant for case of a single-residue repeating unit.
void write_minor_helix_params (
	std::string const &filename,
	utility::vector1 < core::Real > const &r1,
	core::Real const &omega1,
	core::Real const &z1,
	utility::vector1 < core::Real > const &delta_omega1,
	utility::vector1 < core::Real > const &delta_z1
) {
	using namespace utility::io;

	ozstream outfile;
	outfile.open( filename );

	runtime_assert_string_msg( outfile.good(), "In protocols::helical_bundle::write_minor_helix_params: Unable to open file for write!" );

	write_crick_params_file_data( outfile, r1, omega1, z1, delta_omega1, delta_z1 );

	outfile.close();

	return;
}

/// @brief Write out a crick_params file.
/// @details Variant for case of a multi-residue repeating unit.
void write_minor_helix_params (
	std::string const &filename,
	core::Size const &residues_per_repeat,
	utility::vector1 <core::Size> const &atoms_per_residue,
	utility::vector1 < core::Real > const &r1,
	core::Real const &omega1,
	core::Real const &z1,
	utility::vector1 < core::Real > const &delta_omega1,
	utility::vector1 < core::Real > const &delta_z1
) {
	using namespace utility::io;

	ozstream outfile;
	outfile.open( filename );

	runtime_assert_string_msg( outfile.good(), "In protocols::helical_bundle::write_minor_helix_params: Unable to open file for write!" );

	char linebuffer[1024];
	std::string line_out;

	if ( residues_per_repeat > 1 ) {
		sprintf(linebuffer, "residues_per_repeat\t%lu\n", static_cast<unsigned long>( residues_per_repeat ) );
		line_out=linebuffer;
		outfile.write( line_out, line_out.length() );
		runtime_assert_string_msg( atoms_per_residue.size() > 0, "In protocols::helical_bundle::write_minor_helix_params: got an empty vector of residue mainchain atom counts.  This should not happen.  Consult a developer -- it's a programming error, not a user error." );
		for ( core::Size i=1, imax=atoms_per_residue.size(); i<=imax; ++i ) {
			sprintf( linebuffer, "atoms_per_residue\t%lu\t%lu\n", static_cast<unsigned long>( i ), static_cast<unsigned long>( atoms_per_residue[i] ) );
			line_out=linebuffer;
			outfile.write( line_out, line_out.length() );
		}
	}

	write_crick_params_file_data( outfile, r1, omega1, z1, delta_omega1, delta_z1 );

	outfile.close();

	return;
}

/// @brief Read minor helix parameters from a crick_params file.
///
void read_minor_helix_params (
	std::string const &filename,
	utility::vector1 < core::Real > &r1,
	core::Real &omega1,
	core::Real &z1,
	utility::vector1 < core::Real > &delta_omega1,
	utility::vector1 < core::Real > &delta_z1,
	core::Size &residues_per_repeat,
	utility::vector1 <core::Size> &atoms_per_residue
) {
	using namespace utility::io;

	r1.clear();
	delta_omega1.clear();
	delta_z1.clear();
	atoms_per_residue.clear();
	omega1 = 0.0;
	z1 = 0.0;

	std::string filename_formatted = filename;
	if ( utility::file::file_extension(filename_formatted)!="crick_params" ) filename_formatted+= ".crick_params";

	izstream infile;
	infile.open( filename_formatted );
	if ( !infile.good() ) {
		filename_formatted = "protocol_data/crick_parameters/" + utility::file::file_basename(filename_formatted) + ".crick_params";
		basic::database::open( infile, filename_formatted );
		runtime_assert_string_msg( infile.good(), "In protocols::helical_bundle::read_minor_helix_params: Unable to open .crick_params file for read!" );
	}

	if ( TR.Debug.visible() ) TR.Debug << "Reading " << filename_formatted << std::endl;

	bool residues_per_repeat_specified(false);

	while ( true ) {
		std::string current_line = "";
		infile.getline(current_line);//Get the current line and store it in current_line.
		if ( infile.eof() ) break;

		if ( TR.Debug.visible() ) TR.Debug << current_line << std::endl;

		if ( current_line[0] == '#' ) continue; //Ignore lines that start with a pound sign.

		std::stringstream ss(current_line);
		char buffer [25];
		ss.getline(buffer, 25, '\t');
		std::string strbuffer(buffer);
		//TR.Debug << strbuffer << "." << std::endl; //DELETE ME
		if ( strbuffer=="r1" ) {
			ss.getline(buffer, 25);
			r1.push_back( static_cast<core::Real>( atof(buffer) ) );
		} else if ( strbuffer=="delta_omega1" ) {
			ss.getline(buffer, 25);
			delta_omega1.push_back( static_cast<core::Real>( atof(buffer) ) );
		} else if ( strbuffer=="delta_z1" ) {
			ss.getline(buffer, 25);
			delta_z1.push_back( static_cast<core::Real>( atof(buffer) ) );
		} else if ( strbuffer=="omega1" ) {
			ss.getline(buffer, 25);
			omega1=static_cast<core::Real>( atof(buffer) );
		} else if ( strbuffer=="z1" ) {
			ss.getline(buffer, 25);
			z1=static_cast<core::Real>( atof(buffer) );
		} else if ( strbuffer=="residues_per_repeat" ) {
			runtime_assert_string_msg( !residues_per_repeat_specified, "Error in protocols::helical_bundle::read_minor_helix_params(): Only one \"residues_per_repeat\" line should be present in a crick_params file." );
			residues_per_repeat_specified=true;
			ss.getline(buffer, 25);
			residues_per_repeat=static_cast<core::Size>( atoi(buffer) );
			atoms_per_residue.resize( residues_per_repeat, 0 );
		} else if ( strbuffer=="atoms_per_residue" ) {
			runtime_assert_string_msg( residues_per_repeat_specified, "Error in protocols::helical_bundle::read_minor_helix_params():  The \"atoms_per_residue\" lines can only be after a \"residues_per_repeat\" line.  (This is a problem with the crick_params file.)" );
			ss.getline(buffer, 25, '\t');
			auto const curres( static_cast<core::Size>( atoi( buffer ) ) );
			runtime_assert_string_msg( curres > 0 && curres <=residues_per_repeat, "Error in protocols::helical_bundle::read_minor_helix_params():  An \"atoms_per_residue\" line specifies a residue index that's not in the repeating unit.  (This is a problem with the crick_params file.)" );
			ss.getline(buffer, 25);
			auto const curres_atomcount( static_cast<core::Size>( atoi( buffer ) ) );
			runtime_assert_string_msg( curres_atomcount > 0, "Error in protocols::helical_bundle::read_minor_helix_params():  The number of atoms in the current residue must be greater than zero.  (This is a problem with the crick_params file.)" );
			runtime_assert_string_msg( atoms_per_residue[curres] == 0, "Error in protocols::helical_bundle::read_minor_helix_params():  The number of atoms in a residue was specified more than once in the crick_params file.");
			atoms_per_residue[curres] = curres_atomcount;
		} else continue;

	}

	infile.close();

	if ( !residues_per_repeat_specified ) {
		residues_per_repeat=1; //Default to 1 if not specified in the file.
		atoms_per_residue.clear();
		atoms_per_residue.push_back( r1.size() );
		runtime_assert_string_msg( delta_omega1.size() == atoms_per_residue[1] && delta_z1.size() == atoms_per_residue[1],
			"Error in protocols::helical_bundle::read_minor_helix_params():  The number of atoms in the residue is inconsistent, based on the number of r1, delta_omega1, and delta_z1 values provided in the crick_params file." );
	} else {
		core::Size counter(0);
		for ( core::Size i=1, imax=atoms_per_residue.size(); i<=imax; ++i ) {
			runtime_assert_string_msg( atoms_per_residue[i] != 0, "Error in protocols::helical_bundle::read_minor_helix_params():  The number of atoms in at least one of the residues in the repeating unit was not specified in the crick_params file." );
			counter += atoms_per_residue[i];
		}
		runtime_assert_string_msg(
			r1.size() == counter &&
			delta_omega1.size() == counter &&
			delta_z1.size() == counter,
			"Error in protocols::helical_bundle::read_minor_helix_params():  Parameters were not specified for all atoms in all residues in the repeating unit.  Please check the crick_params file provided."
		);
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Finished reading " << filename_formatted << std::endl;
		for ( core::Size i=1, imax=r1.size(); i<=imax; ++i ) {
			TR.Debug << "r1[" << i << "]\t" << r1[i] << std::endl;
		}
		TR.Debug << "omega1\t" << omega1 << std::endl;
		TR.Debug << "z1\t" << z1 << std::endl;
		for ( core::Size i=1, imax=delta_omega1.size(); i<=imax; ++i ) {
			TR.Debug << "delta_omega1[" << i << "]\t" << delta_omega1[i] << std::endl;
		}
		for ( core::Size i=1, imax=delta_z1.size(); i<=imax; ++i ) {
			TR.Debug << "delta_z1[" << i << "]\t" << delta_z1[i] << std::endl;
		}
		TR.Debug << "residues_per_repeat\t" << residues_per_repeat << std::endl;
		for ( core::Size i=1; i<=residues_per_repeat; ++i ) {
			TR.Debug << "atoms in residue " << i << "\t" << atoms_per_residue[i] << std::endl;
		}
	}

	return;
}

/// @brief Generate the x,y,z coordinates of the mainchain atoms using the Crick equations.
/// @details Coordinates will be returned as a vector of vectors of xyzVectors.  The outer
/// index will refer to residue number, and the inner index will refer to atom number.
/// Returns failed=true if coordinates could not be generated, false otherwise.
void generate_atom_positions(
	utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > &outvector,
	core::pose::Pose const &helixpose,
	core::Size const helix_start,
	core::Size const helix_end,
	core::Real const &r0,
	core::Real const &omega0,
	core::Real const &delta_omega0,
	core::Real const &delta_t,
	core::Real const &z1_offset,
	core::Real const &z0_offset,
	core::Real const &epsilon,
	bool const invert_helix,
	utility::vector1 < core::Real > const &r1,
	core::Real const &omega1,
	core::Real const &z1,
	utility::vector1 < core::Real > const &delta_omega1,
	core::Real const &delta_omega1_all,
	utility::vector1 < core::Real > const &delta_z1,
	core::Size const residues_per_repeat,
	utility::vector1 <core::Size> const &atoms_per_residue,
	core::Size const repeating_unit_offset, //0 if the first residue is the first residue in the repeating unit; 1 if we're off by 1, etc.
	bool &failed
) {

	runtime_assert_string_msg( z1 != 0, "Error in protocols::helical_bundle::generate_atom_positions(): The value of z1 cannot be zero!" );

	outvector.clear();
	failed=false;

	core::Size helix_length=helix_end-helix_start+1;

	// The offset along the minor helix axis is converted into delta_t and delta_omega1_all values:
	core::Real const r0_omega0_over_z1( r0*omega0/z1 );
	if ( r0_omega0_over_z1 > 1 ) {
		failed=true;
		outvector.clear();
		return;
	}
	core::Real const z1_offset_prime( z1_offset/sqrt( 1 - r0_omega0_over_z1*r0_omega0_over_z1 ) ); //Correction factor to match the z-offset in Grigoryan and DeGrado (2011) JMB 405(4):1079-1100.
	core::Real const delta_t_prime (delta_t + z1_offset_prime / z1);
	core::Real const delta_omega1_all_prime ( delta_omega1_all - (z1_offset_prime / z1) * omega1 );

	core::Real t = -1.0 * static_cast<core::Real>(helix_length + 2) / 2.0 + delta_t_prime;

	core::Size residue_in_repeating_unit(repeating_unit_offset); //What's the index of the current residue in the repeating unit?
	core::Size curatom_in_repeating_unit(0); //What's the current mainchain atom index in the repeating unit?

	for ( core::Size ir2=helix_start-1; ir2<=helix_end+1; ++ir2 ) { //Loop through all residues in the helix, padded by one on either side
		core::Size ir( ir2 );
		//Repeat start and end residues to pad by one:
		if ( ir2 == helix_start-1 ) ir = helix_start;
		if ( ir2 == helix_end+1 ) ir = helix_end;

		//Figure out which residue we are in the repeating unit:
		++residue_in_repeating_unit;
		if ( ir2==helix_start-1 ) {
			residue_in_repeating_unit=repeating_unit_offset+1;
			curatom_in_repeating_unit=0;
		}
		if ( ir2==helix_end+1 ) { //Repeat the last residue if we're at the end.
			--residue_in_repeating_unit;
			curatom_in_repeating_unit -= helixpose.residue(helix_end).n_mainchain_atoms();
		}
		if ( residue_in_repeating_unit > residues_per_repeat || residue_in_repeating_unit==0 || ir2 == helix_start ) {
			residue_in_repeating_unit=1; //Reset if we're on to the next repeating unit.
			curatom_in_repeating_unit=0;
		}

		utility::vector1 < numeric::xyzVector <core::Real> > innervector;
		for ( core::Size ia=1, iamax=helixpose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia ) { //Loop through all mainchain atoms in the current helix residue

			if ( ia==1 ) {
				runtime_assert_string_msg(
					iamax == atoms_per_residue[ residue_in_repeating_unit ],
					"Error in protocols::helical_bundle::generate_atom_positions():  The number of atoms per residue for one or more residues does not match the expected number given the Crick parameters file."
				);
			}

			//Increment the current atom in the repeating unit:
			++curatom_in_repeating_unit;

			//Actually call the generator function:
			innervector.push_back( numeric::crick_equations::XYZ_BUNDLE(
				t, r0, omega0,
				( invert_helix ? numeric::constants::d::pi - delta_omega0 : delta_omega0 ),
				r1[curatom_in_repeating_unit], omega1-omega0, z1, delta_omega1[curatom_in_repeating_unit]+delta_omega1_all_prime, delta_z1[curatom_in_repeating_unit],
				epsilon,
				failed )
			);
			if ( invert_helix ) {
				innervector[ia].x( -1.0*innervector[ia].x() );
				//innervector[ia].y( -1.0*innervector[ia].y() );
				innervector[ia].z( -1.0*innervector[ia].z() );
			}
			innervector[ia].z( innervector[ia].z() + (invert_helix ? -1.0 : 1.0 ) * z0_offset ); //Offset by z0_offset in direction of helix (inverted for inverted helix)
			if ( failed ) break;

		}
		if ( failed ) break;
		outvector.push_back(innervector);
		t+=1.0;
	}

	if ( failed ) outvector.clear();

	return;
}

/// @brief Place the helix mainchain atoms based on the Crick equations.
///
void place_atom_positions(
	core::pose::Pose &pose,
	utility::vector1 < utility::vector1 < numeric::xyzVector < core::Real >  > > const &atom_positions,
	core::Size const helix_start,
	core::Size const helix_end
) {

	//Index in the outer vector for the current residue.
	core::Size index=2;

	for ( core::Size ir=helix_start; ir<=helix_end; ++ir ) {
		for ( core::Size ia=1, iamax=atom_positions[index].size(); ia<=iamax; ++ia ) {
			pose.set_xyz( core::id::AtomID( pose.residue_type(ir).mainchain_atom(ia), ir ), atom_positions[index][ia] );
		}
		++index;
	}

	pose.update_residue_neighbors();

	return;
}

/// @brief Copy backbone bond length values from one pose, where helix mainchain atom coordinates have been
/// set with the Crick equations, to another with ideal geometry.
void copy_helix_bondlengths(
	core::pose::Pose &pose,
	core::pose::Pose const &ref_pose,
	core::Size const helix_start,
	core::Size const helix_end
) {
	for ( core::Size ir=helix_start; ir<=helix_end; ++ir ) {
		for ( core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia ) {
			if ( ia==1 && ir==helix_start ) continue; //Skip the first atom.
			core::id::AtomID const thisatom( ref_pose.residue_type(ir).mainchain_atom(ia), ir );
			core::id::AtomID const prevatom( (ia==1 ? ref_pose.residue_type(ir - 1).mainchain_atom( ref_pose.residue_type(ir-1).mainchain_atoms().size() ) : ref_pose.residue_type(ir).mainchain_atom( ia - 1 ) ), (ia==1 ? ir-1 : ir) ); //The previous atom is ia-1 in this residue, unless this atom is the first atom, in which case the previous atom is iamax in the previous residue.
			pose.conformation().set_bond_length( thisatom, prevatom, ref_pose.xyz(thisatom).distance( ref_pose.xyz(prevatom) ) );
		}
	}

	//TODO properly handle the first and last residues using the extra residue xyz coordinates that were generated!

	pose.update_residue_neighbors();

	return;
}

/// @brief Copy backbone bond angle values from one pose, where helix mainchain atom coordinates have been
/// set with the Crick equations, to another with ideal geometry.
void copy_helix_bondangles(
	core::pose::Pose &pose,
	core::pose::Pose const &ref_pose,
	core::Size const helix_start,
	core::Size const helix_end
) {
	for ( core::Size ir=helix_start; ir<=helix_end; ++ir ) {
		for ( core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia ) {
			if ( ia==1 && ir==helix_start ) continue; //Skip the first atom.
			if ( ia==iamax && ir==helix_end ) continue; //Skip the last atom.
			core::id::AtomID const thisatom( ref_pose.residue_type(ir).mainchain_atom(ia), ir );
			core::id::AtomID const prevatom( (ia==1 ? ref_pose.residue_type(ir-1).mainchain_atom( ref_pose.residue_type(ir-1).mainchain_atoms().size() ) : ref_pose.residue_type(ir).mainchain_atom(ia - 1) ), (ia==1 ? ir-1 : ir) ); //The previous atom is ia-1 in this residue, unless this atom is the first atom in this residue, in which case the previous atom is iamax in the previous residue.
			core::id::AtomID const nextatom( (ia==iamax ? ref_pose.residue_type(ir+1).mainchain_atom(1) : ref_pose.residue_type(ir).mainchain_atom(ia + 1) ), (ia==iamax ? ir+1 : ir) ); //The next atom is ia+1 in this residue, unless this atom is the last atom in this residue, in which case the next atom is the first atom in the next residue.
			pose.conformation().set_bond_angle( prevatom, thisatom, nextatom, numeric::angle_radians<core::Real>( ref_pose.xyz(prevatom), ref_pose.xyz(thisatom), ref_pose.xyz(nextatom) ) );
		}
	}

	//TODO properly handle the first and last residues using the extra residue xyz coordinates that were generated!

	pose.update_residue_neighbors();

	return;
}

/// @brief Copy backbone dihedral values from one pose, where helix mainchain atom coordinates have been
/// set with the Crick equations, to another with ideal geometry.
void copy_helix_dihedrals(
	core::pose::Pose &pose,
	core::pose::Pose const &ref_pose,
	core::Size const helix_start,
	core::Size const helix_end
) {
	for ( core::Size ir=helix_start; ir<=helix_end; ++ir ) {
		for ( core::Size itors=1, itorsmax=ref_pose.residue(ir).mainchain_torsions().size(); itors<=itorsmax; ++itors ) {
			pose.conformation().set_torsion( core::id::TorsionID( ir, core::id::BB, itors ), ref_pose.conformation().torsion( core::id::TorsionID( ir, core::id::BB, itors ) ) );
		}
	}

	//TODO properly handle the first and last residues using the extra residue xyz coordinates that were generated!

	pose.update_residue_neighbors();

	return;
}

/// @brief Align mainchain atoms of pose to ref_pose mainchain atoms.
///
void align_mainchain_atoms(
	core::pose::Pose &pose,
	core::pose::Pose const &ref_pose,
	core::Size const helix_start,
	core::Size const helix_end
) {
	core::id::AtomID_Map< core::id::AtomID > amap;
	core::pose::initialize_atomid_map(amap, pose, core::id::AtomID::BOGUS_ATOM_ID());

	for ( core::Size ir=helix_start; ir<=helix_end; ++ir ) {
		for ( core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia ) {
			amap[core::id::AtomID(pose.residue_type(ir).mainchain_atom(ia),ir)] = core::id::AtomID(ref_pose.residue_type(ir).mainchain_atom(ia),ir);
		}
	}

	core::scoring::superimpose_pose( pose, ref_pose, amap );

	return;
}

/// @brief Align mainchain atoms of pose to ref_pose mainchain atoms,
/// moving ONLY the residues involved in the alignment.
void align_mainchain_atoms_of_residue_range(
	core::pose::Pose &pose,
	core::pose::Pose const &ref_pose,
	core::Size const helix_start,
	core::Size const helix_end
) {
	core::pose::Pose pose_copy(pose);

	core::id::AtomID_Map< core::id::AtomID > amap;
	core::pose::initialize_atomid_map(amap, pose_copy, core::id::AtomID::BOGUS_ATOM_ID());

	for ( core::Size ir=helix_start; ir<=helix_end; ++ir ) {
		for ( core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia ) {
			amap[core::id::AtomID(pose.residue_type(ir).mainchain_atom(ia),ir)] = core::id::AtomID(ref_pose.residue_type(ir).mainchain_atom(ia),ir);
		}
	}

	//Align pose_copy to ref_pose
	core::scoring::superimpose_pose( pose_copy, ref_pose, amap );

	//Copy only those residues from pose_copy that were used in the alignment to pose:
	for ( core::Size ir=helix_start; ir<=helix_end; ++ir ) {
		pose.replace_residue( ir, pose_copy.residue(ir), false );
	}

	pose.update_residue_neighbors();

	return;
}

/// @brief Given a comma-separated list of residue names, separate these out into a vector of
/// residue names.
/// @details The string_in string is the input; the vect_out vector is the output (which will
/// be reset by this operation).
void parse_resnames (
	std::string const &string_in,
	utility::vector1< std::string > &vect_out
) {
	vect_out.clear();

	//Strip whitespace from the string:
	std::string nowhitespace(string_in);
	for ( signed long /*must be signed*/ i=0; i < static_cast< signed long >(nowhitespace.length()); ++i ) {
		if ( nowhitespace[i]=='\t' || nowhitespace[i]==' ' ) {
			nowhitespace.erase(i,1);
			--i;
		}
	}

	//Split the comma-separated list into the vector:
	std::stringstream ss(nowhitespace);
	while ( ss.good() ) {
		std::string aa_name;
		getline(ss, aa_name, ',');
		vect_out.push_back(aa_name);
	}

	return;
}

} //helical_bundle
} //namespace protocols

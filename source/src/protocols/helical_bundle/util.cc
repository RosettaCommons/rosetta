// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/util.cc
/// @brief  Utility functions for helical bundle construction.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// C++ headers
#include <stdio.h>

// Unit Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
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

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace helical_bundle {

	static thread_local basic::Tracer TR("protocols.helical_bundle.util");

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

		char linebuffer[256];
		std::string line_out;

		for(core::Size i=1, imax=r1.size(); i<=imax; ++i) {
			sprintf(linebuffer, "r1\t%.12f\n", r1[i]);
			line_out=linebuffer;
			outfile.write( line_out, line_out.length() );
		}

		sprintf(linebuffer, "omega1\t%.12f\nz1\t%.12f\n", omega1, z1);
		line_out=linebuffer;
		outfile.write(line_out, line_out.length() );

		for(core::Size i=1, imax=delta_omega1.size(); i<=imax; ++i) {
			sprintf(linebuffer, "delta_omega1\t%.12f\n", delta_omega1[i]);
			line_out=linebuffer;
			outfile.write( line_out, line_out.length() );
		}

		for(core::Size i=1, imax=delta_z1.size(); i<=imax; ++i) {
			sprintf(linebuffer, "delta_z1\t%.12f\n", delta_z1[i]);
			line_out=linebuffer;
			outfile.write( line_out, line_out.length() );
		}

		outfile.close();

		return;
	}

	void read_minor_helix_params (
		std::string const &filename,
		utility::vector1 < core::Real > &r1,
		core::Real &omega1,
		core::Real &z1,
		utility::vector1 < core::Real > &delta_omega1,
		utility::vector1 < core::Real > &delta_z1
	) {
		using namespace utility::io;

		r1.clear();
		delta_omega1.clear();
		delta_z1.clear();
		omega1 = 0.0;
		z1 = 0.0;

		std::string filename_formatted = filename;
		if(utility::file::file_extension(filename_formatted)!="crick_params") filename_formatted+= ".crick_params";

		izstream infile;
		infile.open( filename_formatted );
		if(!infile.good()) {
			filename_formatted = "protocol_data/crick_parameters/" + utility::file::file_basename(filename_formatted) + ".crick_params";
			basic::database::open( infile, filename_formatted );
			runtime_assert_string_msg( infile.good(), "In protocols::helical_bundle::read_minor_helix_params: Unable to open .crick_params file for read!" );
		}

		if(TR.Debug.visible()) TR.Debug << "Reading " << filename_formatted << std::endl;

		while(true) {
			std::string current_line = "";
			infile.getline(current_line);//Get the current line and store it in current_line.
			if(infile.eof()) break;

			if(TR.Debug.visible()) TR.Debug << current_line << std::endl;

			if(current_line[0] == '#') continue; //Ignore lines that start with a pound sign.

			std::stringstream ss(current_line);
			char buffer [25];
			ss.getline(buffer, 25, '\t');
			std::string strbuffer=buffer;
			//TR.Debug << strbuffer << "." << std::endl; //DELETE ME
			if(strbuffer=="r1") {
				ss.getline(buffer, 25);
				r1.push_back( atof(buffer) );
			} else if (strbuffer=="delta_omega1") {
				ss.getline(buffer, 25);
				delta_omega1.push_back( atof(buffer) );
			} else if (strbuffer=="delta_z1") {
				ss.getline(buffer, 25);
				delta_z1.push_back( atof(buffer) );
			} else if (strbuffer=="omega1") {
				ss.getline(buffer, 25);
				omega1=atof(buffer);
			} else if (strbuffer=="z1") {
				ss.getline(buffer, 25);
				z1=atof(buffer);
			} else continue;

		}

		infile.close();

		if(TR.Debug.visible()) {
			TR.Debug << "Finished reading " << filename_formatted << std::endl;
			for(core::Size i=1, imax=r1.size(); i<=imax; ++i) {
				TR.Debug << "r1[" << i << "]\t" << r1[i] << std::endl;
			}
			TR.Debug << "omega1\t" << omega1 << std::endl;
			TR.Debug << "z1\t" << z1 << std::endl;
			for(core::Size i=1, imax=delta_omega1.size(); i<=imax; ++i) {
				TR.Debug << "delta_omega1[" << i << "]\t" << delta_omega1[i] << std::endl;
			}
			for(core::Size i=1, imax=delta_z1.size(); i<=imax; ++i) {
				TR.Debug << "delta_z1[" << i << "]\t" << delta_z1[i] << std::endl;
			}
		}

		runtime_assert_string_msg( (r1.size() == delta_omega1.size()) && ( r1.size() == delta_z1.size() ), "In protocols::helical_bundle::read_minor_helix_params: the r1, delta_omega1, and delta_z1 lists are of different lengths.  Check the .crick_params file, since something is clearly wonky." );

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
		bool const invert_helix,
		utility::vector1 < core::Real > const &r1,
		core::Real const &omega1,
		core::Real const &z1,
		utility::vector1 < core::Real > const &delta_omega1,
		core::Real const &delta_omega1_all,
		utility::vector1 < core::Real > const &delta_z1,
		bool &failed
	) {

		outvector.clear();
		failed=false;

		core::Size helix_length=helix_end-helix_start+1;

		core::Real t = -1.0 * static_cast<core::Real>(helix_length + 2) / 2.0 + delta_t;

		for(core::Size ir2=helix_start-1; ir2<=helix_end+1; ++ir2) { //Loop through all residues in the helix, padded by one on either side
			core::Size ir = ir2;
			//Repeat start and end residues to pad by one:
			if(ir2==helix_start-1) ir=helix_start;
			if(ir2==helix_end+1) ir=helix_end;
			utility::vector1 < numeric::xyzVector <core::Real> > innervector;
			for(core::Size ia=1, iamax=helixpose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia) { //Loop through all mainchain atoms in the current helix residue
				innervector.push_back( numeric::crick_equations::XYZ_BUNDLE(
					t, r0, omega0,
					( invert_helix ? numeric::constants::d::pi - delta_omega0 : delta_omega0 ),
					r1[ia], omega1-omega0, z1, delta_omega1[ia]+delta_omega1_all, delta_z1[ia], failed )
				);
				if(invert_helix) {
					innervector[ia].x( -1.0*innervector[ia].x() );
					//innervector[ia].y( -1.0*innervector[ia].y() );
					innervector[ia].z( -1.0*innervector[ia].z() );
				}
				if(failed) break;
			}
			if(failed) break;
			outvector.push_back(innervector);
			t+=1.0;
		}

		if(failed) outvector.clear();

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

		for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
			for(core::Size ia=1, iamax=atom_positions[index].size(); ia<=iamax; ++ia) {
				pose.set_xyz( core::id::AtomID( ia, ir ), atom_positions[index][ia] );
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
		for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
			for(core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia) {
				if(ia==1 && ir==helix_start) continue; //Skip the first atom.
				core::id::AtomID const thisatom( ia, ir );
				core::id::AtomID const prevatom( (ia==1 ? iamax : ia - 1), (ia==1 ? ir-1 : ir) ); //The previous atom is ia-1 in this residue, unless this atom is the first atom, in which case the previous atom is iamax in the previous residue.
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
		for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
			for(core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia) {
				if(ia==1 && ir==helix_start) continue; //Skip the first atom.
				if(ia==iamax && ir==helix_end) continue; //Skip the last atom.
				core::id::AtomID const thisatom( ia, ir );
				core::id::AtomID const prevatom( (ia==1 ? iamax : ia - 1), (ia==1 ? ir-1 : ir) ); //The previous atom is ia-1 in this residue, unless this atom is the first atom in this residue, in which case the previous atom is iamax in the previous residue.
				core::id::AtomID const nextatom( (ia==iamax ? 1 : ia + 1), (ia==iamax ? ir+1 : ir) ); //The next atom is ia+1 in this residue, unless this atom is the last atom in this residue, in which case the next atom is the first atom in the next residue.
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
		for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
			for(core::Size itors=1, itorsmax=ref_pose.residue(ir).mainchain_torsions().size(); itors<=itorsmax; ++itors) {
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
			core::pose::initialize_atomid_map(amap, pose, core::id::BOGUS_ATOM_ID);

			for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
				for(core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia) {
					amap[core::id::AtomID(ia,ir)] = core::id::AtomID(ia,ir);
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
			core::pose::initialize_atomid_map(amap, pose_copy, core::id::BOGUS_ATOM_ID);

			for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
				for(core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia) {
					amap[core::id::AtomID(ia,ir)] = core::id::AtomID(ia,ir);
				}
			}

			//Align pose_copy to ref_pose
			core::scoring::superimpose_pose( pose_copy, ref_pose, amap );

			//Copy only those residues from pose_copy that were used in the alignment to pose:
			for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
				pose.replace_residue( ir, pose_copy.residue(ir), false );
			}

			pose.update_residue_neighbors();

			return;
	}

} //helical_bundle
} //namespace protocols

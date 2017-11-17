// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


#include <core/id/AtomID_Map.hh>
#include <devel/init.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <fstream>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <protocols/moves/Mover.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/string_util.hh>

#include <apps/pilot/will/will_util.hh>

static basic::Tracer TR( "coiled_coil_search" );


using core::pose::Pose;
using std::string;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using numeric::xyzVector;
using numeric::min;
using namespace ObjexxFCL::format;

typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef vector1<Size> Sizes;

void strip_termini(core::pose::Pose & pose) {
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue(i).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(pose,i);
		if ( pose.residue(i).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(pose,i);
	}
}

Vec
get_helix_axis(core::pose::Pose const & pose, Size start, Size end) {

	Vec summed_axes = Vec(0,0,0);

	for ( Size ir=start; ir < end; ir++ ) {
		for ( Size jr=ir+1; jr <= end; jr++ ) {
			summed_axes += pose.residue(jr).xyz(1) - pose.residue(ir).xyz(1);
		}
	}

	return summed_axes;

}

Vec
get_helix_center(core::pose::Pose const & pose, Size start, Size end) {

	Vec summed_xyzs = Vec(0,0,0);

	for ( Size ir=start; ir<=end; ir++ ) {
		summed_xyzs += pose.residue(ir).xyz(1);
	}

	summed_xyzs /= ((end - start) + 1);

	return summed_xyzs;

}

void
make_pdb_atom(numeric::xyzVector<Real> const & coords, std::ostream & out, Size & atom_id) {

	out << "ATOM  " << I(5,atom_id) << "  CA  ALA A" << I(4,atom_id) << "    " <<
		F(8,3,coords(1)) <<
		F(8,3,coords(2)) <<
		F(8,3,coords(3)) <<
		F(6,2,1.0) << F(6,2,1.0) << '\n';

	atom_id++;

}

int
main (int argc, char *argv[]){
	try{
		using namespace core;
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using namespace core::conformation::symmetry;

		devel::init( argc, argv );

		for ( Size ifile = 1; ifile <= option[in::file::s]().size(); ++ifile ) {
			string fname = option[in::file::s]()[ifile];
			pose::Pose pose;
			import_pose::pose_from_file(pose,fname, core::import_pose::PDB_file);

			// Get the starting and ending positions of the first and last helices
			core::scoring::dssp::Dssp dssp(pose);
			Size nstart=0,nstop=0,cstart=0,cstop=0;
			for ( Size i = 1; i <= pose.size(); ++i ) {
				if ( !nstart && dssp.get_dssp_secstruct(i) == 'H' ) {
					nstart = i;
				} else if (  nstart && dssp.get_dssp_secstruct(i) != 'H' ) {
					nstop = i-1;
					break;
				}
			}
			for ( Size i = pose.size(); i >= 1; --i ) {
				if ( !cstop && dssp.get_dssp_secstruct(i) == 'H' ) {
					cstop = i;
				} else if (  cstop && dssp.get_dssp_secstruct(i) != 'H' ) {
					cstart = i+1;
					break;
				}
			}
			//TR << fname << " " << nstart << " " << nstop << " " << pose.size()-cstart << " " << pose.size()-cstop << std::endl;

			Vec z = Vec(0,0,1);
			for ( Size nc = 0; nc <= 1; nc++ ) {
				Vec h_axis = Vec(0,0,0);
				Vec h_cent = Vec(0,0,0);
				if ( !nc && nstart < 7 && nstop-nstart > 6 ) {
					h_axis = get_helix_axis(pose, nstart, nstop);
					h_cent = get_helix_center(pose, nstart, nstop);
				} else if ( nc && pose.size()-cstop < 7 && cstop-cstart > 6 ) {
					h_axis = get_helix_axis(pose, cstart, cstop);
					h_cent = get_helix_center(pose, cstart, cstop);
				} else {
					continue;
				}
				h_axis.normalize();

				// Angle test
				Real ang = numeric::angle_degrees(h_axis,Vec(0,0,0),Vec(0,0,1));
				if ( fabs(ang-35.2643436) > 10.0 && fabs((180.0-ang)-35.2643436) > 10.0 ) { TR << "ang fail " << ang << std::endl; continue; }

				TR << "HIT: " << fname << " " << ang << std::endl;

				////////////////
				// For those proteins that have helices at the correct angle:
				////////////////

				// Create the symmetric pose
				core::pose::Pose symm_pose = pose;
				core::pose::symmetry::make_symmetric_pose(symm_pose);
				SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(symm_pose);
				std::map<Size,SymDof> dofs = symm_info->get_dofs();
				int sym_jump = 0;
				for ( std::map<Size,SymDof>::iterator i = dofs.begin(); i != dofs.end(); i++ ) {
					Size jump_num = i->first;
					if ( sym_jump == 0 ) {
						sym_jump = jump_num;
					} else {
						utility_exit_with_message("Can only handle one subunit!");
					}
				}
				if ( sym_jump == 0 ) {
					utility_exit_with_message("No jump defined!");
				}

				// Define which subs are part of the oligomeric building block. This will be used to define which
				// subunits need to be checked for helix pairing with chain A.
				Sizes intra_subs;
				if ( !option[in::olig_design::num_subs_building_block]() ) {
					utility_exit_with_message("ERROR: You have not set the required option -in::olig_design::num_subs_building_block");
				} else {
					for ( int intrasub=1; intrasub<=option[in::olig_design::num_subs_building_block](); intrasub++ ) {
						intra_subs.push_back(intrasub);
					}
				}

				// Get starting position
				Vec start_trans = symm_pose.jump(sym_jump).get_translation();
				Mat start_rot   = symm_pose.jump(sym_jump).get_rotation();
				Vec neg_disp = Vec(1.0,0,0);
				if ( ang >= 90.0 ) {
					neg_disp = Vec(-1.0,0,0);
				}

				for ( Size r=30; r<=60; r++ ) {

					// Move the pose to the new radius
					Vec trans = start_trans + r * neg_disp;
					core::kinematics::Jump j = symm_pose.jump(sym_jump);
					j.set_translation(trans);
					symm_pose.set_jump(sym_jump,j);

					// Recalculate the helical axis and position at the new radius.
					Vec h1_axis = Vec(0,0,0);
					Vec h1_cent = Vec(0,0,0);
					if ( !nc ) {
						h1_axis = get_helix_axis(symm_pose, nstart, nstop);
						h1_cent = get_helix_center(symm_pose, nstart, nstop);
					} else {
						h1_axis = get_helix_axis(symm_pose, cstart, cstop);
						h1_cent = get_helix_center(symm_pose, cstart, cstop);
					}
					h1_axis.normalize();

					// Figure out which subunit is next to chain A. The helices in each chain will be antiparallel by construction.
					Size nres_monomer = symm_info->num_independent_residues();
					for ( Size isub = intra_subs.size()+1; isub <= symm_info->subunits(); isub++ ) {
						Size resi_offset = (isub-1)*nres_monomer;
						Vec h2_axis = Vec(0,0,0);
						Vec h2_cent = Vec(0,0,0);
						if ( !nc ) {
							h2_axis = get_helix_axis(symm_pose, nstart + resi_offset, nstop + resi_offset);
							h2_cent = get_helix_center(symm_pose, nstart + resi_offset, nstop + resi_offset);
						} else {
							h2_axis = get_helix_axis(symm_pose, cstart + resi_offset, cstop + resi_offset);
							h2_cent = get_helix_center(symm_pose, cstart + resi_offset, cstop + resi_offset);
						}
						h2_axis.normalize();
						Vec neg_h2_axis = Vec(0,0,0) - h2_axis;
						// The helical vectors cannot be more than 20 degrees off, given the selection above.
						// There will be multiple solutions to this -- equal to the degree of symmetry of the
						// building block. Just take the first and rotate accordingly.
						if ( numeric::angle_degrees(neg_h2_axis, Vec(0,0,0), h1_axis) > 20.0 ) {
							//      TR << "ang fail: " << numeric::angle_degrees(neg_h2_axis, Vec(0,0,0), h1_axis) << " isub: " << isub << std::endl;
							continue;
						} // else, other chain of interest = isub

						// Business between chain A & other chain starts here
						// Sample the rotational degree of freedom
						Real w = 0.0;
						Real dw = 1.0;
						for ( Size iw=0; iw<360; iw++ ) {
							Mat rot = numeric::x_rotation_matrix_degrees(w) * start_rot;
							j.set_rotation(rot);
							symm_pose.set_jump(sym_jump,j);


							// Recalculate the helical axes and positions at the new radius.
							h1_axis = Vec(0,0,0);
							h1_cent = Vec(0,0,0);
							if ( !nc ) {
								h1_axis = get_helix_axis(symm_pose, nstart, nstop);
								h1_cent = get_helix_center(symm_pose, nstart, nstop);
							} else {
								h1_axis = get_helix_axis(symm_pose, cstart, cstop);
								h1_cent = get_helix_center(symm_pose, cstart, cstop);
							}
							h1_axis.normalize();
							h2_axis = Vec(0,0,0);
							h2_cent = Vec(0,0,0);
							if ( !nc ) {
								h2_axis = get_helix_axis(symm_pose, nstart + resi_offset, nstop + resi_offset);
								h2_cent = get_helix_center(symm_pose, nstart + resi_offset, nstop + resi_offset);
							} else {
								h2_axis = get_helix_axis(symm_pose, cstart + resi_offset, cstop + resi_offset);
								h2_cent = get_helix_center(symm_pose, cstart + resi_offset, cstop + resi_offset);
							}
							h2_axis.normalize();


							//      TR << "nstart = " << nstart << " nstop = " << nstop << " cstart = " << cstart << " cstop = " << cstop << " resi_offset = " << resi_offset << " isub = " << isub << std::endl;

							// Find the perpendicular distance between h1_axis and the two-fold
							// Taken from http://mathforum.org/library/drmath/view/54731.html
							Vec P = (h1_cent + h2_cent) / 2;
							Vec AB = h1_cent + h1_axis;
							Vec AP = P - h1_cent;
							Vec x = AB.cross(AP);
							Real a = x.length();
							Real d = a / AB.length();

							/*
							Size atom_id = 1;
							std::ostringstream out;
							make_pdb_atom(h1_cent, out, atom_id);
							make_pdb_atom(h1_cent+h1_axis, out, atom_id);
							make_pdb_atom(h2_cent, out, atom_id);
							make_pdb_atom(h2_cent+h2_axis, out, atom_id);
							make_pdb_atom(P, out, atom_id);
							make_pdb_atom(AB, out, atom_id);
							make_pdb_atom(AP, out, atom_id);
							make_pdb_atom(x, out, atom_id);
							TR << std::endl << out.str() << std::endl;
							*/

							//      TR << "r = " << r << " w = " << w << " d = " << d << std::endl;
							std::string fn = ObjexxFCL::string_of(r)+"_"+ObjexxFCL::string_of(w)+".pdb";
							//      symm_pose.dump_pdb(fn);
							if ( d>4.0 && d<5.71 ) {
								TR << "HIT: r = " << r << " w = " << w << " d = " << d << std::endl;
								//symm_pose.dump_pdb(fn);
								//utility_exit_with_message("Break for now.");
							}
							w = w + dw;
						} // for w

						break;

					} // for isub

					/*
					// Reset the pose to the starting position so that it can be moved in the next cycle
					// IS THIS EVEN NECESSARY, SINCE I AM SETTING THE JUMP EACH TIME INSTEAD OF INCREMENTING IT?
					j.set_translation(start_trans);
					j.set_rotation(start_rot);
					symm_pose.set_jump(sym_jump,j);
					*/

				} // for Size r

			} // for Size nc

		} // ifile
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;


}

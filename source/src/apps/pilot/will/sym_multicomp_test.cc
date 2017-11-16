// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/scoring/rms_util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <devel/init.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

using namespace core;
using namespace core::pose;
using namespace core::kinematics;
using namespace core::id;
using namespace utility;
using namespace numeric;


typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

static basic::Tracer TR( "sym_multicomp_test" );

bool check_coords_match(
	vector1<Vec> const & a,
	vector1<Vec> const & b
){
	if ( a.size() != b.size() ) utility_exit_with_message("ERROR");
	bool err = false;
	for ( int ir = 1; ir <= (int)a.size(); ++ir ) {
		if ( a[ir].distance(b[ir]) > 0.00001 ) {
			// std::cerr << "different at " << ir << std::endl;
			err = true;
		}
	}
	return !err;
}


vector1<Vec>
chain_coords(Pose const & pose, char chain, core::Size nres) {
	vector1<Vec> result;
	std::map<core::Size,char> conf2chain = core::pose::conf2pdb_chain(pose);
	for ( core::Size i = 1; i <= nres; ++i ) {
		if ( conf2chain[pose.chain(i)] == chain ) {
			result.push_back(pose.xyz(AtomID(2,i)));
		}
	}
	return result;
}
vector1<Vec>
non_chain_coords(Pose const & pose, char chain, core::Size nres) {
	vector1<Vec> result;
	std::map<core::Size,char> conf2chain = core::pose::conf2pdb_chain(pose);
	for ( core::Size i = 1; i <= nres; ++i ) {
		if ( conf2chain[pose.chain(i)] != chain ) {
			result.push_back(pose.xyz(AtomID(2,i)));
		}
	}
	return result;
}

void move_jump(core::pose::Pose & pose, int jnum) {
	Jump j = pose.jump(jnum);
	j.set_rotation( x_rotation_matrix_degrees(/*numeric::random::uniform()**/10.0) * j.get_rotation() );
	pose.set_jump(jnum,j);
}

int main (int argc, char *argv[]) {

	try {

		devel::init(argc,argv);
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		utility::vector1<std::string> files = option[in::file::s]();
		for ( int ifile = 1; ifile <= (int)files.size(); ++ifile ) {
			core::pose::Pose pose;
			core::import_pose::pose_from_file(pose,files[ifile], core::import_pose::PDB_file);
			core::Size nres = pose.size();
			Pose init(pose);
			core::pose::symmetry::make_symmetric_pose(pose);
			Pose test;
			core::pose::symmetry::extract_asymmetric_unit(pose,test);
			if ( core::scoring::CA_rmsd(init,test) > 0.001 ) {
				TR << "FAIL " << files[ifile] << " " << option[basic::options::OptionKeys::symmetry::symmetry_definition]() << std::endl;
				continue;
			}
			bool fail = false;
			core::conformation::symmetry::SymmetryInfoCOP syminfo = core::pose::symmetry::symmetry_info(pose);
			std::map<core::Size,core::conformation::symmetry::SymDof> dofs( syminfo->get_dofs() );
			for ( std::map<core::Size,core::conformation::symmetry::SymDof>::const_iterator i = dofs.begin(); i != dofs.end(); ++i ) {
				std::string jname = syminfo->get_jump_name(i->first);
				char chain = jname[jname.size()-1];
				vector1<Vec> preC   =     chain_coords(pose,chain,nres);
				vector1<Vec> preNoC = non_chain_coords(pose,chain,nres);
				move_jump(pose,i->first);
				vector1<Vec> postC   =     chain_coords(pose,chain,nres);
				vector1<Vec> postNoC = non_chain_coords(pose,chain,nres);
				if (  check_coords_match(preC  ,postC)   ) TR.Error << "FAIL: chain " << chain << " not moved" << std::endl;
				if ( !check_coords_match(preNoC,postNoC) ) TR.Error << "FAIL: not chain " << chain << " moved" << std::endl;
				// The following expression result was unused; I have attempted to correct it. ~Labonte
				fail = !check_coords_match(preC,postC) || !check_coords_match(preNoC,postNoC);
			}
			if ( !fail ) TR << "WOOT " << files[ifile] << " " << option[basic::options::OptionKeys::symmetry::symmetry_definition]() << std::endl;
		}
		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}



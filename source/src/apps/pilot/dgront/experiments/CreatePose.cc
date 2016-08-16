// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <core/types.hh>
#include <devel/init.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>


#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <utility/io/izstream.hh>

#include <core/sequence/util.hh>

#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/annotated_sequence.hh>

OPT_1GRP_KEY( File, out, pdbfname )
OPT_1GRP_KEY( File, in, phi_psi_omega )

void register_options() {

	OPT(in::file::s);
	OPT(in::file::native);
	OPT(in::file::fasta);
	OPT(in::file::residue_type_set);
//	OPT(out::pdb);
	NEW_OPT( out::pdbfname, "provides a file name for the output PDB file","out_pose.pdb" );
	NEW_OPT( in::phi_psi_omega, "provides phi,psi,omega angles to create a poly-alanine pose","");
}

int main(int argc, char * argv[]) {
    try {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using std::string;

	register_options();
	devel::init(argc, argv);

	core::pose::Pose extended_pose;
	std::string out_file_name = (option[out::pdbfname]());

//        rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

        core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
	                option[ in::file::residue_type_set ]()
	);
	string sequence;
	if (option[in::file::fasta].user()) {
		sequence = core::sequence::read_fasta_file_return_str(
				option[in::file::fasta]()[1]);

		core::pose::make_pose_from_sequence(extended_pose, sequence, rsd_set);

		// make extended chain
		for (Size pos = 1; pos <= extended_pose.total_residue(); pos++) {
			if (!extended_pose.residue(pos).is_protein())
				continue;
			extended_pose.set_phi(pos, -150);
			extended_pose.set_psi(pos, 150);
			extended_pose.set_omega(pos, 180);
		}
		extended_pose.dump_pdb(out_file_name);
		return 0;
	}

	if (option[in::file::s].user()) {

		core::pose::PoseOP tmp_pose(new core::pose::Pose);
		std::string fn = option[in::file::s](1);
		core::import_pose::pose_from_file(*tmp_pose, fn, core::import_pose::PDB_file);

		sequence = tmp_pose->sequence();
		core::pose::make_pose_from_sequence(extended_pose, sequence,
				*rsd_set);
		for (Size i = 1; i <= sequence.size(); ++i) {
			extended_pose.set_phi(i, tmp_pose->phi(i));
			extended_pose.set_psi(i, tmp_pose->psi(i));
			extended_pose.set_omega(i, tmp_pose->omega(i));
		}
		extended_pose.dump_pdb(out_file_name);
		return 0;
	}

	if (option[in::phi_psi_omega].user()) {

	    utility::io::izstream data(option[in::phi_psi_omega]());

	    std::string line;
	    Real phi,psi,omega;
	    utility::vector1<Real> all_phi;
	    utility::vector1<Real> all_psi;
	    utility::vector1<Real> all_omega;
	    while (!data.fail()) {
		char c = data.peek();
		if (c == '#' || c == '\n') {
			getline(data, line); //comment
			continue;
		}
		data >> phi >> psi >> omega;
		all_phi.push_back(phi);
		all_psi.push_back(psi);
		all_omega.push_back(omega);
	    }
	    std::string seq = "";
	    for(Size i=1;i<=all_phi.size();i++)
		seq += "A";
	    core::pose::make_pose_from_sequence(extended_pose, seq,*rsd_set);

	    // make a chain
	    for (Size pos = 1; pos <= extended_pose.total_residue(); pos++) {
		extended_pose.set_phi(pos, all_phi[pos]);
		extended_pose.set_psi(pos, all_psi[pos]);
		extended_pose.set_omega(pos, all_omega[pos]);
		std::cerr<<all_phi[pos]<<" "<<all_psi[pos]<<" "<<all_omega[pos]<<"\n";
	    }
	    extended_pose.dump_pdb(out_file_name);
	    return 0;
	}

	char aa[] = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
			'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };
	for (int i = 0; i < 20; i++) {
		sequence = "";
		for (int j = 0; j < 10; j++)
			sequence += aa[i];
		std::cerr << sequence << std::endl;
		core::pose::make_pose_from_sequence(extended_pose, sequence,
				*(chemical::ChemicalManager::get_instance()->residue_type_set(
						"fa_standard")));
		for (Size pos = 1; pos <= extended_pose.total_residue(); pos++) {
			if (!extended_pose.residue(pos).is_protein())
				continue;
			extended_pose.set_phi(pos, -150);
			extended_pose.set_psi(pos, 150);
			extended_pose.set_omega(pos, 180);
		}
		std::string s = out_file_name + "-";
		s += aa[i];
		s += ".pdb";
		std::cerr << s << std::endl;
		extended_pose.dump_pdb(s);
	}
    } catch ( utility::excn::EXCN_Base const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
    }
	return 0;
}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file mm_params
/// @brief Test code for reading molecular mechanics bond angle atom types/parameters
/// @details
/// This app will dump all bond angles in all amino acid types. It will then construct a MMBondAngleResidueTypeParam
/// object assuming the use of residue type theta0 angles and dump that as well. This command only requres a single
/// command line flag, -in:path:database. There are two optional flags, -use_residue_type_theta0 and
/// -central_atoms_to_score, whose documentation can be viewed by specifying -help.
/// This app will also sampling all branching atom combinations for N, CA, and C, and write them to the database.


// Core Headers
#include <devel/init.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Atom.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParam.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <core/scoring/ScoringManager.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/Pose.hh>

#include <protocols/branch_angle/BranchAngleOptimizer.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

#include <core/chemical/MMAtomType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


static thread_local basic::Tracer TR( "mm_params" );

OPT_1GRP_KEY(Boolean, mm_params, use_residue_type_theta0)
OPT_1GRP_KEY(StringVector, mm_params, central_atoms_to_score)

int
main
(
	int argc,
	char * argv []
)
{
    try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	OPT(in::path::database);
	NEW_OPT(mm_params::use_residue_type_theta0, "use ResidueType for theta0", false);
	NEW_OPT(mm_params::central_atoms_to_score, "specify central atoms to score", utility::vector1<std::string>());

	devel::init(argc, argv);

	TR << "Initialization Successful" << std::endl;

	core::scoring::mm::MMBondAngleLibrary const & mm_bondangle_library(
		core::scoring::ScoringManager::get_instance()->get_MMBondAngleLibrary()
	);

	//mm_bondangle_library->pretty_print();

	core::chemical::ResidueTypeSetCOP residue_set(
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
	);

	//for (core::chemical::ResidueTypeSet::const_residue_iterator residue_iter(residue_set->all_residues_begin());
	//	residue_iter != residue_set->all_residues_end(); ++residue_iter) {
	//
	//	TR << residue_iter->first << "\n";
	//	TR << " name3: " << residue_iter->second->name3() << " natoms: " << residue_iter->second->natoms() << "\n";
	//}

	core::scoring::mm::MMBondAngleResidueTypeParamSet residue_type_param_set;

	residue_type_param_set.mm_bondangle_library(&mm_bondangle_library);
	residue_type_param_set.use_residue_type_theta0(option[ mm_params::use_residue_type_theta0 ]);
	residue_type_param_set.central_atoms_to_score(option[ mm_params::central_atoms_to_score ]);

	for (std::list< core::chemical::AA >::const_iterator aa_iter(residue_set->aas_defined_begin());
		aa_iter != residue_set->aas_defined_end(); ++aa_iter) {

		TR << *aa_iter << std::endl;

		core::chemical::ResidueTypeCOPs const & aa_caps(residue_set->aa_map(*aa_iter));

		for (core::chemical::ResidueTypeCOPs::const_iterator residue_iter(aa_caps.begin());
		     residue_iter != aa_caps.end(); ++residue_iter) {

			core::chemical::ResidueType const & residue_type(**residue_iter);

			if (!residue_type.is_protein()) continue;
			if (residue_type.name().find("acetylated") != std::string::npos) continue;
			if (residue_type.name().find("carboxylated") != std::string::npos) continue;
			if (residue_type.name().find("monomethylated") != std::string::npos) continue;
			if (residue_type.name().find("dimethylated") != std::string::npos) continue;
			if (residue_type.name().find("trimethylated") != std::string::npos) continue;
			if (residue_type.name().find("hydroxylated") != std::string::npos) continue;
			if (residue_type.name().find("phosphorylated") != std::string::npos) continue;
			if (residue_type.name().find("sulfated") != std::string::npos) continue;
			if (residue_type.name().find("diiodinated") != std::string::npos) continue;
			if (residue_type.name().find("methylamidated") != std::string::npos) continue;

			TR << residue_type.name() << std::endl;
			TR << " name3: " << residue_type.name3() << " natoms: " << residue_type.natoms() << std::endl;
			TR << " is_lower_terminus: " << residue_type.is_lower_terminus() << std::endl;
			TR << " is_upper_terminus: " << residue_type.is_upper_terminus() << std::endl;

			TR << "Intraresidue Bond Angles:" << std::endl;
			for (core::Size i = 1; i <= residue_type.num_bondangles(); ++i) {
				core::chemical::bondangle_atom_set const & atom_set(residue_type.bondangle(i));
				core::Real residue_type_theta0(numeric::angle_radians(residue_type.atom(atom_set.key1()).ideal_xyz(), residue_type.atom(atom_set.key2()).ideal_xyz(), residue_type.atom(atom_set.key3()).ideal_xyz()));

				TR << residue_type.atom_name(atom_set.key1()) << "-"
				   << residue_type.atom_name(atom_set.key2()) << "-"
					 << residue_type.atom_name(atom_set.key3());

				std::string const & type1(residue_type.mm_atom_type(atom_set.key1()).name());
				std::string const & type2(residue_type.mm_atom_type(atom_set.key2()).name());
				std::string const & type3(residue_type.mm_atom_type(atom_set.key3()).name());

				TR << " (" << type1 << "-" << type2 << "-" << type3 << ")";

				core::scoring::mm::mm_bondangle_library_citer_pair mm_pair(mm_bondangle_library.lookup(type1, type2, type3));
				core::Real mm_theta0(0);
				core::Real mm_Ktheta(0);
				bool mm_set(false);

				for ( core::scoring::mm::mm_bondangle_library_citer i = mm_pair.first, e = mm_pair.second; i != e; ++i ) {

					if (mm_set) TR << " (Multiple Parameters Defined)";
					mm_Ktheta = (i->second).key1();
					mm_theta0 = (i->second).key2();
					mm_set = true;
				}

				TR << " mm_Ktheta: " << mm_Ktheta << " mm_theta0: " << numeric::conversions::degrees(mm_theta0)
				   << " rt_theta0: " << numeric::conversions::degrees(residue_type_theta0) << std::endl;
			}

			for (core::Size i = 1; i <= residue_type.n_residue_connections(); ++i) {
				core::chemical::ResidueConnection const & residue_connection(residue_type.residue_connection(i));
				core::Size const connection_atomno(residue_connection.atomno());
				core::Vector external_xyz(residue_connection.icoor().build(residue_type));
				core::chemical::AtomIndices const & bonded_neighbors(residue_type.bonded_neighbor(connection_atomno));

				TR << "Connection " << i << " Bond Angles:" << std::endl;
				for (core::Size j = 1; j <= bonded_neighbors.size(); ++j) {

					core::Real residue_type_theta0(numeric::angle_radians(residue_type.atom(bonded_neighbors[j]).ideal_xyz(), residue_type.atom(connection_atomno).ideal_xyz(), external_xyz));

					TR << residue_type.atom_name(bonded_neighbors[j]) << "-"
					   << residue_type.atom_name(connection_atomno) << "-?";

					std::string const & type1(residue_type.mm_atom_type(bonded_neighbors[j]).name());
					std::string const & type2(residue_type.mm_atom_type(connection_atomno).name());

					TR << " (" << type1 << "-" << type2 << "-?)";
					TR << " rt_theta0: " << numeric::conversions::degrees(residue_type_theta0) << std::endl;
				}
			}

			TR << residue_type_param_set.get(residue_type);
		}
	}

	protocols::branch_angle::BranchAngleOptimizer branchopt;
	//branchopt.read_database();

	std::string all_aas("ACDEFGHIKLMNPQRSTVWY");

	std::string protseq("XXX");

	for (core::Size rttheta0 = 0; rttheta0 <= 1; ++rttheta0) {
		if (rttheta0) {
			std::cout << "ResidueType theta0" << std::endl;
			core::scoring::mm::MMBondAngleResidueTypeParamSetOP param_set( new core::scoring::mm::MMBondAngleResidueTypeParamSet() );
			param_set->use_residue_type_theta0(true);
			branchopt.bond_angle_residue_type_param_set(param_set);
		} else {
			std::cout << "MM theta0" << std::endl;
		}
		for (core::Size aa1 = 1; aa1 <= all_aas.length(); ++aa1) {
			protseq[0] = all_aas[aa1-1];
			std::cout << "Optimizing " << protseq[0] << "XX" << std::endl;
			for (core::Size aa2 = 1; aa2 <= all_aas.length(); ++aa2) {
				protseq[1] = all_aas[aa2-1];
				for (core::Size aa3 = 1; aa3 <= all_aas.length(); ++aa3) {
					protseq[2] = all_aas[aa3-1];
					for (core::Size revft = 0; revft <= 1; ++revft) {
						core::pose::PoseOP pose_op( new core::pose::Pose );
						core::pose::Pose & pose = *pose_op;
						core::pose::make_pose_from_sequence(pose, protseq, "fa_standard");
						if (revft) {
							core::kinematics::FoldTree foldtree;
							foldtree.add_edge(3, 1, core::kinematics::Edge::PEPTIDE);
							pose.fold_tree(foldtree);
						}
						using core::id::AtomID;
						if (protseq[0] != 'P' && revft) branchopt.optimize_angles(pose, pose.atom_tree().atom(AtomID(1, 1)).child(0)->id(), AtomID(1, 1), AtomID(2, 1));
						branchopt.optimize_angles(pose, AtomID(1, 1), AtomID(2, 1), AtomID(3, 1));
						branchopt.optimize_angles(pose, AtomID(2, 1), AtomID(3, 1), AtomID(1, 2));
						if (protseq[1] != 'P') branchopt.optimize_angles(pose, AtomID(3, 1), AtomID(1, 2), AtomID(2, 2));
						branchopt.optimize_angles(pose, AtomID(1, 2), AtomID(2, 2), AtomID(3, 2));
						branchopt.optimize_angles(pose, AtomID(2, 2), AtomID(3, 2), AtomID(1, 3));
						if (protseq[2] != 'P') branchopt.optimize_angles(pose, AtomID(3, 2), AtomID(1, 3), AtomID(2, 3));
						branchopt.optimize_angles(pose, AtomID(1, 3), AtomID(2, 3), AtomID(3, 3));
						if (!revft) branchopt.optimize_angles(pose, AtomID(2, 3), AtomID(3, 3), pose.atom_tree().atom(AtomID(3, 3)).child(0)->id());
					}
				}
			}
		}
	}

	//branchopt.write_database();
	branchopt.write_undefined_coef1("branch_angle_1_undefined.txt");
	branchopt.write_undefined_coef2("branch_angle_2_undefined.txt");

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }

	return 0;
}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/ward_design.cc
/// @brief starting point for xtal contact design with Andrew Ward


#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
//Auto Headers


#include <apps/pilot/will/will_util.ihh>

using numeric::conversions::radians;

static basic::Tracer TR( "ward_design" );

using core::Size;
using core::Real;
using core::pose::Pose;
using core::id::AtomID;
typedef numeric::xyzVector<Real> Vec;
typedef utility::vector1<Vec>    Vecs;
typedef numeric::xyzMatrix<Real> Mat;
using protocols::moves::MoverOP;
using core::scoring::ScoreFunctionOP;
using numeric::random::uniform;
using utility::vector1;

Size INTERFACE_CONTACT_DIS = 9.0; // dist between CBs to call and interface contact


vector1<Size> interface_residues(core::pose::Pose const & pose, vector1<Size> primary_subs) {
	vector1<Size> iface_res;
	core::conformation::symmetry::SymmetryInfoCOP syminfo = core::pose::symmetry::symmetry_info(pose);
	Size nres = syminfo->num_total_residues_without_pseudo();

	for ( Size i = 1; i <= nres ; ++i ) {
		Size isub = syminfo->subunit_index(i);
		if ( std::find(primary_subs.begin(),primary_subs.end(),isub) == primary_subs.end() ) continue; // skip if i not primary
		bool iface_contact_i = false;
		std::string contact_atom_i = "CB"; if ( !pose.residue(i).has(contact_atom_i) ) contact_atom_i = "CA";
		for ( Size j = 1; j <= nres; ++j ) {
			Size jsub = syminfo->subunit_index(j);
			if ( std::find(primary_subs.begin(),primary_subs.end(),jsub) != primary_subs.end() ) continue; // skip if j *is* primary
			// at this point, i is in primary subs and j is not
			std::string contact_atom_j = "CB"; if ( !pose.residue(j).has(contact_atom_j) ) contact_atom_j = "CA";
			if ( pose.residue(i).xyz(contact_atom_i).distance( pose.residue(j).xyz(contact_atom_j) ) <= INTERFACE_CONTACT_DIS ) {
				iface_contact_i = true;
				break;
			}
		}
		if ( iface_contact_i ) {
			iface_res.push_back(i);
		}
	}

	return iface_res;
}

void minimize(core::pose::Pose & pose, ScoreFunctionOP sf, vector1<Size> iface_res ) {
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(false);
	movemap->set_bb(false);
	movemap->set_chi(true);
	// uncomment this to minimize BB in interface res
	// for(vector1<Size>::iterator i = iface_res.begin(); i != iface_res.end(); ++i) movemap->set_bb(*i,true);
	protocols::minimization_packing::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true );
	m.apply(pose);
}


void design(Pose & pose, ScoreFunctionOP sf, vector1<Size> iface_res ) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	vector1< bool > aas(20,true);
	aas[core::chemical::aa_cys] = false;
	aas[core::chemical::aa_pro] = false;
	aas[core::chemical::aa_gly] = false;

	sf->score(pose);
	Real worig = sf->get_weight(core::scoring::res_type_constraint);
	if ( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);
	utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);

	// loop over residues
	for ( Size i = 1; i <= pose.size(); ++i ) {

		// decide what to do with res i
		bool design_this_res = false;
		bool repack_this_res = false;
		if ( pose.residue(i).name3()=="GLY" ||
				pose.residue(i).name3()=="PRO" ||
				pose.residue(i).name3()=="VRT" ||
				pose.residue(i).name3()=="CYS" ||
				pose.residue(i).name3()=="CYD"
				) {
			bool design_this_res = false;
			bool repack_this_res = false;
		} else if ( std::find(iface_res.begin(),iface_res.end(),i) != iface_res.end() ) { // is interface res
			design_this_res = true;
			repack_this_res = false;
		} else {
			design_this_res = false;
			repack_this_res = true;
		}

		if ( design_this_res ) {
			bool tmp = aas[pose.residue(i).aa()];
			aas[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
			aas[pose.residue(i).aa()] = tmp;
			task->nonconst_residue_task(i).or_include_current(true);
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else if ( repack_this_res ) {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_include_current(true);
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else {
			task->nonconst_residue_task(i).prevent_repacking();
		}

	}

	// do the design
	core::pack::make_symmetric_PackerTask_by_truncation(pose,task);
	protocols::minimization_packing::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);


	// clean up
	pose.remove_constraints( res_cst );
	sf->set_weight(core::scoring::res_type_constraint,worig);
}


/////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {


		using namespace core;
		using namespace pose;
		using namespace protocols;
		using namespace moves;
		using namespace ObjexxFCL::format;
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using numeric::random::uniform;
		using ObjexxFCL::lead_zero_string_of;

		devel::init(argc,argv);

		std::string outdir = option[out::file::o]() + "/";

		core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();

		utility::vector1<Size> primary_subs;
		primary_subs.push_back(1);
		primary_subs.push_back(2);

		for ( Size ifile = 1; ifile <= option[in::file::s]().size(); ifile++ ) {

			std::string  in_fname = option[in::file::s]()[ifile];
			std::string out_fname = utility::file_basename(in_fname);

			Pose init;
			import_pose::pose_from_file(init,in_fname, core::import_pose::PDB_file);

			Pose init_sym = init;
			core::pose::symmetry::make_symmetric_pose(init_sym);

			init_sym.dump_pdb( outdir + out_fname + "_init_sym.pdb" );
			vector1<Size> iface_res = interface_residues(init_sym,primary_subs);

			TR << "interface residues " << std::endl;
			for ( vector1<Size>::iterator i = iface_res.begin(); i != iface_res.end(); ++i ) TR << *i << "+";
			TR << std::endl;


			Pose pose = init_sym;

			design(pose,sf,iface_res);
			pose.dump_pdb( outdir + out_fname + "_design_test.pdb" );

			// minimize(pose,sf,iface_res);
			// pose.dump_pdb( outdir + out_fname + "_minimize_test.pdb" );

		}


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}



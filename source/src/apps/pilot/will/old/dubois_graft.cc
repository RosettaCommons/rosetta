// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/dubois_graft.cc
/// @brief graft peptide structs onto dubois catalyst

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <devel/init.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
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
#include <basic/Tracer.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>


#include <apps/pilot/will/will_util.ihh>

using numeric::conversions::radians;

static basic::Tracer TR( "dubois_graft" );

using core::Size;
using core::Real;
using core::id::AtomID;
typedef numeric::xyzVector<Real> Vec;
typedef utility::vector1<Vec>    Vecs;
typedef numeric::xyzMatrix<Real> Mat;
using protocols::moves::MoverOP;
using core::scoring::ScoreFunctionOP;
using numeric::random::uniform;


class DuboisChiMover : public protocols::moves::Mover {
	Real mag12_, mag3_, mag4_;
public:
	DuboisChiMover(Real mag12 = 5.0, Real mag3=60.0, Real mag4=20.0) : mag12_(mag12),mag3_(mag3),mag4_(mag4) {}
	std::string get_name() const { return "DuboisChiMover"; }
	void apply( core::pose::Pose & pose ) {
		// play some games to sample rotation of aro group but keep chi 1/2 pos close
		Real a = mag12_*numeric::random::gaussian();
		Real b = mag12_*numeric::random::gaussian();
		Real n = mag3_ *numeric::random::gaussian();
		Real d = mag4_ *numeric::random::gaussian();
		pose.set_chi(1,1, pose.chi(1,1) + b     );
		pose.set_chi(2,1, pose.chi(2,1) + a - n );
		pose.set_chi(3,1, pose.chi(3,1) + n     );
		pose.set_chi(4,1, pose.chi(4,1) + d     );
		// pose.set_psi(2,pose.psi(2)+10.0*numeric::random::gaussian());
	}
};


void minimize(core::pose::Pose & pose, ScoreFunctionOP sf) {
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(false);
	movemap->set_chi(false);
	movemap->set_bb(false);
	movemap->set_bb(1,true);
	movemap->set_chi(1,true);
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true );
	m.apply(pose);
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

		// init pose from dubois catalyst ncaa
		Pose init;
		import_pose::pose_from_file(init,"input/pdb/NPH_0001.pdb", core::import_pose::PDB_file);

		core::chemical::ResidueTypeSetCAP fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		core::chemical::ResidueType const & alatype( fa_residue_set->name_map("ALA") );
		core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue(alatype);
		core::pose::remove_upper_terminus_type_from_pose_residue(init,1);

		ScoreFunctionOP sf = core::scoring::get_score_function();

		Pose init_sym = init;
		core::pose::symmetry::make_symmetric_pose(init_sym);
		minimize(init_sym,sf);
		Real init_sym_sc = sf->score(init_sym);

		for ( Size ifile = 1; ifile <= option[in::file::s]().size(); ++ifile ) {
			// add pept res
			utility::vector1<Pose> pepts;
			import_pose::pose_from_file(pepts,*fa_residue_set,option[in::file::s]()[ifile], core::import_pose::PDB_file);

			for ( Size imodel = 1; imodel <= pepts.size(); ++imodel ) {
				Pose pept = pepts[imodel];
				core::pose::remove_lower_terminus_type_from_pose_residue(pept,1);
				Real pept_sc = 4.0*sf->score(pept);
				// prepend ala to use for alignment
				pept.prepend_polymer_residue_before_seqpos( *ala, 1, true );
				core::kinematics::Stub s = getxform(pept.residue(1),init.residue(1));
				xform_pose(pept,s);
				// pept.dump_pdb("pept.pdb");

				Pose init_pept = init;
				for ( Size i = 2; i <= pept.size(); ++i ) {
					init_pept.append_polymer_residue_after_seqpos( pept.residue(i), i-1, false );
				}
				init_pept.set_xyz( AtomID(4,1), pept.residue(1).xyz("O") );
				add_lower_terminus_type_to_pose_residue(init_pept,1);
				add_variant_type_to_pose_residue(init_pept,"VIRTUAL_NTERM",1);
				// init_pept.dump_pdb("test0.pdb");

				TR << "init score " << init_sym_sc << " " << pept_sc << " " << pept_sc+init_sym_sc << std::endl;

				core::pose::symmetry::make_symmetric_pose(init_pept);

				std::string fname = utility::file_basename( option[in::file::s]()[ifile] );
				core::io::silent::SilentFileData sfd;
				Pose pose = init_pept;
				Pose best = init_pept;
				Real bestsc = 9e9;
				Size ang_num_step = 36;
				Real ang_per_step = 360.0 / ang_num_step;
				Size ntot=0,nmin=0;
				for ( Size i = 1; i <= ang_num_step; ++i ) {
					TR << "ifile imodel chi2 " << ifile << " " << imodel << " " << i << std::endl;
					Real chi2 = i*ang_per_step;
					for ( Size j = 1; j <= ang_num_step; ++j ) {
						Real chi1 = j*ang_per_step;
						for ( Size k = 1; k <= ang_num_step; ++k ) {
							Real psi = k*ang_per_step;
							pose.set_chi( 2, 1, chi2 );
							pose.set_chi( 1, 1, chi1 );
							pose.set_psi( 1, psi );
							std::string tag = fname+"__"+lead_zero_string_of(imodel,2)+"__"+lead_zero_string_of(ang_per_step*i,3)+"_"+lead_zero_string_of(ang_per_step*j,3)+"_"+lead_zero_string_of(ang_per_step*k,3);
							ntot++;
							sf->score(pose);
							// pose.dump_pdb(tag+".pdb");
							// {
							//  core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
							//  ss_out->fill_struct(pose,tag);
							//            ss_out->add_energy( "scoredelta", sf->score(pose) - init_sym_sc - pept_sc );
							//  sfd.write_silent_struct( *ss_out, fname+"__"+lead_zero_string_of(imodel,2)+".sc" );
							// }
							if ( sf->score(pose) > init_sym_sc+pept_sc+100.0 ) continue;
							nmin++;
							minimize(pose,sf);
							if ( sf->score(pose) < bestsc ) {
								bestsc = sf->score(pose);
								best = pose;
							}
							{
								core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
								ss_out->add_energy( "scoredelta", sf->score(pose) - init_sym_sc - pept_sc );
								ss_out->fill_struct(pose,tag+"_MIN");
								sfd.write_silent_struct( *ss_out, option[out::file::o]+"/"+fname+"__"+lead_zero_string_of(imodel,2)+".sc" );
							}

						}
					}
				}
				if ( sf->score(best) > init_sym_sc+pept_sc+100.0 ) continue;
				TR << "dumping best " << fname << " " << imodel << " ntot " << ntot << " nmin " << nmin << std::endl;
				best.dump_pdb( option[out::file::o]+"/"+fname+"__"+lead_zero_string_of(imodel,2)+"__best.pdb");
			}
		}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

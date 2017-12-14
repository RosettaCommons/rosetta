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
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/types.hh>
#include <fstream>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>
#include <sstream>
#include <utility/io/izstream.hh>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/database/open.hh>


static basic::Tracer TR( "RosettaHoles" );


using std::string;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using numeric::xyzVector;
using namespace core::scoring::packing;
using namespace ObjexxFCL::format;

// HolesParamsRes params("/Users/sheffler/project/holes/holes_params.txt");

std::string tag_from_pose( core::pose::Pose & pose ) {
	std::string tag( "empty_tag" );
	if ( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		tag = static_cast< basic::datacache::CacheableString const & >
			( pose.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	}
	return tag;
}


class RosettaHolesMover : public protocols::moves::Mover {

private:
	core::scoring::packing::HolesParams hp_resl_,hp_dec_,hp_dec15_;
	bool make_pdb_,make_voids_,atom_scores_,residue_scores_;
	core::scoring::ScoreFunctionOP sf_;

public:

	RosettaHolesMover() {
		hp_resl_ .read_data_file(basic::database::full_name("scoring/rosettaholes/resl.params"));
		hp_dec_  .read_data_file(basic::database::full_name("scoring/rosettaholes/decoy25.params"));
		hp_dec15_.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
		make_pdb_       = basic::options::option[ basic::options::OptionKeys::holes::make_pdb ]();
		make_voids_     = basic::options::option[ basic::options::OptionKeys::holes::make_voids ]();
		atom_scores_    = basic::options::option[ basic::options::OptionKeys::holes::atom_scores ]();
		residue_scores_ = basic::options::option[ basic::options::OptionKeys::holes::residue_scores ]();
		sf_             = core::scoring::get_score_function();
		TR  << "HEADER         File/Tag         RosHol    resl   dec25   dec15   AAresl   AAdec25   AAdec15  RosScore" << std::endl;
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) TR << "     RMS";

	}
	void       make_pdb( bool make_pdb       )       { make_pdb_       = make_pdb ;       }
	bool       make_pdb(                     ) const { return            make_pdb_;       }
	void    atom_scores( bool atom_scores    )       { atom_scores_    = atom_scores ;    }
	bool    atom_scores(                     ) const { return            atom_scores_;    }
	void residue_scores( bool residue_scores )       { residue_scores_ = residue_scores ; }
	bool residue_scores(                     ) const { return            residue_scores_; }

	void
	apply(
		core::pose::Pose & pose
	) override {
		using namespace std;
		using namespace core;
		using namespace io::pdb;
		using namespace pose;
		using namespace scoring;
		using namespace packing;
		using namespace basic::options;
		using namespace utility;
		using core::Size;

		TR << "apply to: " << tag_from_pose(pose) << std::endl;

		Size MAX_RES = 5000;
		if ( pose.size() > MAX_RES ) {
			TR << "nres > " << MAX_RES << ", skipping pose " << tag_from_pose(pose) << std::endl;
			return;
		}

		HolesResult result = compute_rosettaholes_score(pose,hp_resl_,hp_dec_,hp_dec15_);

		TR << "got holes result: " << tag_from_pose(pose) << std::endl;


		Real rms = -1;
		if ( basic::options::option[ OptionKeys::in::file::native ].user() ) {
			Pose native;
			core::import_pose::pose_from_file( native, basic::options::option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
			rms = scoring::CA_rmsd( native, pose );
		}

		TR  << RJ(30,tag_from_pose(pose)) << " "
			<< F(7,4,result.score) << " "
			<< F(7,4,result.resl_score) << " "
			<< F(7,4,result.decoy_score) << " "
			<< F(7,4,result.dec15_score);
		TR  << F(9,2,result.resl_score *result.natom) << " "
			<< F(9,2,result.decoy_score*result.natom) << " "
			<< F(9,2,result.dec15_score*result.natom);
		TR << " " << F(9,3,(*sf_)(pose));
		if ( rms != -1 ) TR << " " << F(7,4,rms);
		TR << std::endl;

		if ( residue_scores_ ) {
			// HolesResult dec15result = compute_dec15_score(pose);
			HolesResult dec15result = compute_holes_score(pose,hp_dec15_);
			for ( Size i = 1; i <= pose.size(); ++i ) {
				Real rscore = 0.0;
				for ( Size j = 5; j <= pose.residue(i).nheavyatoms(); ++j ) {
					rscore += dec15result.atom_scores[AtomID(j,i)];
				}
				rscore /= max(1,(int)pose.residue(i).nheavyatoms()-4);
				TR << "residue_score " << I(5,i) << " " << pose.residue(i).name3() << " " << F(10,3,rscore) << std::endl;
			}
		}

		if ( make_pdb_ || make_voids_ ) {
			std::string dir = "";
			if ( basic::options::option[basic::options::OptionKeys::out::file::o].user() ) {
				dir = basic::options::option[basic::options::OptionKeys::out::file::o]() + "/";
			}
			std::string fname = dir  + file_basename(tag_from_pose(pose))+"_cavs.pdb";
			ofstream out(fname.c_str());
			if ( !out.good() ) {
				utility_exit_with_message("can't open file: "+fname);
			}

			// ONE OR THE OTHER OF THESE BLOCKS BELOW!
			TR << "dumping pdb to " << fname << std::endl;

			//JAB -XRW 2016
			core::pose::set_bfactors_from_atom_id_map(pose, result.atom_scores);
			//pose.dump_pdb(out, "NO_MODEL_LINE_IN_OUTPUT");
			pose.dump_pdb(out);

			//dump_bfactor_pdb(pose,result.atom_scores,out,"NO_MODEL_LINE_IN_OUTPUT");

			if ( make_voids_ ) {
				TR << "dumping cavities " << std::endl;
				core::scoring::packstat::output_packstat_pdb( pose, out );
			}
		}

	}

	std::string get_name() const override { return "RosettaHolesMover"; }
}; // end class RosettaHolesMover


int
main (int argc, char *argv[])
{

	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace utility;

		devel::init( argc, argv );

		RosettaHolesMover holes;

		protocols::jobdist::not_universal_main( holes );

		return 0;


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>  //make_pose_from_sequence
#include <core/pose/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/rms_util.tmpl.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <iostream>
#include <string>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>

#include <ObjexxFCL/format.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>

// loop hash closability check
//#include <apps/pilot/rayyrw/closability_measurement.hh>
#include <apps/pilot/rayyrw/overlap_score.hh>
#include <apps/pilot/rayyrw/rms_util.hh>
#include <apps/pilot/rayyrw/util.hh>

static basic::Tracer TR("test_score_methods");

using namespace ObjexxFCL::format;

// options declaration
class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};


//////////////////////////////////////////////////////////////////////////////

OPT_KEY( FileVector, ref_fragfiles )
OPT_KEY( FileVector, fragfiles )
OPT_1GRP_KEY( File, out, score )

// overlap score function options
// if dist < 2:
//     -1*overlap_wt*exp( -1*dist*steepness_wt )
// else:
//     return clash_return_score
OPT_KEY( Real, overlap_dist_cutoff )
//OPT_KEY( Real, overlap_wt )
OPT_KEY( Real, steepness_wt )
OPT_KEY( Integer, seq_sep )
OPT_KEY( Real, clash_dist )
OPT_KEY( Real, clash_return_score )


void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( ref_fragfiles, "[vector] fragment pdb files", "" );
	NEW_OPT( fragfiles, "[vector] fragment pdb files", "" );
	NEW_OPT( out::score, "score_file", "score.sc" );
	NEW_OPT( overlap_dist_cutoff, "cutoff of giving a penalty for a pair", 2.0 );
	//NEW_OPT( overlap_wt, "weight for overlap score (useless)", 1.0 );
	NEW_OPT( steepness_wt, "weight for steepness for the exp function", 1.0 );
	NEW_OPT( seq_sep, "for how many residues separate you would start considering clashes", 5 );
	NEW_OPT( clash_dist, "to what distance you would consider a clash", 3.0 );
	NEW_OPT( clash_return_score, "how much of a penalty per clash", 10 );
}



int main( int argc, char* argv[] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		ThisApplication::register_options();

		devel::init(argc, argv);

		// get pdb names into a vector
		utility::vector1<utility::file::FileName> ref_fragpdbs,
			fragpdbs;

		ref_fragpdbs = option[ ref_fragfiles ]().vector();
		fragpdbs = option[ fragfiles ]().vector();

		// for the for loop speed issue, don't need to walk through vector all the time
		core::Size ref_fragpdbs_size = ref_fragpdbs.size();
		core::Size fragpdbs_size = fragpdbs.size();

		// read poses and positions information from pdbs all at once
		utility::vector1<core::pose::Pose> ref_fragpdb_poses( ref_fragpdbs_size ),
			fragpdb_poses( fragpdbs_size );

		utility::vector1<core::Size> ref_fragpdb_positions( ref_fragpdbs_size ),
			fragpdb_positions( fragpdbs_size );

		utility::vector1<core::Real> ref_fragpdb_rmsds( ref_fragpdbs_size, 0.0 ),
			fragpdb_rmsds( fragpdbs_size, 0.0 );


		// output: poses and positions
		read_pdbs( ref_fragpdbs, ref_fragpdb_poses, ref_fragpdb_positions );
		read_pdbs( fragpdbs, fragpdb_poses, fragpdb_positions );

		if ( option[ in::file::native ].user() ) {
			// native pose
			core::pose::Pose native_pose;
			core::chemical::ResidueTypeSetCOP cen_rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
			core::import_pose::pose_from_pdb( native_pose, *cen_rsd_set, basic::options::option[ basic::options::OptionKeys::in::file::native ]().name() );

			get_frag_rmsd( ref_fragpdb_rmsds, ref_fragpdb_positions, ref_fragpdb_poses, native_pose );
			get_frag_rmsd( fragpdb_rmsds, fragpdb_positions, fragpdb_poses, native_pose );
		}


		// output stream
		std::stringstream outlines;

		for ( core::Size i=1; i<=ref_fragpdbs_size; ++i ) {
			core::Size ref_fragpdb_pos = ref_fragpdb_positions[i];
			core::pose::Pose const &pose1 = ref_fragpdb_poses[i];

			for ( core::Size j=1; j<=fragpdbs_size; ++j ) {
				core::Size fragpdb_pos = fragpdb_positions[j];
				int pos_offset = fragpdb_pos - ref_fragpdb_pos;
				core::pose::Pose const &pose2 = fragpdb_poses[j];


				core::Real score = soft_overlap_score( pose1, ref_fragpdb_pos,
					pose2, fragpdb_pos,
					option[ overlap_dist_cutoff ],
					option[ steepness_wt ],
					option[ seq_sep ],
					option[ clash_dist ],
					option[ clash_return_score ]  );


				outlines << ref_fragpdbs[i].base()
					<< " "
					<< fragpdbs[j].base()
					<< RJ(6, pos_offset)
					<< F(10, 4, score )
					<< F(10, 4, ref_fragpdb_rmsds[i]+fragpdb_rmsds[j] )
					<< std::endl;
			} //ref_frag
		} //frag

		// write to the outfile
		utility::io::ozstream outfn( option[ out::score ] );
		outfn << "#fragA " << "fragB " << "offset " << "overlap_score " << "total_rmsd " << std::endl;
		outfn << outlines.str();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

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

#include <apps/pilot/rayyrw/closability_score.hh>
#include <apps/pilot/rayyrw/clash_score.hh>
#include <apps/pilot/rayyrw/rms_util.hh>
#include <apps/pilot/rayyrw/util.hh>

static THREAD_LOCAL basic::Tracer TR("cal_nonoverlap_score");

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
OPT_KEY( Real, clash_dist_threshold )
OPT_KEY( Integer, radial_search_radius )
OPT_KEY( Boolean, bound )
OPT_KEY( Boolean, bound_stringent )
OPT_KEY( Boolean, ignore_junction_for_gap1 )
OPT_1GRP_KEY( File, out, score )


void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( ref_fragfiles, "[vector] fragment pdb files", "" );
	NEW_OPT( fragfiles, "[vector] fragment pdb files", "" );
	NEW_OPT( clash_dist_threshold, "dist for defining clash", 3.0 );
	NEW_OPT( radial_search_radius, "radial search radius for loophash", 4 );
	NEW_OPT( out::score, "score_file", "score.sc" );
	NEW_OPT( bound, "instead of using loophash - for gap_size within 9, using bound csts", false );
	NEW_OPT( bound_stringent, "instead of using bound score, just return the value", false );
	NEW_OPT( ignore_junction_for_gap1, "when gap_size==1, don't check the terminus residues", false );
}



int main( int argc, char* argv[] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		ThisApplication::register_options();

		devel::init(argc, argv);

		// closab_score
		ClosabilityScore closab_score;

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
			core::import_pose::pose_from_file( native_pose, *cen_rsd_set, basic::options::option[ basic::options::OptionKeys::in::file::native ]().name() , core::import_pose::PDB_file);

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

				core::Real closab(0.0);
				core::Real clash(0.0);
				core::Size counts( 0 );
				int gap_size(0);

				// closab score
				if ( option[ bound ].user() ) {
					// use bound constraints to measure closability
					closab = closab_score.closability_score( pose1, ref_fragpdb_pos,
						pose2, fragpdb_pos,
						gap_size,
						option[ bound_stringent ] );
				} else {
					// use loophash to detect closability
					closab = closab_score.closability_score( pose1, ref_fragpdb_pos,
						pose2, fragpdb_pos,
						counts,
						gap_size,
						option[ radial_search_radius ] );
				}

				// clash score
				if ( option[ ignore_junction_for_gap1 ].user() ) {
					// calling soften_clash_score, which ingores junction residue check
					clash = clash_score( pose1, ref_fragpdb_pos,
						pose2, fragpdb_pos, option[ clash_dist_threshold ] );
				} else {
					clash = clash_score( pose1, pose2, option[ clash_dist_threshold ] );
				}


				outlines << ref_fragpdbs[i].base()
					<< " "
					<< fragpdbs[j].base()
					<< RJ(6, pos_offset )
					<< RJ(6, gap_size )
					<< RJ(6, counts )
					<< F(10, 4, closab )
					<< F(10, 4, clash )
					<< F(10, 4, ref_fragpdb_rmsds[i]+fragpdb_rmsds[j] )
					<< std::endl;
			} //ref_frag
		} //frag

		// write to the outfile
		std::string output_fn = option[ out::score ];
		if ( utility::file::file_exists( output_fn ) ) {
			utility_exit_with_message( output_fn + "has existed!\n" );
		}
		utility::io::ozstream outfn( output_fn );
		outfn << "#fragA " << "fragB " << "offset " << "gap_size " << "counts " << "closab_score " << "clash_score " << "total_rmsd " << std::endl;
		outfn << outlines.str();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/cluster/cluster.hh>
#include <protocols/loops/Loops.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <iostream>
#include <string>
#include <deque>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <ObjexxFCL/format.hh>

using ObjexxFCL::format::F;
using namespace protocols::cluster;


OPT_1GRP_KEY( FileVector, ccd, rmsd_matrix )

void
register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT(ccd::rmsd_matrix,"files to compare","");
}

bool
use_in_rmsd(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size resno,
	core::Size atomno
) {
	if(resno==pose1.n_residue() && atomno==pose1.residue(resno).atom_index( "O")) return false;
	if(resno==pose1.n_residue() && atomno==pose1.residue(resno).atom_index("SG")) return true;
	if(resno==1                 && atomno==pose1.residue(resno).atom_index("SG")) return true;
	if(pose1.residue(resno).has( "N") && pose2.residue(resno).has( "N") && pose1.residue(resno).atom_index( "N")==atomno) return true;
	if(pose1.residue(resno).has("CA") && pose2.residue(resno).has("CA") && pose1.residue(resno).atom_index("CA")==atomno) return true;
	if(pose1.residue(resno).has( "C") && pose2.residue(resno).has( "C") && pose1.residue(resno).atom_index( "C")==atomno) return true;
	if(pose1.residue(resno).has( "O") && pose2.residue(resno).has( "O") && pose1.residue(resno).atom_index( "O")==atomno) return true;
	if(pose1.residue(resno).has("CB") && pose2.residue(resno).has("CB") && pose1.residue(resno).atom_index("CB")==atomno) return true;
	return false;
}


class ClusterCycDisulf: public ClusterPhilStyle {
public:
	ClusterCycDisulf(){}

	virtual ~ClusterCycDisulf() {}
	virtual std::string get_name() const {
		return "ClusterCycDisulf";
	}
	protocols::moves::MoverOP clone() const {
		return new ClusterCycDisulf( *this );
	}

	virtual core::Real
	get_distance_measure(
		const core::pose::Pose & pose1,
		const core::pose::Pose & pose2
	) const {
		return core::scoring::rmsd_with_super( pose1, pose2, use_in_rmsd );
	}


	virtual void
	align_clusters() {
		Pose posecen = poselist[ clusterlist[0][0] ];
		for ( int i=0;i<(int)clusterlist.size();i++ ) {
			for ( int j=0;j<(int)clusterlist[i].size();j++ ) {
				Pose & pose = poselist[ clusterlist[i][j] ];
				core::id::AtomID_Map< core::id::AtomID > amap;
				core::pose::initialize_atomid_map(amap,pose,core::id::BOGUS_ATOM_ID);
				for(int ir = 1; ir <= (int)pose.n_residue(); ++ir) {
					for(int ia = 1; ia <= (int)pose.residue(ir).nheavyatoms(); ++ia) {
						if(use_in_rmsd(posecen,pose,ir,ia)) {
							using core::id::AtomID;
							amap[AtomID(ia,ir)] = AtomID(ia,ir);
						}
					}
				}
				core::scoring::superimpose_pose( pose, posecen, amap );
			}
		}
	}
};
typedef utility::pointer::owning_ptr<ClusterCycDisulf> ClusterCycDisulfOP;

using namespace core;
using namespace ObjexxFCL;
using namespace core::pose;
using namespace protocols;
using namespace basic::options;

int
main( int argc, char * argv [] ) {

	try {


	using namespace protocols;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;
	using namespace protocols::cluster;
	using namespace basic::options::OptionKeys::cluster;

	register_options();
	devel::init(argc, argv);


// def bestalign(s1,s2):
// 	ol1 = [x for x in cmd.get_object_list() if x.startswith(s1)]
// 	ol2 = [x for x in cmd.get_object_list() if x.startswith(s2)]
// 	keepers = []
// 	for o1 in ol1:
// 		mn = 9e9
// 		mo = None
// 		for o2 in ol2:
// 			r = cmd.align(o1,o2)[0]
// 			if r < mn:
// 				mn = r
// 				mo = o2
// 		cmd.align(o1,mo)
// 		keepers.append(mo)
// 		print "%7.3f    %10s %10s"%(mn,o1,mo)
// 	for o2 in ol2:
// 		if o2 not in keepers:
// 			cmd.delete(o2)


	if(option[ccd::rmsd_matrix].user()) {
		utility::vector1<std::string> files1 = option[in::file::s]();
		utility::vector1<std::string> files2 = option[ccd::rmsd_matrix]();
		utility::vector1< core::pose::PoseOP > poses1 = core::import_pose::poseOPs_from_files(files1, core::import_pose::PDB_file);
		utility::vector1< core::pose::PoseOP > poses2 = core::import_pose::poseOPs_from_files(files2, core::import_pose::PDB_file);
		for(int i = 1; i <= (int)poses1.size(); ++i) {
			core::Real mn = 9e9;
			int mj = 0;
			for(int j = 1; j <= (int)poses2.size(); ++j) {
				core::Real r = core::scoring::rmsd_with_super( *poses1[i], *poses2[j], use_in_rmsd );
				if( r < mn ) {
					mn = r;
					mj = j;
				}
			}
			std::cout << F(7,3,mn) << " " << files1[i] << " " << files2[mj] << std::endl;
			core::pose::Pose tmp1(*poses1[i]);
			core::id::AtomID_Map<core::id::AtomID> atom_map;
			core::pose::initialize_atomid_map(atom_map,tmp1,core::id::BOGUS_ATOM_ID);
			core::scoring::superimpose_pose(tmp1,*poses2[mj],atom_map);
		}


		return 0;
	}


	if ( !option[ out::output ].user() ) {
		option[ out::nooutput ].value( true );
	}

	int time_start = time(NULL);

	MoverOP mover;
	core::scoring::ScoreFunctionOP sfxn;
	sfxn = core::scoring::get_score_function();
	if ( option[ basic::options::OptionKeys::symmetry::symmetric_rmsd ]() ) {
		core::scoring::ScoreFunctionOP sfxn_sym =
						new core::scoring::symmetry::SymmetricScoreFunction( *sfxn );
		sfxn = sfxn_sym;
	}

	ClusterCycDisulfOP clustering = new ClusterCycDisulf;

	clustering->set_score_function( sfxn );
	clustering->set_cluster_radius(
		option[ basic::options::OptionKeys::cluster::radius ]()
	);
	clustering->set_population_weight(
		option[ basic::options::OptionKeys::cluster::population_weight ]()
	);


	core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->
		residue_type_set( option[ in::file::fullatom ]() ? "fa_standard" : "centroid");

	// Cluster the first up-to-400 structures by calculating a full rms matrix
	mover = clustering;

	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
	core::Size count = 0;
	while( input.has_another_pose() && (count < 1000 ) ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );
		mover->apply( pose );
		count ++;
	}

	clustering->do_clustering( option[ OptionKeys::cluster::max_total_cluster ]() );
	clustering->do_redistribution();

	std::cout << "Assigning extra structures ... " << std::endl;
	AssignToClustersMoverOP mover_add_structures = new AssignToClustersMover( clustering );
	mover_add_structures->set_score_function( sfxn );

	while( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );
		mover_add_structures->apply( pose );
	}

	clustering->sort_each_group_by_energy();
	if ( option[ sort_groups_by_energy        ].user() ) clustering->sort_groups_by_energy();
	if ( option[ remove_singletons            ].user() ) clustering->remove_singletons();
	if ( option[ limit_cluster_size           ].user() ) clustering->limit_groupsize( option[ limit_cluster_size ] );
	if ( option[ limit_clusters               ].user() ) clustering->limit_groups( option[ limit_clusters ] );
	if ( option[ limit_total_structures       ].user() ) clustering->limit_total_structures( option[ limit_total_structures] );
	if ( option[ remove_highest_energy_member ].user() ) clustering->remove_highest_energy_member_of_each_group();

	if ( option[ export_only_low ].user() ){
		clustering->sort_each_group_by_energy();
		clustering->sort_groups_by_energy( );
		clustering->export_only_low( option[ export_only_low ]() );
	}

	clustering->print_summary();
	clustering->align_clusters();
	clustering->print_cluster_PDBs( option[ out::prefix ]() );

	clustering->print_cluster_assignment();
	std::vector < Cluster >  const & clusterlist=clustering->get_cluster_list();
	std::list < int > sorted_list;
	for ( int i=0; i<(int)clusterlist.size(); i++ ) {
		for ( int j=0; j<(int)clusterlist[i].size(); j++ ) {
			sorted_list.push_back(  clusterlist[i][j]  );
		}
	}

	sorted_list.sort();

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


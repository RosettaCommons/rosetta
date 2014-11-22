/// @file
/// @brief

#include <apps/pilot/frank/spacegroup.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pose/CrystInfo.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/types.hh>
#include <devel/init.hh>

// james' magic
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <fstream>
#include <iostream>
#include <math.h>
#include <cmath>

#include <sstream>
#include <string>
#include <queue>
#include <cstdarg>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::pack;
using namespace core::pack::task;
using namespace core::pack::task::operation;
using namespace protocols::moves;
using namespace scoring;

OPT_1GRP_KEY(String, tors, mode)
OPT_1GRP_KEY(Real, tors, dev)

static thread_local basic::Tracer TR( "torsion.corrections" );

struct FragInfo {
	FragInfo ( ) : weight_(0.0), center_(0) {}

	FragInfo ( core::pose::PoseOP pose, core::Real weight, core::Size center ) :
		pose_(pose), weight_(weight), center_(center) {}

	core::pose::PoseOP pose_;
	core::Real weight_;
	core::Size center_;
};

core::Real distance( FragInfo const &f1, FragInfo const &f2) {
	if (f1.pose_->total_residue() != f2.pose_->total_residue()) return 999.0;
	if (f1.pose_->sequence() != f2.pose_->sequence()) return 999.0;
	if (f1.center_ != f2.center_) return 999.0;

	return core::scoring::CA_rmsd ( *f1.pose_, *f2.pose_ );
}

core::Size getbin(core::Real theta) {
	theta = fmod( theta, 360.0 );
	if (theta>180.0) theta -=360.0; // -180 to 180
	core::Size binnum = std::floor( (theta+180.0)/10.0 + 1.5 );
	if (binnum == 37) binnum = 1;
	return binnum;
}

//
void
mutate_to_ala(core::pose::Pose &pose, int center) {
	using namespace core;
	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	conformation::Residue const replace_res( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA"), true );

	for (int i=1; i<=(int)pose.total_residue(); ++i ) {
		if (i==center) continue;
		chemical::ResidueType const * cur_restype = & pose.residue_type( i );
		if ( cur_restype->aa() == chemical::aa_pro || cur_restype->aa() == chemical::aa_gly ) continue;

		utility::vector1< std::string > current_variants;
		pose.replace_residue( i, replace_res, true);
	}
}


//
void
make_fragments() {
  core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();

	std::map < core::chemical::AA, ObjexxFCL::FArray2D< utility::vector1<FragInfo> > > database;

	for (int i=1; i<=20; ++i) {
		database[ (core::chemical::AA)i ].dimension(36,36);
	}

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

  while( input.has_another_pose() ) {
    core::pose::Pose pose;
    input.fill_pose( pose );

		// dssp
		core::scoring::dssp::Dssp dssp(pose);
		dssp.insert_ss_into_pose(pose);

		// foreach residue
		for (int i=3; i<=((int)pose.total_residue())-2; ++i) {
			// ensure we have a full fragment, no chainbreaks
			char ss = pose.secstruct(i);
	    core::pose::PoseOP frag( new core::pose::Pose);
			core::Real phi, psi;
			core::Size center;
			bool have_connected_frag=true;

			if (ss == 'L') {
				for (int j=i-2; j<i+2; ++j) {
					if ( pose.fold_tree().is_cutpoint(j) ) {
						have_connected_frag = false;
					}
				}
				if (!have_connected_frag) continue;
				for (int j=i-2; j<=i+2; ++j) frag->append_residue_by_bond( pose.residue(j) );
				center = 3;
			}

			if (ss == 'H') {
				if (i-4<=0 || i-4 >= ((int)pose.total_residue())) continue;
				for (int j=i-4; j<i+4; ++j) {
					if ( pose.fold_tree().is_cutpoint(j) ) {
						have_connected_frag = false;
					}
				}
				if (!have_connected_frag) continue;
				for (int j=i-4; j<=i+4; ++j) frag->append_residue_by_bond( pose.residue(j) );

				center = 5;
			}

			if (ss == 'E') {
				for (int j=i-2; j<i+2; ++j) {
					if ( pose.fold_tree().is_cutpoint(j) ) {
						have_connected_frag = false;
					}
				}
				if (!have_connected_frag) continue;
				for (int j=i-2; j<=i+2; ++j) frag->append_residue_by_bond( pose.residue(j) );

				// find bb hbonds
				(*scorefxn)(pose); // to get interaction graph
				utility::vector1<bool> incl(pose.total_residue(), false);
				EnergyGraph const & energy_graph( pose.energies().energy_graph() );
				for (int j=i-2; j<=i+2; ++j) {
					core::conformation::Residue const &rsd_j( pose.residue(j) );

					for ( graph::Graph::EdgeListConstIter
							iru  = energy_graph.get_node(i)->const_edge_list_begin(),
							irue = energy_graph.get_node(i)->const_edge_list_end();
							iru != irue; ++iru ) {
						EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
						int k=(int)edge->get_other_ind(j);

						core::conformation::Residue const &rsd_k( pose.residue(k) );

						if (k>=i-2 && k<=i+2) continue;

						// check hbonding (O--H <= 3Å)
						core::Real dist1 = 999.0;
						if (rsd_k.aa() != core::chemical::aa_pro) dist1 = (rsd_j.atom("O").xyz() - rsd_k.atom("H").xyz()).length_squared();
						core::Real dist2 = 999.0;
						if (rsd_j.aa() != core::chemical::aa_pro) dist2 = (rsd_j.atom("H").xyz() - rsd_k.atom("O").xyz()).length_squared();
						if (dist1<=9 || dist2<=9) incl[k]=true;
					}
				}

				// zap singletons
				for (int j=2; j<=(int)incl.size()-1; ++j) incl[j] = incl[j] || (incl[j-1]&&incl[j+1]);
				for (int j=1; j<=(int)incl.size(); ++j) incl[j] = incl[j] && ( (j>1&&incl[j-1]) || (j<(int)incl.size()&&incl[j+1]));

				for (int j=1; j<=(int)incl.size(); ++j) {
					if (incl[j]) {
						if (j==1 || !incl[j-1]) frag->append_residue_by_jump( pose.residue(j), 1 );
						else frag->append_residue_by_bond( pose.residue(j) );
					}
				}

				center = 3;
			}

			mutate_to_ala( *frag, center);
			phi = frag->phi(center);
			psi = frag->psi(center);

			core::Size phi_bin = getbin(phi);
			core::Size psi_bin = getbin(psi);

			// check for duplicates
			utility::vector1<FragInfo> & fragstore = database[ pose.residue(i).aa() ](phi_bin, psi_bin);
			FragInfo new_fraginfo(frag,1.0,center);

			bool add_to_db=true;
			for (int j=1; j<=(int)fragstore.size(); ++j) {
				if ( distance(new_fraginfo, fragstore[j]) < option[ tors::dev ]() ) {
					fragstore[j].weight_ += 1.0;
					add_to_db = false;
				}
			}

			if (add_to_db) { fragstore.push_back( new_fraginfo ); }

			//frag.dump_pdb( "frag"+utility::to_string(i)+".pdb" );
		}
	}

	for (int i=1; i<=36; ++i)
	for (int j=1; j<=36; ++j) {
		std::cout << i << " " << j << " ";
		for (int a=1; a<=20; ++a) {
			utility::vector1<FragInfo> & fragstore = database[ (core::chemical::AA)a ](i, j);
			core::Size nfrags = fragstore.size();
			core::Real weightsum=0.0;
			for (int x=1; x<=(int)nfrags; ++x) weightsum+= fragstore[x].weight_;
			std::cout << nfrags << "/" << weightsum << " ";
		}
		std::cout << std::endl;
	}
}


int
main( int argc, char * argv [] ) {
try {
	NEW_OPT(tors::mode, "mode", "makefrags");
	NEW_OPT(tors::dev, "dev", 1.0);

	devel::init( argc, argv );

	option[ in::missing_density_to_jump ].value(true); // force this

	std::string mode = option[ tors::mode ]();

	if (mode == "makefrags") {
		make_fragments();
	} else if (mode == "calcscores") {
		;
	}

} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
return 0;

}


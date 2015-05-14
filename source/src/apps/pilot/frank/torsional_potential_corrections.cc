/// @file
/// @brief

#include <apps/pilot/frank/spacegroup.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/optimizeH.hh>
#include <core/pose/CrystInfo.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
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
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

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
#include <basic/datacache/CacheableString.hh>

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

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
OPT_1GRP_KEY(Real, tors, cluster_radius)
OPT_1GRP_KEY(String, tors, fragfile)
OPT_1GRP_KEY(String, tors, scorefile)
OPT_1GRP_KEY(Boolean, tors, rotmin)
OPT_1GRP_KEY(Real, tors, smoothing)
OPT_1GRP_KEY(Real, tors, cap)
OPT_1GRP_KEY(Real, tors, scale)

static thread_local basic::Tracer TR( "torsion.corrections" );

struct FragInfo {
	FragInfo ( ) : weight_(0), center_(0) {}

	FragInfo ( core::pose::PoseOP pose, core::Real weight, core::Size center, std::string tag ) :
		pose_(pose), weight_(weight), center_(center), tag_(tag) {}

	core::pose::PoseOP pose_;
	core::Size weight_;
	core::Size center_;
	std::string tag_;
};

struct ScoreFragInfo {
	ScoreFragInfo() : phibin_(0),psibin_(0),rotidx_(0), weight_(0), fragscore_(0) {}

	ScoreFragInfo(
			std::string id,	int aa, core::Size phibin, core::Size psibin,
			core::Size rotidx, core::Size r1, core::Size r2, core::Size r3, core::Size r4,
			core::Size weight, core::Real fragscore
	) {
		id_=id;
		aa_=aa;
		phibin_=phibin;
		psibin_=psibin;
		rotidx_=rotidx;
		rotid_[0]=r1; rotid_[1]=r2; rotid_[2]=r3; rotid_[3]=r4;
		weight_=weight;
		fragscore_=fragscore;
	}

	std::string id_;
	int aa_;
	core::Size phibin_,psibin_,rotidx_;
	core::Size rotid_[4];
	core::Size weight_;
	core::Real fragscore_;
};

core::Size
getRotID( core::Size r1, core::Size r2=0, core::Size r3=0, core::Size r4=0) {
	core::Size retval = r1;
	if (r2>0) retval += 3*(r2-1);
	if (r3>0) retval += 9*(r3-1);
	if (r4>0) retval += 27*(r4-1);
	return retval;
}

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


bool
is_semi_rot( core::chemical::AA aa ) {
	using namespace core::chemical;
	return (
			aa==aa_asn || aa==aa_asp || aa==aa_gln || aa==aa_glu ||
		  aa==aa_his || aa==aa_phe || aa==aa_trp || aa==aa_tyr
	);
}


void
mutate_to_ala(core::pose::Pose &pose, int center) {
	using namespace core;
	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	conformation::Residue const replace_res( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA"), true );

	// remove disulfides
	for (int i=1; i<=(int)pose.total_residue(); ++i ) {
		if ( pose.residue(i).has_variant_type( chemical::DISULFIDE ) ) {
  		TR << "Reverting disulfide to thiol type at resid " << i << std::endl;
			conformation::change_cys_state( i, "CYS", pose.conformation() );
		}
	}

	for (int i=1; i<=(int)pose.total_residue(); ++i ) {
		if (i==center) continue;
		chemical::ResidueType const * cur_restype = & pose.residue_type( i );
		if ( cur_restype->aa() == chemical::aa_pro || cur_restype->aa() == chemical::aa_gly ) continue;
		pose.replace_residue( i, replace_res, true);
	}
}


void
dump_table( ObjexxFCL::FArray2D<core::Real> table, std::string filename, std::string tag ) {
	std::ofstream outtable( filename.c_str(), std::ofstream::out | std::ofstream::app );

	outtable << tag << " = [ ...\n";
	for (int i=1; i<=table.u1(); ++i) {
		outtable << "   " << table(i,1);
		for (int j=2; j<=table.u2(); ++j) {
			outtable << "," << table(i,j);
		}
		if (i!=table.u1()) outtable << "; ...\n";
		else outtable << "];\n\n";
	}
}


void
correct_dunbrack() {
	using namespace core::chemical;

	// smooth things by assuming this many isoenergetic rotamers in each bin
	//core::Real MEST = 5.0;

	//params
	core::Real cap = option[ tors::cap ]();
	core::Real smoothing = option[ tors::smoothing ]();
	core::Real scale = 1.0 / option[ tors::scale ]();

	/////
	// 1 read scorefile, average in each residue bin and apply cap
	std::ifstream inscore( option[ tors::scorefile ]().c_str(), std::ifstream::in);

	utility::vector1< ScoreFragInfo > currentfrag;
	std::string currtag;
	core::Real bestscore=1e30;
	std::string semitag, currsemitag;

	// store all the tables in memory
	// map (aa,rotidx) -> probs
	std::map< std::pair< int, int > , ObjexxFCL::FArray2D<core::Real> > rotEsum;

	// for semirots, map (aa,rotidx) -> probs over all semirots
	std::map< std::pair< int, int > , utility::vector1< ObjexxFCL::FArray2D< core::Real > > > semirotEsum;

	// map aa -> count
	std::map< int , ObjexxFCL::FArray2D<core::Size> > rotcount;

	ScoreFragInfo fragline;
	inscore >> semitag;
	inscore >> fragline.id_  >> fragline.aa_ >> fragline.rotidx_;
	inscore >> fragline.rotid_[0] >> fragline.rotid_[1] >> fragline.rotid_[2] >> fragline.rotid_[3];
	inscore >> fragline.phibin_ >> fragline.psibin_ >> fragline.weight_ >> fragline.fragscore_;

	while (inscore.good()) {
		if (currtag.length() == 0) { // first fragment
			currtag = fragline.id_;
			currsemitag = semitag;
		}

		if (currtag != fragline.id_ || currsemitag != semitag) { // new fragment (or switch Rot->Semi), aggregate results of previous

			// enforce cap
			core::Real workingcap = bestscore+cap;
			for (int i=1; i<=(int)currentfrag.size(); ++i) {
				//currentfrag[i].fragscore_ -= bestscore;
				currentfrag[i].fragscore_ = std::min( currentfrag[i].fragscore_, workingcap );
			}

			// update tables
			if (currsemitag == "R") {
				//  A: counts (shared across rotamers)
				if (rotcount.find(currentfrag[1].aa_ ) == rotcount.end()) {
					rotcount[currentfrag[1].aa_].dimension( 36, 36 );
					rotcount[currentfrag[1].aa_] = 0;
				}
				rotcount[currentfrag[1].aa_](currentfrag[1].phibin_, currentfrag[1].psibin_) += currentfrag[1].weight_;

				//  B: energies
				for (int i=1; i<=(int)currentfrag.size(); ++i) {
					std::pair< int, int > rottag = std::make_pair( currentfrag[i].aa_ , currentfrag[i].rotidx_ );
					if (rotEsum.find( rottag ) == rotEsum.end()) {
						rotEsum[rottag].dimension( 36, 36 );
						rotEsum[rottag] = 0.0;
					}
					rotEsum[rottag](currentfrag[i].phibin_, currentfrag[i].psibin_) += currentfrag[i].weight_ * currentfrag[i].fragscore_;
				}
			} else {
				// hardcoded: semirot chi
				core::chemical::AA aa = (core::chemical::AA) currentfrag[1].aa_;
				bool semichi_is_2=true;
				if (aa==aa_gln || aa==aa_glu) semichi_is_2=false;

				//  energies only (counts apply to both)
				for (int i=1; i<=(int)currentfrag.size(); ++i) {
					std::pair< int, int > rottag = std::make_pair( currentfrag[i].aa_ , currentfrag[i].rotidx_ );
					if (semirotEsum.find( rottag ) == semirotEsum.end()) {
						semirotEsum[rottag].resize(36);
						for (int j=1;j<=36;++j) {
							semirotEsum[rottag][j].dimension(36,36);
							semirotEsum[rottag][j] = 0.0;
						}
					}
					semirotEsum[rottag][semichi_is_2 ? currentfrag[i].rotid_[1] : currentfrag[i].rotid_[2]](currentfrag[i].phibin_, currentfrag[i].psibin_)
						+= currentfrag[i].weight_ * currentfrag[i].fragscore_;
				}
			}

			bestscore=1e30;
			currentfrag.clear();
			currtag = fragline.id_;
			currsemitag = semitag;
		}

		bestscore = std::min( bestscore , fragline.fragscore_ );
		currentfrag.push_back( fragline );

		inscore >> semitag;
		inscore >> fragline.id_  >> fragline.aa_ >> fragline.rotidx_;
		inscore >> fragline.rotid_[0] >> fragline.rotid_[1] >> fragline.rotid_[2] >> fragline.rotid_[3];
		inscore >> fragline.phibin_ >> fragline.psibin_ >> fragline.weight_ >> fragline.fragscore_;
	}

	/////
	// 2 normalize, apply convolution
	// conv mask
	ObjexxFCL::FArray2D<core::Real> conv;
	core::Real convSum=0.0;
	conv.dimension(36,36); conv=0.0;
	for (int i=1; i<=36; ++i)
	for (int j=1; j<=36; ++j) {
		core::Real i_deg = 10.0*(i-1), j_deg = 10.0*(j-1);
		if (i_deg>=180.0) i_deg -=360.0;
		if (j_deg>=180.0) j_deg -=360.0;
		core::Real dist2 = ( i_deg*i_deg + j_deg*j_deg );
		core::Real x_ij = std::exp( -dist2/(2.0*smoothing*smoothing) );
		conv(i,j) = x_ij;
		convSum += x_ij;
	}
	for (int i=1; i<=36; ++i)
	for (int j=1; j<=36; ++j) {
		conv(i,j) /= convSum;
	}
	ObjexxFCL::FArray2D< std::complex<core::Real> > Fconv;
	numeric::fourier::fft2(conv, Fconv);

	ObjexxFCL::FArray2D< std::complex<core::Real> > Frot;

	/////
	// rotameric smoothing
	utility::vector1< ObjexxFCL::FArray2D<core::Real> > probSums(20, ObjexxFCL::FArray2D<core::Real>(36,36) );
	for (int i=1; i<=20; ++i) { probSums[i] = 0.0; };
	for (std::map< std::pair< int, int >,ObjexxFCL::FArray2D<core::Real> >::iterator it=rotEsum.begin(); it!=rotEsum.end(); ++it) {
		std::pair< int, int > rottag = it->first;
		ObjexxFCL::FArray2D<core::Real> &data = it->second;

		for (int i=1; i<=36; ++i)
		for (int j=1; j<=36; ++j) {
			// strategy 1: m estimates
			//core::Size count = rotcount[rottag.first](i,j) + MEST;
			//data(i,j) /= count;

			// strategy 2: don't correct sparse bins
			core::Size count = rotcount[rottag.first](i,j);
			if (count > 10) {
				data(i,j) /= count;
			} else {
				data(i,j) = 0.0;
			}

			// smooth in probability space
			data(i,j) = std::exp( -scale*data(i,j) );

			// normalize per-AA
			probSums[rottag.first](i,j) += data(i,j);
		}
	}
	for (std::map< std::pair< int, int >,ObjexxFCL::FArray2D<core::Real> >::iterator it=rotEsum.begin(); it!=rotEsum.end(); ++it) {
		std::pair< int, int > rottag = it->first;
		ObjexxFCL::FArray2D<core::Real> &data = it->second;

		for (int i=1; i<=36; ++i)
		for (int j=1; j<=36; ++j) {
			// normalize per-AA
			data(i,j) /= probSums[rottag.first](i,j);
		}

		numeric::fourier::fft2(data, Frot);

		for (int i=1; i<=36; ++i)
		for (int j=1; j<=36; ++j) {
			Frot(i,j) *= std::conj( Fconv(i,j) );
		}
		numeric::fourier::ifft2(Frot, data);

		dump_table ( data, "data_smooth.m",
			"data_" + utility::to_string( (core::chemical::AA)rottag.first ) + "_" + utility::to_string( rottag.second ) );
	}

	/////
	// semirotameric smoothing
	for (int i=1; i<=20; ++i) { probSums[i] = 0.0; };
	for (std::map< std::pair< int, int >, utility::vector1< ObjexxFCL::FArray2D<core::Real> > >::iterator it=semirotEsum.begin(); it!=semirotEsum.end(); ++it) {
		std::pair< int, int > rottag = it->first;
		for (int k=1; k<=36; ++k) {
			ObjexxFCL::FArray2D<core::Real> &data = it->second[k];
			for (int i=1; i<=36; ++i)
			for (int j=1; j<=36; ++j) {
				//core::Size count = rotcount[rottag.first](i,j) + MEST;
				//data(i,j) /= count;

				core::Size count = rotcount[rottag.first](i,j);
				if (count > 10) {
					data(i,j) /= count;
				} else {
					data(i,j) = 0.0;
				}

				// smooth in probability space
				data(i,j) = std::exp( -scale*data(i,j) );

				// normalize per-AA
				probSums[rottag.first](i,j) += data(i,j);
			}
		}
	}
	for (std::map< std::pair< int, int >, utility::vector1< ObjexxFCL::FArray2D<core::Real> > >::iterator it=semirotEsum.begin(); it!=semirotEsum.end(); ++it) {
		for (int k=1; k<=36; ++k) {
			ObjexxFCL::FArray2D<core::Real> &data = it->second[k];
			std::pair< int, int > rottag = it->first;

			for (int i=1; i<=36; ++i)
			for (int j=1; j<=36; ++j) {
				// normalize per-AA
				data(i,j) /= probSums[rottag.first](i,j);
			}

			numeric::fourier::fft2(data, Frot);

			for (int i=1; i<=36; ++i)
			for (int j=1; j<=36; ++j) {
				Frot(i,j) *= std::conj( Fconv(i,j) );
			}
			numeric::fourier::ifft2(Frot, data);

			dump_table ( data, "data_semirot.m",
				"data_" + utility::to_string( (core::chemical::AA)rottag.first ) + "_" + utility::to_string( rottag.second )+"_"+utility::to_string( k ) );
		}
	}


	/////
	// 3 compute corrected probabilities
	for (int aa = 1; aa<=20; aa++) {
		core::chemical::ResidueTypeCOPs restype =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->aa_map((core::chemical::AA)aa);

		core::pack::dunbrack::SingleResidueDunbrackLibraryCOP rotlib =
			utility::pointer::dynamic_pointer_cast<core::pack::dunbrack::SingleResidueDunbrackLibrary const> (
				core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( *(restype[1]) ) );

		if (!rotlib) continue;

		// rotameric...
		for (int x=1; x<=36; ++x)
		for (int y=1; y<=36; ++y) {
			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > allrots=
				rotlib->get_all_rotamer_samples((x-1)*10.0-180.0, (y-1)*10.0-180.0);
			core::Size nrot = allrots.size();
			if (nrot == 0) continue;

			core::Real probsum = 0.0;
			for ( Size ii = 1; ii <= nrot; ++ii ) {
				// use our own rotamer index since get_all_rotamer_samples changes indices based on phi/psi
				core::Size r1=0, r2=0, r3=0, r4=0;
				if (allrots[ii].nchi() > 0) r1 = allrots[ii].rot_well()[1];
				if (allrots[ii].nchi() > 1) r2 = allrots[ii].rot_well()[2];
				if (allrots[ii].nchi() > 2) r3 = allrots[ii].rot_well()[3];
				if (allrots[ii].nchi() > 3) r4 = allrots[ii].rot_well()[4];

				core::Size myrotidx = getRotID( r1,r2,r3,r4 );

				std::pair< int, int > rottag(aa,myrotidx);
				ObjexxFCL::FArray2D<core::Real> &data = rotEsum[rottag];

				// smooth Energy space
				//core::Real newP = allrots[ii].probability() * std::exp( -scale*data(x,y) );

if (aa == aa_ser && x==13 && y==14) {
	TR << "SER " << (x-1)*10.0-180.0 << " " << (y-1)*10.0-180.0 << " " << allrots[ii].probability() << " " << data(x,y) << " " << log( data(x,y)) / -scale << std::endl;
}
				// smooth Prob space
				core::Real newP = allrots[ii].probability() / data(x,y);


				data(x,y) = newP;
				probsum += newP;
			}

			for ( Size ii = 1; ii <= nrot; ++ii ) {
				// use our own rotamer index since get_all_rotamer_samples changes indices based on phi/psi
				core::Size r1=0, r2=0, r3=0, r4=0;
				if (allrots[ii].nchi() > 0) r1 = allrots[ii].rot_well()[1];
				if (allrots[ii].nchi() > 1) r2 = allrots[ii].rot_well()[2];
				if (allrots[ii].nchi() > 2) r3 = allrots[ii].rot_well()[3];
				if (allrots[ii].nchi() > 3) r4 = allrots[ii].rot_well()[4];

				core::Size myrotidx = getRotID( r1,r2,r3,r4 );

				std::pair< int, int > rottag(aa,myrotidx);
				ObjexxFCL::FArray2D<core::Real> &data = rotEsum[rottag];
				data(x,y) /= probsum;
if (aa == aa_ser && x==13 && y==14) {
	TR << "SER " << (x-1)*10.0-180.0 << " " << (y-1)*10.0-180.0 << allrots[ii].probability() << " " << data(x,y) << std::endl;
}
			}
		}

		/////
		// semirot
		if ( !is_semi_rot( (core::chemical::AA)aa ) ) continue;

		// make the residue we need
		core::conformation::ResidueOP working_res ( core::conformation::ResidueFactory::create_residue ( *(restype[1]) ) );

		// semirotameric sampling
		core::Real start=-180.0,step=10.0;
		if (aa == aa_asp || aa == aa_glu) { start = -90.0; step = 5.0; }
		if (aa == aa_phe || aa == aa_tyr) { start = -30.0; step = 5.0; }

		for (int x=1; x<=36; ++x)
		for (int y=1; y<=36; ++y) {
			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > allrots=
				rotlib->get_all_rotamer_samples((x-1)*10.0-180.0, (y-1)*10.0-180.0);
			core::Size nrot = allrots.size();
			if (nrot == 0) continue;

			utility::vector1< core::Real > bbtors(3,180.0);
			bbtors[1] = (x-1)*10.0-180.0;
			bbtors[2] = (y-1)*10.0-180.0;
			working_res->mainchain_torsions( bbtors ); // the only way to set this with an isolated residue ... ?

			core::Real probSum = 0.0;
			for ( Size ii = 1; ii <= nrot; ++ii ) {
				if (allrots[ii].rot_well()[allrots[ii].nchi()] != 1) continue;

				for (Size jj=1; jj<=allrots[ii].nchi()-1; ++jj)
					working_res->set_chi( (int)jj, allrots[ii].chi_mean()[jj] );

				// use our own rotamer index since get_all_rotamer_samples changes indices based on phi/psi
				core::Size r1=0, r2=0;
				if (allrots[ii].nchi() > 1) r1 = allrots[ii].rot_well()[1];
				if (allrots[ii].nchi() > 2) r2 = allrots[ii].rot_well()[2];
				core::Size myrotidx = getRotID( r1,r2 );

				std::pair< int, int > rottag(aa,myrotidx);

				for (Size jj=1; jj<=36; ++jj) {
					ObjexxFCL::FArray2D<core::Real> &data = semirotEsum[rottag][jj];

					dunbrack::RotamerLibraryScratchSpace scratch;
					working_res->set_chi( (int)allrots[ii].nchi(), start + (jj-1.0)*step );
					core::Real dunE = rotlib->rotamer_energy( *working_res, scratch );

					// smooth Energy space
					//core::Real newP = std::exp( -dunE ) * std::exp( -scale*data(x,y) );

					// smooth Prob space
					core::Real newP = std::exp( -dunE ) * data(x,y);

					data(x,y) = newP;
					probSum += newP;
				}
			}

			for ( Size ii = 1; ii <= nrot; ++ii ) {
				if (allrots[ii].rot_well()[allrots[ii].nchi()] != 1) continue;

				// use our own rotamer index since get_all_rotamer_samples changes indices based on phi/psi
				core::Size r1=0, r2=0;
				if (allrots[ii].nchi() > 1) r1 = allrots[ii].rot_well()[1];
				if (allrots[ii].nchi() > 2) r2 = allrots[ii].rot_well()[2];
				core::Size myrotidx = getRotID( r1,r2 );

				std::pair< int, int > rottag(aa,myrotidx);
				for (Size jj=1; jj<=36; ++jj) {
					ObjexxFCL::FArray2D<core::Real> &data = semirotEsum[rottag][jj];
					data(x,y) /= probSum;
				}
			}
		}
	}


	// 4 dump files
	for (std::map< std::pair< int, int >,ObjexxFCL::FArray2D<core::Real> >::iterator it=rotEsum.begin(); it!=rotEsum.end(); ++it) {
		std::pair< int, int > rottag = it->first;
		ObjexxFCL::FArray2D<core::Real> &data = it->second;
		dump_table ( data, "probs.m", "probs_" + utility::to_string( (core::chemical::AA)rottag.first ) + "_" + utility::to_string( rottag.second ) );
	}

	for (std::map< std::pair< int, int >, utility::vector1< ObjexxFCL::FArray2D<core::Real> > >::iterator it=semirotEsum.begin(); it!=semirotEsum.end(); ++it) {
		std::pair< int, int > rottag = it->first;
		for (Size jj=1; jj<=36; ++jj) {
			ObjexxFCL::FArray2D<core::Real> &data = it->second[jj];
			dump_table ( data, "probssemi.m", "probs_semi_"
					+ utility::to_string( (core::chemical::AA)rottag.first ) + "_" + utility::to_string( rottag.second ) + "_semi" + utility::to_string( jj ) );
		}
	}

	for (int aa = 1; aa<=20; aa++) {
		core::chemical::ResidueTypeCOPs restype =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->aa_map((core::chemical::AA)aa);

		core::pack::dunbrack::SingleResidueDunbrackLibraryCOP rotlib =
			utility::pointer::dynamic_pointer_cast<core::pack::dunbrack::SingleResidueDunbrackLibrary const> (
				core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( *(restype[1]) ) );

		if (!rotlib) continue;

		std::string aaname = core::chemical::name_from_aa((core::chemical::AA)aa);
		std::string outfile = aaname+".bbdep.rotamers.lib.gz";
		boost::algorithm::to_lower(outfile);
		utility::io::ozstream outdun(outfile.c_str() );

		// rotameric
		outdun<< "# Copyright (c) 2007-2010\n";
		outdun<< "# Maxim V. Shapovalov and Roland L. Dunbrack Jr.\n";
		outdun<< "# Fox Chase Cancer Center\n";
		outdun<< "# Philadelphia, PA, USA\n";
		//# T  Phi  Psi  Count    r1 r2 r3 r4 Probabil  chi1Val chi2Val chi3Val chi4Val   chi1Sig chi2Sig chi3Sig chi4Sig
		//VAL  -180 -180     1     1  0  0  0  0.707823    62.8     0.0     0.0     0.0       7.5     0.0     0.0     0.0
		for (int x=1; x<=37; ++x)
		for (int y=1; y<=37; ++y) {
			int xp=x, yp=y;
			if (xp==37) xp=1;
			if (yp==37) yp=1;

			int phi = (x-19)*10;
			int psi = (y-19)*10;

			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > allrots=
				rotlib->get_all_rotamer_samples( (core::Real)phi, (core::Real)psi );
			core::Size nrot = allrots.size();
			if (nrot == 0) continue;

			for ( Size ii = 1; ii <= nrot; ++ii ) {
				// use our own rotamer index since get_all_rotamer_samples changes indices based on phi/psi
				core::Size r1=0, r2=0, r3=0, r4=0;
				if (allrots[ii].nchi() > 0) r1 = allrots[ii].rot_well()[1];
				if (allrots[ii].nchi() > 1) r2 = allrots[ii].rot_well()[2];
				if (allrots[ii].nchi() > 2) r3 = allrots[ii].rot_well()[3];
				if (allrots[ii].nchi() > 3) r4 = allrots[ii].rot_well()[4];

				core::Size myrotidx = getRotID( r1,r2,r3,r4 );
				std::pair< int, int > rottag(aa,myrotidx);
				ObjexxFCL::FArray2D<core::Real> &data = rotEsum[rottag];
				core::Real rotprob = data(xp,yp); // rotameric probability

				outdun << aaname << " "
					<< phi << " "
					<< psi << " "
					<< "0" << " " // count ... this info is not stored in Rosetta (as far as I can tell)
					<< allrots[ii].rot_well()[1] << " "
					<< allrots[ii].rot_well()[2] << " "
					<< allrots[ii].rot_well()[3] << " "
					<< allrots[ii].rot_well()[4] << " "
					<< rotprob << " "
					<< allrots[ii].chi_mean()[1] << " "
					<< allrots[ii].chi_mean()[2] << " "
					<< allrots[ii].chi_mean()[3] << " "
					<< allrots[ii].chi_mean()[4] << " "
					<< allrots[ii].chi_sd()[1] << " "
					<< allrots[ii].chi_sd()[2] << " "
					<< allrots[ii].chi_sd()[3] << " "
					<< allrots[ii].chi_sd()[4] << "\n";
			}
		}

		/////
		// semirot
		if ( !is_semi_rot( (core::chemical::AA)aa) ) continue;

		std::string outfilesemi = aaname+".bbdep.densities.lib.gz";
		boost::algorithm::to_lower(outfilesemi);
		utility::io::ozstream outsemi(outfilesemi.c_str() );

		// semirotameric
		outsemi<< "# Copyright (c) 2007-2010\n";
		outsemi<< "# Maxim V. Shapovalov and Roland L. Dunbrack Jr.\n";
		outsemi<< "# Fox Chase Cancer Center\n";
		outsemi<< "# Philadelphia, PA, USA\n";

		bool semichi_is_2 = true;
		if ((core::chemical::AA)aa==aa_gln || (core::chemical::AA)aa==aa_glu) semichi_is_2=false;

		// # T Phi Psi Count r1 Probabil chi1Val chi1Sig -180 -170 -160 -150 -140 -130 -120 -110 -100 -90 -80 -70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170
		for (int x=1; x<=37; ++x)
		for (int y=1; y<=37; ++y) {
			int xp=x, yp=y;
			if (xp==37) xp=1;
			if (yp==37) yp=1;

			int phi = (x-19)*10;
			int psi = (y-19)*10;

			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > allrots=
					rotlib->get_all_rotamer_samples( (core::Real)phi, (core::Real)psi );

			core::Size nrot = allrots.size();
			if (nrot == 0) continue;

			for ( Size ii = 1; ii <= nrot; ++ii ) {
				if (allrots[ii].rot_well()[allrots[ii].nchi()] != 1) continue;  // ugly hack (matching calc_scores logic)

				// use our own rotamer index since get_all_rotamer_samples changes indices based on phi/psi
				core::Size r1=0, r2=0;
				if (allrots[ii].nchi() > 1) r1 = allrots[ii].rot_well()[1];
				if (allrots[ii].nchi() > 2) r2 = allrots[ii].rot_well()[2];

				core::Size myrotidx = getRotID( r1,r2 );

				std::pair< int, int > rottag(aa,myrotidx);

				core::Real rotprob=0.0;

				for ( Size jj = 1; jj <=36 ; ++jj ) {
					ObjexxFCL::FArray2D<core::Real> &data = semirotEsum[rottag][jj];
					core::Real rotprob_j = data(xp,yp);
					rotprob += rotprob_j;
				}
				for ( Size jj = 1; jj <=36 ; ++jj ) {
					ObjexxFCL::FArray2D<core::Real> &data = semirotEsum[rottag][jj];
					data(xp,yp) /= rotprob;
				}

				// rot info
				outsemi << aaname << " " << phi << " " << psi << " "
						<< "0" << " "; // count ... this info is not stored in Rosetta (as far as I can tell)
				if (semichi_is_2) {
					outsemi<< allrots[ii].rot_well()[1] << " "
						<< rotprob << " "
						<< allrots[ii].chi_mean()[1] << " "
						<< allrots[ii].chi_sd()[1] << " ";
				} else {  // semirotameric chi is chi3
					outsemi << allrots[ii].rot_well()[1] << " "
						<< allrots[ii].rot_well()[2] << " "
						<< rotprob << " "
						<< allrots[ii].chi_mean()[1] << " "
						<< allrots[ii].chi_mean()[2] << " "
						<< allrots[ii].chi_sd()[1] << " "
						<< allrots[ii].chi_sd()[2] << " ";
				}
				for ( Size jj = 1; jj <=36 ; ++jj ) {
					ObjexxFCL::FArray2D<core::Real> &data = semirotEsum[rottag][jj];
					outsemi << data(xp,yp) << " ";
				}
				outsemi << "\n";
			}
		}
	}
}


void
calc_scores() {
	using namespace core::chemical;

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

	// zero fa_dun
	if (scorefxn->get_weight( core::scoring::fa_dun ) > 0) {
		TR << "WARNING!!! fa_dun is non-zero in the scorefunction.  Setting it to non-zero" << std::endl;
		TR << "If you are sure you want fa_dun on, set the component terms instead: fa_dun_rot, fa_dun_dev, fa_dun_semi" << std::endl;
		scorefxn->set_weight( core::scoring::fa_dun , 0.0 );
	}

	core::Real fadundev_wt = scorefxn->get_weight( core::scoring::fa_dun_dev );
	scorefxn->set_weight( core::scoring::fa_dun_dev , 0.0 );

	if (fadundev_wt == 0 && option[ tors::rotmin ]() ) {
		TR << "WARNING!!! fa_dun_dev is zero in the scorefunction but minimization is enabled... are you sure this is what you want?" << std::endl;
	}

	// read silent file
	core::io::silent::SilentFileData sfd;
	sfd.read_file( option[ tors::fragfile ]() );

	// packing stuff
	core::pack::task::TaskFactoryOP tf (new core::pack::task::TaskFactory());
	tf->push_back( core::pack::task::operation::TaskOperationOP(new core::pack::task::operation::InitializeFromCommandline) );
	tf->push_back( core::pack::task::operation::TaskOperationOP(new core::pack::task::operation::RestrictToRepacking) );

	std::ofstream outscore( option[ tors::scorefile ]().c_str(), std::ofstream::out );

	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		std::string tag = iter->decoy_tag();
		utility::vector1< std::string > fields = utility::string_split( tag, '_' );

		TR.Debug << "Processing " << tag << std::endl;

		// NOTE: order must match what is written!
		std::string fragtag = fields[1];
		core::Size phibin = atoi( fields[2].c_str() );
		core::Size psibin = atoi( fields[3].c_str() );
		core::chemical::AA aa = (core::chemical::AA) atoi( fields[4].c_str() );
		core::Size weight = atoi( fields[5].c_str() );
		core::Size center = atoi( fields[6].c_str() );

		// get pose
		core::pose::Pose frag;
		iter->fill_pose( frag );

		core::conformation::Residue const &rsd = frag.residue( center );
		runtime_assert( rsd.aa() == aa ); // sanity check!

		core::pack::dunbrack::SingleResidueDunbrackLibraryCOP rotlib =
			utility::pointer::dynamic_pointer_cast<core::pack::dunbrack::SingleResidueDunbrackLibrary const> (
				core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( rsd.type() ) );

		if (!rotlib) continue;

		utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > allrots=
			rotlib->get_all_rotamer_samples(frag.phi(center), frag.psi(center));

		core::Size nrot = allrots.size();
		if (nrot == 0) continue;

		// rotameric samples
		for ( Size ii = 1; ii <= nrot; ++ii ) {
			for (Size jj=1; jj<=allrots[ii].nchi(); ++jj)
				frag.set_chi( (int)jj, center, allrots[ii].chi_mean()[jj] );

			core::pack::optimizeH( frag, *scorefxn );
			core::Real score_ii = (*scorefxn)(frag);
			if ( option[ tors::rotmin ]() ) {
				kinematics::MoveMap mm;
				mm.set_bb  ( false ); mm.set_chi ( false ); mm.set_jump( false );
				mm.set_chi ( center, true );
				core::optimization::MinimizerOptions options( "dfpmin", 0.1, true, false, false );
				core::optimization::AtomTreeMinimizer minimizer;

				scorefxn->set_weight( core::scoring::fa_dun_dev , fadundev_wt );
				minimizer.run( frag, mm, *scorefxn, options );
				scorefxn->set_weight( core::scoring::fa_dun_dev , 0.0 );
				score_ii = (*scorefxn)(frag);
			}

			// use our own rotamer index since get_all_rotamer_samples changes indices based on phi/psi
			core::Size r1=0, r2=0, r3=0, r4=0;
			if (allrots[ii].nchi() > 0) r1 = allrots[ii].rot_well()[1];
			if (allrots[ii].nchi() > 1) r2 = allrots[ii].rot_well()[2];
			if (allrots[ii].nchi() > 2) r3 = allrots[ii].rot_well()[3];
			if (allrots[ii].nchi() > 3) r4 = allrots[ii].rot_well()[4];

			core::Size myrotidx = getRotID( r1,r2,r3,r4 );

			outscore << "R " << fragtag << " " << (int)aa << " " << myrotidx << " ";
			for (int jj=1; jj<=4; ++jj) {
				if (jj<=(int)allrots[ii].nchi()) outscore << allrots[ii].rot_well()[jj]  << " ";
				else outscore << "0 ";
			}
			outscore << phibin << " " << psibin << " " << weight << " " << score_ii << "\n";
		}

		/////
		// semirot
		if ( !is_semi_rot( (core::chemical::AA)aa ) ) continue;

		core::Real start=-180.0,step=10.0;
		if (aa == aa_asp || aa == aa_glu) { start = -90.0; step = 5.0; }
		if (aa == aa_phe || aa == aa_tyr) { start = -30.0; step = 5.0; }

		for ( Size ii = 1; ii <= nrot; ++ii ) {
			if (allrots[ii].rot_well()[allrots[ii].nchi()] != 1) continue;  // ugly hack

			// use our own rotamer index since get_all_rotamer_samples changes indices based on phi/psi
			core::Size r1=0, r2=0;
			if (allrots[ii].nchi() > 1) r1 = allrots[ii].rot_well()[1];
			if (allrots[ii].nchi() > 2) r2 = allrots[ii].rot_well()[2];

			core::Size myrotidx = getRotID( r1,r2 );

			for (Size jj=1; jj<=36; ++jj) {
				// reset rotameric chis
				for (Size kk=1; kk<=allrots[ii].nchi()-1; ++kk) {
					frag.set_chi( (int)kk, center, allrots[ii].chi_mean()[kk] );
				}
				// set semirotameric chi
				core::Real chi_i = start + (jj-1.0)*step;
				frag.set_chi( (int)allrots[ii].nchi(), center, chi_i );

				core::pack::optimizeH( frag, *scorefxn );

				core::Real score_ii = (*scorefxn)(frag);
				if ( option[ tors::rotmin ]() ) {
					kinematics::MoveMap mm;
					mm.set_bb  ( false ); mm.set_chi ( false ); mm.set_jump( false );
					mm.set_chi ( center, true );
					core::optimization::MinimizerOptions options( "dfpmin", 0.1, true, false, false );
					core::optimization::AtomTreeMinimizer minimizer;

					scorefxn->set_weight( core::scoring::fa_dun_dev , fadundev_wt );
					minimizer.run( frag, mm, *scorefxn, options );
					scorefxn->set_weight( core::scoring::fa_dun_dev , 0.0 );

					// check to ensure we did not leave our semirotameric bin ... if we did, reset
					core::Real chi_new = frag.chi( (int)allrots[ii].nchi(), center );
					if (std::fabs( chi_i - chi_new ) > step/2.0 ) {
						frag.set_chi( (int)allrots[ii].nchi(), center, chi_i );
					}

					score_ii = (*scorefxn)(frag);
				}

				outscore << "S " << fragtag << " " << (int)aa << " " << myrotidx << " ";
				for (int kk=1; kk<=4; ++kk) {
					if (kk<(int)allrots[ii].nchi()) outscore << allrots[ii].rot_well()[kk]  << " ";
					else if (kk==(int)allrots[ii].nchi()) outscore << jj << " ";
					else outscore << "0 ";
				}
				outscore << phibin << " " << psibin << " " << weight << " " << score_ii << "\n";
			}
		}

	}
}


void
make_fragments() {
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

  core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();

	std::map < core::chemical::AA, ObjexxFCL::FArray2D< utility::vector1<FragInfo> > > database;

	for (int i=1; i<=20; ++i) {
		database[ (core::chemical::AA)i ].dimension(36,36);
	}

  while( input.has_another_pose() ) {
    core::pose::Pose pose;
    input.fill_pose( pose );

		std::string filetag = pose.data().get_const_ptr< basic::datacache::CacheableString >
				(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG)->str();

		// dssp
		core::scoring::dssp::Dssp dssp(pose);
		dssp.insert_ss_into_pose(pose);

		// foreach residue
		for (int i=3; i<=((int)pose.total_residue())-2; ++i) {
			// ensure we have a full fragment, no chainbreaks
			char ss = pose.secstruct(i);
	    core::pose::PoseOP frag( new core::pose::Pose);
			core::Real phi, psi;
			core::Size center=0;
			bool have_connected_frag=true;

			if (ss == 'L') {
				for (int j=i-2; j<i+2; ++j) {
					if ( pose.fold_tree().is_cutpoint(j) || !pose.residue(j).is_protein()) {
						have_connected_frag = false;
					}
				}
				if (!have_connected_frag || !pose.residue(i+2).is_protein()) continue;
				for (int j=i-2; j<=i+2; ++j) frag->append_residue_by_bond( pose.residue(j) );
				center = 3;
			}

			if (ss == 'H') {
				if (i-4<=0 || i+4 >= ((int)pose.total_residue())) continue;
				for (int j=i-4; j<i+4; ++j) {
					if ( pose.fold_tree().is_cutpoint(j) || !pose.residue(j).is_protein()) {
						have_connected_frag = false;
					}
				}
				if (!have_connected_frag || !pose.residue(i+4).is_protein()) continue;
				for (int j=i-4; j<=i+4; ++j) frag->append_residue_by_bond( pose.residue(j) );

				center = 5;
			}

			if (ss == 'E') {
				for (int j=i-2; j<i+2; ++j) {
					if ( pose.fold_tree().is_cutpoint(j) || !pose.residue(j).is_protein()) {
						have_connected_frag = false;
					}
				}
				if (!have_connected_frag || !pose.residue(i+2).is_protein()) continue;
				for (int j=i-2; j<=i+2; ++j) frag->append_residue_by_bond( pose.residue(j) );

				// find bb hbonds
				(*scorefxn)(pose); // to get interaction graph
				utility::vector1<bool> incl(pose.total_residue(), false);
				EnergyGraph const & energy_graph( pose.energies().energy_graph() );
				for (int j=i-2; j<=i+2; ++j) {
					core::conformation::Residue const &rsd_j( pose.residue(j) );

					for ( graph::Graph::EdgeListConstIter
							iru  = energy_graph.get_node(j)->const_edge_list_begin(),
							irue = energy_graph.get_node(j)->const_edge_list_end();
							iru != irue; ++iru ) {
						EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
						int k=(int)edge->get_other_ind(j);

						core::conformation::Residue const &rsd_k( pose.residue(k) );

						if (!rsd_k.is_protein()) continue;
						if (k>=i-2 && k<=i+2) continue;

						// check hbonding (O--H <= 3Å)
						core::Real dist1 = 999.0;
						if (rsd_k.aa() != core::chemical::aa_pro) {
							if (rsd_k.is_lower_terminus() && !rsd_k.has("H")) {
								dist1 =                  (rsd_j.atom("O").xyz() - rsd_k.atom("1H").xyz()).length_squared();
								dist1 = std::min( dist1, (rsd_j.atom("O").xyz() - rsd_k.atom("2H").xyz()).length_squared() );
								dist1 = std::min( dist1, (rsd_j.atom("O").xyz() - rsd_k.atom("3H").xyz()).length_squared() );
							} else {
								dist1 = (rsd_j.atom("O").xyz() - rsd_k.atom("H").xyz()).length_squared();
							}
						}
						core::Real dist2 = 999.0;
						if (rsd_j.aa() != core::chemical::aa_pro) {
							if (rsd_j.is_lower_terminus() && !rsd_j.has("H")) {
								dist2 =                  (rsd_j.atom("1H").xyz() - rsd_k.atom("O").xyz()).length_squared();
								dist2 = std::min( dist1, (rsd_j.atom("2H").xyz() - rsd_k.atom("O").xyz()).length_squared() );
								dist2 = std::min( dist1, (rsd_j.atom("3H").xyz() - rsd_k.atom("O").xyz()).length_squared() );
							} else {
								dist2 = (rsd_j.atom("H").xyz() - rsd_k.atom("O").xyz()).length_squared();
							}
						}
						if (dist1<=9 || dist2<=9) incl[k]=true;
					}
				}

				// zap singletons
				for (int j=2; j<=(int)incl.size()-1; ++j) incl[j] = incl[j] || (incl[j-1]&&incl[j+1]);
				for (int j=1; j<=(int)incl.size(); ++j) incl[j] = incl[j] && ( (j>1&&incl[j-1]) || (j<(int)incl.size()&&incl[j+1]));

				for (int j=1; j<=(int)incl.size(); ++j) {
					if (incl[j]) {
						if (j==1 || pose.fold_tree().is_cutpoint(j-1) || !incl[j-1]) frag->append_residue_by_jump( pose.residue(j), 1 );
						else frag->append_residue_by_bond( pose.residue(j) );
					}
				}

				center = 3;
			}

			assert( center != 0 );

			mutate_to_ala( *frag, center);
			phi = frag->phi(center);
			psi = frag->psi(center);

			core::Size phi_bin = getbin(phi);
			core::Size psi_bin = getbin(psi);

			// check for duplicates
			utility::vector1<FragInfo> & fragstore = database[ pose.residue(i).aa() ](phi_bin, psi_bin);
			FragInfo new_fraginfo(frag, 1.0, center, filetag);

			bool add_to_db=true;
			for (int j=1; j<=(int)fragstore.size(); ++j) {
				if ( distance(new_fraginfo, fragstore[j]) < option[ tors::cluster_radius ]() ) {
					fragstore[j].weight_ += 1;
					add_to_db = false;
				}
			}

			if (add_to_db) { fragstore.push_back( new_fraginfo ); }

			//frag.dump_pdb( "frag"+utility::to_string(i)+".pdb" );
		}
	}

  core::io::silent::SilentFileData sfd;

	core::Size fragidx=1;
	for (int i=1; i<=36; ++i)
	for (int j=1; j<=36; ++j) {
		for (int a=1; a<=20; ++a) {
			utility::vector1<FragInfo> & fragstore = database[ (core::chemical::AA)a ](i, j);
			core::Size nfrags = fragstore.size();

			for (int x=1; x<=(int)nfrags; ++x) {
				core::io::silent::SilentStructOP ss_out_all( new core::io::silent::BinarySilentStruct );
				std::string outfile =
						"S"+utility::to_string(fragidx++)+"_"
						+utility::to_string(i)+"_"+utility::to_string(j)+"_"+utility::to_string(a)+"_"
						+utility::to_string(fragstore[x].weight_)+"_"+utility::to_string(fragstore[x].center_);
				ss_out_all->fill_struct(*fragstore[x].pose_, outfile);
				ss_out_all->add_energy( "aa" , a );
				ss_out_all->add_energy( "phibin", i );
				ss_out_all->add_energy( "psibin", j );
				ss_out_all->add_energy( "weight", fragstore[x].weight_ );
				ss_out_all->add_energy( "center", fragstore[x].center_ );
				sfd.write_silent_struct( *ss_out_all, option[ tors::fragfile ]() );
			}
		}
	}
}


int
main( int argc, char * argv [] ) {
try {
	NEW_OPT(tors::mode, "mode (makefrags, calcscores)", "makefrags");
	NEW_OPT(tors::fragfile, "name of fragment file", "fragments.out");
	NEW_OPT(tors::scorefile, "name of fragment file", "scores.out");
	NEW_OPT(tors::cluster_radius, "Cluster fragments at this radius", 1.0);
	NEW_OPT(tors::rotmin, "Minimize each rotamer?", false);
	NEW_OPT(tors::smoothing, "Smoothing (in degrees phi/psi)", 10.0);
	NEW_OPT(tors::cap, "Cap on energy", 3.0);
	NEW_OPT(tors::scale, "Scaling factor on corrections", 0.56);

	devel::init( argc, argv );

	option[ in::missing_density_to_jump ].value(true); // force this
	option[ in::file::silent_struct_type ].value("binary"); // force this
	option[ packing::no_optH ].value(false); // force this
	option[ packing::flip_HNQ ].value(true); // force this

	std::string mode = option[ tors::mode ]();

	if (mode == "makefrags") {
		make_fragments();
	} else if (mode == "calcscores") {
		calc_scores();
	} else if (mode == "correctdun") {
		correct_dunbrack();
	}

} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
return 0;

}


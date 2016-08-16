// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/mc_convergence_checks/HPool.cc
/// @brief hierarchical pool
/// @author Yuan Liu (wendao@u.washington.edu)

#include <protocols/canonical_sampling/mc_convergence_checks/HPool.hh>
#include <protocols/toolbox/KCluster.hh>

#include <basic/Tracer.hh>
#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <protocols/toolbox/superimpose.hh>

#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray3P.hh>

#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>
#include <core/io/silent/SilentFileData.hh>

#include <numeric/xyzMatrix.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.hcluster" );

namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

//get desired lib file name from tag
std::string lib_full_path(std::string tag_orign)
{
	using namespace std;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static string subdir_prefix("sub_");
	static string libdir_prefix("c_");
	static string libdir_suffix("_lib");
	static string s("/");
	static int nofsubdir = option[cluster::K_n_sub];
	//static string rootpath = utility::file::FileName(option[mc::known_structures]()).path();
	static string defaultpath = utility::file::FileName(option[mc::known_structures]()).base() + libdir_suffix + s;

	string tag = tag_orign + ".0000";

	Size pos=1;
	Size len=tag.length();
	utility::vector1<Size> id_stack;
	while ( pos<len )
			{
		Size newpos = tag.find(".", pos+1);
		id_stack.push_back(atoi(tag.substr(pos+1,newpos-pos-1).c_str()));
		pos = newpos;
	}

	Size n = id_stack.size();
	if ( n == 1 ) return option[out::file::silent]();
	ostringstream fnstream;
	fnstream << "c_" << setfill ('0') << setw (5) << id_stack[n-1] << ".out";
	string fn(fnstream.str());

	for ( Size i=n-1; i>=1; i-- ) {
		ostringstream pathstream1;
		ostringstream pathstream2;

		//main dir
		if ( i>1 ) pathstream1 << libdir_prefix << setfill ('0') << setw (5) << id_stack[i-1] << libdir_suffix << s ;
		//sub dir
		pathstream2 << subdir_prefix << setfill ('0') << setw (3) << int((id_stack[i]-1)/nofsubdir) << s;

		fn = pathstream1.str()+pathstream2.str()+fn;
	}
	return defaultpath+fn;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
HPool_RMSD::HPool_RMSD(std::string silent_file, core::Size lv)
:Pool_RMSD(silent_file),
	silent_file_(silent_file),
	level_(lv),
	n_extra_(0)
{
	//set radius
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	radius_ = option[cluster::K_radius]()[level_];
	TR << "New cluster loaded: level=" << level_ << " radius=" << radius_ << std::endl;

	//subclusters lib path comes from lib_path_full
	//not this
	//lib_path_ = utility::file::FileName(silent_file).path();
	//lib_path_ += utility::file::FileName(silent_file).base() + "_lib";
	//TR << "Desired lib path: " << lib_path_ << std::endl;

	old_size_ = size();
	//get the tag prefix from the first tag
	//xxx.xxx.xxx.yyy -> xxx.xxx.xxx
	core::Size pos=1,oldpos=1;
	std::string t(tag(1));
	core::Size len=t.length();
	while ( pos<len ) {oldpos=pos; pos = t.find(".", pos+1);}
	tag_prefix_ = t.substr(0,oldpos+1).c_str();
	TR << "tag_prefix=" << tag_prefix_ << std::endl;

	has_child_ = (level_ == core::Size(option[cluster::K_level]())) ? false : true;
	//if (utility::file::file_exists(lib_path)) has_child_ = true;

	if ( has_child_ ) {
		//TR << "with subcluster" << std::endl;
		subpools_.resize(size());
	}

	build_pair_dis_matrix();
}

void HPool_RMSD::build_pair_dis_matrix()
{
	//save pair dist in a 1D vecotr, easy to expand
	pair_dis_.clear();
	pair_dis_.resize(size()*(size()-1)/2);//space

	for ( core::Size i=2; i<=size(); i++ ) {
		for ( core::Size j=1; j<=i-1; j++ ) {
			//cal id in 1D vector
			core::Size ndx = (i-1)*(i-2)/2 + j;
			pair_dis_[ndx] = std::sqrt(dist_square(i,j));
		}
	}
}

core::Real HPool_RMSD::dist_square(ObjexxFCL::FArray2_double &conf1, ObjexxFCL::FArray2_double &conf2)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static ObjexxFCL::FArray1D_double weights( natom(), 1.0 );
	static numeric::xyzMatrix< core::Real > R; //do not care

	//fit
	if ( !option[cluster::K_not_fit_xyz] ) {
		protocols::toolbox::fit_centered_coords(natom(), weights, conf1, conf2, R);
	}

	//cal dist
	core::Real sum=0.0;
	for ( Size i=1; i<=natom(); i++ ) {
		for ( core::Size d=1; d<=3; d++ ) {
			core::Real dx = conf1(d,i) - conf2(d,i);
			sum += dx*dx;
		}
	}
	return sum/natom();
}

core::Real HPool_RMSD::dist_square(core::Size ndx1, core::Size ndx2)
{
	ObjexxFCL::FArray2P_double conf1(coords()(1,1,ndx1), 3, natom());
	ObjexxFCL::FArray2P_double conf2(coords()(1,1,ndx2), 3, natom());
	return dist_square(conf1, conf2);
}
/*
void HPool_RMSD::add(core::pose::Pose& pose, std::string &tag)
{
//fill pose to coords
ObjexxFCL::FArray2D_double coord( 3, natom(), 0.0 );
runtime_assert( natom() == pose.total_residue() );
protocols::toolbox::fill_CA_coords( pose, coord );
add(coord, tag);
}

void HPool_RMSD::add(core::io::silent::SilentStruct& pss, std::string &tag)
{
//fill pose to coords
ObjexxFCL::FArray2D_double coord( pss.get_CA_xyz() );
runtime_assert( pss.nres() == natom() );
add(coord, tag);
}
*/
void HPool_RMSD::add(ObjexxFCL::FArray2D_double &coord, std::string &tag)
{
	//expand the pair dist matrix
	Pool_RMSD::add(coord, natom(), tag);//add new coords
	core::Size nd = size();
	pair_dis_.reserve(pair_dis_.size()+nd-1); //expand space
	for ( core::Size i=1; i<nd; i++ ) {
		//add new pair dist
		pair_dis_.push_back(std::sqrt(dist_square(nd,i)));
	}

	//add new subcluster
	if ( has_child_ ) {
		subpools_.push_back(*(new HPool_RMSD_OP));
	}
}

core::Real HPool_RMSD::get_pair_dist(core::Size i, core::Size j) const
{
	if ( i==j ) return 0.0;
	core::Size col = std::min(i,j);
	core::Size row = std::max(i,j);
	core::Size ndx = (row-1)*(row-2)/2 + col;
	return pair_dis_[ndx];
}

void HPool_RMSD::debug()
{
	for ( core::Size i=1; i<=size(); i++ ) {
		for ( core::Size j=1; j<=size(); j++ ) {
			std::cout << get_pair_dist(i,j) << " ";
		}
		std::cout << std::endl;
	}
}

core::Size HPool_RMSD::evaluate(
	core::pose::Pose& pose,
	core::Real resolution,
	std::string& best_decoy,
	core::Real& best_rmsd )
{
	if ( size() == 0 ) {
		best_decoy = "n/a";
		best_rmsd = 10000;
		return 0;
	}
	ObjexxFCL::FArray2D_double coord( 3, natom(), 0.0 );
	runtime_assert( natom() == pose.total_residue() );
	protocols::toolbox::fill_CA_coords( pose, coord );

	/*
	//TR << "eval: resolustion=" << resolution << std::endl;
	core::Size best_ndx=evaluate_core(coord, best_decoy, best_rmsd, 1);
	//TR << "best " << best_decoy << ":" << best_rmsd << std::endl;
	if (best_rmsd <= resolution)
	{
	//get the right one
	return 0;
	}
	else if ( (best_rmsd<radius_) && load_lib(best_ndx) )
	{
	//search from sub cluster
	subpools_[best_ndx]->evaluate(pose,resolution,best_decoy,best_rmsd);
	return 0;
	}
	else
	{
	//add a new center
	TR << "min_rmsd=" << best_rmsd << " out of R=" << radius_;
	TR << " Adding a new center ..." << std::endl;
	std::ostringstream tag;
	tag << "new." << (++n_extra_);
	best_decoy = tag.str();
	add(pose, best_decoy);
	best_rmsd = 0.0;
	return 0;
	}
	*/

	if ( evaluate(coord, resolution, best_decoy, best_rmsd) == 1 ) {
		//save this strcuture into subcluster's silent file
	}
	return 0;
}

core::Size HPool_RMSD::evaluate(
	core::io::silent::SilentStruct& pss,
	core::Real resolution,
	std::string& best_decoy,
	core::Real& best_rmsd )
{
	if ( size() == 0 ) {
		best_decoy = "n/a";
		best_rmsd = 10000;
		return 0;
	}

	ObjexxFCL::FArray2D_double coord( pss.get_CA_xyz() );
	runtime_assert( pss.nres() == natom() );

	/*
	//TR << "eval: resolustion=" << resolution << std::endl;
	core::Size best_ndx=evaluate_core(coord, best_decoy, best_rmsd, 1);
	//TR << "best " << best_decoy << ":" << best_rmsd << std::endl;
	if (best_rmsd <= resolution)
	{
	//get the right one
	return 0;
	}
	else if ( (best_rmsd<radius_) && load_lib(best_ndx) )
	{
	//search from sub cluster
	subpools_[best_ndx]->evaluate(pss,resolution,best_decoy,best_rmsd);
	return 0;
	}
	else
	{
	//add a new center
	TR << "min_rmsd=" << best_rmsd << " out of R=" << radius_;
	TR << " Adding a new center ..." << std::endl;
	std::ostringstream tag;
	tag << "new." << (++n_extra_);
	best_decoy = tag.str();
	add(pss, best_decoy);
	best_rmsd = 0.0;
	return 0;
	}
	*/

	//hierarchy
	if ( evaluate(coord, resolution, best_decoy, best_rmsd) == 1 ) {
		//save this strcuture into subcluster's silent file
		std::string filename(protocols::toolbox::file_full_path(best_decoy));
		utility::file::create_directory_recursive(utility::file::FileName(filename).path());
		core::io::silent::SilentFileData clusters;
		TR << "tag:" << best_decoy << " -> " << filename << std::endl;
		pss.set_decoy_tag(best_decoy);
		clusters.write_silent_struct(pss,filename);
	}
	return 0;
}

///////////////////////////
//hierarchy evaluate FArray2D
///////////////////////////
core::Size HPool_RMSD::evaluate(
	ObjexxFCL::FArray2D_double& coord,
	core::Real resolution,
	std::string& best_decoy,
	core::Real& best_rmsd )
{
	//TR << "eval: resolustion=" << resolution << std::endl;
	core::Size best_ndx=evaluate_core(coord, best_decoy, best_rmsd, 1);
	//TR << "best " << best_decoy << ":" << best_rmsd << std::endl;
	if ( best_rmsd <= resolution ) {
		//get the right one
		return 0;
	} else if ( (best_rmsd<radius_) && has_child_ ) {
		/*
		else if ( (best_rmsd<radius_) && load_lib(best_ndx) )
		{
		//search from sub cluster
		subpools_[best_ndx]->evaluate(pose,resolution,best_decoy,best_rmsd);
		}
		*/
		//06/20/2010
		//if this pool is not the lowest level but there is no
		//sub cluster for this center, then create a new one
		//
		if ( load_lib(best_ndx) ) {
			//successfully loaded the subcluster
			return subpools_[best_ndx]->evaluate(coord,resolution,best_decoy,best_rmsd);
		} else {
			//there is no element in the subcluster, create the file for it
			best_decoy = tag(best_ndx)+".00001";
			best_rmsd = 0.0;
			return 1;//need to save this strcture
		}
	} else {
		//add a new center
		TR << "min_rmsd=" << best_rmsd << " out of R=" << radius_;
		TR << " Adding a new center ..." << std::endl;

		std::ostringstream tag;
		//tag << "new." << (++n_extra_);
		//new tag, same format as old structure
		tag << tag_prefix_ << std::setfill ('0') << std::setw (5) << (++n_extra_+old_size_);
		best_decoy = tag.str();
		add(coord, best_decoy);
		best_rmsd = 0.0;
		return 1;//need to save this strcture
	}
	return 0;
}

//ref Pool_ConvergenceCheck.cc
//index indicates what index you start evaluation from
core::Size HPool_RMSD::evaluate_core(
	ObjexxFCL::FArray2D_double& coord,
	std::string& best_decoy,
	core::Real& best_rmsd ,
	core::Size index ) const
{
	if ( size() == 0 ) {
		best_decoy = "n/a";
		best_rmsd = 10000;
		return 0;
	}

	ObjexxFCL::FArray1D_double const weights( natom(), 1.0 );
	ObjexxFCL::FArray1D_double transvec( 3 );
	protocols::toolbox::reset_x( natom(), coord, weights, transvec );//center coordinates
	//n_atoms is simply # CA atoms in each "pose"
	core::Real invn( 1.0 / natom() );
	best_rmsd = 1000000;
	core::Size best_index( 0 );
	for ( core::Size i = index; i <= size(); i++ ) {
		//get xyz
		toolbox::Matrix R;
		ObjexxFCL::FArray2P_double xx2( coords()(1,1,i), 3, natom() );
		protocols::toolbox::fit_centered_coords( natom(), weights, xx2, coord, R );

		//triangle inequality
		if ( i>index ) {
			if ( get_pair_dist(best_index,i)/2.0>=best_rmsd ) {
				//TR << "jump" << std::endl;
				continue;
			}
		}

		//calculate rmsd
		core::Real rmsd( 0 );
		for ( core::Size n = index; n <= natom(); n++ ) {
			for ( core::Size d = 1; d<=3; ++d ) {
				rmsd += ( coord( d, n ) -  xx2( d, n ) ) * ( coord( d, n ) - xx2( d, n ) ) * invn;
			}
		}
		rmsd = std::sqrt( rmsd );
		if ( rmsd <= best_rmsd ) {
			best_index = i;
			best_rmsd = rmsd;
		}
	}

	best_decoy = tag(best_index);
	return best_index;
}

bool HPool_RMSD::load_lib(core::Size nc)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !has_child_ ) return false;

	if ( level_<=core::Size(option[cluster::K_deque_level]()) ) {
		//search in the list
		core::Size nlist = sub_ndx_deque_.size();
		core::Size i;
		//remember, deque is from 0 to n-1
		for ( i=0; i<nlist; i++ ) {
			if ( sub_ndx_deque_[i] == nc ) break;
		}

		if ( i>=nlist ) {
			//not in the list
			if ( nlist>=core::Size(option[cluster::K_deque_size]()) ) {
				//reach the limit of the list
				//removing the oldest subcluster in memory
				clear_lib(sub_ndx_deque_[0]);
				sub_ndx_deque_.pop_front();
			}
			//the last one is the newest one
			sub_ndx_deque_.push_back(nc);
		} else {
			//already in the list
			if ( i<nlist-1 ) { //not the last one
				sub_ndx_deque_.erase(sub_ndx_deque_.begin()+i);
				sub_ndx_deque_.push_back(nc);//put nc to the tail
			}
		}
	}//only for top K_deque_level level

	if ( subpools_[nc] ) return true;
	std::string filename = lib_full_path(tag(nc));
	TR << "try to load " << filename <<std::endl;
	if ( utility::file::file_exists(filename) ) {
		TR << "loading new cluster..." << std::endl;
		subpools_[nc] = protocols::canonical_sampling::mc_convergence_checks::HPool_RMSD_OP( new HPool_RMSD(filename, level_+1) );
		return true;
	}
	return false;
}

void HPool_RMSD::clear_lib(core::Size nc)
{
	subpools_[nc].reset(); // to NULL
}

}//mc_conv
}//moves
}//prot

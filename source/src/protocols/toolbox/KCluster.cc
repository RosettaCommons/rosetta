// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	protocols/toolbox/KCluster.hh
/// @brief	Fast clustering algorithm for large silent file
/// @author	Yuan Liu (wendao@u.washington.edu)

#include <protocols/toolbox/KCluster.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <core/io/silent/SilentFileData.hh>
//#include <core/io/silent/SilentStructFactory.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/random/random.hh>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <algorithm>

#include <utility/vector1.hh>


using namespace std;
using namespace core;

static thread_local basic::Tracer TR( "protocols.kcluster" );

namespace protocols {
namespace toolbox {

/// @details Auto-generated virtual destructor
KClusterData::~KClusterData() {}

/// @details Auto-generated virtual destructor
KClusterElement::~KClusterElement() {}

KClusterOP get_K_cluster_engine(const string &style)
{
    if (style == "GKC") return KClusterOP( new GreedyKCenter() );
    else if (style == "KMedoid") return KClusterOP( new KMedoid() );
    else
    {
        utility_exit_with_message("Undefined KCluster type!");
    }

    return KClusterOP( new GreedyKCenter() );
}

//get desired output file name from tag
string file_full_path(string tag)
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

	static string subdir_prefix("sub_");
	static string libdir_prefix("c_");
	static string libdir_suffix("_lib");
	static string s("/");
	static int nofsubdir = option[cluster::K_n_sub];

	//if the mc:hierarchical_pool is specified, then use that dir
	static string rootpath;
	static string defaultpath = utility::file::FileName(option[out::file::silent]()).base() + libdir_suffix + s;
	static string defaultfile;
	static bool flag=true;

	if (flag)
	{
		if (option[mc::hierarchical_pool].user())
		{
			rootpath = option[mc::hierarchical_pool]()+"_dir/sub_000/";
			defaultfile = "c_00001.out";
			defaultpath = "c_00001_lib/";
		}
		else
		{
			rootpath = utility::file::FileName(option[out::file::silent]()).path();
			defaultfile = utility::file::FileName(option[out::file::silent]()).local_name();
			defaultpath = utility::file::FileName(defaultfile).base() + libdir_suffix + s;
		}
		flag=false;
	}

	int pos=tag.find('.');
	int len=tag.length();
	utility::vector1<int> id_stack;
	while(pos<len && pos>0)
	{
		int newpos = tag.find('.', pos+1);
		id_stack.push_back(atoi(tag.substr(pos+1,newpos-pos-1).c_str()));
		pos = newpos;
	}

	Size n = id_stack.size();
	if (n == 1) return rootpath+defaultfile;
	ostringstream fnstream;
	fnstream << "c_" << setfill ('0') << setw (5) << id_stack[n-1] << ".out";
	string fn(fnstream.str());

	for (Size i=n-1; i>=1; i--)
	{
		ostringstream pathstream1;
		ostringstream pathstream2;

		//main dir
		if (i>1) pathstream1 << libdir_prefix << setfill ('0') << setw (5) << id_stack[i-1] << libdir_suffix << s ;
		//sub dir
		pathstream2 << subdir_prefix << setfill ('0') << setw (3) << int((id_stack[i]-1)/nofsubdir) << s;

		fn = pathstream1.str()+pathstream2.str()+fn;
	}
	return rootpath+defaultpath+fn;
}

string fix_tag_suffix(string str)
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
		static Size nlevel = option[cluster::K_level];
		Size n=count(str.begin(), str.end(), '.');
		for (Size i=n; i<nlevel; i++) str+=".1";

		int pos; //a silly hack
		while ((pos = str.find('_'))>0) str.replace(pos,1,1,'.');

		return str;
}

/// @brief assign a data into a cluster
void KClusterElement::assign_type_data(Size ndx_data, Size ndx_cluster, Real d)
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    assert(edit_mode);
    //assign type
    type_list_[ndx_data] = ndx_cluster;
    //save distance
    dist_list_[ndx_data] = d;
    //save the farest one
    if (d>max_distance_)
    {
      max_distance_=d;
      max_dist_ndx_=ndx_data;
    }
    //add into list, real data id
    //???
    //don't save the center into its cluster list
    //???
    if (option[cluster::K_redundant]() || (data_ndx_[ndx_data] != center_ndx_[ndx_cluster]))
    {
        subclusters_[ndx_cluster]->add_new_data(data_ndx_[ndx_data]);
    }
}

//////////////
//KClusterData
//////////////
KClusterData::KClusterData()
:ndata_(0),
natom_(0),
n_ca_atom_(0),
nfile_(0)
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
		using namespace protocols::loops;

    Loops loops( true );
    natom_ = loops.loop_size();
    if (natom_>0)
    {
        //specified loop region
        for( Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it )
        {
            for (core::Size i=it->start(), end=it->stop(); i<=end; i++)
            {
                rmsd_ca_list_.push_back(i);
            }
        }
        runtime_assert(natom_==rmsd_ca_list_.size());
    }

    load_silent_files();
    runtime_assert(n_ca_atom_>0 && ndata_>0);

    //debug
    TR << "Finished loading database:" << endl;
    TR << "Number of data: " << ndata_ << endl;
    TR << "Number of file: " << nfile_ << endl;
}

Real KClusterData::dist_square(FA2d &conf1, FA2d &conf2)
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    static FA1Dd weights( natom_, 1.0 );
    static numeric::xyzMatrix< Real > R; //do not care
    //fit
    if (!option[cluster::K_not_fit_xyz]) protocols::toolbox::fit_centered_coords(natom_, weights, conf1, conf2, R);

    //cal dist
    Real sum=0.0;
    for (Size i=1; i<=natom_; i++)
    {
        for (Size d=1; d<=3; d++)
        {
            Real dx = conf1(d,i) - conf2(d,i);
            sum += dx*dx;
        }
    }
    return sum/natom_;
}

Real KClusterData::dist_square(Size ndx1, Size ndx2)
{
    FA2Pd conf1(dataset_.coords()(1,1,ndx1), 3, natom_);
    FA2Pd conf2(dataset_.coords()(1,1,ndx2), 3, natom_);
    return dist_square(conf1, conf2);
}

void KClusterData::load_silent_files()
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    TR << "Reading Silent Files ..." << endl;

    if (option[in::file::silent].user())
    {
        for (Size i=1, e=option[in::file::silent]().size(); i<=e; i++)
        {
            load_silent_file(option[in::file::silent]()[i], i);
        }
    }
    else if (option[in::file::silent_list].user())
    {
        for (Size i=1, e=option[in::file::silent_list]().size(); i<=e; i++)
        {
            load_silent_file(option[in::file::silent_list]()[i], i);
        }
    }
    else
    {
        //no silent file input
        utility_exit_with_message("Please specify the input silent file or list!");
    }

    if (!option[cluster::K_not_fit_xyz])
    {
        //dataset_.superimpose(); //!!! maybe a new method that just reset_x is better
        FA1Dd transvec( 3 );
        FA1Dd weights(natom_, 1.0);
        for (Size i=1; i<=ndata_; i++)
        {
            FA2Pd conf( dataset_.coords()(1,1,i), 3, natom_ );
            reset_x( natom_, conf, weights, transvec );
        }
    }
}

void KClusterData::load_silent_file(string silent_file, Size nfile)
{
    io::silent::SilentFileData sfd;
    sfd.read_file( silent_file );

    Size extra=sfd.size();

    if ( dataset_.n_decoys_max() < dataset_.n_decoys() + extra + 1 )
    {
        dataset_.reserve( dataset_.n_decoys() + extra );
    }

    Size ncoord=0;
    for ( io::silent::SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it )
    {
				FA2Dd original_xyz(it->get_CA_xyz());
        //check n_ca_atom_
        if ( n_ca_atom_==0 ) n_ca_atom_=it->nres();
        else runtime_assert( n_ca_atom_ == it->nres() );

        //setup final_xyz
        if (natom_==0)
        {
            //no loop specified, use all residues
            natom_ = n_ca_atom_;
            for (core::Size i=1; i<=n_ca_atom_; i++) rmsd_ca_list_.push_back(i);
        }
        //debug
        //TR << "n_ca_atom_ = " << n_ca_atom_ << std::endl;
        //TR << "natom_ = " << natom_ << std::endl;
        FA2Dd final_xyz(3, natom_, 0.0);

        for( core::Size i=1; i<=natom_; i++)
        {
            final_xyz(1,i) = original_xyz(1,rmsd_ca_list_[i]);
            final_xyz(2,i) = original_xyz(2,rmsd_ca_list_[i]);
            final_xyz(3,i) = original_xyz(3,rmsd_ca_list_[i]);
        }

        dataset_.push_back_CA_xyz( final_xyz, natom_ );

        //build tag
        ncoord++;
        ostringstream tag;
        tag << "d."  << setfill ('0') << setw (4) << nfile << "."
        << setfill ('0') << setw (8) << ncoord;
				TagList list;
				list.push_back(tag.str());
        tags_.push_back(list);
        TagList tagvec;
        tagvec.push_back(it->decoy_tag());
        original_tags_.push_back(tagvec);
				//				TagList fn;
				//				fn.push_back( silent_file );
				original_filenames_.push_back( silent_file );
    }

    //check natom_
    runtime_assert(dataset_.n_atoms()==natom_);
    //if (natom!=natom_) {
    //    utility_exit_with_message("Input silent file contains different protein size!");
    //}

    ndata_ += extra;
    assert(ndata_ == dataset_.n_decoys());
    nfile_++;
}

//for saving the cluster center
void KClusterData::mark_tags(KClusterElementOP c, string t)
{
    //mark of the centers in this element
    //if the subcluster's center_ndx_ is not empty, recursive
    Size nc = c->get_cur_ncluster();
    if (nc==0)return;//makesure this element has been clustered
    for (Size i=1; i<=nc; i++)
    {
        //for all center
        ostringstream tag;
        //tag << t << "." << setfill ('0') << setw (5) << i;
				tag << t << "." << i;
        //mark that have not been marked
        //if (tags_[c->get_center_ndx(i)].c_str()[0]=='d') tags_[c->get_center_ndx(i)]=tag.str();
				tags_[c->get_center_ndx(i)].push_back(tag.str());

        //if has subcluster
        mark_tags(c->get_subcluster(i), tag.str());
    }
}

void KClusterData::save_all_in_one()
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    utility::vector1<string> fnlist;

    if (option[in::file::silent].user())
    {
        fnlist=option[in::file::silent]();
    }
    else if (option[in::file::silent_list].user())
    {
        fnlist=option[in::file::silent_list]();
    }

    //save the center
    io::silent::SilentFileData clusters;
    string const silent_outfile(option[out::file::silent]());

    //read from file list
    Size count=0;
    for (Size i=1; i<=nfile_; i++)
    {
        //sort the struct order
        io::silent::SilentFileData sfd;
        sfd.read_file( fnlist[i] );

        for ( io::silent::SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it )
        {
            count++;
            //if (tags_[count].c_str()[0]=='d') continue; //skip data
            //get center that has been marked
            //it->set_decoy_tag(tags_[count]);
            //clusters.write_silent_struct(**it,silent_outfile);
						for(Size i=2, nt=tags_[count].size(); i<=nt; i++)
						{
							//it->set_decoy_tag(tags_[count][i]);
							it->add_string_value("cluster_id",tags_[count][i]);
							clusters.write_silent_struct(**it,silent_outfile);
						}
        }
    }
    assert(count==ndata_);
}

void KClusterData::save_cluster_tree()
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    utility::vector1<string> fnlist;

    if (option[in::file::silent].user())
    {
        fnlist=option[in::file::silent]();
    }
    else if (option[in::file::silent_list].user())
    {
        fnlist=option[in::file::silent_list]();
    }

    //save the center
    io::silent::SilentFileData clusters;

    //read from file list
    Size count=0;
    for (Size i=1; i<=nfile_; i++)
    {
        //sort the struct order
        io::silent::SilentFileData sfd;
        sfd.read_file( fnlist[i] );

        for ( io::silent::SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it )
        {
            count++;
            //if (tags_[count].c_str()[0]=='d') continue; //skip data
            //get center that has been marked
			for(Size i=2, nt=tags_[count].size(); i<=nt; i++)
			{
				string filename(file_full_path(tags_[count][i]));
				utility::file::create_directory_recursive(utility::file::FileName(filename).path());

				if (option[cluster::K_redundant])
				{

					//it->set_decoy_tag(fix_tag_suffix(tags_[count][i]));
					it->add_string_value("cluster_id",fix_tag_suffix(tags_[count][i]));
				}
				else
				{
					//it->set_decoy_tag(tags_[count][i]);
					it->add_string_value("cluster_id",tags_[count][i]);
				}

				if (option[ cluster::K_save_headers ])
				{
					std::ofstream os;
					os.open( filename.c_str(), std::ios::app);
					(*it)->print_header( os ); //this outputs a header to the silent file.
					os.close();
				}
				clusters.write_silent_struct(**it,filename);
			}
        }
    }
    assert(count==ndata_);
}

void
KClusterData::show_cluster_assignments() {
	std:: cout << "outputting cluster assignments..." << std::endl;
	for( core::Size ii = 1; ii <= tags_.size(); ii++ ){
		std::cout << "size of tagslist for index: " << ii << " is " << tags_[ii].size() << std::endl;
		for( core::Size jj = 1; jj <= tags_[ii].size(); jj++ ){
			std::cout << "ndx belonging to cluster  " << ii << " " << tags_[ii][jj] << std::endl;
		}
	}

}

//////////
//KCluster
//////////
KCluster::KCluster()
{
    //using namespace basic::options;
    //using namespace basic::options::OptionKeys;
    //n_cluster_ = option[ cluster::K_n_cluster ];
    n_cluster_ = 0;
}

KCluster::~KCluster(){}

void KCluster::cluster(KClusterElementOP c, KClusterData& d, Size first)
{
    TR << endl;
    TR << "*********** Job ***********"  << endl;
    init(c, first); //init cluster(s), randomly or keep the old cluster
    TR << "Database: " << c->get_ndata() << " structures" << endl;
    TR << "Clustering ..." << endl;

    bool flag(whoami());
    Real old_score = 999999;
    do
    {
        Real score = assign(c,d);
        TR << "Score: " << score << " N_cluster: " << c->get_cur_ncluster() << endl;

        if (flag)
        {
					TR << "in KCenter mode" << std::endl;
            //KCenter
            if (c->get_cur_ncluster()==n_cluster_)
            {
                break;
            }

            if (score<=get_threshold())break;
        }
        else
        {
            //KMedoid
            if (old_score-score<=get_threshold()) break;
            old_score = score;
        }

        update(c,d);
    }
    while (true);
    TR << "Finish!" << endl;
}

//////////
//KMedoid
//////////
KMedoid::KMedoid():KCluster()
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    //threshold_ = option[cluster::K_threshold];
    threshold_ = 0.02;
}

Real KMedoid::get_threshold()
{
    return threshold_;
}

void KMedoid::set_threshold(Real t)
{
    threshold_ = t;
}

bool KMedoid::whoami()
{
    TR << "I am K-Medoid Algorithm!" << endl;
    return false;
}

void KMedoid::init(KClusterElementOP c, Size first)
{
    cur_ncluster_ = c->get_cur_ncluster();
    Size nd = c->get_ndata();

    if (cur_ncluster_ == 0)
    {
        cur_ncluster_ = min(n_cluster_, nd);
        assert(cur_ncluster_ > 0);
        //randomly choose n center
        TR << "Empty cluster, randomly choose center" << endl;
        for (Size i=1; i<=cur_ncluster_; i++)
        {
            Size newcenter;
						if (i==1 && first>0)
						{
								//specify the first one
								newcenter = first;
						}
						else
						{
								newcenter = c->get_data_ndx(static_cast<int>( numeric::random::rg().uniform() * nd + 1 ));
						}

            Size flag = i;
            for (Size j=1; j<i; j++)
            {
                //make sure there's no replica
                if (c->get_center_ndx(j) == newcenter)
                {
                    i--;
                    break;
                }
            }

            if (i==flag) c->add_new_cluster(newcenter);
        }
    }
    else
    {
        TR << "Load clusters with " << cur_ncluster_ << " centers" << endl;
    }
}

void KMedoid::copy_coord(Size len, FA2d &src, FA2d &dst)
{
    for (Size i=1; i<=len; i++)
    {
        for (Size j=1; j<=3; j++)
        {
            dst(j,i) = src(j,i);
        }
    }
}

Real KMedoid::assign(KClusterElementOP c, KClusterData& d)
{
    //assign all data to the nearest center
    //keep the coord that fit the nearest center
    Size nd = c->get_ndata();
    Size na = d.get_natom();
    FA2Dd coord(3, na, 0.0);
//TR << "DEBUG: ncluster" << cur_ncluster_ << endl;
    //DEBUG
    //TR<<"center list:"<< endl;
    //for (Size i=1; i<=cur_ncluster_; i++)
    //{
    //	TR <<c->get_center_ndx(i) << endl;
    //}

    for (Size nc=1; nc<=cur_ncluster_; nc++)
    {
        //for each cluster center

        //build the center dist list
        ObjexxFCL::FArray1D_double center_dis_list(nc,0.0);
        for (Size i=1; i<nc; i++)
        {
            //TR << "Center list: c" << i <<":"<<c->get_center_ndx(i) << ", c" << nc <<":"<<c->get_center_ndx(nc)<< " d: ";
            center_dis_list(i) = sqrt(d.dist_square(
                                          c->get_center_ndx(i),
                                          c->get_center_ndx(nc) ));
            //TR << center_dis_list(i) << endl;
        }

        c->clear();
        for (Size i=1; i<=nd; i++)
        {
            //assign the center to its cluster
            if (c->get_data_ndx(i)==c->get_center_ndx(nc))
            {
                c->assign_type_data(i, nc, 0.0);
                continue;
            }

            //for each data
            if (nc==1)
            {
                //assign all to the the first
                c->assign_type_data(i, nc, sqrt(d.dist_square(
									c->get_center_ndx(nc), c->get_ndx_list()[i]
                                                )));
                continue;
            }

            Real d_old = c->get_distance(i);

            //apply triangle inequality
            if (d_old <= center_dis_list(c->get_type(i))/2.0 )
            {
                //dist to the new center is farer than the old one
                c->assign_type_data(i,c->get_type(i),d_old);
                continue;
            }


            FA2Pd src(d.coords()(1,1,c->get_ndx_list()[i]),3,na);
            copy_coord(na, src, coord);//save old structure
            Real d_new = sqrt(d.dist_square(c->get_center_ndx(nc),c->get_ndx_list()[i]));//it would be align to the new center
            if (d_new<d_old)
            {
                c->assign_type_data(i, nc, d_new);
            }
            else
            {
                copy_coord(na, coord, src);//re, make sure every structure align to their center
                c->assign_type_data(i,c->get_type(i),d_old);
            }
        }
        c->check();
    }

    //calculate score
    Real sum=0.0;
    for (Size i=1; i<=nd; i++)
    {
        sum += c->get_distance(i);
    }
    return sum/nd;
}

void KMedoid::update(KClusterElementOP c, KClusterData& d)
{
    //TR << "Update..." << endl;
    //find the coord most closed to the average of the cluster
    Size na = d.get_natom();
    for (Size i=1; i<=cur_ncluster_; i++)
    {
        //for each cluster
        const utility::vector1<Size> &list(c->get_ndx_list(i));
        Size nd = list.size();
        FA2Dd sum(3, na, 0.0);
        for (Size j=1; j<=nd; j++)
        {
            FA2Pd coord(d.coords()(1,1,list[j]),3,na);
            for (Size m=1; m<=na; m++)
            {
                for (Size n=1; n<=3; n++)
                {
                    sum(n,m)+=coord(n,m);
                }
            }
        }

        /*/debug
        TR << "Cluster: " << i << " -> ";
        for (Size j=1; j<=nd; j++)
        {
        	TR << list[j] << " ";
        }
        TR << endl;*/

        //average
        for (Size m=1; m<=na; m++)
        {
            for (Size n=1; n<=3; n++)
            {
                sum(n,m)/=nd;
            }
        }

        //find the nearest
        Real mindist2 = 999999;
        Size nearest = 0;
        for (Size j=1; j<=nd; j++)
        {
            //DEBUG
            //TR << "U list: " << j << ":" << list[j];
            FA2Pd coord(d.coords()(1,1,list[j]),3,na);
            Real d2_new = d.dist_square(coord,sum);
            //TR << " d_new: " << d2_new << " -- " << mindist2 << endl;
            if (d2_new < mindist2)
            {
                mindist2 = d2_new;
                nearest = list[j];
                //TR << nearest << endl;
            }
        }

        //TR << "final: "<< nearest << endl;
        //re set the center
        assert(nearest>0);
        c->set_cluster(i,nearest);
    }
}

///////////////
//GreedyKCenter
///////////////
GreedyKCenter::GreedyKCenter():KCluster()
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    //radius_ = option[ cluster::K_radius ];
    radius_ = 2.0;
}

//GreedyKCenter::~GreedyKCenter(){}
bool GreedyKCenter::whoami()
{
    TR << "I am Approximate K-Center Algorithm!" << endl;
    return true;
}

Real GreedyKCenter::get_threshold()
{
    return radius_;
}

void GreedyKCenter::set_threshold(Real r)
{
    radius_ = r;
}

void GreedyKCenter::init(KClusterElementOP c, Size first)
{
    //random select a center
    TR << "Initializing ..." << endl;
    assert(c->get_cur_ncluster()==0);//only begin from an empty clusters
    Size center;
	 	if (first>0)
		{
				center = first;
		}
		else
		{
				center = c->get_data_ndx(static_cast< int >( numeric::random::rg().uniform() * c->get_ndata() + 1 ));
		}
    //debug
    //TR << "Rand: " << center << endl;
    c->add_new_cluster(center);
}

Real GreedyKCenter::assign(KClusterElementOP c, KClusterData& d)
{
    //assign all data to the nearest center
    Size nc = c->get_cur_ncluster();

    //save center dist
    ObjexxFCL::FArray1D_double center_dis_list(nc,0.0);
    for (Size i=1; i<nc; i++)
    {
        center_dis_list(i) = sqrt(d.dist_square(
                                      c->get_center_ndx(i),
                                      c->get_center_ndx(nc) ));
        //TR << "Center list: c" << i << ", c" << nc << " d: " << center_dis_list(i) << endl;
    }

    c->clear();
    Size nd = c->get_ndata();
    for (Size i=1; i<=nd; i++)
    {
        //for each data
        Real d_old = c->get_distance(i);
        if (nc>1)
        {
            if (d_old <= center_dis_list(c->get_type(i))/2.0 )
            {
                //dist to the new center is farer than the old one
                c->assign_type_data(i,c->get_type(i),d_old);
                continue;
            }
        }
        Real d_new = sqrt(d.dist_square(c->get_center_ndx(nc),c->get_ndx_list()[i]));
        if (d_new<d_old)
        {
            c->assign_type_data(i, nc, d_new);
        }
        else
        {
            c->assign_type_data(i,c->get_type(i),d_old);
        }
    }
    c->check();

    return c->get_max_distance();
}

void GreedyKCenter::update(KClusterElementOP c, KClusterData&)
{
    if (c->get_cur_ncluster()==n_cluster_)
    {
				assert(false);
        //add a psesudo cluster center, which will be remove later
        c->add_new_cluster(0);
        return;
    }

    //debug
    //TR << "Farest one:" << c->get_max_dist_ndx() << endl;
    c->add_new_cluster(c->get_ndx_list()[c->get_max_dist_ndx()]);
    return;
}

}//toolbox
}//protocols

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

#ifndef INCLUDED_protocols_toolbox_KCluster_hh
#define INCLUDED_protocols_toolbox_KCluster_hh

#include <protocols/toolbox/KCluster.fwd.hh>
#include <protocols/toolbox/DecoySetEvaluation.hh>

#include <core/types.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray3P.hh>

#include <string>

namespace protocols {
namespace toolbox {

////////////////////////////////////////////////////////////////////////
// Hierachical Cluster Data Structure
////////////////////////////////////////////////////////////////////////
class KClusterElement : public utility::pointer::ReferenceCount
{
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~KClusterElement();
    typedef utility::vector1< core::Size > ClusterNdxList;
    typedef utility::vector1< core::Size > ClusterTypList;
    typedef utility::vector1< core::Real > ClusterDisList;

    KClusterElement():
    max_distance_(-1.0),
    max_dist_ndx_(0),
    edit_mode(false){}

    KClusterElement(core::Size nd):
    data_ndx_(nd,0),
    type_list_(nd,0),
    dist_list_(nd,999999),
    max_distance_(-1.0),
    max_dist_ndx_(0),
    edit_mode(false)
    {
    	for(core::Size i=1; i<=nd; i++)data_ndx_[i]=i;
    }

    /// @brief assign a data into a cluster
    void assign_type_data(core::Size, core::Size, core::Real);

    /// @brief add a new struture ndx to the data_ndx_
    void add_new_data(core::Size ndx_data)
    {
    	data_ndx_.push_back(ndx_data);
    	type_list_.push_back(0);
    	dist_list_.push_back(999999);

    }

    /// @brief add a new cluster center's data_ndx
    void add_new_cluster(core::Size ndx_data)
    {
        center_ndx_.push_back(ndx_data);
        subclusters_.push_back(new KClusterElement());
    }

    /// @brief set a cluster center's data_ndx
    void set_cluster(core::Size ndx_cluster, core::Size ndx_data)
    {
    	center_ndx_[ndx_cluster] = ndx_data;
    }

    /// @brief return data's type(local cluster index)
	core::Size get_type(core::Size ndx_data) const
    {
		return type_list_[ndx_data];
    }

	/// @brief return distance between data and center
    core::Real get_distance(core::Size ndx_data) const
    {
        return dist_list_[ndx_data];
    }

    /// @brief return cluster center's data_ndx
    core::Size get_center_ndx(core::Size ndx_cluster) const
    {
        return center_ndx_[ndx_cluster];
    }

    /// @brief return the data ndx list of this cluster
    const utility::vector1<core::Size> &get_ndx_list() const
	{
    	return data_ndx_;
	}

    /// @brief return the ndx list of sub-cluster
    const utility::vector1<core::Size> &get_ndx_list(core::Size c) const
    {
    	return subclusters_[c]->get_ndx_list();
    }

    /// @brief return the subcluster
    KClusterElementOP get_subcluster(core::Size nc) const
    {
    	return subclusters_[nc];
    }

	core::Size ncluster() {
		return center_ndx_.size();
	}

    /// @brief return current cluster number
    core::Size get_cur_ncluster() const
    {
        return subclusters_.size();
    }

    /// @brief return current data number
    core::Size get_ndata() const
    {
        return data_ndx_.size();
    }

    /// @brief return current data number
    core::Size get_data_ndx(core::Size ndx_data) const
    {
        return data_ndx_[ndx_data];
    }

    core::Real get_max_distance() const
    {
    	return max_distance_;
    }

    core::Size get_max_dist_ndx() const
   	{
    	return max_dist_ndx_;
	}

    /// @brief clean the data list
    void clear_data()
    {
    	data_ndx_.clear();
    	type_list_.clear();
    	dist_list_.clear();
    }
    /// @brief clean the subcluster's list, open edit mode
    void clear()
    {
        assert(!edit_mode);
        for (core::Size i=1, e=subclusters_.size(); i<=e; i++)
        {
            subclusters_[i]->clear_data();
        }
        edit_mode=true;
        max_distance_=-1.0;
        max_dist_ndx_=0;
    }

    /// @brief check the list, close edit mode
    void check()
    {
        assert(edit_mode);
        edit_mode=false;
    }

private:
    ClusterNdxList data_ndx_; //the ndx of each data; root empty
    ClusterTypList type_list_; //thet type of each data
    ClusterDisList dist_list_; //the distance of each data to center

    ClusterNdxList center_ndx_; //the data ndx of each cluster center
    utility::vector1< KClusterElementOP > subclusters_; //array of sub-cluster

	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Size n_data_; //number of data
    core::Real max_distance_; //save the farest data's d
    core::Size max_dist_ndx_; //save the farest data's ndx

    bool edit_mode;
};


/// @brief database of a K-style Clustering algorithm
class KClusterData : public utility::pointer::ReferenceCount
{
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~KClusterData();
    typedef ObjexxFCL::FArray1D_double FA1Dd;
    typedef ObjexxFCL::FArray2D_double FA2Dd;
    typedef ObjexxFCL::FArray2_double FA2d;
    typedef ObjexxFCL::FArray2P_double FA2Pd;
    typedef ObjexxFCL::FArray3P_double FA3Pd;
    typedef utility::vector1< std::string > TagList;

    KClusterData();

    void load_silent_files();
    void load_silent_file(std::string, core::Size);
    //void save_cluster();

    core::Size get_ndata() const
    {
        return ndata_;
    }
    core::Size get_natom() const
    {
    	return n_ca_atom_;
    }

    std::string get_tag(core::Size i) {
    	return original_tags_[i][1];
    }

  	std::string source_filename(core::Size i) {
	  	return original_filenames_[i];
	  }

    FA3Pd coords()
    {
    	return dataset_.coords();
    }

    //for saving the cluster center
    void mark_tags(KClusterElementOP, std::string);
    void save_all_in_one();
    void save_cluster_tree();

    core::Real dist_square(FA2d &, FA2d &);
    core::Real dist_square(core::Size, core::Size);

    void show_cluster_assignments(); //ek added 2-4-2011

private:
    protocols::toolbox::DecoySetEvaluation dataset_; //data
    core::Size ndata_; //number of data
    core::Size natom_; //number of atom for calculating rmsd
    core::Size n_ca_atom_; //number of CA in each data
    core::Size nfile_; //number of silent file
    utility::vector1< TagList > tags_; //save the tag for each data
    utility::vector1< TagList > original_tags_; // map the cluster tag back to the original decoy tag
	utility::vector1< std::string > original_filenames_; //map cluster tag back to original filename. (this is because often i have lots of silent-files with decoys of the same name )
    utility::vector1< core::Size > rmsd_ca_list_;
};


/// @brief basic class for performing a K-style Clustering algorithm
class KCluster : public utility::pointer::ReferenceCount
{
public:
    typedef ObjexxFCL::FArray2_double FA2d;
    typedef ObjexxFCL::FArray2P_double FA2Pd;
    typedef ObjexxFCL::FArray2D_double FA2Dd;

    KCluster();
    virtual ~KCluster();
    virtual void init(KClusterElementOP, core::Size first=0)=0;
    virtual void update(KClusterElementOP, KClusterData&)=0;
    virtual bool whoami()=0;
    virtual core::Real get_threshold()=0;
    virtual void set_threshold(core::Real)=0;

    //virtual Real assign(KClusterData&)=0;//typical assign (K-means)
    virtual core::Real assign( KClusterElementOP, KClusterData&)=0;

    void cluster( KClusterElementOP, KClusterData&, core::Size first=0);
    void set_ncluster(core::Size nc){n_cluster_=nc;}

protected:
    core::Size n_cluster_;
};


/// @brief Typical K-Medoids Clustering Algorithm
class KMedoid : public KCluster
{
public:
    KMedoid();

    bool whoami();
    core::Real get_threshold();
    void init(KClusterElementOP, core::Size first=0);
    core::Real assign(KClusterElementOP, KClusterData&);
    void update(KClusterElementOP, KClusterData&);
    void set_threshold(core::Real);
protected:
    core::Size cur_ncluster_;
    core::Real threshold_;
    void copy_coord(core::Size, FA2d &, FA2d &);
};


/// @brief Greedy K-Center Clustering Algorithm
/// @note "A Fast Geometric Clustering Method on Conformation Space of Biomolecules"
/// Jian Sun, Yuan Yao, Xuhui Huang, Vijay Pande, Gunnar Carlsson, Leonidas J. Guibas
class GreedyKCenter : public KCluster
{
public:
    GreedyKCenter();

    bool whoami();
    core::Real get_threshold();
    void init(KClusterElementOP, core::Size first=0);
    core::Real assign(KClusterElementOP, KClusterData&);
    void update(KClusterElementOP, KClusterData&);
    void set_threshold(core::Real);
protected:
    core::Real radius_;
};

//get the type of cluster from options
KClusterOP get_K_cluster_engine(const std::string&);
//parse the tag, get full path of the silent file
std::string file_full_path(std::string);
//fix the missing tag
std::string fix_tag_suffix(std::string);

}//toolbox
}//protocols

#endif

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/toolbox/KCluster.hh>

#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <utility/exit.hh>

#include <iostream>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::toolbox;

static THREAD_LOCAL basic::Tracer TR( "test_cluster" );

void show_cluster_assignment(KClusterElementOP elem, KClusterData dat){
	for(core::Size ii = 1; ii <= elem->get_ndata(); ii++ ){
		TR << ii  << " ";
		TR << "tag is " << dat.get_tag(ii) << " ";
		//std::cerr << " " << elem->get_type(ii) << " ";
		//std::cerr << " center:ndx: " << elem->get_center_ndx(elem->get_type(ii)) << " ";
		TR << "cluster-tag is: " << dat.get_tag(elem->get_center_ndx(elem->get_type(ii))) << " ";
		TR << elem->get_distance(ii) <<std::endl;
	}
	utility::vector1<Size> cluster_ndx = elem->get_ndx_list();


}

void save_cluster_tree_lite(KClusterElementOP elem, KClusterData dat) {
  utility::vector1<string> fnlist;
  if (option[in::file::silent].user())
    {
      fnlist=option[in::file::silent]();
    }
  else if (option[in::file::silent_list].user())
    {
      fnlist=option[in::file::silent_list]();
    }
  utility::vector1<std::string> tags_to_get;
  utility::vector1<std::string> source_files;

  for( Size ii = 1; ii <= elem->ncluster();ii++ ) {
    TR << "cluster " << ii << " has tag: " << dat.get_tag(elem->get_center_ndx(ii)) << " and came from file: " << dat.source_filename( elem->get_center_ndx(ii)) << std::endl;
    tags_to_get.push_back(dat.get_tag(elem->get_center_ndx(ii)));
    source_files.push_back( dat.source_filename( elem->get_center_ndx(ii)) );
  }

  //this will be ok as long as we assume the number of clusters is much less than the amount of data we have to work with.
  core::io::silent::SilentFileData sfd;
  utility::vector1<std::string> tags;
  utility::vector1< std::string > uniq_source_files;
  Size count_tags_found = 0;

  for( Size jj = 1; jj <= source_files.size(); jj++ ) {
    if( std::find( uniq_source_files.begin(), uniq_source_files.end(), source_files[ jj ] ) == uniq_source_files.end() ) {
      uniq_source_files.push_back( source_files[ jj ] );
    }
  }
  TR << " number of unique source files: " << uniq_source_files.size() << std::endl;
  for( Size zz = 1; zz <= uniq_source_files.size(); zz++ ) {
    TR << "one of uniq source files: " << uniq_source_files[ zz ] << std::endl;
  }

  for( Size jj = 1; jj <= uniq_source_files.size(); jj++ ) {
    tags = sfd.read_tags_fast( uniq_source_files[jj] );
    utility::vector1<std::string> tags_to_ext;
    for( Size kk = 1; kk <= tags_to_get.size(); kk++ ) {
      if( std::find( tags.begin(), tags.end(), tags_to_get[kk]) != tags.end() ) {
	TR << "yay! tag: " << tags_to_get[kk] << " was found in file:  " << uniq_source_files[jj] << std::endl;
	count_tags_found++;
	tags_to_ext.push_back(tags_to_get[kk]);
      }
    }
    sfd.read_file( uniq_source_files[ jj ], tags_to_ext );
    TR << "tags: " << count_tags_found << " of expected " << elem->ncluster() << std::endl;
    //sfd.write_all( option[out::file::silent](), false);
    for( core::Size nn = 1; nn <= tags_to_ext.size(); nn++ ) {
      core::io::silent::SilentStructOP ss = sfd[ tags_to_ext[ nn ] ];
      sfd.write_silent_struct( *ss, option[ out::file::silent ](), false );
    }
    sfd.clear();
  }

}

void set_up_engine(KClusterOP &engine, Size ndx)
{
	Size style_ndx = ndx>option[cluster::K_style]().size() ? 1 : ndx; //default: the first style
    engine = get_K_cluster_engine(option[cluster::K_style]()[style_ndx]);
    if (option[cluster::K_n_cluster].user())
	{
        if (ndx>option[cluster::K_n_cluster]().size()) //ncluster should be specified by user
        {
            utility_exit_with_message("Undefined N_cluster for this level!");
        }
        engine->set_ncluster(option[cluster::K_n_cluster]()[ndx]);
    }

	if (option[cluster::K_style]()[ndx] == "GKC" && option[cluster::K_radius].user())
	{
    	if (ndx>option[cluster::K_radius]().size()) //radius should be specified by user
    	{
    		utility_exit_with_message("Undefined Radius for this level!");
    	}
		else
		{
			engine->set_threshold(option[cluster::K_radius]()[ndx]);
		}
	}
}

bool do_clustering( KClusterData &data, KClusterElementOP element, Size level, Size last, Size firstcenter=0)
{
    if (element->get_ndata()<=0) return true; //which means there is only one structure in the cluster and it was already saved in higher level
	KClusterOP engine;
	set_up_engine(engine, level);
	//std::cout << "Debug: level=" << level << std::endl;
	engine->cluster( element, data, firstcenter );

	if (last>0)
	{
		Size nsub = element->get_cur_ncluster();
		for (Size i=1; i<=nsub; i++)
		{
			Size first = 0;
			if ( option[cluster::K_redundant] )
			{
				first = element->get_center_ndx(i);
			}
			do_clustering( data, element->get_subcluster(i), level+1, last-1, first);
		}
	}
	return true;
}

int main(int argc, char *argv[])
{

	try {

	devel::init(argc, argv);

	KClusterData dat;
	KClusterElementOP elem = new KClusterElement(dat.get_ndata());

/* old loop version
    //setup engine
	KClusterOP engine;
	set_up_engine(engine, 1);

	//root
	engine->cluster(elem, dat);
	if (option[cluster::K_level]()>1)
	{
	    set_up_engine(engine, 2);
		Size nsub = elem->get_cur_ncluster();
		for (Size i=1; i<=nsub; i++)
		{
			engine->cluster( elem->get_subcluster(i), dat );
			if (option[cluster::K_level]()>2)
			{
			    set_up_engine(engine, 3);
			    Size nsubsub = elem->get_subcluster(i)->get_cur_ncluster();
				for (Size j=1; j<=nsubsub; j++)
				{
					engine->cluster( elem->get_subcluster(i)->get_subcluster(j), dat );
				}
			}
		}
	}
*/

	//recursive version
	do_clustering(dat, elem, 1, option[cluster::K_level]()-1);

	dat.mark_tags(elem, "c_1");//a silly hack, avoiding the prefix being parse by the parser, and it would be changed back later, by fix_tag_suffix
	//dat.save_all_in_one();

	//SHOULD BE
	//delete all data before read the silent file
	if( option[ out::file::silent ].user() && !option[cluster::lite].user() ) {
	  dat.save_cluster_tree();
	} else if ( option[ cluster::lite ].user() && option[ out::file::silent].user() ) {
	  std::cout << "attempting to store cluster-centers in a lite-weight manner!" << std::endl;
	  save_cluster_tree_lite(elem,dat);
	}else {
	  std::cout << " not saving clusters in silent-file, but i will output the tags so you can figure out the cluster-assignments! " << std::endl;
	  std::cout << " if you want a silent-file at the end, remember to specify -out::file::silent " << std::endl;
	}
	show_cluster_assignment( elem, dat);
	std::cout << "done with clustering! you can now use your clusters!" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


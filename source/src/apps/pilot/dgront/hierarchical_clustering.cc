// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


//#include <core/types.hh>
#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>

#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

// Clusterign stuff
#include <numeric/ClusteringTreeNode.hh>
#include <numeric/agglomerative_hierarchical_clustering.hh>

#include <string>

#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

static thread_local basic::Tracer TR( "cluster_anything" );

OPT_1GRP_KEY( Integer, clustering, n )
OPT_1GRP_KEY( String, clustering, distances )
OPT_1GRP_KEY( Boolean, clustering, single )
OPT_1GRP_KEY( Boolean, clustering, complete )
OPT_1GRP_KEY( Boolean, clustering, average )

//using namespace core;
using namespace numeric;

void register_options() {

    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    NEW_OPT( clustering::n, "set the desired number of clusters",5);
    NEW_OPT( clustering::single, "set single  linkage as the clustering method", "" );
    NEW_OPT( clustering::complete, "set complete  linkage as the clustering method", "" );
    NEW_OPT( clustering::average, "set average linkage as the clustering method","");
    NEW_OPT( clustering::distances, "input file with all pairwise distances", "" );
}

void read_distances(std::istream & input,utility::vector1<utility::vector1<numeric::Real> > & distances) {

    utility::vector1<Size> id1;
    utility::vector1<Size> id2;
    utility::vector1<Real> d;
    Size max_id = 0;
    std::string line;
    while( getline( input, line ) ) {
	if ( line.substr(0,1) == "#" ) continue;
    	std::istringstream line_stream( line );

	Size resi;
	Size resj;
	Real dist;
	line_stream  >> resi >> resj >> dist;
	if(max_id<resi) max_id = resi;
	if(max_id<resj) max_id = resj;
	id1.push_back( resi );
	id2.push_back( resj );
	d.push_back( dist );
    }
    TR << "Got "<<id1.size()<<" data lines, max_id = "<<max_id<<std::endl;
    distances.clear();
    distances.resize(max_id);
    for(Size i=1;i<=max_id;i++)
	distances[i].resize( max_id );

    for(Size i=1;i<=id1.size();i++) {
	distances[ id2[i] ][ id1[i] ] = d[i];
	distances[ id1[i] ][ id2[i] ] = d[i];
    }
}


int main( int argc, char * argv [] ) {
    try {
    using namespace numeric;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    register_options();
    devel::init(argc, argv);

//------------- set up distance matrix -----------------
    utility::vector1<utility::vector1<Real> > dm;
    if(!option[clustering::distances].user())
	utility_exit_with_message("You must give a file that defines distances between clustered objects !");
    read_distances(utility::io::izstream( option[ clustering::distances ]() ), dm );
    utility::vector1<Size> data_in;
    for(Size i=1;i<=dm.size();i++)
	data_in.push_back(i);

    Size n = option[clustering::n]();
    utility::vector1<ClusteringTreeNodeOP> roots;
    if (option[clustering::single].user()) {
		SingleLinkClusterer slc;
		roots = slc.cluster(dm, n);
    }
	if (option[clustering::complete].user()) {
		CompleteLinkClusterer clc;
		roots = clc.cluster(dm, n);
    }
	if (option[clustering::average].user()) {
		AverageLinkClusterer alc;
		roots = alc.cluster(dm, n);
    }

    std::cout<<"# Here are the "<<n<<" clusters you asked for:\n";
    std::cout<<"#   id  size  distance : members\n";
    for(Size i=1;i<=roots.size();i++) {
	std::cout<<std::setw(6)<<roots[i]->id()<<std::setw(6)<<roots[i]->size()<<" "
		<<std::setw(8)<<roots[i]->distance()<< " : ";
	utility::vector1<Size> data_out;
	get_cluster_data(data_in,roots[i],data_out);
	for(Size j=1;j<=data_out.size();j++)
	    std::cout<<" "<<std::setw(5)<<data_out[j];
	std::cout<<"\n";
    }
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }

    return 0;
}



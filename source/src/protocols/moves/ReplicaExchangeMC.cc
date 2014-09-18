// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/ReplicaExchangeMC.cc
/// @brief  implementing a Recplica Exchange Monte Carlo Mover
/// @author Yuan Liu (wendao@u.washington.edu

#ifdef USEMPI
#include <mpi.h>
#endif //USEMPI

#include <cmath>
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>
#include <protocols/moves/ReplicaExchangeMC.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.ReplicaExchangeMC" );

namespace protocols {
namespace moves {

ReplicaExchangeMC::ReplicaExchangeMC(
    Pose const & init_pose, // PoseCOP init_pose,
    ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
    utility::vector1<core::Real> const &tlist,
    core::Size nint
):MonteCarlo(init_pose, scorefxn, 0.0),
rank_(0),
size_(1),
nreplica_frequency_(nint),
//noutput_(200),
ntrials_(0),
Tlist_(tlist),
last_energylist(NULL),
T_tag(NULL),
T_rev(NULL),
T_ndx(0)
{
    init();
}

ReplicaExchangeMC::ReplicaExchangeMC(
    ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
    utility::vector1<core::Real> const &tlist,
    core::Size nint
):MonteCarlo(scorefxn, 0.0),
rank_(0),
size_(1),
nreplica_frequency_(nint),
//noutput_(200),
ntrials_(0),
Tlist_(tlist),
last_energylist(NULL),
T_tag(NULL),
T_rev(NULL),
T_ndx(0)
{
    init();
}

void ReplicaExchangeMC::init()
{
#ifdef USEMPI
    //init mpi parameter
    MPI_Comm_rank( MPI_COMM_WORLD, &rank_ );
    MPI_Comm_size( MPI_COMM_WORLD, &size_ );
#endif

    //make sure the tlist match the nproc
		if( size_ != (int)Tlist_.size() ){
			std::cerr << "size " << size_ << " and Tlist-size: " << Tlist_.size() << " don't match up!" << std::endl;
		}
    runtime_assert(size_ == static_cast<int>(Tlist_.size()));

    //setup init list
    last_energylist = new double[size_];
    T_tag = new int[size_];
    T_rev = new int[size_];

    for (int i=0; i<size_; i++)
    {
        T_tag[i] = i;
        T_rev[i] = i;
    }

    //set temperature
    T_ndx = rank_;
    set_temperature(Tlist_[T_tag[T_ndx]+1]);

    TR << "I am proc:" << rank_ << ", of all " << size_ << " procs!";
    TR << " T=" << Tlist_[T_tag[T_ndx]+1] <<  " ninterval=" << nreplica_frequency_ << std::endl;

    if (rank_>0) return;

    //setup exchange_schedule
    utility::vector1<std::pair<int, int> > list;
    for (int i=0; i<size_-1; i+=2)
    {
        std::pair<int, int> elem(i, i+1);
        list.push_back(elem);
    }
    exchange_schedule.push_back(list);
    list.clear();
    for (int i=1; i<size_-1; i+=2)
    {
        std::pair<int, int> elem(i, i+1);
        list.push_back(elem);
    }
    exchange_schedule.push_back(list);

    //debug
    TR << "Building exchange schedule 1" << std::endl;
    for (Size i=1; i<=exchange_schedule[1].size(); i++)
    {
        TR << exchange_schedule[1][i].first << "<-->" << exchange_schedule[1][i].second << std::endl;
    }
    TR << "Building exchange schedule 2" << std::endl;
    for (Size i=1; i<=exchange_schedule[2].size(); i++)
    {
        TR << exchange_schedule[2][i].first << "<-->" << exchange_schedule[2][i].second << std::endl;
    }
}

ReplicaExchangeMC::~ReplicaExchangeMC()
{
    if (last_energylist!=NULL) delete [] last_energylist;
    if (T_tag!=NULL) delete [] T_tag;
    if (T_rev!=NULL) delete [] T_rev;
}

void ReplicaExchangeMC::build_temperature_list(double *elist)
{
    static int flag=-1;
    int nlist=(3+flag)/2; // 1 or 2, switch

    TR << "Building " << nlist << std::endl;
    for(Size i=1; i<=exchange_schedule[nlist].size(); i++)
    {
        Size node1=T_rev[exchange_schedule[nlist][i].first];
        Size node2=T_rev[exchange_schedule[nlist][i].second];

        //TR << node1 << "<==>" << node2 << std::endl;
        //TR << "proc" << node1 << ": e=" << elist[node1] <<" T=" << Tlist_[T_tag[node1]+1] << std::endl;
        //TR << "proc" << node2 << ": e=" << elist[node2] <<" T=" << Tlist_[T_tag[node2]+1] << std::endl;

        Real delta=(1.0/Tlist_[T_tag[node1]+1]-1.0/Tlist_[T_tag[node2]+1])*(elist[node2]-elist[node1]);
        Real probability = 0;
        //if (delta>0) probability = 1.0;
        //else probability = std::exp( std::max(-40.0,delta) );
        if (delta<0) probability = 1.0;
        else probability = std::exp( std::max(-40.0, -delta) );
        TR << "Try:" << Tlist_[exchange_schedule[nlist][i].first+1]
					 << " <==> " << Tlist_[exchange_schedule[nlist][i].second+1]
					 << std::endl;
				//TR << "Delta=" << delta << " Prob=" << probability << std::endl;

        if (numeric::random::rg().uniform()<probability)
        {
						TR << "Switch:" << Tlist_[exchange_schedule[nlist][i].first+1] << "<==>" << Tlist_[exchange_schedule[nlist][i].second+1] << std::endl;
            //switch
            Size tmp=T_tag[node1];
            T_tag[node1]=T_tag[node2];
            T_tag[node2]=tmp;
            tmp=T_rev[exchange_schedule[nlist][i].first];
            T_rev[exchange_schedule[nlist][i].first]=T_rev[exchange_schedule[nlist][i].second];
            T_rev[exchange_schedule[nlist][i].second]=tmp;
        }
    }
    TR << "Done" << std::endl;

    //switch
    flag = -flag;
}

bool
ReplicaExchangeMC::boltzmann(
        Pose & pose,//PoseOP pose,
        std::string const & move_type,
        core::Real const proposal_density_ratio,
				core::Real const inner_score_temperature_delta)
{
    ntrials_++;

#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
    if (ntrials_ % (nreplica_frequency_/2) == 0)
    {
        double last_energy_ = last_accepted_score();

        //get infomation
        MPI_Gather(&last_energy_, 1, MPI_DOUBLE, last_energylist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //change the T_tag and T_rev at node0
        if (rank_==0) build_temperature_list(last_energylist);
        //public the new T_tag
        MPI_Scatter(T_tag, 1, MPI_INT, &T_ndx, 1, MPI_INT, 0, MPI_COMM_WORLD);
        set_temperature(Tlist_[T_ndx+1]);
    }
#endif

    //if (ntrials_ % noutput_ == 0 ) {
        //this line should be moved to another condition after debug, and controld by another option
        //TR << "proc=" << rank_  << " step=" << ntrials_ << " T=" << temperature() << " E=" << last_accepted_score() << std::endl;
    //}

    return Parent::boltzmann( pose, move_type, proposal_density_ratio, inner_score_temperature_delta );
}

} // moves
} // prot



// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/NV/NVscore.cc
/// @details  Implementation of Neighbor Vector estimation ( from Durham EA, Et. Al. "Solvent Accessible Surface Area Approximations for Protein Structure Prediction"
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#include <utility/vector1.hh>
#include <core/scoring/nv/NVscore.hh>
#include <core/scoring/nv/NVscoreCreator.hh>

// AUTO-REMOVED #include <core/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <numeric/constants.hh>

#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/EnergyMap.hh>


namespace core {
namespace scoring {
namespace nv {


/// @details This must return a fresh instance of the NVscore class,
/// never an instance already in use
methods::EnergyMethodOP
NVscoreCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new NVscore;
}

ScoreTypes
NVscoreCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( neigh_count );
	sts.push_back( neigh_vect );
	sts.push_back( neigh_vect_raw );
	return sts;
}


static basic::Tracer TR("core.scoring.NVscore");


NVscore::NVscore() :
	parent( new NVscoreCreator ),
	lookup_table_(ScoringManager::get_instance()->get_NVLookupTable() )
{}


methods::EnergyMethodOP NVscore::clone() const
{
	return new NVscore(*this);
}


void NVscore::setup_for_scoring(pose::Pose &pose, const ScoreFunction &) const
{
	pose.update_residue_neighbors();
}

void NVscore::setup_for_packing(pose::Pose &pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const
{
	pose.update_residue_neighbors();
}

void NVscore::setup_for_derivatives(pose::Pose &pose, const ScoreFunction &) const
{
	pose.update_residue_neighbors();
}

void NVscore::setup_for_minimizing(pose::Pose & pose, ScoreFunction const & ,kinematics::MinimizerMapBase const &) const
{
	pose.update_residue_neighbors();
}



void NVscore::indicate_required_context_graphs(utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[twelve_A_neighbor_graph] = true;
}

///Calculate the weighted neighbor count given an upper and lower bound
Real NVscore::neighborWeight(Vector::Value  &dist, Real &lBound, Real &uBound) const
{

	if(dist <= lBound)
	{
		//neighbor count score is 1 if less than the lower bound
		return(1);
	}else if(dist >= uBound)
	{
		//neighbor count score is 0 if gerater than upper bound
		return(0);
	}else if( (lBound < dist) && (uBound > dist) )
	{
		//if between upper and lower bound, score follows a smooth function

		Real weight = ( cos( ( (dist-lBound) / (uBound-lBound) ) * numeric::constants::r::pi ) + 1 )/2.0;
		return(weight);
	}
	return(0);
}

void NVscore::residue_energy(const conformation::Residue &rsd, const pose::Pose &pose, EnergyMap & emap) const
{

	//lbound defaults to 3.3 and ubound defaults to 11.1.  If you change these values the lookup table may no longer be accurate
	Real lBound = basic::options::option[ basic::options::OptionKeys::score::NV_lbound]();
	Real uBound = basic::options::option[ basic::options::OptionKeys::score::NV_ubound]();

	Real neighborCount(0);
	Vector neighborVectSum(0,0,0);

	conformation::ResidueOPs::iterator poseIT;
	//use the coordinates of residue neighbor atom for all calcuations
	Vector currentVect(rsd.nbr_atom_xyz());

	//pose::Pose pose(pose);
	//iterate through the pose
	for(core::Size pose_index = 1; pose_index <= pose.total_residue() ; ++pose_index)
	{
		//get the residue to compare to the current ersidue rsd
		conformation::Residue compRes(pose.residue(pose_index));
		//you don't want to compare a residue to itself
		if(rsd.seqpos() == compRes.seqpos()) continue;
		Vector compVect(compRes.nbr_atom_xyz());
		//calculate the distance between the two residues
		Vector::Value dist = currentVect.distance(compVect);
		//get the weighted neighbor count
		Real weight = neighborWeight(dist,lBound, uBound);

		//calculate the weighted neighbor vector for this pair and sum
		Vector weightedVector = ( (compVect-currentVect) / dist) * weight;
		neighborCount += weight;
		neighborVectSum += weightedVector;

	}
	if ( neighborCount == 0.0 ) return; // do not try to divide by zero

	Vector avgSum = neighborVectSum/neighborCount;
	//neighbor vector score is the norm of the average sum of all neighbor vectors
	Real neighborVector = avgSum.norm();

	char single_letter_name = rsd.name1();
	//use the neighbor vector score to look up the potential from the knowledge base
	Real NVPotential = lookup_table_.get_potentials(single_letter_name,neighborVector);

	emap[ neigh_vect ] += NVPotential;
	emap[ neigh_vect_raw ] += neighborVector;
	emap[ neigh_count ] += neighborCount;

}
core::Size
NVscore::version() const
{
	return 1; // Initial versioning
}

} //NV
} //scoring
} //core

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/PoseMetricCalculators/InterfaceVectorDefinitionCalculator.hh
/// @brief  Calculates the residues at an interface between two protein chains.  The calculation is done in the following manner.  First the point graph is used to find all residues within some big cutoff of residues on the other chain.  For these residues near the interface, two metrics are used to decide if they are actually possible interface residues.  The first metric is to itterate through all the side chain atoms in the residue of interest and check to see if their distance is less than the nearby atom cutoff, if so then they are an interface residue.  If a residue does not pass that check, then two vectors are drawn, a CA-CB vector and a vector from CB to a CB atom on the neighboring chain.  The dot product between these two vectors is then found and if the angle between them is less than some cutoff then they are classified as interface.
/// @author Ben Stranges (stranges@unc.edu)



#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_InterfaceVectorDefinitionCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_InterfaceVectorDefinitionCalculator_hh
#include <protocols/toolbox/pose_metric_calculators/InterfaceDefinitionCalculatorBase.hh>
#include <protocols/toolbox/pose_metric_calculators/InterfaceVectorDefinitionCalculator.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
//#include <numeric/xyzVector.fwd.hh>
#include <numeric/HomogeneousTransform.hh>
#include <set>

//Auto Headers
#include <utility/vector1_bool.hh>



namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

class InterfaceVectorDefinitionCalculator : public InterfaceDefinitionCalculator {

public:

	typedef std::pair< std::set<core::Size>,std::set<core::Size> > InterfacePair;
	typedef numeric::HomogeneousTransform< core::Real > HTReal;

	//chain number definition
	InterfaceVectorDefinitionCalculator( core::Size const chain1_number, core::Size const chain2_number );
	//full constructor that takes all of the inputs
	InterfaceVectorDefinitionCalculator( core::Size const chain1_number, core::Size const chain2_number,
																			 core::Real CB_dist_cutoff, core::Real nearby_atom_cutoff,
																			 core::Real vector_angle_cutoff, core::Real vector_dist_cutoff );
	//chain character definition
	InterfaceVectorDefinitionCalculator( char const chain1_letter, char const chain2_letter );
	//full constructor that takes all of the inputs and uses chain characters
	InterfaceVectorDefinitionCalculator( char const chain1_letter, char const chain2_letter,
																			 core::Real CB_dist_cutoff, core::Real nearby_atom_cutoff,
																			 core::Real vector_angle_cutoff, core::Real vector_dist_cutoff );

	virtual core::pose::metrics::PoseMetricCalculatorOP clone() const;

	//setters to make this awesome!
	void set_CB_dist_cutoff( core::Real CB_dist_cutoff){
		CB_dist_cutoff_ = CB_dist_cutoff;
		//make function to empty the data here
	}
	void set_nearby_atom_cutoff(core::Real nearby_atom_cutoff){
		nearby_atom_cutoff_ = nearby_atom_cutoff;
		//make function to empty the data here
	}
	void set_vector_angle_cutoff(core::Real vector_angle_cutoff){
		vector_angle_cutoff_ = vector_angle_cutoff;
	}
	void set_vector_distance_cutoff(core::Real vector_dist_cutoff){
		vector_dist_cutoff_ = vector_dist_cutoff;
	}

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:
	///@brief looks at the big set and figures out what is actually pointing towards the interface
	void find_interface_pointing_residues_from_neighbs(core::pose::Pose pose, InterfacePair interface_pairs );
	///@brief find nearby atoms to other in interface
	bool any_atoms_within_cutoff(core::conformation::Residue & res1,
															 core::conformation::Residue & res2,
															 core::Real & cutoff);
	///@brief neighbors to look for vectors within (big set here)
	InterfacePair find_neighbors_within_CB_cutoff( core::pose::Pose pose, core::Real big_cutoff );
	///@brief the Cbeta vector(s) from on rsd to another
	numeric::xyzVector<core::Real> cbeta_vector( core::conformation::Residue & res);
	///@brief the action coordinate for each residue
	numeric::xyzVector<core::Real> select_coord_for_residue(core::conformation::Residue & res);
	///@brief out if res1 and res2 are pointing at eachother
	bool res1_pointed_at_res2( core::conformation::Residue & res1,
														 core::conformation::Residue & res2,
														 core::Real angle_cutoff /*degrees*/,
														 core::Real dist_cutoff);

	//set by the calculator as metrics
	std::set< core::Size > interface_residues_;
	std::set< core::Size > chain1_interface_residues_;
	std::set< core::Size > chain2_interface_residues_;
	core::Size num_interface_residues_;
	core::Size num_chain1_interface_residues_;
	core::Size num_chain2_interface_residues_;
	//cutoffs for various restrictions
	core::Real CB_dist_cutoff_; //distance for big CB cutoff
	core::Real nearby_atom_cutoff_; // used for finding atoms that are close
	core::Real vector_angle_cutoff_; // used for cutoff for res1 CB to res2 CB angle cutoff
	core::Real vector_dist_cutoff_; // used for cutoff for res1 CB to res2 CB vector distance cutoff

};

} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif

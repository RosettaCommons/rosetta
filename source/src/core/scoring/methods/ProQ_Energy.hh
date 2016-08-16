// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ProQ_Energy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Björn Wallner


#ifndef INCLUDED_core_scoring_methods_ProQ_Energy_hh
#define INCLUDED_core_scoring_methods_ProQ_Energy_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ProQPotential.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/dssp/Dssp.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

#include <ObjexxFCL/FArray2D.hh>
// Utility headers


namespace core {
namespace scoring {
namespace methods {


class ProQ_Energy : public WholeStructureEnergy  {
public:
	typedef WholeStructureEnergy  parent;
public:


	ProQ_Energy();

	ProQ_Energy( ProQ_Energy const & src);

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & scorefxn ) const;

	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const & scorefxn,
		EnergyMap & totals
	) const;
	virtual
	core::Size version() const;


	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {};


private:
	ProQPotential const & potential_;
	MembraneTopology topology_;

	ObjexxFCL::FArray2D< Real > prob_profile_;
	ObjexxFCL::FArray2D< Real > scaled_logodds_profile_;
	ObjexxFCL::FArray1D< Real > entropy_;
	Size nres_;
	ObjexxFCL::FArray2D< Real > ss_pred_;
	ObjexxFCL::FArray1D< char > ss1_;
	ObjexxFCL::FArray1D< Real > z_pred_;
	ObjexxFCL::FArray1D< Real > rsa_pred_; //Only used for ProQM where the prediction is a real value
	//ObjexxFCL::FArray1D< char > rsa_class_; //This is from mpSA currently not used...
	ObjexxFCL::FArray1D< char > rsa_class_pred_; //Only used for ProQ2 where the prediction is buried/exposed

	bool all_inputs_ProQM_;
	bool all_inputs_ProQ2_;


	void initialize();
	void output_local_prediction(pose::Pose & pose, ObjexxFCL::FArray1D< Real > & proq,std::string name) const;
	void calculate_feature_vector(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & feature_vector) const;
	void calculate_feature_vector_proq2(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & feature_vector) const;
	Size read_profiles_and_entropy(std::string profile, std::string mtxfile);
	void read_zpred(std::string zpredfile);
	void read_mpSA(std::string mpSAfile);
	void read_acc(std::string accfile);
	void read_ss2(std::string ss2file);

	//these needs to run on the whole in one shot
	void atom_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,Size index,int windowsize=21) const;
	void res_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,Size index,int windowsize=21) const;
	void surf_feature(pose::Pose & pose, utility::vector1< Real> & rsd_sasa_rel,ObjexxFCL::FArray2D< Real > & vec,Size index,int windowsize=21) const;
	void gss_sc_feature(pose::Pose & pose, dssp::Dssp & ss,ObjexxFCL::FArray2D< Real > & vec,Size index) const;
	void grsa_sc_feature(pose::Pose & pose, utility::vector1< Real> & rsd_sasa_rel,ObjexxFCL::FArray2D< Real > & vec,Size index) const;

	//these are called for each residue
	void stride_feature(pose::Pose & pose, dssp::Dssp & ss,ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index, int windowsize=11) const;
	void ss_feature(pose::Pose & pose, dssp::Dssp & ss,ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index,int windowsize=1) const;
	void ss_sc_feature(pose::Pose & pose, dssp::Dssp & ss,ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index,int windowsize=21) const;
	void topology_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index) const;
	void zpred_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index) const;
	void z_feature(pose::Pose & pose, ObjexxFCL::FArray1D< Real > & Z,ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index) const;
	void entropy_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index,int windowsize=11) const;
	void profile_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index,int windowsize=3) const;
	void rsa_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,utility::vector1< Real> & rsd_sasa_rel,int const pos,Size index,int windowsize=1) const;
	void prsa_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index,int windowsize=1) const;
	void rsa_sc_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,utility::vector1< Real> & rsd_sasa_rel,int const pos,Size index,int windowsize=21) const;
	void termini_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,int const pos,Size index,int windowsize=23) const;

	//closest heavy-atom distance
	Real crd(pose::Pose & pose,Size i,Size j) const;
	int res6(conformation::Residue const & rsd) const;
	int res6(char const aa) const;
	int atom13(conformation::Residue const & rsd, Size atom_i) const;
	int atom13_0(conformation::Residue const & rsd, Size atom_i) const;
	char profile_index_to_aa(int i) const;
	void sum_profile(Size k, ObjexxFCL::FArray1D< Real > & vec) const;
	void calculateZ(pose::Pose & pose,ObjexxFCL::FArray1D< Real > & Z) const;
};


}
}
}

#endif // INCLUDED_core_scoring_methods_ProQ_Energy_HH

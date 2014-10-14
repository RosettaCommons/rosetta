// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, may 2009


#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_MatchConstraintFileInfo_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_MatchConstraintFileInfo_hh

// Unit headers
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.fwd.hh>

// Package headers
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.fwd.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.fwd.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// numeric headers
#include <numeric/HomogeneousTransform.fwd.hh>

// C++ headers
#include <map>
#include <string>
#include <set>
#include <list>

#include <utility/vector1.hh>


namespace protocols{
namespace toolbox{
namespace match_enzdes_util{

/// @brief function to go through a list of restypes and
/// reduce them to chemically identical ones based on the same base_name
/// i.e. this function gets rid of the variant redundancy
void
add_relevant_restypes_to_subset(
	utility::vector1< core::chemical::ResidueTypeCOP > & restype_subset,
	utility::vector1< core::chemical::ResidueTypeCOP > const & restypes,
	core::chemical::ResidueTypeSetCAP restype_set
);


/// @brief class that stores information of one geometric parameter line of the .cst file
/// i.e. angleA or torsionB
class GeomSampleInfo : public utility::pointer::ReferenceCount
{

public:

	GeomSampleInfo(
		std::string tag );

	GeomSampleInfo(
		core::Real ideal_val,
		core::Real tolerance,
		core::Real force_k,
		core::Real periodicity_
	);

	virtual ~GeomSampleInfo();

	/// @brief data reading routine
	bool
	read_data( std::istringstream & line_stream );


	/// @brief creates an explicit sample vector from the input data
	/// i.e. translates ideal value, tolerance, num_steps and periodicity
	/// into a vector of values
	utility::vector1< core::Real >
	create_sample_vector() const;

	std::string
	tag() const {
		return tag_; }

	core::Real
	ideal_val() const {
		return ideal_val_; }

	core::Real
	tolerance() const {
		return tolerance_; }

	core::Real
	periodicity() const {
		return periodicity_; }

	core::Real
	force_const() const {
		return force_const_; }

	core::Size
	num_steps() const {
		return num_steps_; }

	core::Real
	step_size() const {
		return step_size_; }

	std::string
	function_tag() const {
		return function_tag_; }



private:

	std::string tag_, function_tag_;

	// ideal val and the periodicity
	core::Real ideal_val_, tolerance_, periodicity_, force_const_;

	core::Size num_steps_;

	//step_size_ = tolerance / num_steps
	core::Real step_size_;
};


class MatchConstraintFileInfo : public utility::pointer::ReferenceCount
{


public:  //construct / destruct

	typedef numeric::HomogeneousTransform< core::Real > HTReal;

	MatchConstraintFileInfo(
		core::Size index,
		core::chemical::ResidueTypeSetCAP restype_set
	);

	virtual ~MatchConstraintFileInfo();

public:  //atom and residue accessors


	core::Size
	index() const{
		return index_;
	}

	/// @brief all positions where a residue for this geometry
	/// can be placed
	utility::vector1< core::Size > const &
	allowed_seqpos() const {
		return allowed_seqpos_; }


	/// @brief what type of amino acids/ligands make this
	/// constraint
	utility::vector1< std::string > const &
	allowed_res_name3s( core::Size which_cstres ) const {
		return this->enz_cst_template_res( which_cstres )->allowed_res_types(); }

	core::Size
	num_enz_cst_template_res() const{
		return enz_template_res_.size(); }

	/// @brief is this interaction a backbone interaction
	bool
	is_backbone( core::Size which_cstres ) const {
		return this->enz_cst_template_res( which_cstres )->is_backbone(); }


	/// @brief all chemically non-redundant restypes of the given restypes
	utility::vector1< core::chemical::ResidueTypeCOP > const &
	allowed_restypes( core::Size which_cstres ) const;

	/// @brief which one of the residues (1 or 2 ) in this block
	/// is the upstream res. used for classic match algorithm
	/// hardcoded 2 for now
	core::Size
	upstream_res() const {
		return 2; }

	/// @brief which one of the residues (1 or 2 ) in this block
	/// is the upstream res. used for classic match algorithm
	/// hardcoded 1 for now
	core::Size
	downstream_res() const {
		return 1; }

	/// @brief holds information read from ALGORITHM_INFO blocks
	/// in the input file
	std::map< std::string, utility::vector1< std::string > > const &
	algorithm_inputs() const {
		return algorithm_inputs_; }

	bool
	is_covalent() const {
		return is_covalent_; }

	GeomSampleInfoCOP
	dis_U1D1() const {
		return dis_U1D1_; }

	GeomSampleInfoCOP
	ang_U1D2() const {
		return ang_U1D2_; }

	GeomSampleInfoCOP
	ang_U2D1() const {
		return ang_U2D1_; }

	GeomSampleInfoCOP
	tor_U1D3() const {
		return tor_U1D3_; }

	GeomSampleInfoCOP
	tor_U3D1() const {
		return tor_U3D1_; }

	GeomSampleInfoCOP
	tor_U2D2() const {
		return tor_U2D2_; }


	/// @brief all atoms of restype to be used as template_atom
	/// in the matcher/constraints
	utility::vector1< core::Size > const &
	template_atom_inds(
		core::Size which_cstres,
		core::Size which_template_atom,
		core::chemical::ResidueType const & restype ) const;

	EnzCstTemplateResCOP
	enz_cst_template_res( core::Size template_res ) const;

	//Kui 110609 Native
	bool native() const{
		return native_;
	}

public: //geometric sample accessors

	/// @brief returns ExternalGeomSampler only if the user has specified all six degrees of freedom,
	/// otherwise null pointer is returned
	ExternalGeomSamplerOP
	create_exgs() const;

public: //mutators

	/// @brief data reading routine
	bool
	read_data( utility::io::izstream & data );


	/// @brief processes the read data
	/// right now this only generates the template atomnos for every restype
	void
	process_data();

public: //convenience functions

	/// @brief function that takes all rotamers for the ResidueType(s)
	/// that interact with the target residues and places them according
	/// according to the geometry specified
	std::list< core::conformation::ResidueCOP >
	inverse_rotamers_against_residue(
		core::Size const target_template,
		core::conformation::ResidueCOP target_conf ) const;

	std::list< core::conformation::ResidueCOP >
	inverse_rotamers_against_residue(
		core::conformation::Residue const & target_conf,
		core::chemical::ResidueTypeCOP invrot_restype,
		utility::vector1< core::Size > const & target_ats,
		utility::vector1< core::Size > const & invrot_ats,
		bool const flip_exgs_upstream_downstream_samples,
		bool const backbone_interaction
	) const;

	void
	diversify_backbone_only_rotamers(
		utility::vector1< core::conformation::ResidueCOP > & rotamers
	) const;


protected:

	/// @brief reads and stores arbitrary algorithm specific input
	bool
	process_algorithm_info(
		std::string tag,
		utility::io::izstream & data
	);

//data
private:

	//the index of this mcfi in the mcfi list
	core::Size index_;

	//allowed positions
	utility::vector1< core::Size > allowed_seqpos_;

	std::map< core::Size, EnzCstTemplateResOP > enz_template_res_;

	bool is_covalent_;

	// the GeomSampleInfos read from the file
	GeomSampleInfoOP dis_U1D1_, ang_U1D2_, ang_U2D1_, tor_U1D3_, tor_U3D1_, tor_U2D2_;

	//container for arbritrary algorithm specific information
	std::map< std::string, utility::vector1< std::string > > algorithm_inputs_;

	core::chemical::ResidueTypeSetCAP restype_set_;

	//Kui Native 110809
	bool native_;
};


/// @brief a simple container class to contain several MatchConstraintFileInfo
/// instances. this can also query the MatchConstraintFileInfos for common upstream
/// restypes and put all their geomsamples into one list
class MatchConstraintFileInfoList : public utility::pointer::ReferenceCount
{
public:

	MatchConstraintFileInfoList(
		core::chemical::ResidueTypeSetCAP restype_set );

	virtual ~MatchConstraintFileInfoList();

public: //accessors

	utility::vector1< core::chemical::ResidueTypeCOP > const &
	upstream_restypes() const {
		return upstream_restypes_; }

	utility::vector1< MatchConstraintFileInfoCOP > const &
	mcfis_for_upstream_restype( core::chemical::ResidueTypeCOP restype ) const;

	//also functions for downstream builders and launch points

	MatchConstraintFileInfoCOP
	mcfi( core::Size which_mcfi ) const {
		return mcfis_[ which_mcfi ];
	}

	core::Size
	num_mcfis() const {
		return mcfis_.size(); }


public: //mutators

	/// @brief data reading routine
	bool
	read_data( utility::io::izstream & data );

public: //convenience functions

	/// @brief function that takes all rotamers for the ResidueType(s)
	/// that interact with the target residues and places them according
	/// according to the geometry specified
	std::list< core::conformation::ResidueCOP >
	inverse_rotamers_against_residue(
		core::Size const target_template,
		core::conformation::ResidueCOP target_conf ) const;

protected:


	/// @brief goes through all the mcfis and figures out the upstream restypes
	void
	determine_upstream_restypes();

private:

	//core::Size active_mcfi_;

	utility::vector1< MatchConstraintFileInfoOP > mcfis_;
	utility::vector1< core::chemical::ResidueTypeCOP > upstream_restypes_;
	std::map< core::chemical::ResidueTypeCOP, utility::vector1< MatchConstraintFileInfoCOP > > mcfis_for_restype_;
	core::chemical::ResidueTypeSetCAP restype_set_;

};

}
} //namespace enzdes
} //namespace protocols




#endif //

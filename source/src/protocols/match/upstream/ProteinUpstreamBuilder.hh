// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/upstream/ProteinUpstreamBuilder.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_upstream_ProteinUpstreamBuilder_hh
#define INCLUDED_protocols_match_upstream_ProteinUpstreamBuilder_hh

// Unit headers
#include <protocols/match/upstream/ProteinUpstreamBuilder.fwd.hh>

// Package headers
#include <protocols/match/downstream/DownstreamAlgorithm.fwd.hh>
//#include <protocols/match/downstream/ExternalGeomSampler.fwd.hh>

#include <protocols/match/upstream/ProteinSCSampler.fwd.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.fwd.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>


// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/types.hh>

// Numeric headers
#include <numeric/HomogeneousTransform.hh>
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace upstream {

/// @brief This class holds all of the data associated with the logic
/// for generating extra samples for a particular chi angle.  There are
/// tons of ways concievable to build extra rotamers; the data in this class
/// is intended to group all of that data into one place.  This class is
/// not responsible for building extra rotamer samples; that responsibility
/// is given to class FullChiSampleSet.
class SampleStrategyData {
public:
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pack::task::ExtraRotSample ExtraRotSample;

public:
	SampleStrategyData();
	~SampleStrategyData();

	void set_strategy( ChiStrategy );

	void set_sample_level( ExtraRotSample setting );
	void set_step_size( Real setting );
	void set_sd_range(  Real setting );
	void set_n_samples_wi_sd_range( Size setting );
	void set_nrchi_prob_minimum_for_extra_samples( Real setting );
	void set_n_samples_per_side_of_nrchi_bin( Size setting );

	ChiStrategy strategy() const;
	ExtraRotSample sample_level() const;
	Real step_size() const;
	Real sd_range() const;
	Size n_samples_wi_sd_range() const;
	Real nrchi_prob_minimum_for_extra_samples() const;
	Size n_samples_per_side_of_nrchi_bin() const;


private:
	ChiStrategy strategy_;
	ExtraRotSample sample_level_;
	Real step_size_;
	Real sd_range_;
	Size n_samples_wi_sd_range_;
	Real nrchi_prob_minimum_for_extra_samples_;
	Size n_samples_per_side_of_nrchi_bin_;
};

/// @brief A simple class that describes the geometry for a particular
/// residue type.  It describes the coordinate frame geometry for the
/// fourth atom defining each chi dihedral. The fourth atom is called
/// the "chi tip" atom, as it's at the tip of the growing kinematic chain
/// when building chi i. This class also describes the location
/// of the atoms controlled by each chi which are not the chi-tip atoms;
/// it measures their location in the coordinate frame of the chi-tip atom.
///
/// @details To generate the coordinate of the chi-tip atom, the
/// stored coordinate frame is multiplied
/// by the coordinate frame at the third atom after that coordinate frame
/// has been multipled by the chi-angle-z-axis rotation HT.
/*
class UpstreamResTypeGeometry : public utility::pointer::ReferenceCount
{
public:
typedef core::Size                            Size;
typedef core::Real                            Real;
typedef core::Vector                          Vector;
typedef numeric::HomogeneousTransform< Real > HTReal;

public:
UpstreamResTypeGeometry();
UpstreamResTypeGeometry( core::chemical::ResidueType const & );

void initialize_from_residue_type( core::chemical::ResidueType const & );

public:

/// @brief the name of the residue type used to generate this geometry
std::string const & name() const {
return restype_name_;
}

/// @brief the number of atoms in this residue type
Size natoms() const {
return controlling_chi_for_atom_.size();
}

Size nchi() const {
return chitip_atoms_.size();
}

bool atom_controlled_by_any_chi( Size atomno ) const {
return controlling_chi_for_atom_[ atomno ] != 0;
}

bool atom_is_chitip( Size atomno ) const {
return controlling_chi_for_atom_[ atomno ] != 0 && which_point_for_atom_[ atomno ] == 0;
}

utility::vector1< Size > const &
controlling_chi_for_atom() const { return controlling_chi_for_atom_; }

utility::vector1< Size > const &
which_point_for_atom() const { return which_point_for_atom_; }

utility::vector1< Size > const &
chitip_atoms() const { return chitip_atoms_; }

Size
chitip_atom( Size chi ) const {
return chitip_atoms_[ chi ];
}

utility::vector1< HTReal > const &
ht_for_chitip_atoms() const { return ht_for_chitip_atoms_; }

HTReal const &
ht_for_chitip_atom( Size chi ) const {
return ht_for_chitip_atoms_[ chi ];
}

Size
n_nonchitip_atoms_for_chi( Size chi ) const {
return nonchitip_atoms_[ chi ].size();
}

utility::vector1< utility::vector1< Size > > const &
nonchitip_atoms() const { return nonchitip_atoms_; }

Size
nonchitip_atom( Size chi, Size which_nonchitip_atom_for_chi ) const {
return nonchitip_atoms_[ chi ][ which_nonchitip_atom_for_chi ];
}

utility::vector1< utility::vector1< Vector > > const &
points_for_nonchitip_atoms() const { return points_for_nonchitip_atoms_; }

utility::vector1< Vector > const &
points_for_nonchitip_atoms( Size chi ) const {
return points_for_nonchitip_atoms_[ chi ];
}

/// @brief Convenience function: get the coordinate in the chitip frame
/// for a particular atom.  The atom must be a non-chitip atom that is
/// not part of the backbone (it must be controlled by a chi angle).
Vector const &
point_for_nonchitip_atom( Size atom ) {
assert( atom_controlled_by_any_chi( atom ) && !atom_is_chitip( atom ) );
return points_for_nonchitip_atoms_[ controlling_chi_for_atom_[ atom ] ]
[ which_point_for_atom_[ atom ] ];
}


private:
/// Data

std::string restype_name_;

utility::vector1< Size >   controlling_chi_for_atom_;
utility::vector1< Size >   which_point_for_atom_;

utility::vector1< Size >   chitip_atoms_;
utility::vector1< HTReal > ht_for_chitip_atoms_;


utility::vector1< utility::vector1< Size > >   nonchitip_atoms_;
utility::vector1< utility::vector1< Vector > > points_for_nonchitip_atoms_;

};
*/

/// @brief Still sketchy on this class.  It holds all of the data needed
/// for describing the geometry of the downstream partner relative to the upstream
/// partner for a single residue type.  This class holds the data for deciding which
/// base rotamers to consider, how each base rotamer should be expanded to produce a
/// full set of chi rotamers, and how to orient the downstream partner relative to
/// this rotamer.  It also holds the UpstreamResTypeGeometry object for the
/// restype being built.
class BuildSet : public utility::pointer::ReferenceCount
{
public:
	typedef core::Size Size;
	typedef core::Real Real;
	typedef utility::pointer::ReferenceCount parent;

public:

	BuildSet();
	virtual ~BuildSet();
	BuildSet( BuildSet const & );

	BuildSet const & operator = ( BuildSet const & rhs );

public:
	/// initialization
	void set_residue_type( core::chemical::ResidueTypeCOP restype, bool backbone_only = false );

	void set_sample_strategy_for_chi( Size chi, SampleStrategyData const & data );

	void set_downstream_algorithm( downstream::DownstreamAlgorithmOP );

public:
	/// accessors

	bool
	has_restype() const {
		return restype_.get() != 0;
	}

	core::chemical::ResidueType const &
	restype() const {
		return *restype_;
	}

	bool
	backbone_only() const {
		return backbone_only_;
	}

	UpstreamResTypeGeometry const &
	restype_geometry() const {
		return *restype_geom_;
	}

	Real
	probability_limit() const {
		return rot_prob_accumulation_limit_;
	}


	SampleStrategyData const &
	sample_strategy_for_chi( Size chi ) const {
		return sample_strategy_for_chi_[ chi ];
	}

	bool
	has_algorithm() {
		return algorithm_.get() != 0;
	}

	downstream::DownstreamAlgorithm const &
	algorithm() const {
		return *algorithm_;
	}

	downstream::DownstreamAlgorithm &
	algorithm() {
		return *algorithm_;
	}

	Size
	nbonds_from_bb_atom( Size atom_index ) const {
		return nbonds_from_bb_atom_[ atom_index ];
	}

	ProbeRadius
	atom_radius( Size atomno ) const {
		return atom_radii_[ atomno ];
	}

	void
	set_fa_dun_cutoff(
		core::Real cutoff );

	bool
	check_fa_dun() const {
		return check_fa_dun_; }

	core::Real
	fa_dun_cutoff() const{
		return fa_dun_cutoff_; }

public:
	core::chemical::ResidueTypeCOP  restype_;
	bool backbone_only_; // true if only the geometry of the backbone is important

	UpstreamResTypeGeometryOP restype_geom_;
	Real rot_prob_accumulation_limit_;
	utility::vector1< SampleStrategyData >  sample_strategy_for_chi_;

	utility::vector1< Size > nbonds_from_bb_atom_;
	utility::vector1< ProbeRadius >  atom_radii_;

	downstream::DownstreamAlgorithmOP algorithm_;

	//in case we only want to build hits from rotamers below a certain fa_dun
	bool check_fa_dun_;
	core::Real fa_dun_cutoff_;

};

class FullChiSampleSet : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FullChiSampleSet();
	typedef core::Size Size;
	typedef core::Real Real;
	typedef numeric::HomogeneousTransform< Real > HTReal;
	typedef core::pack::task::ExtraRotSample ExtraRotSample;
	typedef core::pack::dunbrack::DunbrackRotamerSampleData DunbrackRotamerSampleData;
	typedef utility::vector1< DunbrackRotamerSampleData > DunbrackRotamerSampleDataVector;

public:

	FullChiSampleSet(
		BuildSet const & build_set,
		core::pack::dunbrack::DunbrackRotamerSampleData const & sample,
		bool dry_run
	);


	utility::vector1< Size >
	n_samples_per_chi() const {
		assert( ! dry_run_ );
		return n_samples_per_chi_;
	}

	Size num_chi_samples_total() const {
		return num_chi_samples_total_;
	}

	utility::vector1< Real > const &
	chi_samples( Size chi ) const {
		assert( ! dry_run_ );
		return chi_samples_[ chi ];
	}


	Real
	chi_sample( Size chi, Size sample_id ) const {
		assert( ! dry_run_ );
		return chi_samples_[ chi ][ sample_id ];
	}

	HTReal const &
	frame( Size chi, Size sample_id ) const {
		assert( ! dry_run_ );
		assert( chi <= frames_.size() );
		return frames_[ chi ][ sample_id ];
	}

private:

	void
	expand_non_dunbrack_chi( Size chi, BuildSet const & build_set );

	void
	expand_samples_by_ex_behavior(
		Size chi,
		ExtraRotSample behavior,
		core::pack::dunbrack::DunbrackRotamerSampleData const & sample
	);

	void
	expand_samples_by_steps_wi_sdrange(
		Size chi,
		SampleStrategyData const & stratdat,
		core::pack::dunbrack::DunbrackRotamerSampleData const & sample
	);

	void
	expand_samples_for_nrchi_wi_nrchi_bin(
		Size chi,
		SampleStrategyData const & stratdat,
		core::pack::dunbrack::DunbrackRotamerSampleData const & sample
	);

	void
	create_hts_for_chi( Size chi );

	/// This doesn't belong in this class -- move to core.
	static
	ExtraRotSample
	ex_level_from_flags( Size chi );

private:
	bool const dry_run_;
	Size num_chi_samples_total_;
	utility::vector1< Size > n_samples_per_chi_;
	utility::vector1< utility::vector1< Real > > chi_samples_;
	utility::vector1< utility::vector1< HTReal > > frames_;

};


class ProteinUpstreamBuilder : public UpstreamBuilder {
public:
	typedef core::Vector Vector;
	typedef core::Real   Real;
	typedef numeric::HomogeneousTransform< Real > HTReal;
	typedef core::pack::dunbrack::DunbrackRotamerSampleData DunbrackRotamerSampleData;
	typedef utility::vector1< DunbrackRotamerSampleData > DunbrackRotamerSampleDataVector;

public:
	ProteinUpstreamBuilder();
	virtual ~ProteinUpstreamBuilder();

	UpstreamBuilderOP
	clone() const;

	/// @brief Iterate across possible conformations for the upstream
	/// half of the hit, and for each (non-coliding) conformation,
	/// sample all external geometries specified by the external_sampler
	/// to construct the three coordinates of the downstream sampler.
	/// Return a list of hits.
	virtual
	std::list< Hit >
	build(
		ScaffoldBuildPoint const & build_point
	) const;

	/// @brief Regenerate the rotamer for a particular hit and give that rotamer
	/// to the UpstreamResidueProcessor.
	virtual
	void
	recover_hit(
		Hit const & hit,
		ScaffoldBuildPoint const & build_point,
		UpstreamResidueProcessor & processor
	) const;

	/// @brief Regenerate a set of rotamers for a subset of hits bound by the
	/// two input hit-list iterators.
	virtual
	void
	recover_hits(
		std::list< Hit >::const_iterator hits_begin,
		std::list< Hit >::const_iterator hits_end,
		ScaffoldBuildPoint const & build_point,
		UpstreamResidueProcessor & processor
	) const;

	virtual
	Size
	n_restypes_to_build() const;

	virtual
	core::chemical::ResidueTypeCOP
	restype( Size which_restype ) const;

	virtual bool compatible(
		Hit const & my_hit,
		ScaffoldBuildPoint const & build_point_mine,
		UpstreamBuilder const & other,
		Hit const & other_hit,
		ScaffoldBuildPoint const & build_point_other,
		bool first_dispatch = true
	) const;

	virtual bool compatible(
		Hit const & my_hit,
		ScaffoldBuildPoint const & build_point_mine,
		ProteinUpstreamBuilder const & other,
		Hit const & other_hit,
		ScaffoldBuildPoint const & build_point_other,
		bool first_dispatch = true
	) const;

	void
	add_build_set(
		BuildSet const & build_set
	);

	Size
	n_build_sets() const {
		return build_sets_.size();
	}

	BuildSet const &
	build_set( core::chemical::ResidueTypeCOP restype ) const;

	BuildSet &
	build_set( core::chemical::ResidueTypeCOP restype );

	void
	set_sampler(
		ProteinSCSamplerCOP sampler
	);

	void set_use_input_sidechain( bool setting );

	//Kui Native 110809
	void set_native_flag (bool native);

private:

	/// @brief Copy the coordinates from the build_point object into the the rescoords
	/// object, compute the coordinate frame at CBeta, copy the CB coordinate into the
	/// rescoords object, and return the CBeta frame.
	HTReal
	initialize_rescoords(
		Size build_set_id,
		core::conformation::Residue & rescoords,
		ScaffoldBuildPoint const & build_point
	) const;

	/// @brief Construct the coordinate frame at CBeta to match the residue type's geometry.
	/// The input coordinates for the build point need not be ideal.
	/// The returned coordinate frame is located at CBeta with the z axis pointing
	/// along the calpha->cbeta bond vector, the y axis is in the N-CA-CB plane,
	/// and the x axis being the crossproduct of y and z.
	HTReal
	compute_cb_frame(
		Size build_set_id,
		ProteinBackboneBuildPoint const & build_point
	) const;

	/// @brief Returns true if a particular atom coordinate for a rotamer rules out
	/// that rotamer as yeilding a potential hit.  Reasons to exclude a rotamer based
	/// on atom placement: collision of that atom with the background, too long of a
	/// distance that atom is from any other hit.
	bool
	atom_coordinate_unacceptable(
		Size build_set_id,
		core::conformation::Residue const & rescoords,
		Size atomno
	) const;

private:
	ProteinSCSamplerCOP sampler_;
	utility::vector1< BuildSet > build_sets_;
	bool use_input_sc_;
	//Kui Native 110809
	bool avoid_building_any_rotamers_dueto_native_;

};

}
}
}

#endif

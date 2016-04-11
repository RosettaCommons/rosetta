// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/simple_moves/bb_sampler/BBDihedralSampler.hh
/// @brief This class functions to hold, access, and set independent and dependent dihedral data.
///   It can act as a base class for particular types of data.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)  and Jason W. Labonte (JWLabonte@jhu.edu)


#ifndef INCLUDED_protocols_simple_moves_bb_sampler_BBDihedralSampler_hh
#define INCLUDED_protocols_simple_moves_bb_sampler_BBDihedralSampler_hh

#include <protocols/simple_moves/bb_sampler/BBDihedralSampler.fwd.hh>

#include <core/id/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace simple_moves {
namespace bb_sampler {

enum BBSampleType {
	minima = 1,
	probability
};


///@brief This class functions to hold, access, sample on, and set independent and dependent dihedral data.
///   It is an abstract base class for particular types of data.
///
///   get_torsion functions should give either the minima on a set of data or sample via the probability
///   If you are subclassing, you do not nessessarily need to use the BBSampleType.
///
///    Feel free to implement more types. See the SugarBBSampler and RangedBBSampler as examples.
///
///   Technically - with now storing the torsion_types as Size (due to waay to many in sugars) - we can now remove the 'BB' part of this whole thing.
///
class BBDihedralSamplerBase : public utility::pointer::ReferenceCount {

public:

	BBDihedralSamplerBase();
	BBDihedralSamplerBase(core::Size torsion_type, BBSampleType sampling_type = probability);

	BBDihedralSamplerBase(BBDihedralSamplerBase const & src);

	virtual ~BBDihedralSamplerBase();

	BBDihedralSamplerBaseOP
	clone() const;


public:

	///@brief Set the torsion type we will be querying.
	void
	set_torsion_type( core::Size torsion_type ) {
		torsion_type_ = torsion_type;
	};

	core::Size
	get_torsion_type( ) const {
		return torsion_type_;
	};

	///@brief Set the sampling type.  Default is to sample probabilistically on the data.
	void
	set_sample_type(BBSampleType sampling_type) {
		sampling_type_ = sampling_type;
	};

	virtual std::string
	name() const = 0;


protected:

	core::Size torsion_type_;
	BBSampleType sampling_type_;


};


///@brief This class functions to hold, access, sample on, and set independent and dependent dihedral data.
///   It is an abstract base class for particular types of data.
///   It should eventually be moved out of here.
///
///   get_torsion functions should give either the minima on a set of data or sample via the probability.
class BBDihedralSampler : public BBDihedralSamplerBase {

public:

	BBDihedralSampler();
	BBDihedralSampler( core::Size torsion_type, BBSampleType sampling_type = probability );

	BBDihedralSampler( BBDihedralSampler const & src );

	virtual ~BBDihedralSampler();

	BBDihedralSamplerOP
	clone() const;

public:

	virtual core::Real
	get_torsion(core::pose::Pose const & pose, core::Size resnum) const = 0;

	///@brief Set torsions to pose
	virtual void
	set_torsion_to_pose(core::pose::Pose & pose, core::Size resnum) const = 0;

	std::string name() const {
		return "BBDihedralSampler";
	};
};



class BBDihedralSampler2D : public BBDihedralSamplerBase {

public:
	BBDihedralSampler2D();
	BBDihedralSampler2D( core::Size torsion_type, BBSampleType sampling_type = probability );

	BBDihedralSampler2D(BBDihedralSampler2D const & src);

	virtual ~BBDihedralSampler2D();

	BBDihedralSampler2DOP
	clone() const;

public:


	///@brief Get a torsion angle dependant on another torsion and torsion angle.
	virtual core::Real
	get_2d_torsion(core::pose::Pose const & pose, core::Size resnum,
		std::pair<core::id::MainchainTorsionType, core::Real > ) const = 0;

	virtual void
	set_2D_torsion_to_pose(core::pose::Pose & pose, core::Size resnum,
		std::pair<core::id::MainchainTorsionType, core::Real > ) const = 0;

	std::string name() const {
		return "BBDihedralSampler2D";
	};

protected:



};


class BBDihedralSampler3D : public BBDihedralSamplerBase {

public:
	BBDihedralSampler3D();
	BBDihedralSampler3D( core::Size torsion_type, BBSampleType sampling_type = probability );

	BBDihedralSampler3D(BBDihedralSampler3D const & src);

	virtual ~BBDihedralSampler3D();

	BBDihedralSampler3DOP
	clone() const;

public:

	///@brief Get a torsion angle dependant on two other torsions and torsion angles.
	/// dependendant types and angles are the std::pair.
	virtual core::Real
	get_3d_torsion(core::pose::Pose const & pose, core::Size resnum,
		std::pair< core::id::MainchainTorsionType, core::Real >,
		std::pair< core::id::MainchainTorsionType, core::Real > ) const = 0;

	virtual void
	set_3D_torsion_to_pose(core::pose::Pose & pose, core::Size resnum,
		std::pair< core::id::MainchainTorsionType, core::Real >,
		std::pair< core::id::MainchainTorsionType, core::Real > ) const = 0;

	std::string name() const {
		return "BBDihedralSampler3D";
	};

};


class BBDihedralSamplerND : public BBDihedralSamplerBase {

public:
	BBDihedralSamplerND();
	BBDihedralSamplerND( core::Size torsion_type, BBSampleType sampling_type = probability );
	BBDihedralSamplerND(BBDihedralSamplerND const & src);

	virtual ~BBDihedralSamplerND();

	BBDihedralSamplerNDOP
	clone() const;

public:

	///@brief Get a torsion angle dependant on n number of other torsion types and angles.
	virtual core::Real
	get_ND_torsion(
		core::pose::Pose const & pose, core::Size resnum,
		utility::vector1< std::pair< core::id::MainchainTorsionType, core::Real > >) const = 0;

	virtual void
	set_ND_torsion_to_pose(
		core::pose::Pose & pose, core::Size resnum,
		utility::vector1< std::pair< core::id::MainchainTorsionType, core::Real > >) const = 0;

	std::string name() const {
		return "BBDihedralSamplerND";
	};

};



} //protocols
} //pose
} //core


#endif //INCLUDED_core_pose_BBDihedralSampler_hh






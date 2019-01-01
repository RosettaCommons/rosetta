// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pre/PREMover.hh
/// @brief   Mover that converts NMR PRE rates into CB-CB atom pair distances and
///          assigns them as atom pair constraints to the pose. In case of degenerate spins
///          (e.g. degenerate protons or symmetric spins) an AmbiguousNMRDistanceConstraint
///          is created.
///          As constraints function a SplineFunc is fitted to a PRE distance histogram.
///          By default, the histogram file is read from the spinlabel database and the generated
///          spline potential converts the HN PRE distances into CB-CB distance constraints.
///          Alternatively, the user can provide a different histogram file to generate CB-CB distance
///          constraints from PRE rates collected for a different type of PRE nucleus e.g. 15N or 1HA
///          or a different spinlabel type.
/// @details last Modified: 12/16/16
///          Note that the calculation of PRE distances is simplified and takes into account only
///          the dipolar and in part Curie relaxation so far. Cross-correlated relaxation is neglected.
///          This approach is reasonable for radicals or paramagnetic ions with a nearly isotropic g-tensor
///          (e.g. nitroxide, Mn2+, Cu2+) which are commonly used in NMR PRE experiments. However, for ions
///          with an anisotropic g-tensor (e.g. lanthanides) Curie and cross-correlated relaxation become
///          more prominent and a direct conversion of PRE rates to distances is not easily possible any longer.
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_pre_PREMover_HH
#define INCLUDED_protocols_nmr_pre_PREMover_HH

// Unit headers
#include <protocols/nmr/pre/PREMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <core/scoring/nmr/pre/PREData.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/id/AtomID.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>
#include <iostream>
#include <map>
#include <set>

namespace protocols {
namespace nmr {
namespace pre {

class PREDistanceRecord : public utility::pointer::ReferenceCount {

public: // Types

	typedef core::Real Real;
	typedef core::Size Size;

public: // Methods

	inline
	PREDistanceRecord() :
		rsds_(),
		dist_(0),
		tol_(0),
		count_(0)
	{}

	inline
	PREDistanceRecord(
		std::set< Size > const & rsds,
		Real d,
		Real t
	) :
		rsds_(rsds),
		dist_(d),
		tol_(t),
		count_(1)
	{}

	inline
	PREDistanceRecord(
		utility::vector1< Size > const & rsds,
		Real d,
		Real t
	) :
		rsds_(rsds.begin(), rsds.end()),
		dist_(d),
		tol_(t),
		count_(1)
	{}

	inline std::set< Size > const & rsds() const { return rsds_; }
	inline std::set< Size > & rsds() { return rsds_; }
	inline Real get_dist() const { return dist_; }
	inline void set_dist(Real d) { dist_ = d; }
	inline Real get_tol() const { return tol_; }
	inline void set_tol(Real t) { tol_ = t; }
	inline Size count() const { return count_; }

	/// @brief Add a new distance and tolerance/weight and update
	///        the class data to hold the current mean values
	inline
	void
	add_and_average(
		Real d,
		Real t
	)
	{
		Real new_dist = ( (dist_ * count_) + d) / (count_ + 1);
		Real new_tol  = ( (tol_  * count_) + t) / (count_ + 1);
		dist_ = new_dist;
		tol_  = new_tol;
		++count_;
	}

private: // Data

	std::set< Size > rsds_;
	Real dist_;
	Real tol_;
	// iteration count, because we want to use the class
	// also for on the fly mean calculation
	Size count_;
};

class PREMover : public protocols::moves::Mover {

public: // Types

	typedef core::Real              Real;
	typedef core::Size              Size;
	typedef core::scoring::nmr::pre::PREData        PREData;
	typedef core::scoring::nmr::pre::PREDataOP         PREDataOP;
	typedef core::scoring::ScoreFunctionOP          ScoreFunctionOP;
	typedef core::pose::Pose            Pose;
	typedef protocols::moves::MoverOP          MoverOP;
	// spinlabel code to histogram file and histogram bin size
	typedef std::map< std::string, std::pair< std::string, core::Real > > SpinlabelHistogramMap;
	// Vector of PRE distance values and errors/weights
	typedef utility::vector1< PREDistanceRecord >       PREDistances;
	// Map spinlabel residue to its PRE distance vector
	typedef std::map< core::Size, PREDistances >        SpinlabelToPREDistances;
	typedef utility::fixedsizearray1< Real, 8 >        Vec8;

public: // Methods

	/// @brief Default constructor
	PREMover();

	/// @brief Construct PREMover from PRE data input file
	PREMover(
		std::string const & pre_data_file,
		Pose const & pose
	);

	/// @brief Copy constructor
	PREMover(PREMover const & other);

	/// @brief Copy assignment
	PREMover &
	operator=(PREMover const & rhs);

	/// @brief destructor
	~PREMover() override;

	/// @brief Get the name of this mover
	std::string get_name() const override;

	static std::string mover_name();

	/// @brief Make a deep copy of this mover
	MoverOP clone() const override;

	/// @brief Create a fresh instance of this mover
	MoverOP fresh_instance() const override;

	/// @brief Calculate CB-CB distances from PRE rates and append
	///        them as atom pair distance constraints to the pose
	void apply(Pose & pose) override;

	void show(std::ostream& TR) const override;

	/// @brief Parse tags of XML script
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & pose
	) override;

	/// @brief Create XML schema definition for PREMover
	static void provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd);

	/// Getter and Setters
	PREDataOP get_pre_data() { return pre_data_; }
	ScoreFunctionOP get_scorefunction() const { return sfxn_; }
	bool weighted_average() const { return weighted_average_; }
	bool minimize_w_pre_csts() const { return minimize_; }

	void set_pre_data(PREDataOP pre_data) { pre_data_ = pre_data; }
	void set_scorefunction(ScoreFunctionOP sfxn) { sfxn_ = sfxn; }
	void set_weighted_average(bool av) { weighted_average_ = av; }
	void set_minimize_w_pre_csts(bool min) {minimize_ = min; }

	void
	add_histogram_file(
		std::string const & spinlabel_name,
		std::string const & histogram_file,
		Real bin_size = 0.5
	);

private: // Methods

	/// @brief Calculate distance from R2 relaxation rate
	/// @details Considers dipolar and Curie relaxation
	/// @params
	/// params[1] = gamma_I:     gyromagnetic ratio of the nuclear spin (must be provided in rad/(s*T), dimension is 10^6)
	/// params[2] = gJ:          electron Lande factor
	/// params[3] = S:           total spin quantum number
	/// params[4] = omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// params[5] = tau_c:       total correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[6] = tau_r:       rotational correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[7] = B0:          magnetic field strength (in Tesla)
	/// params[8] = T:           temperature (in K)
	/// R2:                      R2 relaxation rate (in Hz)
	Real
	R2_to_dist_dd_curie(
		Vec8 const & params,
		Real const R2
	);

	/// @brief Calculate distance from R1 relaxation rate
	/// @details Considers dipolar and Curie relaxation
	/// @params
	/// params[1] = gamma_I:     gyromagnetic ratio of the nuclear spin (must be provided in rad/(s*T), dimension is 10^6)
	/// params[2] = gJ:          electron Lande factor
	/// params[3] = S:           total spin quantum number
	/// params[4] = omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// params[5] = tau_c:       total correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[6] = tau_r:       rotational correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[7] = B0:          magnetic field strength (in Tesla)
	/// params[8] = T:           temperature (in K)
	/// R1:                      R1 relaxation rate (in Hz)
	Real
	R1_to_dist_dd_curie(
		Vec8 const & params,
		Real const R1
	);

	/// @brief Calculate distances from relaxation rates and map them
	///        to their respective spinlabel and protein residue(s)
	void
	pre_data_to_distances(
		core::scoring::nmr::pre::PREData & pre_data,
		SpinlabelToPREDistances & all_sl_distances
	);

private: // Data

	/// @brief Map of histogram files for very spinlabel type that is
	///        used in pre_data_. The histogram is used for instantiation
	///        of a SplineFunc potential for conversion of the measured
	///        PRE distance (e.g. for HN) into a CB-CB atom pair constraint.
	SpinlabelHistogramMap histogram_files_;

	/// @brief collection of all PRE datasets for multiple spinlabel sites
	PREDataOP pre_data_;

	/// @brief scorefunction object
	ScoreFunctionOP sfxn_;

	/// @brief use PRESingleSet weights to calculate an average distance
	///        in case that the same PRE distance was measured multiple times
	///        (e.g. at different field strengths) for the same spinlabel site
	bool weighted_average_;

	/// @brief Do one round of minimization of input pose after PRE distances
	///        constraints are added to the pose
	bool minimize_;

};

} // namespace pre
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_pre_PREMover_HH

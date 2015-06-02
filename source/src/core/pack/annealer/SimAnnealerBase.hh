// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/SimAnnealerBase.hh
/// @brief  Packer's simulated annealing base class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_pack_annealer_SimAnnealerBase_HH
#define INCLUDED_core_pack_annealer_SimAnnealerBase_HH

// Unit Headers
#include <core/pack/annealer/SimAnnealerBase.fwd.hh>

// Package Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers

// STL Headers

#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>


namespace core {
namespace pack {
namespace annealer {

//extern bool annealing_starts_at_low_temperature;

class SimAnnealerBase : public utility::pointer::ReferenceCount
{
public:
	typedef rotamer_set::RotamerSetsBaseCOP RotamerSetsBaseCOP;

public:

	SimAnnealerBase(
		int num_rots_to_pack,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	virtual ~SimAnnealerBase();

	virtual void run() = 0;

	bool pass_metropolis( core::PackerEnergy delta_energy) const;
	bool pass_metropolis( core::PackerEnergy previous_energy, core::PackerEnergy delta_energy ) const;

	void setup_iterations();
	void setup_iterations( const int & num_of_state_changes );
	int get_outeriterations() const;
	int get_inneriterations() const;
	void setup_temperature( const int & nn );
	void setup_temperature( const ObjexxFCL::FArray1D< core::PackerEnergy >& loopenergy,int nn );
	void set_temperature( core::PackerEnergy new_temperature );
	core::PackerEnergy get_temperature() const;
	void set_to_quench();
	void set_not_to_quench();
	bool quench() const;
	bool get_start_with_current() const;
	bool get_calc_rot_freq() const;
	void set_disallow_quench( bool const & setting );

	void set_hightemp( core::PackerEnergy );
	void set_lowtemp( core::PackerEnergy );

	inline void scale_outeriterations( core::PackerEnergy const so )
	{
			outeriterations_scaling_ = so; return;
	}

	inline void scale_inneriterations( core::PackerEnergy const si )
	{
			inneriterations_scaling_ = si; return;
	}

protected:

	static const core::PackerEnergy hightemp;
	static const core::PackerEnergy lowtemp;
	static const int maxouteriterations = 500;
	static const core::PackerEnergy calc_freq_temp;

	Size num_rots_to_pack() const { return num_rots_to_pack_; }
	void num_rots_to_pack( Size setting );

	ObjexxFCL::FArray1D_int& bestrotamer_at_seqpos();
	ObjexxFCL::FArray1D_int const & bestrotamer_at_seqpos() const;
	core::PackerEnergy & bestenergy();
	bool start_with_current() const;
	ObjexxFCL::FArray1_int & current_rot_index();
	ObjexxFCL::FArray1_int const & current_rot_index() const;
	bool calc_rot_freq() const;
	ObjexxFCL::FArray1D< core::PackerEnergy >& rot_freq();
	ObjexxFCL::FArray1D< core::PackerEnergy > const & rot_freq() const;

	core::PackerEnergy get_hightemp() const { return hightemp_; }
	core::PackerEnergy get_lowtemp() const { return lowtemp_; }

	void clear(); // resets counts modified by get_temperature

private:

	Size num_rots_to_pack_;
	ObjexxFCL::FArray1D_int& bestrotamer_at_seqpos_;
	core::PackerEnergy & bestenergy_;
	bool start_with_current_;
	ObjexxFCL::FArray1_int & current_rot_index_; //assert current_rot_index.size() == pose.total_residue()
	bool calc_rot_freq_;
	ObjexxFCL::FArray1D< core::PackerEnergy >& rot_freq_;


	int outeriterations_;
	int inneriterations_;
	bool quench_;
	core::PackerEnergy hightemp_; // initialized at instantiation
	core::PackerEnergy lowtemp_; // initialized at instantiation
	core::PackerEnergy temperature_;
	int jump_;

	core::PackerEnergy outeriterations_scaling_;
	core::PackerEnergy inneriterations_scaling_;

	bool const low_temp_annealing_;
	bool disallow_quench_;

	SimAnnealerBase(const SimAnnealerBase& rhs);
};

} //namespace annealer
} //namespace pack
} //namespace core

#endif

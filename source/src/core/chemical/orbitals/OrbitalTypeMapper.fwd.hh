// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


#ifndef INCLUDED_core_chemical_orbitals_OrbitalTypeMapper_fwd_hh
#define INCLUDED_core_chemical_orbitals_OrbitalTypeMapper_fwd_hh

namespace core {
namespace chemical {
namespace orbitals {
class OrbitalTypeMapper;

enum orbital_type_enum{
	C_pi_sp2=1,
	N_pi_sp2,
	N_p_sp2,
	O_pi_sp2,
	O_p_sp2,
	O_p_sp3,
	S_p_sp3,
	O_pi_sp2_bb,
	O_p_sp2_bb,
	num_orbital_types= O_p_sp2_bb
};


}
}
}

#endif /* ORBITALTYPEMAPPER_FWD_HH_ */

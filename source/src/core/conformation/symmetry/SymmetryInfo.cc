// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Ingemar Andre, Phil Bradley

// Unit Headers
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Package headers
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/id/types.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>

#include <core/conformation/Conformation.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.functions.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// C++ headers
#include <iostream>

// core utilities
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "core.conformation.SymmetryInfo" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace symmetry {


/// @details  helper function for symminfo serialization
template<class S, class T>
void
comma_strings_to_vector_map(
	std::istream & is,
	Size const nbase,
	std::map< S, utility::vector1< T > > & clones,
	std::string const & tag=""
)
{
	bool fail( false );
	std::string tag0;
	if ( tag != "" ) {
		is >> tag0;
		if ( tag0 != tag ) {
			TR << "Input failed: tag mismatch " << tag << " " << tag0 << std::endl;
			return;
		}
	}

	for ( Size i=1; !fail && i<= nbase; ++i ) {
		std::string jump_string;
		is >> jump_string;
		if ( is.fail() ) {
			fail = true;
			break;
		}
		std::replace( jump_string.begin(), jump_string.end(), ',', ' ' );
		std::istringstream l( jump_string );
		S base_jump;
		l >> base_jump;
		if ( l.fail() ) {
			fail = true;
			break;
		}
		while ( true ) {
			T j;
			l >> j;
			if ( l.fail() ) break;
			clones[ base_jump ].push_back( j );
		}
		if ( clones[ base_jump ].size() < 1 ) {
			fail = true;
			break;
		}
	}
	if ( clones.size() != nbase ) {
		is.setstate( std::ios_base::failbit );
	}
}


/// @details  helper function for symminfo serialization
template<class S, class T, class U>
void
comma_strings_to_pair_map(
	std::istream & is,
	Size const nbase,
	std::map< S, std::pair<T,U> > & clones,
	std::string const & tag=""
) {
	bool fail( false );
	std::string tag0;
	if ( tag != "" ) {
		is >> tag0;
		if ( tag0 != tag ) {
			TR << "Input failed: tag mismatch " << tag << " " << tag0 << std::endl;
			return;
		}
	}
	for ( Size i=1; !fail && i<= nbase; ++i ) {
		std::string jump_string;
		is >> jump_string;
		if ( is.fail() ) {
			fail = true;
			break;
		}
		std::replace( jump_string.begin(), jump_string.end(), ',', ' ' );
		std::istringstream l( jump_string );
		S base_jump;
		l >> base_jump;
		if ( l.fail() ) {
			fail = true;
			break;
		}
		T pair1;
		U pair2;
		l >> pair1 >> pair2;
		clones[ base_jump ] = std::pair<T,U>(pair1, pair2);
	}
	if ( clones.size() != nbase ) {
		is.setstate( std::ios_base::failbit );
	}
}

/// @details  helper function for symminfo serialization
template<class S, class T>
void
comma_strings_to_map(
	std::istream & is,
	Size const nbase,
	std::map< S, T > & clones,
	std::string const & tag=""
) {
	bool fail( false );
	std::string tag0;
	if ( tag != "" ) {
		is >> tag0;
		if ( tag0 != tag ) {
			TR << "Input failed: tag mismatch " << tag << " " << tag0 << std::endl;
			return;
		}
	}
	for ( Size i=1; !fail && i<= nbase; ++i ) {
		std::string jump_string;
		is >> jump_string;
		if ( is.fail() ) {
			fail = true;
			break;
		}
		std::replace( jump_string.begin(), jump_string.end(), ',', ' ' );
		std::istringstream l( jump_string );
		S base_jump;
		l >> base_jump;
		if ( l.fail() ) {
			fail = true;
			break;
		}
		l >> clones[ base_jump ];
	}
	if ( clones.size() != nbase ) {
		is.setstate( std::ios_base::failbit );
	}
}


/// @details  helper function for symminfo serialization
template<class S>
void
comma_strings_to_vector(
	std::istream & is,
	Size const nbase,
	utility::vector1< S > & clones,
	std::string const & tag=""
)
{
	bool fail( false );
	std::string tag0;
	if ( tag != "" ) {
		is >> tag0;
		if ( tag0 != tag ) {
			TR << "Input failed: tag mismatch " << tag << " " << tag0 << std::endl;
			return;
		}
	}

	std::string jump_string;
	is >> jump_string;
	if ( is.fail() ) fail = true;
	std::replace( jump_string.begin(), jump_string.end(), ',', ' ' );
	std::istringstream l( jump_string );
	while ( true ) {
		S j;
		l >> j;
		if ( l.fail() ) break;
		clones.push_back( j );
	}
	if ( clones.size() != nbase ) {
		fail = true;
	}

	if ( fail ) {
		is.setstate( std::ios_base::failbit );
	}
}

void
comma_strings_to_map(
	std::istream & is,
	Size const nbase,
	std::map< Size, SymDof > & clones,
	std::string const & tag=""
)
{
	bool fail( false );
	std::string tag0;
	if ( tag != "" ) {
		is >> tag0;
		if ( tag0 != tag ) {
			TR << "Input failed: tag mismatch " << tag << " " << tag0 << std::endl;
			return;
		}
	}
	for ( Size i=1; !fail && i<= nbase; ++i ) {
		std::string jump_string;
		is >> jump_string;
		if ( is.fail() ) {
			fail = true;
			break;
		}
		std::replace( jump_string.begin(), jump_string.end(), ',', ' ' );
		std::istringstream l( jump_string );
		Size base_jump;
		l >> base_jump;
		std::string dof_line;
		l >> dof_line;
		clones[base_jump].read(dof_line);
		if ( l.fail() ) {
			fail = true;
			break;
		}
	}
	if ( fail ) {
		is.setstate( std::ios_base::failbit );
	}
}

/// @details  helper function for symminfo serialization
template<class S, class T>
void vector_map_to_comma_strings(
	std::ostream & out,
	std::map< S, utility::vector1< T > > clones,
	std::string const & tag=""
) {
	typename std::map< S,utility::vector1<T> >::const_iterator it;
	if ( tag != "" ) out << ' ' << tag ;
	for ( it = clones.begin(); it != clones.end(); ++it ) {
		out << ' ' << it->first;
		utility::vector1< T > const & l( it->second );
		for ( Size i=1; i<= l.size(); ++i ) {
			out << ',' << l[i];
		}
	}
}


/// @details  helper function for symminfo serialization
template<class S, class T>
void map_to_comma_strings(
	std::ostream & out,
	std::map< S, T > clones,
	std::string const & tag=""
) {
	typename std::map< S, T >::const_iterator it;
	if ( tag != "" ) out << ' ' << tag ;
	for ( it = clones.begin(); it != clones.end(); ++it ) {
		out << ' ' << it->first << ',' << it->second ;
	}
}

/// @details  helper function for symminfo serialization
template<class S, class T, class U>
void pair_map_to_comma_strings(
	std::ostream & out,
	std::map< S, std::pair<T,U> > clones,
	std::string const & tag=""
) {
	typename std::map< S, std::pair<T,U> >::const_iterator it;
	if ( tag != "" ) out << ' ' << tag ;
	for ( it = clones.begin(); it != clones.end(); ++it ) {
		out << ' ' << it->first << ',' << it->second.first<< ',' << it->second.second;
	}
}


/// @details  helper function for symminfo serialization
template<class S>
void
vector_to_comma_strings(
	std::ostream & out,
	utility::vector1 < S > clones,
	std::string const & tag=""
)
{
	if ( tag != "" ) out << ' ' << tag ;
	for ( Size i=1; i<= clones.size(); ++i ) {
		if ( i == 1 ) {
			out << " " << clones[i] ;
		} else {
			out << ',' << clones[i] ;
		}
	}

}

SymmetryInfo::SymmetryInfo() {
	init_defaults();
}

SymmetryInfo::~SymmetryInfo() = default;

void
SymmetryInfo::init_defaults() {
	use_symmetry_ = false;
	score_multiply_factor_ = 1.0;
	last_indep_residue_ = 0;
	reweight_symm_interactions_ = basic::options::option[ basic::options::OptionKeys::symmetry::reweight_symm_interactions ]();
	contiguous_monomers_ = true;
	torsion_changes_move_other_monomers_ = false;
	nres_monomer_ = 0;
	scoring_subunit_ = 1;
	npseudo_ = 0;
	njump_monomer_ = 0;
	interfaces_ = 0;
	num_components_ = 1;
}

// fpd this is horribly out of date
//bool SymmetryInfo::operator==( SymmetryInfo const & s )
//{
// return ( npseudo_ == s.npseudo_ &&
//           bb_clones_ ==s.bb_clones_ &&
//           chi_clones_ == s.chi_clones_ &&
//           jump_clones_ == s.jump_clones_ );
//}
//bool SymmetryInfo::operator!=( SymmetryInfo const & s ) {
// return !( *this == s );
//}

SymmetryInfo::SymmetryInfo( SymmData const & symm_data, Size const nres_subunit, Size const njump_subunit )
{
	init_defaults();

	Size joff = njump_subunit*symm_data.get_subunits();
	std::map<std::string,Size> const & name2num = symm_data.get_jump_string_to_jump_num();
	for ( auto const & i : name2num ) {
		dofname2jnum_[i.first] = i.second+joff;
		jnum2dofname_[i.second+joff] = i.first;
	}
	if (  symm_data.get_jump_clones().size() > 0 ) {
		initialize( nres_subunit, njump_subunit,
			symm_data.get_subunits(), symm_data.get_num_virtual(),
			symm_data.get_jump_clones(), symm_data.get_dofs(),
			symm_data.get_score_subunit(), symm_data.get_score_multiply_subunit(),
			symm_data.get_slide_info(), symm_data.get_interfaces() );
	} else {
		initialize( nres_subunit, njump_subunit,
			symm_data.get_subunits(), symm_data.get_num_virtual(),
			symm_data.get_dofs(),symm_data.get_score_subunit(),
			symm_data.get_score_multiply_subunit(), symm_data.get_slide_info(),
			symm_data.get_interfaces() );
	}
	TR.Debug << *this << std::endl;
}

SymmetryInfo::SymmetryInfo(
	Size const nres_monomer,
	Size const njump_monomer,
	Size const N,
	std::map< Size, SymDof > dofs,
	Size const score_subunit,
	utility::vector1< Real > score_multiply_subunit,
	SymSlideInfo slide_info,
	Size const num_interfaces,
	std::string const & type
)
{
	contiguous_monomers_ = true;
	torsion_changes_move_other_monomers_ = false;
	initialize( nres_monomer, njump_monomer, N, N, dofs, score_subunit,
		score_multiply_subunit, slide_info, num_interfaces, type );
}

SymmetryInfoOP
SymmetryInfo::clone() const
{
	return SymmetryInfoOP( new SymmetryInfo( *this ) );
}

// Initialize from "de-novo-style" symmdef files (as from make_symmdef_denovo.py)
void
SymmetryInfo::initialize(
	Size const nres_monomer,
	Size const njump_monomer,
	Size const n_subunits,
	Size const num_virtual,
	std::map< Size, SymDof > dofs,
	Size const score_subunit,
	utility::vector1< Real > score_multiply_subunit,
	SymSlideInfo slide_info,
	Size const num_interfaces,
	std::string const & type
)
{
	init_defaults();

	nres_monomer_ = nres_monomer;

	interfaces_ = num_interfaces;
	set_score_multiply_from_subunit_factors(score_multiply_subunit, nres_monomer, n_subunits);
	njump_monomer_ = njump_monomer;
	dofs_ = dofs;
	use_symmetry_ = true;
	slide_info_ = slide_info;
	type_ = type;
	scoring_subunit_ = score_subunit;

	bb_clones_.clear();
	chi_clones_.clear();
	jump_clones_.clear();

	// check that score_monomer makes sense...
	if ( score_subunit > n_subunits || score_subunit < 1 ) {
		utility_exit_with_message("score_subunit must be in the range 1-N");
	}

	// special case of no symmetry
	//fpd   is this used at all?
	if ( type == "c1" ) {
		npseudo_ = num_virtual;

		for ( Size i=1; i<= nres_monomer; ++i ) {
			Clones clones;
			clones.clear();
			bb_clones_.insert( std::make_pair( i, clones ) );
			chi_clones_.insert( std::make_pair( i, clones ) );
		}
		for ( Size i=1; i<= njump_monomer; ++i ) {
			Clones clones;
			clones.clear();
			jump_clones_.insert( std::make_pair( i, clones ) );
		}
		return;
	} //end c1 symmetry

	for ( Size i=1; i<= nres_monomer; ++i ) {
		Clones clones;
		int base ( i + ( score_subunit - 1 ) * nres_monomer );
		for ( Size k=0; k<n_subunits; ++k ) {
			if ( k+1 != score_subunit ) {
				clones.push_back( i + k * nres_monomer );
				add_bb_clone( base, i + k * nres_monomer );
				add_chi_clone( base, i + k * nres_monomer );
			}
		}
		bb_clones_.insert( std::make_pair( base, clones ) );
		chi_clones_.insert( std::make_pair( base, clones ) );
	}

	// the N*njump_monomer internal jumps
	for ( Size i=1; i<= njump_monomer; ++i ) {
		for ( Size k=1; k<n_subunits; ++k ) {
			add_jump_clone( i, i + k * njump_monomer, 0.0 );
		}
	}

	if ( type == "no_pseudo" ) {
		npseudo_ = num_virtual;
	} else if ( type == "simple" ) {
		// 1                 --> N*njump_monomer  : the internal jumps
		// N*njump_monomer+1 --> N*njump_monomer+N: the pseudo-rsd--monomer jumps
		// last N-1 jumps                         : jumps between pseudo-rsds
		npseudo_ = num_virtual;

		// the N jumps from pseudo-residues to monomers
		{
			Size const base_jump( n_subunits*njump_monomer + 1 );
			//Clones clones;
			for ( Size k=1; k<n_subunits; ++k ) {
				add_jump_clone( base_jump, base_jump + k, 0.0 );
			}
		}

	} else {
		std::cerr << "unrecognized type: " << type << std::endl;
		utility_exit();
	}

	// last_independent_residue
	for ( Size i=1; i <=num_total_residues_without_pseudo(); ++i ) {
		if ( bb_is_independent(i) ) {
			last_indep_residue_ = i;
		}
	}

	update_score_multiply_factor();
}

/// @details  This is a helper function for some of the DOF_ID routines below. Really just a best guess...
/// this is a little tricky: the mapping from a DOF_ID to a TorsionID is not straightforward to
/// construct (see kinematics/util.cc:setup_dof_to_torsion_map)
/// So we don't really know whether a dof_id is a bb degree of freedom or a chi degree of freedom...
/// or even which residue it should be attached to, eg the phi of residue i might be dof_id with rsd i-1
/// So we take a guess based on whether id.atomno is a backbone or sidechain atom
id::TorsionType
guess_torsion_type_of_dof_id( id::DOF_ID const & id, Conformation const & conf )
{
	if ( id::RB1 <= id.type() && id.type() <= id::RB6 ) {
		return id::JUMP;
	} else {
		if ( conf.atom_is_backbone_norefold( id.rsd(), id.atomno() ) ) {
			return id::BB;
		} else {
			return id::CHI;
		}
	}
}

void
SymmetryInfo::update_score_multiply_factor()
{
	// compute the score_multiply_factor
	utility::vector1< bool > indep_res = independent_residues();
	for ( int i=1; i<=(int)indep_res.size(); ++i ) {
		if ( indep_res[i] ) {
			score_multiply_factor_ = score_multiply_[ i ];
			break;
		}
	}
}

//fpd force a default score multiply factor
//    return true if any weights have changed
// reasonable default: main subunit = score_multiply_factor
//    all other subunits = score_multiply_factor/2
//    jump derivative weights = 1
bool
SymmetryInfo::reset_score_multiply_to_reasonable_default() {
	bool any_changes = false;
	for ( Size i = 1; i<= num_total_residues_without_pseudo(); ++i ) {
		Real newval = score_multiply_factor_/2.0;
		if ( bb_is_independent(i) ) newval = score_multiply_factor_;
		if ( score_multiply_[i] != newval ) {
			any_changes = true;
			score_multiply_[i] = newval;
		}
	}
	// virtuals can stay unchanged
	return any_changes;
}

// Initialize from explicit VRT symmdef files (as from make_symmdef_file.pl)
void
SymmetryInfo::initialize(
	Size const nres_monomer,
	Size const njump_monomer,
	Size const n_subunits,
	Size const num_virtual,
	std::map< Size, WtedClones > jump_clones,
	std::map< Size, SymDof > dofs,
	Size const score_subunit,
	utility::vector1< Real > score_multiply_subunit,
	SymSlideInfo slide_info,
	Size const num_interfaces,
	std::string const & type
)
{
	contiguous_monomers_ = true;
	torsion_changes_move_other_monomers_ = false;
	nres_monomer_ = nres_monomer;

	interfaces_ = num_interfaces;
	set_score_multiply_from_subunit_factors(score_multiply_subunit, nres_monomer, n_subunits);
	njump_monomer_ = njump_monomer;
	dofs_ = dofs;
	use_symmetry_ = true;
	slide_info_ = slide_info;
	type_ = type;
	scoring_subunit_ = score_subunit;
	bb_clones_.clear();
	chi_clones_.clear();
	jump_clones_.clear();

	// 1                 --> N*njump_monomer  : the internal jumps
	// N*njump_monomer+1 --> N*njump_monomer+N: the pseudo-rsd--monomer jumps
	// last N-1 jumps                         : jumps between pseudo-rsds
	npseudo_ = num_virtual;

	for ( Size i=1; i<= nres_monomer; ++i ) {
		Clones clones;
		int base ( i + ( score_subunit - 1 ) * nres_monomer );
		for ( Size k=0; k<n_subunits; ++k ) {
			if ( k+1 != score_subunit ) {
				clones.push_back( i + k * nres_monomer );
				add_bb_clone( base, i + k * nres_monomer );
				add_chi_clone( base, i + k * nres_monomer );
			}
		}
		bb_clones_.insert( std::make_pair( base, clones ) );
		chi_clones_.insert( std::make_pair( base, clones ) );
	}

	// the N*njump_monomer internal jumps
	for ( Size i=1; i<= njump_monomer; ++i ) {
		for ( Size k=0; k<n_subunits; ++k ) {
			if ( k != ( score_subunit - 1 ) ) {
				add_jump_clone( i + (score_subunit-1)*njump_monomer, i + k*njump_monomer, 0.0 );
			}
		}
	}

	std::map< Size,WtedClones >::const_iterator it, it_start=jump_clones.begin(), it_end=jump_clones.end();
	for ( it=it_start; it != it_end; ++it ) {
		for ( Size i = 1; i<= it->second.size(); ++i ) {
			add_jump_clone( it->first + n_subunits*njump_monomer, it->second[i].first + n_subunits*njump_monomer, it->second[i].second );
		}
	}

	// last_independent_residue
	for ( Size i=1; i <=num_total_residues_without_pseudo(); ++i ) {
		if ( bb_is_independent(i) ) {
			last_indep_residue_ = i;
		}
	}

	// compute the score_multiply_factor
	update_score_multiply_factor();
}

/////////////////////////////////////////////////////////////////////////////
std::istream& operator>> ( std::istream & s, SymmetryInfo & symminfo )
{

	std::string tag;
	Size num_bb_indep, num_chi_indep, num_jump_indep;
	Size num_bb_dep, num_chi_dep, num_jump_dep;
	Size num_dof, num_score_multiply, nnamed_jumps;

	symminfo.set_use_symmetry(true);

	//fpd -- already set in the constructor
	//symminfo.contiguous_monomers_ = true;
	//symminfo.torsion_changes_move_other_monomers_ = false;

	s >> tag ;
	if ( tag == "SYMMETRY_INFO" && !s.fail() ) {
		bool fail=false,
			nres_monomer_set=false,
			scoring_subunit_set=false,
			njump_monomer_set=false,
			npseudo_set=false,
			nnamed_jumps_set=false,
			interfaces_set=false,
			type_set=false,
			num_bb_indep_set=false,
			num_chi_indep_set=false,
			num_jump_indep_set=false,
			num_bb_dep_set=false,
			num_chi_dep_set=false,
			num_jump_dep_set=false,
			num_dof_set=false,
			num_score_multiply_set=false,
			bb_clones_set=false,
			chi_clones_set=false,
			jump_clones_set=false,
			jump_clone_wts_set=false,
			bb_follows_set=false,
			chi_follows_set=false,
			jump_follows_set=false,
			dofs_set=false,
			score_multiply_set=false,
			num_components_set=false,
			jnum2dofname_set=false,
			dofname2jnum_set=false;
		// unused:
		//bool components_set=false,
		//component_bounds_set=false,
		//name2component_set=false,
		//jname2components_set=false,
		//jname2subunits_set=false;

		//fpd  skip first two fields
		s >> tag >> tag;
		while ( s >> tag ) {
			if ( tag == "N_RES_MONOMER" ) {
				s >> symminfo.nres_monomer_;
				nres_monomer_set=true;
			} else if ( tag == "SCORING_SUBUNIT" ) {
				s >> symminfo.scoring_subunit_;
				scoring_subunit_set=true;
			} else if ( tag == "N_JUMP_MONOMER" ) {
				s >> symminfo.njump_monomer_;
				njump_monomer_set=true;
			} else if ( tag == "N_VIRT" ) {
				s >> symminfo.npseudo_;
				npseudo_set=true;
			} else if ( tag == "N_NAMED_JUMPS" ) {
				s >> nnamed_jumps;
				nnamed_jumps_set=true;
			} else if ( tag == "N_INTERFACE" ) {
				s >> symminfo.interfaces_;
				interfaces_set=true;
			} else if ( tag == "TYPE" ) {
				s >> symminfo.type_;
				type_set=true;
			} else if ( tag == "BB_CLONES_SIZE" ) {
				s >> num_bb_indep;
				num_bb_indep_set=true;
			} else if ( tag == "CHI_CLONES_SIZE" ) {
				s >> num_chi_indep;
				num_chi_indep_set=true;
			} else if ( tag == "JUMP_CLONES_SIZE" ) {
				s >> num_jump_indep;
				num_jump_indep_set=true;
			} else if ( tag == "BB_FOLLOWS_SIZE" ) {
				s >> num_bb_dep;
				num_bb_dep_set=true;
			} else if ( tag == "CHI_FOLLOWS_SIZE" ) {
				s >> num_chi_dep;
				num_chi_dep_set=true;
			} else if ( tag == "JUMP_FOLLOWS_SIZE" ) {
				s >> num_jump_dep;
				num_jump_dep_set=true;
			} else if ( tag == "DOFS_SIZE" ) {
				s >> num_dof;
				num_dof_set=true;
			} else if ( tag == "SCORE_MULTIPLY_SIZE" ) {
				s >> num_score_multiply;
				num_score_multiply_set=true;
			} else if ( tag == "BB_CLONES" ) {
				if ( ! num_bb_indep_set ) { utility_exit_with_message("Must have BB_CLONES_SIZE before BB_CLONES in silent file."); }
				comma_strings_to_vector_map( s,   num_bb_indep,  symminfo.bb_clones_, "" );
				bb_clones_set=true;
			} else if ( tag == "CHI_CLONES" ) {
				if ( ! num_chi_indep_set ) { utility_exit_with_message("Must have CHI_CLONES_SIZE before CHI_CLONES in silent file."); }
				comma_strings_to_vector_map( s,  num_chi_indep,  symminfo.chi_clones_, "" );
				chi_clones_set=true;
			} else if ( tag == "JUMP_CLONES" ) {
				if ( ! num_jump_indep_set ) { utility_exit_with_message("Must have JUMP_CLONES_SIZE before JUMP_CLONES in silent file."); }
				comma_strings_to_vector_map( s, num_jump_indep,  symminfo.jump_clones_, "" );
				jump_clones_set=true;
			} else if ( tag == "JUMP_CLONE_WEIGHTS" ) {
				if ( ! num_jump_dep_set || ! num_jump_indep_set ) { utility_exit_with_message("Must have JUMP_CLONES_SIZE and JUMP_FOLLOWS_SIZE before JUMP_CLONE_WEIGHTS in silent file."); }
				comma_strings_to_map( s, num_jump_dep-num_jump_indep,  symminfo.jump_clone_wts_, "" );
				jump_clone_wts_set=true;
			} else if ( tag == "BB_FOLLOWS" ) {
				if ( ! num_bb_dep_set ) { utility_exit_with_message("Must have BB_FOLLOWS_SIZE before BB_FOLLOWS in silent file."); }
				comma_strings_to_map( s,   num_bb_dep,  symminfo.bb_follows_, "" );
				bb_follows_set=true;
			} else if ( tag == "CHI_FOLLOWS" ) {
				if ( ! num_chi_dep_set ) { utility_exit_with_message("Must have CHI_FOLLOWS_SIZE before CHI_FOLLOWS in silent file."); }
				comma_strings_to_map( s,  num_chi_dep,  symminfo.chi_follows_, "" );
				chi_follows_set=true;
			} else if ( tag == "JUMP_FOLLOWS" ) {
				if ( ! num_jump_dep_set ) { utility_exit_with_message("Must have JUMP_FOLLOWS_SIZE before JUMP_FOLLOWS in silent file."); }
				comma_strings_to_map( s, num_jump_dep,  symminfo.jump_follows_, "" );
				jump_follows_set=true;
			} else if ( tag == "DOFS" ) {
				if ( ! num_dof_set ) { utility_exit_with_message("Must have DOFS_SIZE before DOFS in silent file."); }
				comma_strings_to_map( s, num_dof,  symminfo.dofs_, "" );
				dofs_set=true;
			} else if ( tag == "SCORE_MULTIPLY" ) {
				if ( ! num_score_multiply_set ) { utility_exit_with_message("Must have SCORE_MULTIPLY_SIZE before SCORE_MULTIPLY in silent file."); }
				comma_strings_to_vector( s, num_score_multiply,  symminfo.score_multiply_, "" );
				score_multiply_set = true;
			} else if ( tag == "N_COMPONENTS" ) {
				s >> symminfo.num_components_;
				num_components_set = true;
			} else if ( tag == "JUMP2DOFNAME" ) {
				if ( ! nnamed_jumps_set ) { utility_exit_with_message("Must have N_NAMED_JUMPS before JUMP2DOFNAME in silent file."); }
				comma_strings_to_map( s, nnamed_jumps, symminfo.jnum2dofname_, "" );
				jnum2dofname_set = true;
			} else if ( tag == "DOFNAME2JUMP" ) {
				if ( ! nnamed_jumps_set ) { utility_exit_with_message("Must have N_NAMED_JUMPS before DOFNAME2JUMP in silent file."); }
				comma_strings_to_map( s, nnamed_jumps, symminfo.dofname2jnum_, "" );
				dofname2jnum_set = true;
			} else if ( tag == "COMPONENTS" ) {
				if ( ! num_components_set ) { utility_exit_with_message("Must have N_COMPONENTS before COMPONENTS in silent file."); }
				comma_strings_to_vector( s, symminfo.num_components_,  symminfo.components_, "" );
				//components_set = true;
			} else if ( tag == "COMPONENT_BOUNDS" ) {
				if ( ! num_components_set ) { utility_exit_with_message("Must have N_COMPONENTS before COMPONENT_BOUNDS in silent file."); }
				comma_strings_to_pair_map( s, symminfo.num_components_,  symminfo.component_bounds_, "" );
				//component_bounds_set = true;
			} else if ( tag == "NAME2COMPONENT" ) {
				if ( ! num_components_set ) { utility_exit_with_message("Must have N_COMPONENTS before NAME2COMPONENT in silent file."); }
				comma_strings_to_map( s, symminfo.num_components_,  symminfo.name2component_, "" );
				//name2component_set = true;
			} else if ( tag == "JNAME2COMPS" ) {
				if ( ! nnamed_jumps_set ) { utility_exit_with_message("Must have N_NAMED_JUMPS before JNAME2COMPS in silent file."); }
				comma_strings_to_vector_map( s, nnamed_jumps,  symminfo.jname2components_, "" );
				//jname2components_set = true;
			} else if ( tag == "JNAME2SUBS" ) {
				if ( ! nnamed_jumps_set ) { utility_exit_with_message("Must have N_NAMED_JUMPS before JNAME2SUBS in silent file."); }
				comma_strings_to_vector_map( s, nnamed_jumps,  symminfo.jname2subunits_, "" );
				//jname2subunits_set = true;
			} else {
				TR.Debug << "Unknown tag {" << tag << "}" << std::endl;
			}

			if ( s.fail() ) {
				TR << "Error while reading tag {" << tag << "}" << std::endl;
				fail = true;
				break;
			}
		}

		// try to recover from missing fields (that are not taken care of in the default constructor)
		// fail if required field is missing
		if ( !npseudo_set || !interfaces_set || !type_set || !num_bb_indep_set || !num_chi_indep_set || !num_jump_indep_set ||
				!num_bb_dep_set || !num_chi_dep_set || !num_jump_dep_set || !num_dof_set || !num_score_multiply_set ||
				!bb_clones_set || !chi_clones_set || !jump_clones_set  || !bb_follows_set ||
				!chi_follows_set || !jump_follows_set || !dofs_set || !score_multiply_set ) {
			TR.Error << "Failed to read the following fields: ";
			if ( !npseudo_set )  TR.Error << "npseudo ";
			if ( !interfaces_set )  TR.Error << "interfaces ";
			if ( !type_set )  TR.Error << "type ";
			if ( !num_bb_indep_set )  TR.Error << "num_bb_indep ";
			if ( !num_chi_indep_set )  TR.Error << "num_chi_indep ";
			if ( !num_jump_indep_set )  TR.Error << "num_jump_indep ";
			if ( !num_bb_dep_set )  TR.Error << "num_bb_dep ";
			if ( !num_chi_dep_set )  TR.Error << "num_chi_dep ";
			if ( !num_jump_dep_set )  TR.Error << "num_jump_dep ";
			if ( !num_dof_set )  TR.Error << "num_dof ";
			if ( !num_score_multiply_set )  TR.Error << "num_score_multiply ";
			if ( !bb_clones_set )  TR.Error << "bb_clones ";
			if ( !chi_clones_set )  TR.Error << "chi_clones ";
			if ( !jump_clones_set )  TR.Error << "jump_clones ";
			if ( !bb_follows_set )  TR.Error << "bb_follows ";
			if ( !chi_follows_set )  TR.Error << "chi_follows ";
			if ( !jump_follows_set )  TR.Error << "jump_follows ";
			if ( !dofs_set )  TR.Error << "dofs ";
			if ( !score_multiply_set )  TR.Error << "score_multiply ";
			TR.Error << std::endl;
			fail = true;
		}

		if ( !fail ) {
			bool warned = false;

			// now handle any missing fields
			if ( !jump_clone_wts_set ) {
				if ( warned ) {
					TR.Warning << "Symmetric silent file is missing fields; attempting automatic recovery." << std::endl;
					//AMW: cppcheck correctly flags this line as being unreachable
					// every time around the loop warned starts as false and can never reach places where it is made true
					//warned = true;
				}
				// set master jumps to 1; clones to 0
				for ( std::map<Size,SymmetryInfo::Clones>::const_iterator map_it=symminfo.jump_clones_.begin(),
						map_end=symminfo.jump_clones_.end();
						map_it != map_end; ++map_it ) {
					for ( Size i=1; i<=map_it->second.size(); ++i ) {
						symminfo.jump_clone_wts_[ map_it->second[i] ] = 0;
					}
				}
			}

			if ( !nres_monomer_set ) {
				symminfo.nres_monomer_ = num_bb_indep;
			}

			if ( !scoring_subunit_set ) {
				std::map<Size,SymmetryInfo::Clones>::const_iterator first_bb_clone=symminfo.bb_clones_.begin();
				symminfo.scoring_subunit_ = 1 + ((first_bb_clone->first-1) / symminfo.nres_monomer_);
			}

			//fpd this could almost certainly be smarter
			if ( !njump_monomer_set ) {
				symminfo.njump_monomer_ = 0;
			}

			// assume single comp
			if ( !num_components_set ) {
				symminfo.num_components_ = 1;
			}

			// warn on anything else
			if ( !nnamed_jumps_set || !jnum2dofname_set || !dofname2jnum_set ) {
				if ( warned ) {
					TR.Warning << "Symmetric silent file is missing fields; attempting automatic recovery." << std::endl;
					//warned = true;
				}

				for ( std::map<Size,SymmetryInfo::Clones>::const_iterator map_it=symminfo.jump_clones_.begin(),
						map_end=symminfo.jump_clones_.end();
						map_it != map_end; ++map_it ) {
					// add master
					std::ostringstream ossmaster;
					ossmaster << "JUMP" << map_it->first;
					symminfo.dofname2jnum_[ossmaster.str()] = map_it->first;
					symminfo.jnum2dofname_[map_it->first] = ossmaster.str();

					for ( Size i=1; i<=map_it->second.size(); ++i ) {
						std::ostringstream ossclone;
						ossclone << "JUMP" << map_it->second[i];
						symminfo.dofname2jnum_[ossclone.str()] = map_it->second[i];
						symminfo.jnum2dofname_[map_it->second[i]] = ossclone.str();
					}
				}
			}

		} else {
			// fail
			TR.Error << "Symmetry_info failed while reading silent file!" << std::endl;
			s.setstate( std::ios_base::failbit );
			return s;
		}

		// finally, update derived info
		for ( Size i=1; i <=symminfo.num_total_residues_without_pseudo(); ++i ) {
			if ( symminfo.bb_is_independent(i) ) {
				symminfo.last_indep_residue_ = i;
			}
		}
	}

	symminfo.update_score_multiply_factor();

	return s;
}

/////////////////////////////////////////////////////////////////////////////
std::ostream& operator<< ( std::ostream & s, const SymmetryInfo & symminfo )
{
	s << "SYMMETRY_INFO " <<
		"N " << symminfo.subunits() << ' ' <<
		"N_RES_MONOMER " << symminfo.nres_monomer_ << ' ' <<
		"SCORING_SUBUNIT " << symminfo.scoring_subunit_ << ' ' <<
		"N_JUMP_MONOMER " << symminfo.njump_monomer_ << ' ' <<
		"N_NAMED_JUMPS " << symminfo.jnum2dofname_.size() << ' ' <<
		"N_VIRT " << symminfo.npseudo_ << ' ' <<
		"N_INTERFACE " << symminfo.num_interfaces() << ' ' <<
		"N_COMPONENTS " << symminfo.num_components_ << ' ' <<
		"TYPE " << symminfo.type_ << ' ' <<
		"BB_CLONES_SIZE " << symminfo.bb_clones_.size() << ' ' <<
		"CHI_CLONES_SIZE " << symminfo.chi_clones_.size() << ' ' <<
		"JUMP_CLONES_SIZE " << symminfo.jump_clones_.size() << ' ' <<
		"BB_FOLLOWS_SIZE " << symminfo.bb_follows_.size() << ' ' <<
		"CHI_FOLLOWS_SIZE " << symminfo.chi_follows_.size() << ' ' <<
		"JUMP_FOLLOWS_SIZE " << symminfo.jump_follows_.size() << ' ' <<
		"DOFS_SIZE " << symminfo.dofs_.size() << ' ' <<
		"SCORE_MULTIPLY_SIZE " << symminfo.score_multiply_.size() ;

	// clones
	vector_map_to_comma_strings( s, symminfo.bb_clones_, "BB_CLONES" );
	vector_map_to_comma_strings( s, symminfo.chi_clones_, "CHI_CLONES" );
	vector_map_to_comma_strings( s, symminfo.jump_clones_, "JUMP_CLONES" );
	map_to_comma_strings( s, symminfo.jump_clone_wts_, "JUMP_CLONE_WEIGHTS" );

	// follows
	map_to_comma_strings( s, symminfo.bb_follows_, "BB_FOLLOWS" );
	map_to_comma_strings( s, symminfo.chi_follows_, "CHI_FOLLOWS" );
	map_to_comma_strings( s, symminfo.jump_follows_, "JUMP_FOLLOWS" );

	//dof
	map_to_comma_strings( s, symminfo.dofs_, "DOFS" );

	//score_multiply_
	vector_to_comma_strings( s, symminfo.score_multiply_, "SCORE_MULTIPLY" );

	// jumpnames, multicomp stuff
	map_to_comma_strings( s, symminfo.jnum2dofname_, "JUMP2DOFNAME" );
	map_to_comma_strings( s, symminfo.dofname2jnum_, "DOFNAME2JUMP" );

	if ( symminfo.components_.size() > 0 ) {
		vector_to_comma_strings( s, symminfo.components_, "COMPONENTS" );
		pair_map_to_comma_strings( s, symminfo.component_bounds_, "COMPONENT_BOUNDS" );
		map_to_comma_strings( s, symminfo.name2component_, "NAME2COMPONENT" );
		vector_map_to_comma_strings( s, symminfo.jname2components_, "JNAME2COMPS" );
		vector_map_to_comma_strings( s, symminfo.jname2subunits_, "JNAME2SUBS" );
	}

	return s;
}

bool
SymmetryInfo::write_silent_struct(
	std::string const & filename
)
{
	bool success = false;

	utility::io::ozstream output;
	if ( !utility::file::file_exists( filename ) ) {
		output.open( filename );
	} else {
		output.open_append( filename );
	}

	output << *this << '\n';

	output.close();

	success = true;
	return success;
}


bool
SymmetryInfo::read_silent_struct(
	std::string const & filename
)
{
	bool success = false;

	utility::io::izstream input ( filename.c_str() );
	std::istringstream line_stream;
	std::string line("");
	if ( !input ) {
		std::cerr << "ERROR:: Unable to open symmetry info file: "
			<< filename << std::endl;
		return success;
	}

	while ( !input.eof() ) {
		getline(input,line);
		line_stream.clear();
		line_stream.str(line);
		line_stream >> *this;
	}

	input.close();

	success = true;
	return success;
}

bool
SymmetryInfo::is_virtual( Size const seqpos ) const {
	return ( seqpos > num_total_residues_without_pseudo() );
}

Size
SymmetryInfo::bb_follows( Size const seqpos ) const
{
	std::map< Size, Size >::const_iterator it( bb_follows_.find( seqpos ) );
	return ( it == bb_follows_.end() ? 0 : it->second );
}

Size
SymmetryInfo::chi_follows( Size const seqpos ) const
{
	std::map< Size, Size >::const_iterator it( chi_follows_.find( seqpos ) );
	return ( it == chi_follows_.end() ? 0 : it->second );
}

Size
SymmetryInfo::jump_follows( Size const seqpos ) const
{
	std::map< Size, Size >::const_iterator it( jump_follows_.find( seqpos ) );
	return ( it == jump_follows_.end() ? 0 : it->second );
}

std::vector < std::pair < Size, Size > >
SymmetryInfo::map_symmetric_res_pairs( Size res1, Size res2 )
const {
	std::vector < std::pair < Size, Size > > map;
	int delta ( res2 - res1 );
	int mapped_res;
	for ( unsigned long clone : bb_clones( res1 ) ) {
		if ( clone + delta > num_total_residues() ) {
			mapped_res = (clone + delta)%num_total_residues();
		} else {
			mapped_res = clone + delta;
		}
		if ( mapped_res < 0 ) {
			mapped_res += num_independent_residues();
		}
		map.emplace_back( clone, mapped_res );
	}
	return map;
}

bool
SymmetryInfo::bb_is_independent( Size const seqpos ) const
{
	return bb_follows(seqpos) == 0;
}

bool
SymmetryInfo::chi_is_independent( Size const seqpos ) const
{
	return chi_follows(seqpos) == 0;
}

bool
SymmetryInfo::fa_is_independent( Size const seqpos ) const
{
	return ( bb_is_independent(seqpos) && chi_is_independent(seqpos) );
}

bool
SymmetryInfo::jump_is_independent( Size const seqpos ) const
{
	return jump_follows(seqpos) == 0;
}

Size
SymmetryInfo::subunits() const
{
	return num_bb_clones() + 1;
}

utility::vector1< bool >
SymmetryInfo::independent_residues() const
{
	utility::vector1 < bool > residues ( num_total_residues_with_pseudo() );
	for ( Size i=1; i <= num_total_residues_with_pseudo(); ++i ) {
		residues[ i ] = bb_is_independent( i );
	}
	return residues;
}

Size
SymmetryInfo::last_independent_residue() const {
	return last_indep_residue_;
}

Size
SymmetryInfo::num_bb_clones() const
{
	// all these lists have the same size
	if ( bb_clones_.empty() ) {
		return 0;
	}

	return bb_clones_.begin()->second.size();
}

Size
SymmetryInfo::num_chi_clones() const
{
	// all these lists have the same size
	return chi_clones_.begin()->second.size();
}


Size
SymmetryInfo::num_jump_clones() const
{
	// all these lists have the same size
	return jump_clones_.begin()->second.size();
}

Size
SymmetryInfo::num_independent_residues() const
{
	return bb_clones_.size();
}

Size
SymmetryInfo::num_total_residues() const
{
	return num_independent_residues()*( num_bb_clones() + 1 );
}

Size
SymmetryInfo::num_total_residues_with_pseudo() const
{
	return num_independent_residues()*( num_bb_clones() + 1 ) + npseudo_;
}

Size
SymmetryInfo::num_total_residues_without_pseudo() const
{
	return num_independent_residues()*( num_bb_clones() + 1 );
}

Size
SymmetryInfo::num_interfaces() const
{
	return interfaces_;
}

Size
SymmetryInfo::score_multiply_factor() const
{
	return score_multiply_factor_;
}

Size
SymmetryInfo::num_virtuals() const
{
	return npseudo_;
}

void
SymmetryInfo::num_virtuals( Size const setting )
{
	npseudo_ = setting;
}

SymmetryInfo::Clones const &
SymmetryInfo::bb_clones( Size const seqpos ) const
{
	std::map< Size, Clones >::const_iterator it( bb_clones_.find( seqpos ) );
	if ( it == bb_clones_.end() ) {
		return empty_list;
	}
	return it->second;
}

SymmetryInfo::Clones const &
SymmetryInfo::chi_clones( Size const seqpos ) const
{
	std::map< Size, Clones >::const_iterator it( chi_clones_.find( seqpos ) );
	if ( it == chi_clones_.end() ) {
		return empty_list;
	}
	return it->second;
}

SymmetryInfo::Clones const &
SymmetryInfo::jump_clones( Size const seqpos ) const
{
	std::map< Size, Clones >::const_iterator it( jump_clones_.find( seqpos ) );
	if ( it == jump_clones_.end() ) {
		return empty_list;
	}
	return it->second;
}

//fpd remap bb_clones/chi_clones when the ASU size changes
//fpd this recreates the arrays from scratch so it may be somewhat inefficient
//sboyken@gmail.com 15/05/03 updating with Frank's code; properly update score_multiply_factor, etc. when new residue appended
void
SymmetryInfo::resize_asu( Size nres_new ) {
	if ( nres_new == nres_monomer_ ) return; // nothing to do

	Size N = subunits();
	Size nres_monomer_old = nres_monomer_;

	nres_monomer_ = nres_new;
	bb_clones_.clear();
	bb_follows_.clear();
	chi_clones_.clear();
	chi_follows_.clear();

	for ( Size i=1; i<= nres_monomer_; ++i ) {
		Clones clones;
		int base ( i + ( scoring_subunit_ - 1 ) * nres_monomer_ );
		for ( Size k=0; k<N; ++k ) {
			if ( k+1 != scoring_subunit_ ) {
				clones.push_back( i + k * nres_monomer_ );
				add_bb_clone( base, i + k * nres_monomer_ );
				add_chi_clone( base, i + k * nres_monomer_ );
			}
		}
		bb_clones_.insert( std::make_pair( base, clones ) );
		chi_clones_.insert( std::make_pair( base, clones ) );
	}
	// update last_independent_residue
	for ( Size i=1; i <=num_total_residues_without_pseudo(); ++i ) {
		if ( bb_is_independent(i) ) {
			last_indep_residue_ = i;
		}
	}

	// update score multiply
	utility::vector1< Size > score_multiply_new(num_total_residues_with_pseudo(), 1);
	for ( Size i=1; i<= nres_monomer_; ++i ) {
		for ( Size k=0; k<N; ++k ) {
			Size res_old = 1 + k * nres_monomer_old;  // just use 1st res of subunit
			Size res_new = i + k * nres_monomer_;
			score_multiply_new[res_new] = score_multiply_[res_old];
		}
	}
	score_multiply_ = score_multiply_new;
}

//fpd remap jump_clones when the number of monomer jumps changes
void
SymmetryInfo::update_nmonomer_jumps( Size njump_monomer ) {
	if ( njump_monomer == njump_monomer_ ) return; // nothing to do

	//std::cerr << "SymmetryInfo::update_nmonomer_jumps(" << njump_monomer << ")  [old=" << njump_monomer_ << "]\n";
	Size N = subunits();

	// remember previous
	std::map< Size, Clones > old_jump_clones = jump_clones_;
	std::map< Size, Real > old_jump_clone_weights = jump_clone_wts_;
	Size old_njump_monomer = njump_monomer_;

	njump_monomer_ = njump_monomer;
	jump_clones_.clear();
	jump_follows_.clear();
	jump_clone_wts_.clear();

	// make new monomer jumps from scratch
	// the N*njump_monomer internal jumps
	std::map<Size,std::string> jnum2dofname_new;
	std::map<std::string,Size> dofname2jnum_new;

	for ( Size i=1; i<= njump_monomer; ++i ) {
		Size src_intra = (i-1)*N + scoring_subunit_; //i + (scoring_subunit_-1)*njump_monomer_;

		std::string name_new_src = "intra_"+utility::to_string(i)+"_master";
		jnum2dofname_new[src_intra] = name_new_src;
		dofname2jnum_new[name_new_src] = src_intra;

		for ( Size k=0; k<N; ++k ) {
			if ( k != ( scoring_subunit_ - 1 ) ) {
				Size tgt_intra = (i-1)*N + k + 1; //i + k*njump_monomer_;
				add_jump_clone( src_intra, tgt_intra, 0.0 );

				std::string name_new_tgt = "intra_"+utility::to_string(i)+"_"+utility::to_string(k+1);
				jnum2dofname_new[tgt_intra] = name_new_tgt;
				dofname2jnum_new[name_new_tgt] = tgt_intra;
			}
		}
	}

	//fd update jump names
	for ( auto & j_it : jnum2dofname_ ) {
		int i_new = (int)j_it.first + (int)N*((int)njump_monomer-(int)old_njump_monomer);
		std::string old_name = j_it.second;
		if ( i_new > 0 ) {
			jnum2dofname_new[i_new] = old_name;
			dofname2jnum_new[old_name] = i_new;
		}
	}


	// 1                 --> N*njump_monomer  : the internal jumps
	// N*njump_monomer+1 --> N*njump_monomer+N: the pseudo-rsd--monomer jumps
	// last N-1 jumps                         : jumps between pseudo-rsds
	for ( auto it=old_jump_clones.begin(), it_end = old_jump_clones.end(); it != it_end; ++it ) {
		Size source = it->first;
		Clones target = it->second;
		std::string old_source_name = jnum2dofname_[source];

		if ( source > N*old_njump_monomer ) {
			// a symm jump
			Size new_source = source + N*( njump_monomer - old_njump_monomer );

			for ( Size i=1; i<=target.size(); ++i ) {
				Size new_target_i = target[i] + N*(njump_monomer-old_njump_monomer);

				add_jump_clone( new_source, new_target_i, old_jump_clone_weights[target[i]] );
			}
		}
	}

	// dofs
	std::map< Size, SymDof > dofs_new;
	for ( auto & dof : dofs_ ) {
		dofs_new.insert( std::make_pair( dof.first + N*( njump_monomer - old_njump_monomer ), dof.second ) );
	}
	set_dofs( dofs_new );

	jnum2dofname_ = jnum2dofname_new;
	dofname2jnum_ = dofname2jnum_new;
}

/// this actually only tests to see if the independent residues are contiguous...
///
void
SymmetryInfo::update_contiguous_monomers()
{
	utility::vector1< Size > independent_positions;
	for ( std::map< Size, Size >::const_iterator it= bb_follows_.begin(); it != bb_follows_.end(); ++it ) {
		if ( it->second == 0 ) independent_positions.push_back( it->first );
	}

	Size first_indep( utility::min( independent_positions ) ), last_indep( utility::max( independent_positions )),
		n_indep( independent_positions.size() );

	contiguous_monomers_ = ( last_indep - first_indep + 1 == n_indep );

	// TR.Trace << "SymmetryInfo::update_contiguous_monomers: " << first_indep << ' ' << last_indep << ' '<< n_indep << ' ' <<
	//  contiguous_monomers_ << std::endl;

}

void
SymmetryInfo::add_bb_clone( Size const base_pos, Size const clone_pos )
{
	if ( bb_follows_[ base_pos ] != 0 ) {
		std::cerr << "Error: add_bb_clone: base_pos is not independent: " <<
			base_pos << ' ' << bb_follows_[ base_pos ] << std::endl;
		utility_exit();
	}
	if ( bb_follows_[ clone_pos ] != 0 &&
			bb_follows_[ clone_pos ] != base_pos ) {
		std::cerr << "Error: add_bb_clone: clone_pos already a follower: " <<
			clone_pos << ' ' << bb_follows_[ clone_pos ] << ' ' << base_pos <<
			std::endl;
		utility_exit();
	}

	bb_follows_[ clone_pos ] = base_pos;
	bb_clones_[ base_pos ].push_back( clone_pos );
	update_contiguous_monomers();
}

void
SymmetryInfo::add_chi_clone( Size const base_pos, Size const clone_pos )
{
	if ( chi_follows_[ base_pos ] != 0 ) {
		std::cerr << "Error: add_chi_clone: base_pos is not independent: " <<
			base_pos << ' ' << chi_follows_[ base_pos ] << std::endl;
		utility_exit();
	}
	if ( chi_follows_[ clone_pos ] != 0 &&
			chi_follows_[ clone_pos ] != base_pos ) {
		std::cerr << "Error: add_chi_clone: clone_pos already a follower: " <<
			clone_pos << ' ' << chi_follows_[ clone_pos ] << ' ' << base_pos <<
			std::endl;
		utility_exit();
	}

	chi_follows_[ clone_pos ] = base_pos;
	chi_clones_[ base_pos ].push_back( clone_pos );
}

void
SymmetryInfo::add_jump_clone( Size const base_pos, Size const clone_pos, Real const jump_wt )
{
	if ( jump_follows_[ base_pos ] != 0 ) {
		std::cerr << "Error: add_jump_clone: base_pos is not independent: " <<
			base_pos << ' ' << bb_follows_[ base_pos ] << std::endl;
		utility_exit();
	}
	if ( jump_follows_[ clone_pos ] != 0 &&
			jump_follows_[ clone_pos ] != base_pos ) {
		std::cerr << "Error: add_jump_clone: clone_pos already a follower: " <<
			clone_pos << ' ' << jump_follows_[ clone_pos ] << ' ' << base_pos <<
			std::endl;
		utility_exit();
	}

	jump_follows_[ clone_pos ] = base_pos;
	jump_clones_[ base_pos ].push_back( clone_pos );
	jump_clone_wts_[ clone_pos ] = jump_wt;
}

std::map< Size, SymDof > const &
SymmetryInfo::get_dofs() const
{
	return dofs_;
}

void
SymmetryInfo::set_dofs( std::map< Size, SymDof > const & dofs )
{
	dofs_ = dofs;
}

Real
SymmetryInfo::score_multiply( Size const res1, Size const res2 ) const
{
	bool bb1ind = bb_is_independent(res1);
	bool bb2ind = bb_is_independent(res2);
	if ( res1 > num_total_residues_without_pseudo() ) {
		return ( bb2ind ) ? reweight_symm_interactions_*score_multiply_[res2] : 0;
	} else if  ( res2 > num_total_residues_without_pseudo() ) {
		return ( bb1ind ) ? reweight_symm_interactions_*score_multiply_[res1] : 0;
	} else if ( bb1ind ) {
		return (bb2ind?1:reweight_symm_interactions_)*score_multiply_[res2];
	} else if ( bb2ind ) {
		return reweight_symm_interactions_*score_multiply_[res1];
	} else {
		return 0;
	}
}

Real
SymmetryInfo::deriv_multiply( Size const res1, Size const res2 ) const
{
	bool bb1ind = bb_is_independent(res1);
	bool bb2ind = bb_is_independent(res2);
	if ( res1 > num_total_residues_without_pseudo() ) {
		return ( bb2ind ) ? 1 : 0;
	} else if ( res2 > num_total_residues_without_pseudo() ) {
		return ( bb1ind ) ? 1 : 0;
	} else if ( bb1ind && bb2ind ) {
		return 1.0;
	} else if ( bb1ind || bb2ind ) {
		return reweight_symm_interactions_;
	} else {
		return 0.0;
	}
}


Size
SymmetryInfo::interface_number( Size const res1, Size const res2 ) const
{
	return subunit_index(bb_is_independent(res1) ? res2 : res1);
}

void
SymmetryInfo::set_score_multiply_from_subunit_factors( utility::vector1< Real > const & score_multiply_vector_subunit, Size const nres_subunit, Size const n_subunits )
{
	score_multiply_.clear();
	for ( Size i = 1; i<= n_subunits; ++i ) {
		for ( Size j = 1; j<= nres_subunit; ++j ) {
			score_multiply_.push_back( score_multiply_vector_subunit[i] );
		}
	}
	for ( Size i = n_subunits + 1 ; i <= score_multiply_vector_subunit.size() ; ++i ) {
		score_multiply_.push_back( score_multiply_vector_subunit[i] );
	}
}

void
SymmetryInfo::set_score_multiply( Size const res, Size const factor )
{
	debug_assert ( res <= score_multiply_.size() );
	score_multiply_[ res ] = factor;
}

void
SymmetryInfo::set_flat_score_multiply( Size const nres, Size const factor )
{
	score_multiply_.clear();
	score_multiply_.resize( nres, factor );
}

Size
SymmetryInfo::get_nres_subunit() const
{
	return nres_monomer_;
}

void
SymmetryInfo::set_nres_subunit( Size const setting )
{
	nres_monomer_ = setting;
}

Size
SymmetryInfo::get_njumps_subunit() const
{
	return njump_monomer_;
}


bool
SymmetryInfo::get_use_symmetry() const
{
	return use_symmetry_;
}

SymSlideInfo
SymmetryInfo::get_slide_info() const
{
	return slide_info_;
}

void
SymmetryInfo::set_use_symmetry( bool setting )
{
	use_symmetry_ = setting;
}

/// @details  Returns set of DOF_IDs that follow a given one. Inefficient: it creates a list each time
SymmetryInfo::DOF_IDs
SymmetryInfo::dependent_dofs( DOF_ID const & id, Conformation const & conf ) const
{
	if ( !dof_is_independent( id, conf ) ) {
		utility_exit_with_message( "SymmetryInfo::dependent_dofs: dof is not independent!" );
	}

	Size const seqpos( id.rsd() );
	Size const atomno( id.atomno() );
	id::TorsionType const type( guess_torsion_type_of_dof_id( id, conf ) );
	Clones const & clones( type == id::JUMP ? jump_clones( conf.fold_tree().get_jump_that_builds_residue( id.rsd() ) ):
		( type == id::BB ? bb_clones( seqpos ) : chi_clones( seqpos ) ) );

	DOF_IDs dofs;
	for ( unsigned long clone : clones ) {
		if ( type == id::JUMP ) {
			dofs.push_back( DOF_ID( id::AtomID( atomno, conf.fold_tree().downstream_jump_residue( clone )  ), id.type() ) );
		} else {
			dofs.push_back( DOF_ID( id::AtomID( atomno, clone ), id.type() ) );
		}
	}
	return dofs;
}

bool
SymmetryInfo::dof_is_independent( DOF_ID const & id, Conformation const & conf ) const
{
	id::TorsionType const type( guess_torsion_type_of_dof_id( id, conf ) );

	switch ( type ) {
	case id::BB :
		return bb_is_independent( id.rsd() );
	case id::CHI :
		return chi_is_independent( id.rsd() );
	case id::JUMP :
		return jump_is_independent( conf.fold_tree().get_jump_that_builds_residue( id.rsd() ) );
	default :
		utility_exit_with_message("dof_is_independent: unrecognized TorsionType!");
		return false;
	}
}

// get a weight for derivative calculations
// weights are 1 for indep DOFs, 0 for dependent NON-JUMP DOFs
//    and may be any real for dependent jump dofs
core::Real
SymmetryInfo::get_dof_derivative_weight( DOF_ID const & id, Conformation const & conf ) const {
	id::TorsionType const type( guess_torsion_type_of_dof_id( id, conf ) );

	if ( type == id::BB ) {
		return bb_is_independent( id.rsd() ) ? 1. : 0.;
	} else if ( type == id::CHI ) {
		return chi_is_independent( id.rsd() ) ? 1. : 0.;
	} else if ( type == id::JUMP ) {
		int jumpnum = conf.fold_tree().get_jump_that_builds_residue( id.rsd() );
		std::map< Size, Real >::const_iterator it( jump_clone_wts_.find( jumpnum ) );
		return ( it == jump_clone_wts_.end() ? 1. : it->second );
	}

	utility_exit_with_message("get_dof_derivative_weight: unrecognized TorsionType!");
	return false;
}


bool
SymmetryInfo::torsion_is_independent( TorsionID const & id ) const
{
	return ( ( id.type() == id::BB   &&   bb_is_independent( id.rsd() ) ) ||
		( id.type() == id::CHI  &&  chi_is_independent( id.rsd() ) ) ||
		( id.type() == id::JUMP && jump_is_independent( id.rsd() ) ) );
}

bool
SymmetryInfo::atom_is_independent( AtomID const & id ) const
{
	return fa_is_independent( id.rsd() );
}

/// @details  Returns set of TorsionIDs that follow a given one. Inefficient: it creates a list each time
SymmetryInfo::TorsionIDs
SymmetryInfo::dependent_torsions( TorsionID const & id ) const
{
	if ( !torsion_is_independent( id ) ) {
		utility_exit_with_message( "SymmetryInfo::dependent_torsions: torsion is not independent!" );
	}

	Size const seqpos( id.rsd() );
	Clones const & seqpos_clones( id.type() == id::BB ? bb_clones( seqpos ) : chi_clones( seqpos ) );

	//std::cerr << " dependent_torsions( TorsionID const & )  bb=" << (id.type() == id::BB) << "   chi=" << (id.type() == id::CHI) << "   jump=" << (id.type() == id::JUMP) << std::endl;

	TorsionIDs tors;
	for ( unsigned long seqpos_clone : seqpos_clones ) {
		tors.push_back( TorsionID( seqpos_clone, id.type(), id.torsion() ) );
	}
	return tors;
}

/// @details  Returns set of AtomIDs that follow a given one. Inefficient: it creates a list each time
SymmetryInfo::AtomIDs
SymmetryInfo::dependent_atoms( AtomID const & id ) const
{
	if ( !atom_is_independent( id ) ) {
		utility_exit_with_message( "SymmetryInfo::dependent_atoms: atom is not independent!" );
	}
	Size const seqpos( id.rsd() );
	Clones const & seqpos_clones( bb_clones( seqpos ) );

	AtomIDs atoms;
	for ( unsigned long seqpos_clone : seqpos_clones ) {
		atoms.push_back( AtomID( id.atomno(), seqpos_clone ) );
	}
	return atoms;
}

bool
SymmetryInfo::is_asymmetric_seqpos( Size const res ) const
{
	Size nres_monomer = num_independent_residues(), num_monomers = subunits() ;
	return ( !get_use_symmetry() || res > nres_monomer*num_monomers || res <= nres_monomer );
}

Size
SymmetryInfo::get_asymmetric_seqpos( Size const res ) const
{
	if ( res > 0 && get_use_symmetry() ) {
		Size nres_monomer = num_independent_residues(), num_monomers = subunits() ;
		return ( res > nres_monomer*num_monomers ? res - nres_monomer*(num_monomers-1) : (res-1)%nres_monomer + 1 );
	} else {
		return res;
	}
}

Size
SymmetryInfo::subunit_index( Size const seqpos ) const {
	if ( contiguous_monomers_ ) {
		// the old logic
		return ( (seqpos-1) / num_independent_residues() + 1 );
	} else {
		// the new way: this depends critically on the clones all having been added in the same order...
		// slower, obviously
		Size const basepos( bb_follows_.find( seqpos )->second );
		Size index(0);
		if ( basepos ) {
			Clones const & clones( bb_clones_.find( basepos )->second );
			index = std::find( clones.begin(), clones.end(), seqpos ) - clones.begin() + 2;
		} else {
			index = 1; // the independent monomer is numbered 1
		}
		//TR.Trace << "SymmetryInfo::subunit_index: " << seqpos << ' ' << index << std::endl;
		return index;
	}
}

Size
SymmetryInfo::equivalent_residue_on_subunit( Size subunit_index, Size residue_id ) const
{
	if ( contiguous_monomers_ ) {
		return ( subunit_index - 1 ) * num_independent_residues() + ( (residue_id-1) % nres_monomer_ + 1 ) ;
	} else {
		Size basepos = bb_follows_.find( residue_id )->second;
		basepos = basepos == 0 ? residue_id : basepos;
		if ( subunit_index == 1 ) {
			return basepos;
		} else {
			return bb_clones_.find( basepos )->second[ subunit_index - 1 ];
		}
	}
}

std::string
SymmetryInfo::get_jump_name(Size jnum) const {
	if ( 0 == jnum2dofname_.count(jnum) ) utility_exit_with_message("bad jump num");
	return jnum2dofname_.find(jnum)->second;
}

Size
SymmetryInfo::get_jump_num(std::string jname) const {
	if ( 0 == dofname2jnum_.count(jname) ) utility_exit_with_message("bad jump name");
	return dofname2jnum_.find(jname)->second;
}

void
SymmetryInfo::set_jump_name(Size jnum, std::string jname) {
	jnum2dofname_[jnum] = jname;
	dofname2jnum_[jname] = jnum;
}

Size
SymmetryInfo::num_slidablejumps() const {
	Size retval = 0;
	for ( auto const & dof : dofs_ ) {
		if ( dof.second.allow_dof(1) || dof.second.allow_dof(2) || dof.second.allow_dof(3) ) retval++;
	}
	return retval;
}

Size const &
SymmetryInfo::get_num_components() const {
	return num_components_;
}

utility::vector1<char> const &
SymmetryInfo::get_components() const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	return components_;
}

char
SymmetryInfo::get_component(Size i) const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	return components_[i];
}

std::map<char,std::pair<Size,Size> > const &
SymmetryInfo::get_component_bounds() const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	return component_bounds_;
}

std::map<std::string,char> const &
SymmetryInfo::get_subunit_name_to_component() const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	return name2component_;
}

std::map<std::string,utility::vector1<char> > const &
SymmetryInfo::get_jump_name_to_components() const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	return jname2components_;
}

std::map<std::string,utility::vector1<Size> > const &
SymmetryInfo::get_jump_name_to_subunits() const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	return jname2subunits_;
}

std::pair<Size,Size> const &
SymmetryInfo::get_component_bounds(char c) const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	if ( component_bounds_.find(c) == component_bounds_.end() ) {
		utility_exit_with_message(std::string("no symmetry component ")+c);
	}
	return component_bounds_.find(c)->second;
}
Size
SymmetryInfo::get_component_lower_bound(char c) const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	if ( component_bounds_.find(c) == component_bounds_.end() ) {
		utility_exit_with_message(std::string("no symmetry component ")+c);
	}
	return component_bounds_.find(c)->second.first;
}
Size
SymmetryInfo::get_component_upper_bound(char c) const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	if ( component_bounds_.find(c) == component_bounds_.end() ) {
		utility_exit_with_message(std::string("no symmetry component ")+c);
	}
	return component_bounds_.find(c)->second.second;
}
char
SymmetryInfo::get_component_of_residue(Size ir) const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	if ( ir > num_total_residues_without_pseudo() || ir < 1 ) {
		utility_exit_with_message(std::string("no symmetry component for residue "));
	}
	Size irindep = (ir-1)%num_independent_residues()+1;
	for ( auto const & component_bound : component_bounds_ ) {
		char component = component_bound.first;
		Size lower = component_bound.second.first;
		Size upper = component_bound.second.second;
		// std::cerr << component << " " << lower << " " << upper << " " << irindep << std::endl;
		if ( lower <= irindep && irindep <= upper ) return component;
	}
	utility_exit_with_message(std::string("no symmetry component for residue "));
	return ' ';
}
char
SymmetryInfo::get_subunit_name_to_component(std::string const & vname) const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	if ( name2component_.find(vname) == name2component_.end() ) {
		utility_exit_with_message(std::string("no symmetry component for ")+vname);
	}
	return name2component_.find(vname)->second;
}
utility::vector1<char> const &
SymmetryInfo::get_jump_name_to_components(std::string const & jname) const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	if ( jname2components_.find(jname) == jname2components_.end() ) {
		utility_exit_with_message(std::string("no symmetry component for ")+jname);
	}
	return jname2components_.find(jname)->second;
}
utility::vector1<Size> const &
SymmetryInfo::get_jump_name_to_subunits(std::string const & jname) const {
	if ( components_.size()==0 ) utility_exit_with_message("function not for use in single component symmetry");
	if ( jname2subunits_.find(jname) == jname2subunits_.end() ) {
		utility_exit_with_message(std::string("no symmetry component for ")+jname);
	}
	return jname2subunits_.find(jname)->second;
}

void
SymmetryInfo::set_multicomponent_info(
	Size const & num_components,
	utility::vector1<char> const & components,
	std::map<char,std::pair<Size,Size> > const & component_bounds,
	std::map<std::string,char> const & name2component,
	std::map<std::string,utility::vector1<char> > const & jname2component,
	std::map<std::string,utility::vector1<Size> > const & jname2subunits
){
	num_components_ = num_components;
	components_ = components;
	component_bounds_ = component_bounds;
	name2component_ = name2component;
	jname2components_ = jname2component;
	jname2subunits_ = jname2subunits;
}

} // symmetry
} // conformation
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::symmetry::SymmetryInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( bb_clones_ ) ); // std::map<Size, Clones>
	arc( CEREAL_NVP( chi_clones_ ) ); // std::map<Size, Clones>
	arc( CEREAL_NVP( jump_clones_ ) ); // std::map<Size, Clones>
	arc( CEREAL_NVP( jump_clone_wts_ ) ); // std::map<Size, Real>
	arc( CEREAL_NVP( bb_follows_ ) ); // std::map<Size, Size>
	arc( CEREAL_NVP( chi_follows_ ) ); // std::map<Size, Size>
	arc( CEREAL_NVP( jump_follows_ ) ); // std::map<Size, Size>
	arc( CEREAL_NVP( empty_list ) ); // Clones
	arc( CEREAL_NVP( nres_monomer_ ) ); // Size
	arc( CEREAL_NVP( scoring_subunit_ ) ); // Size
	arc( CEREAL_NVP( npseudo_ ) ); // Size
	arc( CEREAL_NVP( njump_monomer_ ) ); // Size
	arc( CEREAL_NVP( last_indep_residue_ ) ); // Size
	arc( CEREAL_NVP( type_ ) ); // std::string
	arc( CEREAL_NVP( interfaces_ ) ); // Size
	arc( CEREAL_NVP( score_multiply_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( reweight_symm_interactions_ ) ); // Real
	arc( CEREAL_NVP( score_multiply_factor_ ) ); // core::Size
	arc( CEREAL_NVP( dofs_ ) ); // std::map<Size, SymDof>
	arc( CEREAL_NVP( use_symmetry_ ) ); // _Bool
	arc( CEREAL_NVP( slide_info_ ) ); // class core::conformation::symmetry::SymSlideInfo
	arc( CEREAL_NVP( jnum2dofname_ ) ); // std::map<Size, std::string>
	arc( CEREAL_NVP( dofname2jnum_ ) ); // std::map<std::string, Size>
	arc( CEREAL_NVP( contiguous_monomers_ ) ); // _Bool
	arc( CEREAL_NVP( torsion_changes_move_other_monomers_ ) ); // _Bool
	arc( CEREAL_NVP( num_components_ ) ); // Size
	arc( CEREAL_NVP( components_ ) ); // utility::vector1<char>
	arc( CEREAL_NVP( component_bounds_ ) ); // std::map<char, std::pair<Size, Size> >
	arc( CEREAL_NVP( name2component_ ) ); // std::map<std::string, char>
	arc( CEREAL_NVP( jname2components_ ) ); // std::map<std::string, utility::vector1<char> >
	arc( CEREAL_NVP( jname2subunits_ ) ); // std::map<std::string, utility::vector1<Size> >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::symmetry::SymmetryInfo::load( Archive & arc ) {
	arc( bb_clones_ ); // std::map<Size, Clones>
	arc( chi_clones_ ); // std::map<Size, Clones>
	arc( jump_clones_ ); // std::map<Size, Clones>
	arc( jump_clone_wts_ ); // std::map<Size, Real>
	arc( bb_follows_ ); // std::map<Size, Size>
	arc( chi_follows_ ); // std::map<Size, Size>
	arc( jump_follows_ ); // std::map<Size, Size>
	arc( empty_list ); // Clones
	arc( nres_monomer_ ); // Size
	arc( scoring_subunit_ ); // Size
	arc( npseudo_ ); // Size
	arc( njump_monomer_ ); // Size
	arc( last_indep_residue_ ); // Size
	arc( type_ ); // std::string
	arc( interfaces_ ); // Size
	arc( score_multiply_ ); // utility::vector1<Size>
	arc( reweight_symm_interactions_ ); // Real
	arc( score_multiply_factor_ ); // core::Size
	arc( dofs_ ); // std::map<Size, SymDof>
	arc( use_symmetry_ ); // _Bool
	arc( slide_info_ ); // class core::conformation::symmetry::SymSlideInfo
	arc( jnum2dofname_ ); // std::map<Size, std::string>
	arc( dofname2jnum_ ); // std::map<std::string, Size>
	arc( contiguous_monomers_ ); // _Bool
	arc( torsion_changes_move_other_monomers_ ); // _Bool
	arc( num_components_ ); // Size
	arc( components_ ); // utility::vector1<char>
	arc( component_bounds_ ); // std::map<char, std::pair<Size, Size> >
	arc( name2component_ ); // std::map<std::string, char>
	arc( jname2components_ ); // std::map<std::string, utility::vector1<char> >
	arc( jname2subunits_ ); // std::map<std::string, utility::vector1<Size> >
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::symmetry::SymmetryInfo );
CEREAL_REGISTER_TYPE( core::conformation::symmetry::SymmetryInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_symmetry_SymmetryInfo )
#endif // SERIALIZATION

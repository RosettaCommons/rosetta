// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/TorsionSampler.cc
///
/// @brief Sample torsions based on experimental distributions
/// @author Guangfeng Zhou and Frank DiMaio

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/database/open.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/scoring/GenericBondedPotential.fwd.hh>
#include <core/scoring/GenericBondedPotential.hh>

#include <utility/io/izstream.hh>

#include <numeric/random/random.hh>


#include <protocols/ligand_docking/GALigandDock/TorsionSampler.hh>

#include <core/chemical/ChemicalManager.hh> // AUTO IWYU For ChemicalManager
#include <core/chemical/AtomTypeSet.hh> // AUTO IWYU For AtomTypeSet
#include <utility/stream_util.hh> // AUTO IWYU For operator<<


namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

static basic::Tracer TR( "protocols.ligand_docking.GALigandDock.TorsionSampler" );
const TorsionDistrParams TorsionSampler::null_tors_distr = TorsionDistrParams("null_torsion", 0);

enum ReadMode {
	rmNONE,
	rmATOM,
	rmBOND,
	rmANGLE,
	rmTORSION,
	rmIMPROPER
};

TorsionDistrParams::TorsionDistrParams():
	n_mode_(0), mult_(1E12), normalized_(false), torsion_type_("")
{ }

TorsionDistrParams::TorsionDistrParams(std::string const& ttype_in,
	core::Size mult_in): n_mode_(0), mult_(mult_in), normalized_(false), torsion_type_(ttype_in)
{ }

void
TorsionDistrParams::add_mode( core::Real A, core::Real mu, core::Real sigma){
	A_vec_.push_back(A);
	mu_vec_.push_back(mu);
	sigma_vec_.push_back(sigma);
	if ( A_cumulative_vec_.empty() ) {
		A_cumulative_vec_.push_back(0.0);
	}
	core::Real A_sum = A_cumulative_vec_.back() + A;
	A_cumulative_vec_.push_back(A_sum);
	n_mode_ = A_vec_.size();
}

void
TorsionDistrParams::get_mode_i(core::Size i, core::Real& A_out,
	core::Real& mu_out, core::Real& sigma_out
) const
{
	assert(i<=n_mode_);
	A_out = A_vec_[i];
	mu_out = mu_vec_[i];
	sigma_out = sigma_vec_[i];
}

void
TorsionDistrParams::normalize()
{
	if ( normalized_ ) return;
	if ( n_mode_ == 0 ) {
		normalized_=true;
		return;
	}
	core::Real mag = 0.0, inv_mag;
	for ( core::Real const& A:A_vec_ ) {
		mag = mag + A;
	}
	if ( mag != 1.0 && mag != 0.0 ) {
		inv_mag = 1.0 / mag;
		for ( core::Size i=1; i<=n_mode_; ++i ) {
			A_vec_[i] = A_vec_[i] * inv_mag;
		}
		TR.Debug << "A_vec_" << A_vec_ << std::endl;
		for ( core::Size i=1; i<=A_cumulative_vec_.size(); ++i ) {
			A_cumulative_vec_[i] = A_cumulative_vec_[i] * inv_mag;
		}
		TR.Debug << "A_cumulative_vec_" << A_cumulative_vec_ << std::endl;
	}
	normalized_ = true;
}

void
TorsionDistrParams::sample_mode (core::Real& A_out,
	core::Real& mu_out,
	core::Real& sigma_out)
const {
	assert(normalized_==true);
	core::Real begin = 0, end = 0 ;
	core::Real p = numeric::random::rg().uniform();
	for ( core::Size i=1; i<=A_cumulative_vec_.size()-1; ++i ) {
		begin = A_cumulative_vec_[i];
		end = A_cumulative_vec_[i+1];
		if ( p >= begin && p < end ) {
			A_out = A_vec_[i];
			mu_out = mu_vec_[i];
			sigma_out = sigma_vec_[i];
			if ( TR.Debug.visible() ) {
				TR.Debug << "p, begin, end, A, mu, sigma: " << p <<", "
					<< begin << ", " << end <<", "
					<< A_out <<", " << mu_out << ", "
					<< sigma_out << std::endl;
			}
			return;
		}
	}
	if ( TR.Debug.visible() ) {
		TR.Debug << " mode not found! p: " << p << std::endl;
	}
}


TorsionSampler::TorsionSampler()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string db_file = option[ score::torsion_distr_params_file ]();
	read_database( db_file );
}


void
TorsionSampler::read_database( std::string const & filename)
{
	using namespace core::scoring;
	utility::io::izstream instream;
	basic::database::open( instream, filename );

	core::chemical::AtomTypeSetCOP ats = core::chemical::ChemicalManager::get_instance()->atom_type_set( "fa_standard" );
	core::Size ntypes = ats->n_atomtypes();
	defined_atom_types_.resize( ntypes, false );
	//runtime_assert( ntypes < 65536 ); // we pack 4 atom indices into 64 bits
	runtime_assert( ntypes < 4096 ); // should be safe?

	// map atom names to ATS indices - handle wildcards
	ReadMode read_mode = rmNONE;
	core::Size linenum( 0 );
	core::Size natomtypes = 0;
	std::string fileline, tag;
	// special treatment of "fallbacks", e.g. 'X' in the params file
	utility::vector1< core::Size > wildcard_vec(1, 0);

	while ( instream ) {
		getline(instream, fileline);
		std::istringstream linestream(fileline);
		linenum ++;

		if ( fileline.length() < 2 ) continue;

		linestream >> tag;
		if ( tag == "ATOM" ) {
			read_mode = rmATOM;
			continue;
		} else if ( tag == "TORSION" ) {
			read_mode = rmTORSION;
			continue;
		} else if ( tag.substr(0,1) == "#" ) {
			continue; // comment or empty line
		}

		if ( read_mode == rmATOM ) {       // in "ATOM" block
			core::Size atm_idx = ats->atom_type_index(tag);
			defined_atom_types_[ atm_idx ] = true;
			TR.Debug << "AtomicIndex: " << atm_idx << " " << tag << std::endl;
			name_index_map[ tag ].push_back( atm_idx );
			std::string alt_tag;
			while ( linestream >> alt_tag ) {
				if ( alt_tag != tag ) name_index_map[ alt_tag ].push_back( atm_idx );
			}
			natomtypes++;

		} else if ( read_mode == rmTORSION ) { // in "TORSION" block
			std::string atm1=tag, atm2, bt, atm3, atm4;
			core::Real A, mu, sigma;
			std::string dummy, temp;

			linestream >> atm2 >> bt >> atm3 >> atm4 >> dummy;
			std::string torsion_type_in = atm1+","+atm2+","+bt+","+atm3+","+atm4;
			if ( TR.Debug.visible() ) {
				TR.Debug << "torsion_type_in: " << torsion_type_in << std::endl;
			}

			utility::vector1< core::Size > indicesBT = bondorders_map(bt);
			utility::vector1< core::Size > const & indices1 = name_index_map[atm1];
			utility::vector1< core::Size > const & indices2 = name_index_map[atm2];
			utility::vector1< core::Size > const & indices3 = name_index_map[atm3];
			utility::vector1< core::Size > const & indices4 = name_index_map[atm4];

			// multiplicity is used for priority ...
			//   the applicable torsion with the lowest multiplicity is selcted
			// we want bondtype to dominate ... if there is a more-specific bt, always prefer that
			//   even if atomtyping is less specific
			// however we need to ensure "fallbacks" aren't favored too much
			utility::vector1< core::Size > const & indicesX = name_index_map["X"];
			core::Size const multi_max( indicesX.size()/2 );
			core::Size multBT = indicesBT.size()*multi_max*multi_max*multi_max*multi_max;

			core::Size multiplicity = multBT + indices1.size() * indices2.size() * indices3.size() * indices4.size();

			TorsionDistrParams tordistparams(torsion_type_in, multiplicity );

			while ( linestream >> A >> mu >> sigma ) {
				tordistparams.add_mode(A, mu, sigma);
				if ( TR.Debug.visible() ) {
					TR.Debug << "distr params: " << A <<", " << mu
						<<", " << sigma << ", " << std::endl;
				}
			}

			tordistparams.normalize();

			while ( linestream >> temp ) {}

			if ( TR.Debug.visible() ) {
				TR.Debug << "---------" << std::endl;
				TR.Debug << "Size of the bondtype: " << bt << ": " << indicesBT.size() << std::endl;
				TR.Debug << "Size of the X: " << indicesX.size() << " half: " << multi_max << std::endl;
				TR.Debug << "Size of the atom type: " << atm1 <<", " << atm2 << ", " << atm3 << ", " << atm4
					<<": " << indices1.size() <<", " << indices2.size() <<", "
					<< indices3.size() <<", " << indices4.size() << std::endl;
				TR.Debug << "multBT: " << multBT << ", multiplicity: " << multiplicity << std::endl;

			}

			utility::vector1< core::Size > const &indices1_loop = (indices1.size() == natomtypes) ? wildcard_vec : indices1;
			utility::vector1< core::Size > const &indices2_loop = (indices2.size() == natomtypes) ? wildcard_vec : indices2;
			utility::vector1< core::Size > const &indices3_loop = (indices3.size() == natomtypes) ? wildcard_vec : indices3;
			utility::vector1< core::Size > const &indices4_loop = (indices4.size() == natomtypes) ? wildcard_vec : indices4;

			// flipping atom order reverses phase -- generate a new constraint if necessary
			tors_distr_.push_back( tordistparams );

			core::Size tgt_pot_idx1 = tors_distr_.size();

			for ( auto i0 : indicesBT ) {
				for ( auto i1 : indices1_loop ) {
					for ( auto i2 : indices2_loop ) {
						for ( auto i3 : indices3_loop ) {
							for ( auto i4 : indices4_loop ) {

								uint64_t hashval = get_parameter_hash( i0, i1, i2, i3, i4 );
								auto it = tors_lookup_.find( hashval );

								// use the MOST SPECIFIC potential
								if ( it == tors_lookup_.end()  ) {
									tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx1) );
									hashval = get_parameter_hash( i0, i4, i3, i2, i1 );
									tors_lookup_.insert( std::make_pair(hashval, tgt_pot_idx1) );
								} else if ( multiplicity < tors_distr_[ it->second ].multiplicity() ) {
									it->second = tgt_pot_idx1;
									hashval = get_parameter_hash( i0, i4, i3, i2, i1 );
									it = tors_lookup_.find( hashval );
									it->second = tgt_pot_idx1;
								}
							}
						}
					}
				}
			}
		}
	} //end while

	TR << "Added total " << tors_distr_.size() << " torsion distribution parameters corresponding to "
		<< tors_lookup_.size() << " unique torsions." << std::endl;
}

core::Real
TorsionSampler::sample(
	core::chemical::BondName bn,
	core::chemical::BondRingness br,
	core::Size type1, core::Size type2,
	core::Size type3, core::Size type4)
const {
	TorsionDistrParams const& tor_distr_params =
		lookup_tors_distr_params(bn, br, type1, type2, type3, type4);
	core::Real A(0.0), mu(0.0), sigma(0.0), retval;
	if ( TR.Debug.visible() ) {
		TR.Debug << "Found tor_distr_params for "<< tor_distr_params.torsion_type()
			<< " has "<< tor_distr_params.n_mode() << " modes." << std::endl;
	}

	if ( tor_distr_params.n_mode() == 0 ) {
		retval = 360.0 * numeric::random::rg().uniform();
		if ( TR.Debug.visible() ) {
			TR.Debug << "fully random sample: " << retval << std::endl;
		}
	} else {
		tor_distr_params.sample_mode(A, mu, sigma);
		retval = mu + numeric::random::rg().gaussian()*sigma;
		if ( retval>180 ) {
			retval =  retval - 360;
		} else if ( retval< -180 ) {
			retval = retval + 360;
		}
	}
	if ( TR.Debug.visible() ) {
		TR.Debug << "TorsionSampler, sample torion: "
			<< tor_distr_params.torsion_type() << " angle: "
			<< retval <<std::endl;
	}
	return retval;
}

TorsionDistrParams const &
TorsionSampler::lookup_tors_distr_params(
	core::chemical::BondName bn,
	core::chemical::BondRingness br,
	core::Size type1, core::Size type2,
	core::Size type3, core::Size type4 ) const
{
	using namespace core::scoring;

	if ( !defined_atom_types_[type1] || !defined_atom_types_[type2] || !defined_atom_types_[type3]|| !defined_atom_types_[type4] ) return null_tors_distr;

	core::Size btidx = bin_from_bond(bn, br);

	//fd look up tgt with best multiplicity
	auto it = tors_lookup_.find( get_parameter_hash(btidx, type1, type2, type3, type4) );

	auto it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, type2, type3, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, 0, type3, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, type2, 0, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, type2, type3, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, 0, type3, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, type2, 0, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, type2, type3, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, 0, 0, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, 0, type3, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, type2, 0, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, type1, 0, 0, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, type2, 0, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, 0, type3, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, 0, 0, type4) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;
	it2 = tors_lookup_.find( get_parameter_hash(btidx, 0, 0, 0, 0) );
	if ( it2 != tors_lookup_.end() &&
			(it == tors_lookup_.end() || tors_distr_[it2->second].multiplicity() < tors_distr_[it->second].multiplicity()) ) it = it2;

	// Final sanity check... (this should never get triggered)
	if ( it == tors_lookup_.end() ) it = tors_lookup_.find( get_parameter_hash(0, 0, 0, 0, 0) );

	return tors_distr_[it->second];
}


}
}
}

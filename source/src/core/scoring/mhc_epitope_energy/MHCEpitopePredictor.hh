// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopePredictor.hh
/// @brief Base class for an MHC epitope predictor, which takes a peptide (string) and returns a score predicting its risk of MHC binding (lower is lower risk, with 0 being none)
/// @details The peptides are of fixed length, specified for the predictor (e.g., for class II MHC: 9 for Propred, 15 for NetMHCII and IEDB)
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#ifndef INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictor_hh
#define INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictor_hh

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictor.fwd.hh>
#include <core/chemical/AA.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

class MHCEpitopePredictor: public utility::pointer::ReferenceCount {

public:
	MHCEpitopePredictor() {}
	virtual ~MHCEpitopePredictor() {}

	/// @brief Is it the same predictor, not just as initialized (e.g., from file), but considering any subsequent modifications?
	/// @details Takes any other epitope predictor so needs to make sure same type.
	virtual bool operator==(MHCEpitopePredictor const & /* other */) = 0;

	/// @brief Get a summary of the data stored in this object
	virtual std::string report() const { return "some kind of epitope predictor"; }

	/// @brief The MHC epitope score, predicting "immunogenicity" of the peptide (lower is lower risk, with 0 being none)
	/// @details The peptide must be of the length specified for this predictor (get_peptide_length()))
	virtual core::Real score(std::string const & /* peptide */) { return 0; }

	/// @brief The length required by this predictor for peptides
	core::Size get_peptide_length() const { return peptide_length_; }

private:
	/// @brief The length required by this predictor for peptides
	/// @details Must be set at some point before using the predictor, but since it might only be determined upon load from disk/database, not required by constructor
	core::Size peptide_length_=0;

protected:
	/// @brief Setter for peptide_length, for derived classes to set upon initialization
	void set_peptide_length(core::Size peptide_length) { peptide_length_ = peptide_length; }

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class MHCEpitopePredictor

} // mhc_epitope_energy
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopePredictor )
#endif // SERIALIZATION


#endif

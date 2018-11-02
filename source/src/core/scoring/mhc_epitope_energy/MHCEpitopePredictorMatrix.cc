// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopePredictorMatrix.cc
/// @brief MHC epitope predictor using a position weight matrix, targeted to Propred though in theory generalizable to others
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <basic/database/open.hh>
#include <iostream>
#include <string>
#include <utility/string_util.hh>
#include <core/chemical/AA.hh>

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorMatrix.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopePredictorMatrix");

AlleleMatrix::AlleleMatrix()
{}

AlleleMatrix::AlleleMatrix(std::string name, utility::vector1< Real > threshes, PWM profile)
: name_(name), threshes_(threshes), profile_(profile)
{}

AlleleMatrix::~AlleleMatrix()
{}

bool AlleleMatrix::operator==(AlleleMatrix const &other)
{
	// TODO: does name matter?
	return profile_ == other.profile_ && threshes_ == other.threshes_;
}

bool AlleleMatrix::is_hit(std::string const &pep, Real thresh)
{
	Real total(0);
	// Loop over the profile of this allele (i.e. each position).
	for ( core::Size p=1; p<=profile_.size(); p++ ) {
		// Look up the appropriate weight in the profile for the peptide's AA at this position, and add it to total.
		total += profile_[p][pep[p-1]];
	}

	// If the total for this peptide is greater than the threshold, it is a hit. Return whether ot not it is.
	return total >= threshes_[(core::Size)thresh];
}

MHCEpitopePredictorMatrix::MHCEpitopePredictorMatrix()
{}

MHCEpitopePredictorMatrix::MHCEpitopePredictorMatrix( std::string const &fn )
{
	load_matrix(fn); // Load the matrix file provided as an argument to the constructor.
}

MHCEpitopePredictorMatrix::~MHCEpitopePredictorMatrix()
{}

bool MHCEpitopePredictorMatrix::operator==(MHCEpitopePredictor const &other)
{
	MHCEpitopePredictorMatrix const *o = dynamic_cast<MHCEpitopePredictorMatrix const *>(&other);
	if ( !o ) return false;

	if ( o->propred_ != propred_ ) return false;
	if ( o->thresh_ != thresh_ ) return false;

	// individually check allele matches, since .mhc file might have subselected
	if ( o->alleles_.size() != alleles_.size() ) return false;
	for ( core::Size a=1; a<=alleles_.size(); a++ ) {
		bool matched = false;
		for ( core::Size oa=1; oa<=alleles_.size() && !matched; oa++ ) {
			matched = alleles_[a] == o->alleles_[oa];
		}
		if ( !matched ) return false;
	}

	return true;
}

std::string MHCEpitopePredictorMatrix::report() const
{
	std::stringstream output("");

	// TODO: more (allele names, ...)?
	output << "Matrix predictor using " << filename_ << "; " << alleles_.size() << " alleles; threshold " << thresh_;

	return output.str();
}

core::Real MHCEpitopePredictorMatrix::score(std::string const &pep)
{
	if ( pep.size() != get_peptide_length() ) {
		TR.Error << "Scoring peptide of size " << pep.size() << " with a matrix expecting peptides of size " << get_peptide_length() << std::endl;
		utility_exit_with_message("MHCEpitopePredictorMatrix is trying to score a peptide of the incorrect size!");
	}

	core::Real total(0);
	// Loop over the allele matrices, considering them separately.
	for ( core::Size a=1; a<=alleles_.size(); a++ ) {
		// Count the total number of allele's saying that the peptide is predicted to be a hit (binder)
		if ( alleles_[a].is_hit(pep, thresh_) ) total++;
	}

	return total;
}

/// @brief Gets the next line that is not blank and not commented out by a #.
///
utility::io::izstream &get_useful_line(utility::io::izstream &stream, std::string &line) {
	while ( getline(stream, line) ) {
		if ( line.size() > 0 && line[0] != '#' ) break;
	}
	return stream;
}

void MHCEpitopePredictorMatrix::load_matrix(std::string const &filename)
{
	filename_ = filename; // TODO: use resolved_fn instead?

	// TODO: This implementation is very fragile
	// TODO: It is also propred-specific, though if there were another ready example of a PWM-based predictor, could easily be generalized

	// The filename allows specifying different sets of matrices, so find the right one.
	utility::file::FileName resolved_fn = basic::database::full_name("scoring/score_functions/mhc_epitope/"+filename+".txt");
	utility::io::izstream input(resolved_fn);
	if ( !input ) utility_exit_with_message("ERROR: Unable to open file " + filename + ", resolved to "+resolved_fn.name());
	TR << "Reading matrix predictor from " << resolved_fn.name() << std::endl;

	// To allow for expansion of format, the first non-blank, non-comment line says the format.
	// Currently that must be propred. Otherwise exit with an error.
	std::string line;
	get_useful_line(input, line);
	if ( line != "propred" ) {
		utility_exit_with_message("ERROR: Unknown epitope predictor " + line);
	} else {
		propred_ = true;
	}

	// The second line should indicate the length of the peptides (9 for Propred)
	get_useful_line(input, line);
	set_peptide_length(utility::string2int(line));

	// The third line gives the amino acid order in the matrix, using one letter codes and no spaces.
	std::string aas;
	get_useful_line(input, aas);
	if ( aas.size() != 20 ) utility_exit_with_message("ERROR: Wrong # AA types for epitope predictor " + line);

	// The fourth line gives the number of alleles being specified in the matrix.
	get_useful_line(input, line);
	Size nallele = utility::string2int(line);
	TR.Debug << nallele << " alleles" << std::endl;
	alleles_.resize(nallele);

	// For the rest of the file, loop over the number of alleles being specified, reading the matrix for each
	for ( core::Size a=1; a<=nallele; a++ ) {
		// Read the name of the allele.
		std::string name;
		get_useful_line(input, name);
		TR.Debug << "reading allele # " << a << " '" << name << "'" << std::endl;
		// The first line in the matrix is the list of thresholds.
		// The thresholds indicate the cutoff in the table that represents the best X% of binders.
		// In a Propred matrix, the "X%" corresponds to 1-10% in increments of 1 percentage point.
		// We store this in the vector threshes.  In a propred matrix, threshes[X] will be the score cutoff
		// for the top X% of binders.
		get_useful_line(input, line);
		std::istringstream thresh_stream( line );
		utility::vector1< Real > threshes;
		core::Real thresh;
		// TODO: We should probably make this a vector of struct (or something) that contain the
		// % threshold and the corresponding score.  The idea of using the index as the percent cutoff
		// seems very fragile to me.
		while ( thresh_stream >> thresh ) threshes.push_back(thresh);

		// Read the PWM
		// There is one line for each position in the peptide (9 for Propred).
		AlleleMatrix::PWM profile;
		profile.resize(get_peptide_length());
		// Loop over the line for each position (9 for Propred).
		// Store the weights in row, and then store each row in profile to get the full matrix.
		for ( core::Size p=1; p<=get_peptide_length(); p++ ) {
			AlleleMatrix::Weights row;
			get_useful_line(input, line);
			std::istringstream row_stream( line );
			for ( core::Size aa=0; aa<20; aa++ ) {
				row_stream >> row[aas[aa]];
			}
			profile[p] = row;
		}
		// Create the matrix and store it in the alleles_ vector.
		alleles_[a] = AlleleMatrix(name, threshes, profile);
	}
}

/// @brief Sets the threshold for what is considered to be an epitope -- top thresh% of peptides in this implementation
/// @details Includes error checking in the propred matrix case to make sure we get a reasonable threshold
void MHCEpitopePredictorMatrix::set_thresh(core::Real thresh) {
	if ( propred_ ) {
		if ( thresh < 1 || thresh > 10 || (core::Size)thresh != thresh ) {
			TR.Error << "A propred matrix must have an integer threshold between 1-10." << std::endl;
			TR.Error << "The score corresponding to these thresholds are explicitly given the matrix." << std::endl;
			utility_exit_with_message( "The threshold of " + std::to_string(thresh) + " that you have selected is invalid." );
		}
	}
	thresh_ = thresh;
}

}//ns mhc_epitope_energy
}//ns scoring
}//ns core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::mhc_epitope_energy::AlleleMatrix::save( Archive & arc ) const {
	arc( CEREAL_NVP( name_ ) ); // std::string
	arc( CEREAL_NVP( threshes_ ) ); // core::vector1
	arc( CEREAL_NVP( profile_ ) ); // PWM
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::mhc_epitope_energy::AlleleMatrix::load( Archive & arc ) {
	arc( name_ ); // std::string
	arc( threshes_ ); // core::vector1
	arc( profile_ ); // PWM
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::mhc_epitope_energy::AlleleMatrix );
CEREAL_REGISTER_TYPE( core::scoring::mhc_epitope_energy::AlleleMatrix )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_AlleleMatrix )

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopePredictorMatrix::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::mhc_epitope_energy::MHCEpitopePredictor >( this ) );
	arc( CEREAL_NVP( filename_ ) ); // std::string
	arc( CEREAL_NVP( alleles_ ) ); // utility::vector1<AlleleMatrix>
	arc( CEREAL_NVP( thresh_ ) ); // core::Real
	arc( CEREAL_NVP( propred_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopePredictorMatrix::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::mhc_epitope_energy::MHCEpitopePredictor >( this ) );
	arc( filename_ ); // std::string
	arc( alleles_ ); // utility::vector1<AlleleMatrix>
	arc( thresh_ ); // core::Real
	arc( propred_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::mhc_epitope_energy::MHCEpitopePredictorMatrix );
CEREAL_REGISTER_TYPE( core::scoring::mhc_epitope_energy::MHCEpitopePredictorMatrix )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopePredictorMatrix )
#endif // SERIALIZATION

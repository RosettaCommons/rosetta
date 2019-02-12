// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/TCRloopRefine.hh
/// @brief Class for parsing and numbering tcr a/b chain sequences
/// @author Ragul Gowthaman (ragul@umd.edu)

#ifndef INCLUDED_protocols_tcr_TCRseqInfo_hh
#define INCLUDED_protocols_tcr_TCRseqInfo_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/antibody/grafting/antibody_sequence.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace tcr {

/// @brief class for parsing/numbering tcr sequences
/// @details TCRseqInfo takes in the alpha/beta chain sequence
/// TCR segments are parsed by regular expressions / Aho numbering / ANARCI / user provided numbering
/// TCR segments: Germline, Framework, CDR1, CDR2, CDR2HV4, CDR3
/// Extended cdr2hv4 segment is used for cdr2. More details can be found in publication/documentation
class TCRseqInfo : public utility::pointer::ReferenceCount {
public:

	TCRseqInfo( std::string const &aseq, std::string const &bseq );
	//TCRseqInfo();

	TCRseqInfoOP clone() const;

	//Default destructor
	~TCRseqInfo() override;

	/// @brief tcrposi holds the start/end position of cdr loops
	/// @details cap has the start/end of the variable domian sequence
	struct tcrposi {
		antibody::grafting::CDR_Bounds cdr1;
		antibody::grafting::CDR_Bounds cdr2hv4;
		antibody::grafting::CDR_Bounds cdr3;
		antibody::grafting::CDR_Bounds cap;
	};

	/// @brief parsed sequence strings of tcr segments
	/// @details truncdomain:variable domain sequence, fr:framework sequence, gm:germline sequence
	/// @details cdr1, cdr2, cdr3 cdr loop sequences
	struct tcrsegs {
		std::string truncdomain;
		std::string fr;
		std::string gm;
		std::string cdr1;
		std::string cdr2;
		std::string cdr3;
	};

	/// @brief returns parsed tcr segments for a/b chains
	tcrsegs atcr() const { return atcr_; };
	tcrsegs btcr() const { return btcr_; };

	/// @brief returns cdr positions for input a/b chains
	tcrposi aseq_posi() const { return aseq_posi_;};
	tcrposi bseq_posi() const { return bseq_posi_;};

	/// @brief returns standard Aho numbering used for cdr positions for a/b chains
	tcrposi aaho_posi() const { return aaho_posi_;};
	tcrposi baho_posi() const { return baho_posi_;};

private:

	/// @brief initialize from options
	void init();

	/// @brief set default values for Aho numbers
	void set_default();

	/// @brief function to parse input alpaha and beta TCR sequences into tcr segments
	/// @details uses regular expression by default for parsing
	void parse_tcrsegs(std::string alpha_chain_sequence, std::string beta_chain_sequence);

private:

	/// @brief parsed tcr segments for tcr a/b chains
	tcrsegs atcr_;
	tcrsegs btcr_;

	/// @brief aho numbers from cdr definitions file for tcr a/b chains
	tcrposi aaho_posi_;
	tcrposi baho_posi_;

	/// @brief standard aho numbers from cdr definitions file for tcr a/b chains
	tcrposi aseq_posi_;
	tcrposi bseq_posi_;

	/// @brief cdr positions for the input a/b chain sequence
	tcrposi ausr_posi_;
	tcrposi busr_posi_;

	/// @brief assign cdr positions by user provided numbers
	bool cdr_posi_by_num_;

	/// @brief path for the ANARCI program used for parsing tcr a/b chains
	std::string anarci_path_;

}; // class TCRseqInfo

} // tcr
} // protocols

#endif // __ANTIBODY_GRAFTING__

#endif

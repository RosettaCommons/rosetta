// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/scs_blast.cc
/// @brief Structural Component Selector (SCS) implementation with NCBI-BLAST+
/// @author Sergey Lyskov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/scs_blast.hh>
#include <protocols/antibody/grafting/exception.hh>
#include <protocols/antibody/grafting/chothia_numberer.hh>
#include <protocols/antibody/grafting/scs_functor.hh>

#include <basic/Tracer.hh>
#include <basic/execute.hh>

#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <utility/string_util.hh>
//#include <utility/stream_util.hh>

#include <algorithm>
#include <fstream>
#include <locale>

namespace protocols {
namespace antibody {
namespace grafting {

using std::string;

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting");




/// @brief Create result set from 'row' element of each SCS_ResultsVector vector
///        if strict is true then throw if no resuls found otherwise use empty OP
///
/// @throw std::out_of_range if for some of the row is not present and strict==true
SCS_ResultSet SCS_Results::get_result_set(uint row, bool strict)
{
	if( h1.size()  <= row  and  h2.size()  <= row  and  h3.size() <= row  and
		l1.size()  <= row  and  l2.size()  <= row  and  l3.size() <= row  and
		frh.size() <= row  and  frl.size() <= row  and  orientation.size() <= row  and  strict) throw std::out_of_range("SCS_Results::get_result_set could not find ALL results for row:" + std::to_string(row));

	SCS_ResultSet r;

	r.h1 = h1.size() > row ? h1[row] : SCS_ResultOP();
	r.h2 = h2.size() > row ? h2[row] : SCS_ResultOP();
	r.h3 = h3.size() > row ? h3[row] : SCS_ResultOP();

	r.l1 = l1.size() > row ? l1[row] : SCS_ResultOP();
	r.l2 = l1.size() > row ? l1[row] : SCS_ResultOP();
	r.l3 = l1.size() > row ? l1[row] : SCS_ResultOP();

	r.frh = frh.size() > row ? frh[row] : SCS_ResultOP();
	r.frl = frl.size() > row ? frl[row] : SCS_ResultOP();

	r.orientation = orientation.size() > row ? orientation[row] : SCS_ResultOP();

	return r;
}



void SCS_Base::add_filter(SCS_FunctorCOP filter)
{
	filters_.push_back(filter);
}

void SCS_Base::set_sorter(SCS_FunctorCOP sorter)
{
	sorter_ = sorter;
}

struct ResultsItem {
	string name;
	SCS_ResultsVector & result;
};




SCS_ResultsOP SCS_Base::select(AntibodySequence const &A)
{
	SCS_ResultsOP r = raw_select(A);

	for(auto &filter : filters_) filter->apply(A, r);

	if(sorter_) sorter_->apply(A, r);

	return r;
}




/// @details Parse plain text file with following format, Replace some of the charactes in legend
/// <legend-prefix> field-1 field-2 field-3 ...
/// f11 f12 f13 ...
/// f21 f22 f23 ...
/// @return utility::vector0< std::map<string field, string value> >
/// @trows _AE_scs_failed_ on unexpeceted formating
utility::vector0< std::map<string, string> > parse_plain_text_with_columns(string file_name, string legend_prefix="# ", char legend_separator=' ', string data_prefix="", char data_separator=' ')
{
	utility::vector0< std::map<string, string> > result;
	utility::vector0<string> legend;

	std::ifstream f(file_name);

	string line;
	while( std::getline(f, line) ) {
		if( utility::startswith(line, legend_prefix) ) {
			if( legend_separator == ' ') legend = utility::split_whitespace( line.substr( legend_prefix.size() ) );
			else legend = utility::string_split( line.substr( legend_prefix.size() ), legend_separator);

			for(auto & l : legend) {
				l = utility::strip(l, ' ');
				l = utility::replace_in(l, ". ", ".");
				l = utility::replace_in(l, " ", "-");
			}

			TR.Trace << "Found legend: " << line << std::endl << "Legend: " << legend << std::endl;
		}
		else {
			if( utility::startswith(line, data_prefix) ) {
				if( !legend.size() ) throw _AE_scs_failed_("File: " + file_name + " is missing legend!");
				else {
					//TR.Trace << "Got line:" << line << std::endl;

					utility::vector0< string > parts;

					//auto parts( utility::string_split(line, '\t') );
					if( data_separator == ' ') parts = utility::split_whitespace(line);
					else parts = utility::string_split(line, data_separator);

					//TR.Trace << "Line parts.size: " << parts.size() << " legend size:" << legend.size() << std::endl;

					if( parts.size() == legend.size() ) {
						std::map<string, string> fields;  int i=0;
						for(auto e : parts) fields[ legend[i++] ] = e;
						result.push_back(fields);
					} else throw _AE_scs_failed_("Number of fileds does not match legend!\nFile: "+file_name+"\nLine: "+line);
				}
			}
		}
	}

	return result;
}


/// @details Parse output of NCBI Blast+
///          Expected input:
///          # BLASTP 2.2.31+
///          # Query: h1.fasta
///          # Database: ./../../tools/antibody/blast_database/database.H1.10
///          # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
///          # 423 hits found
///          h1.fasta	pdb4d9r_chothia.pdb	100.00	10	0	0	1	10	1	10	6e-04	22.7
///          h1.fasta	pdb4d9q_chothia.pdb	100.00	10	0	0	1	10	1	10	6e-04	22.7
///          ...
/// @return utility::vector0< std::map<string field, string value> >
SCS_ResultsVector parse_blastp_output(string file_name, string query, std::map< string, std::map<string, string> > db)
{
	SCS_ResultsVector results;

	auto lines( parse_plain_text_with_columns(file_name, "# Fields: ", ',', query, '\t') );

	for(auto fields : lines ) {
		// Converting text based results filed to SCS_BlastMetric
		SCS_BlastResultOP r(new SCS_BlastResult);
		r->pdb = fields["subject-id"].substr(3,4);  // pdb2adf_chothia.pdb → 2adf

		r->alignment_length = utility::string2int( fields.at("alignment-length") );

		r->identity  = utility::string2Real( fields.at("%-identity") );
		r->bit_score = utility::string2Real( fields.at("bit-score") );

		r->bio_type      = db.at(r->pdb)["bio_type"];
		r->light_type    = db.at(r->pdb)["light_type"];
		r->struct_source = db.at(r->pdb)["struct_source"];

		r->h1 = db[r->pdb]["h1"];  r->h2 = db[r->pdb]["h2"];  r->h3 = db[r->pdb]["h3"];  r->frh = db[r->pdb]["frh"];
		r->l1 = db[r->pdb]["l1"];  r->l2 = db[r->pdb]["l2"];  r->l3 = db[r->pdb]["l3"];  r->frl = db[r->pdb]["frl"];

		//r->sequence = db[r->pdb][results.sequence];
		//TR << "PDB:" << r->pdb << std::endl;

		results.push_back(r);
	}

	return results;
}


void SCS_BlastPlus::init_from_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR << "Setting output prefix from command line-options to: " << TR.bgWhite << TR.Black << option[OptionKeys::antibody::prefix]() << std::endl;
	set_output_prefix( option[OptionKeys::antibody::prefix]() );

	TR << "Setting antibody grafting database path from command line-options to: " << TR.bgWhite << TR.Black << option[OptionKeys::antibody::grafting_database]() << std::endl;
	set_database_path( option[OptionKeys::antibody::grafting_database]() );

	TR << "Setting path to NCBI-Blast+ executable from command line-options to: " << TR.bgWhite << TR.Black << option[OptionKeys::antibody::blastp]() << std::endl;
	set_blastp_executable( option[OptionKeys::antibody::blastp]() );
}


SCS_ResultsOP SCS_BlastPlus::raw_select(AntibodySequence const &A) //, AntibodyNumbering const &N)
{
	//AntibodyFramework const frh( A.heavy_framework() );
	//AntibodyFramework const frl( A.light_framework() );

	FRH_FRL fr( calculate_frh_frl(A) );

	string frh = fr.frh1 + fr.frh2 + fr.frh3 + fr.frh4;
	string frl = fr.frl1 + fr.frl2 + fr.frl3 + fr.frl4;

	string heavy_orientation = fr.frh1 + A.h1_sequence() + fr.frh2 + A.h2_sequence() + fr.frh3 + A.h3_sequence() + fr.frh4;
	string light_orientation = fr.frl1 + A.l1_sequence() + fr.frl2 + A.l2_sequence() + fr.frl3 + A.l3_sequence() + fr.frl4;

	string orientation = heavy_orientation + light_orientation; // Note: this is in reverse order if compared to Python code!!!

	TR.Trace << "frh:" << frh << std::endl;
	TR.Trace << "frl:" << frl << std::endl;

	// _framework_names_ = ['FRL', 'FRH', 'light', 'heavy', 'L1', 'L2', 'L3', 'H1', 'H2', 'H3', 'light_heavy']


	SCS_ResultsOP results(new SCS_Results);
	//results->antibody_sequence = A;

	struct {
		string name, sequence;
		SCS_ResultsVector &results;
	} J[] {
		{"frh", frh, results->frh}, {"h1", A.h1_sequence(), results->h1}, {"h2", A.h2_sequence(), results->h2}, {"h3", A.h3_sequence(), results->h3},
		{"frl", frl, results->frl}, {"l1", A.l1_sequence(), results->l1}, {"l2", A.l2_sequence(), results->l2}, {"l3", A.l3_sequence(), results->l3},
	    {"orientation", orientation, results->orientation},
	};

	string blast_database( database_ + "/blast_database/database" );

	utility::vector0< std::map<string, string> > antibody_info_lines( parse_plain_text_with_columns(database_ + "/info/antibody.info") );
	std::map< string, std::map<string, string> > antibody_info_db;
	for(auto fields : antibody_info_lines) antibody_info_db[fields["pdb"]] = fields;

	//TR << antibody_info_db;

	for(auto &j : J) {
		string fasta_file_name(prefix_ + j.name + ".fasta");
		string align_file_name(prefix_ + j.name + ".align");

		//TR << "Working on: " << j.name << std::endl;

		std::ofstream f(fasta_file_name);  f << "> " << fasta_file_name << '\n' << j.sequence << '\n';  f.close();

		// len_cdr = len(cdr_query[k])
		// if k.count('FR') or k.count('heavy') or k.count('light'): db = blast_database + '/database.%s' % k
		// else: db = blast_database + '/database.%s.%s' % (k, len_cdr)

		//string db_suffix = (j.name.find("fr") < string::npos  or  j.name.find("heavy") < string::npos  or j.name.find("light") < string::npos) ? j.name : j.name + '.' + ( std::stringstream() << j.sequence.size() ).str();
		std::stringstream sq;  sq << j.sequence.size();

		string db_suffix = j.name + '.' + sq.str();

		if (j.name.find("fr") < string::npos) db_suffix = j.name;

		std::transform(db_suffix.begin(), db_suffix.end(), db_suffix.begin(), toupper);

		if( j.name.find("orientation")  < string::npos ) db_suffix = "light_heavy";    // or  j.name.find("heavy") < string::npos  or j.name.find("light") < string::npos

		// wordsize, e_value, matrix = (2, 0.00001, 'BLOSUM62') if k.count('FR') or k.count('heavy') or k.count('light') else (2, 2000, 'PAM30')
		//string extra = (j.name.find("fr") < string::npos  or  j.name.find("heavy") < string::npos  or  j.name.find("light") < string::npos) ? "-evalue 0.00001 -matrix BLOSUM62" : "-evalue 2000 -matrix PAM30";
		string extra = (j.name.find("fr") < string::npos  or  \
						j.name.find("orientation") < string::npos or \
						j.name.find("heavy") < string::npos  or  j.name.find("light") < string::npos) ? "-evalue 0.00001 -matrix BLOSUM62" : "-evalue 2000 -matrix PAM30";

		//string command_line ("cd "+working_dir+" && "+blast+" -db "+blast_database+'.'+db_suffix+" -query "+fasta_file_name+" -out "+align_file_name+" -word_size 2 -outfmt 7 -max_target_seqs 1024 "+extra);
		string command_line (blastp_+" -db "+blast_database+'.'+db_suffix+" -query "+fasta_file_name+" -out "+align_file_name+" -word_size 2 -outfmt 7 -max_target_seqs 1024 "+extra);

		basic::execute("Running blast+ for " + j.name + "...", command_line);

		//j.results.name = j.name;
		//j.results.sequence = j.sequence;
		j.results = parse_blastp_output(align_file_name, fasta_file_name, antibody_info_db);
	}

	return results;
}


void trim_framework(AntibodySequence const &A, AntibodyFramework &heavy_fr, AntibodyFramework &light_fr)
{
	if( heavy_fr.fr1_end - heavy_fr.fr1_begin > 25 ) heavy_fr.fr1_begin = heavy_fr.fr1_end - heavy_fr.fr1_begin - 25;
	heavy_fr.fr4_begin = A.heavy.cdr3.end;  heavy_fr.fr4_end = A.heavy.cdr3.end + 12; // FR_H4 = heavy_chain[H3_end + 1: H3_end + 1 + 12]

	if( heavy_fr.fr4_end > A.heavy.sequence.size() ) heavy_fr.fr4_end = A.heavy.sequence.size();

	heavy_fr.update_sequences(A.heavy.sequence);

	//TR << "heavy_fr heavy_fr.fr4_begin:" << heavy_fr.fr4_begin << " heavy_fr.fr4_end:" << heavy_fr.fr4_end << " frh4:" << heavy_fr.fr4 << ":" << heavy_fr.fr4.size() << std::endl;

	if( light_fr.fr1_end - light_fr.fr1_begin > 23 ) light_fr.fr1_begin = light_fr.fr1_end - light_fr.fr1_begin - 23;
	light_fr.fr2_begin = A.light.cdr1.end;  light_fr.fr2_end = A.light.cdr1.end + 15;  // FR_L2 = light_chain[ L1_end + 1  :  L1_end + 1 + 15]
	light_fr.fr4_begin = A.light.cdr3.end;  light_fr.fr4_end = A.light.cdr3.end + 12;  // FR_L4 = light_chain[ L3_end + 1  :  L3_end + 1 + 12 ]

	if( light_fr.fr4_end > A.light.sequence.size() ) light_fr.fr4_end = A.light.sequence.size();

	light_fr.update_sequences(A.light.sequence);  //TR << "Light framework: " << light_fr << std::endl;
}



FRH_FRL calculate_frh_frl(AntibodySequence const &A)
{
	AntibodyFramework heavy_fr = A.heavy_framework();
	AntibodyFramework light_fr = A.light_framework();

	trim_framework(A, heavy_fr, light_fr);

	AntibodyNumbering an( Chothia_Numberer().number(A, heavy_fr, light_fr) );

	FRH_FRL r;

	struct {
		string sequence;
		AntibodyChainNumbering::NumberingVector numbering;
		uint begin, end;
		string &dest;
	} tasks[] {
		{heavy_fr.fr1, an.heavy.fr1,  10,  26, r.frh1},
		{heavy_fr.fr2, an.heavy.fr2,  36,  40, r.frh2},  {heavy_fr.fr2, an.heavy.fr2, 46, 50, r.frh2},
		{heavy_fr.fr3, an.heavy.fr3,  66,  95, r.frh3},
		{heavy_fr.fr4, an.heavy.fr4, 103, 110, r.frh4},

		{light_fr.fr1, an.light.fr1,  10,  24, r.frl1},
		{light_fr.fr2, an.light.fr2,  35,  39, r.frl2},  {light_fr.fr2, an.light.fr2, 45, 50, r.frl2},
		{light_fr.fr3, an.light.fr3,  57,  67, r.frl3},  {light_fr.fr3, an.light.fr3, 71, 89, r.frl3},
		{light_fr.fr4, an.light.fr4,  98, 105, r.frl4},
	};

	for(auto &t : tasks ) {
		auto s( t.sequence.begin() );  auto n( t.numbering.begin() );

		for(; s!=t.sequence.end(); ++s, ++n) {
			uint i = utility::string2Size( *n );
			if( i >= t.begin and i < t.end ) t.dest.push_back( *s );
		}
	}

	//r.frh = r.frh1 + r.frh2 + r.frh3 + r.frh4;
	//r.frl = r.frl1 + r.frl2 + r.frl3 + r.frl4;

	return r;
}


} // namespace grafting
} // namespace antibody
} // namespace protocols


#endif // __ANTIBODY_GRAFTING__

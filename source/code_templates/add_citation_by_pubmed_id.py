#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  code_templates/add_citation_by_pubmed_id.py
## @brief  Script for adding a citation to the Rosetta citations list in the Rosetta database.
## @details  Requires the pubmed-lookup Python package (pip install pubmed-lookup).
## @author  Vikram K. Mulligan (vmulligan@flatironinstitute.org).

# Enure that we have the pubmed_lookup module:
try:
    from pubmed_lookup import Publication, PubMedLookup
    print("Successfully imported pubmed-lookup package.")
except ModuleNotFoundError:
    print('The pubmed-lookup package was not found.  Please install it with "pip install pubmed-lookup" before running this script.')
    exit()

# Other imports:
import argparse
import sys
import os

## @brief A class for storing an author.
class Author:
    def __init__( self, author_string ) :
        #print( "Processing '" + author_string + "'.")
        names = author_string.split( " " )
        assert len(names) > 0
        self.surname = ""
        if len(names) > 1 :
            self.initials = names[len(names) - 1]
        else :
            self.initials = ""
        for i in range( len(names) - 1 ) :
            self.surname += names[i]
        #print ("surname='" + self.surname + "'\tinitials='" + self.initials + "'")

## @brief Given a list of authors, separate into a list
def get_authors_as_list( author_list ):
    returnlist = []
    for name in author_list :
        returnlist.append( Author(name) )
    return returnlist


## @brief A class for storing a single citation.
class Citation:
    def __init__(self):
        self.lines = []
        self.year = None
        self.doi = None


    def init_from_lines( self, lines ) :
        #Store all lines:
        self.lines = lines
        #Extract DOI and year

        next_is_doi = False
        next_is_year = False
        self.doi = None
        self.year = None
        for line in lines :
            if next_is_year == True :
                next_is_year = False
                self.year = int(line)
                #print( "Found year " + str(self.year))
                continue

            if next_is_doi == True :
                next_is_doi = False
                self.doi = line.strip("\t\n ")
                #print( "Found doi " + str(self.doi))
                continue

            if line.startswith("    [BEGIN_DOI]") :
                next_is_doi = True
                next_is_year = False
            elif line.startswith("    [BEGIN_YEAR]" ):
                next_is_year = True
                next_is_doi = False
            else :
                next_is_doi = False
                next_is_year = False

    def init_from_pmid( self, pmid ):
        publication = Publication( PubMedLookup( pmid, '' ), resolve_doi=False )
        self.year = int(publication.year)
        self.doi = publication.url.replace( "http://dx.doi.org/", "" )
        self.lines = ["[BEGIN_CITATION]\n", "    [BEGIN_PRIMARY_AUTHORS]\n"]
        #print( publication._author_list )
        authors = get_authors_as_list(publication._author_list)
        assert len(authors) > 0
        self.lines.append( '        "" "' + authors[0].surname + '" "' + authors[0].initials + '"\n'   )
        self.lines.append( '    [END_PRIMARY_AUTHORS]\n' )
        if( len(authors) > 2 ) :
            self.lines.append( '    [BEGIN_COAUTHORS]\n' )
            for i in range(1,len(authors)-1) :
                self.lines.append( '        "" "' + authors[i].surname + '" "' + authors[i].initials + '"\n'   )
            self.lines.append( '    [END_COAUTHORS]\n' )
        if( len(authors) > 1 ) :
            self.lines.append( '    [BEGIN_SENIOR_AUTHORS]\n' )
            self.lines.append( '        "" "' + authors[len(authors)-1].surname + '" "' + authors[len(authors)-1].initials + '"\n'   )
            self.lines.append( '    [END_SENIOR_AUTHORS]\n' )
        self.lines.append( '    [BEGIN_YEAR]\n' )
        self.lines.append( '        ' + str(self.year) + "\n" )
        self.lines.append( '    [END_YEAR]\n' )
        self.lines.append( '    [BEGIN_TITLE]\n' )
        self.lines.append( '        ' + publication.title + "\n" )
        self.lines.append( '    [END_TITLE]\n' )
        self.lines.append( '    [BEGIN_JOURNAL]\n' )
        self.lines.append( '        ' + publication.journal + "\n" )  
        self.lines.append( '    [END_JOURNAL]\n' )
        self.lines.append( '    [BEGIN_VOLUME_ISSUE_PAGES]\n' )
        volstring = publication.volume
        if publication.issue != "" :
            volstring += "(" + publication.issue + ")"
        if publication.pages != "" :
            volstring += ":" + publication.pages
        self.lines.append( '        ' + volstring + "\n" )
        self.lines.append( '    [END_VOLUME_ISSUE_PAGES]\n' )
        self.lines.append( '    [BEGIN_DOI]\n' )
        self.lines.append( '        ' + self.doi + "\n" )
        self.lines.append( '    [END_DOI]\n' )
        self.lines.append( '[END_CITATION]\n' )
        print( 'Successfully found article titled \"' + publication.title + '".' )

## @brief Helper function for sorting:
def sortfunc( element ):
    if element.year == None :
        return 0
    return element.year


## @brief A class for storing the list of citations.
class CitationFileContents:
    def __init__( self, filename ) :
        assert os.path.exists( filename ) #Should already be guaranteed true.
        with open( filename ) as filehandle:
            lines = filehandle.readlines()

        self.comments = []
        self.citations = []
        in_citation = False
        for line in lines :
            if line.startswith("#") :
                self.comments.append(line)
                continue
            if in_citation == False:
                if line.startswith("[BEGIN_CITATION]") :
                    in_citation = True
                    citation_accumulator = [line]
                    continue
            else :
                if line.startswith("[END_CITATION]") :
                    in_citation = False
                    citation_accumulator.append(line)
                    newcitation = Citation()
                    newcitation.init_from_lines(citation_accumulator)
                    self.citations.append(newcitation)
                    continue
                citation_accumulator.append(line)

    ## @brief Does this collection of citations contain one with the specified DOI?
    def has( self, doi ):
        for citation in self.citations :
            if citation.doi == doi :
                return True
        return False

    ## @brief Add a new citation:
    def add( self, newcitation ) :
        self.citations.append(newcitation)

    ## @brief Sort the list:
    def sort(self) :
        self.citations.sort( key=sortfunc )
    
    ## @brief Write the citations:
    def write_citations(self, filename) :
        with open(filename, 'w') as f :
            #Write comment lines first:
            for line in self.comments :
                f.write(line)
            f.write('\n')
            #Write each citation:
            for citation in self.citations:
                for line in citation.lines:
                    f.write(line)
                f.write('\n')


## @brief Get the Rosetta/main/database directory.
## @details Throws if can't be found.
def get_database_dir( options ) :
    if( options.database ) :
        assert os.path.exists( options.database ), 'Error in get_database_dir: Could not find path "' + options.database + '".'
        return options.database + "/"
    defaultpath = "../database/"
    assert os.path.exists(defaultpath), 'Error in get_database_dir: Could not find path "../database/".  Are you running this script from the "Rosetta/main/source" directory?  If not, please provide the "--database" option.'
    return defaultpath

## @brief Get the Rosetta citation list file.
def get_citation_file( database_dir ) :
    assert os.path.exists( database_dir + "citations/rosetta_citations.txt" ), 'Error in get_citation_file: Could not find "' + database_dir + 'citations/rosetta_citations.txt".'
    return database_dir + "citations/rosetta_citations.txt"

## @brief Program execution entry point.
def main(args):
    # Parse options:
    parser = argparse.ArgumentParser()
    parser.add_argument( '--database', help="Path to the Rosetta/main/database/ directory.  Not required if this script is run from the Rosetta/main/source/ directory.")
    parser.add_argument('--pmid', help="The PubMed ID of the citation to add to the dataase.")
    options = parser.parse_args(args=args[1:])

    # Check for required options:
    assert options.pmid != None and options.pmid != "", 'Error in main: You must provide a PubMed ID number with the "--pmid" option.'
    new_pmid = options.pmid

    database_dir = get_database_dir( options )
    citation_file = get_citation_file( database_dir )
    citation_file_contents = CitationFileContents( citation_file )
    new_citation = Citation()
    new_citation.init_from_pmid( new_pmid )
    if not citation_file_contents.has(new_citation.doi):
        print( "Adding PubMed ID " + new_pmid + " (DOI " + new_citation.doi + "). to citations file." )
        citation_file_contents.add(new_citation)
    else :
        print( "The PubMed ID " + new_pmid + " was already found in the citations file (DOI " + new_citation.doi + ").  Only sorting citations file." )
    
    citation_file_contents.sort()
    citation_file_contents.write_citations( citation_file )

    print( 'Wrote sorted citations list to "' + citation_file + '".  Please be sure to commit your changes with "git add ' + citation_file + ' && git commit".' )

if __name__ == "__main__" :
    main(sys.argv)
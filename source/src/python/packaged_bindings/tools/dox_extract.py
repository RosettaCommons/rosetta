#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   dox_extract.py
## @brief  Extract Doxygen comments from C++ files
## @author Brian Weitzner, Sergey Lyskov

import re, sys

class dox_extract:
    def __init__( self ):
        self.is_a_function = re.compile( '[;{}]' )
        self.is_a_blank_line = re.compile( '^\s*$' )
        self.line_begins_with_white_space = re.compile( '^\s' )
        self.contains_end_block_dox_comment = re.compile( '\*\*/\s*' )
        self.contains_start_block_dox_comment = re.compile( '\s*/\*\*' )
        self.contains_three_slash_dox_comment = re.compile( '^\s*///' )

        # add doxygen symbols here as needed
        self.dox_symbols = [
            'brief', 'details', 'param', 'return',
            'fn',
            'file',
            'author',
            'detailed',
            'li',
            'note'
        ]
    def clean_line( self, line ):
        line = re.sub( self.contains_three_slash_dox_comment, '', line )
        line = re.sub( self.contains_start_block_dox_comment, '', line )
        return line.rstrip() + '\n'

    def remove_unnecessary_whitespace( self, lines ):
        there_is_whitespace = False
        begins_with_white_space = self.line_begins_with_white_space.search( lines[ 0 ] )
        if begins_with_white_space is not None:
            there_is_whitespace = True
        while there_is_whitespace:
            for i in range( len( lines ) ):
                lines[ i ] = re.sub( self.line_begins_with_white_space, '', lines[ i ] )
            begins_with_white_space = self.line_begins_with_white_space.search( lines[ 0 ] )
            if begins_with_white_space is not None:
                there_is_whitespace = True
            else:
                there_is_whitespace = False
        return lines

    def remove_dox_symbols( self, lines ):
        for dox_symbol in self.dox_symbols:
            for i in range( len( lines ) ):
                lines[ i ] = re.sub( re.compile( '@' + dox_symbol + '\s*' ), '', lines[ i ] )
        return lines

    def getDoxygenComment( self, file_name, line_number ):
        ''' Search and return doxygen comments at given file/line position and return them
            Note: - line is a number in string form here.
                  - if no appropriate comments could be found (or result is empty) function should return  file_ + ':' + line  string
        '''

        line_number = int( line_number )

        i_say_so = True
        block_comment = False

        file = open( file_name, 'r' )
        lines = file.readlines()
        file.close()

        comment_block = []
        current_line_number = line_number - 1

        while i_say_so:
            current_line_number -= 1

            three_slashes = self.contains_three_slash_dox_comment.search( lines[ current_line_number ] )

            end_block = self.contains_end_block_dox_comment.search( lines[ current_line_number ] )
            start_block = self.contains_start_block_dox_comment.search( lines[ current_line_number ] )

            blank_line = self.is_a_blank_line.search( lines[ current_line_number ] )

            # if we decide to keep blank lines, remove the last condition below and then uncomment the first elif
            if three_slashes is not None or block_comment:
                comment_block.append( self.clean_line( lines[ current_line_number ] ) )
            elif blank_line is not None:
                continue
            elif end_block is not None:
                block_comment = True
            else:
                function = self.is_a_function.search( lines[ current_line_number ] )
                if function is not None:
                     i_say_so = False

            if block_comment and start_block is not None:
                block_comment = False

        comment_block.reverse()
        if len( comment_block ):
            comment_block = self.remove_unnecessary_whitespace( comment_block )
            comment_block = self.remove_dox_symbols( comment_block )
        return ''.join(comment_block).__repr__()[1:-1].replace('"', '\\"') or (file_name + ':' + str(line_number))

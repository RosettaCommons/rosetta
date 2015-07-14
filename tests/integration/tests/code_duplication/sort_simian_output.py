#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @brief  Script to sort Simian output so test results are stable and does not depend on file order
## @author Sergey Lyskov


import re, types


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = 'NT: |'
        for i in dir(self):
            if not i.startswith('__') and not isinstance(getattr(self, i), types.MethodType): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'


with file('log') as f:
    _ = 'ATTENTION: For raw Simian output containing line numbers and block lengths please see raw_simian_log.ignore\n\n' + f.read()
    blocks = re.split('(Found \d+ duplicate lines in the following files:\n)', _)
    #blocks = re.split('Found ', _)


#blocks = [b for b in blocks if b]
#blocks = [ NT(header=blocks[i*2], text=blocks[i*2+1]) for i in range( len(blocks)/2 )]

nb = [];  i=0
while i<len(blocks):
    if blocks[i].startswith('Found ')  and  i+1<len(blocks):
        if blocks[i+1].startswith(' Between lines '): nb.append( NT(header=blocks[i], text=blocks[i+1]) ); i+=1  # also filter empty block, sometime simian produce two 'Found ...' in a row etc.
    elif blocks[i]: nb.append( NT(header=blocks[i], text='') )
    i += 1

blocks = nb


# for b in blocks:
#     if not b.text: print '______EMPTY:', b
#     #else: print b

def replace_line_numbers(line):
    ''' Expecting: Between lines 120 and 134 in ROSETTA_MAIN/source/src/utility/string_util.cc
        Output:    Between lines XXXX and XXXX in ROSETTA_MAIN/source/src/utility/string_util.cc
    '''
    a = line.split(' ')
    if len(a) > 5:
        a[3] = 'XXXX'
        a[5] = 'XXXX'
    return ' '.join(a)


def sort_block(block):
    ''' Sort block of text and return sorting (key, new-text) for this block '''

    #text = re.split('Found \d+ duplicate lines in the following files:\n', text)[0]  # removing 'Found...' header if any
    #print 'text:', text

    if block.text.startswith(' Between lines '):
        lines, block_key = block.text.split('\n'), ''
        lines = [l for l in lines if l]
        if lines:
            lines.sort(key=lambda l: l.split()[-1])
            block_key =  ' '.join( [l.split()[-1] for l in lines] ) + '\n'.join(lines)
        return block_key, block.header+'\n'.join( [ replace_line_numbers(l) for l in lines] )

    else: return '', block.header+block.text  # This is probably Simian header, we ignore it


blocks = [ sort_block(b) for b in blocks]

blocks.sort(key=lambda k_l: k_l[0])

#for b in blocks: print '_____', b, '_____\n'

blocks = [ b for k, b in blocks]

#text = '\n\nFound duplicated code in the following files:\n'.join(blocks) + '\n'
text = '\n\n'.join(blocks) + '\n'

with file('log', 'w') as f: f.write(text)

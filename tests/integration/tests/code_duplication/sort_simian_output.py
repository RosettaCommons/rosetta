#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @brief  Script to sort Simian output so test results are stable and does not depend on file order
## @author Sergey Lyskov


import re

with file('log') as f:
    _ = f.read()
    blocks = re.split('Found \d+ duplicate lines in the following files:\n', _)

blocks = [b for b in blocks if b]

def sort_block(text):
    ''' Sort block of text and return sorting (key, new-text) for this block '''

    if text.startswith(' Between lines '):
        lines, block_key = text.split('\n'), ''
        lines = [l for l in lines if l]
        if lines:
            lines.sort(key=lambda l: l.split()[-1])
            block_key =  ' '.join( [l.split()[-1] for l in lines] ) + '\n'.join(lines)
        return block_key, '\n'.join(lines)

    else: return '', text  # This is probably Simian header, we ignore it


blocks = [ sort_block(b) for b in blocks]

blocks.sort(key=lambda k_l: k_l[0])

#for b in blocks: print '_____', b, '_____\n'

blocks = [ b for k, b in blocks]

text = '\n\nFound duplicated code in the following files:\n'.join(blocks) + '\n'

with file('log', 'w') as f: f.write(text)

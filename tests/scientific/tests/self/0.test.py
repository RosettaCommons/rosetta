#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  self/0.test.py
## @brief scientific self-test script
## @author Sergey Lyskov

import json

import benchmark
from benchmark import _StateKey_, _ResultsKey_, _LogKey_, _S_passed_, _S_failed_, _S_build_failed_, _S_script_failed_, _multi_step_result_

config = benchmark.config()

_index_html_template_ = '''\
<html>
<body>
<h1>Scientific test: self</h1>
test!
</body>
</html>
'''

working_dir = config['working_dir']

with open(f'{working_dir}/index.html', 'w') as f: f.write(_index_html_template_)

with open(_multi_step_result_, 'w') as f:
    r = {
        _StateKey_  : _S_passed_,
        _ResultsKey_ : {},
        _LogKey_ : 'Done!',
    }

    json.dump(r, f, sort_keys=True, indent=2)

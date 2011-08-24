PATH_TO_ROOT = '../../'

import os, string, sys
sys.path.append(PATH_TO_ROOT + 'script/')
sys.path.append(PATH_TO_ROOT + 'script/novice/')
sys.path.append(PATH_TO_ROOT + 'script/novice/components/')
import components

components.make_component(PATH_TO_ROOT, 'resources', False, False)
components.make_component(PATH_TO_ROOT, 'database', False, False)

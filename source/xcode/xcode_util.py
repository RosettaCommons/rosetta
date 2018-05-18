PATH_TO_ROOT = '../../'

import hashlib, os, string, sys
#sys.path.append(PATH_TO_ROOT + 'script/')
import build_util


#HOW TO UPDATE XCODE PROJECT WITH NEW LIBRARIES
#==============================================

#NOTE: This is ONLY for adding new Rosetta libraries - aka protocols.4.src.settings. 
#If you add, move, delete files and directories, no change needs to happen for you and re-running make_project.py will 
#update xcode with your new code. 

#You will need to first add the new library to your Xcode project.  Automation will happen at some point, but is not trivial.
#
# Steps:
#        1) Double Click on Rosetta.  See the Target list.  Right click on the closest library to yours and hit duplicate
#        2) Click on the target like you do on Mac finder window and rename.  Drag the target to where it should be.
#        3) Double click on your target.  To the right select build phases.  Add your dependencies that are missing to target dependancies
#             and link binary to library
#        4) Right click sources, hit create group.  Name the group the library.  Drag the group to the appropriate place.
#        5) Now, you have to open the src.setting file and find a directory within it.  Right click on your new library in sources.
#            Add the file into the group.  Select options (lower left).  Make sure the target is selected and new group is checked
#        6) Finally, meander back to targets. Go to build phases.  If you see a compile sources section, click the upper right x to clear it.
#           Next, select the + button on the upper left.  Select 'new compile sources' tag.  Under the section, click the + button.
#           Select a .cc file corresponding in the directory named the same as this target.
#        7) Run make_project.py - make sure everything works.  Copy your .pbxproj file into .pbxproj.template.
#        8) Use the grep commands below to get the keys.  Source key is right before it says /* source * /

# group keys grep:
# grep -A 50 'BE8A17150CA8365000D67A6F \/\* Sources \*\/ =' Rosetta.xcodeproj/project.pbxproj
# protocols sources keys grep:
# grep -A 2000 '^\/\* Begin PBXNativeTarget' Rosetta.xcodeproj/project.pbxproj | grep -A 6 '\/\* protocols' | grep -E '(\/\* Sources|\/\* protocols)'


PROJECT_KEYS = {  #        group                       sources
    'basic'       :      ('D1FC4A1A13687F50006C102D', 'D1FC49F813687F42006C102D'),
    'utility'     :      ('BEA64D500DB2E11C003495F8', 'BEDC3FB80CA84F63000FAD97'),
    'numeric'     :      ('BEA64B640DB2E09B003495F8', 'BEDC40E50CA8592F000FAD97'),
    'ObjexxFCL'   :      ('BE8A1AF00CA8367600D67A6F', 'BEDC3F960CA84E43000FAD97'),
    'core.1'      :      ('BEA646040DB2DF4E003495F8', 'BEDC400C0CA8522A000FAD97'),
    'core.2'      :      ('D1FC4A1B13687F65006C102D', 'D139B5B113687DEE0030829C'),
    'core.3'      :      ('D1FC4A1C13687F7A006C102D', 'D139B7E413687E0F0030829C'),
    'core.4'      :      ('D1FC4A1D13687F7E006C102D', 'D139BA0C13687E180030829C'),
    'core.5'      :      ('D1FC4A1E13687F83006C102D', 'D139BC3413687E200030829C'),
    'protocols.1'   :    ('3D5BA02A1492D65F00524FBB', 'D1FD035D0E2851B800B3C702'),
    'protocols_a.2' :    ('3D5BA02E1492D69900524FBB', '3B945C85148C9D5200C2FBE5'),
    'protocols_b.2' :    ('3D5BA0311492D6BB00524FBB', '3D5BA0861492D95800524FBB'),
    'protocols.3'   :    ('3D5BA02B1492D66A00524FBB', '3B94613C148CA62A00C2FBE5'),
    'protocols.4'   :    ('7092E3B920AA0DC000495DD9', '7092E45120AA126100495DD9'),
    'protocols_a.5' :    ('3D5BA0301492D6AD00524FBB', '3B9464A6148CBD9D00C2FBE5'),
    'protocols_b.5' :    ('3D5BA0331492D6D200524FBB', '3B9464EE148CBE0A00C2FBE5'),
    'protocols_c.5' :    ('3D5BA0351492D6E700524FBB', '3B946529148CBE0C00C2FBE5'),
    'protocols_d.5' :    ('3D5BA0371492D6FC00524FBB', '3D5BA06A1492D94500524FBB'),
    'protocols_e.5' :    ('3D5BA0391492D71200524FBB', '3B94656E148CBF6800C2FBE5'),
    'protocols_f.5' :    ('3D5BA03B1492D72800524FBB', '3B9465A0148CBF6E00C2FBE5'),
    'protocols_a.6' :    ('7092E67220AA3D9C00495DD9', '7092EB3420AA59CD00495DD9'),
    'protocols_b.6' :    ('7092E67320AA3DB600495DD9', '7092EB3620AA5A7700495DD9'),
    'protocols_c.6' :    ('7092E88A20AA575A00495DD9', '7092EB3820AA5B8500495DD9'),
    'protocols_d.6' :    ('7092E88B20AA576500495DD9', '7092EB3A20AA5BC700495DD9'),
    'protocols_e.6' :    ('7092E8F120AA57B600495DD9', '7092EB3C20AA5BF400495DD9'),
    'protocols.7'   :    ('3D5BA02D1492D68900524FBB', '3B94661A148CBFBF00C2FBE5'),
    'protocols.8':       ('7092EA7B20AA58BA00495DD9', '7092EB3E20AA5C2700495DD9'),

    #'interactive' :      ('BE8A199B0CA8367500D67A6F', 'BEDC41120CA85A65000FAD97'),
    #'game'        :      ('D1CA2C0E110BFEC30085C48B', 'D1CA2C08110BFE910085C48B'),
    #'interactive.test' : ('D18BC31A0FCDF26D000AE673', 'D18BC3030FCDF1C6000AE673'),
    #'game.test'        : ('D1CA2C26110BFFE40085C48B', 'D1CA2C10110BFEDC0085C48B'),
    'devel'       :      ('3D3A9A1413ECE8330081959A', '3DDF62B913EE13B800401374'),
}

TEST_KEYS = {  #        group                       sources
    'basic.test'       :      ('D1FC4A1A13687F50006C102D', 'D1FC49F813687F42006C102D'),
    'apps.test'        :      ('BEA64D500DB2E11C003495F8', 'BEDC3FB80CA84F63000FAD97'),
    'core.test'        :      ('BEA64B640DB2E09B003495F8', 'BEDC40E50CA8592F000FAD97'),
    'demo.test'        :      ('BE8A1AF00CA8367600D67A6F', 'BEDC3F960CA84E43000FAD97'),
    'devel.test'       :      ('BEA646040DB2DF4E003495F8', 'BEDC400C0CA8522A000FAD97'),
    'numeric.test'     :      ('D1FC4A1B13687F65006C102D', 'D139B5B113687DEE0030829C'),
    'protocols.test'   :      ('D1FC4A1C13687F7A006C102D', 'D139B7E413687E0F0030829C'),
    'utility.test'     :      ('D1FC4A1D13687F7E006C102D', 'D139BA0C13687E180030829C'),
}


def find_line(s, lines, start):
    for ii in range(start, len(lines)):
        if lines[ii].startswith(s):
            return ii
    raise RuntimeError('Can\'t find line starting with: "' + s + '".')

def get_key(line):
    key = line.split()[0]
    if len(key) != 24:
        raise RuntimeError('Bad key: "' + key + '".')
    return key

def remove_groups_and_file_refs(group_key, group_name, group_lines, file_ref_lines):
    is_group = False
    try:
        group_begin_ln = find_line('\t\t' + group_key + ' /* ' + group_name + ' */ = {', group_lines, 0)
        is_group = True
    except:
        pass

    if is_group:
        children_begin_line = find_line('\t\t\tchildren = (', group_lines, group_begin_ln)
        children_end_line = find_line('\t\t\t);', group_lines, children_begin_line)
        group_end_ln = find_line('\t\t};', group_lines, children_end_line)

        children_lines = group_lines[children_begin_line + 1:children_end_line]
        del group_lines[group_begin_ln:group_end_ln + 1]

        for ln in children_lines:
            key = get_key(ln)
            name = ln.split()[2]
            remove_groups_and_file_refs(key, name, group_lines, file_ref_lines);
    else:
        file_ref_ln = find_line('\t\t' + group_key + ' /* ' + group_name + ' */ = {', file_ref_lines, 0)
        del file_ref_lines[file_ref_ln]

def remove_build_files_and_sources(source_files_lines, build_file_lines):
    for l in source_files_lines:
        key = get_key(l)

        bl = find_line('\t\t' + key, build_file_lines, 0)
        del build_file_lines[bl]

    del source_files_lines[:]

def genkey(s):
    return hashlib.md5(s.encode('utf-8')).hexdigest().upper()[0:24]

def to_dir_tpl(d):
    dl = d.strip('/').split('/')
    if dl == ['']:
        dl = []
    return tuple(dl)

def from_dir_tpl(dl):
    # import string
    return '/'.join(dl) + '/'

def make_groups_and_file_refs(root_group_key, project, project_path, project_files):
    node_map = {}

    for d, fs in sorted(project_files):
        d_dir_tpl = to_dir_tpl(d)
        d_dir_top = d_dir_tpl[0]
        d_dir_tpl = d_dir_tpl[1:]
        for ii in range(len(d_dir_tpl) + 1):
            part_dir_tpl = d_dir_tpl[0:ii]
            if not part_dir_tpl in node_map:
                full_path = d_dir_top + '/' + from_dir_tpl(part_dir_tpl)

                new_group = {}
                if len(part_dir_tpl) == 0:
                    new_group['name'] = project
                    new_group['key'] = root_group_key
                    parent = None
                else:
                    new_group['name'] = part_dir_tpl[-1]
                    new_group['key'] = genkey(full_path + ':group:' + project)
                    parent = node_map[part_dir_tpl[0:-1]]
                new_group['path'] = full_path
                new_group['children_dir'] = []
                new_group['children_file'] = []

                if parent != None:
                    parent['children_dir'].append(new_group)

                node_map[part_dir_tpl] = new_group

        group_node = node_map[d_dir_tpl]

        for f in sorted(fs):
            full_path = d_dir_top + '/' + from_dir_tpl(d_dir_tpl) + f

            new_file_ref = {}
            new_file_ref['name'] = f
            new_file_ref['path'] = full_path
            new_file_ref['key_build'] = genkey(full_path + ':build:' + project)
            new_file_ref['key_ref'] = genkey(full_path + ':file_ref:' + project)
            if f.endswith('.cc'):
                new_file_ref['sourcetype'] = 'sourcecode.cpp.cpp'
                new_file_ref['is_source'] = True
            elif f.endswith('.hh'):
                new_file_ref['sourcetype'] = 'sourcecode.cpp.h'
                new_file_ref['is_source'] = False
            elif f.endswith('.c'):
                new_file_ref['sourcetype'] = 'sourcecode.c.c'
                new_file_ref['is_source'] = True
            elif f.endswith('.h'):
                new_file_ref['sourcetype'] = 'sourcecode.c.h'
                new_file_ref['is_source'] = False
            else:
                raise RuntimeError('Unknown extension: ' + f)

            group_node['children_file'].append(new_file_ref)

    return node_map[()]

def add_new_lines(path_to_files, group, group_lines, file_ref_lines, build_file_lines, source_files_lines):
    group_lines.append('\t\t' + group['key'] + ' /* ' + group['name'] + ' */ = {\n')
    group_lines.append('\t\t\tisa = PBXGroup;\n')
    group_lines.append('\t\t\tchildren = (\n')
    for child in group['children_dir']:
        group_lines.append('\t\t\t\t' + child['key'] + ' /* ' + child['name'] + ' */,\n')
    for child in group['children_file']:
        group_lines.append('\t\t\t\t' + child['key_ref'] + ' /* ' + child['name'] + ' */,\n')
    group_lines.append('\t\t\t);\n')
    group_lines.append('\t\t\tname = ' + group['name'] + ';\n')
    group_lines.append('\t\t\tpath = ' + path_to_files + group['path'] + ';\n')
    group_lines.append('\t\t\tsourceTree = SOURCE_ROOT;\n')
    group_lines.append('\t\t};\n')

    for child in group['children_file']:
        file_ref_lines.append('\t\t' + child['key_ref'] + ' /* ' + child['name'] + ' */ = {isa = PBXFileReference; fileEncoding = 30; lastKnownFileType = ' + child['sourcetype'] + '; name = ' + child['name'] + '; path = ' + path_to_files + child['path'] + '; sourceTree = SOURCE_ROOT; };\n')

        if child['is_source']:
            build_file_lines.append('\t\t' + child['key_build'] + ' /* ' + child['name'] + ' in Sources */ = {isa = PBXBuildFile; fileRef = ' + child['key_ref'] + ' /* ' + child['name'] + ' */; };\n')
            source_files_lines.append('\t\t\t\t' + child['key_build'] + ' /* ' + child['name'] + ' in Sources */,\n')

    for child in group['children_dir']:
        add_new_lines(path_to_files, child, group_lines, file_ref_lines, build_file_lines, source_files_lines)

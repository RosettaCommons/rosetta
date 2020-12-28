import os, shutil, distutils.dir_util

def move(source, destination):
    ''' move file/directory
    '''
    assert source[0] != '/'
    assert destination[0] != '/'

    assert '..' not in source
    assert '..' not in destination

    print( f'moving {source} --> {destination}' )
    shutil.move(source, destination)


def copy(source, destination):
    ''' move file/directory
    '''
    assert source[0] != '/'
    assert destination[0] != '/'

    assert '..' not in source
    assert '..' not in destination

    print( f'copying {source} --> {destination}' )

    if os.path.isfile(source) and not os.path.isdir(destination): shutil.copyfile(source, destination)
    else: distutils.dir_util.copy_tree(source, destination)


def mkdir(path):
    ''' create directory
    '''
    assert ' ' not in path
    assert '..' not in path
    assert path[0] != '/'

    print( f'mkdir {path!r}' )

    path = os.path.realpath( os.getcwd() + '/' + path )
    assert path.startswith( os.getcwd() + '/' )

    if not os.path.isdir(path): os.makedirs(path)


def register():
    return dict(
        move = move,
        copy = copy,
        mkdir = mkdir,
    )

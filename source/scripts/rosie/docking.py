
import os, os.path


def create_ensemble_lists(*names):
    ''' create docking-ensemble-list files by collecting all name*.pdb files and storing them in name-ensemble.txt files
    '''
    for name in names:
        for c in ' /': assert c not in name

        with open(f'{name}-ensemble.txt', 'w') as f:
            for file in os.listdir():
                if os.path.isfile(file) \
                   and file.startswith(name) \
                   and file.endswith('.pdb') \
                   and not file.endswith('_low.pdb') \
                   and not file.endswith('_last.pdb') \
                   and len(file) > len(name) + len('.pdb'): f.write(file + '\n')

def register():
    return {
        'docking.create-ensemble-lists' : create_ensemble_lists,
    }

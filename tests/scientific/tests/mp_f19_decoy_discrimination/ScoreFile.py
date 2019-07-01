#!/usr/bin/env python2.7
import re
import string
import operator
import bz2
import gzip

open_func = {
                'zip': gzip.GzipFile,
                'bz2': bz2.BZ2File,
                None : open
}

def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

class BadRegister( Exception ):
    pass

class ScoreData:
  def __init__(self, verbose=True, compression=None):
    self.data = {}
    self.compression = compression
    self.verbose = verbose
    self.headers = [ ]

  def add_register(self,score_register):
    get_score = re.compile("^SCORE:\s*(.+)")
    skipped = []
    if not score_register.startswith( 'SCORE' ):
        raise BadRegister( score_register )
    else:
      scores = get_score.findall(score_register.strip())[0].split()
      if len(self.headers) != len(scores):
          raise BadRegister( scores[-1] )
      data = dict(zip(self.headers,scores))
      self.data[scores[-1]] = data


  def write_score_file(self,filename = None, handle = None, sorted_by = "score"):
    assert sorted_by in self.headers
    out = None
    if filename:
      out = open_func[ self.compression ](filename,'w')
    else:
      out = handle
    assert(out)
    out.write("SCORE: " + string.join(self.headers,'\t') + '\n')
    for  keys, values in sorted(self.data.iteritems(),key= lambda x: x[1][sorted_by]):
      scores = [values[header] for header in self.headers]
      out.write("SCORE: " + string.join(map(lambda x: str(x), scores),'\t') + '\n')

  def get_top(self, sorted_by = 'score', top = 1, obj_return = 'register', key='lower'):
    assert sorted_by in self.headers
    assert obj_return in self.headers or obj_return == 'register'
    assert top < len(self.data)
    assert key == 'lower' or key == 'upper'
    rev = True if key == 'lower' else False
    data = []
    for  keys, values in sorted(self.data.iteritems(),key= lambda x: float(x[1][sorted_by]) , reverse = rev)[0:top+1]:
      if obj_return == 'register':
        scores = [values[header] for header in self.headers]
        data.append("SCORE: " + string.join(map(lambda x:str(x),  scores),'\t') + '\n')
      elif obj_return == 'data':
        scores = [values[header] for header in self.headers]
        data.append(scores)
      else:
        data.append(values[obj_return])
    return data

  def get_col(self, header, rtype=float):
    values = []
    for id in self.data.keys():
      try:
        values.append(rtype(self.data[id][header]))
      except:
        if self.verbose: print(id + ' does not seem to have a ' + header + ' field!')
    return values

  def append_col(self, header, default_value = '0.0'):
    if header not in self.headers:
      self.headers.insert(-1,header)
      for values in self.data.itervalues():
        values[header] = default_value
    else:
      if self.verbose: print(header + " already in score file. new column will not be appended")

  def set_value(self, header, id, value, overwrite = False):
    if overwrite:
      self.data[id][header] = value
    elif float(self.data[id][header]) != 0.0:
       pass
    else:
      self.data[id][header] = value

  def get_value(self, id, term):
    return self.data[id][term]


  def read_score(self, filename ):
    input =  open_func[ self.compression ](filename,'r')
    col_headers = input.readline()
    if col_headers.startswith('SEQUENCE'):
      col_headers = input.readline()
    get_col_headers = re.compile("^SCORE:\s*(.+)")
    self.headers = get_col_headers.findall(col_headers.strip())[0].split()
    skipped = []
    for (lnum, line) in enumerate(input):
        if not is_number(line.strip().split()[1]):
            continue
        try:
            self.add_register(line.strip())
        except BadRegister as e:
            skipped.append( str(e) )
    if len(skipped):
            print('Skipped', len(skipped), '/', lnum, 'bad registers')

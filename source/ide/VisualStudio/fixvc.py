import sys

if __name__ == '__main__':
  with open(sys.argv[1]) as solution:
    for line in solution:
      line = line.replace('vcxproj', 'vcproj')
      print line

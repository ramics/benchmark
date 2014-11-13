
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy

def main(filename):

    with open(filename) as f:
        for cut in [2,3]:
            floats = [float(line.split(";")[cut].split("=")[1].strip()) for line in f]
            lists = [floats[i::5] for i in range(4)]
            print lists
            averages = [(sum(i) / float(len(i))) for i in lists]
            temp = averages[0]
            averages[0] = averages[1]
            averages[1] = temp
            stddevs = [numpy.std(i) for i in lists]
            temp = stddevs[0]
            stddevs[0] = stddevs[1]
            stddevs[1] = temp
            for i in range(len(averages)):
                print "%.3f +/- %.3f" % (averages[i], stddevs[i]/2)
            print ""
            f.seek(0)

        

if __name__ == "__main__":
  main(sys.argv[1])
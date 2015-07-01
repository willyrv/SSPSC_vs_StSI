#!/usr/bin/env python

import lib.model
import sys
from lib.model import *
import time
CONFIG_FILE = 'experiment_settings.txt'
OUTPUT_FILE = 'experiment_OUT.txt'

if __name__ == "__main__":
    #OUTPUT_FILE = sys.argv[1]    
    start = time.time()
    p = Parser(CONFIG_FILE)
    options = p.get_array_of_parameters()
    c = Comparator()
    [output_text, observation_text] = c.multiple_experiments(options)
    end = time.time()
    time_elapsed = end - start
    output_text += '\n TIME (in seconds)  ' + time_elapsed.__str__()
    a = open(OUTPUT_FILE, 'w')
    a.write(output_text)
    a.close()
    a = open(OUTPUT_FILE+'_obs', 'w')
    a.write(observation_text)
    a.close()
    print 'Done'


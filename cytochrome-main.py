#This is a main file for the cytochrome calculations. This file simply interact with a user and access chosen subroutines. 
import sys
import os
sys.path.insert(1,os.getcwd()+'\\libraries')

import numpy as np
import cytochrome_lib as cyt
reload(cyt)


config_dict = cyt.read_config_file('config.cyt')
command_dict = cyt.read_config_file('config.cyt',method = 'commands')
#print command_dict
while 1:
    text_input = str.lower(raw_input("Please enter the command (type 'get help' for list of commands:)"))
    print "you have entered:", text_input
    if text_input == 'exit':
        break
    if text_input[0:3] == 'run':
        print 'running:' + text_input[4:]
        try:
            exec(open(text_input[4:]).read(), globals())
        except (ValueError,RuntimeError, TypeError, NameError,SyntaxError, IOError):
            print 'failed to run the script' + text_input[4:]
    elif text_input[0:3] == 'get':
        print 'getting:' + text_input[4:]
        try:
            print(command_dict[text_input[4:]])
        except (ValueError,RuntimeError, TypeError, NameError,SyntaxError, IOError):
            print(command_dict['other'])
    else:
        print(command_dict['other'])
    
    




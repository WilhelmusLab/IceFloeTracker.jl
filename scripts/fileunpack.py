import os, sys

unpack = sys.argv
#dir, _ = os.path.split(unpack[1]) 
f = open(os.path.join(unpack[1], "myfile.txt"), "x") 
f.close()
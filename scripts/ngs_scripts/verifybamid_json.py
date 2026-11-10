import sys
import json

input_file = sys.argv[1]
out_file = sys.argv[2]
dict = {}
with open(input_file) as verfiybamid_file:
	for line in verfiybamid_file:
		if line.startswith("#"):
			continue
		info = line.split()
		sample_name = info[0]
		contaminate = info[6]
		dict[sample_name] = {"contaminate":contaminate}

ofile = open(out_file, 'w')
jstr = json.dumps(dict)
ofile.write(jstr)
ofile.close()

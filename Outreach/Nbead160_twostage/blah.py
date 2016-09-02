#!/usr/bin/env python

line = raw_input()
while line != "Bonds":
	print line
	line = raw_input()

print line
print raw_input()

line = raw_input()
while line:
	spline = line.split()
	if int(spline[2])%10==0:
		line = line[:10] + "2" + line[11:]
	print line
	line = raw_input()

while True:
	print line
	try:
		line = raw_input()
	except:
		break

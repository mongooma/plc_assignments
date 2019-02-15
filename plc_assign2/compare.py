import sys

file1 = sys.argv[1];
file2 = sys.argv[2];

c = 0
with open(file1) as f1, open(file2) as f2:
	c1 = f1.read(1)
	c2 = f2.read(1)
	if (c1 != c2):
		c += 1

print("Different digits: %s" % str(c))




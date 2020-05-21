f = open('all_genes.txt','r')

other = 0
count = 0
for line in f.readlines():
	if line.split()[0] != '1':
		count += 1
	else:
		other += 1
print count, other

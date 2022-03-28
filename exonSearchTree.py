'''
Making an exon search tree which will hopefully speed up the whole get_variants section.
Then again, it might slow it down. Maybe I'll test it both ways.
Or maybe I won't. Who knows what I'll do? I'm a WILD CARD

I tested it on the setup step and it's several times faster (30x speedup in system time).
Next, search...

Another thing to test might be whether it's faster with two chunk layers. But for now
I'm just gonna keep on trucking through the main code.

It's a several-layer tree. There's a chromosome node, which has daughter "chunk" nodes (one chunk for each million
BP in a chromosome), which have leaf nodes for each exon. (An exon that spans multiple chunk nodes will appear
in all, since we'll be searching for single positions.)
'''

class ExonSearchTree:
	'''The root of the tree, so to speak. Everything gets called through this class and propagates "down".'''
	def __init__(self, chunk_size = 1000000):
		self.chr_nodes = {}
		self.chunk_size = chunk_size
		self.total_exons = 0
	
	def print_tree(self):
		'''Used this instead of __str__ for formatting reasons, I guess.''' 
		for node in self.chr_nodes:
			self.chr_nodes[node].print_node()
	
	def add_exon(self, chr, start, end, pos_in_gene, name = "."):
		'''Adds an exon to the tree. Call the one from this class, it'll call the rest appropriately.'''
		chr = str(chr)
		try:
			self.chr_nodes[chr].add_node(start, end, name, pos_in_gene)
		except KeyError:
			self.chr_nodes[chr] = ESTChrNode(chr, self.chunk_size)
			self.chr_nodes[chr].add_node(start, end, name, pos_in_gene)
		self.total_exons += 1
			
	def find_exon(self, chr, pos):
		'''Finds all exons containing a query variant.'''
		chr = str(chr)
		try:
			exon = self.chr_nodes[chr].find_exon(pos)
		except KeyError:
			exon = []
		return exon

class ESTChrNode:
	# A layer of chromosome nodes
	def __init__(self, chr_num, chunk_size):
		self.chr_num = chr_num
		self.chunk_nodes = {}
		self.chunk_size = chunk_size
		self.total_exons = 0
	
	def add_node(self, start, end, name, pos_in_gene):
		start_chunk = start//self.chunk_size
		end_chunk = end//self.chunk_size
		for chunk in range(start_chunk, end_chunk+1):
			try:
				self.chunk_nodes[chunk].add_node(start, end, name, pos_in_gene)
			except KeyError:
				self.chunk_nodes[chunk] = ESTChunkNode(self.chr_num, chunk)
				self.chunk_nodes[chunk].add_node(start, end, name, pos_in_gene)
		self.total_exons += 1
	
	def find_exon(self, pos):
		chunk = pos//self.chunk_size
		try:
			exon = self.chunk_nodes[chunk].find_exon(pos)
		except KeyError:
			exon = []
		return exon

	def print_node(self):
		for node in self.chunk_nodes:
			self.chunk_nodes[node].print_node()

class ESTChunkNode:
	# A layer of chunk nodes (for now, every million BP gets its own chunk node)
	def __init__(self, chr, number):
		self.chr = chr
		self.chunk_num = number
		self.exons = {}
	
	def add_node(self, start, end, name, pos_in_gene):
		ident = "%d#%s", (start, name)
		self.exons[ident] = ESTExonLeaf(self.chr, start, end, name, pos_in_gene)
	
	def find_exon(self, pos):
		names = []
		for ex in self.exons:
			if pos >= self.exons[ex].start:
				exon = self.exons[ex]
				if pos <= exon.end:
					 names.append([exon.name, exon.pos_in_gene+(pos-exon.start)])
		return names

	def print_node(self):
		for node in self.exons:
			self.exons[node].print_exon()

class ESTExonLeaf:
	# The leaf node that contains the exon.
	def __init__(self, chr, start, end, name, pos_in_gene):
		self.chr = chr
		self.start = start
		self.end = end
		self.name = name
		self.pos_in_gene = pos_in_gene
		self.id = "%d#%s" % (start, name)
	
	def print_exon(self):
		print(self.chr, self.start, self.end, self.name, self.pos_in_gene)
		
if __name__=="__main__":
	# This section only gets run if you call exonSearchTree.py directly from the command line.
	# So I can test the code without running it in the larger context if I want.
	est = ExonSearchTree(1000)
	est.add_exon(1, 100000000, 100002000, 1500, "NP_5")
	est.add_exon(1, 99999900, 100000500, 1499, "NP_7")
	est.add_exon('M', 999990, 1000500, 280, "NP_6")
	est.add_exon('Q', 3939, 8858, 0)
	est.print_tree()
	print(est.find_exon('M', 999995))
	print(est.find_exon('1', 100000001))
	print(est.find_exon('1', 99999901))
	print(est.find_exon('2', 949034))

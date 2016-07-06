class Domain():
	"""
	Represents a single domain
	"""
	def __init__(self):
		self.dom_id = ""
		self.dom_name = ""
		self.coords = ()
		self.pre_post_coords = ()

	def is_within(self, mut_pos):
		"""
		Return True if the mutation is within (+-10 aa) a given domain but NOT inside.
		"""
		# Check if BEFORE the domain
		before = False
		after = False
		if (mut_pos < self.coords[0] and mut_pos >= self.pre_post_coords[0]):
		 	before = True
		 	print "Mutation is located just BEFORE the domain ", self.dom_name
		
		# Check if AFTER the domain		
		elif (mut_pos > self.coords[1] and mut_pos <= self.pre_post_coords[1]):
			after = True
			print "Mutation is located just AFTER the domain ", self.dom_name

		if before or after:
			return True
		else:
			return False

	def is_inside(self, mut_pos):
		"""
		Return True if the mutation is inside a given domain.
		"""
		if mut_pos >= self.coords[0] and mut_pos <= self.coords[1]:
			print "Mutation is INSIDE the domain ", self.dom_name
			return True
		else:
			return False

class Mutation():
	"""
	Represents a single mutation
	"""
	def __init__(self):
		self.aa_code = ""
		self.aa_mut = ""
		self.pos = 0
		self.mut_type = ""
		self.diseases= {}
		self.pubmed_ids = []
	
class Protein():
	"""
	Czyta pliki z mutacja i domenami i buduje listy obiektow z domenami i mutacjami.
	"""
	def __init__(self, domains_file, mut_file):
		self.domains_file = domains_file
		self.mut_file = mut_file
		self.uniprot_id = ""
		self.seq_name = ""
		self.domains = []
		self.sequence = ""
		self.mutations = []
		self.seq_len = 0

	def read_mutations(self):
		"""
		Create mutations objects
		"""
		mut_file = open(self.mut_file, "r").read().splitlines()

		for line in mut_file:
			if line.startswith("KEY"):
				continue
			if line.startswith("UNIPROT_ACC"):
				continue
			if line.startswith("UNIPROT_ID"):
				continue
			if line.startswith("SEQ_NO"):
				mutation = Mutation()
				mutation.pos = float(line.split()[1])
				continue
			if line.startswith("AA_CODE"):
				mutation.aa_code = line.split()[1]
				continue
			if line.startswith("AA_MUT"):
				mutation.aa_mut = line.split()[1]
				continue
			if line.startswith("MUT_TYPE"):
				mutation.mut_type = line.split()[1]
				continue
			if line.startswith("DISEASE_ID"):
				id = line.split()[1]
				mutation.diseases[id] = ""
				continue
			if line.startswith("DISEASE_NAME"):
				mutation.diseases[id] = " ".join(line.split()[1:])
				continue

			self.mutations.append(mutation)

		def locate_mutations(self):
			pass

	def read_domains(self):
		"""
		Create domains objects
		"""
		dom_file = open(self.domains_file, "r").read().splitlines()

		new_domain = False
		for line in dom_file:
			if line.startswith("ID"):
				self.uniprot_id = line.split()[1]
				continue
			elif line.startswith("SEQ_NAME"):
				self.seq_name = " ".join(line.split()[1:])
				continue
			elif line.startswith("SEQUENCE"):
				self.sequence = line.split()[1]
				continue
			elif line.startswith("DOMAIN_ID"):
				new_domain = True
				dom_id = line.split()[1]
			elif line.startswith("DOMAIN_NAME"):
				dom_name = line.split()[1]
				continue
			elif line.startswith("DOMAIN_START"):
				domain = Domain()
				domain.dom_name = dom_name
				domain.dom_id = dom_id
				domain.coords = (float(line.split()[1]),)
				continue
			elif line.startswith("DOMAIN_PRESTART"):
				domain.pre_post_coords = (float(line.split()[1]),)
				continue
			elif line.startswith("DOMAIN_END"):
				domain.coords += (float(line.split()[1]),)
				continue
			elif line.startswith("DOMAIN_POSTEND"):
				domain.pre_post_coords += (float(line.split()[1]),)
				self.domains.append(domain)
			

		for el in self.domains:
			print el.pre_post_coords


if __name__ == '__main__':
	mutfile = "umutations.dat"
	domfile = "domains_Q9NQV7.dat"

	protein = Protein(domfile, mutfile)
	protein.read_mutations()
	protein.read_domains()


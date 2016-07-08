import csv, re

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
		Return True if the mutation/natural variant is within (+-10 aa) a given domain but NOT inside.
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

	def natural_variants_enrichment(self):
		"""
		"""
		pass

class Mutation():
	"""
	Represents a single mutation.
	"""
	def __init__(self):
		self.aa_code = ""
		self.aa_mut = ""
		self.pos = 0
		self.mut_type = ""
		self.diseases= {}
		self.pubmed_ids = []


class Natural_variant():
	"""
	Represents a single natural variant extracted from ExAC database.
	"""
	def __init__(self):
		self.protein_loc = 0
		self.aa_code = "" # source aa
		self.aa_mut = "" # mutated
		self.annotation = "" # missennse, 
		self.flags = "" # is it LoF
		self.allele_count = 0
		self.allele_frequency = 0.0
		self.homo_number = 0

	def is_disease_assoc(self, mut_pos, mut_aa):
		"""
		Is the variant in the uniprot disease associated list?
		To make sure, check the aa code of the mutated aa.
		"""

		assert(mut_aa==aa_code)

		if mut_pos == self.protein_loc:
			return True
		else:
			return False

		

class Alignment():
	"""
	Represents a single alignment of target (our) sequence with template.
	"""
	def __init__(self):
		self.pdf_code = ""
		self.chain = ""
		self.identity = 0.0
		self.evalue = 0.0
		self.pdb_start = 0
		self.pdb_end = 0
		self.ali_start = 0 #first aligned residue in source seq
		self.ali_end = 0
		self.ref_seq = ""
		self.pdb_seq = ""
	
class Protein():
	"""
	Czyta pliki z mutacja i domenami i buduje listy obiektow z domenami i mutacjami.
	"""
	def __init__(self, domains_file, mut_file, variants_file, ali_file):
		self.domains_file = domains_file
		self.mut_file = mut_file
		self.variants_file = variants_file
		self.ali_file = ali_file
		self.uniprot_id = ""
		self.seq_name = ""
		self.domains = []
		self.sequence = ""
		self.mutations = []
		self.variants = []
		self.alignments = []
		self.seq_len = 0

	def read_mutations(self):
		"""
		Create mutations objects
		"""
		mut_file = open(self.mut_file, "r").read().splitlines()

		for line in mut_file:
			if line.startswith("KEY"):
				continue
			elif line.startswith("UNIPROT_ACC"):
				continue
			elif line.startswith("UNIPROT_ID"):
				continue
			elif line.startswith("SEQ_NO"):
				mutation = Mutation()
				mutation.pos = float(line.split()[1])
				continue
			elif line.startswith("AA_CODE"):
				mutation.aa_code = line.split()[1]
				continue
			elif line.startswith("AA_MUT"):
				mutation.aa_mut = line.split()[1]
				continue
			elif line.startswith("MUT_TYPE"):
				mutation.mut_type = line.split()[1]
				continue
			elif line.startswith("DISEASE_ID"):
				id = line.split()[1]
				mutation.diseases[id] = ""
				continue
			elif line.startswith("DISEASE_NAME"):
				mutation.diseases[id] = " ".join(line.split()[1:])
				continue
			elif line.startswith("NREFS"):
				self.mutations.append(mutation)


	def locate_mutations(self):
		"""
		Locates the mutations and answers the question if they lay in
		the close proximity to any functional domain or inside it.
		"""

		in_domain = 0
		in_prox = 0

		for mutation in self.mutations:
			for domain in self.domains:
				if domain.is_within(mutation.pos):
					in_prox+=1
					break
				elif domain.is_inside(mutation.pos):
					in_domain+=1
					break
				else: continue
			else: break


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
				self.seq_len = float(len(self.sequence))
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


	def read_natural_variants(self):
		"""
		Creates list of nv objects.
		For now - all of them, even those inside introns and 3'UTR
		"""
		csvfile = open(self.variants_file, 'rb')
		variants = csv.reader(csvfile, delimiter=',', quotechar='"')
		next(variants) # omit header

		for line in variants:
			variant = Natural_variant()
			if line[6]:
				aa_exch = line[6].split(".")[1]
				parts = re.findall("([A-Za-z]{3})([0-9]+)([A-Za-z]{3}|\?)",aa_exch)
				if parts:
					variant.protein_loc = parts[0][1]
					variant.aa_code = parts[0][0]
					variant.aa_mut = parts[0][2]		
				else:pass

			variant.annotation = line[9] # missennse etc
			variant.flags = line[10] # is it LoF
			variant.allele_count = float(line[11])
			variant.allele_frequency = float(line[14])
			variant.homo_number = float(line[13])

			self.variants.append(variant)		

	def read_alignments(self):
		"""
		Creates list of alignment objects.
		"""
		alifile = open(self.ali_file, "r").read().splitlines()

		for line in alifile:
			if line.startswith("ALIGN_PDB_CODE"):
				alignment = Alignment()
				alignment.pdf_code = line.split()[1]
				continue
			elif line.startswith("ALIGN_CHAIN"):
				alignment.chain = line.split()[1]
				continue
			elif line.startswith("ALIGN_E_VALUE"):
				alignment.evalue  = float(line.split()[1])
				continue
			elif line.startswith("ALIGN_OVERLAP_FROM"):
				alignment.pdb_start = float(line.split()[1])
				continue
			elif line.startswith("ALIGN_OVERLAP_TO"):
				alignment.pdb_end = float(line.split()[1])
				continue
			elif line.startswith("UNIPROT_SEQ"):
				alignment.ref_seq = line.split()[1]
				continue
			elif line.startswith("PDB_SEQ"):
				alignment.pdb_seq = line.split()[1]
				continue
			elif line.startswith("ALIGN_START_RES"):
				alignment.ali_start = float(line.split()[1])
				continue
			elif line.startswith("ALIGN_END_RES"):
				alignment.ali_end = float(line.split()[1])
				self.alignments.append(alignment)


	def get_natural_variant_stats(self):#, region_start=0, region_end=self.seq_len):
		"""
		Calculates basic stats for the variants list.
		"""
		n_variants = float(len(self.variants))
		n_missense = 0
		n_synonymous = 0
		n_stops = 0
		n_introns = 0
		n_utr = 0
		n_frameshift = 0
		n_deletions = 0

		number_of_variants_per_seq_len = n_variants/self.seq_len

		for variant in self.variants:
			if variant.annotation == "5' UTR":
				n_utr+=1
			elif variant.annotation == "synonymous":
				n_synonymous+=1
			elif variant.annotation == "missense":
				n_missense+=1
			elif variant.annotation == "stop gained":
				n_stops+=1
			elif variant.annotation == "frameshift":
				n_frameshift+=1
			elif variant.annotation == "inframe deletion":
				n_deletions+=1
			elif variant.annotation == "intron":
				n_introns+=1

		stats_list = {"n_variants":n_variants, "n_missense":n_missense, "n_synonymous":n_synonymous, "n_stop":n_stops,
						"n_frameshift":n_frameshift, "n_introns":n_introns, "n_utr": n_utr, "n_deletions": n_deletions, 
						"num_per_len": number_of_variants_per_seq_len}

		return stats_list

if __name__ == '__main__':
	mutfile = "umutations.dat"
	domfile = "domains_Q9NQV7.dat"
	variantfile = "exac_ENSG00000164256.csv"
	alifile = "alignments.dat"

	protein = Protein(domfile, mutfile, variantfile, alifile)
	protein.read_mutations()
	protein.read_domains()
	protein.read_natural_variants()
	protein.read_alignments()
	#protein.locate_mutations()
	stats = protein.get_natural_variant_stats()
	print stats


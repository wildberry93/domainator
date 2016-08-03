import csv, re, glob, os
import collections

class Domain():
	"""
	Represents a single domain
	"""
	def __init__(self):
		self.dom_id = ""
		self.dom_name = ""
		self.coords = ()
		self.pre_post_coord = ()

	def is_within(self, mut_pos):
		"""
		Return True if the mutation/natural variant is 
		within (+-10 aa) a given domain but NOT inside.
		"""
		# Check if BEFORE the domain
		before = False
		after = False
		if (mut_pos < self.coords[0] and mut_pos >= self.pre_post_coords[0]):
		 	before = True
		 	print "Mutation is located just BEFORE the domain ", self.dom_id
		
		# Check if AFTER the domain		
		elif (mut_pos > self.coords[1] and mut_pos <= self.pre_post_coords[1]):
			after = True
			print "Mutation is located just AFTER the domain ", self.dom_id

		if before or after:
			return True
		else:
			return False

	def is_inside(self, mut_pos):
		"""
		Return True if the mutation is inside a given domain.
		"""
		if mut_pos >= self.coords[0] and mut_pos <= self.coords[1]:
			print "Mutation is INSIDE the domain ", self.dom_id
			return True
		else:
			return False
		

class Model():
	"""
	Necessary for pdb start and pdb end.
	"""
	def __init__(self):
		self.mapped_pdb = ""
		self.pdb_chain = ""
		self.start = 0
		self.end = 0
		self.start_aa = ""
		self.end_aa = ""

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

	def is_natural_variant(self, variants):
		"""
		Check if the mutation is in natural variants set.
		"""

		positions = [variant.protein_loc for variant in variants]

		if self.pos in positions:
			return True
		else: 
			return False

class Natural_variant():
	"""
	Represents a single natural variant extracted from ExAC database.
	"""
	def __init__(self):
		self.protein_loc = 0
		self.aa_code = "" # source aa
		self.aa_mut = "" # mutated
		self.annotation = "" # missennse
		self.flags = "" # is it LoF
		self.allele_count = 0
		self.allele_frequency = 0.0
		self.homo_number = 0

	def is_disease_assoc(self, mutations):
		"""
		Is the variant in the uniprot disease associated list?
		To make sure, check the aa code of the mutated aa.
		"""

		positions = [mutation.pos for mutation in mutations]

		if self.protein_loc in positions:
			return True
		else:
			return False

	def is_interesting(self):
		"""
		Check if variant is interesting
		"""
		af_threshold = 10e-6
		if self.annotation not in ["neutral","intron","synonymous","5' UTR","start lost"] and\
			self.allele_frequency <= af_threshold:
			return True
		else:
			return False

class Alignment():
	"""
	Represents a single alignment of target (our) sequence with template.
	"""
	def __init__(self):
		self.pdb_code = ""
		self.chain = ""
		self.identity = 0.0
		self.evalue = 0.0
		self.pdb_start = 0
		self.pdb_end = 0
		self.ali_start = 0 #first aligned residue in source seq
		self.ali_end = 0
		self.ref_seq = ""
		self.pdb_seq = ""
		self.ali_len = 0
		self.ali_dict = {}#collections.OrderedDict()
		self.model = 0

	def make_ali_dict(self):
		ali_dict = {}

		start_ali = self.ali_start #1 based
		print "ali_start", start_ali
		print "ali_end", self.ali_end
		start_pdb = self.pdb_start
		print "pdb_start", start_pdb
		print "pdb_end",self.pdb_end

		print "r",self.ref_seq
		print "p", self.pdb_seq

		for i in range(0, len(self.ref_seq)):
			print str(start_ali) + "\t" + str(start_pdb)
			if self.ref_seq[i] != '-' and self.pdb_seq[i] == '-':
				ali_dict[start_ali] = '-'
				start_ali+=1
			elif self.ref_seq[i] == '-' and self.pdb_seq[i] != '-':
				start_pdb+=1
			elif self.ref_seq[i] != '-' and self.pdb_seq[i] != '-':
				ali_dict[start_ali] = start_pdb
				start_ali+=1
				
			else: 
				print "cos jest nie tak - i -"

			
			print self.ref_seq[i]+ "\t" + self.pdb_seq[i]


		return ali_dict
	
class Alignment2():
	"""
	Represents a single alignment from Romans' alignments of
	target (our) sequence with template.
	"""
	def __init__(self):
		self.pdb_code = ""
		self.chain = ""
		self.identity = 0.0
		self.evalue = 0.0	
		self.ali_start = 0 
		self.ali_end = 0
	
class Protein():
	"""
	Does all the work.
	"""
	def __init__(self, domains_file, mut_file, variants_file, ali_file, model_file, ali_file2):
		self.domains_file = domains_file
		self.mut_file = mut_file
		self.variants_file = variants_file
		self.ali_file = ali_file
		self.ali_file2 = ali_file2
		self.model_file = model_file
		self.uniprot_id = ""
		self.seq_name = ""
		self.domains = []
		self.sequence = ""
		self.mutations = []
		self.variants = []
		self.alignments = []
		self.alignments2 = []
		self.models = []
		self.seq_len = 0

	def read_mutations(self):
		"""
		Create mutations objects.
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
		outside = 0

		domains_file = open("mutations2domains_"+self.uniprot_id+".txt", "w")

		for mutation in self.mutations:
			for domain in self.domains:
				if domain.is_within(mutation.pos):
					in_prox+=1

					domains_file.write(domain.dom_id+"\t"+str(mutation.pos)+"\t"+
						mutation.aa_code+"\t"+mutation.aa_mut+"\t"+
						",".join(mutation.diseases.values())+"\n")
					break
				elif domain.is_inside(mutation.pos):
					in_domain+=1

					domains_file.write(domain.dom_id+"\t"+str(mutation.pos)+"\t"+
						mutation.aa_code+"\t"+mutation.aa_mut+"\t"+
						",".join(mutation.diseases.values())+"\n")
					break
				else:
					outside+=1 
					continue
			else: break

		domains_file.close()

		return (in_domain, in_prox, outside)


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
					variant.protein_loc = float(parts[0][1])
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
			if line.startswith("ID"):
				alignment = Alignment()
				self.alignments.append(alignment)
				id = line.split(":")[1]
				if len(id) == 5:
					alignment.chain = id[-1]
				alignment.pdb_code = id[:4]
			else:
				ref = line.split("\t")[0]
				pdb = line.split("\t")[1]
				ref_pos = ref.split()[0]
				ref_aa = ref.split()[1]
				pdb_pos = pdb.split()[0]
				pdb_aa = pdb.split()[1]
				alignment.ali_dict[int(ref_pos)] = (ref_aa, int(pdb_pos), pdb_aa)

	def read_alignments2(self):
		"""
		Creates list of alignment objects.
		"""
		alifile = open(self.ali_file2, "r").read().splitlines()

		for line in alifile:
			if line.startswith("ALIGN_PDB_CODE"):
				alignment = Alignment2()
				self.alignments2.append(alignment)
				alignment.pdb_code = line.split()[1]
				continue
			elif line.startswith("ALIGN_CHAIN"):
				alignment.chain = line.split()[1]
				continue
			elif line.startswith("ALIGN_E_VALUE"):
				alignment.evalue  = float(line.split()[1])
				continue
			
	def read_models(self):
		"""
		Reads models.dat for every protein.
		"""
		modfile = open(self.model_file, "r").read().splitlines()

		for line in modfile:
			if line.startswith("MAPPED_PDB_CODE"):
				model = Model()
				model.mapped_pdb = line.split()[1]
				continue
			elif line.startswith("MAPPED_PDB_CHAIN"):
				model.pdb_chain = line.split()[1]
				continue
			elif line.startswith("MAPPED_START_RES"):
				model.start_aa = line.split()[1]
				model.start = float(line.split()[2])
				continue
			elif line.startswith("MAPPED_END_RES"):
				model.end_aa = line.split()[1]
				model.end = float(line.split()[2])
				self.models.append(model)


	def get_variants_from_region(self, region_start=0, region_end=0):
		"""
		Get the subset of variants that are inside or in 
		the proximity to the domain.
		"""
		if not region_end:  # if the end is not given, then take the whole sequence
			region_end = self.seq_len

		reg_variants = []
		
		for variant in self.variants:
			if float(variant.protein_loc) >= region_start and\
				 float(variant.protein_loc) <= region_end:
				if variant.annotation not in ["neutral","intron","synonymous","5' UTR","start lost"]: 
					reg_variants.append(variant)

		return reg_variants

	def get_mutations_from_region(self, region_start=0, region_end=0):
		"""
		Get the subset of mutations that are inside or in the proximity to the domain.
		"""
		if not region_end:  # if the end is not given, then take the whole sequence
			region_end = self.seq_len

		reg_mutations = []
		
		for mutation in self.mutations:
			if mutation.pos >= region_start and mutation.pos <= region_end:
				reg_mutations.append(mutation)

		return reg_mutations

	def get_natural_variant_stats(self, region_start=0, region_end=0):
		"""
		Calculates basic stats for the variants list.
		"""

		if not region_end:  # if the end is not given, then take the whole sequence
			region_end = self.seq_len

		n_all_variants = float(len(self.variants)) #number of variants across entire protein
		n_missense = 0
		n_synonymous = 0
		n_stops = 0
		n_introns = 0
		n_utr = 0
		n_frameshift = 0
		n_deletions = 0
		n_variants = 0

		region_length = region_end-region_start

		for variant in self.variants:
			if variant.protein_loc>=region_start and variant.protein_loc<=region_end:
				n_variants+=1
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


		number_of_variants_per_seq_len = n_variants/region_length

		stats_list = {"n_all_variants":n_all_variants, "n_variants_in_region":n_variants, 
					  "n_missense":n_missense, "n_synonymous":n_synonymous, "n_stop":n_stops,
					  "n_frameshift":n_frameshift, "n_introns":n_introns, "n_utr": n_utr, 
					  "n_deletions": n_deletions,  "num_per_len": number_of_variants_per_seq_len}


		return stats_list

	def get_af_distribution():
		"""
		Return the allele frequency threshold all alleles 
		below this level are assumed important and later on 
		draw a plot (later - I hate plotting in Python)
		"""
		frequencies = [variant.allele_frequency for variant in self.variants]
		mean_freq = mean(frequencies)

		return mean_freq

	def natural_variants_enrichment(self):
		"""
		Calculate the natural variants enrichment inside the annotated domains.
		Basically, just run get_natural_variant_stats() on domain coordinates.
		Filter by those with low Allele Frequency that indicates they are rare.
		Proposes the new disease causing mutations.

		Which allele frequency is considered low?
		http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3207674/
		"""

		af_threshold = 10e-6 #self.get_af_distribution()
	
		#enrichment_file = open("enrichment_"+self.uniprot_id+".txt", "w")
		interesting_file = open("interesting_variants_"+self.uniprot_id+".txt", "w")

		print self.uniprot_id

		for domain in self.domains:
			if domain.dom_id == "UNASS": continue
			interesting_variants = []
			dom_variants = self.get_variants_from_region(domain.pre_post_coords[0], 
															domain.pre_post_coords[1])
			dom_stats = self.get_natural_variant_stats(domain.pre_post_coords[0], 
														domain.pre_post_coords[1])	
			#enrichment_file.write(domain.dom_id+"\t"+",".join([":".join((x, str(dom_stats[x]))) for x in dom_stats])+"\n")		

			for variant in dom_variants:
				if variant.allele_frequency <= af_threshold:
					if not variant.is_disease_assoc(self.mutations):
						if variant.annotation not in ["neutral","intron","synonymous","5' UTR","start lost"]: 
							interesting_variants.append(variant)
							interesting_file.write(domain.dom_id+"\t"+str(variant.protein_loc)+
													"\t"+variant.aa_code+"\t"+variant.aa_mut+"\t"+
													variant.annotation+"\t"+variant.flags+"\t"+
													str(variant.allele_frequency)+"\n")
					else:
						print "WARNING! Low frequency variant but already in UniProt mutations!"
		#enrichment_file.close()
		interesting_file.close()

	

def do_single(dir):
	"""
	Perform whole procedure over a single protein entry. 
	The input directory should contain mutations and domains files.
	The alignment file is not required but recommended since it's 
	necessary for full analysis.
	"""
	infiles = glob.glob(dir+"\*")
	#dis_file = open("diseases_mapping.txt", "a")
	strs_file = open("proteins2structures.txt", "a")
	muts_pdbs = open("mutations_in_pdbs.txt", "a")
	var_pdbs = open("variants_in_pdbs.txt", "a")

	#if alifile does not exist we dont need models and the other way round
	if os.path.join(dir,"alignments.dat") not in infiles:
		print "WATCH OUT! The alignment in unavaliable!"
		alifile = None 
		alifile2 = None
		modelfile = None
	else:
		alifile = os.path.join(dir,"all_alignments.dat")
		alifile2 = os.path.join(dir, "alignments.dat")
		modelfile = os.path.join(dir, "models.dat")


	prot_name = os.path.basename(dir)
	domfile = os.path.join(dir, "domains_"+prot_name+".dat")
	mutfile = os.path.join(dir, "umutations.dat")
	varfile = os.path.join(dir, prot_name+".csv")
	

	protein = Protein(domfile, mutfile, varfile, alifile, modelfile, alifile2)
	protein.read_mutations()
	protein.read_domains()
	
	#dis_file.write(protein.uniprot_id+"\t"+muts+"\n")

	try:
		protein.read_natural_variants()
	except:
		varfile = None
		pass

	if alifile and alifile2: 
		protein.read_alignments() #if not none
		protein.read_alignments2()
		#protein.read_models()

	not_dbd = ["BTB","SCAN", "PAS", "KRAB","MH2","PF07710","AXH", 
			   "PF01166","PF15784","SANT", "PF00104","PINIT","LIM", 
			   "ELM2","MYST-type","AWS","KID", "ANK","PF04814","OTU",
			   "PF11717", "JmjN","PF12284","SH2","PF02865","PF12188",
			   "PF02518","PF16534","PF11620","PF03531","Bromo","KIX",
			   "Tudor","Pre-SET","PF08169","SWIRM","BAH","PF16159",
			   "PF11569","PF01167","CVC","ELM2","PF14048","Helicase",
			   "WHOLE"] #"PF01167"

	#print prot_name
	#filter the alignments
	red = [] # redundant
	new_alignments = [] # filtered alignments
	for ali2 in protein.alignments2:
		for ali in protein.alignments:
			if ali2.pdb_code == ali.pdb_code and ali2.chain == ali.chain and\
				(ali2.pdb_code, ali2.chain) not in red:
				new_alignments.append(ali)
				red.append((ali2.pdb_code, ali2.chain))


	#get natural variants and mutations from alignments
	muts_mapped = {}
	vars_mapped = {}
	for ali in new_alignments:
		ali.ali_start = min(ali.ali_dict.keys())
		ali.ali_end = max(ali.ali_dict.keys())
		muts = protein.get_mutations_from_region(ali.ali_start, ali.ali_end)
		vars = protein.get_variants_from_region(ali.ali_start, ali.ali_end)

		for domain in protein.domains:
			ol = overlap(domain.coords, (ali.ali_start, ali.ali_end))
			if ol >= 30 and domain.dom_id and domain.dom_id not in not_dbd:

				if ali.pdb_code+"_"+ali.chain not in muts_mapped:
					muts_mapped[ali.pdb_code+"_"+ali.chain] = []

				if ali.pdb_code not in vars_mapped:
					vars_mapped[ali.pdb_code+"_"+ali.chain] = []

				muts_mapped[ali.pdb_code+"_"+ali.chain]+=[ali.ali_dict[int(mut.pos)][1]\
							for mut in muts if int(mut.pos) in ali.ali_dict]
				vars_mapped[ali.pdb_code+"_"+ali2.chain]+=[ali.ali_dict[int(var.protein_loc)][1]\
							for var in vars if var.is_interesting()]
				break
	'''if prot_name == "P78411":
		print "vars_mapped", muts_mapped
		print "muts", [x.pos for x in muts]
		print "alis", [a.ali_dict for a in new_alignments]'''

	print muts_mapped

	if muts_mapped:
		muts_pdbs.write(prot_name+"\t"+ali.pdb_code+"\t"+
						ali.chain+"\t"+",".join(muts_mapped)+"\n")
	if vars_mapped:
		var_pdbs.write(prot_name+"\t"+ali.pdb_code+"\t"+
						ali.chain+"\t"+",".join(vars_mapped)+"\n")


	#if varfile: protein.natural_variants_enrichment()

	#protein.locate_mutations()


def do_all(up_dir, selected_set=None):
	"""
	Perform for all given proteins.
	"""

	all_proteins = glob.glob(up_dir+"\*")
	
	if selected_set:
		prot_seq = [x for x in all_proteins if os.path.basename(x) in selected_set]
		all_proteins = prot_seq


	for prot in all_proteins:
		do_single(prot)


def overlap( (s1, e1), (s2, e2) ):
	""" Return length of the overlapping fragment of two segments. """
	
	(s1, e1), (s2, e2) = sorted( [(s1, e1), (s2, e2)] )
	if s2 <= e1:
		return e1 - s2 + 1
	return 0

def get_start_end(string):
	stru_rev = string[::-1]

	parts = re.findall("[0-9]+",string)

	pre_cnt = 0
	for char in string:
		if char.isdigit():
			break
		else: 
			pre_cnt+=1

	stringi = str(parts[0])

	a = 0
	if len(stringi) == 3:
		a = 2
	if len(stringi) == 2:
		a = 1
	if len(stringi) == 1:
		a = 1

	pre_add = float(parts[0])-pre_cnt-a

	post_cnt = 0

	for char in stru_rev:
		if char.isdigit():
			break
		else: 
			post_cnt+=1

	post_add = post_cnt+float(parts[-1])

	return pre_add, post_add


if __name__ == '__main__':
	#do_single("C:\Users\jagoda\Desktop\domainator\processed_entries\O00327")

	#sel_set = open("C:\Users\jagoda\Desktop\domainator\/not_mutated_list.txt", "r").read().splitlines()
	do_all("C:\Users\jagoda\Desktop\domainator\processed_entries")#, sel_set)
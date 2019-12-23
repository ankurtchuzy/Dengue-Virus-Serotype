
# <--------------------DEFENITIONS-------------------->
#Reading fasta sequence from file and storing sequence without header in a string
def read_fasta(file1):
	lines = []
	for a in open (file1):
		a = a.strip()
		if a != "":
			lines.append(a)
	
	header = lines.pop(0)
	header = header.replace('>','').split(' ')
	seq = ''.join(lines)
	return seq,header[0]


#Reading a CSV file and creating dictionary of SNPs position and variant
def read_poly (file2):
	records = []
	for n in open (file2):
		n = n.strip().replace('-','0')
		records.append(n)
	records.pop(0)

	pairs = {}
	for b in records:
		if b != "":
			fields = b.split(',')
			if fields[1] != '0':
				variant = 'A'
				var_depth = int(fields[7])
			elif fields[2] != '0':
				variant = 'C'
				var_depth = int(fields[8])
			elif fields[3] != '0':
				variant = 'G'
				var_depth = int(fields[9])
			elif fields[4] != '0':
				variant = 'T'
				var_depth = int(fields[10])
		
			if fields[5] == 'a':
				target = 'A'
				tar_depth = int(fields[7])
			elif fields[5] == 'c':
				target = 'C'
				tar_depth = int(fields[8])
			elif fields[5] == 'g':
				target = 'G'
				tar_depth = int(fields[9])
			elif fields[5] == 't':
				target = 'T'
				tar_depth = int(fields[10])
		pairs[fields[0]] = [variant,var_depth,target,tar_depth]
	return pairs


#translating genomic region using standard codon table
def comp_translation (wild,mut):
	codons = {'TTT':'F','TCT':'S','TAT':'Y','TGT':'C','TTC':'F','TCC':'S','TAC':'Y','TGC':'C','TTA':'L','TCA':'S','TAA':'*','TGA':'*','TTG':'L','TCG':'S','TAG':'*','TGG':'W',
	'CTT':'L','CCT':'P','CAT':'H','CGT':'R','CTC':'L','CCC':'P','CAC':'H','CGC':'R','CTA':'L','CCA':'P','CAA':'Q','CGA':'R','CTG':'L','CCG':'P','CAG':'Q','CGG':'R',
	'ATT':'I','ACT':'T','AAT':'N','AGT':'S','ATC':'I','ACC':'T','AAC':'N','AGC':'S','ATA':'I','ACA':'T','AAA':'K','AGA':'R','ATG':'M','ACG':'T','AAG':'K','AGG':'R',
	'GTT':'V','GCT':'A','GAT':'D','GGT':'G','GTC':'V','GCC':'A','GAC':'D','GGC':'G','GTA':'V','GCA':'A','GAA':'E','GGA':'G','GTG':'V','GCG':'A','GAG':'E','GGG':'G'}
	wild_peptide=[]
	mut_peptide=[]
	#translating wild and mutant region
	count = 1
	for i in range(0,len(mut),3):
		if (''.join(wild[i:i+3])) != (''.join(mut[i:i+3])):
			if codons[''.join(wild[i:i+3])] == codons[''.join(mut[i:i+3])]:
				replacement = codons[''.join(wild[i:i+3])]+" == "+codons[''.join(mut[i:i+3])]+"\t"+str(count)+"\t"+str(len(wild)//3)+"\tsynm"
			else:
				replacement = codons[''.join(wild[i:i+3])]+" -> "+codons[''.join(mut[i:i+3])]+"\t"+str(count)+"\t"+str(len(wild)//3)+"\tnon-synm"
		
		wild_peptide.append(codons[''.join(wild[i:i+3])])
		mut_peptide.append(codons[''.join(mut[i:i+3])])
		count = count + 1
	
	#comparing peptides and determining position of difference
	#count = 1
	#for w,m in zip(wild_peptide,mut_peptide):
		#if w != m:
			#replacement = w+" -> "+m+" at "+str(count)+" ("+str(len(wild_peptide))+"aa long)"
			#break
		#else:
			#replacement = "Synonymous"
		#count = count + 1
	return replacement
	

#getting regions for translation with SNP information and reporting changes in peptide
def trans_region (start,end,pos,var):
	mut_genome = list(genome)
	mut_genome[pos] = var
	wild_region = list(genome [start:end])      #region from the wild type genome
	mut_region = mut_genome [start:end]         #region from the mutant genome
	change = comp_translation(wild_region,mut_region)
	return change 
	
	
# <----------------------PROGRAM---------------------->


#getting genome from fasta file
fasta_file_name = "ncbi_dengu4.fasta"
genome,accn = read_fasta (fasta_file_name)

#getting dictionary of SNPs from CSV file
poly_file_name = "poly_dengu4.csv"
poly = read_poly (poly_file_name)

#mature peptide coding regions
#start,end,gene,locus_tag,gene_synonym,product,protein_id

peptide_regions = [[102,440,"POLY","DV4_gp1","polyprotein gene","anchored capsid protein ancC","NP_740314.1"],
[102,398,"POLY","DV4_gp1","polyprotein gene","capsid protein C","NP_740313.1"],
[441,938,"POLY","DV4_gp1","polyprotein gene","membrane glycoprotein precursor prM","NP_740315.1"],
[441,713,"POLY","DV4_gp1","polyprotein gene","protein pr","YP_009164957.1"],
[714,938,"POLY","DV4_gp1","polyprotein gene","membrane glycoprotein M","NP_740316.1"],
[939,2423,"POLY","DV4_gp1","polyprotein gene","envelope protein E","NP_740317.1"],
[2424,3479,"POLY","DV4_gp1","polyprotein gene","nonstructural protein NS1","NP_740318.1"],
[3480,4133,"POLY","DV4_gp1","polyprotein gene","nonstructural protein NS2A","NP_740319.1"],
[4134,4523,"POLY","DV4_gp1","polyprotein gene","nonstructural protein NS2B","NP_740320.1"],
[4524,6377,"POLY","DV4_gp1","polyprotein gene","nonstructural protein NS3","NP_740321.1"],
[6378,6758,"POLY","DV4_gp1","polyprotein gene","nonstructural protein NS4A","NP_740322.1"],
[6759,6827,"POLY","DV4_gp1","polyprotein gene","protein 2K","NP_740323.1"],
[6828,7562,"POLY","DV4_gp1","polyprotein gene","nonstructural protein NS4B","NP_740324.1"],
[7563,10262,"POLY","DV4_gp1","polyprotein gene","RNA-dependent RNA polymerase NS5","NP_740325.1"]]

#iterating over peptide region and cheking SNPs in those region

print ("genome accession","genomic position","Ref nucleotide","Alt nucleotide","Ref depth","Alt depth","Protein accession","Product","amino change","position","peptide length","type",sep = "\t")
for reg in peptide_regions:
	flag = 0 #flag if no SNP found in region
	for k,v in poly.items():
		if int(k)>=(reg[0]-1) and int(k)<=(reg[1]-1):
			pep_match_result = trans_region(reg[0]-1,reg[1],int(k)-1,v[0])
			print (accn,k,genome[int(k)-1],v[0],v[3],v[1],reg[6],reg[5],pep_match_result,sep = "\t")
			flag = 1
	#if flag == 0:
		#print (accn,'-','-','-','-','-',reg[6],reg[5],'-','-','-',"No SNPs in this region",sep ="\t")


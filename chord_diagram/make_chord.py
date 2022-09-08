import csv

collapsing_file = 'data/Collapsing_MASTER_Olink_n47k_comb_0.0001.csv'
pvalue_threshold = 1e-8
collapsing_models = ['ptv']

output_file = 'circos/collapsing_links.txt'


ensembl_gtf = 'data/Homo_sapiens.GRCh38.107.gtf'
hgnc_file = 'data/HGNC_names.txt'
olink_file = 'data/olink_protein_map_1.5k_v1-Manifest.txt'

color = {'pos': 'positive', 'neg': 'negative', 'both': 'both'}
color_opacity = 0.4
gene_width = 1000

ensg2pos = {}
for line in open(ensembl_gtf):
	if line[:1] == '#':
		continue
	row = line.split()
	if row[2] == 'gene':
		start, end = int(row[3]), int(row[4])
		chrom = row[0]
		ensg = line.split('gene_id "')[1].split('"')[0]
		ensg2pos[ensg] = (chrom, start, end)

hgnc2ensg = {}
for x in csv.DictReader(open(hgnc_file), delimiter='\t'):
	hgnc2ensg[x['Approved symbol']] = x['Ensembl gene ID']
	for y in x['Previous symbols'].split(', '):
		if y not in hgnc2ensg:
			hgnc2ensg[y] = x['Ensembl gene ID']

prot2genes = {}
for row in csv.DictReader(open(olink_file), delimiter='\t'):
	acc = row['UKBPPP_ProteinID'].split(':')[0]
	prot2genes[acc] = prot2genes.get(acc, []) + [row['HGNC.symbol']]
	hgnc2ensg[row['HGNC.symbol']] = row['ensembl_id']

complexes = set()
interactions = {}
for row in csv.DictReader(open(collapsing_file), delimiter=';'):
	gene = row['Gene']
	
	proteins = list(set(prot2genes[row['Protein']]))
	if len(proteins) > 1: # skip complexes
		complexes.add(row['Protein'])
		continue

	prot = proteins[0]
	if gene != prot:
		p = float(row['p'].replace(',','.'))
		model = row['model']
		beta = float(row['beta'].replace(',','.'))
		inter = 'neg' if beta < 0 else 'pos'

		#if model in collapsing_models:
		if p <= pvalue_threshold:
			if (gene, prot) in interactions:
				if interactions[(gene, prot)] != inter:
					interactions[(gene, prot)] = 'both'
			else:
				interactions[(gene, prot)] = inter

print(len(complexes), 'complexes skipped', complexes)
print(len(interactions), 'total interactions')
print()


with open(output_file, 'w') as f:
	for gene, prot in interactions:
		gene_pos, prot_pos = ensg2pos[hgnc2ensg[gene]], ensg2pos[hgnc2ensg[prot]]
		
		source_chrom = 'hs{}'.format(gene_pos[0])
		source_start = gene_pos[1]
		source_end = gene_pos[2]#source_start + gene_width

		target_chrom = 'hs{}'.format(prot_pos[0])
		target_start = prot_pos[1]
		target_end = prot_pos[2]#target_start + gene_width

		print(source_chrom, source_start, source_end, target_chrom, target_start, target_end, 
			'color='+color[interactions[(gene, prot)]], 
			sep='\t', file=f)
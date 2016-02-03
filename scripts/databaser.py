#!/usr/local/bin/python

import os, sys, re, datetime, subprocess, urllib, ftplib, gzip, tarfile

print '\n'.join((
	'',
	'Databaser: generate taxonomically balanced, non-redundant proteome database.',
	'Author: Qiyun Zhu (qzhu@jcvi.org)',
	'Affiliation: J. Craig Venter Institute',
	'Last updated: Jan 18, 2016',
	''
	))

if len(sys.argv) == 1:
	print '\n'.join((
	'Usage:',
	'  python scripts/databaser.py [-format=none|blast] [-range=all|microbe|etc] [-represent=no|auto|filename] [-subsample=0|1|2...] [-out=dbname]',
	'Note:',
	'  By default, the script will download NCBI RefSeq proteomic data of microbial organisms, keep one proteome per species, plus all NCBI-defined representative proteomes.',
	'  Use "-help" to display details of program parameters.'
	))
	raw_input('Press Enter to continue or Ctrl+C to exit...')
elif sys.argv[1] == '-h' or sys.argv[1].endswith('-help'):
	print '\n'.join((
	'Details of parameters:',
	'  -format=<none, blast (default)>',
	'    Format of the output database. If "none", only a master fasta file will be generated, based on which one may manually create customized databases.',
	'  -range=<taxonomic categories>',
	'    The NCBI RefSeq genome database has the following categories: archaea, bacteria, fungi, invertebrate, plant, protozoa, vertebrate_mammalian, vertebrate_other, and viral.',
	'    Type one or more desired categories, separated by comma (,).',
	'    Type "all" for all categories, or "microbe" for archaea, bacteria, fungi, and protozoa (default).',
	'  -represent=<no, auto (default) or filename>',
	'    By default, the script will download a list of NCBI-defined representative genomes, and include all of them in the resulting database.',
	'    One may also provide a user-defined genome list (format: one taxid per line), or switch off this function',
	'  -subsample=<0, 1 (default), 2, 3...>',
	'    By default, the script will keep up to one proteome per species.',
	'    Change the number for more proteomes, or "0" for all (no subsampling).',
	'  -out=<database name>',
	'    Please use a simple string without spaces or symbols. The default name is "stdb".',
	'',
	'Examples:',
	'  Default behavior:',
	'    python script/databaser.py -format=blast -range=microbe -subsample=1 -represent=auto',
	'  Download microbial eukaryotic proteomes:',
	'    python script/databaser.py -format=none -range=protozoa,fungi -subsample=0 -represent=no -out="MyEukarie"',
	''
	))
	sys.exit()

today = datetime.date.today()
fout = open('date.txt', 'w')
fout.write(str(today))
fout.close()

(format, range, repre, subsp, outname) = ('blast', 'microbe', 'auto', 1, 'stdb')
for argv in sys.argv[1:]:
	l = argv.split('=')
	if len(l) < 2: continue
	if l[0] == '-format': format = l[1]
	elif l[0] == '-range': range = l[1]
	elif l[0] == '-represent': repre = l[1]
	elif l[0] == '-subsample': subsp = int(l[1])
	elif l[0] == '-out': outname = l[1]

print 'Connecting to the NCBI FTP server...',
sys.stdout.flush()
ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov', 'anonymous', '')
# ftp.set_pasv(False)
print 'done.'

print 'Downloading the NCBI taxonomy database...',
sys.stdout.flush()
if os.path.isdir('taxdump'):
	if os.path.isfile('taxdump/names.dmp'): os.remove('taxdump/names.dmp')
	if os.path.isfile('taxdump/nodes.dmp'): os.remove('taxdump/nodes.dmp')
else:
	os.makedirs('taxdump')
ftp.cwd('/pub/taxonomy')
fout = open('taxdump.tar.gz', 'wb')
ftp.retrbinary('RETR ' + 'taxdump.tar.gz', fout.write, 1024)
fout.close()
tar = tarfile.open('taxdump.tar.gz', 'r:gz')
tar.extract('names.dmp', 'taxdump')
tar.extract('nodes.dmp', 'taxdump')
tar.close()
os.remove('taxdump.tar.gz')
print 'done.'

name2id = {}
nname = {} # number of names per TaxID
if repre == 'auto' or subsp:
	print 'Reading the NCBI taxonomy database...',
	sys.stdout.flush()
	fin = open('taxdump/names.dmp', 'r')
	for line in fin:
		l = re.split('\s+\|\s+', line)
		if nname.has_key(l[0]): nname[l[0]] += 1
		else: nname[l[0]] = 1
		if re.search('scientific name\s*\|', line):
			if not name2id.has_key(l[1]):
				name2id[l[1]] = l[0]
	fin.close()
	print 'done.'
	print ' ', str(len(nname)), 'taxa read.'

represents = {}
if repre == 'auto':
	print 'Downloading the NCBI representative genome list...',
	sys.stdout.flush()
	if os.path.exists('representative_genomes.txt'): os.remove('representative_genomes.txt')
	try:
		urllib.urlretrieve('http://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=refgenomes&download=on', 'representative_genomes.txt')
		print 'done.'
	except:
		print 'failed.'
	if os.path.exists('representative_genomes.txt'):
		print 'Reading the NCBI representative genome list...',
		sys.stdout.flush()
		fin = open('representative_genomes.txt', 'r')
		for line in fin:
			if line[0] == '#': continue
			line = line.rstrip('\r\n')
			if not line: continue
			name = line.split('\t')[2]
			if name2id.has_key(name):
				represents[name2id[name]] = 1
		fin.close()
		print 'done.'
		print ' ', str(len(represents)), ' genomes read.'
elif repre != 'no' and os.path.exists(repre):
	print 'Reading user-defined representative genome list...',
	sys.stdout.flush()
	fin = open(repre, 'r')
	for line in fin:
		if line[0] == '#': continue
		line = line.rstrip('\r\n')
		if not line: continue
		l = line.split('\t')
		represents[l[0]] = 1
	fin.close()
	print ' ', str(len(represents)), ' genomes read.'

print 'Reading RefSeq genome list...',
sys.stdout.flush()
gcfs = {}
downfiles = {}
if range == 'all':
	downfiles['assembly_summary_refseq.txt'] = 'assembly_summary_refseq.txt'
elif range == 'microbe':
	for i in ('archaea', 'bacteria', 'fungi', 'protozoa'):
		downfiles['assembly_summary_' + i + '.txt'] = i + '/assembly_summary.txt'
else:
	for i in range.split(','):
		downfiles['assembly_summary_' + i + '.txt'] = i + '/assembly_summary.txt'
for downfile in downfiles:
	if os.path.isfile(downfile): os.remove(downfile)
	if '/' in downfiles[downfile]: ftp.cwd('/genomes/refseq/' + downfiles[downfile].split('/')[0])
	else: ftp.cwd('/genomes/refseq')
	fout = open(downfile, 'wb')
	if '/' in downfiles[downfile]: ftp.retrbinary('RETR ' + downfiles[downfile].split('/')[1], fout.write, 1024)
	else: ftp.retrbinary('RETR ' + downfiles[downfile], fout.write, 1024)
	fout.close()
	fin = open(downfile, 'r')
	for line in fin:
		line = line.rstrip('\r\n')
		if line[0] == '#': continue
		l = line.split('\t')
		if l[10] == 'latest':
			gcfs[l[0]] = {'wgs_master': l[3], 'taxid': l[5], 'species_taxid': l[6], 'organism_name': l[7], 'assembly_level': l[11], 'seqs': 0, 'length': 0}
	fin.close()
	os.remove(downfile)
print 'done.'
print ' ', str(len(gcfs)), 'genomes read.'

if subsp:
	print 'Subsampling genomes...',
	sys.stdout.flush()
	sp2gcfs = {}
	for gcf in gcfs:
		sp = gcfs[gcf]['species_taxid']
		if sp2gcfs.has_key(sp): sp2gcfs[sp] += ',' + gcf
		else: sp2gcfs[sp] = gcf
	for sp in sp2gcfs:
		l = sp2gcfs[sp].split(',')
		if len(l) <= subsp:
			for gcf in l:
				gcfs[gcf]['in'] = 1
		else:
			n = 0
			if represents: # get representative organisms as priority
				repid2gcfs = {} # one representative TaxID may correspond to multiple GCFs
				for gcf in l:
					id = gcfs[gcf]['taxid']
					if represents.has_key(id):
						if repid2gcfs.has_key(id): repid2gcfs[id] += ',' + gcf
						else: repid2gcfs[id] = gcf
				for id in repid2gcfs:
					if ',' in repid2gcfs[id]:
						l = repid2gcfs[id].split(',')
						thisgcf = ''
						for gcf in l:
							if gcfs[gcf]['assembly_level'] == 'Complete Genome' or gcfs[gcf]['assembly_level'] == 'Chromosome':
								thisgcf = gcf
								break
						if not thisgcf:
							for gcf in l:
								if gcfs[gcf]['assembly_level'] == 'Scaffold':
									thisgcf = gcf
									break
							if not thisgcf:
								thisgcf = l[0]
						gcfs[thisgcf]['in'] = 1
					else:
						gcfs[repid2gcfs[id]]['in'] = 1
					n += 1
			if n >= subsp: continue
			names = {} # GCF => number of names
			for gcf in l:
				if gcfs[gcf]['assembly_level'] == 'Complete Genome' or gcfs[gcf]['assembly_level'] == 'Chromosome':
					id = gcfs[gcf]['taxid']
					if nname.has_key(id): names[gcf] = nname[id]
					else: names[gcf] = 0
			if names:
				for (gcf, number) in sorted(names.items(), reverse=True, key=lambda x: x[1]):
					gcfs[gcf]['in'] = 1
					n += 1
					if n >= subsp: break
			if n >= subsp: continue
			names = {}
			for gcf in l:
				if gcfs[gcf]['assembly_level'] == 'Scaffold':
					id = gcfs[gcf]['taxid']
					if nname.has_key(id): names[gcf] = nname[id]
					else: names[gcf] = 0
			if names:
				for (gcf, number) in sorted(names.items(), reverse=True, key=lambda x: x[1]):
					gcfs[gcf]['in'] = 1
					n += 1
					if n >= subsp: break
			if n >= subsp: continue
			names = {}
			for gcf in l:
				if gcfs[gcf]['assembly_level'] == 'Contig':
					id = gcfs[gcf]['taxid']
					if nname.has_key(id): names[gcf] = nname[id]
					else: names[gcf] = 0
			if names:
				for (gcf, number) in sorted(names.items(), reverse=True, key=lambda x: x[1]):
					gcfs[gcf]['in'] = 1
					n += 1
					if n >= subsp: break
	n = 0
	for gcf in gcfs:
		if gcfs[gcf].has_key('in'): n += 1
	print 'done.'
	print ' ', str(n), 'genomes retained.'

print 'Reading RefSeq genomic data...',
sys.stdout.flush()
ftp.cwd('/genomes/all')
files = []
try: files = ftp.nlst()
except ftplib.error_perm, resp:
	if str(resp) == '550 No files found': sys.exit('Error: Data source not available.')
	else: raise
dirs = []
for file in files:
	if not file.startswith('GCF_'): continue
	l = file.split('_')
	if len(l) < 3: continue
	gcf = l[0] + '_' + l[1]
	if gcfs.has_key(gcf):
		if subsp:
			if not gcfs[gcf].has_key('in'): continue
		dirs.append(file)
print 'done.'
print ' ', str(len(dirs)), 'genomes to download.'

print 'Extracting non-redundant proteomic data from NCBI RefSeq...'
sys.stdout.flush()
(gilist, gi2taxids, sharedgis) = ([], {}, {})
(ngcf, nseq, naa) = (0, 0, 0)
(nodata, haserror) = ([], [])
print ' ',
fout = open('sampled_genomes.txt', 'w')
fout.write('#' + '\t'.join(('assembly_accession', 'wgs_master', 'taxid', 'species_taxid', 'organism_name', 'assembly_level', 'proteins', 'residues')) + '\n')
fout.close()
for dir in dirs:
	l = dir.split('_')
	gcf = l[0] + '_' + l[1]
	ftp.cwd('/genomes/all/' + dir)
	files = []
	retry = 0
	while True:
		try:
			files = ftp.nlst()
			break
		except ftplib.error_perm, resp:
			if str(resp) == '550 No files found':
				nodata.append(gcf)
				break
			else:
				if retry >= 5:
					haserror.append(gcf)
					print str(resp)
					break
				else:  
					retry += 1
					continue
	downfile = dir + '_protein.gpff.gz'
	if not downfile in files:
		nodata.append(gcf)
		continue
	fout = open(downfile, 'wb')
	ftp.retrbinary('RETR ' + downfile, fout.write, 1024)
	fout.close()
	fout = open(outname + '.faa', 'a')
	fin = gzip.open(downfile, 'r')
	taxid = gcfs[gcf]['taxid']
	(accn, gi, name, seq) = ('', '', '', '')
	for line in fin:
		line = line.rstrip('\r\n')
		if line.startswith('DEFINITION  '):
			name = line[12:].rstrip('.')
		elif line.startswith('VERSION     '):
			l = line[12:].split(' ')
			(accn, gi) = (l[0], l[-1].split(':')[-1])
		elif len(line) > 10 and re.search(r'^\s+\d+\s$', line[:10]):
			seq += line[10:].replace(' ', '').upper()
		elif line == '//':
			if gi2taxids.has_key(gi):
				gi2taxids[gi] += ',' + taxid
				if not sharedgis.has_key(gi): sharedgis[gi] = '1'
			else:
				gi2taxids[gi] = taxid
				gilist.append(gi)
				fout.write('>gi|' + gi + '|ref|' + accn + '| ' + name + '\n' + seq +'\n')
			gcfs[gcf]['seqs'] += 1
			gcfs[gcf]['length'] += len(seq)
			(accn, gi, name, seq) = ('', '', '', '')
	fin.close()
	fout.close()
	nseq += gcfs[gcf]['seqs']
	naa += gcfs[gcf]['length']
	fout = open('sampled_genomes.txt', 'a')
	fout.write('\t'.join((gcf, gcfs[gcf]['wgs_master'], taxid, gcfs[gcf]['species_taxid'], gcfs[gcf]['organism_name'], gcfs[gcf]['assembly_level'], str(gcfs[gcf]['seqs']), str(gcfs[gcf]['length']))) + '\n')
	fout.close()
	os.remove(downfile)
	ngcf += 1
	if not ngcf % 25:
		print str(ngcf),
		sys.stdout.flush()
print 'done.'

print ' ', str(ngcf), 'genomes,', str(nseq), 'protein sequences, and', str(naa), 'residues downloaded.'
if nodata:
	print ' ', len(nodata), 'genomes have no data.'
	fout = open('genomes_without_data.txt', 'w')
	for i in nodata: fout.write(i + '\n')
	fout.close()
if haserror:
	print ' ', len(haserror), 'genomes cannot be downloaded.'
	fout = open('cannot_download.txt', 'w')
	for i in haserror: fout.write(i + '\n')
	fout.close()
fout = open('gi2taxid.txt', 'w')
for gi in gilist:
	fout.write(gi + '\t' + gi2taxids[gi] + '\n')
fout.close()
if sharedgis:
	print ' ', str(len(gi2taxids)), 'proteins extracted, including', str(len(sharedgis)), 'proteins shared by multiple organisms.'
else:
	print ' ', str(len(gi2taxids)), 'proteins extracted.'

if format == 'blast':
	print 'Building BLAST database...'
	sys.stdout.flush()
	if os.path.isdir('blast'):
		filelist = os.listdir('blast')
		if filelist:
			for filename in filelist:
				os.remove('blast/' + filename)
	else:
		os.makedirs('blast')
	p = subprocess.Popen('makeblastdb -in ' + outname + '.faa -parse_seqids -dbtype prot -out blast/' + outname + ' -title "Auto-generated standard BLAST database ' + str(today) + '"', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
	out = p.stdout.read()
	if 'command not found' in out:
		print 'failed. Check if ncbi-blast+ has been correctly installed.'
	else:
		for line in out.split('\n'):
			if not line: continue
			if line[0] == '\t': line = '  ' + line[1:]
			print '  ' + line
	p.stdout.close()
	print 'done.'
print 'Task completed.'


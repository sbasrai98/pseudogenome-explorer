import os
import pandas as pd
from Bio import SeqIO
from PIL import Image, ImageDraw, ImageFont
font = ImageFont.truetype("OpenSans-Regular.ttf", 60)
import pickle
with open('pseud_db.pickled', 'rb') as filein:
    pseud_db = pickle.load(filein)
megapath = '/mnt/c/Program\\ Files/MEGA-X/megacc.exe'

def search(query, pseud_db):
    '''
    Returns a DataFrame of all pseudogene families that contain the query in\
     the gene name or description (or top BLAST hits for sequence queries)
    query: string
    pseud_db: pd.DataFrame
    '''
    if query.startswith('-file '):
        filename = query.split(' ')[1]
        os.system('blastn -db pseudogenes_blastdb/pseudogenes -query %s \
            -out blasted.tsv -outfmt 6' % filename)
        os.system(
            "sed -i '1s/^/qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\
\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\\n/' blasted.tsv")
        results = pd.read_csv('blasted.tsv', sep='\t')['sseqid']
        results = list(map(lambda x: x.split('-')[1], results))
        families = []
        for id in results:
            df = pseud_db[pseud_db['gene'] == id]
            if len(df) > 0: 
                families.append(list(df['family'])[0])
        os.system("rm blasted.tsv")
        return list(dict.fromkeys(families))
    else:
        matches = []      # list of indices that contain a match
        for i in range(len(pseud_db)):
            row = pseud_db.loc[i]
            if query in row['gene'] or query in row['description']:
                matches.append(i)
        result = pseud_db.iloc[matches]
        return list(dict.fromkeys(result['family']))

def chrom_map(features, original):
    '''
    Returns a PIL image object representing a chromosomal map of all\
     pseudogene locations
    features: list of Bio.SeqFeature
    original: Bio.SeqFeature or None
    '''
    chrom_lens = [2489, 2421, 1982, 1902, 1815, 1708, 1593, 1451, 
                  1383, 1337, 1350, 1332, 1143, 1070, 1019, 903, 
                  832, 803, 586, 644, 467, 508, 1560, 572]
    im = Image.new('RGB', (2490+200, (25*140)), (255, 255, 255))
    draw = ImageDraw.Draw(im)   
    # group features by chromosome
    chromosomes = {}
    for i in range(24): 
        chromosomes[i] = [ ]
    for ft in features:
        chr = int(ft.qualifiers['from'][0].split('.')[0][-2:]) - 1
        chromosomes[chr].append(ft)
    # draw all chromosomes
    ystart = 70
    for i in range(24):
        draw.rectangle([150, ystart, 150+chrom_lens[i], ystart+50], 
                       outline=(0,0,0), fill=(245,205,105))
        draw.rectangle([150, ystart+50, 150+chrom_lens[i], ystart+100], 
                       outline=(0,0,0), fill=(245,205,105))
        for ft in chromosomes[i]:
            if ft.location.strand > 0: strand = ystart # positive sense
            if ft.location.strand < 0: strand = ystart+50 # negative sense
            if ft.qualifiers['gene_biotype'][0] == 'transcribed_pseudogene':
                color = (0,0,255)
            else:
                color = (255,0,0)
            pos = int(round((ft.location.start + ft.location.end) / 200000)) 
            draw.rectangle([150+(pos-4), strand, 150+(pos+4), strand+50], 
                        outline=(0,0,0), fill=color)
        label = str(i+1)
        if i < 9: label = ' '+label
        if i == 22: label = ' X'
        if i == 23: label = ' Y'
        draw.text([50, ystart+10], label, fill=(0,0,0),
                  font=font)
        ystart += 140
    if original:  # if an original gene has been specified
        chr = int(original.qualifiers['from'][0].split('.')[0][-2:]) - 1
        pos = int(round(
            (original.location.start + original.location.end) / 200000))
        ystart = 70 + (chr*140)
        if original.location.strand < 0: ystart += 50 # negative sense 
        draw.rectangle([150+(pos-4), ystart, 150+(pos+4), ystart+50], 
                    outline=(0,0,0), fill=(0,255,0))
    return im

def get_phylo(family):
    '''
    Executes shell commands to align sequences and construct phylogeny\
     using MEGA X
    family: string
    '''
    os.system(megapath+' -a muscle_align.mao -d "'+family+' sequences.fa" -o "'+family+'"')
    os.system(megapath+' -a ML_phylo.mao -d "'+family+'.meg" -o "'+family+'"')

def get_info(features, original):
    '''
    Returns a DataFrame containing information about all pseudogenes provided
    features: list of Bio.SeqFeature
    original: Bio.SeqFeature or None
    '''
    d = {'gene':[], 'description':[], 'source':[], 'location':[],
         'biotype':[]}
    all_entries = list(features)
    if original:
        all_entries = [original] + list(features)
    for ft in all_entries:
        d['gene'].append(ft.qualifiers['gene'][0])
        d['description'].append(ft.qualifiers['description'][0])
        d['source'].append(ft.qualifiers['from'][0])
        d['location'].append(ft.location)
        d['biotype'].append(ft.qualifiers['gene_biotype'][0])
    df = pd.DataFrame(d)
    return df

def get_fasta(family_df, orig_feature):
    '''
    Returns all FASTA sequences for the given pseudogene family
    family_df: pd.DataFrame
    orig_feature: Bio.SeqFeature or None
    '''
    fam_seqs = []
    fam_ids = list(map(lambda x: 
                family_df.loc[x]['pseud_feature'].qualifiers['from'][0] + \
                '-' + family_df.loc[x]['gene'], family_df.index))
    seqs = SeqIO.parse('all_pseudogenes.fa', 'fasta')
    for seq in seqs:
        if seq.id in fam_ids:
            fam_seqs.append(seq)
        if len(fam_seqs) == len(fam_ids):
            break
    if orig_feature:
        orig_id = orig_feature.qualifiers['from'][0] + '-' + \
                orig_feature.qualifiers['gene'][0] 
        seqs = SeqIO.parse('all_genes.fa', 'fasta')
        for seq in seqs:
            if seq.id == orig_id:
                fam_seqs = [seq] + fam_seqs
                break
    return fam_seqs

def main():
    '''
    Executes the main program
    '''
    query = input('Search for a pseudogene family by gene name (ex. DDX11L1),\
 description (ex. helicase), or FASTA sequence (ex. -file gene.fa).\nQuery: ')
    results = search(query, pseud_db)
    if len(results) == 0:
        print('Sorry, no matches found.')
        main()
    else:
        for i, fam in enumerate(results):
            print(i, fam)
        valid = list(map(str, range(len(results)))) + ['q']
        select = 'a'
        while not(select in valid):
            select = input('Select a pseudogene (0-%i) or q to return: ' \
                            % (len(results)-1))
        if select == 'q':
            main()
        else:
            family = results[int(select)]
            family_df = pseud_db[pseud_db['family'] == family]
            pseudogenes = family_df['pseud_feature']
            original = list(family_df['orig_feature'])[0]
            print('Selection:', family)
            print(' - %i pseudogene(s) found' % len(pseudogenes))
            if original:
                print(' - Original gene found')
            else:
                print(' - Original gene not found')
            ## What data would you like to view?
            valid_options = ['m', 'p', 'i', 'f', 'q']
            options = ['a']
            while sum(map(lambda x: not(x in valid_options), options)) != 0:
                options = input('Options:\n\
 - m: Map pseudogene locations on chromosomes\n\
 - p: Construct a phylogeny with pseudogene sequences\n\
 - i: Save information about this pseudogene family\n\
 - f: Save original FASTA sequences of this pseudogene family\n\
Select actions to perform (ex. mpf) or q to return: ')
                options = list(options)
            if 'q' in options:
                main()
            if 'm' in options:
                print('Mapping pseudogene locations...')
                im = chrom_map(pseudogenes, original)
                im.save(family+' map.png', quality=100)
            if 'i' in options:
                print('Saving pseudogene information...')
                output = get_info(pseudogenes, original)
                output.to_csv(family+' info.tsv', sep='\t', index=False)
            if 'f' in options or 'p' in options:
                print('Saving FASTA sequences...')
                SeqIO.write(get_fasta(family_df, original), 
                            family+' sequences.fa', 'fasta')
            if 'p' in options:
                get_phylo(family)


if __name__ == "__main__":
    main()


# returns an object containing all the proteins + metadata in the fasta file,  
# accessed by the stable protein ID or accession in case of artificial identifier
# all_proteins[proteinID] = {'tag': tag, 'accession': accession, 'description': description, 'sequence': sequence}
# if fasta headers are stored in a separate file, provide the data frame with the headers as the second argument, with the protein accession as index
def read_fasta(filename, headers=None):
    proteinDB_file = open(filename, 'r')

    # read protein db
    metadata = proteinDB_file.readline()    # line starting with '>'
    sequence = ""                           # all the other lines are considered protein sequences (considering also a multi-line format)

    all_proteins = {}

    # read the reference sequence database
    while metadata != "":

        line = proteinDB_file.readline()
        while not line.startswith('>') and not line == "":
            sequence += line[:-1]
            line = proteinDB_file.readline()

        tag = ''
        accession = ''
        description = ''

        if "|" in metadata:					# the header is at least partially formated
            metadata_parsed = metadata[1:].split('|')
            if (headers is not None):
                tag = headers.loc[accession]['tag']
            elif 'generic' in metadata_parsed[0]:
                tag = metadata_parsed[0]
            else:
                tag = 'generic_' + metadata_parsed[0]

            if len(metadata_parsed) == 2:					# accession and description are potentially merged -> separate them
                if " " in metadata_parsed[1]:
                    accession = metadata_parsed[1].split(' ')[0]
                    description = metadata_parsed[1].split(' ', 1)[1]
                else:
                    accession = metadata_parsed[1]				# no description -> keep accession as it is
            elif len(metadata_parsed) == 3:					# descripton and accesson already separated -> keep
                accession = metadata_parsed[1]
                description = metadata_parsed[2]

        else:						                    # the header is not formated
            accession = metadata[1:-1].split(" ")[0]
            if " " in metadata:
                description = metadata.split(" ", 1)[1]
            if (headers is not None):
                tag = headers.loc[accession]['tag']

        proteinID = accession.split('.')[0]

        matching_proteins = []
        seq_positions = []
        reading_frames = []

        if 'position_within_protein:' in description:
            seq_positions = [ int(pos) for pos in description.split('position_within_protein:', 1)[1].split(' ', 1)[0].split(';') ]
        elif (headers is not None):
            seq_positions = [ int(pos) for pos in headers.loc[accession]['position_within_protein'].split(';') ]

        if 'matching_proteins:' in description:
            proteinLists = description.split('matching_proteins:', 1)[1].split(' ', 1)[0].split(';')
            matching_proteins = [ l.split(',') for l in proteinLists ]
        elif (headers is not None):
            matching_proteins = [ l.split(',') for l in headers.loc[accession]['matching_proteins'].split(';') ]

        if 'reading_frame:' in description:
            rfLists = description.split('reading_frame:', 1)[1].split(maxsplit=1)[0].split(';')
            reading_frames = [ l.split(',') for l in rfLists ]
        elif (headers is not None):
            reading_frames = [ l.split(',') for l in headers.loc[accession]['reading_frame'].split(';') ]

        all_proteins[proteinID] = {'tag': tag, 'accession': accession, 'description': description.replace('\n', ''), 'sequence': sequence, 'seq_positions': seq_positions, 'matching_proteins': matching_proteins, 'reading_frames': reading_frames}

        metadata = line
        sequence = ""

    proteinDB_file.close()

    return all_proteins

# returns an object containing all the protein metadata in the fasta file, 
# accessed by the stable protein ID or accession in case of artificial identifier
# if fasta headers are stored in a separate file, provide the data frame with the headers as the second argument, with the protein accession as index
def read_fasta_metadata(filename, headers=None):
    proteinDB_file = open(filename, 'r')

    # read protein db
    line = proteinDB_file.readline()
    all_proteins = {}

    while line != '':
        if (not line.startswith('>')):
            line = proteinDB_file.readline()  
        else:
            tag = ''
            accession = ''
            description = ''

            if "|" in line:					# the header is at least partially formated
                metadata_parsed = line[1:].split('|')
                if (headers is not None):
                    tag = headers.loc[accession]['tag']
                elif 'generic' in metadata_parsed[0]:
                    tag = metadata_parsed[0]
                else:
                    tag = 'generic_' + metadata_parsed[0]

                if len(metadata_parsed) == 2:					# accession and description are potentially merged -> separate them
                    if " " in metadata_parsed[1]:
                        accession = metadata_parsed[1].split(' ')[0]
                        description = metadata_parsed[1].split(' ', 1)[1]
                    else:
                        accession = metadata_parsed[1]				# no description -> keep accession as it is
                elif len(metadata_parsed) == 3:					# descripton and accesson already separated -> keep
                    accession = metadata_parsed[1]
                    description = metadata_parsed[2]

            else:						                    # the header is not formated
                accession = line[1:].split(" ")[0]
                if " " in line:
                    description = line.split(" ", 1)[1]
            if (headers is not None):
                tag = headers.loc[accession]['tag']

            proteinID = accession.split('.')[0]
            all_proteins[proteinID] = {'tag': tag, 'accession': accession, 'description': description.replace('\n', '')}

            line = proteinDB_file.readline()

    proteinDB_file.close()
    return all_proteins

def get_protein_name_dict(fasta_file):
    all_proteins = read_fasta(fasta_file)
    name_dict = {}

    for proteinID in all_proteins:
        if (proteinID.startswith('ENSP')):
            name_dict[proteinID] = proteinID + ':REF'
        elif (proteinID.startswith('enshap')):
            ref_id = all_proteins[proteinID]['description'].split('ref:')[1].split('.')[0]
            hap = all_proteins[proteinID]['description'].split('haplotype:')[1].split()[0]
            name_dict[proteinID] = ref_id + ':' + hap
        else:
            name_dict[proteinID] = proteinID
    
    return name_dict, all_proteins

# returns a list of tryptic peptides, allowed number of missed cleavages provided in argument
# peptide length threshold in argument
def digest(seq, missed_c, min_length, max_length, cleavage_pattern):
    buffer = ""
    queue = []          # queue of tryptic peptides of any length
    position_queue = [] # starting positions of these peptides within given sequence

    peptides = []
    positions = []
    missed_cl = []

    prev_cleavage = 0

    for aa_idx in range(len(seq)-1):
        # check if there's a residue preventing cleavage
        if (seq[aa_idx+1] in cleavage_pattern['restrictionAfter']) or ((aa_idx > 0) and (seq[aa_idx-1] in cleavage_pattern['restrictionBefore'])):
            buffer += seq[aa_idx]
            continue

        # check whether we shouldn't cleave before adding the current residue
        if not (seq[aa_idx] in cleavage_pattern['aminoAcidAfter']):
            buffer += seq[aa_idx]

        if (seq[aa_idx] in cleavage_pattern['aminoAcidBefore']) or (seq[aa_idx] in cleavage_pattern['aminoAcidAfter']):
            queue.append(buffer)
            position_queue.append(prev_cleavage)

            if (len(queue) > (missed_c+1)):
                queue = queue[-(missed_c+1):]
                position_queue = position_queue[-(missed_c+1):]

            if (len(buffer) >= min_length and len(buffer) <= max_length):
                peptides.append(buffer)
                positions.append(prev_cleavage)
                missed_cl.append(0)

            for no_missed in range(1, missed_c+1):
                if (len(queue) > no_missed):
                    pep_missed = ''.join(queue[-(no_missed+1):])
                    if (len(pep_missed) >= min_length and len(pep_missed) <= max_length):
                        peptides.append(pep_missed)
                        positions.append(position_queue[-(no_missed+1)])
                        missed_cl.append(no_missed)

            buffer = "" if not (seq[aa_idx] in cleavage_pattern['aminoAcidAfter']) else seq[aa_idx]
            prev_cleavage = aa_idx+1 if not (seq[aa_idx] in cleavage_pattern['aminoAcidAfter']) else aa_idx

    # store the remaining peptide on the c-terminal
    queue.append(buffer)
    position_queue.append(prev_cleavage)

    if (len(queue) > (missed_c+1)):
        queue = queue[-(missed_c+1):]
        position_queue = position_queue[-(missed_c+1):]

    if (len(buffer) >= min_length and len(buffer) <= max_length):
        peptides.append(buffer)
        positions.append(prev_cleavage)
        missed_cl.append(0)

    for no_missed in range(1, missed_c+1):
        if (len(queue) > no_missed):
            pep_missed = ''.join(queue[-(no_missed+1):])
            if (len(pep_missed) >= min_length and len(pep_missed) <= max_length):
                peptides.append(pep_missed)
                positions.append(position_queue[-(no_missed+1)])
                missed_cl.append(no_missed)

    return peptides, positions, missed_cl

def get_protein_coverage(data):
    regions = []
    current_pos = 0
    region_start = 0
    current_type = -1

    active_regions = []
    eventQ = [] if (data[0][0] == 0) else [{'pos': 0, 'event': 'start', 'type': -1, 'pep_idx': -1}, {'pos': data[0][0], 'event': 'end', 'type': -1, 'pep_idx': -1}]

    for idx, pep in enumerate(data):
        eventQ.append({ 'pos': pep[0], 'event': 'start', 'type': pep[2], 'pep_idx': idx })
        eventQ.append({ 'pos': pep[0] + pep[1], 'event': 'end', 'type': pep[2], 'pep_idx': idx })

    eventQ = sorted(eventQ, key=lambda x: x['pos'])

    for evt in eventQ:
        current_pos = evt['pos']

        if evt['event'] == 'start':
            if evt['type'] > current_type:
                if current_pos > region_start:
                    regions.append([region_start, current_pos, current_type])
                current_type = evt['type']
                region_start = current_pos
            active_regions.append([evt['pep_idx'], evt['type']])

        if evt['event'] == 'end':
            active_regions = list(filter(lambda x: x[0] != evt['pep_idx'], active_regions))
            new_type = -1
            if (len(active_regions) > 0):
                new_type = max(list(map(lambda x: x[1], active_regions)))
            if (new_type < current_type):
                if current_pos > region_start:
                    regions.append([region_start, current_pos, current_type])
                region_start = current_pos
                current_type = new_type

    return regions


def ORF_gff_reader(filename):
    """
    i0_c100027/f6p76/1257|scaffold_2:3387555-3389419        .       CDS     2       613     .       +       .       ID=cds.m.31872;Parent=m.31872
    """
    last_id = None
    last_list = []
    with open(filename) as f:
        for line in f:
            raw = line.strip()
            if len(raw)==0: continue
            raw = raw.split()
            if raw[2]=='CDS':
                id = raw[0]
                start = int(raw[3])
                end = int(raw[4])
                if last_id == id: last_list.append((start, end))
                else:
                    if last_id is not None: yield last_id, last_list
                    last_id = id
                    last_list = [(start, end)]
    yield last_id, last_list

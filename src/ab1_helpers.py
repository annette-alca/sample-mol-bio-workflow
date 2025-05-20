"""These are several helper functions I wrote to be used in MolBio and Bioinfo workflows.
These are imported by the Jupyter notebooks (.ipynb) files. 

Do not remove this file from the src folder. Do not rename the src folder.
-Annette"""

from Bio import SeqIO
import csv, os, glob, time
from Bio.Blast import NCBIWWW, NCBIXML


def ab1_reader(ab1file):
    abrecord = SeqIO.read(ab1file,'abi')
    return abrecord

def mult_ab1_to_dict(ab1_subfolder):
    seqdict = {}
    ab1_path = os.getcwd() + '/seq_data/' + ab1_subfolder
    try:
        ab1files = glob.glob(ab1_path + "/*.ab1")
    except FileNotFoundError:
        print("""Cannot find ab1 subfolder, please make sure {ab1_subfolder} is in:
            \n{path[:len(path)-len(ab1_subfolder)]}.""")
        return
    for ab1file in ab1files:
        sample_id = ab1file[len(ab1_path)+1:-4] #get file name
        record = ab1_reader(ab1file)
        seqdict[sample_id] = record.seq
    return seqdict

def mult_ab1_to_single_fasta(ab1_subfolder, filename = "my_ab1s.fasta"):
    '''saves one fasta from ab1 folder'''
    if ".fa" not in filename:
        filename += ".fasta"
    output_folder = f"python_outputs/{ab1_subfolder}"
    try:
        os.mkdir(output_folder)
    except FileExistsError:
        pass

    fasta_path = f"python_outputs/{ab1_subfolder}/{filename}"
    seqdict = mult_ab1_to_dict(ab1_subfolder)
    with open (fasta_path,'w', encoding='utf-8') as f:
        for sample_id, sequence in seqdict.items():
            f.write (f">{sample_id}\n")
            f.write(str(sequence) + '\n')
    print(f'{fasta_path} saved with {len(seqdict)} sequences.')
    return seqdict, output_folder



def blast_seq(seq_obj):
    '''uses qblast from Bio.Blast.NCBIWWW'''
    result = NCBIWWW.qblast("blastn","nt", seq_obj, format_type="XML")
    return result

def blast_seqdict(seqdict,output_folder):
    try:
        os.listdir(output_folder)
    except FileNotFoundError:
        print(f"""folder {output_folder} does not exist. 
              Please create it or name a different folder.""")
        return
    for sample_id, sequence in seqdict.items():
        try:
            result = blast_seq(sequence)
            xml_path = f"{output_folder}/{sample_id}_blast.xml"
            with open(xml_path, "w") as out_handle:
                out_handle.write(result.read())
            
            print(f"{sample_id} qblast complete, result saved to {xml_path}")

        except Exception as e:
            print(f"Failed to BLAST {sample_id}: {e}")
        time.sleep(3)

def xml_to_table(seqdict,output_folder):
    '''convert XML to single csv file'''
    csvfile = output_folder + "/summary_of_blast.csv"
    xmls = glob.glob(output_folder + "/*.xml")

    with open(csvfile, "w") as out_file:
        out_file.write("Sample_ID,Sample_DNA_Length,Accession_of_Top_Hit,Title,Identity(%),Alignment Length,NCBI_Link\n")

        for xml_file in xmls:
            sample_id = xml_file[len(output_folder)+1:-10]
            dna_length = len(seqdict[sample_id])
            with open(xml_file) as handle:
                blast_record = NCBIXML.read(handle)

                if blast_record.alignments:
                    alignment = blast_record.alignments[0]
                    hsp = alignment.hsps[0]
                    accession = alignment.accession
                    title = alignment.title.replace(",", ";")  # prevent CSV breaking
                    identity = (hsp.identities / hsp.align_length) * 100
                    align_len = hsp.align_length
                    link = f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"

                    out_file.write(f"{sample_id},{dna_length},{accession},{title},{identity:.2f},{align_len},{link}\n")
                else:
                    out_file.write(f"{sample_id},No hits found,NA,NA,NA,NA\n")
    print (f"{csvfile} is saved with a summary of Blast results.")
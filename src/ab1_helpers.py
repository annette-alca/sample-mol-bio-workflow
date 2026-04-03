"""These are several helper functions I wrote to be used in MolBio and Bioinfo workflows.
These are imported by the Jupyter notebooks (.ipynb) files. 

Do not remove this file from the src folder. Do not rename the src folder.
-Annette"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv, os, glob, time, io
from Bio.Blast import NCBIWWW, NCBIXML


def ab1_reader(ab1file: str) -> SeqRecord:
    """Read a single .ab1 Sanger sequencing file.

    Args:
        ab1file: Path to the .ab1 file.

    Returns:
        A BioPython SeqRecord containing the sequence and metadata.
    """
    abrecord = SeqIO.read(ab1file,'abi')
    return abrecord

def mult_ab1_to_dict(ab1_subfolder: str) -> dict[str, Seq] | None:
    """Read all .ab1 files in a subfolder and return a dictionary of sequences.

    Args:
        ab1_subfolder: Name of the subfolder within seq_data/ containing .ab1 files.

    Returns:
        A dictionary mapping sample IDs (filenames without extension) to Seq objects,
        or None if the subfolder is not found.
    """
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

def mult_ab1_to_single_fasta(ab1_subfolder: str, filename: str = "my_ab1s.fasta") -> tuple[dict[str, Seq], str]:
    """Convert all .ab1 files in a subfolder to a single FASTA file.

    Args:
        ab1_subfolder: Name of the subfolder within seq_data/ containing .ab1 files.
        filename: Name for the output FASTA file. Defaults to "my_ab1s.fasta".

    Returns:
        A tuple of (sequence dictionary, output folder path).
    """
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



def blast_seq(seq_obj: Seq) -> io.StringIO:
    """Run a BLAST search for a single sequence against NCBI's nucleotide database.

    Args:
        seq_obj: A BioPython Seq object to search.

    Returns:
        A file-like handle containing XML-formatted BLAST results.
    """
    result = NCBIWWW.qblast("blastn","nt", seq_obj, format_type="XML")
    return result

def blast_seqdict(seqdict: dict[str, Seq], output_folder: str) -> None:
    """Run BLAST searches for all sequences in a dictionary and save results as XML.

    Includes a 3-second delay between queries to respect NCBI rate limits.

    Args:
        seqdict: Dictionary mapping sample IDs to Seq objects.
        output_folder: Path to the folder where XML result files will be saved.
    """
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

def xml_to_table(seqdict: dict[str, Seq], output_folder: str) -> None:
    """Parse BLAST XML results and produce a CSV summary of top hits.

    Reads all XML files in the output folder and writes a summary_of_blast.csv
    with the top hit for each sequence (accession, title, identity %, etc.).

    Args:
        seqdict: Dictionary mapping sample IDs to Seq objects (used for DNA length).
        output_folder: Path to the folder containing XML result files.
    """
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
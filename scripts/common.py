# This file is part of StaG-mwc
import csv
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
 
class UserMessages():
    """
    Deferred user messages for delayed printout
    """

    def __init__(self):
        self.messages = {
            "info": set(),
            "warning": set(),
        }

    def print_messages(self):
        for level, messages in self.messages.items():
            for message in messages:
                print(level.upper()+":", message)

    def print_info(self):
        for message in self.info:
            print(level.upper()+":", message)

    def print_warnings(self):
        for message in self.warn:
            print(level.upper()+":", message)

    def info(self, message):
        self.messages["info"].add(message)

    def warn(self, message):
        self.messages["warning"].add(message)

        
class SampleSheet():
    """Parse three-column samplesheet with header and column names:
    sample_id -- Unique sample ID 
    fastq_1   -- (Local/Remote) path to fastq_1
    fastq_2   -- (Local/Remote) path to fastq_2

    Example:
    sample_id\tfastq_1\t_fastq_2
    UniqueSampleName1\tinput/sample1_1.fq.gz\tinput/sample1_2.fq.gz
    OtherSample2\ts3://bucket/QWERTY12345_1.fq.gz\ts3://bucket/QWERTY12345_2.fq.gz
    """
    def __init__(self, samplesheet, endpoint_url="http://s3.amazonaws.com"):
        self.sample_info = {}
        self.S3 = S3RemoteProvider(host=endpoint_url)

        with open(samplesheet) as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter="\t")
            for row in reader:
                fq1_source = "local"  # For debugging
                fq2_source = "local"  # For debugging
                fq1 = row["fastq_1"]
                fq2 = row["fastq_2"]
                if fq1.startswith("s3://"):                   
                    fq1 = self.S3.remote(fq1.lstrip("s3://"))      
                    fq1_source = "S3"                         
                if fq2.startswith("s3://"):                   
                    fq2 = self.S3.remote(fq2.lstrip("s3://"))      
                    fq2_source = "S3"                         
                                                              
                self.sample_info[row["sample_id"]] = {                 
                        "read1": fq1, "read1_src": fq1_source,
                        "read2": fq2, "read2_src": fq2_source,
                }                                             
        self.samples = list(self.sample_info.keys())



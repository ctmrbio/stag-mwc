# This file is part of StaG-mwc
import csv
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
 
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
    def __init__(self, samplesheet, keep_local=False, endpoint_url="http://s3.amazonaws.com"):
        self.sample_info = {}
        self.S3 = S3RemoteProvider(host=endpoint_url)
        self.HTTP = HTTPRemoteProvider()
        self.required_columns = ("sample_id", "fastq_1", "fastq_2")

        def create_provider(uri):
            if uri.startswith("s3://"):
                return self.S3.remote(uri.lstrip("s3://"), keep_local=keep_local)
            if uri.startswith(("http://", "https://")):
                return self.HTTP.remote(uri.split("//", maxsplit=1)[1], keep_local=keep_local)
            else:
                return uri

        with open(samplesheet) as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter="\t")
            for line, row in enumerate(reader, start=1):
                if not all(col in row for col in self.required_columns):
                    raise ValueError(f"Missing column! {self.required_columns} are required.")
                if row["sample_id"] in self.sample_info:
                    raise ValueError(f"{row['sample_id']} exists more than once!")

                try:
                    fq1 = create_provider(row["fastq_1"])
                    fq2 = create_provider(row["fastq_2"])
                except AttributeError as e:
                    raise ValueError(f"Cannot parse line {line} in {samplesheet}")
                                                              
                self.sample_info[row["sample_id"]] = {                 
                        "read1": fq1,
                        "read2": fq2,
                }                                             
        self.samples = list(self.sample_info.keys())

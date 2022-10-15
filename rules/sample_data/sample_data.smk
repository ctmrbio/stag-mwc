# vim: syntax=python expandtab
# Downloads sample data to .input as to keep separate from other input files

import csv
import os

sampledatadir = ".input"
samples_ = []
if os.path.exists('.input/sample_data.csv'):
    with open('.input/sample_data.csv', 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            sample_ = "_".join([row["sample_id"], row["read"]]) + "." + row["suffix"]
            samples_.append(f"{sampledatadir}/{sample_}")

    rule all:
        input:
            samples_

    rule curl_sampledata:
        input: ".input/sample_data.csv",
        output: samples_,
        shell:
            """
            wget --no-check-certificate --content-disposition https://github.com/boulund/stag-mwc_test_data/archive/master.zip -P {sampledatadir}
            cd {sampledatadir}
            unzip -j stag-mwc_test_data-master.zip
            rm README.md
            rm stag-mwc_test_data-master.zip
            """

else:
        print(f"no sample data sheet in ./input")

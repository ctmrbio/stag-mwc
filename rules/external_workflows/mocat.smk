# vim: syntax=python expandtab
# Rules for MOCAT2 used inside of mwc

rule download_mocat2_db:
    output:
        config["dbdir"]+"/mocat/"
    shadow:
        "shallow"
    shell:
        """
        wget http://vm-lux.embl.de/~kultima/share/MOCAT/data/1506MG.tar.gz \
        && \
        tar -xf 1506MG.tar.gz \
        """

rule mocat2:
    input:
        sample_file=config["mocat2_sample_file"],
        config_file=config["mocat2_config_file"],
    output:
    shadow:
        "shallow"
    threads:
        20
    singularity:
        "docker://boulund/mocat2"
    shell:
        """
        export EMAIL={config[mocat2_email]}
        runMOCAT.sh \
            -sf {sample_file} \
            -cfg {config_file}
        """

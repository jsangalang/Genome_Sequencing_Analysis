include: "mutect2.smk"
include: "filterMC.smk"
include: "../Common/collectSeqAM.smk"
include: "../Common/estiContamination.smk"
include: "filterOB.smk"

if config["pipeline"]["SEQ"] == "humain":
    include: "Split/oncotator.smk"
    include: "Exom/oncotator.smk"
else:
    include: "Annovar/annovar.smk"

include: "mutect2.smk"
include: "filterVC.smk"
include: "../Common/collectSeqAM.smk"
include: "../Common/estiContamination.smk"
include: "filterOB.smk"

if config["pipeline"]["SEQ"] == "humain":
    include: "Split/main.smk"
    include: "Exom/main.smk"

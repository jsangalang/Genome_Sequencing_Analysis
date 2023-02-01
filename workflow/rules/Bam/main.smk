if config["pipeline"]["PAIRED"]:
    include: "paired.smk"
else:
    include: "single.smk"

include: "bam.smk"


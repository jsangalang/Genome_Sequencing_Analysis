if config["pipeline"]["PAIR"] == "PE":
    include: "paired.smk"
else:
    include: "single.smk"

"bam.smk"


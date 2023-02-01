if config["pipeline"]["PAIRED"]:
    include: "fastpPE.smk"
else:
    include: "fastpSE.smk"

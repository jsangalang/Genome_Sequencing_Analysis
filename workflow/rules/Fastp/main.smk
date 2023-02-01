if config["pipeline"]["PAIR"] == "PE":
    include: "fastpPE.smk"
else:
    include: "fastpSE.smk"

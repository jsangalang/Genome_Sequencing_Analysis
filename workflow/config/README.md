#### Note 1: default values for target interval

TARGET_INTERVAL_GATK = ""

TARGET_INTERVAL_BQSR = " -L 1 "

#### Note 2: any configuration of DIR should not be terminated by "/"

#### Note 3: settings of pipeline

samples: ["human"|"mouse"]

seq_type: ["WGS"|"WES"]

paired: [True|False]

remove_duplicates: [True|False]

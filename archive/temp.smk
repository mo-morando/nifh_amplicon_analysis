# Import the datetime module
import datetime


# Function to generate a timestamp
def get_timestamp():
    return datetime.datetime.now().strftime("%Y%m%d")

configfile: "config/config_blast_db_ucyna.json"
configfile: "config/config_blast_db_arb.json"

configfile: "config/config_blast_ucyna.json"
configfile: "config/config_blast_arb.json"

configfile: "config/config_filt_blast_ucyna.json"
configfile: "config/config_filt_blast_arb.json"

# List of databases or config files
databases = ["arb", "ucyna"]


########################


configfile: "./config_change_csv_header_ucyna.json"
configfile: "./config_change_csv_header_arb.json"


# Define a rule to collect all outputs
rule all:
    input:
        expand(
            "{blast_output_dir}{output}",
            blast_output_dir=[config[f"blast_output_dir_{db}"] for db in databases],
            output=[config[f"output_header_csv_{db}"] for db in databases]
        )


# Dynamically create rules for each database
for database in databases:
    # configfile: f"config/config_change_csv_{database}.json"
    blast_output_dir = config[f"blast_output_dir_{database}"]
    # print(blast_output_dir)

    rule:
        name: f"database_headers_{database}"
        input:
            # config["input_csv"],
            expand(
                "{blast_output_dir}{output}",
                blast_output_dir=blast_output_dir,
                output=config[f"filtered_blast_output_{database}"],
            ),
        params:
            header_prefix=config[f"header_prefix_{database}"],
        output:
            expand(
                "{blast_output_dir}{output}",
                blast_output_dir=blast_output_dir,
                output=config[f"output_header_csv_{database}"],
            ),
        log:
            expand(
                "{blast_output_dir}logs/headers_{timestamp}.log",
                timestamp=get_timestamp(),
                blast_output_dir=blast_output_dir,
            ),
        shell:
            """
            bash ../../annotations/scripts/change_csv_headers.sh {input} {output} {params.header_prefix} > {log} 2>&1
            """







# ## Rule to add AUID key
# ## At some point, new AUIDs were generated but all my files and analysis need these old AUID IDs, so this key needs to be added

configfile: "./config_add_auid_key_ucyna.json"
configfile: "./config_add_auid_key_arb.json"

# Define a rule to collect all outputs
rule all_add_auid_key:
        input:
            expand(
                "{blast_output_dir}{output}",
                blast_output_dir=[config[f"blast_output_dir_{db}"] for db in databases],
                output=[config[f"keyed_csv_{db}"] for db in databases]
            )
            # expand(
            #     "{blast_output_dir}{output}",
            #     blast_output_dir=config[f"blast_output_dir{database}"],
            #     output=config[f"keyed_csv_{database}"],
            #     database=databases
            # )

# Dynamically create rules for each database
for database in databases:
    # configfile: f"config/config_change_csv_{database}.json"
    blast_output_dir = config[f"blast_output_dir_{database}"]
    # print(blast_output_dir)
    
    rule:
        name:f"add_auid_key{database}"
        input:
            expand(
                "{blast_output_dir}{output}",
                blast_output_dir=blast_output_dir,
                output=config[f"output_header_csv_{database}"],
            ),
        params:
            left_on_key=config[f"left_on_key_{database}"],
            right_on_key=config[f"right_on_key_{database}"],
            how_key=config[f"how_key_{database}"],
        output:
            expand(
                "{blast_output_dir}{output}",
                blast_output_dir=blast_output_dir,
                output=config[f"keyed_csv_{database}"],
            ),
        log:
            # "logs/filter_blast{timestamp}.log",
            expand(
                "{blast_output_dir}logs/auid_key_{timestamp}.log",
                timestamp=get_timestamp(),
                blast_output_dir=blast_output_dir,
            ),
        shell:
            """
            python3 ../../annotations/scripts/add_auid_key.py --left_csv {input} --merged_csv {output} --left_on {params.left_on_key} --right_on {params.right_on_key} --how {params.how_key} >> {log} 2>&1
            """




######

## Rule to merge csv files
configfile: "./config/config_merge_arb_genome879.json"
configfile: "./config/config_merge_cart.json"

# Define a rule to collect all outputs
rule all_merge:
    input:
        expand(
            "{merged_dir}{output}",
            merged_dir = [config[f"merged_dir_{db}"] for db in databases],
            output=[config[f"merged_csv_{db}"] for db in databases],
        ),

# Define merging order with merge_key
merge_key = config["merge_key"]
# Dynamically creat rules for each database
for database in merge_key:
    blast_output_dir = config[f"blast_output_dir_{database}"]
    ## Define output directory
    merged_dir = config[f"merged_dir_{database}"]


    rule: 
        name: f"merge{database}"
        input:
            expand(
                "{blast_output_dir}{output}",
                blast_output_dir=blast_output_dir,
                output=config[f"keyed_csv_{database}"],
            ),
            # config["right_csv"],
        params:
            left_csv=config[f"left_csv_{database}"],
            # left_on=config[f"left_on_{database}"],
            left_on=config[f"left_on_key_{database}"],
            right_on=config[f"right_on_{database}"],
            how=config[f"how_{database}"],
        output:
            # config["merged_csv"],
            expand(
                "{merged_dir}{output}",
                merged_dir=merged_dir,
                output=config[f"merged_csv_{database}"],
            ),
        log:
            expand(
                "{merged_dir}logs/merge_{timestamp}.log",
                timestamp=get_timestamp(),
                merged_dir=merged_dir,
            ),
        shell:
            """
            python3 ../scripts/merge_csv.py --left_csv {params.left_csv} --right_csv {input} --merged_csv {output} --left_on {params.left_on} --right_on {params.right_on} --how {params.how} >> {log} 2>&1
            """


####


# CART - assigning nifh clusters

# Addtional confige file for blast
configfile: "./config_cart.json"

# input_fasta_dir_cart = config["input_fasta_dir_ucyna"]
# blast_database_dir_ucyna = expand(
#     "{input_fasta_dir}{blast_database_dir}",
#     input_fasta_dir=input_fasta_dir_ucyna,
#     blast_database_dir=config["blast_database_dir_ucyna"],
# )


## assing global variable
# db1=config["database_suffix_ucyna"]
blast_output_dir=config["blast_output_dir_cart"]

rule cart:
    input:
        config[f"input_fasta_cart"],
    output:
        # config["final_blst_csv"],
        expand(
            "{blast_output_dir}{output}",
            blast_output_dir=blast_output_dir,
            output=config[f"final_blst_csv_cart"],
        )
    params:
        # blast_database_blast=expand(
        #     "{blast_database_dir}{blast_database}",
        #     blast_database_dir=blast_database_dir,
        #     blast_database=config[f"blast_database_{db1}"],
        # ),
        # output_csv=expand(
        #     "{blast_database_dir}{output_csv}",
        #     blast_database_dir=blast_database_dir,
        #     output_csv=config[f"output_csv_{db1}"],
        # ),
        # # blast_database_blast=config["blast_database_blast"],
        # blast_type=config["blast_type"],
        # num_threads=config[f"num_threads_{db1}"],
        # output_csv=config["output_csv"],
    log:
        # "logs/filter_blast{timestamp}.log",
        expand(
            # "{UCYNA_oligoreps_dir}logs/blast_{timestamp}.log",
            "{blast_output_dir}logs/{blast_database}_blast_{timestamp}.log",
            timestamp=get_timestamp(),
            blast_database=config[f"database_suffix_cart"],
            blast_output_dir=blast_output_dir,
        ),
    shell:
        """
        bash ../../databases/NifHClustersFrank2016/assignNifHclustersToNuclSeqs.sh {input} {output} > {log} 2>&1
        """

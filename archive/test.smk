# Import the datetime module
import datetime


# Function to generate a timestamp
def get_timestamp():
    return datetime.datetime.now().strftime("%Y%m%d")

configfile: "config/config_blast_db_ucyna.json"
configfile: "config/config_blast_db_arb.json"

configfile: "config/config_blast_ucyna.json"
configfile: "config/config_blast_arb.json"

# List of databases or config files
databases = ["arb", "ucyna"]

# Dictionary to store blast database directories
full_blast_database_dir_ = {}

## Define the full path for the blast database with for loop:
for database in databases:
    full_blast_database_dir = expand(
        "{input_fasta_dir}{blast_database_dir}",
        input_fasta_dir=config[f"input_fasta_dir_{database}"],
        blast_database_dir=config[f"blast_database_dir_{database}"],
    )
    full_blast_database_dir_[database] = full_blast_database_dir
    # print(full_blast_database_dir_[database])
    

## Blasting
## Blasting
## Blasting
## Blasting

# Define a rule to collect all outputs
rule all:
    input:
        expand(
            "{blast_output_dir}{output}",
            blast_output_dir=[config[f"blast_output_dir_{db}"] for db in databases],
            output=[config[f"final_blst_csv_{db}"] for db in databases]
        )


# Dynamically create rules for each database
for database in databases:
    # configfile: f"config/config_blast_{database}.json"
    blast_output_dir = config[f"blast_output_dir_{database}"]
    # print(full_blast_database_dir_[database])

    rule:
        name: f"blast_{database}"
        input:
            fasta=config[f"input_fasta_blast_{database}"]
        output:
            # csv="{blast_output_dir}final_blst.csv"
            expand(
            "{blast_output_dir}{output}",
            blast_output_dir=blast_output_dir,
            output=config[f"final_blst_csv_{database}"],
        )
        params:
            # blast_database="{blast_database_dir}{database}"
            blast_database_blast=expand(
            "{full_blast_database_dir}/{blast_database}",
            full_blast_database_dir=full_blast_database_dir_[database],
            # blast_database_dir=blast_database_dir,
            blast_database=config[f"blast_database_{database}"],
            
        ),
            output_csv=expand(
            "{full_blast_database_dir}/{output_csv}",
            full_blast_database_dir=full_blast_database_dir_[database],
            # blast_database_dir=blast_database_dir,
            output_csv=config[f"output_csv_{database}"],
        ),
            blast_type=config[f"blast_type_{database}"],
            num_threads=config[f"num_threads_{database}"]
        log:
            # "{blast_output_dir}logs/{database}_blast.log"
            expand(
            # "{UCYNA_oligoreps_dir}logs/blast_{timestamp}.log",
            "{blast_output_dir}logs/{blast_database}_blast_{timestamp}.log",
            timestamp=get_timestamp(),
            blast_database=config[f"blast_database_{database}"],
            blast_output_dir=blast_output_dir,
        ),
        shell:
            """
            bash ../scripts/blast_function.sh {input} {params.blast_database_blast} {params.blast_type} {params.num_threads} {params.output_csv} {output} > {log} 2>&1
            """


#####

# filtering blast out to keep only columns we want

configfile:"./config/config_filt_blastcols_arb.json"


rule:
    name: filter_blast_cols
    input:
        config["col_filt_blast_output"]
        # expand(
        #     "{blast_output_dir}{output}",
        #     blast_output_dir=blast_output_dir,
        #     output=config[f"cols_filt_blast_output_{database}"],
        # )
    params:
        headers_to_keep=config["headers_to_keep"]
    output:
        config["col_filt_blast_output"],
        # expand(
        #     "{blast_output_dir}{output}",
        #     blast_output_dir=blast_output_dir,
        #     output=config[f"col_filt_blast_output{database}"],
        # ),
    shell:
        """
        python3 ../scripts/filter_csv_cols_by_header.py --input_csv {input} --filtered_csv {output} --headers_to_keep {params.headers_to_keep}
        """





####
# Addtional confige file for blast
configfile: "./config/config_cart.json"

## Rule to merge csv files
configfile: "./config/config_add_auid_key_ucyna.json"
configfile: "./config/config_add_auid_key_arb.json"
configfile: "./config/config_add_auid_key_cart.json"
configfile: "./config/config_add_auid_key_genome879.json"
# configfile: "./config_merge_ucyna.json"
configfile: "./config/config_merge_cart_arb.json"
configfile: "./config/config_merge_genome879.json"
configfile: "./config/config_merge_ncd_cyano.json"
configfile: "./config/config_merge_ucyna.json"

## Merging files 


## define a merging key
merge_key = config["merge_key"]

##- rule to merge them all
rule all_merge_databases:
    input:
        # config["merged_csv"],
        expand(
            "{merged_dir}{output}",
            merged_dir=config[f"merged_dir_{merge_key[4]}"],
            output=config[f"merged_csv_{merge_key[4]}"],
        ),

## merge 1 cart with arb 
merged_dir = config[f"merged_dir_{merge_key[0]}"]
# blast_output_dir = config[f"blast_output_dir_{merge_key[0]}"]

rule: 
    name: f"merge_databases_1"
    input:
        expand(
            "{blast_output_dir}{output}",
            blast_output_dir=config[f"blast_output_dir_{merge_key[0]}"],
            output=config[f"keyed_csv_{merge_key[0]}"],
        ),
        # config["right_csv"],
    params:
        left_csv=config[f"left_csv_{merge_key[0]}"],
        # left_on=config[f"left_on_{merge_key}"],
        left_on=config[f"left_on_{merge_key[0]}"],
        right_on=config[f"right_on_{merge_key[0]}"],
        how=config[f"how_{merge_key[0]}"],
    output:
        # config["merged_csv"],
        expand(
            "{merged_dir}{output}",
            merged_dir=merged_dir,
            output=config[f"merged_csv_{merge_key[0]}"],
        ),
    log:
        expand(
            "{merged_dir}logs/merge_{timestamp}.log",
            timestamp=get_timestamp(),
            merged_dir=merged_dir
        ),
    shell:
        """
        python3 ../scripts/merge_csv.py --left_csv {params.left_csv} --right_csv {input} --merged_csv {output} --left_on {params.left_on} --right_on {params.right_on} --how {params.how} >> {log} 2>&1
        """



### merge 2 merging output with genomes879

merged_dir = config[f"merged_dir_{merge_key[2]}"]
# blast_output_dir = config[f"blast_output_dir_{merge_key[0]}"]

rule: 
    name: f"merge_databases_2"
    input:
        # config["merged_csv"],
        expand(
            "{merged_dir}{output}",
            merged_dir=config[f"merged_dir_{merge_key[0]}"],
            output=config[f"merged_csv_{merge_key[0]}"],
        ),
    params:
        left_csv=config[f"left_csv_{merge_key[2]}"],
        # left_on=config[f"left_on_{merge_key}"],
        left_on=config[f"left_on_{merge_key[2]}"],
        right_on=config[f"right_on_{merge_key[2]}"],
        how=config[f"how_{merge_key[2]}"],
    output:
        # config["merged_csv"],
        expand(
            "{merged_dir}{output}",
            merged_dir=merged_dir,
            output=config[f"merged_csv_{merge_key[2]}"],
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



### merge 3 merging output with ncd_cyanos

merged_dir = config[f"merged_dir_{merge_key[3]}"]
# blast_output_dir = config[f"blast_output_dir_{merge_key[0]}"]

rule: 
    name: f"merge_databases_3"
    input:
        # config["merged_csv"],
        expand(
            "{merged_dir}{output}",
            merged_dir=config[f"merged_dir_{merge_key[0]}"],
            output=config[f"merged_csv_{merge_key[2]}"],
        ),
    params:
        left_csv=config[f"left_csv_{merge_key[3]}"],
        # left_on=config[f"left_on_{merge_key}"],
        left_on=config[f"left_on_{merge_key[3]}"],
        right_on=config[f"right_on_{merge_key[3]}"],
        how=config[f"how_{merge_key[3]}"],
    output:
        # config["merged_csv"],
        expand(
            "{merged_dir}{output}",
            merged_dir=merged_dir,
            output=config[f"merged_csv_{merge_key[3]}"],
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



### merge 4 merge output with ucyna_oligos

merged_dir = config[f"merged_dir_{merge_key[4]}"]
# blast_output_dir = config[f"blast_output_dir_{merge_key[0]}"]

rule: 
    name: f"merge_databases_4"
    input:
        # config["merged_csv"],
        expand(
            "{merged_dir}{output}",
            merged_dir=config[f"merged_dir_{merge_key[0]}"],
            output=config[f"merged_csv_{merge_key[3]}"],
        ),
    params:
        left_csv=config[f"left_csv_{merge_key[4]}"],
        # left_on=config[f"left_on_{merge_key}"],
        left_on=config[f"left_on_{merge_key[4]}"],
        right_on=config[f"right_on_{merge_key[4]}"],
        how=config[f"how_{merge_key[4]}"],
    output:
        # config["merged_csv"],
        expand(
            "{merged_dir}{output}",
            merged_dir=merged_dir,
            output=config[f"merged_csv_{merge_key[4]}"],
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




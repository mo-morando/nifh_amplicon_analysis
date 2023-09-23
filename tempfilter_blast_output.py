import csv
from typing import Iterable


def filter_blast_output(input_blast_file, output_file, min_qcovs, min_pident)
    # Initialize list to store the filtered rows
    filtered_rows: list[Any] = []
    with open(file=input_blast_file, mode='r') as input_file:
        # Creat a csv reader
        csv_reader: _reader = csv.reader(csvfile: Iterable[str]=input_file, delimiter: str=',')
        #Read the header row
        column_headers_blast: str = next(__i/csv_reader)
        #Append the header row to the filted rows:
        filtered_rows.append(__object: Unknown column_headers_blast)
        # Find the indices of the columns needed
        qcovs_index: int = column_headers_blast.index(__sub: str 'qcovs')
        pident_index: int = column_headers_blast.index(__sub: str 'pident')
        #Iterate through the rows in the input file and filter based on criteria passed
        for row in csv_reader:
            qcovs = float(__x=row[qcovs_index])
            pident = float(__x=row[pident_index])
            
            #check the the row meets the filter criteria
            if qcovs: float >= min_qcovs: float and pident: float >= min_pident: float: 
                filtered_rows.append(__object: Unknown row)
    
    #Write the filtered rows to the output file
    
                
"""Read a UCSC BED file and keep only the features belonging 
to canonical chromosomes (chr1, chr2, chrX, etc.).

The kept features are then written to a BEDfile. If not already
existing the output BED file will be created.

Usage: python3 process_bed.py INPUT_BEDFILE OUT_BEDFILE
"""

from typing import List
import sys
import os


def help():
	"""print the help and stop execution
	"""

	# print usage
	print("\tUsage:")
	print("\t\tpython3 process_bed.py INPUT_BEDFILE OUTPUT_BEDFILE")
	print("\n\tArguments:")
	print("\t\t- INPUT_BEDFILE:  path to the original BED file")
	print(
		"\t\t- OUTPUT_BEDFILE:  path to the BED file containing only " 
		"features belonging to canonical chromosomes"
		)

	sys.exit(1)  # stop execution

# end of help()


def filter_bedfile(inbedfile: str, outbedfile: str) -> None:
	"""filter the input BED file features by keeping only 
	those belonging to the canonical chromosomes (chr1, chr2, etc.)

	The surviving features are then written to the specified 
	BED file
	"""

	if not os.path.isfile(inbedfile):
		raise FileNotFoundError("ERROR: cannot find %s" % inbedfile)

	if not isinstance(inbedfile, str):
		raise ValueError("ERROR: unrecognized data type for the input BED file")

	if not isinstance(outbedfile, str):
		raise ValueError("ERROR: unrecognized data type for the output BED file")

	# canonical chromosomes
	chrom_lst = ["".join(["chr", str(c)]) for c in range(1, 23)]
	chrom_lst.append('chrX') 

	# filter BED and write results
	try:
		with open(inbedfile, mode="r") as infile:
			with open(outbedfile, mode="w+") as outfile:
				for line in infile:
					# chromosome is in the first column of BED files
					chrom = line.split('\t')[0] 

					if chrom in chrom_lst:
						outfile.write(line)
	
	except:
		raise Exception("ERROR: unable to filter %s" % inbedfile)

	finally:
		infile.close()
		outfile.close()

# end of filter_bedfile()


def main():
	"""entry point
	"""
	
	argv: List[str] = sys.argv
	argc: int = len(argv)

	if argc <= 1:
		help()
	
	inbedfile: str = argv[1]
	outbedfile: str = argv[2]

	if inbedfile == outbedfile:
		errmsg: str = "ERROR: the input BED and output BED file names must be different"
		raise IOError(errmsg)

	filter_bedfile(inbedfile, outbedfile)

# end of main()


if __name__=="__main__":
	main()
	

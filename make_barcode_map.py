""""""

def main(config, threshold):
	
	#open config file, create SeqLib object
	c = json.load(open(config, "U"))
	lib = seqlib.overlap.OverlapSeqLib(c["libraries"][0])
	
	conflicts = open(lib.forward.strip("fq")+"conflicts", "w")
	map_file = open(lib.forward.strip("fq")+"mapping", "w")
	
	#dictionary to store all barcodes
	barcode_dict = {}
	
	#read FASTQ files 
	for fwd, rev, bc in read_fastq_multi([lib.forward, lib.reverse, c["libraries"][0]["fastq"]["barcode"]]):

		#fuse fwd and rev
		fused = lib.fuse_reads(fwd, rev)
		
		#store data
		if fused != None:
			if bc.sequence in barcode_dict:
				barcode_dict[bc.sequence].append(fused)
			else:
				barcode_dict[bc.sequence] = [fused]

	#for each barcode, define consensus sequence
	for b in barcode_dict:
		
		#how many reads had this barcode?
		len_b = float(len(barcode_dict[b]))
		
		#count the unique sequences of that barcode
		count_b = Counter([x.sequence for x in barcode_dict[b]])

		#if there's only one variant, then we've got our consensus
		if len(count_b) == 1:
			print >> map_file, "\t".join([b, barcode_dict[b][0].sequence])
		
		#otherwise see if there's a consensus sequence
		else:
			conc_b = count_b.most_common(10)
			
			#if there's a conflict (based on a threshold), then print to conflicts file
			if conc_b[1][1]/len_b > threshold:
				#print all the sequences that are above the conflict threshold
				all_conflicts = []
				for seq in conc_b:
					if seq[1]/len_b > threshold: all_conflicts.append(seq[0]+":"+str(seq[1]))
				print >> conflicts, "\t".join([b, ",".join(all_conflicts)])
			
			#if there isn't a conflict, print to the mapping file
			else:
				print >> map_file, "\t".join([b, conc_b[0][0]])
	
	conflicts.close()
	map_file.close()
		
				
if __name__ == "__main__":
	import seqlib.overlap
	import json
	from fastq_util import read_fastq_multi
	from optparse import OptionParser
	from collections import *
			
	parser = OptionParser()
	parser.add_option('--config', action='store', type = 'string', dest = 'config', help = 'json config file (as for Enrich)')
	parser.add_option('--threshold', action = 'store', type = 'float', dest = 'threshold', help = 'frequency of next highest sequence for a given barcode to define conflict', default = .1)
	(option, args) = parser.parse_args()
	
	main(option.config, option.threshold)

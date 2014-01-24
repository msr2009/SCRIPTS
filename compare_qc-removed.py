#/usr/bin/python

"""
compare_qc-removed.py

Script takes input of two qc_removed files, finds either intersection or difference of reads between them, 
then outputs a list of readIDs for those reads

python compare_qc_removed.py --f1 FILE1 --f2 FILE2 --mode MODE
"""


def main():
	
	f1_set = set()
	f2_set = set()
	
	#create sets of readIDs for each file
	for line in open(option.infile1, 'rU'):
		l = line.strip().split('\t')
		if 'local' in l[1]:
			f1_set.add(l[0])
	
	for line in open(option.infile2, 'rU'):
		l = line.strip().split('\t')
		if 'local' in l[1]:
			f2_set.add(l[0])
		
	#take difference and print 
	if option.intersection == 'difference':
		f_out_name = '/'.join(option.infile1.split('/')[0:-1]) +'/'+ option.infile1.split('/')[-1].split('_index_')[0] + '-' + option.infile2.split('/')[-1].split('_index_')[0] + '-qc_removed_difference'
		f_out = open(f_out_name, 'w')
		for read in f1_set.symmetric_difference(f2_set):
			print >> f_out, '\t'.join([read, '1'])
#		for read in f2_set.difference(f1_set):
#			print >> f_out, '\t'.join([read, '2'])	
		f_out.close()
			
	elif option.intersection == 'intersection':
		f_out_name = '/'.join(option.infile1.split('/')[0:-1]) + '/' + option.infile1.split('/')[-1].split('_index_')[0] + '-' + option.infile2.split('/')[-1].split('_index_')[0] + '-qc_removed_intersection'
		f_out = open(f_out_name, 'w')
		for read in f1_set.intersection(f2_set):
			print >> f_out, '\t'.join([read, '1'])
		f_out.close()

if __name__ == '__main__':
	from optparse import OptionParser
	import sys, time

	parser = OptionParser()
	parser.add_option('--f1', action = 'store', type = 'string', dest = 'infile1', help = 'path to qc_removed file')
	parser.add_option('--f2', action = 'store', type = 'string', dest = 'infile2', help = 'path to qc_removed file')
	parser.add_option('--mode', action = 'store', type = 'string', dest = 'intersection', help = 'mode for comparison: intersection,difference,both')
	(option, args) = parser.parse_args()
	
	print time.asctime()
	main()
	print time.asctime()
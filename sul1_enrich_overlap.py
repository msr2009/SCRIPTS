{
  "libraries" : [
    {
      "name" : "20131105_sul1p",
      "description" : "SUL1p library -- 500cycle miseq run",
      "wild type" : {
        "sequence" : "TAACGAACAAGTTGTCAAATTAGACCCATAATAATTTTGAACACTTCTACCTGTTCATGTCTTTTCTCGAACACTGTCATTTGAAATTATGCACTGTGAAAAAAAGAAACAAAGACCAAAAGAATAATATAAATAGTGAAGTAAAATGTGTTGTAATGCACATGGATCTTGTACTGCTCAAACTTAATATTTCTATTGTAGAAAAATTTTCGATTTAAAATTGTGAAACCGATTATATAAAAGTATATTAGCTGACATTAACGTCTCAAAACCAGGTCAATAGCTTTAAAAATAAAAATAAAT",
        "coding" : false
      },
      "fastq" : {
        "forward" : "/net/fields/vol1/home/mattrich/DUNHAM/SUL1/20131105/20131105_sul1p12N_1_reads.q30.fq",
        "reverse" : "/net/fields/vol1/home/mattrich/DUNHAM/SUL1/20131105/20131105_sul1p12N_2.q30.fq"
      },
      "filters" : {
        "remove unresolvable" : true,
        "chastity" : true,
        "max mutations" : 10,
        "min quality" : 15
      },
      "overlap" : {
        "forward start" : 54,
        "reverse start" : 66,
        "length" : 185,
        "overlap only" : false,
        "max mismatches" : 10
      }
    }
  ]
}

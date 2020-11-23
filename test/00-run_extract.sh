clear
# python3.7 ../01-extract_geno.py -v rs.file -o /data/programs/scripts/UKBB/test_extract_geno/out1/out1 -m imp -e exclude.list $1
# python3.7 ../01-extract_geno.py -v rs.file -o /data/programs/scripts/UKBB/test_extract_geno/out2/out2 -m geno -e exclude.list $1
python3.7 ../01-extract_geno.py -v rs.file -o out3 -m wes $1


 #  487352 out1/out1
 #  488320 out2/out2
 #   49903 out3/out3
 # 1025575 total

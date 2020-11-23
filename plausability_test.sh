


# head -n 1 /data/studies/06_UKBB/UKBB_150k/01_Raw_Data/ukb8974.csv | sed 's/"//g' | tr "," "\n" | grep -n "21000-0.0"


#### Ethnic background
# cat /data/studies/06_UKBB/UKBB_150k/01_Raw_Data/ukb8974.csv | cut -d "," -f 1095 | sed 's/"//g' | awk '{ if ($1 == 1001) print }' | wc -l


# cat /data/studies/06_UKBB/UKBB_150k/01_Raw_Data/ukb8974.csv | cut -d "," -f 1111 > test.txt
# sed -i 's/"//g' test.txt

cat test.txt | sort -u

echo "-------------------"

cauc=$(cat test.txt | awk '{if ($1 == 1) print}' | wc -l)

nonCauc=$(cat test.txt | awk '{if ($1 == "") print}' | wc -l)

echo "Caucasian: $cauc"
echo "Non-Caucasian: $nonCauc"

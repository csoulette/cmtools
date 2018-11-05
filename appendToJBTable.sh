FILE1=$1
FILE2=$2


if [ -z "$FILE1" ]
then
   echo "Missing FILE1 argument. Correct usage syntax: bash appendToJBTable.sh FILE1 (drim_out.tsv) FILE2 (amended_jb_Table.tsv)"
   exit 1
fi

if [ -z "$FILE2" ]
then
   echo "Missing FILE2 argument. Correct usage syntax: bash appendToJBTable.sh FILE1 (drim_out.tsv) FILE2 (amended_jb_Table.tsv)"
   exit 1
fi


awk 'FNR==NR{a[$2]=$(NF-1)"\t"$NF; next} {if ($1 in a) { print $0,a[$1] } else if ($1=="gene_id") { print }}' $FILE1 $FILE2

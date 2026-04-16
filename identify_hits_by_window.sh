

### identify_hits_by_window_Mar6.sh:

WINLIST="CnanuMar6_90pc_csv"   # your CSV list
TSV="C1_13_Feb26_input.tsv"

while IFS=, read -r chr s e; do
  awk -v chr="$chr" -v s="$s" -v e="$e" '
    BEGIN{OFS="\t"}
    {
      # reorder to match your one-off output exactly
      print $5, $3, $1, $4, $2
    }
  ' "$TSV" \
  | awk -v chr="$chr" -v s="$s" -v e="$e" '
      $1==chr && $2>=s && $2<=e
    ' \
  | grep "Cnanu"
done < "$WINLIST"

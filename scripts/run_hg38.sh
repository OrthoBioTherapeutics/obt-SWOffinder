mm=4
w=4
b=2

ref_path='../misc/references/hg38.fa'
sgrna_list_path='./guides_human.txt'
output_path="./out/human/${mm}mm${w}w${b}b"
maxE=$mm
maxM=$mm
maxMB=3
maxB=$b
threads=4
best_window=true
best_window_size=$w
PAM='NGG'
PAM_edits=true

java -cp bin \
SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign \
$ref_path \
$sgrna_list_path \
$output_path \
$maxE \
$maxM \
$maxMB \
$maxB \
$threads \
$best_window \
$best_window_size \
$PAM \
$PAM_edits

#2. **sgRNA list**: The path of a text file containing a list of sgRNAs with their PAM (see sgRNAs.txt file for example).
#4. **maxE**: Max edits allowed (integer).
#5. **maxM**: Max mismatches allowed without bulges (integer).
#6. **maxMB**: Max mismatches allowed with bulges (integer).
#7. **maxB**: Max bulges allowed (integer).
#9. **Best in-a-window**: Flag whether to choose the best off-target site in a window or not (true or false).
#10. **Best in-a-window size**: The window size for choosing the best in a window (integer). Please insert even if **Best in-a-window** is false.
#12. **Allow PAM edits**: Flag whether to allow PAM edits or not (true or false).


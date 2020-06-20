for ii in SFig/Fig*.pdf;do 
  base=`basename $ii`
  lab=`echo $base|sed 's/Fig.\?_\(S[0-9]\+\).*/Supplemental Figure \1/'`
  echo $base $lab
  bash resizeTab.bash $ii adjusted/$base 1
  bash addHeader.bash "$lab" adjusted/$base labeled/$base
done
pdftk labeled/Fig._S*.pdf output out/SupplementalFigures.pdf


for ii in fig/Fig*.pdf;do 
  base=`basename $ii`
  lab=`echo $base|sed 's/Fig.\?_\([0-9]\+\).*/Figure \1/'`
  echo $base $lab
  bash resizeTab.bash $ii adjusted/$base 1
  bash addHeader.bash "$lab" adjusted/$base labeled/$base
done
pdftk labeled/Fig._[0-9]*.pdf output out/Figures.pdf


for ii in tab/Table*.pdf;do 
  base=`basename $ii`
  lab=`echo $base|sed 's/Table_\(S[0-9]\).*/Table \1/'`
  bash resizeTab.bash $ii adjusted/$base
  bash addHeader.bash "$lab" adjusted/`basename $ii` labeled/$base 0
done
pdftk labeled/Table*.pdf output out/SupplementaryTables.pdf



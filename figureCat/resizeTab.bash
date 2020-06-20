targetWidth=8.5
#scale=`Rscript -e "message(round($targetWidth/$width,5))"
pdfcrop $1 tmp.pdf
width=$(bash pdfScale.sh -i tmp.pdf|grep Inch|sed 's/[^0-9]\+\([0-9.]\+\) x.*/\1/')
height=$(bash pdfScale.sh -i tmp.pdf|grep Inch|sed 's/[^0-9]\+\([0-9.]\+\) x \([0-9.]\+\)/\2/')
scale=`echo "scale=5;$targetWidth/$width"|bc --mathlib`
newHeight=`echo "scale=5;(${height}*$scale)+1" |bc --mathlib`
newWidth=`echo "scale=5;${width}*$scale" |bc --mathlib`
if $4;then
  targetWidth=$width
  newHeight=$(echo $height+1.0|bc --mathlib)
fi
echo $width $height $scale $newWidth $newHeight
bash pdfScale.sh -v -r "custom in $targetWidth $newHeight" --vert-align top -a none tmp.pdf $2
#pdfcrop tmp2.pdf $2
bash pdfScale.sh -i $2


#echo ""|enscript -L1 --header="||$1" --output - |ps2pdf ->tmp.pdf
#pdfcrop tmp.pdf tmp2.pdf
#width=$(bash pdfScale.sh -i tmp2.pdf|grep Inch|sed 's/[^0-9]\+\([0-9.]\+\) x.*/\1/')
width=$(bash pdfScale.sh -i $2|grep Inch|sed 's/[^0-9]\+\([0-9.]\+\) x.*/\1/')
height=$(bash pdfScale.sh -i $2|grep Inch|sed 's/[^0-9]\+\([0-9.]\+\) x \([0-9.]\+\)/\2/')
echo "\documentclass{article}
\usepackage{helvet}
\renewcommand*\familydefault{\sfdefault}
\usepackage[paperheight=${height}in,paperwidth=${width}in,margin=0.1in]{geometry}
\begin{document}
\noindent
\textbf{\Huge $1}
\end{document}
">tmp.tex
pdflatex tmp.tex 1>/dev/null
pdftk $2 stamp tmp.pdf output $3


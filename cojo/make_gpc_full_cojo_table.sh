gpc=2

cat header.txt > gpc"$gpc".txt
tail -n +2 -q gpc"$gpc"/*jma* >> gpc"$gpc".txt

chmod -R 777 .

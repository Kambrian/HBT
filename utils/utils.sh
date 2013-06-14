#to get a list of snapshot numbers
ls groups_* -d|grep -o [0-9]+ -E

for i in {0..1023};
 do 
 name=snapdir_`printf %03d $i`; 
 gname=groups_`printf %03d $i`;
 if [ -d $name -a -d $gname ]; 
 then 
 n=`ls $name |wc -l`; 
 m=`ls $gname |wc -l`;
 if [ $n -gt 0 -a $m -gt 0 ]; then
 echo $i $n $m; 
 fi
 fi; 
done

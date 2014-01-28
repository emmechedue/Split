cat ./1/ensambleN.txt >> ensambleN.txt
cat ./2/ensambleN.txt >> ensambleN.txt
cat ./3/ensambleN.txt >> ensambleN.txt
cat ./4/ensambleN.txt >> ensambleN.txt
cat ./5/ensambleN.txt >> ensambleN.txt
cat ./6/ensambleN.txt >> ensambleN.txt

cat ./1/ensamblex.txt >> ensamblex.txt
cat ./2/ensamblex.txt >> ensamblex.txt
cat ./3/ensamblex.txt >> ensamblex.txt
cat ./4/ensamblex.txt >> ensamblex.txt
cat ./5/ensamblex.txt >> ensamblex.txt
cat ./6/ensamblex.txt >> ensamblex.txt

cp ./1/time.txt ./
sed "s|N_loop= 50|N_loop= 300|" ./1/parameters.txt >./parameters.txt

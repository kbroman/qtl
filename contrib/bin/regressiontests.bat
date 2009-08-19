echo "Regressiontests for standalone MQM"
echo "Dataset 1 F2"
mqm.exe -v -pTest/std/phenotypes1.txt -gTest/std/genotypes1.txt -mTest/std/markers1.txt -sTest/t11/settings.txt > test/t11out.txt
mqm.exe -v -pTest/std/phenotypes1.txt -gTest/std/genotypes1.txt -mTest/std/markers1.txt -sTest/t12/settings.txt -cTest/t12/cofactors.txt  > test/t12out.txt

echo "Dataset 2 BC=~ Hyper (genotypes NA -> 0)"
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2.txt -mTest/std/markers2.txt -sTest/t21/settings.txt > test/t21out.txt
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2.txt -mTest/std/markers2.txt -sTest/t22/settings.txt -cTest/t22/cofactors.txt > test/t22out.txt
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2.txt -mTest/std/markers2.txt -sTest/t23/settings.txt -cTest/t23/cofactors.txt > test/t23out.txt
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2m.txt -mTest/std/markers2.txt -sTest/t24/settings.txt > test/t24out.txt
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2m.txt -mTest/std/markers2.txt -sTest/t25/settings.txt > test/t25out.txt

echo "Dataset 3 F2=~ Listeria"
mqm.exe -v -pTest/std/phenotypes3.txt -gTest/std/genotypes3.txt -mTest/std/markers3.txt -sTest/t31/settings.txt > test/t31out.txt
mqm.exe -v -pTest/std/phenotypes3.txt -gTest/std/genotypes3.txt -mTest/std/markers3m.txt -sTest/t32/settings.txt > test/t32out.txt
mqm.exe -v -pTest/std/phenotypes3.txt -gTest/std/genotypes3.txt -mTest/std/markers3.txt -sTest/t33/settings.txt -cTest/t33/cofactors.txt > test/t33out.txt

echo "Output done"
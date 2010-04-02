echo "Regressiontests for standalone MQM"
echo "Dataset 1 F2"
mqm.exe -v -pTest/std/phenotypes1.txt -gTest/std/genotypes1.txt -mTest/std/markers1.txt -sTest/std/settings1.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --neglect=1 > test/t11out.txt
mqm.exe -v -pTest/std/phenotypes1.txt -gTest/std/genotypes1.txt -mTest/std/markers1.txt -sTest/std/settings1.txt -cTest/t12/cofactors.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --neglect=1 > test/t12out.txt

echo "Dataset 2 BC=~ Hyper (genotypes NA -> 0)"
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2.txt -mTest/std/markers2.txt -sTest/std/settings2.txt  --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --neglect=1 > test/t21out.txt
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2.txt -mTest/std/markers2.txt -sTest/std/settings2.txt -cTest/t22/cofactors.txt  --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --neglect=1 > test/t22out.txt
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2.txt -mTest/std/markers2.txt -sTest/std/settings2.txt -cTest/t23/cofactors.txt  --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --neglect=1 > test/t23out.txt
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2m.txt -mTest/std/markers2.txt -sTest/std/settings2.txt  --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=y --maugment=10000 --miaugment=250 --neglect=1 > test/t24out.txt
mqm.exe -v -pTest/std/phenotypes2.txt -gTest/std/genotypes2m.txt -mTest/std/markers2.txt -sTest/std/settings2.txt  --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=y --maugment=10000 --miaugment=250 --neglect=1 > test/t25out.txt

echo "Dataset 3 F2=~ Listeria"
mqm.exe -v -pTest/std/phenotypes3.txt -gTest/std/genotypes3.txt -mTest/std/markers3.txt -sTest/std/settings3.txt  --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --neglect=1> test/t31out.txt
mqm.exe -v -pTest/std/phenotypes3.txt -gTest/std/genotypes3.txt -mTest/std/markers3m.txt -sTest/std/settings3.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=y --maugment=10000 --miaugment=250 --neglect=1 > test/t32out.txt
mqm.exe -v -pTest/std/phenotypes3.txt -gTest/std/genotypes3.txt -mTest/std/markers3.txt -sTest/std/settings3.txt -cTest/t33/cofactors.txt  --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --neglect=1 > test/t33out.txt

echo "Output done"
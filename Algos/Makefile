exp:
	export PATH="/scratch/shared/cs6219/dependencies/GCC/bin:/scratch/shared/cs6219/dependencies/squashfs-root/usr/bin:$PATH"
	export LD_LIBRARY_PATH=/scratch/shared/cs6219/dependencies/GCC/lib64/:$LD_LIBRARY_PATH

run_our_bma: exp
	python3.6 BMA-Our/recon.py --i ../Data/cov5/nanopore-sim-cond-ld/clusters.txt  --o  ../Data/cov5/nanopore-sim-cond-ld/recons-our-bma.txt --path /dev/shm --L 110

build_sab_bma: exp
	g++ -std=c++11 BMA-Sabary/BMALookaheadMod.cpp -o BMA

run_sab_bma: build_sab_bma
	./BMA ../Data/cov5/nanopore-sim-cond-ld/clusters.txt ../Data/cov5/nanopore-sim-cond-ld/recons-sab-bma.txt

build_iterative: exp
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o LCS2.o Iterative/LCS2.cpp
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o EditDistance.o Iterative/EditDistance.cpp
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Clone.o Iterative/Clone.cpp
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o Cluster2.o Iterative/Cluster2.cpp
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o LongestPath.o Iterative/LongestPath.cpp
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o CommonSubstring2.o Iterative/CommonSubstring2.cpp
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o DividerBMA.o Iterative/DividerBMA.cpp
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o DNA.o Iterative/DNA.cpp
	g++ -o DNA *.o

b_iterative: exp
	g++ -std=c++0x -O3 -g3 -Wall -c -fmessage-length=0 -o DNA.o Iterative/DNA.cpp
	g++ -o DNA *.o

run_iterative: 
	./DNA ${clusters} ${refs} ${strandLen} ${output}

compare:
	wc -l $(refs) $(recons)
	python3 compare.py $(refs) $(recons) ${strandLen}
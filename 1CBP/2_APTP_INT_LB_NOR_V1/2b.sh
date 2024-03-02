#!/bin/bash

for i in 01 02 03
do
	./2_APTP_BIN_LB/MODEL "./_INPUT/CGCUT/" "CGCUT${i}.TXT" "./_OUTPUT/2_APTP_BIN_LB/Gen.txt" > "./_OUTPUT/2_APTP_BIN_LB/garbage.txt"
done

for i in 01 02 03 04 05 06 07 08 09
do
	./2_APTP_BIN_LB/MODEL "./_INPUT/HT/" "HT${i}.TXT" "./_OUTPUT/2_APTP_BIN_LB/Gen.txt" > "./_OUTPUT/2_APTP_BIN_LB/garbage.txt"
done

for i in 01 02 03 04 05 06 07 08 09 10
do
	./2_APTP_BIN_LB/MODEL "./_INPUT/BENG/" "BENG${i}.TXT" "./_OUTPUT/2_APTP_BIN_LB/Gen.txt" > "./_OUTPUT/2_APTP_BIN_LB/garbage.txt"
done

for i in 01 02 03 04 05 06 07 08 09 10 11 12
do
	./2_APTP_BIN_LB/MODEL "./_INPUT/NGCUT/" "NGCUT${i}.TXT" "./_OUTPUT/2_APTP_BIN_LB/Gen.txt" > "./_OUTPUT/2_APTP_BIN_LB/garbage.txt"
done

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13
do
	./2_APTP_BIN_LB/MODEL "./_INPUT/BKW/" "BKW${i}.TXT" "./_OUTPUT/2_APTP_BIN_LB/Gen.txt" > "./_OUTPUT/2_APTP_BIN_LB/garbage.txt"
done

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13
do
	./2_APTP_BIN_LB/MODEL "./_INPUT/GCUT/" "GCUT${i}.TXT" "./_OUTPUT/2_APTP_BIN_LB/Gen.txt" > "./_OUTPUT/2_APTP_BIN_LB/garbage.txt"
done

for i in 01 02 03 05 07 09 10
do
	for j in 020 040 060 080 100
	do
		for k in 01 02 03 04 05 06 07 08 09 10
		do
			./2_APTP_BIN_LB/MODEL "./_INPUT/CLASS/" "cl_${i}_${j}_${k}.TXT" "./_OUTPUT/2_APTP_BIN_LB/Gen.txt" > "./_OUTPUT/2_APTP_BIN_LB/garbage.txt"
		done
	done
done



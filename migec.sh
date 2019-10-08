cd /home/users/kjyi/Projects/cocult/cocult2
cat sample.txt | while read num b1 b2 b3 b4 am as bm bs r1 r2 name; do
	cat << EOF | qsub -q day -l nodes=1:ppn=1 -e /dev/null -o /dev/null	
	cd /home/users/kjyi/Projects/cocult/cocult2
	mkdir -p migec/$name
	sh /home/users/kjyi/Projects/cocult/pipe.sh \
		-o migec/$name \
		--a1 $am \
		--a2 $as \
		--b1 $bm \
		--b2 $bs \
		--f1 $r1 \
		--f2 $r2 &> migec/$name/migec_pipeline.log
EOF
done

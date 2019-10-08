
cat << EOF 1>&2
--------------
migec pipeline
--------------

`date`

args: $@


EOF

. ~kjyi/src/parse2 $@ << EOF

-o		odir	"."	output directory
--a1	a1		""	master barcode sequence for TCRA
--a2	a2		""	slave barcode sequence for TCRA
--b1	b1		""	master barcode sequence for TCRB
--b2	b2		""	slave barcode sequence for TCRB
--f1	f1		""	input fastq 1
--f2	f2		""	input fastq 2

EOF


if [ -f $f1 ] && [ -f $f2 ] ; then
	echo Input files present: $f1, $f2
else
	echo Input files are not exist;	exit 1
fi

mkdir -p $odir
cat << EOF > $odir/barcodes.txt
TCRB	$b1	$b2
TCRA	$a1	$a2
EOF

echo "[[ Checkout ]]"
mkdir -p $odir/checkout
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Checkout -cute $odir/barcodes.txt $f1 $f2 $odir/checkout

echo "[[ Histogram ]]"
mkdir -p $odir/histogram
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Histogram $odir/checkout $odir/histogram

pushd $odir/histogram
R --slave << EOF

#
# Plots a fancy oversequencing histogram from tables pre-computed by Histogram util
#
# usage:
# \$cd histogram_output_dir/
# \$RScript histogram.R
#

require(ggplot2); require(reshape)

build_df <- function(stat, units, sweep_df = NULL) {
	df<- read.table(paste(stat, units, ".txt", sep= ""), header=FALSE, comment="", sep = "\\t")
  x<-df[1,3:ncol(df)]
	id<-df[2:nrow(df), 1] 
	if (is.null(sweep_df)) {
		sweep_df <- df
	}
	df[2:nrow(df),3:ncol(df)] <- sweep(df[2:nrow(df),3:ncol(df)], 1, rowSums(sweep_df[2:nrow(df),3:ncol(df)]), FUN = "/")
	      
	df.m <- melt(df[2:nrow(df), c(1,3:ncol(df))])
	df.m\$x <- as.vector(t(x[df.m\$variable]))
	df.m\$s <- rep(stat, nrow(df.m))
	list(sweep_df, df.m)						   
}
   
for (units in c("", "-units")) {  
   dfo <- build_df("overseq", units)
   dfc <- build_df("collision1", units, dfo[[1]])     
   df <- rbind(dfo[[2]], dfc[[2]])
   
   pdf(paste("histogram", units, ".pdf", sep= ""))
   
	 print(
		ggplot(df, aes(x=x,y=value))+geom_smooth()+
			scale_x_log10(name = "MIG size, reads", expand=c(0,0), limits=c(1, 10000), breaks = c(1:10,100,1000,10000), labels = c("1", rep("", 8), 10, 100, 1000, 10000), oob=scales::rescale_none)+
			scale_y_continuous(name = "",expand=c(0,0), limits = c(0, max(df\$value)), oob=scales::rescale_none)+theme_bw()+theme(panel.grid.minor = element_blank()) + facet_grid(s~.)
	 )
   
   dev.off()  
}
EOF

R --slave << EOF

#
# Plots position-weight matrices for UMI sequences based on data precomputed by Histogram util
#
# usage:
# \$cd histogram_output_dir/
# \$RScript pwm.R
#

require(seqLogo) # available @ bioconductor

# read in
logo <- function(prefix) {
	df<-read.table(paste(prefix,".txt",sep=""))
	rownames(df)<-df[,1]
	df<-df[,2:ncol(df)]
	df[, ] <- apply(df[, ], 2, as.numeric)
	
	# build seqlogo
	df.p <- makePWM(df)
	pdf(paste(prefix,".pdf",sep=""))
	seqLogo(df.p, ic.scale = F)
	dev.off()
	
	con <- file(paste(prefix,"-stats.txt",sep=""))
	sink(con, append=TRUE)
	sink(con, append=TRUE, type="output")
	
	# stats
	print(paste("IC =",sum(df.p@ic)))
	
	#entropy
	h<-sum(2-df.p@ic)
	print(paste("H =",h))
	
	#correlation with pos
	print(paste("R(Hi,i) =", cor(-df.p@ic,1:length(df.p@ic))))
	
	#number of variants
	print(paste("N_obs =",2^h))
	
	#theoretical number of variants
	print(paste("N_exp =",2^(2*ncol(df))))
	sink()
}

logo("pwm-summary")
logo("pwm-summary-units")
EOF

popd

echo "[[ Assemble ]]"
mkdir -p $odir/assemble
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar AssembleBatch --force-collision-filter --force-overseq 5 $odir/checkout $odir/histogram $odir/assemble

echo "[[ CdrBlast A ]]"
mkdir -p $odir/cdr_A
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar CdrBlastBatch -R TRA $odir/checkout $odir/assemble $odir/cdr_A

echo "[[ CdrBlast B ]]"
mkdir -p $odir/cdr_B
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar CdrBlastBatch -R TRB $odir/checkout $odir/assemble $odir/cdr_B

echo "[[ filter TRA ]]"
mkdir -p $odir/cdr_A_filter
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar FilterCdrBlastResultsBatch $odir/cdr_A $odir/cdr_A_filter

echo "[[ filter TRB ]]"
mkdir -p $odir/cdr_B_filter
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar FilterCdrBlastResultsBatch $odir/cdr_B $odir/cdr_B_filter

echo "[[ reporting ]]"
mkdir -p $odir/report/A
mkdir -p $odir/report/B
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Report -c $odir/checkout -i $odir/histogram -a $odir/assemble -b $odir/cdr_A -f $odir/cdr_A_filter $odir/report/A
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Report -c $odir/checkout -i $odir/histogram -a $odir/assemble -b $odir/cdr_B -f $odir/cdr_B_filter $odir/report/B
sed -i '1{s,^,Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc");,}' $odir/report/A/*.R
sed -i '1{s,^,Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc");,}' $odir/report/B/*.R
sed -i '2{s!)!, knit_root_dir = "'$(pwd)'")!}' $odir/report/A/*.R
sed -i '2{s!)!, knit_root_dir = "'$(pwd)'")!}' $odir/report/B/*.R
sed -i '200,230{s,df$READS_DROPPED_WITHIN_MIG / df$READS_TOTAL,(df$READS_DROPPED_WITHIN_MIG_1 + df$READS_DROPPED_WITHIN_MIG_2) / df$READS_TOTAL,}' $odir/report/A/*.Rmd
sed -i '200,230{s,df$READS_DROPPED_WITHIN_MIG / df$READS_TOTAL,(df$READS_DROPPED_WITHIN_MIG_1 + df$READS_DROPPED_WITHIN_MIG_2) / df$READS_TOTAL,}' $odir/report/B/*.Rmd
sed -i '168{s!levels=\(.*\)!levels=unique(\1)!}' $odir/report/A/*.Rmd
sed -i '168{s!levels=\(.*\)!levels=unique(\1)!}' $odir/report/B/*.Rmd
Rscript $odir/report/A/*.R
Rscript $odir/report/B/*.R

echo "[[ done ]]"




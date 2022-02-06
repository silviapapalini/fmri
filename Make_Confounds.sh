subjects="sub-FG01 sub-FG02"
task="Extinction"

for subj in $subjects; do
	prefix="derivatives/afni/$subj"
	all_confounds="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_desc-confounds_timeseries.tsv"

	Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
n <- names(X)

#CONFOUNDS 1
S <- startsWith(n, "trans") | startsWith(n, "rot") | startsWith(n, "white_matter") | startsWith(n, "csf") | startsWith(n, "motion_outlier")


#CONFOUNDS 2
#S <- startsWith(n, "trans") | startsWith(n, "rot") | n == "white_matter" | n == "csf" | startsWith(n, "motion_outlier")
X <- lapply(X[S], function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
# tried on FG01: these confounds produced activations outside the brain and in the crf,(visible from the f-stat contrast)


write.table(X, file=args[2], row.names=F, col.names=F, sep=" ")
' $all_confounds $prefix/confounds
done

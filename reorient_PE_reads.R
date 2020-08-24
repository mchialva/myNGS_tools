#!/usr/bin/env Rscript
start.time <- proc.time()

args <- commandArgs(TRUE)

R1_path<-args[[1]]
R2_path<-args[[2]]
fw_primer<-args[[3]]
rev_primer<-args[[4]]
primer_mismatch_fw<-as.integer(args[[5]])
primer_mismatch_rev<-as.integer(args[[6]])
output_dir<-args[[7]]
barcode_length<-as.integer(args[[8]])
nthreads<-as.integer(args[[9]])

errQuit <- function(mesg, status=1) {
  message("Error: ", mesg)
  q(status=status)
}

### VALIDATE ARGUMENTS ###
# Input directory is expected to contain .fastq file(s)
if(!(file.exists(R1_path) && file.exists(R2_path))) {
  errQuit("Input files do not exist.")
}
if(R1_path == R2_path) {
  errQuit("Duplicate .fastq file in input")
}
# Output path is to be a directory name
if(!(dir.exists(output_dir))) {
  if (substr(output_dir, nchar(output_dir), nchar(output_dir)) =="/"){
    suppressWarnings(dir.create(output_dir))
  } else {
    suppressWarnings(dir.create(paste(output_dir, "/", sep="")))
    output_dir<-paste(output_dir, "/", sep="")
  } 
}

if(!(dir.exists(output_dir))) { 
  errQuit("Output path does not exist.")
}

# Check if fqgrep is installed
if(!system('fqgrep -V', ignore.stdout = T)==0) {
  errQuit("fqgrep is not installed in the system")
}

# Generate temporary folder
dir.create(file.path(output_dir, "temp"))
temp_dir<-paste(output_dir, "temp/", sep="")

# Check if packages are installed and load them
if(suppressWarnings(suppressMessages(!require(stringr))))
{install.packages("stringr", dependencies = T, repos='https://cloud.r-project.org/')
  suppressMessages(require(stringr))}

if(suppressWarnings(suppressMessages(!require(dplyr))))
{install.packages("dplyr", dependencies = T, repos='https://cloud.r-project.org/')
  suppressMessages(require(dplyr))}

if(suppressWarnings(suppressMessages(!require(parallel))))
{install.packages("parallel", dependencies = T, repos='https://cloud.r-project.org/')
  suppressMessages(require(parallel))}

#if(suppressWarnings(suppressMessages(!require(ShortRead))))
#{source("http://www.bioconductor.org/biocLite.R")
#  biocLite("ShortRead", suppressUpdates = T)
#  suppressMessages(require(ShortRead))}

#### Display Multithreading Warning #### 
message("The script runs at maximum 4 cores:")
message("Multithreading enabled on ", ifelse(nthreads>4,4, nthreads), " cores")

# Check format and number of reads
if(all(str_split_fixed(R1_path, pattern="\\.", 2)[,2]=="fastq" && str_split_fixed(R2_path, pattern="\\.", 2)[,2] == "fastq")){
  message("fastq reads detected\n", "R1: ", R1_path, "\nR2: ", R2_path, "\noutput set to .fastq")
  output_extension<-"fastq"
  fqF<-system(paste("echo $(cat ", R1_path, "|wc -l)/4|bc", paste=""), intern=T)
  fqR<-system(paste("echo $(cat ", R2_path, "|wc -l)/4|bc", paste=""), intern=T)
} else{
  if(all(str_split_fixed(R1_path, pattern="\\.", 2)[,2]=="fastq.gz" && str_split_fixed(R2_path, pattern="\\.", 2)[,2] == "fastq.gz")){
    message("fastq.gz reads detected\n", "R1: ", R1_path, "\nR2: ", R2_path, "\noutput set to .fastq.gz")
    output_extension<-"gz"
    fqF<-system(paste("echo $(zcat ", R1_path, "|wc -l)/4|bc", paste=""), intern=T)
    fqR<-system(paste("echo $(zcat ", R2_path, "|wc -l)/4|bc", paste=""), intern=T)
  } else{
    errQuit("input read files extension unknown")}
  }

# Print number of Reads in input and check if reads in R1 == reads in R2
# Stream on input .fastq files
#fqF <- FastqStreamer(R1_path, 1e6)
#fqR <- FastqStreamer(R2_path, 1e6)

# Check length
#if((R1<-length(yield(fqF))) != (R2<-length(yield(fqR)))) {
#  errQuit("R1 and R2 contain a different number of reads")
#} else {
#  message("Reads in input\n", "R1: ", R1, "\nR2: ", R2)
#}
#close(fqF, fqR)
if(fqF != fqR) {
  errQuit("R1 and R2 contain a different number of reads")
  } else {
    message("Reads in input\n", "R1: ", fqF, "\nR2: ", fqR)
  }

######### Re-orient paired-end reads and  ########
######### trim variable Ns before barcode ########
#### fqgrep primer sequence on R1 and R2 reads allowing mismatches and generate report ####
# R1
R1_fw_fq<-paste("fqgrep -c -e -r -p ", "'", fw_primer, "' -m ", primer_mismatch_fw, " ", R1_path, ">", paste(temp_dir, "FW_R1_report.txt", sep=""), sep="")
# R2
R2_rev_fq<-paste("fqgrep -c -e -r -p ", "'", rev_primer, "' -m ", primer_mismatch_rev, " ", R2_path, ">", paste(temp_dir, "REV_R2_report.txt", sep=""), sep="")

#### fqgrep primer sequence on R1 and R2 reads inverting FW and REV ####
# R2
R2_fw_fq<-paste("fqgrep -c -e -r -p ", "'", fw_primer, "' -m ", primer_mismatch_fw, " ", R2_path, ">", paste(temp_dir, "FW_R2_report.txt", sep=""), sep="")
# R1
R1_rev_fq<-paste("fqgrep -c -e -r -p ", "'", rev_primer, "' -m ", primer_mismatch_rev, " ", R1_path, ">", paste(temp_dir, "REV_R1_report.txt", sep=""), sep="")

# generate reports with multicore fqgrep commands
message("Matching primers...")
mc_fqgrep<-mclapply(list(R1_fw_fq, R2_rev_fq, R2_fw_fq, R1_rev_fq), function(list) system(list), mc.cores=nthreads, mc.preschedule=T, mc.cleanup=T)

# check report files dimension
report_dim<-unlist(lapply(list.files("/home/matteo/Scrivania/Work_on_Server/Olivo/Trimmed/temp/", full.names = T), function(x) file.size(x)))
if(sum(report_dim[2],report_dim[3])==0){
  errQuit("Reads are already in correct orientation")
}
# import report files
reports_list<- mclapply(list.files(temp_dir, full.names=T), function(x) read.delim(x, header=T), mc.cores=nthreads, mc.preschedule=T, mc.cleanup=T)

# Set trim_Ns function
# x is a list of report files
trim_Ns<-function(x) {
  # Discard spurious matches >30)
  report<- filter(x, start.position<=30)
  # set trim position before first barcode nucleotide (start position - barcode_length)
  report<-mutate(report, barcode.start=start.position-barcode_length)
  # isolate adapter (primer start-8)
  report<-mutate(report, barcode.seq=substr(substr(sequence, 1, start.position), barcode.start+1, start.position))
  # isolate the remaining part of the read
  report<- mutate(report, seq=as.character(str_split_fixed(sequence, pattern = "0m", n=2)[,2]))
  # compose trimmed reads in sequence field
  report<-mutate(report, sequence=paste(barcode.seq, match.string, seq, sep=""))
  # adjust qualities
  report<-mutate(report, quality=substr(quality, barcode.start+1, str_length(quality)))
  return(report)
}

# Launch trimming command
message("Trimming...")                        
trimmed_list<-mclapply(reports_list, function(x)
  trim_Ns(x), mc.cores=nthreads, mc.preschedule=T, mc.cleanup=T)

# merge FW and REV
message("Coupling R1 with R2...")                       
FW_REV_R1<-merge(trimmed_list[[1]], trimmed_list[[4]], by="read.name")
FW_REV_R1[1]<-paste(FW_REV_R1$read.name,FW_REV_R1$read.comments.x, sep=" ")
FW_REV_R1[28]<-paste(substr(FW_REV_R1$read.name, 1, str_length(FW_REV_R1$read.name)-str_length(FW_REV_R1$read.comments.x)),FW_REV_R1$read.comments.y, sep="")

FW_REV_R2<-merge(trimmed_list[[2]], trimmed_list[[3]], by="read.name")
FW_REV_R2[1]<-paste(FW_REV_R2$read.name,FW_REV_R2$read.comments.x, sep=" ")
FW_REV_R2[28]<-paste(substr(FW_REV_R2$read.name, 1, str_length(FW_REV_R2$read.name)-str_length(FW_REV_R2$read.comments.x)),FW_REV_R2$read.comments.y, sep="")

#### Remove paired reads without any barcodes
FW_REV_R1<-filter(FW_REV_R1, str_length(barcode.seq.x)>=barcode_length | str_length(barcode.seq.y)>=barcode_length)
FW_REV_R2<-filter(FW_REV_R2, str_length(barcode.seq.x)>=barcode_length | str_length(barcode.seq.y)>=barcode_length)

#### Assemble .fastq files tables ####
FW_REV_R1$read.name<-paste("@", FW_REV_R1$read.name, sep="")
FW_REV_R1$V28<-paste("@", FW_REV_R1$V28, sep="")
FW_REV_R1$sep<-"+"

FW_REV_R2$read.name<-paste("@", FW_REV_R2$read.name, sep="")
FW_REV_R2$V28<-paste("@", FW_REV_R2$V28, sep="")
FW_REV_R2$sep<-"+"

R1<- rbind(
  do.call(rbind, lapply(seq(nrow(FW_REV_R1[,c(1,10,29,11)])), function(i) t(FW_REV_R1[,c(1,10,29,11)][i, ]))),
  do.call(rbind, lapply(seq(nrow(FW_REV_R2[,c(1,10,29,11)])), function(i) t(FW_REV_R2[,c(1,10,29,11)][i, ])))
  )

R2<-rbind(
  do.call(rbind, lapply(seq(nrow(FW_REV_R1[,c(28,23,29,24)])), function(i) t(FW_REV_R1[,c(28,23,29,24)][i, ]))),
  do.call(rbind, lapply(seq(nrow(FW_REV_R2[,c(28,23,29,24)])), function(i) t(FW_REV_R2[,c(28,23,29,24)][i, ])))
  )

# Write Output
message("Writing output files...")
if(output_extension=="fastq"){
# Write re-oriented trimmed as .fastq
write.table(R1, file=paste(output_dir, str_replace(basename(R1_path), ".fastq", "_trimmed.fastq"), sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(R2, file=paste(output_dir, str_replace(basename(R2_path), ".fastq", "_trimmed.fastq"), sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE)
fqF_o<-system(paste("echo $(cat ", paste(output_dir, str_replace(basename(R1_path), ".fastq", "_trimmed.fastq"), sep=""), "|wc -l)/4|bc", paste=""), intern=T)
fqR_o<-system(paste("echo $(cat ", paste(output_dir, str_replace(basename(R2_path), ".fastq", "_trimmed.fastq"), sep=""), "|wc -l)/4|bc", paste=""), intern=T)
} else {
# Write re-oriented trimmed as .fastq.gz files
gz1 = gzfile(paste(output_dir, str_replace(basename(R1_path), ".fastq.gz", "_trimmed.fastq.gz"), sep=""),"w")
write(R1, gz1)
close(gz1)

gz2 = gzfile(paste(output_dir, str_replace(basename(R2_path), ".fastq.gz", "_trimmed.fastq.gz"), sep=""),"w")
write(R2, gz2)
close(gz2)
fqF_o<-system(paste("echo $(zcat ", paste(output_dir, str_replace(basename(R1_path), ".fastq", "_trimmed.fastq"), sep=""), "|wc -l)/4|bc", paste=""), intern=T)
fqR_o<-system(paste("echo $(zcat ", paste(output_dir, str_replace(basename(R2_path), ".fastq", "_trimmed.fastq"), sep=""), "|wc -l)/4|bc", paste=""), intern=T)
}

# Print number of Reads in outputs and check if reads in R1 == reads in R2
#fqF_o <- FastqStreamer(paste(output_dir, str_replace(basename(R1_path), ".fastq", "_trimmed.fastq"), sep=""), 1e6)
#fqR_o <- FastqStreamer(paste(output_dir, str_replace(basename(R2_path), ".fastq", "_trimmed.fastq"), sep=""), 1e6)

#if((R1<-length(yield(fqF_o))) != (R2<-length(yield(fqR_o)))) {
#  errQuit("An error occurred, R1/R2 output contain a different number of reads")
#} else {
#  message("Reads in output\n", "R1: ", R1, "\nR2: ", R2)
#}
#close(fqF_o, fqR_o)

#fqF_o<-system(paste("echo $(zcat ", paste(output_dir, str_replace(basename(R1_path), ".fastq", "_trimmed.fastq"), sep=""), "|wc -l)/4|bc", paste=""))
#fqR_o<-system(paste("echo $(zcat ", paste(output_dir, str_replace(basename(R2_path), ".fastq", "_trimmed.fastq"), sep=""), "|wc -l)/4|bc", paste=""))

if(fqF_o != fqR_o) {
  errQuit("An error occurred, R1 and R2 output contain a different number of reads")
} else {
  message("Reads in output\n", "R1: ", fqF_o, "\nR2: ", fqR_o)
}

stop.time<-proc.time()
message("Elapsed time (min):")
print(round((stop.time-start.time)/60,3))

# delete the temporary folder
unlink(file.path(output_dir, "temp"), recursive = T)

  # Save the Global Enviroment in the output directory for debug purposes
  #save.image(paste(output_dir, "debug.RData", sep=""))
q(status=0)

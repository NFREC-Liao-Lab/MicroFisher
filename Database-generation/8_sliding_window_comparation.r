#################################################################################
#sliding window calculation :https://cloud.tencent.com/developer/article/1511065
#################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


if (!require('seqinr')) install.packages('seqinr')
if (!require('ape')) install.packages('ape')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require('Biostrings')) BiocManager::install("Biostrings")
if (!require('ggplot2')) install.packages('ggplot2')
if (!require('ggseqlogo')) install.packages('ggseqlogo')
if (!require('seqinr')) install.packages('seqinr')
if (!require('gridExtra')) install.packages('gridExtra')
if (!require('seqLogo')) BiocManager::install("seqLogo")



#read the sequence
my_fasta <- read.fasta(file = "order_seq_aligned.aln", seqtype = "DNA", as.string = T)

#
my_fasta_string = vector(mode = 'character')
for (i in 1:length(my_fasta)){
  my_fasta_string[i] = toupper(c2s(my_fasta[[i]]))
}


seq_string <- DNAStringSet(my_fasta_string)


#generate the consensus matrix
consensusMatrix <- consensusMatrix(seq_string)
pwm_dash <- consensusMatrix[c("A","C","G","T","-"), ]
pwm_DNA_BASES <- consensusMatrix[DNA_BASES, ]


#calculate the prop table
#prop_DNA_BASES <-  prop.table(pwm_DNA_BASES, 2)
prop_dash <- prop.table(pwm_dash, 2)[5,]
prop_DNA_BASES <- prop.table(prop.table(pwm_dash, 2)[1:4,],2)


prop_DNA_BASES[is.na(prop_DNA_BASES)] <- 0.000000001
pwm <- makePWM(prop_DNA_BASES)




slotNames(pwm)
pwm(pwm)
ic(pwm)
consensus(pwm)

bit_score <- ic(pwm)


#1. trim the double end gaps
for (i in 1:200){
  if (prop_dash[i] <= 0.75){
    begin_base_site=i
    break
  }
}

for (i in length(prop_dash):(length(prop_dash)-200)){
  if (prop_dash[i] <= 0.75){
    finish_base_site=i
    break
  }
}


#assign the values for gaps
for (i in begin_base_site:finish_base_site){
  if (prop_dash[i] > 0.90) {bit_score[i] <- 0.8}
}


#sliding window calculate the best trimming position
######################


#calculate the trimming position of D1 region

D1D2_mid=begin_base_site+((finish_base_site-begin_base_site)/2)
if (as.integer(D1D2_mid) != as.numeric(D1D2_mid)  ){
  D1D2_mid=D1D2_mid+0.5
}


short_base_length_D1=140
base_site_D1 <- begin_base_site:D1D2_mid
sum_min=10000000
start_min=0
end_min=0
for (i in base_site_D1) {
  start=i
  end=start+short_base_length_D1+sum(bit_score[start:(start+140)]==0.8)
  sum_i = sum(bit_score[start:end])
  
  if(sum_i < sum_min){
    sum_min = sum_i
    start_min=start
    end_min=end
  } else {
    sum_min = sum_min
    start_min=start_min
    end_min=end_min
  }
}

cut_start_1=start_min
cut_end_1=start_min+140
Bits_mean_1=sum(bit_score[cut_start_1:cut_end_1])/140
length_D1=cut_end_1-cut_start_1

#calculate the trimming position of D2 region
D1D2_mid=begin_base_site+((finish_base_site-begin_base_site)/2)
if (as.integer(D1D2_mid) != as.numeric(D1D2_mid)  ){
  D1D2_mid=D1D2_mid+0.5
}

short_base_length=240
seq_length=finish_base_site-begin_base_site

base_site <- finish_base_site:D1D2_mid
sum_min=10000000
start_min=0
end_min=0


for (i in base_site) {
  end=i
  start=end-short_base_length-sum(bit_score[start:(start+240)]==0.8)
  sum_i = sum(bit_score[start:end])
  
  if(sum_i < sum_min){
    sum_min = sum_i
    start_min=start
    end_min=end
  } else {
    sum_min = sum_min
    start_min=start_min
    end_min=end_min
  }
}


#write the start and end positions
cut_start_2=start_min
cut_end_2=start_min+140
Bits_mean_2=sum(bit_score[cut_start_2:cut_end_2])/140
length_D2=cut_end_2-cut_start_2

cut_start_3=end_min-140
cut_end_3=end_min
Bits_mean_3=sum(bit_score[cut_start_3:cut_end_3])/140
length_D3=cut_end_3-cut_start_3

write.table(cut_start_1, file="cut_start_1.txt")
write.table(cut_end_1, file="cut_end_1.txt")
write.table(Bits_mean_1, file="Bits_mean_1.txt")
write.table(length_D1, file="length_D1.txt")
write.table(cut_start_2, file="cut_start_2.txt")
write.table(cut_end_2, file="cut_end_2.txt")
write.table(Bits_mean_2, file="Bits_mean_2.txt")
write.table(length_D2, file="length_D2.txt")
write.table(cut_start_3, file="cut_start_3.txt")
write.table(cut_end_3, file="cut_end_3.txt")
write.table(Bits_mean_3, file="Bits_mean_3.txt")
write.table(length_D3, file="length_D3.txt")
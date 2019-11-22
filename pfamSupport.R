
library(data.table)
library(Biostrings)

pfamA31 = readAAStringSet("~/Downloads/Pfam-A31.fasta")
pfamA32 = readAAStringSet("~/Downloads/Pfam-A32.fasta")


nms31 = names(pfamA31)
nms32 = names(pfamA32)

# set the names of PFAM-A entries read in with readAAStringset() to PFAM accession
nms31 =  gsub(".* .* |;.*","",nms31)
nms32 =  gsub(".* |;.*","",nms32)

# get the version of pfam in which the family was added
vers31 = gsub(".*\\.","",nms31)


nms31 = gsub("\\..*","",nms31)


pfamTable32 = data.table(sequence = as.data.table(pfamA32),width = width(pfamA32),name = names(pfamA32),family = nms32)
setkey(pfamTable32,width)


res = list()
unms32 = unique(pfamTable32$family[pfamTable32$width >= 18])
print(paste("anzahl gesamt: ",length(unms32)))
for(l in 1:length(unms32)){

    fam = list()
    nm = list()
    print(paste("gesamt zaehler: ",l))
    j = 1
    print(paste("anzahl Familie: ",length(which(pfamTable32$family == unms32[l]))))

    for(i in which(pfamTable32$family == unms32[l])){
        print(paste("fam zaehler: ",j))
        tmp = (grepl(pfamTable32$sequence.x[i],pfamTable32$sequence.x[i+1:length(pfamTable32$width)]))

        fam[[j]] = (pfamTable32$family[i+1:length(pfamTable32$width)])[tmp]
        nm[[j]] = pfamTable32$name[i]
        j = j + 1
    }

    res[[l]] = list(fam,nm)
}

res2 = list()

for(i in 1:length(res)){
    partOfFam = c()
    famm = c()
    z = rle(sort(res[[i]][[1]]))
    for(j in 1:length(z)){
        if(z$value[j] != unms32[i]){
            partOfFam[j] = z$length[j]/length(pfamTable32$name[pfamTable$name == z$value[j]])
            famm[j] = z$value[j]
        }
    }
    res2[[i]] = list(partOfFam,famm)
}

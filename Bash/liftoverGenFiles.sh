#Script to lift over gen files per chromosome. Note requires specific naming of files.

intputDir=""
outputDir=""
tmpDir=""

mkdir -p ${outputDir}
mkdir -p ${tmpDir}

for chr in {1..22} "X" "XY" "Y"
do
    
    awk '
    
        BEGIN {
            OFS="\t"
        }
    
        {
            
            if($1=="XY"){
                $1="Y"
            }
            
            print "chr"$1,$3,$3+1,$2
            
        }
    
    ' < ${intputDir}/chr${chr} > ${tmpDir}/chr${chr}_b36.bed
    
    
    ./liftOver \
	-bedPlus=4  ${tmpDir}/chr${chr}_b36.bed \
	./hg18ToHg19.over.chain.gz \
	${tmpDir}/chr${chr}_b37.bed \
	${tmpDir}/chr${chr}_b37.unmapped.txt
    
    awk -v mappingFile=${tmpDir}/chr${chr}_b37.bed '
    
        BEGIN{
            OFS=" "
            split("", mappings)
            while( (getline < mappingFile) > 0){
                mappings[$4] = $2
            }
        }
        
        $2 in mappings{
            $3 = mappings[$2]
            print $0
        }

    
    ' < ${intputDir}/chr${chr} > ${outputDir}/chr${chr}
    
    cp ${intputDir}/chr${chr}.sample ${outputDir}/chr${chr}.sample
    
done

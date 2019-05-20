
#shell.prefix("source activate snakemake_py2")
shell.prefix(" module load bioinfo; module load samtools; module load r")
import os
if not os.path.exists('logs'):
    os.makedirs('logs')

    
localrules: cordinate, outerMostUTRCor, merge_samples, largest_extension, merge_counts, DESeq2


## Get the sample names provided in the sample.txt and make a list
smlst = []
sampleFile = config["sample"]

with open(sampleFile) as f:
	for line in f:
		smlst.append(line.strip())



rule all:
    input:#"utr_OuterCor.txt",
            #expand(config["geneicDir"] + "{sample}.geneic.sorted.bam", sample=smlst),
	    #expand(config["totalReadsDir"] + "{sample}.totalReads.txt", sample=smlst),	
            #expand(config["outputDir"]+ "{sample}_five_prime_extension_plus.txt", sample=smlst),
            #expand(config["outputDir"] + "{sample}_five_prime_extension_minus.txt", sample=smlst),
            #expand(config["outputDir"] + "{sample}_three_prime_extension_plus.txt", sample=smlst),
            #expand(config["outputDir"] + "{sample}_three_prime_extension_minus.txt", sample=smlst),
            #config["outputDir"] + "five_prime_final_plus.txt",
            #config["outputDir"] + "five_prime_final_minus.txt",
            #config["outputDir"]  + "three_prime_final_plus.txt",
            #config["outputDir"]  + "three_prime_final_minus.txt"
            #config["UTRFile"]
            #expand(config["countsDir"] + "/{sample}_counts.txt", sample=smlst),
            #config["finalUTRFile"],
            config["finalUTR"]


## Extract the cordinates of UTRs and genes          
rule cordinate:
	input:
		gtf = config["reference"]["gtf"]
	output:
		utr = "utrCor.txt",
		gene = "gene.bed"
	shell:"""
                grep UTR {input.gtf} > {output.utr}
		grep 'gene'  {input.gtf} | awk '{{print $1"\t"$4"\t"$5}}' > {output.gene}

            """
	
## Get the outer most cordinates of UTRs as one gene has multiple UTRs
	
rule outerMostUTRCor:
    input:
        rules.cordinate.output.utr
    output:
        "utr_OuterCor.txt"

    run:            
        def outer1(utr, strand):
            
            d = {}
            with open(str(input)) as f:
                for line in f:
                    ln = line.strip().split('\t')
                    gene = ln[-1].split(";")[-1]
                    if utr in line and strand in line:
                        gene1 =  ln[8].split(';')[1].split('.')[1]
                        lst = [ln[0], ln[2], ln[3], ln[6]]
                        strF = '\t'.join(lst)
                        d[strF] = gene1
            return d
        import os

        utrfile = str(output)
        if os.path.exists(utrfile):
            os.remove(utrfile)

        def outer2(utr, strand):
            
            d = {}
            dd = {}
            with open(str(input)) as f:
                for line in f:
                    ln = line.strip().split('\t')
                    gene = ln[-1].split(";")[-1]
                    if utr in ln and strand in line:
                        gene1 = ln[8].split(';')[1].split('.')[1]
                        lst = [ln[0], ln[2], ln[3], ln[4], ln[6], gene1]
                        if gene in d.keys():
                            d[gene].append(lst)
                        else:
                            d[gene] = [lst]
                for k, v in d.items():
                    vv = v[-1]
                    strF = '\t'.join([vv[0], vv[1], vv[3], vv[4]])
                    dd[strF] = vv[5]
            return dd

        myfile =  open(utrfile, 'a')
            
        for k, v in sorted(outer1("five_prime_UTR.1", "+").items()):
            myfile.write(k + '\t' + str(v) + '\n')
            
        for k, v in sorted(outer1("three_prime_UTR.1", "-").items()):
            myfile.write(k + '\t' + str(v) + '\n')
                
        for k, v in sorted(outer2("five_prime_UTR", "-").items()):
            myfile.write(k + '\t' + str(v) + '\n')                

        for k, v in sorted(outer2("three_prime_UTR", "+").items()):
            myfile.write(k + '\t' + str(v) + '\n')
            

## Generate the gene bam file with reads mapping only to geneic regions using gene.bed
            
rule gene_bam:
	input:	
		sample = config["inputDir"]+"/{sample}.bam",
		genebed = rules.cordinate.output.gene

	output:
		first = config["geneicDir"]+"/{sample}.geneic.bam",
		second = config["geneicDir"]+"/{sample}.geneic.sorted.bam"
	shell:"""
		ml samtools
		samtools view -hL {input.genebed} {input.sample} -@ 5 > {output.first} 
		samtools sort {output.first} -o {output.second} -@ 5
		samtools index {output.second}
		"""

## Get the total number of input reads for calculating the RPKM values
	
rule total_reads:
	input:
		sample = config["inputDir"]+"/{sample}.bam"
	output:
		totalReads = config["totalReadsDir"]+"/{sample}.totalReads.txt"
	shell:"""
		ml samtools
		samtools view -F 4 {input.sample} | grep 'NH:i:1' |awk '$5 == 50 {{print $5}}' | wc -l > {output.totalReads}
		"""

## Get the five prime UTR extensions on the minus strand
	
rule five_prime_UTR_minus:
    input:
        bam = config["inputDir"]+"/{sample}.bam",
        geneic = config["geneicDir"]+"/{sample}.geneic.sorted.bam",
        reads = config["totalReadsDir"]+"/{sample}.totalReads.txt",
    	cor =  rules.outerMostUTRCor.output
    output:
        out5minus = config["outputDir"]+ "/{sample}_five_prime_extension_minus.txt"
    run:
        def five_prime_UTR_minus(utr, strand, outputfile):
            import pysam
            import inspect
            entire = pysam.AlignmentFile(str(input.bam), "rb")
            geneic = pysam.AlignmentFile(str(input.geneic), "rb")
            cordinate = str(input.cor)
            #outfile = str(output.out)
            totalReads = open(str(input.reads))
            total = int(totalReads.readline().strip())
            myfile = open(outputfile,'w')
            with open(cordinate) as f:
                count = 0
                d = {}
                lst = []
                biglst = []
                name = 0
                for line in f:
                    if utr in line and strand in line:
                        #print(line)
                        ln = line.strip().split('\t')
                        bin  = 10
                        start = int(ln[2])
                        st = int(ln[2])
                        end = int(ln[2]) + bin
                        chr = ln[0]
                        gene = ln[4]
                        i = 0
                        tup = ('NH', 1)
                        if start > 10:
                            while i<1000:
                                for read in geneic.fetch(chr, start, end):
                                    #gene_obj =  dir(read)
                                    read_name_gene = read.query_name # getattr(read, gene_obj[90])
                                    flags_gene = read.tags #getattr(read, gene_obj[105])
                                    mapq_gene = read.mapping_quality #getattr(read, gene_obj[64])
                                    d[read_name_gene] = ''
                                for read1 in entire.fetch(chr, start, end):
                                    #entire_obj =  dir(read1)
                                    read_name_entire = read1.query_name
                                    flags_entire = read1.tags
                                    mapq_entire = read1.mapping_quality 
                                    if (read_name_entire not in d):
                                        if mapq_entire == 50 and tup in flags_entire:
                                            count+=1
                                d ={}
                                
                                RPM = (float(count) * 10**9)/(total *10)
                                if RPM >= 1:
                                    #print(chr, start, end,count, RPM)
                                    lst.append(end)
                                    start = start + bin
                                    end = end + bin
                                    biglst.append(lst)
                                else:
                                    break
                        
                                i+=1
                                count = 0
                    
                         
                            if len(lst)>= 1 and start > 10:
                                name+=1
                                dist = lst[-1] - st
                                final = [chr, str(st), str(lst[-1]), str(dist), 'five','-', gene]
                                strF =  '\t'.join(final)
                                #print(strF)
                                myfile.write(strF+ '\n')                
                                lst = []
            myfile.close()
        five_prime_UTR_minus("five_prime_UTR", "-", str(output.out5minus))


## Get the three prime UTR extensions on the plus strand
           
rule three_prime_UTR_plus:
        input:
                bam = config["inputDir"]+"/{sample}.bam",
                geneic = config["geneicDir"]+"/{sample}.geneic.sorted.bam",
		reads = config["totalReadsDir"]+"/{sample}.totalReads.txt",
		cor =  rules.outerMostUTRCor.output
        output:
                out3plus = config["outputDir"]+"/{sample}_three_prime_extension_plus.txt"
           
                
        run:
            def three_prime_UTR_plus(utr, strand, outputfile): 
                import pysam
                import inspect
                entire = pysam.AlignmentFile(str(input.bam), "rb")
                geneic = pysam.AlignmentFile(str(input.geneic), "rb")
                cordinate = str(input.cor)
                #outfile = str(output.out)
                totalReads = open(str(input.reads))
                total = int(totalReads.readline().strip())
                myfile = open(outputfile,'w')
                with open(cordinate) as f:
                    count = 0
                    d = {}
                    lst = []
                    biglst = []
                    name = 0
                    for line in f:
                        if utr in line and strand in line:
                            ln = line.strip().split('\t')
                            bin  = 10
                            start = int(ln[2])
                            st = int(ln[2])
                            end = int(ln[2]) + bin
                            chr = ln[0]
                            gene = ln[4]
                            i = 0 #print(line)
                            
                            tup = ('NH', 1)
                            if start > 10:
                                while i<1000:
                                    for read in geneic.fetch(chr, start, end):
                                        #gene_obj =  dir(read)
                                        read_name_gene = read.query_name # getattr(read, gene_obj[90])
                                        flags_gene = read.tags #getattr(read, gene_obj[105])
                                        mapq_gene = read.mapping_quality #getattr(read, gene_obj[64])
                                        d[read_name_gene] = ''
                                    for read1 in entire.fetch(chr, start, end):
                                        #entire_obj =  dir(read1)
                                        read_name_entire = read1.query_name
                                        flags_entire = read1.tags
                                        mapq_entire = read1.mapping_quality
                                        if (read_name_entire not in d):
                                            if mapq_entire == 50 and tup in flags_entire:
                                                count+=1
                                    d ={}
                        
                                    RPM = (float(count) * 10**9)/(total *10)
                                    if RPM >= 1:
                                        #print(chr, start, end,count, RPM)
                                        lst.append(end)
                                        start = start + bin
                                        end = end + bin
                                        biglst.append(lst)
                                    else:
                                        break
                        
                                    i+=1
                                    count = 0
                    
                         
                                if len(lst)>= 1 and start > 10:
                                    name+=1
                                    dist = lst[-1] - st
                                    final = [chr, str(st), str(lst[-1]), str(dist), 'three','+', gene]
                                    strF =  '\t'.join(final)
                                    #print(strF)
                                    myfile.write(strF+ '\n')
                
                                lst = []
                myfile.close()
            
            three_prime_UTR_plus("three_prime_UTR", "+", str(output.out3plus))


## Get the five prime UTR extensions on the plus strand

rule five_prime_UTR_plus:
        input:
                bam = config["inputDir"]+"/{sample}.bam",
                geneic = config["geneicDir"]+"/{sample}.geneic.sorted.bam",
		reads = config["totalReadsDir"]+"/{sample}.totalReads.txt",
		cor =  rules.outerMostUTRCor.output
        output:
                out5plus =  config["outputDir"]+"/{sample}_five_prime_extension_plus.txt",
                
        run:
            def five_prime_UTR_plus(utr, strand, outputfile): 
                import pysam
                import inspect
                entire = pysam.AlignmentFile(str(input.bam), "rb")
                geneic = pysam.AlignmentFile(str(input.geneic), "rb")
                cordinate = str(input.cor)
                #outfile = str(output.out)
                totalReads = open(str(input.reads))
                total = int(totalReads.readline().strip())
                myfile = open(outputfile,'w')
                with open(cordinate) as f:
                    count = 0
                    d = {}
                    lst = []
                    biglst = []
                    name = 0
                    for line in f:
                        if utr in line and strand in line:
                            #print(line)
                            ln = line.strip().split('\t')
                            bin  = 10
                            end = int(ln[2])
                            en = int(ln[2])
                            start = int(ln[2]) - bin
                            chr = ln[0]
                            gene = ln[4]
                            i = 0
                            tup = ('NH', 1)
                            if start > 10:
                                while i<1000:
                                    for read in geneic.fetch(chr, start, end):
                                        #gene_obj =  dir(read)
                                        read_name_gene = read.query_name # getattr(read, gene_obj[90])
                                        flags_gene = read.tags #getattr(read, gene_obj[105])
                                        mapq_gene = read.mapping_quality #getattr(read, gene_obj[64])
                                        d[read_name_gene] = ''
                                    for read1 in entire.fetch(chr, start, end):
                                        #entire_obj =  dir(read1)
                                        read_name_entire = read1.query_name
                                        flags_entire = read1.tags
                                        mapq_entire = read1.mapping_quality
                                        if (read_name_entire not in d):
                                            if mapq_entire == 50 and tup in flags_entire:
                                                count+=1
                                    d ={}
                        
                                    RPM = (float(count) * 10**9)/(total *10)
                                    if RPM >= 1 and start > 10:
                                        #print(chr, start, end,count, RPM)
                                        lst.append(end)
                                        start = start - bin
                                        end = end - bin
                                        biglst.append(lst)
                                    else:
                                        break
                        
                                    i+=1
                                    count = 0
                    
                         
                                if len(lst)>= 1:
                                    name+=1
                                    dist = en - lst[-1]
                                    final = [chr, str(lst[-1]), str(en),  str(dist), 'five', '+', gene]
                                    strF =  '\t'.join(final)
                                    myfile.write(strF+ '\n')
                
                                lst = []
                myfile.close()
            five_prime_UTR_plus("five_prime_UTR", "+", str(output.out5plus))


## Get the three prime UTR extensions on the minus strand

rule three_prime_UTR_minus:
        input:
                bam = config["inputDir"]+"/{sample}.bam",
                geneic = config["geneicDir"]+"/{sample}.geneic.sorted.bam",
		reads = config["totalReadsDir"]+"/{sample}.totalReads.txt",
		cor =  rules.outerMostUTRCor.output
        output:
                out3minus =  config["outputDir"]+"/{sample}_three_prime_extension_minus.txt"
                
        run:
            def three_prime_UTR_minus(utr, strand, outputfile): 
                import pysam
                import inspect
                entire = pysam.AlignmentFile(str(input.bam), "rb")
                geneic = pysam.AlignmentFile(str(input.geneic), "rb")
                cordinate = str(input.cor)
                #outfile = str(output.out)
                totalReads = open(str(input.reads))
                total = int(totalReads.readline().strip())
                myfile = open(outputfile,'w')
                with open(cordinate) as f:
                    count = 0
                    d = {}
                    lst = []
                    biglst = []
                    name = 0
                    for line in f:
                        if utr in line and strand in line:
                            #print(line)
                            ln = line.strip().split('\t')
                            bin  = 10
                            end = int(ln[2])
                            en = int(ln[2])
                            start = int(ln[2]) - bin
                            chr = ln[0]
                            gene = ln[4]
                            i = 0
                            tup = ('NH', 1)
                            if start > 10:
                                while i<1000:
                                    for read in geneic.fetch(chr, start, end):
                                        #gene_obj =  dir(read)
                                        read_name_gene = read.query_name # getattr(read, gene_obj[90])
                                        flags_gene = read.tags #getattr(read, gene_obj[105])
                                        mapq_gene = read.mapping_quality #getattr(read, gene_obj[64])
                                        d[read_name_gene] = ''
                                    for read1 in entire.fetch(chr, start, end):
                                        #entire_obj =  dir(read1)
                                        read_name_entire = read1.query_name
                                        flags_entire = read1.tags
                                        mapq_entire = read1.mapping_quality
                                        if (read_name_entire not in d):
                                            if mapq_entire == 50 and tup in flags_entire:
                                                count+=1
                                    d ={}
                        
                                    RPM = (float(count) * 10**9)/(total *10)
                                    if RPM >= 1  and start > 10:
                                        #print(chr, start, end,count, RPM)
                                        lst.append(end)
                                        start = start - bin
                                        end = end - bin
                                        biglst.append(lst)
                                    else:
                                        break
                        
                                    i+=1
                                    count = 0
                    
                         
                                if len(lst)>= 1 and start > 10:
                                    name+=1
                                    dist = en - lst[-1]
                                    final = [chr, str(lst[-1]), str(en),  str(dist), 'three', '-', gene]
                                    strF =  '\t'.join(final)
                                    myfile.write(strF+ '\n')
                
                                lst = []
                myfile.close()
            three_prime_UTR_minus("three_prime_UTR", "-", str(output.out3minus))


## Merge all the UTRs into one file 
rule merge_samples:
                    input:
                        utr5plus = expand(config["outputDir"] + "/{sample}_five_prime_extension_plus.txt", sample=smlst),
                        utr5minus = expand(config["outputDir"] + "/{sample}_five_prime_extension_minus.txt", sample=smlst),
                        utr3plus = expand(config["outputDir"] +  "/{sample}_three_prime_extension_plus.txt", sample=smlst),
                        utr3minus = expand(config["outputDir"] +  "/{sample}_three_prime_extension_minus.txt", sample=smlst)
                    output:
                            merged5plus = config["outputDir"] + "/five_prime_extension_plus.txt",
                            merged5minus = config["outputDir"] + "/five_prime_extension_minus.txt",
                            merged3plus = config["outputDir"] + "/three_prime_extension_plus.txt",
                            merged3minus = config["outputDir"] + "/three_prime_extension_minus.txt"
                            
                    shell:"""
                            cat {input.utr5plus} | sort -u > {output.merged5plus}
                            cat {input.utr5minus}  | sort -u > {output.merged5minus}
                            cat {input.utr3plus}  | sort -u > {output.merged3plus}
                            cat {input.utr3minus} | sort -u > {output.merged3minus}
                            """

## get the largest Extension for each gene    
rule largest_extension:
                        input:
                            merged5plus = config["outputDir"] +  "/five_prime_extension_plus.txt",
                            merged5minus = config["outputDir"] + "/five_prime_extension_minus.txt",
                            merged3plus = config["outputDir"] +  "/three_prime_extension_plus.txt",
                            merged3minus = config["outputDir"] +  "/three_prime_extension_minus.txt"

                        output: config["UTRFile"]
                        run:
                            finalfile = "{output}"
                            if os.path.exists(finalfile):
                                os.remove(finalfile)
                            def merged(infile, outfile, utr):
                                d = {}
                                myfile = open(outfile, 'a')
                                #myfile.write('UTRID\tChr\tStart\tEnd\tLength\tUTR\tStrand\n')
                                with open(infile) as f:
                                    for line in f:
                                        ln = line.strip().split('\t')
                                        key = ln[6]
                                        val = line.strip().split('\t')
                                        if key in d:
                                            d[key].append(val)
                                        else:
                                            d[key] = [val]
                                    

                                    for k, v in d.items():
                                        dist = [int(x[3]) for x in v]
                                        maxD =  max(dist)
                                        maxCor = v[dist.index(maxD)]
                                        name = k + utr
                                        lst = [name] + maxCor[:-1]
                                        strF = '\t'.join(lst)
                                        myfile.write(strF+'\n')
                                        
                                myfile.close()
                            merged(str(input.merged5plus), str(output), '.5UTR')
                            merged(str(input.merged5minus), str(output), '.5UTR')
                            merged(str(input.merged3plus), str(output), '.3UTR')
                            merged(str(input.merged3minus), str(output), '.3UTR')


## Get the counts for each extension

rule get_counts_UTR:
                    input:
                        bam = config["inputDir"] + "/{sample}.bam",
                        geneic = config["geneicDir"] + "/{sample}.geneic.sorted.bam",
                        cordinate = rules.largest_extension.output
                    output:
                        counts= config["countsDir"] + "/{sample}_counts.txt"
                    run:
                        import pysam
                        import inspect
                        bam = pysam.AlignmentFile(str(input.bam), "rb")
                        geneic = pysam.AlignmentFile(str(input.geneic), "rb")
                       
                        def counts():
                            myfile = open(str(output.counts),'w')
                            myfile.write('UTRID\tChr\tStart\tEnd\tLength\tReads\n')
                            with open(str(input.cordinate)) as f:
                                #f.readline()
                                for line in f:
                                    ln = line.strip().split('\t')
                                    start = int(ln[2])
                                    end = int(ln[3])
                                    chr = ln[1]
                                    id = ln[0]
                                    length = ln[4]
                                    count = 0
                                    tup = ('NH', 1)
                                    dgene = {}
                                    for read1 in geneic.fetch(chr, start, end):
                                        read_name_gene = read1.query_name
                                        flags_gene = read1.tags
                                        mapq_gene = read1.mapping_quality
                                        dgene[read_name_gene] = ''

                                    for read in bam.fetch(chr, start, end):
                                        read_name_utr = read.query_name
                                        flags_utr = read.tags
                                        mapq_utr = read.mapping_quality
                                        if read_name_utr not in dgene and mapq_utr == 50 and tup in flags_utr:
                                            #print mapq_gene, flags_gene
                                            count+=1
                                    lst = [id, chr, str(start), str(end), str(length), str(count)]
                                    strF = '\t'.join(lst)
                                    myfile.write(strF+'\n')
                                    count = 0   
                            myfile.close()

    
                        counts()
                        

## Merge the all the samples counts into one file

rule merge_counts:
                input:
                    cnt = list(expand(config["countsDir"] + "/{sample}_counts.txt", sample=smlst)),
                    utrfile = rules.largest_extension.output
                    
                output: config["finalUTRFile"]
                run:
                    import os
                    import pandas as pd
                    import functools
                    files = input.cnt
                    names =  ["UTRID"] + [file.replace('Counts/','').split('_')[0] for file in files]
                    fields = ['UTRID', 'Reads']
                    dataframes = [ pd.read_csv( f, sep ='\t', usecols=fields) for f in files ] 
                    merged = functools.reduce(lambda left,right: pd.merge(left,right,on='UTRID', how='outer'), dataframes)
                    merged.columns = names
                    colnames = ['UTRID', 'Chr', 'Start','End','Length','UTR','Strand']
                    utr = pd.read_csv(str(input.utrfile), sep ='\t', encoding='latin-1',names=colnames)
                    utrFrame = pd.DataFrame(utr, columns = colnames)
                    #print (utrFrame.head())
                    #print (merged.head())
                    newLst = [utrFrame, merged]
                    newMerged = pd.merge(utr, merged, on='UTRID', how='outer')
                    newMerged.to_csv(str(output), sep='\t', index= False)



## Perform statistical test using DESeq2 to get the treatment responsice UTRs that would be the real UTRs.
                    
rule DESeq2:
        input:
                infile = rules.merge_counts.output,
                meta = config["metadata"]
        output:
                config["finalUTR"]
        shell:"""
                Rscript DESEq2.r {input.meta} {input.infile} {output}
                """

        

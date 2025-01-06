DNA双端测序原理illumina：   
(https://blog.csdn.net/keepaware/article/details/114486031)   
(https://blog.csdn.net/weixin_43843918/article/details/135899404)   
RNA-seq:(https://www.bilibili.com/video/BV1XJ411r7bJ/?from=search&seid=9478135992416820142&spm_id_from=333.788.comment.all.click&vd_source=3720075e67a7b3080234abc490244540)

# 工具
sratoolkit：3.1.1   
fastqc：0.12.1   
multiqc：1.26   
cutadapt：5.0   
trimgalore：0.6.3   
fastp：0.24.0   
trimmomatic：0.39   
hisat2：2.2.1   
sortmerna：2.1  
samtools：1.21   
HTseq：2.0.9   
R：4.4.2   
parallel：202411122   
stringtie：2.2.3   
R包：   
    ballgown：2.38.0    
    DESeq2：1.46.0     
    pheatmap：1.0.12   
    biomaRt：2.62.0    
    org.Rn.eg.db：3.20.0    
    clusterProfiler：4.14.4    
    ggplot2:3.5.1    
    

## 3.1.参考数据下载步骤中：
大鼠ensembl下载ref部分信息：
```
>1 dna:primary_assembly primary_assembly:mRatBN7.2:1:1:260522016:1 REF
>2 dna:primary_assembly primary_assembly:mRatBN7.2:2:1:249053267:1 REF
>3 dna:primary_assembly primary_assembly:mRatBN7.2:3:1:169034231:1 REF
>4 dna:primary_assembly primary_assembly:mRatBN7.2:4:1:182687754:1 REF
>5 dna:primary_assembly primary_assembly:mRatBN7.2:5:1:166875058:1 REF
```
染色体编号后描述信息（组装版本、长度等等），影响后续脚本统计或分析，需要先去掉
```
# 去除染色体编号后的描述信息
cat rn6.raw.fa | perl -n -e 'if(m/^>(.+?)(?:\s|$)/){ print ">$1\n";}else{print}' > rn6.fa
```
代码解释：   
`m/^>(.+?)(?:\s|$)/`：`m/../`match匹配操作符（后接正则表达式），`^>`以>开头，捕获`>`后面的内容直到空白字符或行尾，存储在`$1`   
`(.+?)`：`.`匹配任意字符（除换行符）、`+`前一字符至少匹配一次、`？`非贪婪修饰符（尽可能少地匹配字符），`(..)`捕获组（匹配的内容被存储，后续`$1`引用）   
`(?:\s|$)`中，`(?:)`非捕获组（匹配内容分组但不存储），`\s`空白字符，`｜`逻辑或（匹配空白字符/字符串结束符），`$`字符串结束符   
`print ">$1\n";`如果前面匹配成功，将捕获的内容前面`$1`加上`>`、后面加上换行符`\n`输出
`else{print}`没匹配成功，则原样输出    

`(.+?)`逐字符捕获匹配（`？`非贪婪模式使得一旦遇到后面非捕获组成功匹配到空白/行尾，匹配就停止），`(?:\s|$)`限制`(.+?)`匹配到空格/行尾结束
输出结果示例：
```
>seq1  #描述信息删除了
ATGCATGCATGC
>seq2
GCTAGCTAGCTA
```

**统计染色体长度的脚本**
```
cat rn6.fa | perl -n -e '
    s/\r?\n//;
    if(m/^>(.+?)\s*$/){
        $title = $1;
        push @t, $title;
    }elsif(defined $title){
        $title_len{$title} += length($_);
    }
    END{
        for my $title (@t){
            print "$title","\t","$title_len{$title}","\n";
        }
    }
'
```
解释：   
- `s/\r?\n//`替换字符串内容。具体为删除回车和换行符
    `s/../../`替换操作符。`\r?`匹配可选的回车符`\r`。`\n`换行符
- `m/^>(.+?)\s*$/`匹配、提取标题行并存储在`$1`中
    `\s*`匹配可选的空白字符，`$`行尾
- `$title = $1`捕获的`$1`赋值给`$title1`
- `push @t, $title`将`@title`存入数组`@t`
    `push`perl的数组操作，将元素追加到数组末尾。
- `elsif(defined $title)`
    `elsif`是if的扩展，“否则”。`defined`判断变量是否已经定义，
- `$title_len{$title} += length($_)`
    `$title_len{...}`哈希的键值对，`+=`累加，`length($_)`返回当前行$_的长度



## 3.2.测试数据/实验数据下载步骤
用SRAtoolkit中prefetch下载：
```
# 后台下载
$ nohup prefetch SRR2190795 SRR224018{2..7} SRR2240228 --output-directory . &
```
代码解释：   
`nohup`命令，让任务在后台运行且不受终端关闭的影响。`prefectch`从SRA下载序列数据。`SRR2190795 SRR224018{2..7} SRR2240228`SRA ID号。    
`--output-directory`output输出目录，`.`当前路径。`{2..7}`大括号扩展，一系列连续数字。`&`后台运行符（放在整个命令的最后）   
下载结果是.sra文件，转换格式fastq-dump
```
parallel -j 4 "
    fastq-dump --split-3 --gzip {1}
" ::: $(ls *.sra)
```
解释：   
`parallel`并行化工具，同时运行多个任务。`-j 4`jobs指定任务数为4个。   
`"fastq-dump --split-3 --gzip {1}"`parallel要执行的每个任务。   
`fastq-dump`从.sra文件提取fastq格式序列数据。`--split-3`双端数据输出两个文件、单端数据输出一个文件。`--gzip`生成的fastq进行gzip压缩   
`{1}`parallel的占位符，表示传递给当前任务的第一个输入参数（.sra文件路径）  
`:::`分隔命令和输入参数（`$(ls *.sra)`）

**FASTQ格式介绍**
```
@SRR2190795.1 HWI-ST1147:240:C5NY7ACXX:1:1101:1320:2244 length=100 
ATGCTGGGGGCATTAGCATTGGGTACTGAATTATTTTCAGTAAGAGGGAAAGAATCCATCTCCNNNNNNNNNNNNNNNNNNNNNNAAANAAAAATAAAAT
+SRR2190795.1 HWI-ST1147:240:C5NY7ACXX:1:1101:1320:2244 length=100 
CCCFFFFFHHHHHJIJJJJJJJJDHHJJJIJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJHH#####################################
```
- 第一行：
    `@SRR2190795.1`序列ID，`HWI-ST1147:240:C5NY7ACXX:1:1101:1320:2244`测序信息（仪器ID、通道号、行号、列号），`length=100`序列长度
- 第二行：DNA序列，ACTG碱基+N（无法确定的碱基）
- 第三行：`+`序列标识符，与第一行对应
- 第四行：
    质量分数字符串，每个字符对应第二行的一个碱基。质量分数衡量测序准确性，值越高测序越可靠，`C``F``H`映射质量分数再计算得出

## 4.1.质量评估
创建目录、fastqc质量评估：
```
mkdir -p ../output/fastqc
fastqc -t 6 -o ../output/fastqc *.gz
```
分析报告在fastqc目录下，多个合并为一个：
```
multiqc .
```
- 平均GC含量
- 序列质量：出现低于30需剔除
- 平均质量值的read数量：平均质量低于20需剔除
- 接头含量：部分序列包含接头，剔除

## 4.2.剔除接头和质量差的
### 使用`cutadapt`剔除接头：
```
cd ～/project/rat/sequence
for i in $(ls *.fastq.gz);
do
    cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
    --minimum-length 30 --overlap 4 --trim-n \
    -o ../output/adapter/${i}  ${i}
done
```
- `-a`去除正向read1的接头adapter P7   
- `--minimum-length`如果剔除接头后read长度低于30，这条read将会被丢弃    
- `--overlap` 如果两端的序列与接头有4个碱基的匹配将会被剔除    
- `--trim-n`剔除两端的N    

### 使用`trimmomatic`去除低质量区域：
```
cd ~/project/rat/output/adapter
mkdir trim

parallel -j 4 "    trimmomatic  SE -phred33 {1} ../trim/{1} \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:5:15 MINLEN:30 \
" ::: $(ls *.gz)
```
- `SE`数据为单端测序single-end，双端PE   
- `-phred33`质量分数编码格式。`{1}`占位符（输入文件）。`../trim/{1}`输出路径和名称（与输入同名）   

- trimmomatic的剪切参数：
    `LEADING:20`从序列*开头*开始去掉质量值小于20的碱基
    `TRAILING:20`从序列*末尾*开始去掉质量值小于20的碱基
    `SLIDINGWINDOW:5:15`去除中间质量差的区域。从5'端开始以 5bp 的窗口计算碱基平均质量，若该值小于15，从这个位置截断read。
    `MINLEN:30`若read长度小于 30 bp 则丢去整条read

### 再次质量评估
```
cd ~/project/rat/output/trim
mkdir ../fastqc_trim
parallel -j 4 "
    fastqc -t 4 -o ../fastqc_trim {1}
" ::: $( ls *.gz)

cd ../fastqc_trim
multiqc .
```

# 5.去除rRNA序列
如果在提取RNA过程中没有对RNA进行筛选的情况下，那么得到的大部分将会是rRNA，用`sortmerna`去除rRNA序列    
ps：测序文件是为压缩格式
```
# sortmrna数据库中，定义数据库变量
sortmerna_ref_data=$(pwd)/rRNA_databases/silva-bac-16s-id90.fasta,$(pwd)/index/silva-bac-16s-db:\
$(pwd)/rRNA_databases/silva-bac-23s-id98.fasta,$(pwd)/index/silva-bac-23s-db:\
$(pwd)/rRNA_databases/silva-arc-16s-id95.fasta,$(pwd)/index/silva-arc-16s-db:\
$(pwd)/rRNA_databases/silva-arc-23s-id98.fasta,$(pwd)/index/silva-arc-23s-db:\
$(pwd)/rRNA_databases/silva-euk-18s-id95.fasta,$(pwd)/index/silva-euk-18s-db:\
$(pwd)/rRNA_databases/silva-euk-28s-id98.fasta,$(pwd)/index/silva-euk-28s-db:\
$(pwd)/rRNA_databases/rfam-5s-database-id98.fasta,$(pwd)/index/rfam-5s-db:\
$(pwd)/rRNA_databases/rfam-5.8s-database-id98.fasta,$(pwd)/index/rfam-5.8s-db

# 真核生物的rRNA
euk_rNRA_ref_data=$(pwd)/rRNA_databases/silva-euk-18s-id95.fasta,$(pwd)/index/silva-euk-18s-db:\
$(pwd)/rRNA_databases/silva-euk-28s-id98.fasta,$(pwd)/index/silva-euk-28s-db:\
$(pwd)/rRNA_databases/rfam-5s-database-id98.fasta,$(pwd)/index/rfam-5s-db:\
$(pwd)/rRNA_databases/rfam-5.8s-database-id98.fasta,$(pwd)/index/rfam-5.8s-db
```
```
cd ~/project/rat/output
mkdir -p ./rRNA/discard

cd trim

parallel -j 4 "
    gzip -d {1}*.fastq.gz
    sortmerna \
    --ref $euk_rNRA_ref_data \
    --reads {1}*.fastq \
    --aligned ../rRNA/discard/{1} \
    --other ../rRNA/{1} \
    --fastx \
    --log \
    -a 4 \
    -v
    gzip ../rRNA/{1}.fastq
    gzip ../rRNA/discard/{1}.fastq 
" ::: $( ls *.fastq.gz | perl -n -e 'print $1."\n" if m/(.+?)_/' ) 
```
解释：  
- `--aligned`与rRNA数据库能比对上的序列(说明是rRNA，后续不需要分析)
- `--other`与rRNA数据库不能比对上的序列(非rRNA，后续需要分析)
- `-a`线程数。`-v`verbose详细信息
ps：`--aligned`和`--other`后面接文件名的前缀就ok，不用加后缀.fastq，后缀在`--fastx`已经指定好了
- `print $1."\n" if m/(.+?)_/`处理.fastq.gz文件名，提取文件名中第一个`_`前的部分，可改成一下提取文件名前缀
    > print $1."\n" if m/(.+?)\./

# 6.序列比对
read比对定位到参考基因组的位置，确定read属于哪个基因，RNA-seq会出现跨范围比对（内含子/外显子差别）
## 6.1.建立索引
使用`hisat2`中的`hisat2-build`
```
hisat2-build [options] <参考基因组> <索引文件的前缀名称>
```
```
cd ~/project/rat/genome
mkdir index
cd index

hisat2-build -p 6 ../rn6.chr1.fa rn6.chr1
```
`-p`线程数
生成八个文件，**序列比对时直接用这八个文件**，不用rn6.chr1.fa文件   
## 6.2.序列比对
使用`hisat2`
```
hisat2 [options] -x <索引文件> < -1 1测序文件 -2 2测序文件 -U 未成对测序文件 > < -S 输出的sam文件>
```
```
cd ~/project/rat/output
mkdir align
cd rRNA

parallel -k -j 4 "
    hisat2 -t -x ../../genome/index/rn6.chr1 \
      -U {1}.fastq.gz -S ../align/{1}.sam \
      2>../align/{1}.log
" ::: $( ls *.gz | perl -p -e 's/.fastq.gz$//' )
```
`-p`逐行处理并打印结果，`-n`逐行处理不打印结果   
`s/.fastq.gz$//`中第一个`.`表示任意单个字符，整体表示匹配任意字符后面跟着fastq.gz并删除这个部分`fastq.gz`

- 总结比对情况：
```
cd ~/project/rat/output/align
file_list=($(ls *.log))

echo -e "sample\tratio\ttime"
for i in ${file_list[@]};
do
    prefix=$(echo ${i} | perl -p -e 's/\.log//')
    echo -n -e "${prefix}\t"
    cat ${i} |
      grep -E "(overall alignment rate)|(Overall time)" |
      perl -n -e '
        if(m/alignment/){
          $hash{percent} = $1 if m/([\d.]+)%/;
        }elsif(m/time/){
          if(m/(\d\d):(\d\d):(\d\d)/){
            my $time = $1 * 60 + $2 + $3 / 60;
            $hash{time} = $time;
          }
        }
        END{
          $hash{percent} = "NA" if not exists $hash{precent};
          $hash{time} = "NA" if not exists $hash{time};
          printf "%.2f\t%.2f\n", $hash{precent}, $hash{time};
        }
      '
done
```
解释：
`file_list=($(ls *.log))`定义**数组变量**   
`ehco -e`中`-e`enable启动转义字符的解析，识别换行\n、制表符\t等特殊字符    
`echo -n`中`-n`do not output the trailing newline不输出结尾的换行符   
`-grep -E`中`-E`extended扩展正则表达式   
`(overall alignment rate)|(Overall time)`中，`(...)`用于分组，将表达式分为不同部分，两个匹配模式   
perl里面片段：读取比对率和比对时间，并控制输出结果为两位小数
    `$hash{percent} = $1 if m/([\d.]+)%/`哈希中，通过`percent`这个键向哈希中添加/更新值，`$hash{percent}`被赋值为`$1`    
    `m/([\d.]+)%/`，`([\d.]+)`捕获组`$1`，匹配一个或多个`+`数字/小数点`[\d.]`，即一个有小数的数字。`%`匹配百分号符号，数字后接百分号才匹配成功。   
    `m/(\d\d):(\d\d):(\d\d)/`，有三个两位数数字的字符串，`hh:mm:ss`，时间格式。`(\d\d)`两个数字，`:`分隔   
    `$hash{percent} = "NA" if not exists $hash{precent}`检查hash中是否存在percent键，不存在赋值NA，确保程序不出错   
    `%.2f`占位符，格式化输出，指定输出的数字格式。`%`格式控制符的开始标志，`.2`控制小数点后的位数为2位，`f`表示格式化输出的值是浮动点数float   

- 格式转换和排序
sam格式存储核酸比对结果，将sam压缩得到bam或cram格式
```
cd ~/project/rat/output/align

# 排序并转bam格式，再建立索引
parallel -k -j 4 "
    samtools sort -@ 4 {1}.sam > {1}.sort.bam;
    samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')

rm *.sam
ls
```
`-k`keep same order。`-@ 4`4线程

# 7.表达量统计
使用`HTseq-count`判断read属于哪个基因，**union**、intersection_strict、intersection_nonempty三种模型，一般用union   
```
htseq-count [options] <比对后的sam/bam文件> <注释的gff文件>
```
```
cd ~/project/rat/output
mkdir HTseq

cd align
parallel -j 4 "
    htseq-count -s no -f bam {1}.sort.bam ../../annotation/rn6.gff \
        >../HTseq/{1}.count 2>../HTseq/{1}.log
" ::: $( ls *.sort.bam | perl -p -e 's/\.sort\.bam$//' )
```

# 8.合并表达矩阵、标准化
## 8.1.合并
使用`R`中merge合并：
```
rm(list=ls())
setwd("~/project/rat/output/HTseq")

# 得到文件样本编号，提取编号和对应文件
files <- list.files(".","*.count")
f_lists <- list()
for(i in files){
    prefix = gsub("(_\\w+)?\\.count", "", i, perl=TRUE)
    f_lists[[prefix]] = i
}

id_list <- names(f_lists)
data <- list()
count <- 0
for(i in id_list){
    count <- count + 1
    a <- read.table(f_lists[[i]], sep="\t", col.names = c("gene_id",i))
    data[[count]] <- a
}

# 合并文件
data_merge <- data[[1]]
for(i in seq(2, length(id_list))){
    data_merge <- merge(data_merge, data[[1]],by="gene_id")
}

write.csv(data_merge, "merge.csv", quote = FALSE, row.name = FALSE)
```
解释
- `rm(list=ls())`，`ls()`函数，列出当前R环境中所有对象名称（返回**字符向量**）。`rm()`函数，接受**对象名称**作为参数并删除。
    `list=ls()`将ls()返回的对象列表作为list参数传给rm()进行删除
- `setwd("~/project/rat/output/HTseq")`，设置工作目录
- `files <- list.files(".", "*.count")`
    `list.files(路径, 文件名匹配模式)`函数，列出某个路径下的文件。`<-`list.files函数结果赋值给`files`变量
- `f_lists <- list()`，`list()`函数，创建列表（数据结构），可以包含不同类型元素（数字、字符、数据框、矩阵等）
- `prefix = gsub("(_\\w+)?\\.count", "", i, perl=TRUE)`，提取前缀
    `gsub()`函数，字符串替换。gsub(匹配模式, 替换后的内容, 要处理的字符串, 是否使用perl风格正则表达式)。
    `(_\\w+)?\\.count`， _xxx.count或者.count
    `\\w`任何字母/数字/字符，`+`匹配一个或多个字符，`(_\\w+)`一个下划线后跟一个/多个字母数字字符，`?`匹配0次/1次。`\\.`转义。
- `f_lists[[prefix]] = i`，`[[..]]`访问列表中的元素。将提取到的前缀`prefix`作为键，文件名`i`作为值，存储到列表`f_lists`中。
- `id_list <- names(f_lists)`，`names()`函数，提取或设置对象名称（键，prefix）。
- `count <- count + 1`,`count`记录当前是第几次循环，在循环外必须**先初始化为0**。
- `a <- read.table(f_lists[[i]], sep="\t", col.names = c("gene_id",i))`，提取前面align对比完的每个文件内的结果
    `read.table()`函数，读取文件并将文件内容加载为数据框。`sep="\t"`指定文件中的列有制表符分隔。`col.names`指定文件列名，gene_id和i
- `data[[count]] <- a`，`data[[count]]`访问data列表中第count个位置，并用a赋值添加。
- `data_merge <- data[[1]]`，初始化一个合并结果的数据框`data_merge`
- `seq(2, length(id_list))`，`seq()`函数，生成序列，从2到`id_list`的长度的整数序列。`length()`函数，统计长度。
- `merge`函数，合并数据框，`by`指定基于哪个列合并。

## 8.2.标准化
[http://www.360doc.com/content/18/0112/02/50153987_721216719.shtml]
RPKM：基因counts数/总reads数/基因长度。单端测序，先测序深度标准化，再基因长度标准化   
FPKM：双端测序   
TPM：先基因长度标准化，再测序深度标准化。每个样本TPM总和相同，可样本间比较   
CPM：基因counts数/总reads数   
统计基因长度
```
library(GenomicFeatures)
# 构建Granges对象
txdb <- txdbmaker::makeTxDbFromGFF("rn6.gff" )
# 查找基因的外显子
exons_gene <- exonsBy(txdb, by = "gene")
# 计算总长度
# reduce()、width()是Irange对象的方法
gene_len <- list()
for(i in names(exons_gene)){
    range_info = reduce(exons_gene[[i]])
    width_info = width(range_info)
    sum_len    = sum(width_info)
    gene_len[[i]] = sum_len
}

# 或者写为lapply的形式(快很多)
gene_len <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

data <- t(as.data.frame(gene_len))
# 写入文件
write.table(data, file = "rn6_gene_len.tsv", row.names = TRUE, sep="\t", quote = FALSE, col.names = FALSE)
```
- CPM
```
CPM=10^6*基因counts数/比对到基因组的总reads数
```
- RPKM
```
RPKM=10^6*基因counts数/比对到基因组的总reads数/基因长度
# 基因长度单位kb，外显子长度/1000
```
```
gene_len_file <- "rn6_gene_len.tsv"
count_file <- "samples.count"

gene_len <- read.table(gene_len_file, header = FALSE, row.name = 1)
colnames(gene_len) <- c("length")

count <- read.table(count_file, header = FALSE, row.name = 1)
colnames(count) <- c("count")
# all read number
all_count <- sum(count["count"])

RPKM <- c()
for(i in row.names(count)){
    count_ = 0
    exon_kb = 1
    rpkm = 0
    count_ = count[i, ]
    exon_kb  = gene_len[i, ] / 1000
    rpkm    = (10 ^ 6 * count_ ) / (exon_kb * all_count )
    RPKM = c(RPKM, rpkm)
}
```
- TPM
```
TPM=（10^6*基因reads数/外显子长度之和）/ ∑(ni / gi)
∑(ni / gi)每个基因（reads/外显子长度和）之和
```
```
# 首先得到总的结果
sum_ <- 0
for(i in row.names(count)){
    count_ = 0
    exon = 1
    count_ = count[i, ]
    exon  = gene_len[i, ]
    value = count_ / exon
    if(is.na(value)){
        print(paste(i, " is error! please check"))
    }else{
        sum_ = sum_ + value
    }
}

TPM <- c()
for(i in row.names(count)){
    count_ = 0
    exon = 1
    count_ = count[i, ]
    exon  = gene_len[i, ]
    tpm = (10 ^ 6 * count_ / exon ) / sum_
    TPM = c(TPM, tpm)
}

count["RPKM"] <- RPKM
count["TPM"] <- TPM       
           
write.table(count, "123.normalize.count", col.names = TRUE, row.names = TRUE, sep="\t", quote = FALSE)
```

# 9.差异表达分析
## 9.1.数据前处理
### 删除HTseq-count结果的总结行:
```
dataframe <- read.csv("merge.csv", header=TRUE, row.names = 1)
# 删除前五行
countdata <- dataframe[-(1:5),]

# 查看数据
head(countdata)
```
### 删除基因id的版本号
基因名后面会有`.1`、`.2`标号，需删除
```
# 得到行的名
row_names <- row.names(countdata)

# 开始替换
name_replace <- gsub("\\.\\w+","", row.names(countdata))

row.names(countdata) <- name_replace
```
### 去除低表达的基因
```
countdata <- countdata[rowSums(countdata) > 0,]
```
## 9.2.差异分析
### 安装R包
- R包
```
# 使用bioconductor安装
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# R包
BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
BiocManager::install("biomaRt")
BiocManager::install("org.Rn.eg.db")
BiocManager::install("clusterProfiler")

# 加载
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Rn.eg.db)
library(clusterProfiler)
```
### 构建对象
用法：
```
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~ batch + condition)
```
测序样本信息：
```
cat <<EOF >./phenotype/phenotype.csv
"ids","state","condition","treatment"
"SRR2240185","Liver cirrhosis","DEN","treatment"
"SRR2240186","Liver cirrhosis","DEN","treatment"
"SRR2240187","Healthy control","PBS","control"
"SRR2240228","Healthy control","PBS","control"
EOF
```
导入R中
```
# 刚才countdata已经得到
countdata

# 读取样本分组信息(注意，需要加上row.names = 1, header = TRUE，将行列名需要看好)
coldata <- read.table("../phenotype/phenotype.csv", row.names = 1, header = TRUE, sep = "," )
# 确认一下行列名是否有（不是简单的数值）
head(coldata)
# 调整数据顺序
countdata <- countdata[row.names(coldata)]

# 构建dds对象
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)

# 查看dds
dds
```
### 样本相关性
- PCA分析
(解读)[http://blog.genesino.com/2016/10/PCA/]
```
# 接续着上面的构建得到的dds对象
# DEseq2包提供了相应的函数、归一
vsdata <- rlog(dds, blind=FALSE)
# intgroup 分组
plotPCA(vsdata, intgroup="treatment") + ylim(-10, 10)
```
- sample-to-sample distances热图
```
# 颜色管理包（不是必须）
library("RColorBrewer")
# 得到数据对象中基因的计数的转化值
gene_data_transform <- assay(vsdata)
# 使用t()进行转置
# 使用dist方法求样本之间的距离
sampleDists <- dist(t(gene_data_transform))
# 转化为矩阵用于后续pheatmap()方法的输入
sampleDistMatrix <- as.matrix(sampleDists)
# 将矩阵的名称进行修改
# rownames(sampleDistMatrix) <- paste(vsdata$treatment, vsdata$condition, vsdata$ids, sep="-")
# colnames(sampleDistMatrix) <- paste(vsdata$treatment, vsdata$condition, vsdata$ids, sep="-")
# 设置色盘
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# 绘制热图与聚类
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```
### 差异基因
- `DESeq()`计算不同组别间的基因的表达差异，输入是之前构建的`dds`对象：
```
# 改变样本组别顺序
dds$treatment <- factor(as.vector(dds$treatment), levels = c("control","treatment"))

# 基于统计学方法进行计算
dds <- DESeq(dds)

# 查看实验组和对照组对比结果
result <- results(dds, pAdjustMethod = "fdr", alpha = 0.05)
head(result)

# 结果按照p-value排序
result_order <- result[order(result$pvalue),]
head(result_order)
```
输出：
```
log2 fold change (MLE): treatment treatment vs control 
Wald test p-value: treatment treatment vs control 
DataFrame with 6 rows and 6 columns
                    baseMean log2FoldChange     lfcSE      stat       pvalue         padj
                   <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
ENSRNOG00000061890  1795.683       -4.27345  0.188101  -22.7189 2.91415e-114 5.95070e-111
ENSRNOG00000054181   969.487       -5.32981  0.242085  -22.0163 2.01216e-107 2.05441e-104
ENSRNOG00000018086  1878.796       -3.59390  0.190641  -18.8516  2.85134e-79  1.94081e-76
ENSRNOG00000013552 30917.241       -3.89394  0.218504  -17.8209  4.86172e-71  2.48191e-68
ENSRNOG00000016807  2598.426       -3.17581  0.183086  -17.3460  2.11439e-67  8.63516e-65
ENSRNOG00000020560   517.413        4.61549  0.273691   16.8639  8.29356e-64  2.82258e-61
```
其中`log2 fold change (MLE): treatment treatment vs control `这行很重要    
- 总结基因上下调情况
```
summary(result_order)
```
输出：
```
out of 2473 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 340, 14%
LFC < 0 (down)     : 317, 13%
outliers [1]       : 0, 0%
low counts [2]     : 431, 17%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```
- 查看显著的基因数量
```
table(result_order$padj<0.05)
```
- 保存数据
```
# 新建文件夹
dir.create("../DESeq2")
# 不用按照padj排序的结果，就保存按照基因名排序的
write.csv(result, file="../DESeq2/results.csv", quote = F)
```
# 10.提取差异表达基因与注释
## 10.1.提取差异基因
```
# padj 小于 0.05 并且 Log2FC 大于 1 或者小于 -1
diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)

# 查看数据框的大小
dim(diff_gene)

# 把差异基因写入文件
dir.create("../DESeq2/")
write.csv(diff_gene, file="../DESeq2/difference.csv", quote = F)
```
## 10.2.转化基因ID
使用`ClusterProfiler`
```
# 安装clusterProfiler包
BiocManager::install("clusterProfiler")
# 这里我们分析的是大鼠，安装大鼠的数据库
BiocManager::install("org.Rn.eg.db")

# 加载包
library(clusterProfiler)
library(org.Rn.eg.db)

# 得到基因ID(这个ID是Ensembl数据库的编号)
ensembl_gene_id <- row.names(diff_gene)

# 转换函数
ensembl_id_transform <- function(ENSEMBL_ID){
    # geneID是输入的基因ID，fromType是输入的ID类型，toType是输出的ID类型，OrgDb注释的db文件，drop表示是否剔除NA数据
    a = bitr(ENSEMBL_ID, fromType="ENSEMBL", toType=c("SYMBOL","ENTREZID"), OrgDb="org.Rn.eg.db")
    return(a)
}

# 开始转化
ensembl_id_transform(ensembl_gene_id)
```
## 10.3.注释
使用`biomaRt`
```
BiocManager::install("biomaRt")
library("biomaRt")

# 选择数据库
mart <- useDataset("rnorvegicus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))

# 得到基因ID(这个ID是Ensembl数据库的编号)
ensembl_gene_id <- row.names(diff_gene)
rat_symbols <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"), filters = 'ensembl_gene_id', values = ensembl_gene_id, mart = mart)
```
- 把基因矩阵和`symbols`合并数据框
```
# 生成用于合并的列
diff_gene$ensembl_gene_id <- ensembl_gene_id
# 将DESeq2对象转换为数据库
diff_gene_dataframe <- as.data.frame(diff_gene)
# 合并
diff_gene_symbols <- merge(diff_gene_dataframe, rat_symbols, by = c("ensembl_gene_id"))
```
- 保存结果
```
write.table(result, "../stat/all_gene.tsv", sep="\t", quote = FALSE)
write.table(diff_gene_symbols, "../stat/diff_gene.tsv", row.names = F,sep="\t", quote = FALSE)
```
- 统计样本的差异基因***???**
```bash
cd ~/project/rat/output/stat
echo -e "sample\tnum" > all_samples.tsv
for i in $(ls);
do
    if [ -d ${i} ];
    then
        prefix=$i
        diff_num=$(cat $i/diff_gene.tsv | tail -n+2 | wc -l)
        echo -e "${prefix}\t${diff_num}" >> all_samples.tsv
    fi
done
```

# 11.可视化
- MA图
```
plotMA(result_order, ylim=c(-10,10))
```
- 热图

# 12.富集分析
使用`clusterProfiler`中`enrichGO`
```R
# 接续着上面的结果
ensembl_gene_id <- row.names(diff_gene)

# 得到symbol
rat_symbols <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"), filters = 'ensembl_gene_id', values = ensembl_gene_id, mart = mart)
diff_gene_ensembl_id <- rat_symbols$ensembl_gene_id
```
`enrichGO`参数：
```
gene	        差异基因对应的向量
keyType	        指定的gene的ID类型，一般都用ENTREZID，该参数的取值可以参考keytypes(org.Hs.eg.db)的结果
OrgDb	        该物种对应的org包的名字
ont	            代表GO的3大类别，BP, CC, MF
pAdjustMethod	指定多重假设检验矫正的方法
pvalueCutoff	对应的阈值
qvalueCutoff	对应的阈值
```
### GO分析：
```
for(i in c("MF", "BP", "CC")){
    ego <- enrichGO(gene       = rat_symbols$ensembl_gene_id,
                    OrgDb      = org.Rn.eg.db,
                    keyType    = 'ENSEMBL',
                    ont        = i,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05)
    dotplot(ego, showCategory = 30, title = paste("The GO ", i, " enrichment analysis", sep = ""))
}
```

### KEGG分析
```
kk <- enrichKEGG(gene = rat_symbols$entrezgene_id, 
                 organism ='rno',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data = FALSE)

dotplot(kk, showCategory = 30, title = paste("The KEGG ", " enrichment analysis", sep = ""))
```

### GSEA分析
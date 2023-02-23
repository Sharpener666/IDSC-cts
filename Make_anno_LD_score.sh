####计算分层LD分数####
#Step 1: Creating an annot file#
python make_annot.py \
		--gene-set-file GTEx_Cortex.GeneSet \#前排10%富集上调的geneset
		--gene-coord-file ENSG_coord.txt \ #基因坐标文件
		--windowsize 100000 \ #扩大基因区域前后100K的区域
		--bimfile 1000G.EUR.QC.22.bim \ #bim文件计算富集基因的LD分数
		--annot-file GTEx_Cortex.annot.gz #最后输出的注释文件

#Step 2: Computing LD scores with an annot file#
python ldsc.py\
		--l2\ #计算ld分数
		--bfile 1000G.EUR.QC.22\  #指定plink文件路径及前缀
		--ld-wind-cm 1\ #指定窗口的大小，单位为cM
		--annot Brain_DPC_H3K27ac.annot.gz\ #前面生成的注释文件
		--thin-annot
		--out Brain_DPC_H3K27ac\
		--print-snps hm.22.snp

#conda create -n iannot
#conda activate iannot
#conda install python=3.10
#pip install iAnnotateSV
#pip install pysam
#git clone https://github.com/rhshah/iAnnotateSV.git
#manually download files from github.com/rhshah/iAnnotateSV/iAnnotateSV/data

module load conda3
source activate res
source config.cfg
#source activate iannot

#python prep_input.py TB-22-17692_S4_extracted.tsv
#python iAnnotateSV/iAnnotateSV.py -r hg19 -ofp TB-22-17692_S4_out -o /home/mjani/aml/iAnnotateSV -i TB-22-17692_S4_input.tsv

mkdir ${outdir}
#echo ${basedir}SV_output/24P-010IP0021*.sv.vcf
#if false; then
for d in ${basedir}SV_output/*.sv.vcf; do
 name=$(basename $d .sv.vcf)
 cmd="python extract.py $d ${svbed} ${outdir}/${name}_extracted.tsv"
 echo $cmd
 bsub -M 10000 -J ${name}_ex -eo logs/${name}_ex.err -oo logs/${name}_ex.oo.out "$cmd"

 cmd="python prep_input.py ${outdir}/${name}_extracted.tsv ${outdir}/"
 echo $cmd
 bsub -M 10000 -J ${name}_prep -w "done(${name}_ex)" -eo logs/${name}_prep.err -oo logs/${name}_prep.oo.out "$cmd"

 cmd="python ${src} -r hg19 -ofp ${name} -o ${outdir} -i ${outdir}/${name}_input.tsv"
 echo $cmd
 bsub -M 10000 -J ${name}_ann -w "done(${name}_prep)" -eo logs/${name}.eo.err -oo logs/${name}.oo.out "$cmd"

 cmd="python annotate_iann.py ${outdir}/${name}_Annotated.txt ${outdir}/${name}_extracted.tsv ${basedir}varcall/${name}/${name}_tumor.bam ${outdir}/${name}_annotated.group.tsv"
 bsub -M 10000 -J ${name}_iann -w "done(${name}_ann)" -eo logs/eoann.err -oo logs/ooann.out "$cmd"
 echo $cmd
 #break
done
#fi

agfusion(){
 outdir=$1
 name=$2
 cmd="bash run_agfusion.sh ${outdir}/${name}_annotated.group.tsv ${outdir} ${name} $3 $4"
 echo $cmd
 bsub -M 90000 -J ${name}_gen -w "done(${name}_iann)" -eo logs/${name}.sag.err -oo logs/${name}.sag.out "$cmd"
 #bsub -M 90000 -J ${name}_gen -eo logs/${name}.sag.err -oo logs/${name}.sag.out "$cmd"
}

for d in ${basedir}SV_output/*.sv.vcf; do
  name=$(basename $d .sv.vcf)
  agfusion ${outdir} ${name} ${hgref} ${mittledb}
  #break
done

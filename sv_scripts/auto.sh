#conda create -n iannot
#conda activate iannot
#conda instal python=3.10
#pip install iAnnotateSV
#pip install pysam
#git clone https://github.com/rhshah/iAnnotateSV.git
#manually download files from github.com/rhshah/iAnnotateSV/iAnnotateSV/data

#module load conda3
#source activate iannot

#python prep_input.py TB-22-17692_S4_extracted.tsv
#python iAnnotateSV/iAnnotateSV.py -r hg19 -ofp TB-22-17692_S4_out -o /home/mjani/aml/iAnnotateSV -i TB-22-17692_S4_input.tsv

if false; then
basedir="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/"
outdir="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/SV_output"
fi

basedir="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/normals/cfdna/Output_combined/"
outdir="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/normals/cfdna/Output_combined/SV_output/annotated"

src="/home/mjani/aml/iAnnotateSV/iAnnotateSV/iAnnotateSV.py"
svbed="/clinical/ccs01/home/clusterHome/mjani/aml/scripts/aml_mrd_sv_genes.bed"
hgref="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/SV_output/self_annot/hg19.refGene.gtf"
mittledb="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/SV_output/combined_scripts/mittleman_genes.tsv"

if false; then
for d in ${basedir}SV_output/*.sv.vcf; do
 name=$(basename $d .sv.vcf)
 cmd="python extract.py $d ${svbed} ${outdir}/${name}_extracted.tsv"
 $cmd
 cmd="python prep_input.py ${outdir}/${name}_extracted.tsv ${outdir}/"
 $cmd
 cmd="python ${src} -r hg19 -ofp ${name} -o ${outdir} -i ${outdir}/${name}_input.tsv"
 bsub -M 1000 -J ${name}_ann -eo logs/eo.err -oo logs/oo.out "$cmd"
 #echo $cmd
 cmd="python annotate_iann.py ${outdir}/${name}_Annotated.txt ${outdir}/${name}_extracted.tsv ${basedir}varcall/${name}/${name}_tumor.bam ${outdir}/${name}_annotated.group.tsv"
 bsub -M 1000 -J ${name}_iann -w "done(${name}_ann)" -eo logs/eoann.err -oo logs/ooann.out "$cmd"
 #bsub -M 1000 -J ${name}_iann -eo logs/eoann.err -oo logs/ooann.out "$cmd"
done
fi

agfusion(){
 outdir=$1
 name=$2
 source activate agfusion
 cmd="python driver_post_annotate.py ${outdir}/${name}_annotated.group.tsv ${outdir} ${name} $3 $4"
 echo $cmd
 bsub -M 90000 -J ${name}_gen -eo logs/${name}.sag.err -oo logs/${name}.sag.out "$cmd"
}

for d in ${basedir}SV_output/*.sv.vcf; do
 name=$(basename $d .sv.vcf)
 agfusion ${outdir} ${name} ${hgref} ${mittledb}
done

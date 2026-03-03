outdir="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/normals/cfdna/dragen_systematic_noise"
source="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/normals/cfdna/Output_combined"
bed="/home/mjani/aml/bedfile_SJPedPan_combined/Annotated_Target_MRD_Panel.bed"
LOG="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/normals/cfdna/dragen_systematic_noise/logs"
dragen_ref_dir="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Validation/ref/ref_dragen42"
mkdir ${LOG}

for bam in ${source}/varcall/*/*.bam; do
 pref=$(basename $bam _tumor.bam)
 mkdir "${outdir}/${pref}"
 cmd="dragen \
 --ref-dir ${dragen_ref_dir} \
 --tumor-bam-input ${bam} \
 --output-directory ${outdir}/${pref} \
 --output-file-prefix ${pref} \
 --enable-map-align false \
 --enable-variant-annotation true \
 --variant-annotation-data /staging/Data/ \
 --variant-annotation-assembly GRCh37 \
 --enable-maf-output false \
 --enable-variant-caller true \
 --vc-detect-systematic-noise true \
 --vc-enable-germline-tagging true \
 --vc-enable-umi-liquid true \
 --vc-target-bed ${bed} \
 --vc-enable-triallelic-filter false"
 echo $cmd
 bsub -app dragen-500g -q dragen_test -J ${pref} -eo $LOG/${pref}.eo.vc.txt -oo $LOG/${pref}.oo.vc.txt "$cmd"
 #break
done

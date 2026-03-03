dragen_ref_dir="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Validation/ref/ref_dragen42"
outdir="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/normals/cfdna/dragen_systematic_noise"
normvcf="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/normals/cfdna/dragen_systematic_noise/normal_vcfs.txt"
outfile="aml_mrd_cfdna_systematic_noise.bed.gz"
LOG="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/normals/cfdna/dragen_systematic_noise/logs"

cmd="dragen \
 --ref-dir ${dragen_ref_dir} \
 --output-directory ${outdir} \
 --output-file-prefix aml_mrd_cfdna_systematic_noise \
 --enable-map-align false \
 --build-sys-noise-vcfs-list ${normvcf} \
 --build-sys-noise-use-germline-tag true \
 --build-sys-noise-min-supporting-samples 2"

echo $cmd
bsub -app dragen-500g -q dragen_dev -J map -eo $LOG/eo.txt -oo $LOG/oo.txt "$cmd";

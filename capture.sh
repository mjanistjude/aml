module load conda3

bash runqc.sh capture
source activate res

source config.cfg

sample=$(basename "source")
LOG="${output}/piplogs"

if [[ ${sample_type} == "TWIST" ]]; then
 ctstart=7
 blen=5
 whitelist=${twist_whitelist}
else
 ctstart=8
 blen=8
 whitelist=${idt_whitelist}
fi


OUTDIR=${output}/varcall
temp="sample_name"
fastq="Temp_fastq_name"
svout=${output}/SV_output


cmd2="dragen -f \
--ref-dir ${dragen_ref_dir} \
--tumor-fastq1 ${source}/${fastq}_R1_001.fastq.gz \
--tumor-fastq2 ${source}/${fastq}_R2_001.fastq.gz \
--output-directory ${OUTDIR}/${temp} \
--output-file-prefix ${temp} \
--RGSM-tumor ${temp} \
--RGID-tumor ${temp} \
--enable-map-align true \
--enable-umi true \
--umi-source qname \
--umi-library-type nonrandom-duplex \
--umi-nonrandom-whitelist ${whitelist} \
--umi-min-supporting-reads 2 \
--umi-emit-multiplicity both"

cmd3="dragen \
--ref-dir ${dragen_ref_dir} \
--tumor-bam-input BAMFILE \
--output-directory OUTPUTDIR \
--output-file-prefix OUTPUTPRE \
--enable-map-align false \
--enable-variant-annotation true \
--variant-annotation-data /staging/Data/ \
--variant-annotation-assembly GRCh37 \
--enable-maf-output true \
--maf-transcript-source RefSeq \
--vc-systematic-noise ${sysnoise} \
--enable-variant-caller true \
--vc-enable-germline-tagging true  \
--vc-enable-umi-liquid true \
--vc-target-bed ${bed} \
--vc-target-vaf 0.001 \
--vc-enable-triallelic-filter false"


setup_dir() {
 mkdir ${output}
 mkdir ${LOG}
 mkdir ${OUTDIR}
 mkdir ${OUTDIR}/vc_calls
 mkdir ${output}/SV_output
}

call_sv() {
  mate=$(basename "$1" _R1_001.fastq.gz)
  cmd="java -jar ${apjar} call -1 $4 -1BL ${blen} -1BS 0 -1S ${ctstart} -2 ${source}/${mate}_R2_001.fastq.gz -2BL ${blen} -2BS 0 -2S ${ctstart} -D ${svout} -I ${apertureRef}/aperture_hg19 -P $3 -T 300"
  echo $cmd
  bsub -J SV${3} -M 50000 -eo ${LOG}/${1}.sv.err -oo ${LOG}/${1}.sv.out "$cmd"
  cmd="gunzip ${svout}/${3}.sv.vcf.gz"
  bsub -J SVgz${3} -w "done(SV${3})" -M 100000 -eo ${LOG}/${1}.sv.err -oo ${LOG}/${1}.sv.out "$cmd"
}

dragaln() {
  cmd="${cmd2//$temp/$3}"
  cmd="${cmd//$fastq/$3}"
  mkdir ${OUTDIR}/${3}
  echo $cmd
  bsub -app dragen-500g -q dragen_dev -J map${3} -eo $LOG/${3}.eo.map.txt -oo $LOG/${3}.oo.map.txt "$cmd";
}

dragvc() {
  bam=${OUTDIR}/${3}/${3}_tumor.bam
  vcout=${OUTDIR}/vc_calls/$3
  cmd="${cmd3//BAMFILE/$bam}"
  cmd="${cmd//OUTPUTDIR/$vcout}"
  cmd="${cmd//OUTPUTPRE/$3}"

  mkdir -p ${vcout}
  echo $cmd
  bsub -app dragen-500g -q dragen_dev -J dvc${3} -w "done(map${3})" -eo $LOG/${3}.eo.vc.txt -oo $LOG/${3}.oo.vc.txt "$cmd"

  cmd="gunzip ${vcout}/*.gz"
  bsub -M 10000 -J gz${3} -w "done(dvc${3})" -eo $LOG/${3}.eo.vcgz.txt -oo $LOG/${3}.oo.vcgz.txt "$cmd"
}

annotate() {
 cmd="python gnomad_annot.py $1"
 echo $cmd
 bsub -M 10000 -J annot${3} -w "done(gz${3})" -eo $LOG/${3}.eo.annot.txt -oo $LOG/${3}.oo.annot.txt "$cmd"
}

annotsv(){
 cmd="python extract.py ${svout}/${3}.sv.vcf"
 echo $cmd
 bsub -M 100000 -J extSV${3} -w "done(SVgz${3})" -eo $LOG/${3}.eo.extsv.txt -oo $LOG/${3}.oo.extsv.txt "$cmd"

 cmd="bash annot_sv.sh ${svout}/${3}_extracted.vcf ${svout} ${3}_annotated.tsv ${3}"
 echo $cmd
 bsub -M 100000 -J stSV${3} -w "done(extSV${3})" -eo $LOG/${3}.eo.runannotsv.txt -oo $LOG/${3}.oo.runannotsv.txt "$cmd"
}

postqc(){
 cmd="python mappingqc.py ${OUTDIR} ${output}/QC ${bed}"
 echo $cmd
 echo "bsub -M 10000 -J mapqc -w "$1""
 echo $1
 bsub -M 10000 -J mapqc -w "$1" -eo $LOG/mapqc.eo.txt -oo $LOG/mapqc.oo.txt "$cmd"
}

#if false;
#then
setup_dir
jobid=""
for f in "$source"/*.fastq.gz; do
  fname=$(basename "$f")
  outname=${fname//.fastq.gz/}
  bname=$(basename "$f" _R1_001.fastq.gz)

  if [[ $fname == *"Undetermined"* || $fname == *"WGS"* || $fname == *"R2_001.fastq.gz"* ]]; then
   continue
  fi
  call_sv ${fname} ${outname} ${bname} $f
  dragaln ${fname} ${outname} ${bname}
  dragvc ${fname} ${outname} ${bname}
  annotate ${OUTDIR}/vc_calls/${bname} ${outname} ${bname}
  jobid=${jobid}"done(gz${bname}) && "
  #break
done
#fi

jobid="${jobid:0:-4}"
#jobid="${jobid//(/(\"}"
#jobid="${jobid//)/\")}"
echo $jobid
postqc "$jobid"

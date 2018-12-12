# main work dir
export project=/projects/nick_matthews
export scripts=${project}/chlamy_locus_map_github/Key_Coding_Scripts/

# linking relvant bam files from pipeline into project folder
bash ${scripts}/generate_bw_files.sh

# generate bw from bams:
cd ${project}/bams
echo $project
for file in $(ls *.bam); do
  echo $file
  bash ${scripts}/bam2bw.sh $file ${project}/resources/Creinhardtii_236.fa.chrsizes sort
done

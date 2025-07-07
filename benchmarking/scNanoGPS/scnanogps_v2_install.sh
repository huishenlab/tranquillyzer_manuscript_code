
#  Python 3.10

# May 28, 2025


#--------------- Notes
# v1 & v2 requires 3.11; does not work with >3.11
# On 5/28/25, github for scnanogps says to use pysam=0.19.0, but this requires python >=3.9,<3.10.0a0

#--------------- End Notes

#~ https://github.com/gaolabtools/scNanoGPS?tab=readme-ov-file#installation
#~ conda create -n scNanoGPS_v2 python=3.11 numpy scipy
#                                                                                                                                                                                                                  
# To activate this environment, use                                                                                                                                                                                
#                                                                                                                                                                                                                  
#     $ conda activate scNanoGPS_v2                                                                                                                                                                                   
#                                                                                                                                                                                                                  
# To deactivate an active environment, use                                                                                                                                                                         
#                                                                                                                                                                                                                  
#     $ conda deactivate  

# then install more dependencies
conda install -c conda-forge biopython=1.80 # this requires python >=3.10,<3.11.0a0
conda install -c conda-forge distance=0.1.3
conda install -c conda-forge matplotlib=3.8.2
conda install -c conda-forge pandas=2.1.4
#~ conda install -c bioconda pysam=0.19.0 # this is where Python>3.11 gives depend error with v1 (was 0.23.0 with v1 so v2 uses older version)
conda install -c bioconda pysam=0.23.0 # this is where Python>3.11 gives depend error with v1 (was 0.23.0 with v1 so v2 uses older version)
conda install -c bioconda seaborn=0.13.1 # switching down to 0.13.1 from 0.13.2 in v1

#~ git clone https://github.com/gaolabtools/scNanoGPS/
#~ cd scNanoGPS
#~ pip3 install -r requirements.txt # ran this despite installing most of these, no version numbers specified here

# install other tools
conda install -c bioconda minimap2
conda install -c bioconda samtools
conda install -c bioconda tabix
conda install -c bioconda spoa
conda install -c bioconda subread
conda install -c bioconda longshot
conda install -c conda-forge -c bioconda isoquant # this is new to v2 and not working b/c version conflicts... :)
conda install -c bioconda bcftools

# downloaded annovar, then tar -xvf annovar.latest.tar.gz
# then add this to my path
export PATH="2025_03_12_tranquilizer_benchmarking/scnanogps/annovar:$PATH"
# and qualimap, wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip
export PATH="2025_03_12_tranquilizer_benchmarking/scnanogps/qualimap_v2.3:$PATH"

# linked in ENCODE references from the BBC
# build index for minimap2 specific for 10x3p
#~ minimap2 -x splice -d hg38_gencode.mmi ../sequence/hg38_gencode.fa

# configure ANNOVAR
#~ perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene hg38db/
#~ perl annotate_variation.pl -buildver hg38 -downdb cytoBand hg38db/
#~ perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome hg38db/
#~ perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 hg38db/


The files in this repository represent the code used to 
process data and carry out analysis in "Predicting gene 
expression in the human malaria parasite Plasmodium falciparum"

Shell scripts run within the "dataPreProcessing" subdirectory used the command
echo "source ~/myenv.sh" >> $(a).job;
or similar commands. The file ~/myenv.sh contained instructions
for which software/versions to load before running the specified job
file. These instructions are specific to your computing environment.
Anytime "source ~/myenv.sh" is run, you need to have the following
software/versions loaded:
python 2.7.3
numpy 1.11.0
pandas 0.18.1
samtools 0.1.19
pysam 0.8.1
cython 0.25.2
matplotlib 1.5.1
bwa 0.7.3
fastqc 0.11.3
fastx-toolkit 0.0.14
pysam 0.8.4
curl 7.48.0

A number of the scripts contained here require reference datasets
that we have not directly included in this repository. To run
these scripts locally, you will need:

1) An indexed Plasmodium genome for use in BWA alignment
Use BWA version XXX
Replace all instances of /net/noble/vol1/data/bwa-indices/pfal3D7-PlasmoDBv29/pfal3D7-PlasmoDBv29
with the appropriate path to your local indexed Plasmodium genome.



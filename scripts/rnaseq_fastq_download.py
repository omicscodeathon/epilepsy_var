import subprocess

# samples correspond to 17 neocortex, 17 hippocampus, and 2 control without Epilepsy neocortex  
sra_numbers = [
    "SRR9733947","SRR9733948","SRR9733949","SRR9733950","SRR9733951","SRR9733952",\
        "SRR9733953","SRR9733954","SRR9733955","SRR9733956","SRR9733957","SRR9733958",\
            "SRR9733959","SRR9733960","SRR9733961","SRR9733962","SRR9733963","SRR9733964",\
                "SRR9733965","SRR9733966","SRR9733967","SRR9733968","SRR9733969","SRR9733970",\
                    "SRR9733971","SRR9733972","SRR9733973","SRR9733974","SRR9733975","SRR9733976",\
                        "SRR9733977","SRR9733978","SRR9733979","SRR9733980","SRR9733981","SRR9733982"
    ]

# This will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = "prefetch " + sra_id
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)

# This will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)
# Input files:
# ~/d/artic-4.1.fastq
# ~/d/midnight-v2.fastq
f := "~/d/artic-4.1.fastq"
test1:
    cargo run -r -- {{f}}
    

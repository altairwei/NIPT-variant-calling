rule Module_2_Calling_Step_1:
    input:
        bam=expand(os.path.join(OUTPUT_DIR, "alignments", 
            "{sample_id}.sorted.rmdup.realign.BQSR.bam"), sample_id=[s[3] for s in SAMPLES]),
        bai=expand(os.path.join(OUTPUT_DIR, "alignments",
            "{sample_id}.sorted.rmdup.realign.BQSR.bam.bai"), sample_id=[s[3] for s in SAMPLES])
    output:
        os.path.join(OUTPUT_DIR, "calls", "all.bam.list")
    run:
        with open(output[0], "w") as f:
            for bam_file in input.bam:
                f.write(bam_file + "\n")

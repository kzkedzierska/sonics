logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s %(levelname)s %(message)s"
)
block_name = "Block 21313"
sample = "sample_55454"

logging.warning(("Partial repetition of the motif in"
                 "block: %s sample: %s. Excluding"
                 "allele"), block_name, sample)

logging.info("Initiating simulation for %s %s", sample, block_name)

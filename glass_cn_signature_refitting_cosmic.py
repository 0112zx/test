

import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze

samples = '/boot3/bio_liaojl/GLASS_copy_number_signature_work/Data/Processed/for_SigProfiler/ASCAT_all440.CNV48.matrix.tsv'
output = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_CN_signature_refitting'
refit_sigs = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/glioma_cosmic_cnsig_definition.txt'

Analyze.denovo_fit( samples,
                    output, 
                    signatures=refit_sigs,
                    signature_database=None,
                    genome_build="GRCh37", 
                    verbose=False)








from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

# SBS96 matrix generation

project = 'Glioma_SBS'
genome = 'GRCh37'
vcfFiles = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler_SBS'

matrices = matGen.SigProfilerMatrixGeneratorFunc(project, genome, vcfFiles, exome=True, bed_file=None, chrom_based=False, plot=False, tsb_stat=False, seqInfo=False)

# COSMIC refitting

import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze

samples = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler_SBS/output/SBS/Glioma_SBS.SBS96.exome'
output = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/SBS_signature_refitting/'
sigs = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cosmic_sig_ref.txt'

Analyze.cosmic_fit( samples, 
                    output, 
                    signatures=None,
                    signature_database=sigs,
                    genome_build="GRCh37",
                    cosmic_version=3.3,
                    verbose=False,
                    collapse_to_SBS96=False,
                    make_plots=True,
                    exclude_signature_subgroups=None,
                    exome=True)





















import os
from SigProfilerExtractor import sigpro as sig
from SigProfilerExtractor import estimate_best_solution as ebs
from SigProfilerMatrixGenerator.scripts import CNVMatrixGenerator as scna

os.chdir('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification')
# os.getcwd()

################ ABSOLUTE

# Copy number matrix generation

file_type = "ASCAT"
output_path = "/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler"
input_file = "/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler/glioma_seg_ascat_form.tsv"
project = "tcga_glioma_ABSOLUTE"

scna.generateCNVMatrix(file_type, input_file, project, output_path)

# copy number signature extraction

gold_matrix = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler/tcga_glioma_ABSOLUTE.CNV48.matrix.tsv'
sig.sigProfilerExtractor("matrix", "sigprofiler_output", gold_matrix, maximum_signatures = 12)

# os.chdir('/home/xudahua/temp')
# gold_matrix = '/home/xudahua/tcga_glioma_ABSOLUTE.CNV48.matrix.tsv'
# sig.sigProfilerExtractor("matrix", "sigprofiler_output", gold_matrix, maximum_signatures = 12)


################ ASCAT

### glioma

# Copy number matrix generation

file_type = "ASCAT"
output_path = "/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler"
input_file = "/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler/glioma_cn_ascat.tsv"
project = "tcga_glioma_ASCAT"

scna.generateCNVMatrix(file_type, input_file, project, output_path)

# copy number signature extraction

gold_matrix = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler/tcga_glioma_ASCAT.CNV48.matrix.tsv'
sig.sigProfilerExtractor("matrix", "sigprofiler_output", gold_matrix, maximum_signatures = 12)

# os.chdir('/home/xudahua/temp')
# gold_matrix = '/home/xudahua/tcga_glioma_ASCAT.CNV48.matrix.tsv'
# sig.sigProfilerExtractor("matrix", "sigprofiler_output", gold_matrix, maximum_signatures = 12)


### LGG

# Copy number matrix generation

file_type = "ASCAT"
output_path = "/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler"
input_file = "/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler/lgg_cn_ascat.tsv"
project = "tcga_lgg_ASCAT"

scna.generateCNVMatrix(file_type, input_file, project, output_path)

# copy number signature extraction

gold_matrix = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler/tcga_lgg_ASCAT.CNV48.matrix.tsv'
sig.sigProfilerExtractor("matrix", "sigprofiler_output", gold_matrix, maximum_signatures = 25)

# os.chdir('/home/xudahua/temp')
# gold_matrix = '/home/xudahua/tcga_lgg_ASCAT.CNV48.matrix.tsv'
# sig.sigProfilerExtractor("matrix", "sigprofiler_output", gold_matrix, maximum_signatures = 25)



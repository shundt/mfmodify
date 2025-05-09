# IMPORT
import os
import zipfile

# FUNCTIONS
# def extract_if_needed(zip_file_path, extract_to=None):
#     # get an extract_to path if not specified 
#     if extract_to is None:
#         extract_to = os.path.dirname(zip_file_path)
#     basename = os.path.basename(zip_file_path)
#     filename, _ = os.path.splitext(basename)
#     extract_final_name = os.path.join(extract_to, filename)
#     # check if extracted
#     if not os.path.exists(extract_final_name):
#         print(f'Extracting {zip_file_path} to {extract_final_name}...')
#         with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
#             zip_ref.extractall(extract_to)
#     else:
#         print(f'{zip_file_path} already extracted to {extract_final_name}')
# Python tools for Interoperable Neuromorphometry Reporting
Python tools for Interoperable Neuromorphometry Reporting (`pinr`) allows for the creation of DICOM Structured Reporting Documents and FHIR Diagnostic Report resources from a T1w image processed by FreeSurfer.

## Installation
`pip install` is coming soon!

For now, installation requires [poetry](https://python-poetry.org/).
Run the following command to install:
```
poetry install
```
Requirements: python >= 3.8

## Usage
Currently, pinr exposes a class called `FreeSurferSeg`, which stores volumetric data and measurements.
A CLI is coming soon!

To initiate:
```python
t1w_dicom_file = "/path/to/file.dcm"
aseg_file = "/path/to/mri/aseg.mgz"
aseg_stats_file = "/path/to/aseg.stats"
# Choose where to save output files
dicom_seg_output_file = "/path/to/aseg.dcm"
dicom_sr_output_file = "/path/to/aseg_sr.dcm"
fhir_output_file = "/path/to/aseg.fhir.json"

aseg = FreeSurferSeg(
    t1w_dicom_file=t1w_dicom_file,
    aseg_file=aseg_file,
    aseg_stats_file=aseg_stats_file,
    dicom_seg_output_file=dicom_seg_output_file,
    dicom_sr_output_file=dicom_sr_output_file,
    fhir_output_file=fhir_output_file,
)
```
Data is stored as properties of the object, and is written to the given file name if one is provided.

```python
seg = aseg.seg
sr = aseg.sr
fdr = aseg.fhir
```


----
Work here was supported by the National Institute Of Biomedical Imaging And Bioengineering of the
National Institutes of Health under Award Number R43EB030910.
